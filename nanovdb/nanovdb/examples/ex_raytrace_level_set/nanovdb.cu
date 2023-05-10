// Copyright Contributors to the OpenVDB Project
// SPDX-License-Identifier: MPL-2.0

#include <cmath>
#include <chrono>

#include <nanovdb/util/IO.h>
#include <nanovdb/util/CudaDeviceBuffer.h>
#include <nanovdb/util/Ray.h>
#include <nanovdb/util/HDDA.h>
#include <nanovdb/util/GridBuilder.h>

#include "common.h"
#include <nanovdb/nanovtt/NanoVTT.h>
#include <nanovdb/nanovtt/util/GridHandle.h>
#include <nanovdb/nanovtt/util/HDDA.h>
#include <nanovdb/nanovtt/util/SampleFromVoxels.h>

#define PNANOVDB_C
#define PNANOVDB_HDDA
#include <nanovdb/nanovtt/PNanoVTT.h>

#if defined(NANOVDB_USE_CUDA)
using BufferT = nanovdb::CudaDeviceBuffer;
#else
using BufferT = nanovdb::HostBuffer;
#endif


void runNanoVDB(nanovdb::GridHandle<BufferT>& handle, int numCPUIterations, int numGPUIterations, std::string& imagePostFix, int width, int height, BufferT& imageBuffer)
{
    using GridT = nanovdb::FloatGrid;
    using CoordT = nanovdb::Coord;
    using RealT = float;
    using Vec3T = nanovdb::Vec3<RealT>;
    using RayT = nanovdb::Ray<RealT>;

    auto* h_grid = handle.grid<float>();
    if (!h_grid)
        throw std::runtime_error("GridHandle does not contain a valid host grid");

    std::cout << "NanoVDB buffer size = " << static_cast<float>(handle.size()) / 1024 / 1024 << "MB" << std::endl;
    std::cout << "Number of grids = " << handle.gridCount() << std::endl;

    float* h_outImage = reinterpret_cast<float*>(imageBuffer.data());

    float              wBBoxDimZ = (float)h_grid->worldBBox().dim()[2] * 2;
    Vec3T              wBBoxCenter = Vec3T(h_grid->worldBBox().min() + h_grid->worldBBox().dim() * 0.5f);
    nanovdb::CoordBBox treeIndexBbox = h_grid->tree().bbox();
    std::cout << "Bounds: "
              << "[" << treeIndexBbox.min()[0] << "," << treeIndexBbox.min()[1] << "," << treeIndexBbox.min()[2] << "] -> ["
              << treeIndexBbox.max()[0] << "," << treeIndexBbox.max()[1] << "," << treeIndexBbox.max()[2] << "]" << std::endl;

    RayGenOp<Vec3T> rayGenOp(wBBoxDimZ, wBBoxCenter);
    CompositeOp     compositeOp;

    auto renderOp = [width, height, rayGenOp, compositeOp, treeIndexBbox, wBBoxDimZ] __hostdev__(int start, int end, float* image, const GridT* grid) {
        // get an accessor.
        auto acc = grid->tree().getAccessor();

        for (int i = start; i < end; ++i) {
            Vec3T rayEye;
            Vec3T rayDir;
            rayGenOp(i, width, height, rayEye, rayDir);
            // generate ray.
            RayT wRay(rayEye, rayDir);
            // transform the ray to the grid's index-space.
            RayT iRay = wRay.worldToIndexF(*grid);
            // intersect...
            float  t0;
            CoordT ijk;
            float  v;
            if (nanovdb::ZeroCrossing(iRay, acc, ijk, v, t0)) {
                // write distance to surface. (we assume it is a uniform voxel)
                float wT0 = t0 * float(grid->voxelSize()[0]);
                compositeOp(image, i, width, height, wT0 / (wBBoxDimZ * 2), 1.0f);
            } else {
                // write background value.
                compositeOp(image, i, width, height, 0.0f, 0.0f);
            }
        }
    };

    if (numCPUIterations > 0)
    {
        float durationAvg = 0;
        for (int i = 0; i < numCPUIterations; ++i) {
            float duration = renderImage(false, renderOp, width, height, h_outImage, h_grid);
            //std::cout << "Duration(NanoVDB-Host) = " << duration << " ms" << std::endl;
            durationAvg += duration;
        }
        durationAvg /= numCPUIterations;
        std::cout << "Average Duration(NanoVDB-Host) = " << durationAvg << " ms" << std::endl;

		saveImage(std::string("raytrace_level_set") + imagePostFix + "-nanovdb-host.pfm", width, height, (float*)imageBuffer.data());
    }

#if defined(NANOVDB_USE_CUDA)
    if (numGPUIterations > 0)
    {
        handle.deviceUpload();

        auto* d_grid = handle.deviceGrid<float>();
        if (!d_grid)
            throw std::runtime_error("GridHandle does not contain a valid device grid");

        imageBuffer.deviceUpload();
        float* d_outImage = reinterpret_cast<float*>(imageBuffer.deviceData());

        {
            float durationAvg = 0;
            for (int i = 0; i < numGPUIterations; ++i) {
                float duration = renderImage(true, renderOp, width, height, d_outImage, d_grid);
                //std::cout << "Duration(NanoVDB-Cuda) = " << duration << " ms" << std::endl;
                durationAvg += duration;
            }
            durationAvg /= numGPUIterations;
            std::cout << "Average Duration(NanoVDB-Cuda) = " << durationAvg << " ms" << std::endl;

            imageBuffer.deviceDownload();
            saveImage(std::string("raytrace_level_set") + imagePostFix + "-nanovdb-cuda.pfm", width, height, (float*)imageBuffer.data());
        }
    }
#endif
}


void runNanoVTT(nanovtt::GridHandle<BufferT>& handle, int numCPUIterations, int numGPUIterations, bool useSubVoxelAccuracy, std::string& imagePostFix, int width, int height, BufferT& imageBuffer)
{
    using GridT = nanovtt::FloatGrid;
    using CoordT = nanovdb::Coord;
    using RealT = float;
    using Vec3T = nanovdb::Vec3<RealT>;
    using RayT = nanovdb::Ray<RealT>;

    auto* h_grid = handle.template grid<float>();
    if (!h_grid)
        throw std::runtime_error("GridHandle does not contain a valid host grid");

    std::cout << "NanoVTT buffer size = " << static_cast<float>(handle.size()) / 1024 / 1024 << "MB" << std::endl;
    std::cout << "Number of grids = " << handle.gridCount() << std::endl;

    float* h_outImage = reinterpret_cast<float*>(imageBuffer.data());

    float              wBBoxDimZ = (float)h_grid->worldBBox().dim()[2] * 2;
    Vec3T              wBBoxCenter = Vec3T(h_grid->worldBBox().min() + h_grid->worldBBox().dim() * 0.5f);
    nanovdb::CoordBBox treeIndexBbox = h_grid->tree().bbox();
    std::cout << "Bounds: "
        << "[" << treeIndexBbox.min()[0] << "," << treeIndexBbox.min()[1] << "," << treeIndexBbox.min()[2] << "] -> ["
        << treeIndexBbox.max()[0] << "," << treeIndexBbox.max()[1] << "," << treeIndexBbox.max()[2] << "]" << std::endl;

    RayGenOp<Vec3T> rayGenOp(wBBoxDimZ, wBBoxCenter);
    CompositeOp     compositeOp;

    // HDDA with sampler and sub-voxel accuracy
    auto renderOpSubVoxel = [width, height, rayGenOp, compositeOp, treeIndexBbox, wBBoxDimZ] __hostdev__(int start, int end, float* image, const GridT * grid) {
        auto acc = grid->tree().getAccessor();
        auto sampler = nanovtt::createSampler<1>(acc);
        constexpr float isoValue = 0.f;
        for (int i = start; i < end; ++i) {
            Vec3T rayEye;
            Vec3T rayDir;
            rayGenOp(i, width, height, rayEye, rayDir);
            // generate ray.
            RayT wRay(rayEye, rayDir);
            // transform the ray to the grid's index-space.
            RayT iRay = wRay.worldToIndexF(*grid);
            // intersect...
            float  t0 = 0.f;
            if (nanovtt::zeroCrossingSubVoxel(
                isoValue, iRay, sampler, t0)) {
                // write distance to surface. (we assume it is a uniform voxel)
                float wT0 = t0 * float(grid->voxelSize()[0]);
                compositeOp(image, i, width, height, wT0 / (wBBoxDimZ * 2), 1.0f);
            }
            else {
                // write background value.
                compositeOp(image, i, width, height, 0.0f, 0.0f);
            }
        }
    };

    // HDDA with sampler and sub-voxel accuracy for PNanoVTT
    auto renderOpSubVoxelPNanoVTT = [width, height, rayGenOp, compositeOp, treeIndexBbox, wBBoxDimZ] (int start, int end, float* image, const GridT * grid) {

        const pnanovdb_buf_t buf = pnanovdb_make_buf(const_cast<uint32_t*>(
            reinterpret_cast<uint32_t const*>(grid->tree().getContainer())),
            grid->tree().getContainer()->m_bufferSize);
        const pnanovdb_grid_handle_t pGridHandle{ pnanovdb_address_offset64(pnanovdb_address_null(), grid->tree().getContainer()->m_containerSize) };
        const pnanovtt_tree_handle_t pTreeHandle{ pnanovdb_address_offset64(pnanovdb_address_null(), grid->tree().getContainer()->m_containerSize + PNANOVDB_GRID_SIZE) };

        pnanovtt_cached_sampler_t psampler;
        pnanovtt_cached_sampler_init(buf, PNANOVDB_REF(psampler), pGridHandle, pTreeHandle);

        constexpr float isoValue = 0.f;
        for (int i = start; i < end; ++i) {
            Vec3T rayEye;
            Vec3T rayDir;
            rayGenOp(i, width, height, rayEye, rayDir);
            // generate ray.
            RayT wRay(rayEye, rayDir);
            // transform the ray to the grid's index-space.
            RayT iRay = wRay.worldToIndexF(*grid);
            // intersect...
            float  t0 = 0.f;
            float vHit = 0.f;
            pnanovdb_vec3_t pHit;
            const pnanovdb_vec3_t origin{ iRay.eye()[0], iRay.eye()[1], iRay.eye()[2] };
            const pnanovdb_vec3_t direction{ iRay.dir()[0], iRay.dir()[1], iRay.dir()[2] };
            const float tMin = iRay.t0();
            const float tMax = iRay.t1();

            const pnanovdb_bool_t foundZC = pnanovtt_hdda_zero_crossing_sub_voxel(
                buf,
                PNANOVDB_REF(psampler),
                PNANOVDB_REF(origin),
                tMin,
                PNANOVDB_REF(direction),
                tMax,
                isoValue,
                PNANOVDB_REF(t0),
                PNANOVDB_REF(vHit),
                PNANOVDB_REF(pHit),
                PNANOVTT_TRAITS_MAX_DEPTH);
            if (foundZC) {
                // write distance to surface. (we assume it is a uniform voxel)
                float wT0 = t0 * float(grid->voxelSize()[0]);
                compositeOp(image, i, width, height, wT0 / (wBBoxDimZ * 2), 1.0f);
            }
            else {
                // write background value.
                compositeOp(image, i, width, height, 0.0f, 0.0f);
            }
        }
    };

    // Original HDDA approach
    auto renderOp = [width, height, rayGenOp, compositeOp, treeIndexBbox, wBBoxDimZ] __hostdev__(int start, int end, float* image, const GridT * grid) {
        auto acc = grid->tree().getAccessor();

        for (int i = start; i < end; ++i) {
            Vec3T rayEye;
            Vec3T rayDir;
            rayGenOp(i, width, height, rayEye, rayDir);
            // generate ray.
            RayT wRay(rayEye, rayDir);
            // transform the ray to the grid's index-space.
            RayT iRay = wRay.worldToIndexF(*grid);
            // intersect...
            float  t0;
            CoordT ijk;
            float  v;
            // Note: A few changes were required to nanovdb::ZeroCrossing to make it work with adaptive volumes.
            //       Hence we use nanovtt::ZeroCrossing here.
            if (nanovtt::ZeroCrossing(iRay, acc, ijk, v, t0)) {
                // write distance to surface. (we assume it is a uniform voxel)
                float wT0 = t0 * float(grid->voxelSize()[0]);
                compositeOp(image, i, width, height, wT0 / (wBBoxDimZ * 2), 1.0f);
            }
            else {
                // write background value.
                compositeOp(image, i, width, height, 0.0f, 0.0f);
            }
        }
    };

    // Original HDDA approach for PNanoVTT
    auto renderOpPNanoVTT = [width, height, rayGenOp, compositeOp, treeIndexBbox, wBBoxDimZ] (int start, int end, float* image, const GridT * grid) {

        const pnanovdb_grid_type_t grid_type = PNANOVDB_GRID_TYPE_FLOAT;
        const pnanovdb_buf_t buf = pnanovdb_make_buf(const_cast<uint32_t*>(
            reinterpret_cast<uint32_t const*>(grid->tree().getContainer())),
            grid->tree().getContainer()->m_bufferSize);
        const pnanovdb_grid_handle_t pGridHandle{ pnanovdb_address_offset64(pnanovdb_address_null(), grid->tree().getContainer()->m_containerSize) };
        const pnanovtt_tree_handle_t pTreeHandle{ pnanovdb_address_offset64(pnanovdb_address_null(), grid->tree().getContainer()->m_containerSize + PNANOVDB_GRID_SIZE) };

        pnanovtt_readaccessor_t pacc;
        pnanovtt_readaccessor_init(buf, PNANOVDB_REF(pacc), pGridHandle, pTreeHandle);

        for (int i = start; i < end; ++i) {
            Vec3T rayEye;
            Vec3T rayDir;
            rayGenOp(i, width, height, rayEye, rayDir);
            // generate ray.
            RayT wRay(rayEye, rayDir);
            // transform the ray to the grid's index-space.
            RayT iRay = wRay.worldToIndexF(*grid);
            // intersect...
            float  t0 = 0.f;
            float vHit = 0.f;
            pnanovdb_coord_t ijk;
            const pnanovdb_vec3_t origin{ iRay.eye()[0], iRay.eye()[1], iRay.eye()[2] };
            const pnanovdb_vec3_t direction{ iRay.dir()[0], iRay.dir()[1], iRay.dir()[2] };
            const float tMin = iRay.t0();
            const float tMax = iRay.t1();

            const pnanovdb_bool_t foundZC = pnanovtt_hdda_zero_crossing(
                grid_type,
                buf,
                PNANOVDB_REF(pacc),
                PNANOVDB_REF(origin),
                tMin,
                PNANOVDB_REF(direction),
                tMax,
                PNANOVDB_REF(t0),
                PNANOVDB_REF(vHit),
                PNANOVDB_REF(ijk),
                PNANOVTT_TRAITS_MAX_DEPTH);
            if (foundZC) {
                // write distance to surface. (we assume it is a uniform voxel)
                float wT0 = t0 * float(grid->voxelSize()[0]);
                compositeOp(image, i, width, height, wT0 / (wBBoxDimZ * 2), 1.0f);
            }
            else {
                // write background value.
                compositeOp(image, i, width, height, 0.0f, 0.0f);
            }
        }
    };

    if (numCPUIterations > 0)
    {
        float durationAvg = 0;
        for (int i = 0; i < numCPUIterations; ++i) {
            float duration = 0;
            if (useSubVoxelAccuracy)
            {
                duration = renderImage(false, renderOpSubVoxel, width, height, h_outImage, h_grid);
            }
            else {
                duration = renderImage(false, renderOp, width, height, h_outImage, h_grid);
            }
            //std::cout << "Duration(NanoVTT-Host) = " << duration << " ms" << std::endl;
            durationAvg += duration;
        }
        durationAvg /= numCPUIterations;
        std::cout << "Average Duration(NanoVTT-Host) = " << durationAvg << " ms" << std::endl;
        saveImage(std::string("raytrace_level_set") + imagePostFix + "-nanovtt-host.pfm", width, height, (float*)imageBuffer.data());
    }

    if (numCPUIterations > 0)
    {
        float durationAvg = 0;
        for (int i = 0; i < numCPUIterations; ++i) {
            float duration = 0;
            if (useSubVoxelAccuracy)
            {
                duration = renderImageHost(renderOpSubVoxelPNanoVTT, width, height, h_outImage, h_grid);
            }
            else {
                duration = renderImageHost(renderOpPNanoVTT, width, height, h_outImage, h_grid);
            }
            //std::cout << "Duration(PNanoVTT-Host) = " << duration << " ms" << std::endl;
            durationAvg += duration;
        }
        durationAvg /= numCPUIterations;
        std::cout << "Average Duration(PNanoVTT-Host) = " << durationAvg << " ms" << std::endl;
        saveImage(std::string("raytrace_level_set") + imagePostFix + "-pnanovtt-host.pfm", width, height, (float*)imageBuffer.data());
    }

#if defined(NANOVDB_USE_CUDA)
    if (numGPUIterations > 0)
    {
        handle.deviceUpload();

        auto* d_grid = handle.template deviceGrid<float>();
        if (!d_grid)
            throw std::runtime_error("GridHandle does not contain a valid device grid");

        imageBuffer.deviceUpload();
        float* d_outImage = reinterpret_cast<float*>(imageBuffer.deviceData());

        {
            float durationAvg = 0;
            for (int i = 0; i < numGPUIterations; ++i) {
                float duration = 0;
                if (useSubVoxelAccuracy)
                {
                    duration = renderImage(true, renderOpSubVoxel, width, height, d_outImage, d_grid);
                }
                else {
                    duration = renderImage(true, renderOp, width, height, d_outImage, d_grid);
                }
                //std::cout << "Duration(NanoVTT-Cuda) = " << duration << " ms" << std::endl;
                durationAvg += duration;
            }
            durationAvg /= numGPUIterations;
            std::cout << "Average Duration(NanoVTT-Cuda) = " << durationAvg << " ms" << std::endl;
            imageBuffer.deviceDownload();
            saveImage(std::string("raytrace_level_set") + imagePostFix + "-nanovtt-cuda.pfm", width, height, (float*)imageBuffer.data());
        }
    }
#endif
}
