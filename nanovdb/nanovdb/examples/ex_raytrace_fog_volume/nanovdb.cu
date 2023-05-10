// Copyright Contributors to the OpenVDB Project
// SPDX-License-Identifier: MPL-2.0

#include <cmath>
#include <chrono>

#include <cuda_profiler_api.h>

#include <nanovdb/util/IO.h>
#include <nanovdb/util/CudaDeviceBuffer.h>
#include <nanovdb/util/Ray.h>
#include <nanovdb/util/HDDA.h>
#include <nanovdb/util/SampleFromVoxels.h>

#include "common.h"

#include <nanovdb/nanovtt/NanoVTT.h>
#include <nanovdb/nanovtt/util/GridHandle.h>
#include <nanovdb/nanovtt/util/SampleFromVoxels.h>

#if defined(NANOVDB_USE_CUDA)
using BufferT = nanovdb::CudaDeviceBuffer;
#else
using BufferT = nanovdb::HostBuffer;
#endif

namespace nanovdb {

template<typename RayT, typename AccT>
__hostdev__ inline bool firstActive(RayT& ray, AccT& acc, float& t)
{
    using namespace nanovdb;
    static const float Delta = 1.0001f;        // forward step-size along the ray to avoid getting stuck
    t = ray.t0();                              // initiate time
    Coord ijk = RoundDown<Coord>(ray.start()); // first voxel inside bbox
    for (HDDA<RayT, Coord> hdda(ray, acc.getDim(ijk, ray)); !acc.isActive(ijk); hdda.update(ray, acc.getDim(ijk, ray)))
    {
        if (!hdda.step())
            return false;                        // leap-frog HDDA and exit if ray bound is exceeded
        t = hdda.time() + Delta;                // update time
        ijk = RoundDown<Coord>(ray(t));       // update ijk
    }
    return true;
}

}

template<typename VoxelT>
void runNanoVDBInternal(nanovdb::GridHandle<BufferT>& handle, int numCPUIterations, int numGPUIterations, bool useHDDAAndSampler, std::string& imagePostFix, int width, int height, BufferT& imageBuffer)
{
    using GridT = nanovdb::NanoGrid<VoxelT>;
    using CoordT = nanovdb::Coord;
    using RealT = float;
    using Vec3T = nanovdb::Vec3<RealT>;
    using RayT = nanovdb::Ray<RealT>;

    auto* h_grid = handle.grid<VoxelT>();
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

    auto renderOpHDDA = [width, height, rayGenOp, compositeOp, treeIndexBbox] __hostdev__(int start, int end, float* image, const GridT * grid) {
        // get an accessor.
        auto acc = grid->tree().getAccessor();
        auto sampler = nanovdb::createSampler<1, decltype(acc), /*cache*/true>(acc);

        for (int i = start; i < end; ++i) {
            Vec3T rayEye;
            Vec3T rayDir;
            rayGenOp(i, width, height, rayEye, rayDir);
            // generate ray.
            RayT wRay(rayEye, rayDir);
            // transform the ray to the grid's index-space.
            RayT iRay = wRay.worldToIndexF(*grid);
            // clip to bounds.
            if (iRay.clip(treeIndexBbox) == false) {
                compositeOp(image, i, width, height, 0.0f, 0.0f);
#ifdef __CUDA_ARCH__
                return;
#else
                continue;
#endif
            }
            float ti0 = iRay.t0(); // index-space hit params of the grid bbox
            float ti1 = iRay.t1();
            float ti01 = fabs(ti1 - ti0);
            {
                RayT fRay(iRay(ti0), iRay.dir(), 0.0, ti01);
                float tf;
                if (firstActive(fRay, acc, tf))
                    ti0 = ti0 + tf; // update ti0 to the first voxel hit
                else {
                    compositeOp(image, i, width, height, 0.0f, 0.0f);
#ifdef __CUDA_ARCH__
                    return;
#else
                    continue;
#endif
                }
            }
            {
                RayT bRay(iRay(ti1), -iRay.dir(), 0.0, ti01);
                float tb;
                if (firstActive(bRay, acc, tb))
                    ti1 = ti1 - tb; // update ti1 to the last voxel hit
                else {
                    compositeOp(image, i, width, height, 0.0f, 0.0f);
#ifdef __CUDA_ARCH__
                    return;
#else
                    continue;
#endif
                }
            }
            if (ti1 <= ti0) {
                compositeOp(image, i, width, height, 0.0f, 0.0f);
#ifdef __CUDA_ARCH__
                return;
#else
                continue;
#endif
            }
            ti0 -= 1.f; // expand by 1 voxel to ensure DDA march covers entirety of touched voxels
            ti1 += 1.f;

            // integrate...
            const float dt = 0.5f;
            float       transmittance = 1.0f;
            for (float t = ti0; t < ti1; t += dt) {
                float sigma = sampler(iRay(t)) * 0.1f;
                transmittance *= 1.0f - sigma * dt;
            }
            // write transmittance.
            compositeOp(image, i, width, height, 0.0f, 1.0f - transmittance);
        }
    };


    auto renderOp = [width, height, rayGenOp, compositeOp, treeIndexBbox] __hostdev__(int start, int end, float* image, const GridT* grid) {
        auto acc = grid->tree().getAccessor();

        for (int i = start; i < end; ++i) {
            Vec3T rayEye;
            Vec3T rayDir;
            rayGenOp(i, width, height, rayEye, rayDir);
            // generate ray.
            RayT wRay(rayEye, rayDir);
            // transform the ray to the grid's index-space.
            RayT iRay = wRay.worldToIndexF(*grid);
            // clip to bounds.
            if (iRay.clip(treeIndexBbox) == false) {
                compositeOp(image, i, width, height, 0.0f, 0.0f);
#ifdef __CUDA_ARCH__
                return;
#else
                continue;
#endif
            }
            // integrate...
            const float dt = 0.5f;
            float       transmittance = 1.0f;
            for (float t = iRay.t0(); t < iRay.t1(); t += dt) {
                float sigma = acc.getValue(CoordT::Floor(iRay(t))) * 0.1f;
                transmittance *= 1.0f - sigma * dt;
            }
            // write transmittance.
            compositeOp(image, i, width, height, 0.0f, 1.0f - transmittance);
        }
    };

    if (numCPUIterations > 0)
    {
        float durationAvg = 0;
        for (int i = 0; i < numCPUIterations; ++i) 
        {
            float duration = 0;
            if (useHDDAAndSampler) {
                duration = renderImage(false, renderOpHDDA, width, height, h_outImage, h_grid);
            }
            else {
                duration = renderImage(false, renderOp, width, height, h_outImage, h_grid);
            }
            //std::cout << "Duration(NanoVDB-Host) = " << duration << " ms" << std::endl;
            durationAvg += duration;
        }
        durationAvg /= numCPUIterations;
        std::cout << "Average Duration(NanoVDB-Host) = " << durationAvg << " ms" << std::endl;
		std::string fileName = std::string("raytrace_fog_volume") + imagePostFix + "-nanovdb-host.pfm";
		saveImage(fileName, width, height, (float*)imageBuffer.data());
    }

#if defined(NANOVDB_USE_CUDA)
    if (numGPUIterations > 0) {
        handle.deviceUpload();

        auto* d_grid = handle.deviceGrid<VoxelT>();
        if (!d_grid)
            throw std::runtime_error("GridHandle does not contain a valid device grid");

        imageBuffer.deviceUpload();
        float* d_outImage = reinterpret_cast<float*>(imageBuffer.deviceData());

        {
            float durationAvg = 0;
            for (int i = 0; i < numGPUIterations; ++i)
            {
                float duration = 0;
                if (useHDDAAndSampler) {
                    duration = renderImage(true, renderOpHDDA, width, height, d_outImage, d_grid);
                }
                else {
                    duration = renderImage(true, renderOp, width, height, d_outImage, d_grid);
                }
                //std::cout << i << ": Duration(NanoVDB-Cuda) = " << duration << " ms" << std::endl;
                durationAvg += duration;
            }
            durationAvg /= numGPUIterations;
            std::cout << "Average Duration(NanoVDB-Cuda) = " << durationAvg << " ms" << std::endl;
            imageBuffer.deviceDownload();
			std::string fileName = std::string("raytrace_fog_volume") + imagePostFix + "-nanovdb-cuda.pfm";
			saveImage(fileName, width, height, (float*)imageBuffer.data());
        }
    }
#endif
}

void runNanoVDB(nanovdb::GridHandle<BufferT>& handle, int numCPUIterations, int numGPUIterations, bool useHDDAAndSampler, std::string& imagePostFix, int width, int height, BufferT& imageBuffer)
{
    if (handle.grid<float>()) {
        runNanoVDBInternal<float>(handle, numCPUIterations, numGPUIterations, useHDDAAndSampler, imagePostFix, width, height, imageBuffer);
    }
    else if (handle.grid<nanovdb::Fp4>()) {
        runNanoVDBInternal<nanovdb::Fp4>(handle, numCPUIterations, numGPUIterations, useHDDAAndSampler, imagePostFix, width, height, imageBuffer);
    }
    else if (handle.grid<nanovdb::Fp8>()) {
        runNanoVDBInternal<nanovdb::Fp8>(handle, numCPUIterations, numGPUIterations, useHDDAAndSampler, imagePostFix, width, height, imageBuffer);
    }
    else if (handle.grid<nanovdb::FpN>()) {
        runNanoVDBInternal<nanovdb::FpN>(handle, numCPUIterations, numGPUIterations, useHDDAAndSampler, imagePostFix, width, height, imageBuffer);
    }
    else
        throw std::runtime_error("GridHandle does not contain a valid host grid");
}


namespace nanovtt {

constexpr bool useCachedSampler = false;

template<typename RayT, typename AccT>
__hostdev__ inline bool firstActive(RayT& ray, AccT& acc, float& t)
{
    using namespace nanovdb;
    static const float Delta = 1.0001f;        // forward step-size along the ray to avoid getting stuck
    t = ray.t0();                              // initiate time
    Coord ijk = RoundDown<Coord>(ray.start()); // first voxel inside bbox
    for (HDDA<RayT, Coord> hdda(ray, acc.getDim(ijk, ray)); !acc.isActive(ijk); hdda.update(ray, acc.getDim()))
    {
        if (!hdda.step())
            return false;                        // leap-frog HDDA and exit if ray bound is exceeded
        t = hdda.time() + Delta;                // update time
        ijk = RoundDown<Coord>(ray(t));       // update ijk
    }
    return true;
}

}

template<typename VoxelT>
void runNanoVTTInternal(nanovtt::GridHandle<BufferT>& handle, int numCPUIterations, int numGPUIterations, bool useHDDAAndSampler, std::string& imagePostFix, int width, int height, BufferT& imageBuffer)
{
    using GridT = nanovtt::Grid<nanovtt::Tree<VoxelT>>;
    using CoordT = nanovdb::Coord;
    using RealT = float;
    using Vec3T = nanovdb::Vec3<RealT>;
    using RayT = nanovdb::Ray<RealT>;

    auto* h_grid = handle.grid<VoxelT>();
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

    auto renderOpHDDA = [width, height, rayGenOp, compositeOp, treeIndexBbox] __hostdev__(int start, int end, float* image, const GridT * grid) {
        // get an accessor.
        auto acc = grid->tree().getAccessor();
        auto sampler = nanovtt::createSampler<1, decltype(acc), /*cache*/nanovtt::useCachedSampler>(acc);

        for (int i = start; i < end; ++i) {
            Vec3T rayEye;
            Vec3T rayDir;
            rayGenOp(i, width, height, rayEye, rayDir);
            // generate ray.
            RayT wRay(rayEye, rayDir);
            // transform the ray to the grid's index-space.
            RayT iRay = wRay.worldToIndexF(*grid);
            // clip to bounds.
            if (iRay.clip(treeIndexBbox) == false) {
                compositeOp(image, i, width, height, 0.0f, 0.0f);
#ifdef __CUDA_ARCH__
                return;
#else
                continue;
#endif
            }
            float ti0 = iRay.t0(); // index-space hit params of the grid bbox
            float ti1 = iRay.t1();
            float ti01 = fabs(ti1 - ti0);
            {
                RayT fRay(iRay(ti0), iRay.dir(), 0.0, ti01);
                float tf;
                if (nanovtt::firstActive(fRay, acc, tf))
                    ti0 = ti0 + tf; // update ti0 to the first voxel hit
                else {
                    compositeOp(image, i, width, height, 0.0f, 0.0f);
#ifdef __CUDA_ARCH__
                    return;
#else
                    continue;
#endif
                }
            }
            {
                RayT bRay(iRay(ti1), -iRay.dir(), 0.0, ti01);
                float tb;
                if (nanovtt::firstActive(bRay, acc, tb))
                    ti1 = ti1 - tb; // update ti1 to the last voxel hit
                else {
                    compositeOp(image, i, width, height, 0.0f, 0.0f);
#ifdef __CUDA_ARCH__
                    return;
#else
                    continue;
#endif
                }
            }
            if (ti1 <= ti0) {
                compositeOp(image, i, width, height, 0.0f, 0.0f);
#ifdef __CUDA_ARCH__
                return;
#else
                continue;
#endif
            }
            ti0 -= 1.f; // expand by 1 voxel to ensure DDA march covers entirety of touched voxels
            ti1 += 1.f;

            // integrate...
            const float dt = 0.5f;
            float       transmittance = 1.0f;
            for (float t = ti0; t < ti1; t += dt) {
                float sigma = sampler(iRay(t)) * 0.1f;
                transmittance *= 1.0f - sigma * dt;
            }
            // write transmittance.
            compositeOp(image, i, width, height, 0.0f, 1.0f - transmittance);
        }
    };


    auto renderOp = [width, height, rayGenOp, compositeOp, treeIndexBbox] __hostdev__(int start, int end, float* image, const GridT * grid) {
        auto acc = grid->tree().getAccessor();
        for (int i = start; i < end; ++i) {
            Vec3T rayEye;
            Vec3T rayDir;
            rayGenOp(i, width, height, rayEye, rayDir);
            // generate ray.
            RayT wRay(rayEye, rayDir);
            // transform the ray to the grid's index-space.
            RayT iRay = wRay.worldToIndexF(*grid);
            // clip to bounds.
            if (iRay.clip(treeIndexBbox) == false) {
                compositeOp(image, i, width, height, 0.0f, 0.0f);
#ifdef __CUDA_ARCH__
                return;
#else
                continue;
#endif
            }
            // integrate...
            const float dt = 0.5f;
            float       transmittance = 1.0f;          
            for (float t = iRay.t0(); t < iRay.t1(); t += dt) {
                float sigma = acc.getValue(CoordT::Floor(iRay(t))) * 0.1f;
                transmittance *= 1.0f - sigma * dt;
        }
            // write transmittance.
            compositeOp(image, i, width, height, 0.0f, 1.0f - transmittance);
        }
    };

    if (numCPUIterations > 0)
    {
        float durationAvg = 0;
        for (int i = 0; i < numCPUIterations; ++i) 
        {
            float duration = 0;
            if (useHDDAAndSampler) {
                duration = renderImage(false, renderOpHDDA, width, height, h_outImage, h_grid);
            }
            else {
                duration = renderImage(false, renderOp, width, height, h_outImage, h_grid);
            }
            //std::cout << "Duration(NanoVTT-Host) = " << duration << " ms" << std::endl;
            durationAvg += duration;
        }
        durationAvg /= numCPUIterations;
        std::cout << "Average Duration(NanoVTT-Host) = " << durationAvg << " ms" << std::endl;
		std::string fileName = std::string("raytrace_fog_volume") + imagePostFix + "-nanovtt-host.pfm";
		saveImage(fileName, width, height, (float*)imageBuffer.data());
    }

#if defined(NANOVDB_USE_CUDA)
    if (numGPUIterations > 0) {
        handle.deviceUpload();

        auto* d_grid = handle.deviceGrid<VoxelT>();
        if (!d_grid)
            throw std::runtime_error("GridHandle does not contain a valid device grid");

        imageBuffer.deviceUpload();
        float* d_outImage = reinterpret_cast<float*>(imageBuffer.deviceData());

		float durationAvg = 0;
		for (int i = 0; i < numGPUIterations; ++i)
		{
			float duration = 0;
			if (useHDDAAndSampler) {
				duration = renderImage(true, renderOpHDDA, width, height, d_outImage, d_grid);
			}
			else {
				duration = renderImage(true, renderOp, width, height, d_outImage, d_grid);
			}
			//std::cout << i << ": Duration(NanoVTT-Cuda) = " << duration << " ms" << std::endl;
			durationAvg += duration;
		}
		durationAvg /= numGPUIterations;
		std::cout << "Average Duration(NanoVTT-Cuda) = " << durationAvg << " ms" << std::endl;
		imageBuffer.deviceDownload();
		std::string fileName = std::string("raytrace_fog_volume") + imagePostFix + "-nanovtt-cuda.pfm";
		saveImage(fileName, width, height, (float*)imageBuffer.data());
    }
#endif
}

void runNanoVTT(nanovtt::GridHandle<BufferT>& handle, int numCPUIterations, int numGPUIterations, bool useHDDAAndSampler, std::string& imagePostFix, int width, int height, BufferT& imageBuffer)
{
    if (handle.grid<float>()) {
        runNanoVTTInternal<float>(handle, numCPUIterations, numGPUIterations, useHDDAAndSampler, imagePostFix, width, height, imageBuffer);
    }
    else if (handle.grid<nanovdb::Fp4>()) {
        runNanoVTTInternal<nanovdb::Fp4>(handle, numCPUIterations, numGPUIterations, useHDDAAndSampler, imagePostFix, width, height, imageBuffer);
    }
    else if (handle.grid<nanovdb::Fp8>()) {
        runNanoVTTInternal<nanovdb::Fp8>(handle, numCPUIterations, numGPUIterations, useHDDAAndSampler, imagePostFix, width, height, imageBuffer);
    }
    else if (handle.grid<nanovdb::FpN>()) {
        runNanoVTTInternal<nanovdb::FpN>(handle, numCPUIterations, numGPUIterations, useHDDAAndSampler, imagePostFix, width, height, imageBuffer);
    }
    else
        throw std::runtime_error("GridHandle does not contain a valid host grid");
}
