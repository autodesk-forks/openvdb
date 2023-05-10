// Copyright Contributors to the OpenVDB Project
// SPDX-License-Identifier: MPL-2.0

#include <algorithm>
#include <iostream>
#include <fstream>
#include <nanovdb/util/IO.h>
#include <nanovdb/util/Primitives.h>
#include <nanovdb/util/CudaDeviceBuffer.h>
#include <nanovdb/nanovtt/NanoVTT.h>
#include <nanovdb/nanovtt/util/GridHandle.h>
#include <nanovdb/nanovtt/util/Primitives.h>
#include <nanovdb/nanovtt/util/IO.h>

#if defined(NANOVDB_USE_CUDA)
using BufferT = nanovdb::CudaDeviceBuffer;
#else
using BufferT = nanovdb::HostBuffer;
#endif

extern void runNanoVDB(nanovdb::GridHandle<BufferT>& handle, int numCPUIterations, int numGPUIterations, std::string& imagePostFix, int width, int height, BufferT& imageBuffer);
#if defined(NANOVDB_USE_OPENVDB)
extern void runOpenVDB(nanovdb::GridHandle<BufferT>& handle, int numIterations, int width, int height, BufferT& imageBuffer);
#endif
extern void runNanoVTT(nanovtt::GridHandle<BufferT>& handle, int numCPUIterations, int numGPUIterations, bool useSubVoxelAccuracy, std::string& imagePostFix, int width, int height, BufferT& imageBuffer);


void usage [[noreturn]] (const std::string& progName, int exitStatus = EXIT_FAILURE)
{
    std::cerr << "\nUsage: " << progName << " [options] vdb [vtt]\n"
        << "Options:\n"
        << "-h,--help\tPrints this message\n"
        << "-subvoxel_accuracy\t Intersects the level set using sub-voxel accuracy\n"
        << "-image_post_fix string\t Adds this postfix to the saved image file name\n"
        << "-num_cpu_iter int\t Sets the number of CPU iterations\n"
        << "-num_gpu_iter int\t Sets the number of GPU iterations\n"
        << "-num_vdb_iter int\t Sets the number of VDB CPU iterations\n";
    exit(exitStatus);
}

int main(int argc, char** argv)
{

    int numCPUIterations = 50;
    int numGPUIterations = 50;
    int numVDBIterations = 1;
    const int width = 1024;
    const int height = 1024;
    bool useSubVoxelAccuracy = false;
    std::string imagePostFix = "";

    std::vector<std::string> fileNames;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg[0] == '-') {
            if (arg == "-h" || arg == "--help") {
                usage(argv[0], EXIT_SUCCESS);
            }
            else if (arg == "-subvoxel_accuracy") {
                useSubVoxelAccuracy = true;
            }
            else if (arg == "-image_post_fix") {
                imagePostFix = argv[++i];
            }
            else if (arg == "-num_cpu_iter") {
                if (i + 1 < argc)
                    numCPUIterations = atoi(argv[++i]);
            }
            else if (arg == "-num_gpu_iter") {
                if (i + 1 < argc)
                    numGPUIterations = atoi(argv[++i]);
            }
            else if (arg == "-num_vdb_iter") {
                numVDBIterations = atoi(argv[++i]);
            }
            else {
                std::cerr << "\nIllegal option: \"" << arg << "\"\n";
                usage(argv[0]);
            }
        }
        else if (!arg.empty()) {
            fileNames.push_back(arg);
        }
    }
    if (!useSubVoxelAccuracy) {
        try {
            nanovdb::GridHandle<BufferT> handle;
            if (!fileNames.empty()) {
                handle = nanovdb::io::readGrid<BufferT>(fileNames[0]);
                std::cout << "Loaded NanoVDB grid[" << handle.gridMetaData()->shortGridName() << "] from file '" << fileNames[0] << "'...\n";
            }
            else {
                handle = nanovdb::createLevelSetSphere<float, float, BufferT>(100.0f, nanovdb::Vec3f(-20, 0, 0), 1.0, 3.0, nanovdb::Vec3d(0), "sphere");
            }

            if (handle.gridMetaData()->isLevelSet() == false) {
                throw std::runtime_error("Grid must be a level set");
            }

            {
                BufferT   imageBuffer;
                imageBuffer.init(width * height * sizeof(float));
                runNanoVDB(handle, numCPUIterations, numGPUIterations, imagePostFix, width, height, imageBuffer);
            }

#if defined(NANOVDB_USE_OPENVDB)
            if (numVDBIterations > 0)
            {
                BufferT   imageBuffer;
                imageBuffer.init(width * height * sizeof(float));
                runOpenVDB(handle, numVDBIterations, width, height, imageBuffer);
            }
#endif
        }
        catch (const std::exception& e) {
            std::cerr << "An exception occurred: \"" << e.what() << "\"" << std::endl;
        }
    }

    // NANOVTT
    if (fileNames.size() > 1) {
        try {
            nanovtt::GridHandle<BufferT> handle;
            handle = nanovtt::io::readGrid<BufferT>(fileNames[1]);
            std::cout << "Loaded NanoVTT grid[" << handle.gridMetaData()->shortGridName() << "] from file '" << fileNames[1] << "'...\n";

            if (handle.gridMetaData()->isLevelSet() == false) {
                throw std::runtime_error("Grid must be a level set");
            }

            BufferT   imageBuffer;
            imageBuffer.init(width * height * sizeof(float));

            runNanoVTT(handle, numCPUIterations, numGPUIterations, useSubVoxelAccuracy, imagePostFix, width, height, imageBuffer);
        }
        catch (const std::exception& e) {
            std::cerr << "An exception occurred: \"" << e.what() << "\"" << std::endl;
        }
    }

    return 0;
}
