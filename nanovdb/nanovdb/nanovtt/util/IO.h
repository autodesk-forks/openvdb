// Copyright Contributors to the OpenVDB Project
// SPDX-License-Identifier: MPL-2.0

/*!
    \file IO.h

    \authors Ken Museth and Autodesk
*/

#ifndef NANOVTT_IO_H_HAS_BEEN_INCLUDED
#define NANOVTT_IO_H_HAS_BEEN_INCLUDED

#include <fstream>
#include <iostream>
#include <nanovdb/nanovtt/util/GridHandle.h>
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-parameter"
#pragma clang diagnostic ignored "-Wextra-semi"
#pragma clang diagnostic ignored "-Wold-style-cast"
#pragma clang diagnostic ignored "-Wimplicit-float-conversion"
#pragma clang diagnostic ignored "-Wmissing-braces"
#pragma clang diagnostic ignored "-Wimplicit-int-conversion"
#pragma clang diagnostic ignored "-Wreserved-id-macro"
#pragma clang diagnostic ignored "-Wcast-qual"
#pragma clang diagnostic ignored "-Wdouble-promotion"
#pragma clang diagnostic ignored "-Wnewline-eof"
#pragma clang diagnostic ignored "-Wshadow"
#pragma clang diagnostic ignored "-Wshorten-64-to-32"
#pragma clang diagnostic ignored "-Wunused-template"
#endif
#include <nanovdb/util/IO.h>
#if defined(__clang__)
#pragma clang diagnostic pop
#endif

namespace nanovtt
{
namespace io
{

/// @brief Write a single grid-handle to file (over-writing existing content of the file)
template<typename BufferT>
void writeGrid(const std::string& in_fileName, const GridHandle<BufferT>& in_handle, nanovdb::io::Codec in_codec = nanovdb::io::Codec::NONE, const int in_verbose = 0)
{
    if (in_verbose > 0) {
        std::cout << "Saving NanoVTT as '" << in_fileName << "' begin" << std::endl;
		std::cout << "Codec = " << (in_codec == nanovdb::io::Codec::NONE ? "NONE" : "Not NONE") << std::endl;
    }
    const uint64_t bufferSize = *reinterpret_cast<const uint64_t*>(in_handle.data());
    std::ofstream os(in_fileName, std::ios::binary);
    os.write(reinterpret_cast<const char*>(in_handle.data()), bufferSize);
    if (in_verbose > 0) {
        std::cout << "Saving NanoVTT as '" << in_fileName << "' end" << std::endl;
    }
}

/// @brief Read a single grid-handle from file. The grid-handle can contain several grids.
template<typename BufferT = nanovdb::HostBuffer>
GridHandle<BufferT> readGrid(const std::string& in_fileName, const int in_verbose = 0, const BufferT& in_buffer = BufferT())
{
    nanovtt::GridHandle<BufferT> handle;
    if (in_verbose > 0) {
        std::cout << "Reading NanoVTT grid-handle from file '" << in_fileName << "' begin" << std::endl;
    }
    std::ifstream is(in_fileName, std::ios::binary);
    if (!is.is_open()) {
        if (in_verbose > 0)
            std::cerr << "failed to open nvtt file" << std::endl;
        return handle;
    }
    uint64_t bufferSize = 0;
    std::streamoff so = 0;
    // Chorus: here we cannot overflow buffer 'bufferSize' since we read sizeof(uint64_t)
    is.read(reinterpret_cast<char*>(&bufferSize), sizeof(decltype(bufferSize))); /* Flawfinder: ignore */
    // Chorus: here we make sure that bufferSize is not larger than the file itself
    const uint64_t fileSize = [&]() {
        std::ifstream is2(in_fileName, std::ifstream::ate | std::ifstream::binary);
        if (!is2.is_open() || is2.tellg() < 0) {
            return static_cast<uint64_t>(0);
        }
        else {
            return static_cast<uint64_t>(is2.tellg());
        }
    }();
    if (bufferSize > fileSize) {
        // Not possible, returning
        return handle;
    }
    is.seekg(so, is.beg);
    // Chorus: here the buffer is created using bufferSize and below we read at most bufferSize.
    handle = BufferT::create(bufferSize, &in_buffer);
    if (handle.size() >= bufferSize) {
        is.read(reinterpret_cast<char*>(handle.data()), bufferSize); /* Flawfinder: ignore */
    }
    return handle;
}

using GridMetaData = nanovdb::GridData;

/// @brief Read the meta data for all grids stored in the file.
inline std::vector<GridMetaData> readGridMetaData(const std::string& in_fileName, const int in_verbose = 0)
{
    std::vector<GridMetaData> metas;
    std::ifstream is(in_fileName, std::ios::in | std::ios::binary);
    if (!is.is_open()) {
        if (in_verbose > 0)
            std::cerr << "Failed to open nvtt file \"" << in_fileName << "\"" << std::endl;
        return metas;
    }
    uint64_t bufferSize = 0;
    uint64_t containerSize = 0;
    // Chorus: here the buffers cannot overflow since we read sizeof(uint64_t)
    is.read(reinterpret_cast<char*>(&bufferSize), sizeof(decltype(bufferSize))); /* Flawfinder: ignore */
    is.read(reinterpret_cast<char*>(&containerSize), sizeof(decltype(containerSize))); /* Flawfinder: ignore */
    // Chorus: here we make sure that bufferSize and containerSize are not larger than the file itself
    const uint64_t fileSize = [&]() {
        std::ifstream is2(in_fileName, std::ifstream::ate | std::ifstream::binary);
        if (!is2.is_open() || is2.tellg() < 0) {
            return static_cast<uint64_t>(0);
        }
        else {
            return static_cast<uint64_t>(is2.tellg());
        }
    }();
    if (bufferSize > fileSize || containerSize > fileSize) {
        // Not possible, returning
        return metas;
    }
    GridMetaData meta;
    std::streamoff so = static_cast<std::streamoff>(containerSize);
    if (sizeof(decltype(meta)) + static_cast<uint64_t>(so) > bufferSize) {
        return metas;
    }
    is.seekg(so, is.beg);
    // Chorus: here we cannot overflow buffer 'meta' since we read sizeof(meta)
    is.read(reinterpret_cast<char*>(&meta), sizeof(decltype(meta))); /* Flawfinder: ignore */
    metas.push_back(meta);
    so += meta.mGridSize;
    uint32_t gridCount = meta.mGridCount;
    for (uint32_t gridId = 1; gridId < gridCount; ++gridId) {
        if (sizeof(decltype(meta)) + static_cast<uint64_t>(so) > bufferSize)
            return metas;
        is.seekg(so, is.beg);
        // Chorus: here we cannot overflow buffer 'meta' since we read sizeof(meta)
        is.read(reinterpret_cast<char*>(&meta), sizeof(decltype(meta))); /* Flawfinder: ignore */
        metas.push_back(meta);
        so += meta.mGridSize;
    }
    return metas;
}

} // namespace io
} // namespace nanovtt

#endif // NANOVTT_IO_H_HAS_BEEN_INCLUDED
