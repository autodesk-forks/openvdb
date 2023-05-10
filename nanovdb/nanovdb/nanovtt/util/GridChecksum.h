// Copyright Contributors to the OpenVDB Project
// SPDX-License-Identifier: MPL-2.0

/*!
    \file GridChecksum.h

    \authors Ken Museth and Autodesk
*/

#ifndef NANOVTT_GRIDCHECKSUM_H_HAS_BEEN_INCLUDED
#define NANOVTT_GRIDCHECKSUM_H_HAS_BEEN_INCLUDED

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
#pragma clang diagnostic ignored "-Wextra-semi-stmt"
#endif
#include <nanovdb/util/GridChecksum.h>
#if defined(__clang__)
#pragma clang diagnostic pop
#endif

namespace nanovtt {

/// @brief Return the (2 x CRC32) checksum of the specified @a grid
///
/// @param grid Grid from which the checksum is computed.
/// @param mode Defines the mode of computation for the checksum.
template <typename ValueT>
uint64_t checksum(const NanoVTTGrid<ValueT> &grid, ::nanovdb::ChecksumMode mode = ::nanovdb::ChecksumMode::Default);

/// @brief Return true if the checksum of the @a grid matches the expected
///        value already encoded into the grid's meta data.
///
/// @param grid Grid whose checksum is validated.
/// @param mode Defines the mode of computation for the checksum.
template <typename ValueT>
bool validateChecksum(const NanoVTTGrid<ValueT> &grid, ::nanovdb::ChecksumMode mode = ::nanovdb::ChecksumMode::Default);

/// @brief Updates the checksum of a grid
///
/// @param grid Grid whose checksum will be updated.
/// @param mode Defines the mode of computation for the checksum.
template <typename ValueT>
void updateChecksum(NanoVTTGrid<ValueT> &grid, ::nanovdb::ChecksumMode mode = ::nanovdb::ChecksumMode::Default);

template <typename ValueT>
nanovdb::GridChecksum gridChecksum(const NanoVTTGrid<ValueT> &grid, ::nanovdb::ChecksumMode mode)
{
    static const size_t offset = 16;

    if (mode == ::nanovdb::ChecksumMode::Disable) return nanovdb::GridChecksum();

    const auto &tree = grid.tree();
    nanovdb::CRC32 crc;
    uint8_t const* begin;
    uint8_t const* end;
    uint32_t mCRC[2] = { nanovdb::CRC32::EMPTY, nanovdb::CRC32::EMPTY };

    // process Container + (Grid + TreeBase), but exclude mMagic and mChecksum
    std::vector<std::uint_fast32_t> checksums;
    begin = reinterpret_cast<const uint8_t*>(tree.getContainer());
    end = begin + tree.getContainer()->memUsage();
    crc(begin, end);
    checksums.push_back(crc.checksum());
    crc.reset();
    begin = reinterpret_cast<const uint8_t*>(&grid);
    end = begin + grid.memUsage() + sizeof(nanovtt::TreeBase<ValueT>);
    crc(begin + offset, end);
    checksums.push_back(crc.checksum());
    crc.reset();
    crc(checksums.data(), sizeof(std::uint_fast32_t) * checksums.size());
    mCRC[0] = static_cast<uint32_t>(crc.checksum());
    crc.reset();

    if (mode == ::nanovdb::ChecksumMode::Partial || tree.isEmpty()) return nanovdb::GridChecksum(mCRC[0], mCRC[1]);

    // process Tree Data
    begin = reinterpret_cast<const uint8_t*>(&tree);
    end = begin + tree.memUsage();
    crc(begin + sizeof(nanovtt::TreeBase<ValueT>), end);
    mCRC[1] = static_cast<uint32_t>(crc.checksum());

    return nanovdb::GridChecksum(mCRC[0], mCRC[1]);
}

template <typename ValueT>
uint64_t checksum(const NanoVTTGrid<ValueT> &grid, ::nanovdb::ChecksumMode mode)
{
    nanovdb::GridChecksum cs = gridChecksum(grid, mode);
    return cs.checksum();
}

template <typename ValueT>
bool validateChecksum(const NanoVTTGrid<ValueT> &grid, ::nanovdb::ChecksumMode mode)
{
    nanovdb::GridChecksum cs1(grid.checksum(), mode), cs2;
    cs2 = gridChecksum(grid, cs1.mode());
    return cs1 == cs2;
}

template <typename ValueT>
void updateChecksum(NanoVTTGrid<ValueT> &grid, ::nanovdb::ChecksumMode mode)
{
    nanovdb::GridChecksum cs = gridChecksum(grid, mode);
    grid.data()->mChecksum = cs.checksum();
}

} // namespace nanovtt

#endif // NANOVTT_GRIDCHECKSUM_H_HAS_BEEN_INCLUDED
