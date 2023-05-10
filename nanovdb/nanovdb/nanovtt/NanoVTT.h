// Copyright Contributors to the OpenVDB Project
// SPDX-License-Identifier: MPL-2.0

/*!
    \file   NanoVTT.h

    \authors Ken Museth and Autodesk

    \brief NanoVTT.h contains modified code from NanoVDB.h and implements a light-weight
           self-contained spatially adaptive Voxel Tile Tree (VTT) data structure named
           NanoVTT.

    \details 
           # Overview

           NanoVTT supports the same high-level API as NanoVDB enabling
           third-party applications to use the same code to manipulate NanoVTT and
           NanoVDB. NanoVTT is spatially adaptive and supports level-of-detail by
           storing data at progressively coarser levels. For fog and level set volumes,
           the adaptivity can be used to represent data at coarser resolution if it is
           sufficiently smooth. Thereby NanoVTT can preserve visual quality and in some
           cases incur a lower memory footprint. A single NanoVTT instance stores a
           number of typed Grids. All grids share the same grid topology in the form of
           a NanoVTT Container.
           
           By default, data is stored in tiles of size 4^3 or 5^3 at each
           level of the data structure, and resolution jumps by a factor of two between
           each level. A tile of size 4^3 is named a non-extended tile, and a tile of size
           5^3 is named an extended tile. Extended tiles take up more storage, but offer
           performance advantages when sampling as only a single lookup per sample point
           is required in the case of linear interpolation. NanoVTT supports storing
           a mix of extended and non-extended tiles thereby offering a tradeoff between
           storage requirements and sampling performance.

           When NANOVTT_SUPPORT_COALESCED_TILES is defined, NanoVTT also supports tiles
           coalesced into larger tiles. For example, coalescing 8 child tiles of size 4^3
           into a single larger tile of size 8^3 etc. This leads to a smaller memory footprint.
           However, it is slower in practice if the problem is not sufficiently memory bounded.
           On the other hand it can lead to speedups for memory bounded problems.
                      
           Spatial linear interpolators ensure a continuous interpolant
           across resolution jumps. The interpolators assume that the data is sampled at
           voxel corners to guarantee this. However, NanoVTT itself does not make any
           assumptions about where data is sampled.
           

           # Limitations

           \li Packing the data as a texture is currently only supported by PNanoVTT.
           \li NANOVTT_SUPPORT_COALESCED_TILES cannot currently be used in combination with quantization.
*/

#pragma once


// Supports coalescing tile data into as large 3D arrays as possible: 8^3, 16^3 etc..
//#define NANOVTT_SUPPORT_COALESCED_TILES


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
#endif
#if defined(_MSC_VER)
#include <intrin.h>
#endif
#include <nanovdb/NanoVDB.h>
#if defined(__clang__)
#pragma clang diagnostic pop
#endif

#define NANOVTT_MAGIC_NUMBER 0x4e616e6f56545430UL // "NanoVTT0" in hex - little endian (uint64_t)

#define NANOVTT_MAJOR_VERSION_NUMBER 1 // reflects changes to the ABI and hence also the file format
#define NANOVTT_MINOR_VERSION_NUMBER 0 //  reflects changes to the API but not ABI
#define NANOVTT_PATCH_VERSION_NUMBER 0 //  reflects changes that does not affect the ABI or API

namespace nanovtt
{
    template<typename T>
    __hostdev__  inline T Max(T a, T b) { return a > b ? a : b; }

    template<typename T>
    __hostdev__  inline T Min(T a, T b) { return a < b ? a : b; }


    // The data type used to index tiles.
    // uint32_t allows for 256GB of pure data (ie not including topology) in a single grid.
    using TileIndex = uint32_t;

#define kTILESIDELUT(in_name)                                            \
    static const int in_name[] = {                                       \
        1073741824, 536870912, 268435456, 134217728, 67108864, 33554432, \
        16777216,   8388608,   4194304,   2097152,   1048576,  524288,   \
        262144,     131072,    65536,     32768,     16384,    8192,     \
        4096,       2048,      1024,      512,       256,      128,      \
        64,         32,        16,        8,         4,        2};
#ifdef __CUDACC__
    __constant__ kTILESIDELUT(kTileSideLUTDevice)
#endif
    kTILESIDELUT(kTileSideLUT)

#define kINVTILESIDELUT(in_name)               \
    float in_name[] = {9.313225746154785e-10f, \
                       1.862645149230957e-09f, \
                       3.725290298461914e-09f, \
                       7.450580596923828e-09f, \
                       1.490116119384766e-08f, \
                       2.980232238769531e-08f, \
                       5.960464477539063e-08f, \
                       1.192092895507813e-07f, \
                       2.384185791015625e-07f, \
                       4.76837158203125e-07f,  \
                       9.5367431640625e-07f,   \
                       1.9073486328125e-06f,   \
                       3.814697265625e-06f,    \
                       7.62939453125e-06f,     \
                       1.52587890625e-05f,     \
                       3.0517578125e-05f,      \
                       6.103515625e-05f,       \
                       0.0001220703125f,       \
                       0.000244140625f,        \
                       0.00048828125f,         \
                       0.0009765625f,          \
                       0.001953125f,           \
                       0.00390625f,            \
                       0.0078125f,             \
                       0.015625f,              \
                       0.03125f,               \
                       0.0625f,                \
                       0.125f,                 \
                       0.25f,                  \
                       0.5f};
#ifdef __CUDACC__
        __constant__ static const kINVTILESIDELUT(kInvTileSideLUTDevice)
#endif
        static constexpr kINVTILESIDELUT(kInvTileSideLUT)

    // Traits<2, 4>
    struct Traits
    {
        static constexpr int kTileWidth = 2;
        static constexpr int kTileWidth2 = kTileWidth * kTileWidth;
        static constexpr int kTileSize = 8;
        static constexpr int kTileDataWidth = 4;
        static constexpr int kTileDataWidth2 = kTileDataWidth * kTileDataWidth;
        static constexpr int kTileDataSize = kTileDataWidth2 * kTileDataWidth;
        static constexpr int kExtendedTileDataWidth = 5;
        static constexpr int kExtendedTileDataWidth2 = kExtendedTileDataWidth * kExtendedTileDataWidth;
        static constexpr int kExtendedTileDataSize = kExtendedTileDataWidth2 * kExtendedTileDataWidth;
        static constexpr int kMaxDepth = 28;
        static constexpr bool kBinaryCompatible = false;
        static constexpr int kTileDomainWidth = 1073741824;
        static constexpr int kTileOffset = -536870912;
        static constexpr int kMinTileCoord = kTileOffset;
        static constexpr int kMaxTileCoord = kMinTileCoord + kTileDomainWidth - 1;
        static constexpr int kVoxelWidthRoot = 536870912;
        static constexpr int kVoxelDataWidthRoot = 268435456;
        static constexpr uint32_t kLog2Dim = 2;
        static constexpr uint32_t kLog2Dim2 = 2*kLog2Dim;
        static constexpr uint32_t kMask = (1u << kLog2Dim) - 1u; // mask for bit operations
        // kMaxCoarsenedPath = 3 means that a 4^3 tileData can at most be coalesced into a 32^3 tileData.
        static constexpr int kMaxCoalescedPath = 3;

        static constexpr uint64_t kTileDataStrideShift = 48;
        static constexpr uint64_t kTileDataShiftShift = kTileDataStrideShift + 8;
        static constexpr uint64_t kTileDataOffsetMask = ((1ull << kTileDataStrideShift) - 1ul);
        static constexpr uint64_t kTileDataShiftMask = ((~kTileDataOffsetMask) << 8ul);
        static constexpr uint64_t kTileDataStrideMask = kTileDataShiftMask >> 8ul;

        static constexpr uint64_t k3DTextureShift = 16;
        static constexpr uint64_t k3DTextureMask = ((1ull << k3DTextureShift) - 1ul);

        __hostdev__ static float invTileSideLUT(const int in_depth) {
#ifdef __CUDA_ARCH__
            return kInvTileSideLUTDevice[in_depth];
#else
            return kInvTileSideLUT[in_depth];
#endif
        }

        __hostdev__ static int tileSideLUT(const int in_depth) {
#ifdef __CUDA_ARCH__
            return kTileSideLUTDevice[in_depth];
#else
            return kTileSideLUT[in_depth];
#endif
        }

        __hostdev__ static int voxelDataWidth(const int in_depth) {
            return tileSideLUT(in_depth) / kTileDataWidth;
        }

        __hostdev__ static int voxelWidth(const int in_depth)
        {
            return tileSideLUT(in_depth + 1);
        }

        __hostdev__ static int tileVoxelWidth(const int in_depth)
        {
            return tileSideLUT(in_depth);
        }
    };

    static constexpr TileIndex kInvalidTileIndex = static_cast<TileIndex>(-1);
    static constexpr int kInvalidDepth = -1;
    static constexpr int kTileStrideTableSize = 4;

#define kTILESTRIDEINDEX(in_name)                                        \
    int8_t in_name[] = {Traits::kTileDataWidth, Traits::kTileDataWidth2, \
                        Traits::kExtendedTileDataWidth,                      \
                        Traits::kExtendedTileDataWidth2};
#ifdef __CUDACC__
    __constant__ static const kTILESTRIDEINDEX(kTileStrideDevice)
#endif
        static constexpr kTILESTRIDEINDEX(kTileStride)

        __hostdev__ static int8_t tileStride(uint64_t in_index)
    {
#ifdef __CUDA_ARCH__
        return kTileStrideDevice[in_index];
#else
        return kTileStride[in_index];
#endif
    }

    template<typename BuildT, typename dummy = void> class Tree;
    template<typename TreeT> class Grid;

    // Explicit alignment not necessary (will be aligned to 8 or 16 bytes)
    struct TreeIndex
    {
        TileIndex m_tile;
        int m_depth;

        __hostdev__ TreeIndex(const TileIndex in_tile = kInvalidTileIndex, const int in_depth = kInvalidDepth) : m_tile(in_tile), m_depth(in_depth) {}
        __hostdev__ bool valid() const { return m_tile != kInvalidTileIndex && m_depth != kInvalidDepth; }
        __hostdev__ bool operator==(const TreeIndex& in) const { return m_tile == in.m_tile && m_depth == in.m_depth; }
        __hostdev__ bool operator!=(const TreeIndex& in) const { return m_tile != in.m_tile || m_depth != in.m_depth; }
    };


    // Explicit alignment not necessary (will be aligned to 16 bytes)
    struct VoxelCoord
    {
        nanovdb::Coord m_ijk;
        int m_depth;

        __hostdev__ VoxelCoord() : m_ijk({ 0,0,0 }), m_depth(kInvalidDepth) {}
        __hostdev__ VoxelCoord(const nanovdb::Coord& in_ijk, const int in_depth)
            : m_ijk(in_ijk), m_depth(in_depth) {}
    };

    struct VoxelCoordLessThan
    {
        /// \return true if in_t1 < in_t2
        __hostdev__ bool operator()(const VoxelCoord& in_t1, const VoxelCoord& in_t2) const
        {
            return (in_t1.m_ijk[0] < in_t2.m_ijk[0] || (in_t1.m_ijk[0] == in_t2.m_ijk[0] && (in_t1.m_ijk[1] < in_t2.m_ijk[1] || (in_t1.m_ijk[1] == in_t2.m_ijk[1] && in_t1.m_ijk[2] < in_t2.m_ijk[2]))))
                || (in_t1.m_ijk[0] == in_t2.m_ijk[0] && in_t1.m_ijk[1] == in_t2.m_ijk[1] && in_t1.m_ijk[2] == in_t2.m_ijk[2] && in_t1.m_depth < in_t2.m_depth);
        }
    };

    struct Tile
    {
        nanovdb::Coord m_tijk; // 12 bytes
        int m_depth;           // 4 bytes
        TileIndex m_tile;        // 4 bytes
        TileIndex m_parent;      //  4 bytes
        TileIndex m_children[Traits::kTileSize]; // 8 * 8 = 64 bytes

        __hostdev__ bool valid() const { return m_tile != kInvalidTileIndex && m_depth != kInvalidDepth; }

        __hostdev__ nanovdb::Coord const& tijk() const { return m_tijk; }

        __hostdev__ int depth() const { return m_depth; }

        __hostdev__ bool atMaxDepth() const { return m_depth == Traits::kMaxDepth; }

        __hostdev__ TileIndex tile() const { return m_tile; }

        __hostdev__ TileIndex parent() const { return m_parent; }

        __hostdev__ nanovdb::Coord childIJK(
            const nanovdb::Coord& in_tijk) const {
            const int vw = Traits::voxelWidth(m_depth);
            return nanovdb::Coord{(in_tijk[0] - m_tijk[0]) / vw,
                                  (in_tijk[1] - m_tijk[1]) / vw,
                                  (in_tijk[2] - m_tijk[2]) / vw};
        }

        __hostdev__ TileIndex childIndex(const nanovdb::Coord& in_ijk) const {
            return m_children[in_ijk[0] + in_ijk[1] * Traits::kTileWidth +
                              in_ijk[2] * Traits::kTileWidth2];
        }

        __hostdev__ TileIndex childIndex(const int in_i) const {
            return m_children[in_i];
        }

        __hostdev__ bool isLeaf() const {
            for (int i = 0; i < Traits::kTileSize; i++) {
                if (m_children[i] != kInvalidTileIndex) {
                    return false;
                }
            }
            return true;
        }
    };

    // The container is stored first in the buffer
    //
    // Limitations
    // - There's just a single container associated with all grids which means that
    //   all grids must have the same topology.
    class NANOVDB_ALIGN(NANOVDB_DATA_ALIGNMENT) Container
    {
	public:
        // Size of the whole NanoVTT buffer in bytes
        uint64_t m_bufferSize;
        // Size of this container in bytes
        uint64_t m_containerSize;
        // The tileCount at a certain depth is computed as m_tileCountAccum[d+1] - m_tileCountAccum[d]
        // The total tile count is computed as m_tileCountAccum[Traits::kMaxDepth + 1]
        TileIndex m_tileCountAccum[Traits::kMaxDepth + 2];

        // Smaller table sizes typically result in worse performance.
        static constexpr int kTableBits = 6;
        static constexpr int kTableBits2 = 2 * kTableBits;
        static constexpr int kTableWidth = (1<<kTableBits);
        static constexpr int kTableSize = kTableWidth * kTableWidth * kTableWidth;
        TileIndex m_table[kTableSize];
        nanovdb::Coord m_tableMin;
        int m_rootDepth;
        uint64_t m_parentOffset;
        uint8_t m_tiles[];

        ///////////////////////////////////////////////////////////
        // NON-STANDARD API BEGIN
        ///////////////////////////////////////////////////////////

        __hostdev__ TileIndex activeTileCount() const
        {
            return m_tileCountAccum[Traits::kMaxDepth + 1];
        }

        __hostdev__ TileIndex activeTileCount(uint32_t in_level) const
		{
			NANOVDB_ASSERT(in_level <= Traits::kMaxDepth);
			return (m_tileCountAccum[in_level + 1] - m_tileCountAccum[in_level]);
		}
	
		__hostdev__ uint64_t activeVoxelCount(uint32_t in_level) const
        {
            NANOVDB_ASSERT(in_level <= Traits::kMaxDepth);
            return static_cast<uint64_t>(m_tileCountAccum[in_level + 1] - m_tileCountAccum[in_level]) * static_cast<uint64_t>(Traits::kTileDataSize);
        }

        __hostdev__ uint64_t activeVoxelCount() const
        {
            return static_cast<uint64_t>(m_tileCountAccum[Traits::kMaxDepth + 1] - m_tileCountAccum[0]) * static_cast<uint64_t>(Traits::kTileDataSize);
        }

        __hostdev__ TileIndex nodeCount(uint32_t in_level) const
        {
            NANOVDB_ASSERT(in_level <= Traits::kMaxDepth);
            return (m_tileCountAccum[in_level + 1] - m_tileCountAccum[in_level]);
        }

        __hostdev__ TileIndex nodeCount() const
        {
            return m_tileCountAccum[Traits::kMaxDepth + 1];
        }

        __hostdev__  TileIndex root(nanovdb::Coord const& in_ijk, nanovdb::Coord& io_tijk) const
        {
            assert(in_ijk[0] >= 0 && in_ijk[0] < kTableWidth && in_ijk[1] >= 0 && in_ijk[1] < kTableWidth && in_ijk[2] >= 0 && in_ijk[2] < kTableWidth);
            const int shift = Traits::kMaxDepth - m_rootDepth + 2;
            io_tijk = { m_tableMin[0] + (in_ijk[0] << shift),
                       m_tableMin[1] + (in_ijk[1] << shift),
                       m_tableMin[2] + (in_ijk[2] << shift) };
            return  m_table[(in_ijk[2] << kTableBits2) + (in_ijk[1] << kTableBits) + in_ijk[0]];
        }

        __hostdev__  TileIndex root(nanovdb::Coord const& in_ijk) const
        {
            assert(in_ijk[0] >= 0 && in_ijk[0] < kTableWidth && in_ijk[1] >= 0 && in_ijk[1] < kTableWidth && in_ijk[2] >= 0 && in_ijk[2] < kTableWidth);
            return  m_table[(in_ijk[2] << kTableBits2) + (in_ijk[1] << kTableBits) + in_ijk[0]];
        }

        __hostdev__  void addRoot(nanovdb::Coord const& in_key, TileIndex in_index)
        {
            // Clamp in_key to tile space at the root depth.
            const int shift = Traits::kMaxDepth - m_rootDepth + 2;
            const int i = (in_key[0] - m_tableMin[0]) >> shift;
            const int j = (in_key[1] - m_tableMin[1]) >> shift;
            const int k = (in_key[2] - m_tableMin[2]) >> shift;
            assert(i >= 0 && j >= 0 && k >= 0 && i < kTableWidth&& j < kTableWidth&& k < kTableWidth);
            m_table[(k << kTableBits2) + (j << kTableBits) + i] = in_index;
        }

        __hostdev__  TileIndex findRoot(nanovdb::Coord const& in_key, nanovdb::Coord& io_tijk) const
        {
            // Clamp in_key to tile space at the root depth.
            const int shift = Traits::kMaxDepth - m_rootDepth + 2;
            const int i = (in_key[0] - m_tableMin[0]) >> shift;
            const int j = (in_key[1] - m_tableMin[1]) >> shift;
            const int k = (in_key[2] - m_tableMin[2]) >> shift;
            const TileIndex hash = static_cast<TileIndex>(i) |
                static_cast<TileIndex>(j) |
                static_cast<TileIndex>(k);
            io_tijk = { m_tableMin[0] + (i << shift),
                       m_tableMin[1] + (j << shift),
                       m_tableMin[2] + (k << shift) };
            return (hash < kTableWidth
                ? m_table[(k << kTableBits2) + (j << kTableBits) + i]
                : kInvalidTileIndex);
        }

        __hostdev__  TileIndex findRoot(nanovdb::Coord const& in_key) const
        {
            // Clamp in_key to tile space at the root depth.
            const int shift = Traits::kMaxDepth - m_rootDepth + 2;
            const int i = (in_key[0] - m_tableMin[0]) >> shift;
            const int j = (in_key[1] - m_tableMin[1]) >> shift;
            const int k = (in_key[2] - m_tableMin[2]) >> shift;
            const TileIndex       hash = static_cast<TileIndex>(i) |
                static_cast<TileIndex>(j) |
                static_cast<TileIndex>(k);
            return (hash < kTableWidth ? m_table[(k << kTableBits2) + (j << kTableBits) + i] : kInvalidTileIndex);
        }

        __hostdev__ int64_t memUsage() const { return m_containerSize; }

        __hostdev__ int64_t toLinearIndex(const TreeIndex& in_treeIndex) const
        {
            return in_treeIndex.m_tile;
        }

        __hostdev__  int rootDepth() const
        {
            return m_rootDepth;
        }

        __hostdev__  void setTileOffset(const TileIndex in_index, const uint64_t in_offset)
        {
            reinterpret_cast<uint64_t*>(m_tiles)[in_index] = in_offset;
        }

        __hostdev__ void setTileShift(const TileIndex in_tile, const uint8_t in_shift)
        {
            reinterpret_cast<uint64_t*>(m_tiles)[in_tile] |= (static_cast<uint64_t>(in_shift) << Traits::kTileDataShiftShift);
        }

#ifdef NANOVTT_SUPPORT_COALESCED_TILES
        __hostdev__ uint8_t tileShift(const TileIndex in_tile) const
        {
            assert(in_tile <= m_tileCountAccum[Traits::kMaxDepth]);
            return static_cast<uint8_t>(1 + ((reinterpret_cast<uint64_t const*>(m_tiles)[in_tile]) >> Traits::kTileDataShiftShift));
        }
#else
        __hostdev__ uint8_t tileShift(const TileIndex /*in_tile*/) const
        {
            return 1;
        }
#endif

        __hostdev__  TileIndex childIndex(const TileIndex in_index, const nanovdb::Coord& in_ijk) const
        {
            // To avoid branching here we store one extra tile at the end with invalid child indices.
#ifdef __CUDACC__
            const uint64_t i1 = reinterpret_cast<const uint64_t*>(m_tiles)[Min(in_index, m_tileCountAccum[Traits::kMaxDepth])];
#else
            const uint64_t i1 = reinterpret_cast<const uint64_t*>(m_tiles)[Min(in_index, m_tileCountAccum[Traits::kMaxDepth])];
#endif

#ifdef NANOVTT_SUPPORT_COALESCED_TILES
            const int shift = static_cast<int>(1 + (i1 >> Traits::kTileDataShiftShift));
            return reinterpret_cast<const TileIndex*>(&(m_tiles[i1 & Traits::kTileDataOffsetMask]))[in_ijk[0] + (in_ijk[1] << shift) + (in_ijk[2] << (shift << 1))];
#else
            return reinterpret_cast<const TileIndex*>(&(m_tiles[i1]))[in_ijk[0] + (in_ijk[1] << 1) + (in_ijk[2] << 2)];
#endif
        }

        __hostdev__  int childIndexAndShift(const TileIndex in_index, const nanovdb::Coord& in_diff, int& io_vw, nanovdb::Coord& io_ijk, TileIndex& out_index) const
        {
            // To avoid branching here we store one extra tile at the end with invalid child indices.
#ifdef __CUDACC__
            const uint64_t i1 = reinterpret_cast<const uint64_t*>(m_tiles)[Min(in_index, m_tileCountAccum[Traits::kMaxDepth])];
#else
            const uint64_t i1 = reinterpret_cast<const uint64_t*>(m_tiles)[Min(in_index, m_tileCountAccum[Traits::kMaxDepth])];
#endif
#ifdef NANOVTT_SUPPORT_COALESCED_TILES
            const int shift = 1 + static_cast<int>(i1 >> Traits::kTileDataShiftShift);
            io_vw += 1 - shift;
            io_ijk[0] = in_diff[0] >> io_vw;
            io_ijk[1] = in_diff[1] >> io_vw;
            io_ijk[2] = in_diff[2] >> io_vw;
            out_index = reinterpret_cast<const TileIndex*>(&(m_tiles[i1 & Traits::kTileDataOffsetMask]))[io_ijk[0] + (io_ijk[1] << shift) + (io_ijk[2] << (shift << 1))];
            return shift;
#else
            io_ijk[0] = in_diff[0] >> io_vw;
            io_ijk[1] = in_diff[1] >> io_vw;
            io_ijk[2] = in_diff[2] >> io_vw;
            out_index = reinterpret_cast<const TileIndex*>(&(m_tiles[i1]))[io_ijk[0] + (io_ijk[1] << 1) + (io_ijk[2] << 2)];
            return 1;
#endif
        }

        __hostdev__  TileIndex* children(const TileIndex in_index)
        {
            return reinterpret_cast<TileIndex*>(&(m_tiles[reinterpret_cast<uint64_t*>(m_tiles)[in_index] & Traits::kTileDataOffsetMask]));
        }

        __hostdev__  TileIndex const* children(const TileIndex in_index) const
        {
            return reinterpret_cast<const TileIndex*>(&(m_tiles[reinterpret_cast<const uint64_t*>(m_tiles)[in_index] & Traits::kTileDataOffsetMask]));
        }

        __hostdev__  TileIndex parent(const TileIndex in_index) const
        {
            return reinterpret_cast<const TileIndex*>(&(m_tiles[m_parentOffset]))[in_index];
        }

        __hostdev__  int parentAndShift(TileIndex& io_index) const
        {
            io_index = reinterpret_cast<const TileIndex*>(&(m_tiles[m_parentOffset]))[io_index];
            return tileShift(io_index);
        }

        __hostdev__  void setParent(const TileIndex in_index, const TileIndex in_parent)
        {
            reinterpret_cast<TileIndex*>(&(m_tiles[m_parentOffset]))[in_index] = in_parent;
        }

        ///////////////////////////////////////////////////////////
        // NON-STANDARD API END
        ///////////////////////////////////////////////////////////
    };

    template<typename BuildT>
    class ReadAccessor
    {
    public:
        using BuildType = BuildT;
        using ValueType = typename nanovdb::BuildToValueMap<BuildType>::Type;
        using CoordType = nanovdb::Coord;
        using CoordValueType = typename CoordType::ValueType;
        using RootType = Tree<BuildType>;

    private:
        // Note: here it's ok to store a pointer to the grid because the accessor is instantiated on either the device or host.
        Container const* const m_container;
        Tree<BuildType> const* const m_tree;
        const float m_voxelSize;
        const int m_rootDepth;
        const bool m_isLevelSet;

        // mutable variables
        // Note: we can reduce this if we impose a lower limit on the depth
        // Implicitly stores the parent.
        // We allocate an extra entry to avoid branching in getValue().
        mutable int m_depth;
        mutable int m_vw;
        mutable CoordType m_tijk;
        mutable CoordType m_ijk;
        mutable bool m_found;
        mutable TileIndex m_tile;
        mutable int m_tileDataStride;

    public:

        __hostdev__ ReadAccessor(Container const* const in_container,
            Tree<BuildType> const* const in_tree);

        /// @brief Return the voxel value at the given coordinate.
        __hostdev__ ValueType getValue(const CoordType& in_ijk) const;

        template<typename RayT>
        __hostdev__ uint32_t getDim(const CoordType& in_ijk, const RayT& /*ray*/, const int in_maxDepth = Traits::kMaxDepth) const;

        /// @brief Return the active state of the given voxel (regardless of state or location in the tree.)
        __hostdev__ inline bool isActive(const CoordType& in_ijk) const;

        __hostdev__ const RootType& root() const;

        ///////////////////////////////////////////////////////////
        // NON-STANDARD API BEGIN
        ///////////////////////////////////////////////////////////

        __hostdev__ uint32_t getDim() const;

        __hostdev__ CoordType& ijk() const { return m_ijk; }

        /// \return The shift corresponding to the depth of the tiledata. The depth of the tiledata is Traits::kMaxDepth - vw().
        __hostdev__ int vw() const { return m_vw; }

        __hostdev__ TileIndex tileIndex() const { return m_tile; }

        __hostdev__ bool isExtended() const { return m_tree->isExtended(m_tile); }

        __hostdev__ CoordType const& tijk() const { return m_tijk; }

        ///  \return The depth of the tile (the tiledata may be have a resolution corresponding to a finer depth).
        __hostdev__ int depth() const { return m_depth; }

        __hostdev__ int tileDataStride() const { return m_tileDataStride; }

        __hostdev__ bool found() const { return m_found; }

        __hostdev__ ValueType getOffsetExtendedValue(const CoordType& in_off) const;

        __hostdev__ void getStencilExtendedValues(ValueType(&v)[2][2][2]) const;

        __hostdev__ ValueType getExtendedValue(const CoordType& in_ijk, const CoordType& in_off, const int in_maxDepth) const;

        __hostdev__ void cache(const CoordType& in_ijk, const int in_maxDepth) const;

        __hostdev__ void cacheTile(const CoordType& in_ijk, const int in_maxDepth) const;

        /// @brief Return the voxel value at the given coordinate.
        __hostdev__ ValueType getValue(const CoordType& in_ijk, const int in_maxDepth) const;

        /// @brief Return the voxel value at the given coordinate.
        __hostdev__ ValueType getExistingValueSameDepth(const CoordType& in_ijk) const;

        __hostdev__ Tree<BuildType> const* getTree() const {
            return m_tree;
        }

    protected:

        __hostdev__ bool isActiveLevelSetOrFog(const CoordType& in_ijk) const;

        ///////////////////////////////////////////////////////////
        // NON-STANDARD API END
        ///////////////////////////////////////////////////////////
    };

    /// \brief The Tree stores the data of a grid.
    /// Note: does not suppport getValue() and isActive(). These are only supported on the accessor.
    template<typename BuildT>
    class NANOVDB_ALIGN(NANOVDB_DATA_ALIGNMENT) TreeBase
    {
	public:
        using BuildType = BuildT;
        using ValueType = typename nanovdb::BuildToValueMap<BuildType>::Type;
        using CoordType = nanovdb::Coord;
        using AccessorType = ReadAccessor<BuildType>;

        nanovdb::BBox<CoordType> m_indexBBox;
        // The container is offset by 'm_containerOffset' from this Tree
        uint64_t m_containerOffset;
        uint64_t m_treeSize;
        ValueType m_defaultValue;

        /// @brief Return true if this tree is empty, i.e. contains no values or nodes
        __hostdev__ bool isEmpty() const { return getContainer()->activeVoxelCount() == 0; }

        /// @brief Return the total number of active voxels in this tree.
        __hostdev__ uint64_t activeVoxelCount() const { return getContainer()->activeVoxelCount(); }

        // @brief Return the memory footprint in bytes of this Tree
        __hostdev__ uint64_t memUsage() const { return m_treeSize; }

        /// @brief   Return the total number of active tiles at the specified level of the tree.
         ///
         /// @details n = 0 corresponds to leaf level tiles.
        __hostdev__ TileIndex activeTileCount(uint32_t in_level) const
        {
            return getContainer()->activeTileCount(in_level);
        }

        __hostdev__ uint32_t nodeCount(uint32_t in_level) const
        {
            return getContainer()->nodeCount(in_level);
        }

        /// @brief Return a const reference to the index bounding box of all the active values in this tree, i.e. in all nodes of the tree
        __hostdev__ const nanovdb::BBox<CoordType>& bbox() const { return m_indexBBox; }

        ///////////////////////////////////////////////////////////
        // NON-STANDARD API BEGIN
        ///////////////////////////////////////////////////////////

        __hostdev__ TileIndex activeTileCount() const
        {
            return getContainer()->activeTileCount();
        }

        __hostdev__ int toLinearIndex(const nanovdb::Coord& in_ijk) const
        {
            return in_ijk[0] + in_ijk[1] * Traits::kTileDataWidth + in_ijk[2] * Traits::kTileDataWidth2;
        }

        __hostdev__ const Container* getContainer() const
        {
            return reinterpret_cast<Container const*>(reinterpret_cast<uint8_t const*>(this) - m_containerOffset);
        }

        ///////////////////////////////////////////////////////////
        // NON-STANDARD API END
        ///////////////////////////////////////////////////////////
    };

    struct /*NANOVDB_ALIGN(NANOVDB_DATA_ALIGNMENT)*/ TileDataFnBase
    {
        float m_minimum; //  4B - minimum of ALL values in this node
        float m_quantum; //  = (max - min)/15 4B

        void init(const float in_min, const float in_max, const uint8_t in_bitWidth)
        {
            m_minimum = in_min;
            m_quantum = (in_max - in_min) / static_cast<float>((1 << in_bitWidth) - 1);
        }
    };

    template<typename ArrayT, uint32_t bitWidthT>
    struct TileDataFixed : public TileDataFnBase
    {
        static constexpr uint32_t kBitWidth = bitWidthT;
        static constexpr uint64_t kTileDataBytes = Traits::kTileDataSize * kBitWidth / (sizeof(ArrayT) * 8) * sizeof(ArrayT);
        // Here we add one because the number of bits is not divisible by 8.
        static constexpr uint64_t kExtendedTileDataBytes = (Traits::kExtendedTileDataSize * kBitWidth / (sizeof(ArrayT) * 8) + 1) * sizeof(ArrayT);
        using ArrayType = ArrayT;
        // Note that this array is stored after TileDataFixed because it's not
        // possible to store an array of structs with zero-sized arrays.
        //ArrayT m_code[tileDataSizeT * kBitWidth / (sizeof(ArrayT) * 8)];  // quantized data

        __hostdev__ ArrayT* data()
        {
            return reinterpret_cast<ArrayT*>(this + 1);
        }

        template<typename CoordType>
        __hostdev__ float getValue(const CoordType& in_ijk, const uint64_t in_index) const
        {
            const int i = in_ijk[0] +
                in_ijk[1] * tileStride(in_index) +
                in_ijk[2] * tileStride(in_index + 1);

            // Ifs optimized out at compile time
            if (kBitWidth == 4) {
                return ((reinterpret_cast<ArrayT const*>(this + 1)[i >> 1] >> ((i & 1) << 2)) & uint8_t(15)) * TileDataFnBase::m_quantum + TileDataFnBase::m_minimum;
            }
            else if (kBitWidth == 8) {
                return reinterpret_cast<ArrayT const*>(this + 1)[i] * TileDataFnBase::m_quantum + TileDataFnBase::m_minimum;
            }
            else /*if (kBitWidth == 16)*/ {
                return reinterpret_cast<ArrayT const*>(this + 1)[i] * TileDataFnBase::m_quantum + TileDataFnBase::m_minimum;
            }
        }

        template<typename CoordType>
        __hostdev__ void getValues(const CoordType& in_ijk, const uint64_t in_index, float(&v)[2][2][2]) const
        {
            const int tds = tileStride(in_index);
            const int tds2 = tileStride(in_index + 1);
            const int i = in_ijk[0] +
                in_ijk[1] * tds +
                in_ijk[2] * tds2;
            ArrayT const* const td = reinterpret_cast<ArrayT const*>(this + 1);

            // Ifs optimized out at compile time
            if (kBitWidth == 4) {
                auto decode = [&](const int in) -> float {
                    return ((td[in >> 1] >> ((in & 1) << 2)) & uint8_t(15)) *
                               TileDataFnBase::m_quantum +
                           TileDataFnBase::m_minimum;
                };
                v[0][0][0] = decode(i);
                v[1][0][0] = decode(i + 1);
                v[0][1][0] = decode(i + tds);
                v[1][1][0] = decode(i + 1 + tds);
                v[0][0][1] = decode(i + tds2);
                v[1][0][1] = decode(i + 1 + tds2);
                v[0][1][1] = decode(i + tds + tds2);
                v[1][1][1] = decode(i + 1 + tds + tds2);
            } else if (kBitWidth == 8) {
                auto decode = [&](const int in) -> float {
                    return td[in] * TileDataFnBase::m_quantum +
                           TileDataFnBase::m_minimum;
                };
                v[0][0][0] = decode(i);
                v[1][0][0] = decode(i + 1);
                v[0][1][0] = decode(i + tds);
                v[1][1][0] = decode(i + 1 + tds);
                v[0][0][1] = decode(i + tds2);
                v[1][0][1] = decode(i + 1 + tds2);
                v[0][1][1] = decode(i + tds + tds2);
                v[1][1][1] = decode(i + 1 + tds + tds2);
            } else /*if (kBitWidth == 16)*/ {
                auto decode = [&](const int in) -> float {
                    return td[in] * TileDataFnBase::m_quantum +
                           TileDataFnBase::m_minimum;
                };
                v[0][0][0] = decode(i);
                v[1][0][0] = decode(i + 1);
                v[0][1][0] = decode(i + tds);
                v[1][1][0] = decode(i + 1 + tds);
                v[0][0][1] = decode(i + tds2);
                v[1][0][1] = decode(i + 1 + tds2);
                v[0][1][1] = decode(i + tds + tds2);
                v[1][1][1] = decode(i + 1 + tds + tds2);
            }
        }

        __hostdev__ uint32_t bitWidth() const
        {
            return kBitWidth;
        }

    };

    struct TileDataVariable : public TileDataFnBase
    {
        int16_t m_logBitWidth;
        int16_t m_isExtended;
        // Note that this array is stored after TileDataVariable because it's not
        // possible to store an array of structs with zero-sized arrays.
        //uint32_t m_code[];

        __hostdev__ uint8_t* data()
        {
            return reinterpret_cast<uint8_t*>(this + 1);
        }

        template<typename CoordType>
        __hostdev__ float getValue(const CoordType& in_ijk) const
        {
            const int16_t index = static_cast<int16_t>(m_isExtended << 1);
            const int i = in_ijk[0] +
                in_ijk[1] * tileStride(index) +
                in_ijk[2] * tileStride(index + 1);

            const int b = m_logBitWidth; // b = 0, 1, 2, 3, 4 corresponding to 1, 2, 4, 8, 16 bits
            uint32_t code = reinterpret_cast<uint32_t const*>(this + 1)[i >> (5 - b)];
            code >>= (i & ((32 >> b) - 1)) << b;
            code &= (1 << (1 << b)) - 1;
            return float(code) * TileDataFnBase::m_quantum + TileDataFnBase::m_minimum;// code * (max-min)/UNITS + min
        }

        template<typename CoordType>
        __hostdev__ void getValues(const CoordType& in_ijk, float (&v)[2][2][2]) const
        {
            const int16_t index = static_cast<int16_t>(m_isExtended << 1);
            const int tds = tileStride(index);
            const int tds2 = tileStride(index + 1);
            const int i = in_ijk[0] +
                in_ijk[1] * tds +
                in_ijk[2] * tds2;
            const int b = m_logBitWidth; // b = 0, 1, 2, 3, 4 corresponding to 1, 2, 4, 8, 16 bits
            const int mask1 = ((32 >> b) - 1);
            const int mask2 = (1 << (1 << b)) - 1;
            const int shift = (5 - b);
            uint32_t const* const td = reinterpret_cast<uint32_t const*>(this + 1);

            auto decode = [&](const int in) -> float
            {
                uint32_t code = td[in >> shift];
                code >>= (in & mask1) << b;
                code &= mask2;
                return float(code) * TileDataFnBase::m_quantum + TileDataFnBase::m_minimum;// code * (max-min)/UNITS + min
            };

            v[0][0][0] = decode(i);
            v[1][0][0] = decode(i + 1);
            v[0][1][0] = decode(i + tds);
            v[1][1][0] = decode(i + 1 + tds);
            v[0][0][1] = decode(i + tds2);
            v[1][0][1] = decode(i + 1 + tds2);
            v[0][1][1] = decode(i + tds + tds2);
            v[1][1][1] = decode(i + 1 + tds + tds2);
        }

        __hostdev__ uint32_t bitWidth() const
        {
            return (1 << m_logBitWidth);
        }
    };

    template <typename BuildT>
    class NANOVDB_ALIGN(NANOVDB_DATA_ALIGNMENT)
        Tree<BuildT,
             typename std::enable_if<
                 !std::is_same<BuildT, nanovdb::Fp4>::value &&
                     !std::is_same<BuildT, nanovdb::Fp8>::value &&
                     !std::is_same<BuildT, nanovdb::Fp16>::value &&
                     !std::is_same<BuildT, nanovdb::FpN>::value,
                 void>::type> : public TreeBase<BuildT> {
	public:
        using BuildType = typename TreeBase<BuildT>::BuildType;
        using ValueType = typename TreeBase<BuildT>::ValueType;
        using CoordType = typename TreeBase<BuildT>::CoordType;
        using AccessorType = typename TreeBase<BuildT>::AccessorType;
        static constexpr uint64_t kTileDataBytes = Traits::kTileDataSize * sizeof(ValueType);
        static constexpr uint64_t kExtendedTileDataBytes = Traits::kExtendedTileDataSize * sizeof(ValueType);

        // Since we can only store a single zero-sized array we pack the data as follows:
        // uint64_t tileDataOffset[numTiles + 1]
        // ValueType m_data[]

        __hostdev__ uint8_t* data() { return reinterpret_cast<uint8_t*>(this) + sizeof(TreeBase<BuildT>); }

        __hostdev__ const uint8_t* data() const { return reinterpret_cast<const uint8_t*>(this) + sizeof(TreeBase<BuildT>); }

        /// @brief Return a new instance of a ReadAccessor used to access values in this grid
        __hostdev__ AccessorType getAccessor() const
        {
            return AccessorType(this->getContainer(), this);
        }

        ///////////////////////////////////////////////////////////
        // NON-STANDARD API BEGIN
        ///////////////////////////////////////////////////////////

        __hostdev__ bool isExtended(const TileIndex in_tile) const
        {
            return (tileDataStride(in_tile) & 1) > 0;
        }

        __hostdev__ void setTileDataOffset(const TileIndex in_tile, const uint64_t in_offset)
        {
            reinterpret_cast<uint64_t*>(data())[in_tile] = in_offset;
        }

        __hostdev__ uint64_t tileDataOffset(const TileIndex in_tile) const
        {
            return (reinterpret_cast<uint64_t const*>(data())[in_tile]) & Traits::kTileDataOffsetMask;
        }

        __hostdev__ void setTileDataStride(const TileIndex in_tile, const uint8_t in_stride)
        {
            // To ensure all encodings fit within 64 bits
            assert(in_stride <= 33);
            reinterpret_cast<uint64_t*>(data())[in_tile] |= (static_cast<uint64_t>(in_stride) << Traits::kTileDataStrideShift);
        }

        __hostdev__ uint8_t tileDataStride(const TileIndex in_tile) const
        {
            return static_cast<uint8_t>(((reinterpret_cast<uint64_t const*>(data())[in_tile]) & Traits::kTileDataStrideMask) >> Traits::kTileDataStrideShift);
        }

        __hostdev__ void setTileDataShift(const TileIndex in_tile, const uint8_t in_shift)
        {
            reinterpret_cast<uint64_t*>(data())[in_tile] |= (static_cast<uint64_t>(in_shift) << Traits::kTileDataShiftShift);
        }

        __hostdev__ uint8_t tileDataShift(const TileIndex in_tile) const
        {
            // Masking not necessary here because it's the most significant bits
            return static_cast<uint8_t>(((reinterpret_cast<uint64_t const*>(data())[in_tile]) /* & Traits::kTileDataShiftMask*/) >> Traits::kTileDataShiftShift);
        }

        /// @brief Return the base address of the tile data
        __hostdev__ ValueType const* tileData() const
        {
            return reinterpret_cast<ValueType const*>(&data()[(this->activeTileCount() + 1) * sizeof(uint64_t)]);
        }

        __hostdev__ ValueType const* tileData(const TileIndex in_tile) const
        {
            return reinterpret_cast<ValueType const*>(&data()[tileDataOffset(in_tile)]);
        }

        __hostdev__ ValueType* tileData(const TileIndex in_tile)
        {
            return reinterpret_cast<ValueType*>(&data()[tileDataOffset(in_tile)]);
        }

        __hostdev__ ValueType const* tileDataFromOffset(const uint64_t in_offset) const
        {
            return reinterpret_cast<ValueType const*>(&data()[in_offset]);
        }

        __hostdev__ ValueType* tileDataFromOffset(const uint64_t in_offset)
        {
            return reinterpret_cast<ValueType*>(&data()[in_offset]);
        }

        __hostdev__ void tileDataOffsetAndStrideAndShift(const TileIndex in_tile, uint64_t& out_offset, int& out_stride, int& out_shift) const
        {
            const uint64_t v = reinterpret_cast<uint64_t const*>(data())[in_tile];
            out_offset = v & Traits::kTileDataOffsetMask;
            out_stride = static_cast<int>((v & Traits::kTileDataStrideMask) >> Traits::kTileDataStrideShift);
            out_shift = static_cast<int>((v /* & Traits::kTileDataShiftMask*/) >> Traits::kTileDataShiftShift);
        }

        __hostdev__ void tileDataStrideAndShift(const TileIndex in_tile, int& out_stride, int& out_shift) const
        {
            const uint64_t v = reinterpret_cast<uint64_t const*>(data())[in_tile];
            out_stride = static_cast<int>((v & Traits::kTileDataStrideMask) >> Traits::kTileDataStrideShift);
            out_shift = static_cast<int>((v /* & Traits::kTileDataShiftMask*/) >> Traits::kTileDataShiftShift);
        }

        __hostdev__ void tileDataOffsetAndStride(const TileIndex in_tile, uint64_t& out_offset, int& out_stride) const
        {
            const uint64_t v = reinterpret_cast<uint64_t const*>(data())[in_tile];
            out_offset = v & Traits::kTileDataOffsetMask;
            out_stride = static_cast<int>((v & Traits::kTileDataStrideMask) >> Traits::kTileDataStrideShift);
        }

        __hostdev__ ValueType const& getValue(const TileIndex in_tile, const CoordType& in_ijk, int& io_shift, int& io_stride, CoordType& io_ijk) const
        {
            uint64_t offset;
            int tileDataStride, tileDataShift;
            tileDataOffsetAndStrideAndShift(in_tile, offset, tileDataStride, tileDataShift);
            tileDataShift = io_shift - tileDataShift;
            io_ijk[0] = (in_ijk[0] >> tileDataShift);
            io_ijk[1] = (in_ijk[1] >> tileDataShift);
            io_ijk[2] = (in_ijk[2] >> tileDataShift);
            io_shift = tileDataShift;
            io_stride = tileDataStride;
            return tileDataFromOffset(offset)[io_ijk[0] +
                io_ijk[1] * tileDataStride +
                io_ijk[2] * tileDataStride * tileDataStride];
        }

        __hostdev__ ValueType const& getValue(const TileIndex in_tile, const CoordType& in_ijk) const
        {
            uint64_t offset;
            int tileDataStride;
            tileDataOffsetAndStride(in_tile, offset, tileDataStride);
            return tileDataFromOffset(offset)[in_ijk[0] +
                in_ijk[1] * tileDataStride +
                in_ijk[2] * tileDataStride * tileDataStride];
        }

        __hostdev__ ValueType& getValue(const TileIndex in_tile, const CoordType& in_ijk)
        {
            uint64_t offset;
            int tileDataStride;
            tileDataOffsetAndStride(in_tile, offset, tileDataStride);
            return tileDataFromOffset(offset)[in_ijk[0] +
                in_ijk[1] * tileDataStride +
                in_ijk[2] * tileDataStride * tileDataStride];
        }

        __hostdev__ nanovdb::Coord getTileTextureCoord(const TileIndex in_tile) const
        {
            uint64_t offset;
            int tileDataStride;
            tileDataOffsetAndStride(in_tile, offset, tileDataStride);
            offset = offset & Traits::kTileDataOffsetMask;
            nanovdb::Coord texCoord;
            texCoord.x() = static_cast<int>(offset & Traits::k3DTextureMask);
            offset = offset >> Traits::k3DTextureShift;
            texCoord.y() = static_cast<int>(offset & Traits::k3DTextureMask);
            offset = offset >> Traits::k3DTextureShift;
            texCoord.z() = static_cast<int>(offset & Traits::k3DTextureMask);
            return texCoord;
        }

        __hostdev__ void getValues(const TileIndex in_tile, const CoordType& in_ijk, ValueType(&v)[2][2][2]) const
        {
            uint64_t offset;
            int tileDataStride;
            tileDataOffsetAndStride(in_tile, offset, tileDataStride);
            const int tileDataStride2 = tileDataStride * tileDataStride;
            ValueType const* const td = tileDataFromOffset(offset);
            const int index = in_ijk[0] +
                in_ijk[1] * tileDataStride +
                in_ijk[2] * tileDataStride2;
            v[0][0][0] = td[index];
            v[1][0][0] = td[index + 1];
            v[0][1][0] = td[index + tileDataStride];
            v[1][1][0] = td[index + 1 + tileDataStride];
            v[0][0][1] = td[index + tileDataStride2];
            v[1][0][1] = td[index + 1 + tileDataStride2];
            v[0][1][1] = td[index + tileDataStride + tileDataStride2];
            v[1][1][1] = td[index + 1 + tileDataStride + tileDataStride2];
        }

        __hostdev__ const Grid<Tree<BuildType>>* getGrid() const
        {
            return reinterpret_cast<Grid<Tree<BuildType>> const*>(reinterpret_cast<uint8_t const*>(this) - sizeof(nanovdb::GridData));
        }

        ///////////////////////////////////////////////////////////
        // NON-STANDARD API END
        ///////////////////////////////////////////////////////////
    };

    template <typename BuildT>
    class NANOVDB_ALIGN(NANOVDB_DATA_ALIGNMENT) Tree<
        BuildT,
        typename std::enable_if<std::is_same<BuildT, nanovdb::Fp4>::value ||
                                    std::is_same<BuildT, nanovdb::Fp8>::value ||
                                    std::is_same<BuildT, nanovdb::Fp16>::value,
                                void>::type> : public TreeBase<BuildT>
    {
	public:
        using BuildType = typename TreeBase<BuildT>::BuildType;
        using ValueType = typename TreeBase<BuildT>::ValueType;
        using CoordType = typename TreeBase<BuildT>::CoordType;
        using AccessorType = typename TreeBase<BuildT>::AccessorType;
        using TileDataType = TileDataFixed<
            typename std::conditional<
                std::is_same<BuildT, nanovdb::Fp16>::value,
                uint16_t,
                uint8_t>::type,
            (std::is_same<BuildT, nanovdb::Fp16>::value
                 ? 16
                 : (std::is_same<BuildT, nanovdb::Fp8>::value ? 8 : 4))>;
        // Here we ensure 4-byte alignment
        static constexpr uint64_t kTileDataBytes = TileDataType::kTileDataBytes + sizeof(TileDataType) + 4 - (TileDataType::kTileDataBytes + sizeof(TileDataType)) % 4;
        static constexpr uint64_t kExtendedTileDataBytes = TileDataType::kExtendedTileDataBytes + sizeof(TileDataType) + 4 - (TileDataType::kExtendedTileDataBytes + sizeof(TileDataType)) % 4;

        // Since we can only store a single zero-sized array we pack the data as follows:
        // uint64_t tileDataOffset[numTiles + 1]
        // struct {
        //   TileDataFixed;
        //   uint8_t[numCodes];
        // }[numTiles]

        /// @brief Return a new instance of a ReadAccessor used to access values in this grid
        __hostdev__ AccessorType getAccessor() const
        {
            return AccessorType(this->getContainer(), this);
        }

        __hostdev__ uint8_t* data() { return reinterpret_cast<uint8_t*>(this) + sizeof(TreeBase<BuildT>); }

        __hostdev__ const uint8_t* data() const { return reinterpret_cast<const uint8_t*>(this) + sizeof(TreeBase<BuildT>); }

        ///////////////////////////////////////////////////////////
        // NON-STANDARD API BEGIN
        ///////////////////////////////////////////////////////////

        __hostdev__ bool isExtended(const TileIndex in_tile) const
        {
            return (tileDataOffset(in_tile + 1) - tileDataOffset(in_tile)) == kExtendedTileDataBytes;
        }

        __hostdev__ void setTileDataOffset(const TileIndex in_tile, const uint64_t in_offset)
        {
            reinterpret_cast<uint64_t*>(data())[in_tile] = in_offset;
        }

        __hostdev__ uint64_t tileDataOffset(const TileIndex in_tile) const
        {
            return reinterpret_cast<uint64_t const*>(data())[in_tile];
        }

        __hostdev__ TileDataType const& tileData(const TileIndex in_tile) const
        {
            return *reinterpret_cast<TileDataType const*>(&data()[tileDataOffset(in_tile)]);
        }

        __hostdev__ TileDataType& tileData(const TileIndex in_tile)
        {
            return *reinterpret_cast<TileDataType*>(&data()[tileDataOffset(in_tile)]);
        }

        __hostdev__ ValueType getValue(const TileIndex in_tile, const CoordType& in_ijk, const int in_shift, int& io_stride, CoordType& io_ijk) const
        {
            io_ijk[0] = in_ijk[0] >> in_shift;
            io_ijk[1] = in_ijk[1] >> in_shift;
            io_ijk[2] = in_ijk[2] >> in_shift;
            const uint64_t index = ((tileDataOffset(in_tile + 1) - tileDataOffset(in_tile)) / kExtendedTileDataBytes) << 1;
            io_stride = tileStride(index);
            return tileData(in_tile).getValue(io_ijk, index);
        }

        __hostdev__ void getValues(const TileIndex in_tile, const CoordType& in_ijk, ValueType(&v)[2][2][2]) const
        {
            const uint64_t index = ((tileDataOffset(in_tile + 1) - tileDataOffset(in_tile)) / kExtendedTileDataBytes) << 1;
            tileData(in_tile).getValues(in_ijk, index, v);
        }

        __hostdev__ ValueType getValue(const TileIndex in_tile, const CoordType& in_ijk) const
        {
            const uint64_t index = ((tileDataOffset(in_tile + 1) - tileDataOffset(in_tile)) / kExtendedTileDataBytes) << 1;
            return tileData(in_tile).getValue(in_ijk, index);
        }

        __hostdev__ const Grid<Tree<BuildType>>* getGrid() const
        {
            return reinterpret_cast<Grid<Tree<BuildType>> const*>(reinterpret_cast<uint8_t const*>(this) - sizeof(nanovdb::GridData));
        }

        __hostdev__ void setTileDataStride(const TileIndex /*in_tile*/, const uint8_t /*in_stride*/)
        {
            // NOP
        }

        __hostdev__ void setTileDataShift(const TileIndex /*in_tile*/, const uint8_t /*in_shift*/)
        {
            // NOP
        }

        __hostdev__ ValueType* tileDataFromOffset(const uint64_t /*in_offset*/)
        {
            // NOP
            return nullptr;
        }

        ///////////////////////////////////////////////////////////
        // NON-STANDARD API END
        ///////////////////////////////////////////////////////////
    };

    template<>
    class NANOVDB_ALIGN(NANOVDB_DATA_ALIGNMENT) Tree<nanovdb::FpN, void> : public TreeBase<nanovdb::FpN>
    {
	public:
        using BuildType = typename TreeBase<nanovdb::FpN>::BuildType;
        using ValueType = typename TreeBase<nanovdb::FpN>::ValueType;
        using CoordType = typename TreeBase<nanovdb::FpN>::CoordType;
        using AccessorType = typename TreeBase<nanovdb::FpN>::AccessorType;
        using TileDataType = TileDataVariable;
        static constexpr uint64_t kTileDataBytes = sizeof(TileDataType); // Does not include the storage of codes
        static constexpr uint64_t kExtendedTileDataBytes = sizeof(TileDataType); // Does not include the storage of codes

        // Since we can only store a single zero-sized array we pack the data as follows:
        // uint64_t tileDataOffset[numTiles]
        // struct {
        //   TileDataVariable;
        //   uint8_t[numCodes];
        // }[numTiles]
        __hostdev__ uint8_t* data() { return reinterpret_cast<uint8_t*>(this) + sizeof(TreeBase<nanovdb::FpN>); }

        __hostdev__ const uint8_t* data() const { return reinterpret_cast<const uint8_t*>(this) + sizeof(TreeBase<nanovdb::FpN>); }


        /// @brief Return a new instance of a ReadAccessor used to access values in this grid
        __hostdev__ AccessorType getAccessor() const
        {
            return AccessorType(this->getContainer(), this);
        }

        ///////////////////////////////////////////////////////////
        // NON-STANDARD API BEGIN
        ///////////////////////////////////////////////////////////

        __hostdev__ bool isExtended(const TileIndex in_tile) const
        {
            return tileData(in_tile).m_isExtended > 0;
        }

        __hostdev__ void setTileDataOffset(const TileIndex in_tile, const uint64_t in_offset)
        {
            reinterpret_cast<uint64_t*>(data())[in_tile] = in_offset;
        }

        __hostdev__ uint64_t tileDataOffset(const TileIndex in_tile) const
        {
            return reinterpret_cast<uint64_t const*>(data())[in_tile];
        }

        __hostdev__ TileDataType const& tileData(const TileIndex in_tile) const
        {
            return *reinterpret_cast<TileDataType const*>(&data()[tileDataOffset(in_tile)]);
        }

        __hostdev__ TileDataType& tileData(const TileIndex in_tile)
        {
            return *reinterpret_cast<TileDataType*>(&data()[tileDataOffset(in_tile)]);
        }

        __hostdev__ ValueType getValue(const TileIndex in_tile, const CoordType& in_ijk, const int in_shift, int& io_stride, CoordType& io_ijk) const
        {
            io_ijk[0] = in_ijk[0] >> in_shift;
            io_ijk[1] = in_ijk[1] >> in_shift;
            io_ijk[2] = in_ijk[2] >> in_shift;
            const auto& td = tileData(in_tile);
            io_stride = Traits::kTileDataWidth + td.m_isExtended;
            return td.getValue(io_ijk);
        }

        __hostdev__ void getValues(const TileIndex in_tile, const CoordType& in_ijk, ValueType(&v)[2][2][2]) const
        {
            tileData(in_tile).getValues(in_ijk, v);
        }

        __hostdev__ ValueType getValue(const TileIndex in_tile, const CoordType& in_ijk) const
        {
            return tileData(in_tile).getValue(in_ijk);
        }

        __hostdev__ const Grid<Tree<BuildType>>* getGrid() const
        {
            return reinterpret_cast<Grid<Tree<BuildType>> const*>(reinterpret_cast<uint8_t const*>(this) - sizeof(nanovdb::GridData));
        }

        __hostdev__ void setTileDataStride(const TileIndex /*in_tile*/, const uint8_t /*in_stride*/)
        {
            // NOP
        }

        __hostdev__ void setTileDataShift(const TileIndex /*in_tile*/, const uint8_t /*in_shift*/)
        {
            // NOP
        }

        __hostdev__ ValueType* tileDataFromOffset(const uint64_t /*in_offset*/)
        {
            // NOP
            return nullptr;
        }

        ///////////////////////////////////////////////////////////
        // NON-STANDARD API END
        ///////////////////////////////////////////////////////////
    };

    template<typename TreeT>
    class Grid : private nanovdb::GridData
    {
    public:
        using TreeType = TreeT;
        using RootType = TreeT;
        using DataType = nanovdb::GridData;
        using ValueType = typename TreeT::ValueType;
        using BuildType = typename TreeT::BuildType;// in rare cases BuildType != ValueType, e.g. then BuildType = ValueMask and ValueType = bool
        using CoordType = typename TreeT::CoordType;
        using AccessorType = ReadAccessor<BuildType>;

        /// @brief Disallow constructions, copy and assignment
        ///
        /// @note Only a Serializer, defined elsewhere, can instantiate this class
        Grid(const Grid&) = delete;
        Grid& operator=(const Grid&) = delete;
        ~Grid() = delete;

        __hostdev__ nanovdb::Version version() const { return DataType::mVersion; }

        __hostdev__ DataType* data() { return reinterpret_cast<DataType*>(this); }

        __hostdev__ const DataType* data() const { return reinterpret_cast<const DataType*>(this); }

        /// @brief Return memory usage in bytes for this class only.
        __hostdev__ static uint64_t memUsage() { return sizeof(GridData); }

        /// @brief Return the memory footprint of the entire grid, i.e. including all nodes and blind data
        __hostdev__ uint64_t gridSize() const { return DataType::mGridSize; }

        /// @brief Return index of this grid in the buffer
        __hostdev__ uint32_t gridIndex() const { return DataType::mGridIndex; }

        /// @brief Return total number of grids in the buffer
        __hostdev__ uint32_t gridCount() const { return DataType::mGridCount; }

        /// @brief Return a const reference to the tree
        __hostdev__ const TreeT& tree() const { return *reinterpret_cast<const TreeT*>(this->treePtr()); }

        /// @brief Return a non-const reference to the tree
        __hostdev__ TreeT& tree() { return *reinterpret_cast<TreeT*>(this->treePtr()); }

        /// @brief Return a new instance of a ReadAccessor used to access values in this grid
        __hostdev__ AccessorType getAccessor() const { return AccessorType(this->tree().getContainer(), &this->tree()); }

        /// @brief Return a const reference to the size of a voxel in world units
        __hostdev__ const nanovdb::Vec3R& voxelSize() const { return DataType::mVoxelSize; }

        /// @brief Return a const reference to the Map for this grid
        __hostdev__ const nanovdb::Map& map() const { return DataType::mMap; }

        /// @brief world to index space transformation
        template<typename Vec3T>
        __hostdev__ Vec3T worldToIndex(const Vec3T& xyz) const { return this->applyInverseMap(xyz); }

        /// @brief index to world space transformation
        template<typename Vec3T>
        __hostdev__ Vec3T indexToWorld(const Vec3T& xyz) const { return this->applyMap(xyz); }

        /// @brief transformation from index space direction to world space direction
        /// @warning assumes dir to be normalized
        template<typename Vec3T>
        __hostdev__ Vec3T indexToWorldDir(const Vec3T& dir) const { return this->applyJacobian(dir); }

        /// @brief transformation from world space direction to index space direction
        /// @warning assumes dir to be normalized
        template<typename Vec3T>
        __hostdev__ Vec3T worldToIndexDir(const Vec3T& dir) const { return this->applyInverseJacobian(dir); }

        /// @brief transform the gradient from index space to world space.
        /// @details Applies the inverse jacobian transform map.
        template<typename Vec3T>
        __hostdev__ Vec3T indexToWorldGrad(const Vec3T& grad) const { return this->applyIJT(grad); }

        /// @brief world to index space transformation
        template<typename Vec3T>
        __hostdev__ Vec3T worldToIndexF(const Vec3T& xyz) const { return this->applyInverseMapF(xyz); }

        /// @brief index to world space transformation
        template<typename Vec3T>
        __hostdev__ Vec3T indexToWorldF(const Vec3T& xyz) const { return this->applyMapF(xyz); }

        /// @brief transformation from index space direction to world space direction
        /// @warning assumes dir to be normalized
        template<typename Vec3T>
        __hostdev__ Vec3T indexToWorldDirF(const Vec3T& dir) const { return this->applyJacobianF(dir); }

        /// @brief transformation from world space direction to index space direction
        /// @warning assumes dir to be normalized
        template<typename Vec3T>
        __hostdev__ Vec3T worldToIndexDirF(const Vec3T& dir) const { return this->applyInverseJacobianF(dir); }

        /// @brief Transforms the gradient from index space to world space.
        /// @details Applies the inverse jacobian transform map.
        template<typename Vec3T>
        __hostdev__ Vec3T indexToWorldGradF(const Vec3T& grad) const { return DataType::applyIJTF(grad); }

        /// @brief Computes a AABB of active values in world space
        __hostdev__ const nanovdb::BBox<nanovdb::Vec3R>& worldBBox() const { return DataType::mWorldBBox; }

        /// @brief Computes a AABB of active values in index space
        ///
        /// @note This method is returning a floating point bounding box and not a CoordBBox. This makes
        ///       it more useful for clipping rays.
        __hostdev__ const nanovdb::BBox<CoordType>& indexBBox() const { return this->tree().bbox(); }

        /// @brief Return the total number of active voxels in this tree.
        __hostdev__ uint64_t activeVoxelCount() const { return this->tree().activeVoxelCount(); }

        /// @brief Methods related to the classification of this grid
        __hostdev__ bool  isValid() const { return DataType::mMagic == NANOVTT_MAGIC_NUMBER; }
        __hostdev__ const nanovdb::GridType& gridType() const { return DataType::mGridType; }
        __hostdev__ const nanovdb::GridClass& gridClass() const { return DataType::mGridClass; }
        __hostdev__ bool             isLevelSet() const { return DataType::mGridClass == nanovdb::GridClass::LevelSet; }
        __hostdev__ bool             isFogVolume() const { return DataType::mGridClass == nanovdb::GridClass::FogVolume; }
        __hostdev__ bool             isStaggered() const { return DataType::mGridClass == nanovdb::GridClass::Staggered; }
        __hostdev__ bool             isPointIndex() const { return DataType::mGridClass == nanovdb::GridClass::PointIndex; }
        __hostdev__ bool             isPointData() const { return DataType::mGridClass == nanovdb::GridClass::PointData; }
        __hostdev__ bool             isMask() const { return DataType::mGridClass == nanovdb::GridClass::Topology; }
        __hostdev__ bool             isUnknown() const { return DataType::mGridClass == nanovdb::GridClass::Unknown; }
        __hostdev__ bool             hasMinMax() const { return DataType::mFlags & static_cast<uint32_t>(nanovdb::GridFlags::HasMinMax); }
        __hostdev__ bool             hasBBox() const { return DataType::mFlags & static_cast<uint32_t>(nanovdb::GridFlags::HasBBox); }
        __hostdev__ bool             hasLongGridName() const { return DataType::mFlags & static_cast<uint32_t>(nanovdb::GridFlags::HasLongGridName); }
        __hostdev__ bool             hasAverage() const { return DataType::mFlags & static_cast<uint32_t>(nanovdb::GridFlags::HasAverage); }
        __hostdev__ bool             hasStdDeviation() const { return DataType::mFlags & static_cast<uint32_t>(nanovdb::GridFlags::HasStdDeviation); }
        __hostdev__ bool             isBreadthFirst() const { return DataType::mFlags & static_cast<uint32_t>(nanovdb::GridFlags::IsBreadthFirst); }

        /// @brief return true if the specified node type is layed out breadth-first in memory and has a fixed size.
        ///        This allows for sequential access to the nodes.
        template <typename NodeT>
        __hostdev__ bool isSequential() const { return true; }

        /// @brief return true if the specified node level is layed out breadth-first in memory and has a fixed size.
        ///        This allows for sequential access to the nodes.
        template <int LEVEL>
        __hostdev__ bool isSequential() const { return true; }

        /// @brief Return a c-string with the name of this grid
        __hostdev__ const char* gridName() const
        {
            if (this->hasLongGridName()) {
                const auto& metaData = this->blindMetaData(DataType::mBlindMetadataCount - 1);// always the last
                NANOVDB_ASSERT(metaData.mDataClass == nanovdb::GridBlindDataClass::GridName);
                return metaData.template getBlindData<const char>();
            }
            return DataType::mGridName;
        }

        /// @brief Return a c-string with the name of this grid, truncated to 255 characters
        __hostdev__ const char* shortGridName() const { return DataType::mGridName; }

        /// @brief Return checksum of the grid buffer.
        __hostdev__ uint64_t checksum() const { return DataType::mChecksum; }

        /// @brief Return true if this grid is empty, i.e. contains no values or nodes.
        __hostdev__ bool isEmpty() const { return this->tree().isEmpty(); }

        /// @brief Return the count of blind-data encoded in this grid
        __hostdev__ int blindDataCount() const { return DataType::mBlindMetadataCount; }

        /// @brief Return the index of the blind data with specified semantic if found, otherwise -1.
        __hostdev__ int findBlindDataForSemantic(nanovdb::GridBlindDataSemantic semantic) const;

        /// @brief Returns a const pointer to the blindData at the specified linear offset.
        ///
        /// @warning Point might be NULL and the linear offset is assumed to be in the valid range
        __hostdev__ const void* blindData(uint32_t n) const
        {
            if (DataType::mBlindMetadataCount == 0) {
                return nullptr;
            }
            NANOVDB_ASSERT(n < DataType::mBlindMetadataCount);
            return this->blindMetaData(n).template getBlindData<void>();
        }

        __hostdev__ const nanovdb::GridBlindMetaData& blindMetaData(int n) const { return *DataType::blindMetaData(n); }

    private:
        static_assert(sizeof(GridData) % NANOVDB_DATA_ALIGNMENT == 0, "sizeof(GridData) is misaligned");
    };

    template <typename BuildT>
    __hostdev__ ReadAccessor<BuildT>::ReadAccessor(
        Container const* const       in_container,
        Tree<BuildType> const* const in_tree)
        : m_container(in_container),
          m_tree(in_tree),
         m_voxelSize(static_cast<float>(in_tree->getGrid()->voxelSize()[0])),
         m_rootDepth(in_container->m_rootDepth),
         m_isLevelSet(in_tree->getGrid()->isLevelSet()),
         m_depth(in_container->m_rootDepth),
         m_vw(Traits::kMaxDepth - in_container->m_rootDepth),
         m_tijk(Traits::kTileDomainWidth, Traits::kTileDomainWidth, Traits::kTileDomainWidth),
         m_found(false),
         m_tile(0)
        , m_tileDataStride(Traits::kTileDataWidth)
    {
    }

        /// @brief Return the voxel value at the given coordinate.
    template <typename BuildT>
    __hostdev__ typename ReadAccessor<BuildT>::ValueType ReadAccessor<BuildT>::getValue(const CoordType& in_ijk) const {
        return getValue(in_ijk, Traits::kMaxDepth);
    }

    template <typename BuildT>
    __hostdev__ typename ReadAccessor<BuildT>::ValueType ReadAccessor<BuildT>::getOffsetExtendedValue(const CoordType& in_off) const
    {
        // Assumptions:
        // - We've called getValue() on the accessor such that the tile is correct.
        return m_tree->getValue(m_tile, m_ijk + in_off);
    }

    template <typename BuildT>
    __hostdev__ void ReadAccessor<BuildT>::getStencilExtendedValues(ValueType(&v)[2][2][2]) const
    {
        // Assumptions:
        // - We've called getValue() on the accessor such that the tile is correct.
        m_tree->getValues(m_tile, m_ijk, v);
    }

    template <typename BuildT>
    __hostdev__ typename ReadAccessor<BuildT>::ValueType ReadAccessor<BuildT>::getExtendedValue(const CoordType& in_ijk, const CoordType& in_off, const int in_maxDepth) const
    {
        // Assumptions:
        // - The tile exists (found is true).
        // - All non-graded cases have an extended tile
        // - The tile (in_ijk) for which we're calling getExtendedValue is not extended.

        /* Algorithm:
         - If sample point is inside an extended tile we just assemble directly from the extended tile.
         - If sample point is interior to tile of 3x3x3 voxels then just use tile itself to interpolate,
           otherwise we need to assemble extended points.
         - For each extended point :
            -Find containing tile.
            - If tile is at same depth, then just return the corresponding value.
            - If tile is at depth - 1 :
              -If extended point is inside interior of 3x3x3 voxels or extended tile, then just interpolate from that tile alone.
             - If not, then find all extended points assuming they are all at the same depth.
            - Once all extended points are assembled, interpolate stencil.
        */

        // The code below has been optimized for this case
        static_assert(Traits::kTileWidth == 2, "getExtendedValue requires Traits::kTileWidth == 2");

        ValueType value = getValue(in_ijk + in_off, in_maxDepth);
        // Here we need to take the shift into account when checking for depth.
        if (Traits::kMaxDepth - m_vw == in_maxDepth) {
            // Neighboring tile is at same depth, so no need to interpolate, we can just return value found.
            return value;
        }
        else {
            // Here we assume a graded tree
            // We know that resolution differs by two between levels so in the case of an extended point,
            // the offset will be 0 or 0.5. At even local indices it will be 0, and at odd indices it will be 0.5.
            // Note: the offsets of 0 or 0.5 could be used to avoid some lookups at the cost of more branching.
            const CoordType baseIJK = m_tijk + nanovdb::Coord(m_ijk[0] << m_vw, m_ijk[1] << m_vw, m_ijk[2] << m_vw);
            const nanovdb::Vec3f uvw(((baseIJK[0] >> (m_vw - 1)) & 1) * 0.5f, ((baseIJK[1] >> (m_vw - 1)) & 1) * 0.5f, ((baseIJK[2] >> (m_vw - 1)) & 1) * 0.5f);
            auto lerp = [](ValueType a, ValueType b, float w) {
                return a + w * (b - a);
            };
            // Linearly interpolate and return.
            // Here we must take the stride into account when checking upper bounds
            if (isExtended() || (m_ijk[0] < m_tileDataStride - 1 &&
                m_ijk[1] < m_tileDataStride - 1 &&
                m_ijk[2] < m_tileDataStride - 1))
            {
                // In this case all values can be looked up from the extended tile
                // itself. The selective lookup optimization used in the non-extended case does
                // not speed up the extended case.
                return lerp(
                    lerp(
                        lerp(value, getOffsetExtendedValue(nanovdb::Coord(0, 0, 1)),
                            uvw[2]),
                        lerp(getOffsetExtendedValue(nanovdb::Coord(0, 1, 0)),
                            getOffsetExtendedValue(nanovdb::Coord(0, 1, 1)),
                            uvw[2]),
                        uvw[1]),
                    lerp(lerp(getOffsetExtendedValue(nanovdb::Coord(1, 0, 0)),
                        getOffsetExtendedValue(nanovdb::Coord(1, 0, 1)),
                        uvw[2]),
                        lerp(getOffsetExtendedValue(nanovdb::Coord(1, 1, 0)),
                            getOffsetExtendedValue(nanovdb::Coord(1, 1, 1)),
                            uvw[2]),
                        uvw[1]),
                    uvw[0]);
            }
            else {
                // Here we optimize by only fetching the tiles we actually need.
                // The accessor does perform some internal caching as well.
                // Precompute weights to avoid further lookups.
                // Note: Optimized method for getting stencil values that only fetches the tiles actually required might result in speedups.
                const int vdw = (1 << m_vw);
                value *= (1 - uvw[0]) * (1 - uvw[1]) * (1 - uvw[2]);
                float w = (1 - uvw[0]) * (1 - uvw[1]) * uvw[2];
                value += (w > 0 ? w * getExistingValueSameDepth(baseIJK + nanovdb::Coord(0, 0, vdw)) : ValueType{}); // i, j, k + 1
                w = (1 - uvw[0]) * uvw[1] * uvw[2];
                value += (w > 0 ? w * getExistingValueSameDepth(baseIJK + nanovdb::Coord(0, vdw, vdw)) : ValueType{}); // i, j+1, k + 1
                w = (1 - uvw[0]) * (uvw[1]) * (1 - uvw[2]);
                value += (w > 0 ? w * getExistingValueSameDepth(baseIJK + nanovdb::Coord(0, vdw, 0)) : ValueType{}); // i, j+1, k
                w = (uvw[0]) * (uvw[1]) * (1 - uvw[2]);
                value += (w > 0 ? w * getExistingValueSameDepth(baseIJK + nanovdb::Coord(vdw, vdw, 0)) : ValueType{}); // i+1, j+1, k
                w = (uvw[0]) * (1 - uvw[1]) * (1 - uvw[2]);
                value += (w > 0 ? w * getExistingValueSameDepth(baseIJK + nanovdb::Coord(vdw, 0, 0)) : ValueType{}); // i+1, j, k
                w = (uvw[0]) * (1 - uvw[1]) * (uvw[2]);
                value += (w > 0 ? w * getExistingValueSameDepth(baseIJK + nanovdb::Coord(vdw, 0, vdw)) : ValueType{}); // i+1, j, k + 1
                w = (uvw[0]) * (uvw[1]) * (uvw[2]);
                value += (w > 0 ? w * getExistingValueSameDepth(baseIJK + nanovdb::Coord(vdw, vdw, vdw)) : ValueType{}); // i+1, j+1, k + 1
                return value;
            }
        }
    }

    template <typename BuildT>
    template <typename RayT>
    __hostdev__ uint32_t
    ReadAccessor<BuildT>::getDim(const CoordType& in_ijk,
                         const RayT& /*ray*/,
                         const int in_maxDepth) const {
        cache(in_ijk, in_maxDepth);
        return 1 << m_vw;
    }

    template <typename BuildT>
    __hostdev__ uint32_t
        ReadAccessor<BuildT>::getDim() const {
        return 1 << m_vw;
    }

    /// @brief Handles SDF and fog.
    template <typename BuildT>
    __hostdev__ bool ReadAccessor<BuildT>::isActiveLevelSetOrFog(const CoordType& in_ijk) const {
        // When using the HDDA + nanovtt::ZeroCrossing method it's faster to always return true if m_isLevelSet is true.
        return m_isLevelSet || (!m_isLevelSet && getValue(in_ijk, Traits::kMaxDepth) > 0);
    }

    /// @brief Return the active state of the given voxel (regardless of state
    /// or location in the tree.)
    template <typename BuildT>
    __hostdev__ inline bool ReadAccessor<BuildT>::isActive(const CoordType& in_ijk) const {
        const ValueType v = getValue(in_ijk, Traits::kMaxDepth);
        return v != 0;
    }

    template <>
    __hostdev__ inline bool ReadAccessor<float>::isActive(const CoordType& in_ijk) const {
        return isActiveLevelSetOrFog(in_ijk);
    }

    template <>
    __hostdev__ inline bool ReadAccessor<double>::isActive(const CoordType& in_ijk) const {
        return isActiveLevelSetOrFog(in_ijk);
    }

    template <>
    __hostdev__ inline bool ReadAccessor<nanovdb::Fp4>::isActive(const CoordType& in_ijk) const {
        return isActiveLevelSetOrFog(in_ijk);
    }

    template <>
    __hostdev__ inline bool ReadAccessor<nanovdb::Fp8>::isActive(const CoordType& in_ijk) const {
        return isActiveLevelSetOrFog(in_ijk);
    }

    template <>
    __hostdev__ inline bool ReadAccessor<nanovdb::Fp16>::isActive(const CoordType& in_ijk) const {
        return isActiveLevelSetOrFog(in_ijk);
    }

    template <>
    __hostdev__ inline bool ReadAccessor<nanovdb::FpN>::isActive(const CoordType& in_ijk) const {
        return isActiveLevelSetOrFog(in_ijk);
    }

    template <>
    __hostdev__ inline bool ReadAccessor<nanovdb::Vec3f>::isActive(const CoordType& in_ijk) const {
        const ValueType v = getValue(in_ijk, Traits::kMaxDepth);
        return nanovdb::Abs(v[0]) > 0 || nanovdb::Abs(v[1]) > 0 || nanovdb::Abs(v[2]) > 0;
    }

    template <>
    __hostdev__ inline bool ReadAccessor<nanovdb::Vec3d>::isActive(const CoordType& in_ijk) const {
        const ValueType v = getValue(in_ijk, Traits::kMaxDepth);
        return nanovdb::Abs(v[0]) > 0 || nanovdb::Abs(v[1]) > 0 || nanovdb::Abs(v[2]) > 0;
    }

    template <>
    __hostdev__ inline bool ReadAccessor<nanovdb::Vec4f>::isActive(const CoordType& in_ijk) const {
        const ValueType v = getValue(in_ijk, Traits::kMaxDepth);
        return nanovdb::Abs(v[0]) > 0 || nanovdb::Abs(v[1]) > 0 || nanovdb::Abs(v[2]) > 0 || nanovdb::Abs(v[3]) > 0;
    }

    template <typename BuildT>
    __hostdev__ const typename ReadAccessor<BuildT>::RootType& ReadAccessor<BuildT>::root() const {
        return *m_tree;
    }

    template <typename BuildT>
    __hostdev__ void ReadAccessor<BuildT>::cache(const CoordType& in_ijk,
                                         const int        in_maxDepth) const {
        // Note: there's a slight overhead with getValue because value is not
        // used, but
        //       avoiding the get of the value did not seem to speed things up.
        getValue(in_ijk, in_maxDepth);
    }

    /// @brief Return the voxel value at the given coordinate.
    template <typename BuildT>
    __hostdev__ typename ReadAccessor<BuildT>::ValueType ReadAccessor<BuildT>::getValue(const CoordType& in_ijk,
        const int in_maxDepth) const
    {
        // The code below has been optimized for this case
        static_assert(Traits::kTileWidth == 2, "getValue requires Traits::kTileWidth == 2");
        assert(in_maxDepth >= 0 && in_maxDepth <= Traits::kMaxDepth);

        // This version traverses bottom up from the last tile
        m_found = true;
        nanovdb::Coord diff{ in_ijk[0] - m_tijk[0], in_ijk[1] - m_tijk[1],
                            in_ijk[2] - m_tijk[2] };
        uint32_t       hash = static_cast<uint32_t>(diff[0]) |
            static_cast<uint32_t>(diff[1]) |
            static_cast<uint32_t>(diff[2]);

        // m_depth is the depth of the tile (also for coalesced tiles).
        int tvw = Traits::tileVoxelWidth(m_depth);
        m_vw = Traits::kMaxDepth - m_depth + 2;
        uint32_t mask = (0xFFFFFFFF << static_cast<uint32_t>(Traits::kMaxDepth - m_depth + 2));
        while ((hash >= static_cast<uint32_t>(tvw) || m_depth > in_maxDepth) &&
            m_depth > m_rootDepth) {
            const int tShift = m_container->parentAndShift(m_tile);
            m_depth -= tShift;
            m_vw += tShift;
            mask = (mask << static_cast<uint32_t>(tShift));
            m_tijk[0] = m_tijk[0] & static_cast<CoordValueType>(mask);
            m_tijk[1] = m_tijk[1] & static_cast<CoordValueType>(mask);
            m_tijk[2] = m_tijk[2] & static_cast<CoordValueType>(mask);
            diff = { in_ijk[0] - m_tijk[0], in_ijk[1] - m_tijk[1],
                    in_ijk[2] - m_tijk[2] };
            hash = static_cast<uint32_t>(diff[0]) |
                static_cast<uint32_t>(diff[1]) |
                static_cast<uint32_t>(diff[2]);
            tvw = (tvw << tShift);
        }

        // At this point we've found a containing tile or we're at the root depth.
        if ((hash >= static_cast<uint32_t>(tvw) || m_tile == kInvalidTileIndex) &&
            m_depth == m_rootDepth) {
            // The ijk is not inside the current root
            m_tile = m_container->findRoot(in_ijk, m_tijk);
            if (m_tile == kInvalidTileIndex) {
                m_vw = Traits::kMaxDepth - m_rootDepth;
                m_ijk = {};
                m_tijk = { Traits::kTileDomainWidth, Traits::kTileDomainWidth, Traits::kTileDomainWidth };
                m_found = false;
                return m_tree->m_defaultValue;
            }
            else {
                diff = { in_ijk[0] - m_tijk[0], in_ijk[1] - m_tijk[1],
                        in_ijk[2] - m_tijk[2] };
            }
        }

        // The performance of this if-statement is ok because there's typically no branch-divergence.
        if (m_depth < in_maxDepth) {
            m_vw = Traits::kMaxDepth - m_depth + 1;
            TileIndex childTile = kInvalidTileIndex;
            int tShift = m_container->childIndexAndShift(m_tile, diff, m_vw, m_ijk, childTile);
            bool done = (childTile == kInvalidTileIndex);
            while (!done) {
                m_tile = childTile;
                m_depth += tShift;
                m_tijk[0] += (m_ijk[0] << m_vw);
                m_tijk[1] += (m_ijk[1] << m_vw);
                m_tijk[2] += (m_ijk[2] << m_vw);
                diff[0] = in_ijk[0] - m_tijk[0];
                diff[1] = in_ijk[1] - m_tijk[1];
                diff[2] = in_ijk[2] - m_tijk[2];
                m_vw -= 1;
                tShift = m_container->childIndexAndShift(m_tile, diff, m_vw, m_ijk, childTile);
                done = (childTile == kInvalidTileIndex ||
                    m_depth >= in_maxDepth);
            }
        }

        m_vw = Traits::kMaxDepth - m_depth;
        return m_tree->getValue(m_tile, diff, m_vw, m_tileDataStride, m_ijk);
    } // getValue

    /// @brief Caches the tile at the given coordinate.
    template <typename BuildT>
    __hostdev__ void ReadAccessor<BuildT>::cacheTile(const CoordType& in_ijk,
        const int in_maxDepth) const
    {
        // The code below has been optimized for this case
        static_assert(Traits::kTileWidth == 2, "getValue requires Traits::kTileWidth == 2");
        assert(in_maxDepth >= 0 && in_maxDepth <= Traits::kMaxDepth);

        // This version traverses bottom up from the last tile
        m_found = true;
        nanovdb::Coord diff{ in_ijk[0] - m_tijk[0], in_ijk[1] - m_tijk[1],
                            in_ijk[2] - m_tijk[2] };
        uint32_t       hash = static_cast<uint32_t>(diff[0]) |
            static_cast<uint32_t>(diff[1]) |
            static_cast<uint32_t>(diff[2]);

        // m_depth is the depth of the tile (also for coalesced tiles).
        int tvw = Traits::tileVoxelWidth(m_depth);
        m_vw = Traits::kMaxDepth - m_depth + 2;
        uint32_t mask = (0xFFFFFFFF << static_cast<uint32_t>(Traits::kMaxDepth - m_depth + 2));
        while ((hash >= static_cast<uint32_t>(tvw) || m_depth > in_maxDepth) &&
            m_depth > m_rootDepth) {
            const int tShift = m_container->parentAndShift(m_tile);
            m_depth -= tShift;
            m_vw += tShift;
            mask = (mask << static_cast<uint32_t>(tShift));
            m_tijk[0] = m_tijk[0] & static_cast<CoordValueType>(mask);
            m_tijk[1] = m_tijk[1] & static_cast<CoordValueType>(mask);
            m_tijk[2] = m_tijk[2] & static_cast<CoordValueType>(mask);
            diff = { in_ijk[0] - m_tijk[0], in_ijk[1] - m_tijk[1],
                    in_ijk[2] - m_tijk[2] };
            hash = static_cast<uint32_t>(diff[0]) |
                static_cast<uint32_t>(diff[1]) |
                static_cast<uint32_t>(diff[2]);
            tvw = (tvw << tShift);
        }

        // At this point we've found a containing tile or we're at the root depth.
        if ((hash >= static_cast<uint32_t>(tvw) || m_tile == kInvalidTileIndex) &&
            m_depth == m_rootDepth) {
            // The ijk is not inside the current root
            m_tile = m_container->findRoot(in_ijk, m_tijk);
            if (m_tile == kInvalidTileIndex) {
                m_vw = Traits::kMaxDepth - m_rootDepth;
                m_ijk = {};
                m_tijk = { Traits::kTileDomainWidth, Traits::kTileDomainWidth, Traits::kTileDomainWidth };
                m_found = false;
                return;
            }
            else {
                diff = { in_ijk[0] - m_tijk[0], in_ijk[1] - m_tijk[1],
                        in_ijk[2] - m_tijk[2] };
            }
        }

        // The performance of this if-statement is ok because there's typically no branch-divergence.
        if (m_depth < in_maxDepth) {
            m_vw = Traits::kMaxDepth - m_depth + 1;
            TileIndex childTile = kInvalidTileIndex;
            int tShift = m_container->childIndexAndShift(m_tile, diff, m_vw, m_ijk, childTile);
            bool done = (childTile == kInvalidTileIndex);
            while (!done) {
                m_tile = childTile;
                m_depth += tShift;
                m_tijk[0] += (m_ijk[0] << m_vw);
                m_tijk[1] += (m_ijk[1] << m_vw);
                m_tijk[2] += (m_ijk[2] << m_vw);
                diff[0] = in_ijk[0] - m_tijk[0];
                diff[1] = in_ijk[1] - m_tijk[1];
                diff[2] = in_ijk[2] - m_tijk[2];
                m_vw -= 1;
                tShift = m_container->childIndexAndShift(m_tile, diff, m_vw, m_ijk, childTile);
                done = (childTile == kInvalidTileIndex ||
                    m_depth >= in_maxDepth);
            }
        }

        m_vw = Traits::kMaxDepth - m_depth;
        int tileDataShift;
        m_tree->tileDataStrideAndShift(m_tile, m_tileDataStride, tileDataShift);
        m_vw -= tileDataShift;
        m_ijk[0] = (m_ijk[0] >> m_vw);
        m_ijk[1] = (m_ijk[1] >> m_vw);
        m_ijk[2] = (m_ijk[2] >> m_vw);
    } // cacheTile

    /// @brief Return the voxel value at the given coordinate.
    template <typename BuildT>
    __hostdev__ typename ReadAccessor<BuildT>::ValueType ReadAccessor<BuildT>::getExistingValueSameDepth(const CoordType& in_ijk) const {
        // Assumptions:
        // - Tile containing in_ijk exists and is at same depth as the accessor currently points to.
        // - The tile is NOT at the root (because at the root all tiles are assumed to be extended tiles).

        return getValue(in_ijk, Traits::kMaxDepth - m_vw);
    }

    template<typename BuildT> using NanoVTTGrid = Grid<Tree<BuildT>>;

    using FloatGrid = NanoVTTGrid<float>;
    using DoubleGrid = NanoVTTGrid<double>;
    using Int32Grid = NanoVTTGrid<int32_t>;
    using UInt32Grid = NanoVTTGrid<uint32_t>;
    using Int64Grid = NanoVTTGrid<int64_t>;
    using Vec3fGrid = NanoVTTGrid<nanovdb::Vec3f>;
    using Vec3dGrid = NanoVTTGrid<nanovdb::Vec3d>;
    using Vec4fGrid = NanoVTTGrid<nanovdb::Vec4f>;
    using Vec4dGrid = NanoVTTGrid<nanovdb::Vec4d>;
    using Vec3IGrid = NanoVTTGrid<nanovdb::Vec3i>;
    using MaskGrid = NanoVTTGrid<nanovdb::ValueMask>;
    using BoolGrid = NanoVTTGrid<bool>;

    //////////////////////////////////////////////////

    /// @brief This is a convenient class that allows for access to grid meta-data
    ///        that are independent of the value type of a grid. That is, this class
    ///        can be used to get information about a grid without actually knowing
    ///        its ValueType.
    class GridMetaData
    {
        // We cast to a grid templated on a dummy ValueType which is safe because we are very
        // careful only to call certain methods which are known to be invariant to the ValueType!
        // In other words, don't use this technique unless you are intimately familiar with the
        // memory-layout of the data structure and the reasons why certain methods are safe
        // to call and others are not!
        using GridT = NanoVTTGrid<int>;
        __hostdev__ const GridT& grid() const { return *reinterpret_cast<const GridT*>(this); }

    public:
        __hostdev__ bool        isValid() const { return this->grid().isValid(); }
        __hostdev__ uint64_t    gridSize() const { return this->grid().gridSize(); }
        __hostdev__ uint32_t    gridIndex() const { return this->grid().gridIndex(); }
        __hostdev__ uint32_t    gridCount() const { return this->grid().gridCount(); }
        __hostdev__ const char* shortGridName() const { return this->grid().shortGridName(); }
        __hostdev__ nanovdb::GridType    gridType() const { return this->grid().gridType(); }
        __hostdev__ nanovdb::GridClass   gridClass() const { return this->grid().gridClass(); }
        __hostdev__ bool        isLevelSet() const { return this->grid().isLevelSet(); }
        __hostdev__ bool        isFogVolume() const { return this->grid().isFogVolume(); }
        __hostdev__ bool        isPointIndex() const { return this->grid().isPointIndex(); }
        __hostdev__ bool        isPointData() const { return this->grid().isPointData(); }
        __hostdev__ bool        isMask() const { return this->grid().isMask(); }
        __hostdev__ bool        isStaggered() const { return this->grid().isStaggered(); }
        __hostdev__ bool        isUnknown() const { return this->grid().isUnknown(); }
        __hostdev__ const nanovdb::Map& map() const { return this->grid().map(); }
        __hostdev__ const nanovdb::BBox<nanovdb::Vec3R>& worldBBox() const { return this->grid().worldBBox(); }
        __hostdev__ const nanovdb::BBox<nanovdb::Coord>& indexBBox() const { return this->grid().indexBBox(); }
        __hostdev__ nanovdb::Vec3R              voxelSize() const { return this->grid().voxelSize(); }
        __hostdev__ int                blindDataCount() const { return this->grid().blindDataCount(); }
        __hostdev__ const nanovdb::GridBlindMetaData& blindMetaData(int n) const { return this->grid().blindMetaData(n); }
        __hostdev__ uint64_t                 activeVoxelCount() const { return this->grid().activeVoxelCount(); }
        __hostdev__ uint32_t                 activeTileCount(uint32_t n) const { return this->grid().tree().activeTileCount(n); }
        __hostdev__ uint32_t                 nodeCount(uint32_t level) const { return this->grid().tree().nodeCount(level); }
        __hostdev__ uint64_t                 checksum() const { return this->grid().checksum(); }
        __hostdev__ bool                     isEmpty() const { return this->grid().isEmpty(); }
        __hostdev__ nanovdb::Version                  version() const { return this->grid().version(); }
    }; // GridMetaData

}
