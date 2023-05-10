// Copyright Contributors to the OpenVDB Project
// SPDX-License-Identifier: MPL-2.0

/*!
    \file Primitives.h

    \authors Ken Museth and Autodesk
*/

#ifndef NANOVTT_PRIMITIVES_H_HAS_BEEN_INCLUDED
#define NANOVTT_PRIMITIVES_H_HAS_BEEN_INCLUDED

#include <cstring>
#include <string.h>
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
#pragma clang diagnostic ignored "-Wunused-template"
#endif
#include <nanovdb/util/Primitives.h>
#if defined(__clang__)
#pragma clang diagnostic pop
#endif
#include <list>
#include <vector>
#include <cstdint>
#include <tbb/parallel_for.h>
#include <nanovdb/nanovtt/util/GridHandle.h>
#include <nanovdb/nanovtt/util/Quantization.h>
#include <nanovdb/nanovtt/util/GridChecksum.h>

namespace nanovtt {

/// <summary>
/// Create an empty grid+tree handle of a given type that can be used in
/// conjunction with the container stored in in_topology. All tile data will
/// have a extended to be used for faster sampling. Can be used to generate
/// standalone grids that can re-use an existing topology/container. The grid is
/// initialized to the default value.
///
/// For now only works if in_coalesceTileData was false in the vttToNanoVTT
/// conversion for the topology.
/// </summary>
/// <typeparam name="ValueT">The value type.</typeparam>
/// <typeparam name="BufferT">The buffer type.</typeparam>
/// <param name="in_topology">The GridHandle containing the container and at
/// least one grid.</param> <param name="in_defaultValue">The default
/// value.</param> <param name="in_name">The name of the output grid.</param>
/// <param name="in_gridClass">The grid class.</param>
/// <param name="in_padBufferForVolumetricTexture">If true, will pad the buffer
/// such that it can be used as a 3D texture of dimensions N cubed.</param>
/// <param name="in_volumetricTexureElementSize">The size of each element in the
/// 3D texture in bytes.</param> <param name="in_enable3DTextureTileLookup">If
/// enabled, the buffer allocated will contain a 3D texture of tile data and no
/// additional information.</param> <param name="in_tileDataWidth">The tile data
/// width of the tiles stored in the 3D texture. Only relevant when
/// in_enable3DTextureTileLookup is true.</param> <param
/// name="in_discardFinestLevels">The number of finest levels that are
/// discarded. Can be used to reduce memory requirements and increase
/// performance, at the cost of quality.</param> <param name="io_buffer">The
/// buffer used to construct the data.</param>
template<typename BufferT = nanovdb::HostBuffer, typename ValueT = float>
void createEmptyExtendedTileGrid(
    const GridHandle<BufferT>& in_topology,
    ValueT                       in_defaultValue,
    const std::string& in_name = "",
    const nanovdb::GridClass     in_gridClass = nanovdb::GridClass::Unknown,
    const bool                   in_padBufferForVolumetricTexture = false,
    const int                    in_volumetricTexureElementSize = 4,
    const bool                   in_enable3DTextureTileLookup = false,
    const int                    in_tileDataWidth = 0,
    const int                    in_discardFinestLevels = 0,
    BufferT&                     io_buffer    = BufferT())
{
    assert(in_discardFinestLevels >= 0);
    using BuildType = ValueT;
    constexpr uint64_t valueTypeSize = sizeof(typename Tree<BuildType>::ValueType);
    const Container* container = in_topology.container();
    if (container == nullptr)
        return;

    const uint64_t numNanoVTTTiles = container->activeTileCount();

    if (in_enable3DTextureTileLookup) {
        // In this case we only store the tile data itself
        const uint64_t dim = static_cast<uint64_t>(
            // Note: we add 1 to account for a tile of default values
            std::ceil(std::cbrt(static_cast<double>(numNanoVTTTiles + 1))) * in_tileDataWidth);
        const uint64_t bufferSize = [&]() {
            if (in_discardFinestLevels == 0) {
                return dim * dim * dim * in_volumetricTexureElementSize;
            }
            else
            {
                uint64_t numNanoVTTTiles2 = 0;
                for (int d = container->rootDepth(); d <= Traits::kMaxDepth - in_discardFinestLevels; ++d) {
                    numNanoVTTTiles2 += container->activeTileCount(d);
                }
                assert(dim % in_tileDataWidth == 0);
                const uint64_t dimTD = dim / in_tileDataWidth;
                const uint64_t dimZ = static_cast<uint64_t>(std::ceil(static_cast<double>(numNanoVTTTiles2 + 1) / (dimTD * dimTD)) * in_tileDataWidth);
                return dim * dim * dimZ * in_volumetricTexureElementSize;
            }
        }();

        io_buffer = BufferT::create(bufferSize);
        tbb::parallel_for(
            tbb::blocked_range<uint64_t>(0, bufferSize / in_volumetricTexureElementSize),
            [&](tbb::blocked_range<uint64_t>& in_r) {
                ValueT* td = reinterpret_cast<ValueT*>(io_buffer.data());
                for (uint64_t t = in_r.begin(); t != in_r.end(); ++t) {
                    td[t] = in_defaultValue;
                }
            });
    }
    else {
        assert(in_discardFinestLevels == 0);
        const uint64_t tdOffsetConstant = (numNanoVTTTiles + 1) * sizeof(uint64_t);
        const uint64_t treeSize = sizeof(Tree<ValueT>) + Traits::kExtendedTileDataSize * numNanoVTTTiles * valueTypeSize + tdOffsetConstant;
        uint64_t bufferSize = sizeof(Grid<Tree<ValueT>>) + treeSize;
        bufferSize += NANOVDB_DATA_ALIGNMENT - (bufferSize % NANOVDB_DATA_ALIGNMENT);
        if (in_padBufferForVolumetricTexture) {
            unsigned int dim =
                static_cast<unsigned int>(std::ceil(std::cbrt(static_cast<double>(bufferSize) / in_volumetricTexureElementSize)));
            bufferSize = dim * dim * dim * in_volumetricTexureElementSize;
        }
        io_buffer = BufferT::create(bufferSize);

        if (io_buffer.data() == nullptr)
            return;

        // Make sure buffer contains no uninitialized memory which can otherwise cause uninitialized
        // memory in data required for alignment.
        std::memset(io_buffer.data(), 0, bufferSize);

        // Grid and Tree
        Grid<Tree<BuildType>>* gridVTT =
            reinterpret_cast<Grid<Tree<BuildType>>*>(io_buffer.data());
        Tree<BuildType>* treeVTT = &gridVTT->tree();

        // Grid
        const GridMetaData* gridMetaData = in_topology.gridMetaData();
        const double voxelSize = gridMetaData->voxelSize()[0];
        nanovdb::GridData* gridDataVTT = reinterpret_cast<nanovdb::GridData*>(gridVTT);
        gridDataVTT->setFlagsOff();
        gridDataVTT->mMagic = NANOVTT_MAGIC_NUMBER;
        gridDataVTT->mVersion =
            nanovdb::Version(NANOVTT_MAJOR_VERSION_NUMBER,
                NANOVTT_MINOR_VERSION_NUMBER,
                NANOVTT_PATCH_VERSION_NUMBER);
        gridDataVTT->mGridClass = in_gridClass;
        gridDataVTT->mGridIndex = 0;
        gridDataVTT->mGridCount = 1;
        gridDataVTT->mGridSize = bufferSize;
#if defined(_MSC_VER)
        strncpy_s(gridDataVTT->mGridName, nanovdb::GridData::MaxNameSize - 1,
                  in_name.c_str(), in_name.size());
#else
        // Chorus: strncpy_s not available on other platforms.
        // Since mGridName is allocated as a fixed size array, the copy is ok.
        assert(gridDataVTT != nullptr);
        strncpy(gridDataVTT->mGridName, in_name.c_str(), nanovdb::GridData::MaxNameSize - 1); /* Flawfinder: ignore */
#endif
        gridDataVTT->mGridName
            [nanovdb::GridData::MaxNameSize - 1] = 0;
        gridDataVTT->mGridType = nanovdb::mapToGridType<BuildType>();
        gridDataVTT->mVoxelSize[0] = voxelSize;
        gridDataVTT->mVoxelSize[1] = voxelSize;
        gridDataVTT->mVoxelSize[2] = voxelSize;
        gridDataVTT->setBreadthFirstOn();

        gridDataVTT->mBlindMetadataOffset = 0;
        gridDataVTT->mBlindMetadataCount = 0;
        gridDataVTT->setLongGridNameOn(false);

        const double Sx = static_cast<double>(voxelSize), Sy = static_cast<double>(voxelSize), Sz = static_cast<double>(voxelSize);
        constexpr double Tx = 0, Ty = 0, Tz = 0;
        const double SxInv = 1 / Sx, SyInv = 1 / Sy, SzInv = 1 / Sz;
        const double mat[4][4] = {
            {Sx, 0.0, 0.0, 0.0}, // row 0
            {0.0, Sy, 0.0, 0.0}, // row 1
            {0.0, 0.0, Sz, 0.0}, // row 2
            {Tx, Ty, Tz, 1.0},          // row 3
        };
        const double invMat[4][4] = {
            {SxInv , 0.0, 0.0, 0.0}, // row 0
            {0.0, SyInv, 0.0, 0.0}, // row 1
            {0.0, 0.0, SzInv, 0.0}, // row 2
            {-Tx, -Ty, -Tz, 1.0},           // row 3
        };
        gridDataVTT->mMap.set(mat, invMat, 1.0);

        // By default, the grid is empty
        gridDataVTT->mWorldBBox[0][0] = std::numeric_limits<double>::max();
        gridDataVTT->mWorldBBox[0][1] = std::numeric_limits<double>::max();
        gridDataVTT->mWorldBBox[0][2] = std::numeric_limits<double>::max();
        gridDataVTT->mWorldBBox[1][0] = -std::numeric_limits<double>::max();
        gridDataVTT->mWorldBBox[1][1] = -std::numeric_limits<double>::max();
        gridDataVTT->mWorldBBox[1][2] = -std::numeric_limits<double>::max();

        gridDataVTT->mChecksum = 0;

        // Tree

        // Tree
        // Note: for a standalone grid, the container offset is not valid
        treeVTT->m_containerOffset = std::numeric_limits<uint64_t>::max();
        treeVTT->m_treeSize = treeSize;
        treeVTT->m_defaultValue = in_defaultValue;

        // By default, the tree is empty
        treeVTT->m_indexBBox[0][0] = std::numeric_limits<int32_t>::max();
        treeVTT->m_indexBBox[0][1] = std::numeric_limits<int32_t>::max();
        treeVTT->m_indexBBox[0][2] = std::numeric_limits<int32_t>::max();
        treeVTT->m_indexBBox[1][0] = -std::numeric_limits<int32_t>::max();
        treeVTT->m_indexBBox[1][1] = -std::numeric_limits<int32_t>::max();
        treeVTT->m_indexBBox[1][2] = -std::numeric_limits<int32_t>::max();

        tbb::parallel_for(
            tbb::blocked_range<uint64_t>(0, numNanoVTTTiles),
            [&](tbb::blocked_range<uint64_t>& in_r) {
                for (uint64_t t = in_r.begin(); t != in_r.end(); ++t) {
                    treeVTT->setTileDataOffset(static_cast<TileIndex>(t), t * Tree<BuildType>::kExtendedTileDataBytes + tdOffsetConstant);
                    treeVTT->setTileDataStride(static_cast<TileIndex>(t), static_cast<uint8_t>(Traits::kExtendedTileDataWidth));
                    treeVTT->setTileDataShift(static_cast<TileIndex>(t), 0);
                    ValueT* td = treeVTT->tileData(static_cast<TileIndex>(t));
                    for (int i = 0; i < Traits::kExtendedTileDataSize; ++i) {
                        td[i] = in_defaultValue;
                    }
                }
            });

        // Set entry 'numNanoVTTTiles' to facilitate
        // computation of sizes of tiles.
        treeVTT->setTileDataOffset(
            static_cast<TileIndex>(numNanoVTTTiles),
            treeSize);
    }
}



/// @brief Returns a handle to a narrow-band level set of a sphere
///
/// @param radius    Radius of sphere in world units
/// @param center    Center of sphere in world units
/// @param voxelSize Size of a voxel in world units
/// @param origin    Origin of grid in world units
/// @param name      Name of the grid
/// @param sMode     Mode of computation for the statistics.
/// @param cMode     Mode of computation for the checksum.
/// @param tolerance Global error tolerance use when VoxelT = FpN
/// @param ditherOn  If true dithering will be applied when VoxelT = {Fp4,Fp8,Fp16,FpN}
/// @param buffer    Buffer used for memory allocation by the handle
///
/// @details The @c ValueT template parameter must be float (default) or double.
///          The @c VoxelT template parameter must be one of the following:
///          float (default), double, Fp4, Fp8, Fp16 or FpN. The @c tolerance
///          argument is only used when VoxelT is set to FpN.
template<typename ValueT = float,
         typename VoxelT = ValueT,
         typename BufferT = nanovdb::HostBuffer>
GridHandle<BufferT> createLevelSetSphere(
    ValueT                       radius    = 100,
    const nanovdb::Vec3<ValueT>& center    = {},
    double                       voxelSize = 1.0,
	const nanovdb::Vec3<ValueT>& origin = {},
    const std::string&           name      = "sphere_ls",
    nanovdb::StatsMode           sMode     = nanovdb::StatsMode::Default,
    nanovdb::ChecksumMode        cMode     = nanovdb::ChecksumMode::Default,
    float                        tolerance = -1.0f,
    bool                         ditherOn  = false,
    const BufferT&               buffer    = BufferT()) {
    static_assert(nanovdb::is_floating_point<ValueT>::value, "createLevelSetSphere: expect floating point");
    static_assert(nanovdb::is_floating_point<typename nanovdb::BuildToValueMap<VoxelT>::Type>::value, "createLevelSetSphere: expect floating point");
    if (!(radius > 0))
        throw std::runtime_error("Sphere: radius must be positive!");
    if (!(voxelSize > 0))
        throw std::runtime_error("Sphere: voxelSize must be positive!");

    constexpr bool verbose = false;

	if (verbose) {
		std::cout << "sMode = " << static_cast<int>(sMode) << std::endl;
	}

    // Define radius of sphere with narrow-band in voxel units
    const ValueT r0 = radius / static_cast<ValueT>(voxelSize), rmax = r0 + 2;

    // Define center of sphere in voxel units
    const nanovdb::Vec3<ValueT> c(ValueT(center[0] - origin[0]) / static_cast<ValueT>(voxelSize),
        ValueT(center[1] - origin[1]) / static_cast<ValueT>(voxelSize),
        ValueT(center[2] - origin[2]) / static_cast<ValueT>(voxelSize));

    // Define bounds of the voxel coordinates
    const int imin = nanovdb::Floor(c[0] - rmax), imax = nanovdb::Ceil(c[0] + rmax);
    const int jmin = nanovdb::Floor(c[1] - rmax), jmax = nanovdb::Ceil(c[1] + rmax);
    const int kmin = nanovdb::Floor(c[2] - rmax), kmax = nanovdb::Ceil(c[2] + rmax);

    int64_t tileCountAccum[Traits::kMaxDepth + 2] = {};
    std::vector<Tile *> tiles;
    std::vector<float> tileDatas;
    uint64_t numNanoVTTTilesAtCoarserLevels = 0;
    uint64_t numLeafTiles = 0;
    std::vector<uint64_t> tOffset;
    tOffset.push_back(0);
    std::list<std::pair<TileIndex, VoxelCoord>> queue;
    queue.push_back({ kInvalidTileIndex, VoxelCoord(nanovdb::Coord{ Traits::kMinTileCoord, Traits::kMinTileCoord, Traits::kMinTileCoord }, 0) });
    tileCountAccum[1] += 1;
    const ValueT defaultValue = static_cast<ValueT>(Traits::tileVoxelWidth(0) * static_cast<ValueT>(voxelSize) * static_cast<ValueT>(1.732051));

    // Breadth first serial grid construction using a distance function. In this case a sphere.
    while (!queue.empty()) {
        const auto qf = queue.front();
        // The index of the parent tile.
        const TileIndex parentIndex = qf.first;
        const VoxelCoord& vc = qf.second;
        // The index of the current tile, and the parent of the child.
        const TileIndex index = static_cast<TileIndex>(tiles.size());
        queue.pop_front();

        Tile* tile = new Tile();
        tiles.push_back(tile);
        tile->m_tijk = vc.m_ijk;
        tile->m_parent = parentIndex;
        tile->m_tile = index;
        tile->m_depth = vc.m_depth;

        // First determine if we should refine any children further.
        // We sample the distance function at the center and check if its less than the refinement threshold.
        const int vw = Traits::voxelWidth(vc.m_depth);
        constexpr ValueT kSqrt3Half = static_cast<ValueT>(0.8660254037844386);
        bool isLeaf = true;
        // The 1.1 is to make the estimate conservative
        const ValueT refinementThreshold = static_cast<ValueT>(vw * kSqrt3Half * static_cast<ValueT>(1.1));
        if (vc.m_depth < Traits::kMaxDepth) {
            for (int k = 0, ci = 0; k < Traits::kTileWidth; ++k) for (int j = 0; j < Traits::kTileWidth; ++j) for (int i = 0; i < Traits::kTileWidth; ++i, ++ci) {
                const auto x2 = nanovdb::Pow2(static_cast<ValueT>((i + 0.5) * vw + vc.m_ijk[0]) - c[0]);
                const auto y2 = nanovdb::Pow2(static_cast<ValueT>((j + 0.5) * vw + vc.m_ijk[1]) - c[1]);
                const auto z2 = nanovdb::Pow2(static_cast<ValueT>((k + 0.5) * vw + vc.m_ijk[2]) - c[2]);
                const auto v = static_cast<ValueT>(nanovdb::Sqrt(x2 + y2 + z2)) - r0;
                if (std::abs(v) < refinementThreshold) {
                    // This child should be refined
                    isLeaf = false;
                    queue.push_back({ index,
                        VoxelCoord(nanovdb::Coord{vc.m_ijk[0] + i * vw, vc.m_ijk[1] + j * vw, vc.m_ijk[2] + k * vw}, vc.m_depth + 1) });
                    // Set child index in the parent
                    tile->m_children[ci] = static_cast<TileIndex>(tiles.size() + queue.size() - 1);
                    if (tileCountAccum[vc.m_depth + 2] == 0) {
                        // This is the first node at this depth, so we propagation the accumulation to this depth
                        tileCountAccum[vc.m_depth + 2] = tileCountAccum[vc.m_depth + 1];
                    }
                    tileCountAccum[vc.m_depth + 2] += 1;
                }
                else {
                    // The child is not refined
                    tile->m_children[ci] = kInvalidTileIndex;
                }
            }
        }
        else {
            for (int k = 0, ci = 0; k < Traits::kTileWidth; ++k) for (int j = 0; j < Traits::kTileWidth; ++j) for (int i = 0; i < Traits::kTileWidth; ++i, ++ci) {
                tile->m_children[ci] = kInvalidTileIndex;
            }
        }
        if (vc.m_depth < Traits::kMaxDepth) {
            ++numNanoVTTTilesAtCoarserLevels;
            if (isLeaf) {
                tOffset.push_back(0);
                ++numLeafTiles;
            }
            else {
                tOffset.push_back(sizeof(TileIndex)* Traits::kTileSize);
            }
        }

        // Next voxelize values of this tile
        const int vdw = Traits::voxelDataWidth(vc.m_depth);
        for (int k = 0; k < Traits::kTileDataWidth; ++k) for (int j = 0; j < Traits::kTileDataWidth; ++j) for (int i = 0; i < Traits::kTileDataWidth; ++i) {
            const auto x2 = nanovdb::Pow2(static_cast<ValueT>(i * vdw + vc.m_ijk[0]) - c[0]);
            const auto y2 = nanovdb::Pow2(static_cast<ValueT>(j * vdw + vc.m_ijk[1]) - c[1]);
            const auto z2 = nanovdb::Pow2(static_cast<ValueT>(k * vdw + vc.m_ijk[2]) - c[2]);
            const auto v = static_cast<ValueT>((static_cast<ValueT>(nanovdb::Sqrt(x2 + y2 + z2)) - r0) * static_cast<ValueT>(voxelSize));
            tileDatas.push_back(v);
        }
    }
    const uint64_t numNanoVTTTiles = (tileDatas.size() / Traits::kTileDataSize);
    const uint64_t tOffsetConstant = (numNanoVTTTilesAtCoarserLevels + 1) * sizeof(uint64_t) + numNanoVTTTiles * sizeof(TileIndex);
    tOffset[0] = tOffsetConstant;
    // Note: here we could use a parallel prefix sum
    for (uint64_t i = 1; i <= numNanoVTTTilesAtCoarserLevels; ++i) {
        tOffset[i] += tOffset[i - 1];
    }
    const uint64_t invalidTileOffset = (numNanoVTTTilesAtCoarserLevels - numLeafTiles) * sizeof(TileIndex) * Traits::kTileSize + tOffsetConstant;
    for (uint64_t i = 0; i < numNanoVTTTilesAtCoarserLevels; ++i) {
        if (tOffset[i] == tOffset[i + 1])
            tOffset[i] = invalidTileOffset;
    }
    if (verbose) {
        std::cout << "numNanoVTTTiles = " << numNanoVTTTiles << std::endl;
        std::cout << "numNanoVTTTilesAtCoarserLevels = " << numNanoVTTTilesAtCoarserLevels << std::endl;
        std::cout << "numLeafTiles = " << numLeafTiles << std::endl;
        std::cout << "invalidTileOffset = " << invalidTileOffset << std::endl;
    }
    // For FpN we need to call the compression routine here to know the size of the
    // resulting data.
    uint64_t treeSize = 0;
    std::vector<Codec> codecs;
    std::vector<uint64_t> tdOffset;
    uint64_t tdOffsetConstant = (numNanoVTTTiles + 1) * sizeof(uint64_t);
    Codec codecDummy{};
	nanovdb::AbsDiff oracle(tolerance);
    nanovdb::DitherLUT lut(ditherOn);
    if (std::is_same<VoxelT, nanovdb::FpN>::value) {
        treeSize = numNanoVTTTiles * sizeof(uint64_t);
        codecs.resize(numNanoVTTTiles);
        tdOffset.resize(numNanoVTTTiles);
        initOracle<VoxelT>(nanovdb::GridClass::LevelSet,
            voxelSize,
            oracle);
        for (int d = 0; d <= Traits::kMaxDepth; ++d) {
            const int64_t accum = tileCountAccum[d];
            const int64_t size = tileCountAccum[d+1] - accum;
            for (int64_t i = 0; i < size; ++i) {
                const int64_t index = accum + i;
                compression<VoxelT>(
                    &tileDatas[index * Traits::kTileDataSize],
                    lut, oracle, codecs[index]);
                treeSize += codecs[index].size + Tree<VoxelT>::kTileDataBytes;
            }
        }
        // Note: could be optimized by a parallel prefix sum
        tdOffset[0] = numNanoVTTTiles * sizeof(uint64_t);
        for (uint64_t i = 1; i < numNanoVTTTiles; ++i) {
            tdOffset[i] = codecs[i - 1].size + Tree<VoxelT>::kTileDataBytes + tdOffset[i - 1];
        }
    }
    else if (std::is_same<VoxelT, void>::value) {
        treeSize = Tree<VoxelT>::kTileDataBytes * numNanoVTTTiles + tdOffsetConstant;
    }
    else {
        // Fp4, Fp8, Fp16
        const uint64_t numExtendedTiles = 0;
        const uint64_t numNonExtendedTiles = numNanoVTTTiles;
        treeSize =
            numExtendedTiles *
            Tree<VoxelT>::kExtendedTileDataBytes +
            numNonExtendedTiles *
            Tree<VoxelT>::kTileDataBytes +
            tdOffsetConstant;
    }

    uint64_t containerSize =
        sizeof(Container) + numNanoVTTTiles * sizeof(nanovtt::Tile) +
        tOffsetConstant +
        (numNanoVTTTilesAtCoarserLevels + 1 - numLeafTiles) *
            sizeof(TileIndex) * Traits::kTileSize;
    // Make sure following data is aligned
    containerSize += NANOVDB_DATA_ALIGNMENT - (containerSize % NANOVDB_DATA_ALIGNMENT);
    uint64_t gridAndTreeSize = sizeof(Grid<Tree<VoxelT>>) +
        sizeof(Tree<VoxelT>) + treeSize;
    gridAndTreeSize += NANOVDB_DATA_ALIGNMENT - (gridAndTreeSize % NANOVDB_DATA_ALIGNMENT);
    uint64_t bufferSize = containerSize + gridAndTreeSize;
    bufferSize += NANOVDB_DATA_ALIGNMENT - (bufferSize % NANOVDB_DATA_ALIGNMENT);
    GridHandle<BufferT> handle(BufferT::create(bufferSize, &buffer));

    if (handle.data() == nullptr)
        return handle;

    // Make sure buffer contains no uninitialized memory which can otherwise cause uninitialized
    // memory in data required for alignment.
    std::memset(handle.data(), 0, bufferSize);

    // Build data structure:
    // Container
    {
        Container* container = reinterpret_cast<Container*>(handle.data());
        container->m_bufferSize = bufferSize;
        container->m_containerSize = containerSize;
        container->m_parentOffset = (numNanoVTTTilesAtCoarserLevels + 1) * sizeof(uint64_t);
        container->m_rootDepth = 0;
        for (uint32_t i = 0; i < Container::kTableSize; ++i) {
            container->m_table[i] =
                kInvalidTileIndex;
        }
        container->m_tableMin = { Traits::kMinTileCoord, Traits::kMinTileCoord, Traits::kMinTileCoord };
        container->addRoot({ 0, 0, 0 }, 0);
        for (int d = 0; d <= Traits::kMaxDepth + 1; ++d) {
            container->m_tileCountAccum[d] = static_cast<TileIndex>(tileCountAccum[d]);
        }
        for (uint64_t i = 0; i < tiles.size(); ++i) {
            const Tile& tile = *tiles[i];
            container->setParent(static_cast<TileIndex>(i), tile.parent());

            if (!tile.atMaxDepth()) {
                container->setTileOffset(static_cast<TileIndex>(i), tOffset[i]);
#ifdef NANOVTT_SUPPORT_COALESCED_TILES
                container->setTileShift(static_cast<TileIndex>(i), 0);
#endif
                if (!tile.isLeaf()) {
                    TileIndex* tileVTT = container->children(static_cast<TileIndex>(i));
                    for (int i2 = 0; i2 < Traits::kTileSize; ++i2) {
                        if (tile.childIndex(i2) != kInvalidTileIndex) {
                            tileVTT[i2] = tile.childIndex(i2);
                        } else {
                            tileVTT[i2] = kInvalidTileIndex;
                        }
                    }
                }
            }
            delete tiles[i];
        }

        {
            container->setTileOffset(static_cast<TileIndex>(tOffset.size() - 1), tOffset.back());
#ifdef NANOVTT_SUPPORT_COALESCED_TILES
            container->setTileShift(static_cast<TileIndex>(tOffset.size() - 1), 0);
#endif
            TileIndex* tileVTT = container->children(static_cast<TileIndex>(tOffset.size() - 1));
            if (verbose) {
                std::cout << "tOffset.size" << tOffset.size() << std::endl;
                std::cout << "invalid tile = " << tileVTT << std::endl;
            }
            for (int i = 0; i < Traits::kTileSize; ++i) {
                tileVTT[i] = nanovtt::kInvalidTileIndex;
            }
        }

        // Grid
        Grid<Tree<VoxelT>>* grid = reinterpret_cast<Grid<Tree<VoxelT>>*>(reinterpret_cast<uint8_t*>(handle.data()) + containerSize);

        // Grid Data
        nanovdb::GridData* gridData = reinterpret_cast<nanovdb::GridData*>(grid);
        gridData->setFlagsOff();
        gridData->mMagic = NANOVTT_MAGIC_NUMBER;
        gridData->mVersion =
            nanovdb::Version(NANOVTT_MAJOR_VERSION_NUMBER,
                NANOVTT_MINOR_VERSION_NUMBER,
                NANOVTT_PATCH_VERSION_NUMBER);
        gridData->mGridClass = nanovdb::GridClass::LevelSet;
        gridData->mGridIndex = 0;
        gridData->mGridCount = 1;
        gridData->mGridSize = gridAndTreeSize;
#if defined(_MSC_VER)
        strncpy_s(gridData->mGridName, nanovdb::GridData::MaxNameSize - 1,
            name.c_str(), name.size());
#else
        // Chorus: strncpy_s not available on other platforms.
        // Since mGridName is allocated as a fixed size array, the copy is ok.
        assert(gridData != nullptr);
        strncpy(gridData->mGridName, name.c_str(), nanovdb::GridData::MaxNameSize - 1); /* Flawfinder: ignore */
#endif

        gridData->mGridName[nanovdb::GridData::MaxNameSize - 1] = 0;
        gridData->mGridType = nanovdb::mapToGridType<VoxelT>();
        gridData->mWorldBBox[0][0] = imin * voxelSize;
        gridData->mWorldBBox[0][1] = jmin * voxelSize;
        gridData->mWorldBBox[0][2] = kmin * voxelSize;
        gridData->mWorldBBox[1][0] = imax * voxelSize;
        gridData->mWorldBBox[1][1] = jmax * voxelSize;
        gridData->mWorldBBox[1][2] = kmax * voxelSize;
        gridData->mVoxelSize[0] = voxelSize;
        gridData->mVoxelSize[1] = voxelSize;
        gridData->mVoxelSize[2] = voxelSize;
        gridData->setBreadthFirstOn();
        // Map
        const double Tx = 0, Ty = 0, Tz = 0;
        const double mat[4][4] = {
            {voxelSize, 0.0, 0.0, 0.0}, // row 0
            {0.0, voxelSize, 0.0, 0.0}, // row 1
            {0.0, 0.0, voxelSize, 0.0}, // row 2
            {Tx, Ty, Tz, 1.0}, // row 3
        };
        const double invMat[4][4] = {
            {1 / voxelSize, 0.0, 0.0, 0.0}, // row 0
            {0.0, 1 / voxelSize, 0.0, 0.0}, // row 1
            {0.0, 0.0, 1 / voxelSize, 0.0}, // row 2
            {-Tx, -Ty, -Tz, 1.0}, // row 3
        };
        gridData->mMap.set(mat, invMat, 1.0);
        gridData->mBlindMetadataOffset = 0;
        gridData->mBlindMetadataCount = 0;
        gridData->setLongGridNameOn(false);

        // Tree
        Tree<VoxelT>& tree = grid->tree();
        tree.m_containerOffset = reinterpret_cast<uint8_t*>(&tree) -
            reinterpret_cast<uint8_t*>(container);
        tree.m_treeSize = sizeof(Tree<VoxelT>) + treeSize;
        tree.m_defaultValue = defaultValue;
        tree.m_indexBBox[0][0] = imin;
        tree.m_indexBBox[0][1] = jmin;
        tree.m_indexBBox[0][2] = kmin;
        tree.m_indexBBox[1][0] = imax;
        tree.m_indexBBox[1][1] = jmax;
        tree.m_indexBBox[1][2] = kmax;
        tbb::parallel_for(
            tbb::blocked_range<uint64_t>(0, numNanoVTTTiles),
            [&](tbb::blocked_range<uint64_t>& in_r) {
                for (uint64_t i = in_r.begin(); i != in_r.end(); ++i) {
                    const float* tdData = &tileDatas[i * Traits::kTileDataSize];
                    if (std::is_same<VoxelT, nanovdb::FpN>::value) {
                        // For FpN we need to add data per tile using the
                        // per-tile codec.
                        addDataToTile(tdData, lut, static_cast<TileIndex>(i), codecs[i], tdOffset[i],
                                      tree);
                    } else if (std::is_same<VoxelT, void>::value) {
                        addDataToTile(
                            tdData, lut, static_cast<TileIndex>(i), codecDummy,
                            tdOffsetConstant + i * Tree<VoxelT>::kTileDataBytes,
                            tree);
                    } else {
                        // Fp4, Fp8, Fp16
                        addDataToTile(
                            tdData, lut, static_cast<TileIndex>(i), codecDummy,
                            tdOffsetConstant + i * Tree<VoxelT>::kTileDataBytes,
                            tree);
                    }
                }
            });

        if (!std::is_same<VoxelT, nanovdb::FpN>::value) {
            // Set entry 'numNanoVTTTiles' to facilitate
            // computation of sizes of tiles.
            tree.setTileDataOffset(
                static_cast<TileIndex>(numNanoVTTTiles), tdOffsetConstant + numNanoVTTTiles * Tree<VoxelT>::kTileDataBytes);
        }

        // We compute the checksum as the last step as all data may be required here.
        gridData->mChecksum = checksum(*grid, cMode);
    }

    return handle;
}

} // namespace nanovtt

#endif // NANOVTT_PRIMITIVES_H_HAS_BEEN_INCLUDED
