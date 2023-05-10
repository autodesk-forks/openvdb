// Copyright Contributors to the OpenVDB Project
// SPDX-License-Identifier: MPL-2.0

/*!
    \file Quantization.h

    \authors Ken Museth and Autodesk
*/

#ifndef NANOVTT_QUANTIZATION_H_HAS_BEEN_INCLUDED
#define NANOVTT_QUANTIZATION_H_HAS_BEEN_INCLUDED

#include <type_traits>
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
#pragma clang diagnostic ignored "-Wdeprecated-copy"
#if defined(__clang_major__) && __clang_major__ >= 14
#pragma clang diagnostic ignored "-Wnull-pointer-subtraction"
#endif
#endif
#include <nanovdb/util/DitherLUT.h>
#include <nanovdb/util/GridBuilder.h>
#if defined(__clang__)
#pragma clang diagnostic pop
#endif
#include <nanovdb/nanovtt/NanoVTT.h>

namespace nanovtt
{

struct Codec { float min, max; uint16_t log2, size; };// used for adaptive bit-rate quantization

template <typename BuildT>
typename std::enable_if<!std::is_same<nanovdb::Fp4, BuildT>::value &&
    !std::is_same<nanovdb::Fp8, BuildT>::value &&
    !std::is_same<nanovdb::Fp16, BuildT>::value &&
    !std::is_same<nanovdb::FpN, BuildT>::value>::type
    initAddDataToTile(const TileIndex in_i,
        const uint64_t in_offset,
        const int in_tileDataStride,
        const int in_tileDataShift,
        Tree<BuildT>& io_tree)
{
    io_tree.setTileDataOffset(in_i, in_offset);
    io_tree.setTileDataStride(in_i, static_cast<uint8_t>(in_tileDataStride));
    io_tree.setTileDataShift(in_i, static_cast<uint8_t>(in_tileDataShift));
}

template <typename BuildT, typename SrcT>
typename std::enable_if<!std::is_same<nanovdb::Fp4, BuildT>::value &&
    !std::is_same<nanovdb::Fp8, BuildT>::value &&
    !std::is_same<nanovdb::Fp16, BuildT>::value &&
    !std::is_same<nanovdb::FpN, BuildT>::value>::type
    addDataToTile(const SrcT* in_td,
        nanovdb::DitherLUT& /*io_lut*/,
        const TileIndex in_i,
        const Codec& /*in_codec*/,
        const int in_tdStrideSrc,
        const int in_tdStrideDst,
        const nanovdb::Coord& in_ijk,
        const nanovdb::Coord& in_tdw,
        Tree<BuildT>& io_tree)
{
    const int tdStrideSrc2 = in_tdStrideSrc * in_tdStrideSrc;
    const int tdStrideDst2 = in_tdStrideDst * in_tdStrideDst;
    using ValueType = typename Tree<BuildT>::ValueType;
    ValueType* io_td = io_tree.tileData(in_i);
    for (int k = in_ijk[2], kk = 0; k < in_ijk[2] + in_tdw[2]; ++k, ++kk)
        for (int j = in_ijk[1], jj = 0; j < in_ijk[1] + in_tdw[1]; ++j, ++jj)
            for (int i = in_ijk[0], ii = 0; i < in_ijk[0] + in_tdw[0]; ++i, ++ii)
    {
        io_td[i + j * in_tdStrideDst + k * tdStrideDst2] = *reinterpret_cast<ValueType const*>(&in_td[ii + jj * in_tdStrideSrc + kk * tdStrideSrc2]);
    }
}

template <typename BuildT, typename SrcT>
typename std::enable_if<!std::is_same<nanovdb::Fp4, BuildT>::value &&
                        !std::is_same<nanovdb::Fp8, BuildT>::value &&
                        !std::is_same<nanovdb::Fp16, BuildT>::value &&
                        !std::is_same<nanovdb::FpN, BuildT>::value>::type
addDataToTile(const SrcT* in_td,
              nanovdb::DitherLUT& /*io_lut*/,
              const TileIndex in_i,
              const Codec& /*in_codec*/,
              uint64_t in_offset,
              Tree<BuildT>& io_tree,
              const int     in_size = Traits::kTileDataSize) {
    io_tree.setTileDataOffset(in_i, in_offset);
    using ValueType = typename Tree<BuildT>::ValueType;
    ValueType* io_td = io_tree.tileData(in_i);
    if (in_size == Traits::kTileDataSize)
    {
        io_tree.setTileDataStride(in_i, static_cast<uint8_t>(Traits::kTileDataWidth));
        io_tree.setTileDataShift(in_i, 0);
        // Non-extended
        for (int i = 0; i < in_size; ++i) {
            io_td[i] = *reinterpret_cast<ValueType const*>(&in_td[i]);
        }
    }
    else {
        io_tree.setTileDataStride(in_i, static_cast<uint8_t>(Traits::kExtendedTileDataWidth));
        io_tree.setTileDataShift(in_i, 0);
        // Extended
        for (int i = 0; i < Traits::kExtendedTileDataSize; ++i) {
            io_td[i] = *reinterpret_cast<ValueType const*>(&in_td[i]);
        }
    }
}

template <typename BuildT>
typename std::enable_if<std::is_same<nanovdb::FpN, BuildT>::value>::type
initAddDataToTile(const TileIndex /*in_i*/,
    const uint64_t /*in_offset*/ ,
    const int /*in_tileDataStride*/ ,
    const int /*in_tileDataShift*/ ,
    Tree<BuildT>& /*io_tree*/ )
{
    // \todo BIFROST-TBD - We may support coalesced tile data for quantized values in the future.
}

template <typename BuildT>
typename std::enable_if<std::is_same<nanovdb::FpN, BuildT>::value>::type
addDataToTile(const float* /*in_td*/,
    nanovdb::DitherLUT& /*io_lut*/,
    const TileIndex /*in_i*/,
    const Codec& /*in_codec*/,
    const int /*in_tdStrideSrc*/,
    const int /*in_tdStrideDst*/,
    const nanovdb::Coord& /*in_ijk*/,
    const nanovdb::Coord& /*in_tdw*/,
    Tree<BuildT>& /*io_tree*/)
{
    // \todo BIFROST-TBD - We may support coalesced tile data for quantized values in the future.
}

template <typename BuildT>
typename std::enable_if<std::is_same<nanovdb::FpN, BuildT>::value>::type
addDataToTile(const float*        in_td,
              nanovdb::DitherLUT& in_lut,
              const TileIndex     in_i,
              const Codec&        in_codec,
              uint64_t            in_offset,
              Tree<BuildT>&       io_tree,
              const int           in_size = Traits::kTileDataSize) {
    io_tree.setTileDataOffset(in_i, in_offset);
    auto& io_td = io_tree.tileData(in_i);
    const uint8_t logBitWidth = static_cast<uint8_t>(in_codec.log2);
    io_td.m_logBitWidth = logBitWidth;
    io_td.m_isExtended = (in_size == Traits::kTileDataSize ? 0 : 1);
    const float min = in_codec.min, max = in_codec.max;
    io_td.init(min, max, static_cast<uint8_t>(static_cast<uint8_t>(1) << logBitWidth));
    // perform quantization relative to the values in the current leaf node
    int offset = 0;
    switch (logBitWidth) {
        case 0u: { // 1 bit
            auto*       dst    = io_td.data();
            const float encode = 1.0f / (max > min ? (max - min) : 1.f);
            if (in_size == Traits::kTileDataSize) {
                // non-extended case
                for (int j = 0; j < in_size / 8; ++j) {
                    uint8_t a = 0;
                    for (int k = 0; k < 8; ++k) {
                        a |= static_cast<uint8_t>(encode * (*in_td++ - min) + in_lut(offset++)) << k;
                    }
                    *dst++ = a;
                }
            }
            else {
                // extended case
                for (int j = 0; j < in_size; ) {
                    uint8_t a = 0;
                    for (int k = 0; k < 8; ++k, ++j) {
                        const bool ok = (j < Traits::kExtendedTileDataSize);
                        a |= static_cast<uint8_t>(std::min(encode * std::max(((ok ? in_td[j] : 0) - min), 0.f) + (ok ? in_lut(j) : 0.f), 1.f)) << k;
                    }
                    *dst++ = a;
                }
            }
        } break;
        case 1u: { // 2 bits
            auto*       dst    = io_td.data();
            const float encode = 3.0f / (max - min);
            if (in_size == Traits::kTileDataSize) {
                // non-extended case
                for (int j = 0; j < in_size / 4; ++j) {
                    auto a = static_cast<uint8_t>(encode * (*in_td++ - min) +
                        in_lut(offset++));
                    a |= static_cast<uint8_t>(static_cast<uint8_t>(encode * (*in_td++ - min) +
                        in_lut(offset++))
                        << 2);
                    a |= static_cast<uint8_t>(static_cast<uint8_t>(encode * (*in_td++ - min) +
                        in_lut(offset++))
                        << 4);
                    *dst++ = static_cast<uint8_t>(static_cast<uint8_t>(encode * (*in_td++ - min) +
                        in_lut(offset++))
                        << 6) |
                        a;
                }
            }
            else {
                // extended case
                for (int j = 0; j < in_size; ) {
                    auto a = static_cast<uint8_t>(encode * (in_td[j] - min) +
                        in_lut(j));
                    ++j;
                    bool ok = j < Traits::kExtendedTileDataSize;
                    a |= static_cast<uint8_t>(static_cast<uint8_t>(std::min(encode * std::max(((ok ? in_td[j] : 0.f) - min), 0.f) +
                        (ok ? in_lut(j) : 0.f), 3.f))
                        << 2);
                    ++j;
                    ok = j < Traits::kExtendedTileDataSize;
                    a |= static_cast<uint8_t>(static_cast<uint8_t>(std::min(encode * std::max(((ok ? in_td[j] : 0.f) - min), 0.f) +
                        (ok ? in_lut(j) : 0.f), 3.f))
                        << 4);
                    ++j;
                    ok = j < Traits::kExtendedTileDataSize;
                    *dst++ = static_cast<uint8_t>(static_cast<uint8_t>(std::min(encode * std::max(((ok ? in_td[j] : 0.f) - min), 0.f) +
                        (ok ? in_lut(j) : 0.f), 3.f))
                        << 6) |
                        a;
                    ++j;
                }
            }
        } break;
        case 2u: { // 4 bits
            auto*       dst    = io_td.data();
            const float encode = 15.0f / (max - min);
            if (in_size == Traits::kTileDataSize) {
                // non-extended case
                for (int j = 0; j < in_size / 2; ++j) {
                    auto a = static_cast<uint8_t>(encode * (*in_td++ - min) +
                        in_lut(offset++));
                    *dst++ = static_cast<uint8_t>(static_cast<uint8_t>(encode * (*in_td++ - min) +
                        in_lut(offset++))
                        << 4) |
                        a;
                }
            }
            else {
                // extended case
                for (int j = 0; j < in_size; ) {
                    auto a = static_cast<uint8_t>(encode * (in_td[j] - min) +
                        in_lut(j));
                    ++j;
                    const bool ok = j < Traits::kExtendedTileDataSize;
                    *dst++ = static_cast<uint8_t>(static_cast<uint8_t>(std::min(encode * std::max(((ok ? in_td[j] : 0.f) - min), 0.f) +
                        (ok ? in_lut(j) : 0.f), 15.f))
                        << 4) |
                        a;
                    ++j;
                }
            }
        } break;
        case 3u: { // 8 bits
            auto*       dst    = io_td.data();
            const float encode = 255.0f / (max - min);
            if (in_size == Traits::kTileDataSize) {
                // non-extended case
                for (int j = 0; j < in_size; ++j) {
                    *dst++ = static_cast<uint8_t>(encode * (*in_td++ - min) +
                        in_lut(offset++));
                }
            }
            else {
                // extended case
                for (int j = 0; j < in_size; ++j) {
                    *dst++ = static_cast<uint8_t>(std::min(encode * std::max((in_td[j] - min), 0.f) +
                        in_lut(j), 255.f));
                }
            }
        } break;
        default: { // 16 bits
            auto*        dst = reinterpret_cast<uint16_t*>(io_td.data());
			const double mind = static_cast<double>(min);
			const double maxd = static_cast<double>(max);
            const double encode =
                65535.0 / (maxd - mind); // note that double is required!
            if (in_size == Traits::kTileDataSize) {
                // non-extended case
                for (int j = 0; j < in_size; ++j) {
                    *dst++ = static_cast<uint16_t>(encode * (static_cast<double>(*in_td++) - mind) +
                        static_cast<double>(in_lut(offset++)));
                }
            }
            else {
                // extended case
                for (int j = 0; j < in_size; ++j) {
                    *dst++ = static_cast<uint16_t>(std::min(encode * (std::max(static_cast<double>(in_td[j]) - mind, 0.0)) +
                        static_cast<double>(in_lut(j)), 65535.0));
                }
            }
        }
    } // end switch
}

template <typename BuildT>
typename std::enable_if<std::is_same<nanovdb::Fp4, BuildT>::value ||
    std::is_same<nanovdb::Fp8, BuildT>::value ||
    std::is_same<nanovdb::Fp16, BuildT>::value>::type
initAddDataToTile(const TileIndex /*in_i*/,
    const uint64_t /*in_offset*/ ,
    const int /*in_tileDataStride*/ ,
    const int /*in_tileDataShift*/ ,
    Tree<BuildT>& /*io_tree*/)
{
    // \todo BIFROST-TBD - We may support coalesced tile data for quantized values in the future.
}

template <typename BuildT>
typename std::enable_if<std::is_same<nanovdb::Fp4, BuildT>::value ||
    std::is_same<nanovdb::Fp8, BuildT>::value ||
    std::is_same<nanovdb::Fp16, BuildT>::value>::type
addDataToTile(const float* /*in_td*/,
    nanovdb::DitherLUT& /*io_lut*/,
    const TileIndex /*in_i*/,
    const Codec& /*in_codec*/,
    const int /*in_tdStrideSrc*/,
    const int /*in_tdStrideDst*/,
    const nanovdb::Coord& /*in_ijk*/,
    const nanovdb::Coord& /*in_tdw*/,
    Tree<BuildT>& /*io_tree*/)
{
    // \todo BIFROST-TBD - We may support coalesced tile data for quantized values in the future.
}

template <typename BuildT>
typename std::enable_if<std::is_same<nanovdb::Fp4, BuildT>::value ||
    std::is_same<nanovdb::Fp8, BuildT>::value ||
    std::is_same<nanovdb::Fp16, BuildT>::value>::type
addDataToTile(const float* in_td, nanovdb::DitherLUT& in_lut, const TileIndex in_i, const Codec& /*in_codec*/, uint64_t in_offset, Tree<BuildT>& io_tree, const int in_size = Traits::kTileDataSize)
{
    io_tree.setTileDataOffset(in_i, in_offset);
    using FloatT = typename std::conditional<Tree<BuildT>::TileDataType::kBitWidth >= 16, double, float>::type;// 16 compression and higher requires double
    using ArrayT = typename Tree<BuildT>::TileDataType::ArrayType;
    static constexpr FloatT UNITS = FloatT((1 << Tree<BuildT>::TileDataType::kBitWidth) - 1);// # of unique non-zero values
    auto& io_td = io_tree.tileData(in_i);
    // compute extrema values
    float min = std::numeric_limits<float>::max(), max = -min;
    for (int i = 0; i < in_size; ++i) {
        const float v = in_td[i];
        if (v < min) min = v;
        if (v > max) max = v;
    }
    io_td.init(min, max, Tree<BuildT>::TileDataType::kBitWidth);
    // perform quantization relative to the values in the current leaf node
    const FloatT encode = (max > min ? (UNITS / static_cast<FloatT>(max - min)) : static_cast<FloatT>(1));
    auto* code = io_td.data();
    int offset = 0;
    if (std::is_same<nanovdb::Fp4, BuildT>::value) {// resolved at compile-time
        if (in_size == Traits::kTileDataSize) {
            // non-extended
            for (int j = 0; j < in_size / 2; ++j) {
                ArrayT tmp = static_cast<ArrayT>(encode * static_cast<FloatT>(*in_td++ - min) + static_cast<FloatT>(in_lut(offset++)));
                *code++ = static_cast<ArrayT>(static_cast<ArrayT>(encode * static_cast<FloatT>(*in_td++ - min) + static_cast<FloatT>(in_lut(offset++))) << 4) | tmp;
            }
        }
        else {
            // extended
            for (int j = 0; j < in_size; ) {
                ArrayT tmp = static_cast<ArrayT>(encode * static_cast<FloatT>(in_td[j] - min) + static_cast<FloatT>(in_lut(offset++)));
                ++j;
                const bool ok = j < Traits::kExtendedTileDataSize;
                *code++ = static_cast<ArrayT>(ok ? (static_cast<ArrayT>(encode * static_cast<FloatT>(in_td[j] - min) + static_cast<FloatT>(in_lut(offset++))) << 4) : 0) | tmp;
                ++j;
            }
        }
    }
    else {
        if (in_size == Traits::kTileDataSize) {
            // non-extended
            for (int j = 0; j < in_size; ++j) {
                *code++ = static_cast<ArrayT>(encode * static_cast<FloatT>(*in_td++ - min) + static_cast<FloatT>(in_lut(offset++)));
            }
        }
        else {
            for (int j = 0; j < in_size; ++j) {
                // extended
                code[j] = static_cast<ArrayT>(
                    encode * static_cast<FloatT>(*in_td++ - min) + static_cast<FloatT>(in_lut(offset++)));
            }
        }
    }
}


template <typename BuildT, typename ValueT, typename OracleT>
inline typename std::enable_if<!std::is_same<BuildT, nanovdb::FpN>::value>::type
initOracle(const nanovdb::GridClass /*in_gridClass*/,
    const ValueT /*in_defaultValue*/,
    OracleT& /*io_oracle*/)
{
}

template <typename BuildT, typename ValueT, typename OracleT>
inline typename std::enable_if<std::is_same<BuildT, nanovdb::FpN>::value>::type
initOracle(const nanovdb::GridClass in_gridClass,
    const ValueT in_voxelSize,
    OracleT& io_oracle)
{
    if (std::is_same<nanovdb::AbsDiff, OracleT>::value && io_oracle.getTolerance() < 0.0f) {// default tolerance for level set and fog volumes
        if (in_gridClass == nanovdb::GridClass::LevelSet) {
            io_oracle.setTolerance(0.1f * static_cast<float>(in_voxelSize));
        }
        else if (in_gridClass == nanovdb::GridClass::FogVolume) {
            io_oracle.setTolerance(0.01f);// range of FOG volumes: [0;1]
        }
        else {
            io_oracle.setTolerance(0.0f);
        }
    }
}

template <typename BuildT, typename ValueT, typename OracleT>
inline typename std::enable_if<!std::is_same<BuildT, nanovdb::FpN>::value>::type
compression(const ValueT* /*in_td*/,
    nanovdb::DitherLUT& /*in_lut*/,
    const OracleT& /*in_oracle*/,
    Codec& /*io_codec*/,
    const int in_size = Traits::kTileDataSize,
    const int in_tdSize = Traits::kTileDataSize)
{
    // Force usage of in_size to be able to set default arguments
    if (in_size == 0 || in_tdSize) {
    }
}


template <typename BuildT, typename OracleT>
inline typename std::enable_if<std::is_same<BuildT, nanovdb::FpN>::value>::type
compression(const float* in_td,
    nanovdb::DitherLUT& in_lut,
    const OracleT& in_oracle,
    Codec& io_codec,
    const int in_size = Traits::kTileDataSize,
    const int in_tdSize = Traits::kTileDataSize)
{
    float min = std::numeric_limits<float>::max(), max = -min;
    for (int j = 0; j < in_size; ++j) {
        float v = in_td[j];
        if (v < min) min = v;
        if (v > max) max = v;
    }
    io_codec.min            = min;
    io_codec.max            = max;
    const float range       = max - min;
    uint16_t    logBitWidth = 0; // 0,1,2,3,4 => 1,2,4,8,16 bits
    while (range > 0.0f && logBitWidth < 4u) {
        const uint32_t mask =
            (static_cast<uint32_t>(1) << (static_cast<uint32_t>(1) << logBitWidth)) - 1u;
        const float encode = mask / range;
        const float decode = range / mask;
        int         j      = 0;
        do {
            const float    exact  = in_td[j]; // exact value
            const uint32_t code   = static_cast<uint32_t>(encode * (exact - min) + in_lut(j));
            const float    approx = code * decode + min; // approximate value
            j += in_oracle(exact, approx) ? 1 : (in_size+1);
        } while (j < in_size);
        if (j == in_size) break;
        ++logBitWidth;
    }
    io_codec.log2 = logBitWidth;
    const uint32_t numCodeBits = (1u << logBitWidth) * in_tdSize;
    // Here we compute the number of bits required by the codes where the codes are stored with 32 bits each.
    // Note that this ensures 4-byte alignment.
    io_codec.size = static_cast<uint16_t>(sizeof(TileDataVariable) + 4 * (numCodeBits / 32 + (numCodeBits % 32 == 0 ? 0 : 1)));
} // compression


} // namespace nanovtt

#endif // NANOVTT_QUANTIZATION_H_HAS_BEEN_INCLUDED
