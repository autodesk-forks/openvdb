// Copyright Contributors to the OpenVDB Project
// SPDX-License-Identifier: MPL-2.0

//////////////////////////////////////////////////////////////////////////
///
/// @file SampleFromVoxels.h
///
/// @authors Ken Museth and Autodesk
///
///////////////////////////////////////////////////////////////////////////

#ifndef NANOVTT_SAMPLE_FROM_VOXELS_H_HAS_BEEN_INCLUDED
#define NANOVTT_SAMPLE_FROM_VOXELS_H_HAS_BEEN_INCLUDED

/// SIMD Intrinsic Headers
#if defined(NANOVDB_USE_INTRINSICS)
#if defined(_WIN32)
#include <intrin.h>
#elif defined(__GNUC__)
#if defined(__x86_64__) || defined(__i386__)
#include <x86intrin.h>
#elif defined(__ARM_NEON__)
#include <arm_neon.h>
#endif
#endif
#endif

#include <nanovdb/nanovtt/NanoVTT.h>
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
#include <nanovdb/util/SampleFromVoxels.h>
#if defined(__clang__)
#pragma clang diagnostic pop
#endif

namespace nanovtt {

// Forward declaration of sampler with specific polynomial orders
template<typename TreeT, int Order, bool UseCache = true>
class SampleFromVoxels;

/// @brief Factory free-function for a sampler of specific polynomial orders
///
/// @details This allows for the compact syntax:
/// @code
///   auto acc = grid.getAccessor();
///   auto smp = nanovdb::createSampler<1>( acc );
/// @endcode
template<int Order, typename AccT, bool UseCache = true>
__hostdev__ SampleFromVoxels<AccT, Order, UseCache> createSampler(const AccT& acc)
{
    return SampleFromVoxels<AccT, Order, UseCache>(acc);
}

// ------------------------------> TrilinearSampler <--------------------------------------

/// @brief Tri-linear sampler, i.e. first order, interpolator
template<typename AccT>
class TrilinearSampler
{
public:
    using ValueT = typename AccT::ValueType;
    using CoordT = typename AccT::CoordType;
    static const int ORDER = 1;

    const AccT& mAcc;

    /// @brief Protected constructor from a ReadAccessor
    __hostdev__ TrilinearSampler(const AccT& acc) : mAcc(acc) {}

    __hostdev__ const AccT& accessor() const { return mAcc; }

    /// @brief Extract the stencil of 8 values
    inline __hostdev__ void stencil(const CoordT& ijk, ValueT (&v)[2][2][2], nanovdb::Vec3f& vIJK, CoordT& tijk, int& shift, float& factor, const int in_maxDepth = Traits::kMaxDepth) const;

    template<typename RealT, template<typename...> class Vec3T>
    static inline __hostdev__ ValueT sample(const Vec3T<RealT> &uvw, const ValueT (&v)[2][2][2]);

    template<typename RealT, template<typename...> class Vec3T>
    static inline __hostdev__ Vec3T<ValueT> gradient(const Vec3T<RealT> &uvw, const ValueT (&v)[2][2][2], const int in_maxDepth = Traits::kMaxDepth);

    static inline __hostdev__ bool zeroCrossing(const ValueT (&v)[2][2][2]);
}; // TrilinearSamplerBase

template<typename AccT>
void TrilinearSampler<AccT>::stencil(const CoordT& ijk, ValueT (&v)[2][2][2], nanovdb::Vec3f& vIJK, CoordT& tijk, int& shift, float &factor, const int in_maxDepth) const
{
    v[0][0][0] = mAcc.getValue(ijk, in_maxDepth); // i, j, k

    if (!mAcc.found())
        return;

    // Account for varying voxelsize in the uvw offset
    tijk = mAcc.tijk();
    shift = mAcc.vw();
    const int depth = Traits::kMaxDepth - shift;
    const auto& ijk2 = mAcc.ijk();
    const float icw = Traits::invTileSideLUT(depth + 1);
    constexpr float ratio = static_cast<float>(Traits::kTileDataWidth) / Traits::kTileWidth;
    factor = ratio * icw;
    vIJK[0] = static_cast<float>(tijk[0] + (ijk2[0] << shift));
    vIJK[1] = static_cast<float>(tijk[1] + (ijk2[1] << shift));
    vIJK[2] = static_cast<float>(tijk[2] + (ijk2[2] << shift));
    // Here we must take the tile stride into account in the check for upper bounds.
    const int size = mAcc.tileDataStride() - 1;
    if (ijk2[0] < size && ijk2[1] < size && ijk2[2] < size) {
        // In this case all values can be looked up from the extended itself
        mAcc.getStencilExtendedValues(v);
    }
    else
    {
        const int vdw = (1 << shift);
        const nanovdb::Coord baseIJK = tijk + nanovdb::Coord(ijk2[0] << shift, ijk2[1] << shift, ijk2[2] << shift);
        v[0][0][1] = mAcc.getExtendedValue(baseIJK, { 0,0,vdw }, depth); // i, j, k + 1
        v[0][1][1] = mAcc.getExtendedValue(baseIJK, { 0,vdw,vdw }, depth); // i, j+1, k + 1
        v[0][1][0] = mAcc.getExtendedValue(baseIJK, { 0,vdw,0 }, depth); // i, j+1, k
        v[1][1][0] = mAcc.getExtendedValue(baseIJK, { vdw,vdw,0 }, depth); // i+1, j+1, k
        v[1][0][0] = mAcc.getExtendedValue(baseIJK, { vdw,0,0 }, depth); // i+1, j, k
        v[1][0][1] = mAcc.getExtendedValue(baseIJK, { vdw,0,vdw }, depth); // i+1, j, k + 1
        v[1][1][1] = mAcc.getExtendedValue(baseIJK, { vdw,vdw,vdw }, depth); // i+1, j+1, k + 1
    }
}

template<typename AccT>
template<typename RealT, template<typename...> class Vec3T>
typename AccT::ValueType TrilinearSampler<AccT>::sample(const Vec3T<RealT> &uvw, const ValueT (&v)[2][2][2])
{
    assert(uvw[0] >= 0.f);
    assert(uvw[1] >= 0.f);
    assert(uvw[2] >= 0.f);

    auto lerp = [](ValueT a, ValueT b, RealT w) { return a + ValueT(w) * (b - a); };
    return lerp(lerp(lerp(v[0][0][0], v[0][0][1], uvw[2]), lerp(v[0][1][0], v[0][1][1], uvw[2]), uvw[1]),
        lerp(lerp(v[1][0][0], v[1][0][1], uvw[2]), lerp(v[1][1][0], v[1][1][1], uvw[2]), uvw[1]),
        uvw[0]);
}

template<typename AccT>
template<typename RealT, template<typename...> class Vec3T>
Vec3T<typename AccT::ValueType> TrilinearSampler<AccT>::gradient(const Vec3T<RealT> &uvw, const ValueT (&v)[2][2][2], const int in_maxDepth)
{
    static_assert(nanovdb::is_floating_point<ValueT>::value, "TrilinearSampler::gradient requires a floating-point type");
    auto lerp = [](ValueT a, ValueT b, RealT w) { return a + ValueT(w) * (b - a); };

    ValueT D[4] = {v[0][0][1] - v[0][0][0], v[0][1][1] - v[0][1][0], v[1][0][1] - v[1][0][0], v[1][1][1] - v[1][1][0]};

    // Z component
    Vec3T<ValueT> grad(0, 0, lerp(lerp(D[0], D[1], uvw[1]), lerp(D[2], D[3], uvw[1]), uvw[0]));

    const ValueT w = ValueT(uvw[2]);
    D[0] = v[0][0][0] + D[0] * w;
    D[1] = v[0][1][0] + D[1] * w;
    D[2] = v[1][0][0] + D[2] * w;
    D[3] = v[1][1][0] + D[3] * w;

    // X component
    grad[0] = lerp(D[2], D[3], uvw[1]) - lerp(D[0], D[1], uvw[1]);

    // Y component
    grad[1] = lerp(D[1] - D[0], D[3] - D[2], uvw[0]);

    grad /= static_cast<typename AccT::ValueType>(Traits::voxelDataWidth(in_maxDepth));

    return grad;
}

template<typename AccT>
bool TrilinearSampler<AccT>::zeroCrossing(const ValueT (&v)[2][2][2])
{
    static_assert(nanovdb::is_floating_point<ValueT>::value, "TrilinearSampler::zeroCrossing requires a floating-point type");
    const bool less = v[0][0][0] < ValueT(0);
    return (less ^ (v[0][0][1] < ValueT(0))) ||
           (less ^ (v[0][1][1] < ValueT(0))) ||
           (less ^ (v[0][1][0] < ValueT(0))) ||
           (less ^ (v[1][0][0] < ValueT(0))) ||
           (less ^ (v[1][0][1] < ValueT(0))) ||
           (less ^ (v[1][1][1] < ValueT(0))) ||
           (less ^ (v[1][1][0] < ValueT(0)));
}

/// @brief Template specialization that does not use caching of stencil points
template<typename AccT>
class SampleFromVoxels<AccT, 1, false> : public TrilinearSampler<AccT>
{
    using BaseT = TrilinearSampler<AccT>;
    using ValueT = typename AccT::ValueType;
    using CoordT = typename AccT::CoordType;

public:

    /// @brief Construction from a ReadAccessor
    __hostdev__ SampleFromVoxels(const AccT& acc) : BaseT(acc) {}

    /// @note xyz is in index space space
    template<typename RealT, template<typename...> class Vec3T>
    inline __hostdev__ ValueT operator()(Vec3T<RealT> xyz, const int in_maxDepth = Traits::kMaxDepth) const;

    /// @note ijk is in index space space
    __hostdev__ ValueT operator()(const CoordT &ijk, const int in_maxDepth = Traits::kMaxDepth) const {return BaseT::mAcc.getValue(ijk, in_maxDepth);}

    /// @brief Return the gradient in index space.
    ///
    /// @warning Will only compile with floating point value types
    template<typename RealT, template<typename...> class Vec3T>
    inline __hostdev__ Vec3T<ValueT> gradient(Vec3T<RealT> xyz, const int in_maxDepth = Traits::kMaxDepth) const;

    /// @brief Return true if the tr-linear stencil has a zero crossing at the specified index position.
    ///
    /// @warning Will only compile with floating point value types
    template<typename RealT, template<typename...> class Vec3T>
    inline __hostdev__ bool zeroCrossing(Vec3T<RealT> xyz, const int in_maxDepth = Traits::kMaxDepth) const;

}; // SampleFromVoxels<AccT, 1, false>

/// @brief Template specialization with caching of stencil values
template<typename AccT>
class SampleFromVoxels<AccT, 1, true> : public TrilinearSampler<AccT>
{
    using BaseT = TrilinearSampler<AccT>;
    using ValueT = typename AccT::ValueType;
    using CoordT = typename AccT::CoordType;

    mutable CoordT mPos;
    mutable nanovdb::Vec3f mIJK;
    mutable CoordT mTIJK;
    mutable int mShift;
    mutable float mFactor;
    mutable int mMaxDepth;
    mutable ValueT mVal[2][2][2];

public:

    template<typename RealT, template<typename...> class Vec3T>
    __hostdev__ void cache(Vec3T<RealT>& xyz, const int in_maxDepth = Traits::kMaxDepth) const;

    /// @brief Construction from a ReadAccessor
    __hostdev__ SampleFromVoxels(const AccT& acc) : BaseT(acc), mPos(CoordT::max()), mIJK({}), mTIJK(), mShift(0), mFactor(0.f), mMaxDepth(0) {}

    /// @note xyz is in index space space
    template<typename RealT, template<typename...> class Vec3T>
    inline __hostdev__ ValueT operator()(Vec3T<RealT> xyz, const int in_maxDepth = Traits::kMaxDepth) const;

    // @note ijk is in index space space
    __hostdev__ ValueT operator()(const CoordT &ijk, const int in_maxDepth = Traits::kMaxDepth) const;

    /// @brief Return the gradient in index space.
    ///
    /// @warning Will only compile with floating point value types
    template<typename RealT, template<typename...> class Vec3T>
    inline __hostdev__ Vec3T<ValueT> gradient(Vec3T<RealT> xyz, const int in_maxDepth = Traits::kMaxDepth) const;

    /// @brief Return true if the tri-linear stencil has a zero crossing at the specified index position.
    ///
    /// @warning Will only compile with floating point value types
    template<typename RealT, template<typename...> class Vec3T>
    inline __hostdev__ bool zeroCrossing(Vec3T<RealT> xyz, const int in_maxDepth = Traits::kMaxDepth) const;

    /// @brief Return true if the cached tri-linear stencil has a zero crossing.
    ///
    /// @warning Will only compile with floating point value types
    __hostdev__ bool zeroCrossing() const { return BaseT::zeroCrossing(mVal); }

    ///////////////////////////////////////////////////////////
    // NON-STANDARD API BEGIN
    ///////////////////////////////////////////////////////////

    /// @brief Sample inside the already cached voxel.
    template<typename RealT, template<typename...> class Vec3T>
    inline __hostdev__ ValueT sampleCachedVoxel(Vec3T<RealT> xyz) const;

    inline __hostdev__ CoordT const& getCachedCoord() const { return mPos; }

    ///////////////////////////////////////////////////////////
    // NON-STANDARD API END
    ///////////////////////////////////////////////////////////

}; // SampleFromVoxels<AccT, 1, true>

template<typename AccT>
template<typename RealT, template<typename...> class Vec3T>
typename AccT::ValueType SampleFromVoxels<AccT, 1, true>::operator()(Vec3T<RealT> xyz, const int in_maxDepth) const
{
    // Appears to be faster to check for found here than always calling sample.
    this->cache(xyz, in_maxDepth);
    return BaseT::mAcc.found() ? BaseT::sample(xyz, mVal) : mVal[0][0][0];
}

template<typename AccT>
template<typename RealT, template<typename...> class Vec3T>
typename AccT::ValueType SampleFromVoxels<AccT, 1, true>::sampleCachedVoxel(Vec3T<RealT> xyz) const
{
    // Account for varying voxelsize in the uvw offset
    xyz[0] = (xyz[0] - mIJK[0]) * mFactor;
    xyz[1] = (xyz[1] - mIJK[1]) * mFactor;
    xyz[2] = (xyz[2] - mIJK[2]) * mFactor;
    return BaseT::sample(xyz, mVal);
}

template<typename AccT>
typename AccT::ValueType SampleFromVoxels<AccT, 1, true>::operator()(const CoordT &ijk, const int in_maxDepth) const
{
    return  (ijk == mPos && mMaxDepth == in_maxDepth) ? mVal[0][0][0] : BaseT::mAcc.getValue(ijk, in_maxDepth);
}

template<typename AccT>
template<typename RealT, template<typename...> class Vec3T>
Vec3T<typename AccT::ValueType> SampleFromVoxels<AccT, 1, true>::gradient(Vec3T<RealT> xyz, const int in_maxDepth) const
{
    this->cache(xyz, in_maxDepth);
    return BaseT::mAcc.found() ? BaseT::gradient(xyz, mVal, Traits::kMaxDepth - mShift) : Vec3T<typename AccT::ValueType>{};
}

template<typename AccT>
template<typename RealT, template<typename...> class Vec3T>
__hostdev__ bool SampleFromVoxels<AccT, 1, true>::zeroCrossing(Vec3T<RealT> xyz, const int in_maxDepth) const
{
    this->cache(xyz, in_maxDepth);
    return BaseT::mAcc.found() ? BaseT::zeroCrossing(mVal) : false;
}

template<typename AccT>
template<typename RealT, template<typename...> class Vec3T>
void SampleFromVoxels<AccT, 1, true>::cache(Vec3T<RealT>& xyz, const int in_maxDepth) const
{
    using CVT = typename CoordT::ValueType;
    CoordT ijk;
    // Note: if kMinTileCoord could be represented in floating point
    //       we could avoid the floorf here by subtracting, casting and adding.

#if defined(NANOVDB_USE_INTRINSICS) && !defined( __CUDA_ARCH__) && (defined(_MSC_VER) || defined(__GNUC__) || defined(__clang__))
    int32_t ijktemp[4];
    __m128 packedval = _mm_set_ps(xyz[0], xyz[1], xyz[2], xyz[0]);
    __m128 floorval = _mm_floor_ps(packedval);
    __m128i intcvt = _mm_cvtps_epi32(floorval);
    _mm_storeu_si128(reinterpret_cast<__m128i*>(ijktemp), intcvt);
    ijk[0] = static_cast<CVT>(ijktemp[3]);
    ijk[1] = static_cast<CVT>(ijktemp[2]);
    ijk[2] = static_cast<CVT>(ijktemp[1]);
#else
    for (int i = 0; i < 3; ++i) {
        ijk[i] = static_cast<CVT>(nanovdb::Floor(xyz[i]));
    }
#endif
    // Currently it seems that cached lookups are slightly slower than non-cached lookups.
#ifdef NANOVTT_USE_CACHED_SAMPLED_LOOKUPS
    const int mask = static_cast<int>(~((1u << mShift) - 1));
    const CoordT ijk2(
        ijk[0] & mask,
        ijk[1] & mask,
        ijk[2] & mask);

    if (!BaseT::mAcc.found() || ijk2 != mPos || in_maxDepth != mMaxDepth)
    {
        mMaxDepth = in_maxDepth;
        BaseT::stencil(ijk, mVal, mIJK, mTIJK, mShift, mFactor, in_maxDepth);
        mPos = { static_cast<int>(mIJK[0]), static_cast<int>(mIJK[1]), static_cast<int>(mIJK[2]) };
    }
#else
    mPos = ijk;
    mMaxDepth = in_maxDepth;
    BaseT::stencil(ijk, mVal, mIJK, mTIJK, mShift, mFactor, in_maxDepth);
#endif
    // Account for varying voxelsize in the uvw offset
    xyz[0] = (xyz[0] - mIJK[0]) * mFactor;
    xyz[1] = (xyz[1] - mIJK[1]) * mFactor;
    xyz[2] = (xyz[2] - mIJK[2]) * mFactor;
}

template<typename AccT>
template<typename RealT, template<typename...> class Vec3T>
typename AccT::ValueType SampleFromVoxels<AccT, 1, false>::operator()(Vec3T<RealT> xyz, const int in_maxDepth) const
{
    using CVT = typename CoordT::ValueType;
    nanovdb::Vec3f ijkOut;
    CoordT tijk;
    int shift;
    float factor;
    ValueT val[2][2][2];
    CoordT ijk;
#if defined(NANOVDB_USE_INTRINSICS) && !defined( __CUDA_ARCH__) && (defined(_MSC_VER) || defined(__GNUC__) || defined(__clang__))
    int32_t ijktemp[4];
    __m128 packedval = _mm_set_ps(xyz[0], xyz[1], xyz[2], xyz[0]);
    __m128 floorval = _mm_floor_ps(packedval);
    __m128i intcvt = _mm_cvtps_epi32(floorval);
    _mm_storeu_si128(reinterpret_cast<__m128i*>(ijktemp), intcvt);
    ijk[0] = static_cast<CVT>(ijktemp[3]);
    ijk[1] = static_cast<CVT>(ijktemp[2]);
    ijk[2] = static_cast<CVT>(ijktemp[1]);
#else
    for (int i = 0; i < 3; ++i) {
        ijk[i] = static_cast<CVT>(nanovdb::Floor(xyz[i]));
    }
#endif
    BaseT::stencil(ijk, val, ijkOut, tijk, shift, factor, in_maxDepth);
    // Account for varying voxelsize in the uvw offset
    xyz[0] = (xyz[0] - ijkOut[0]) * factor;
    xyz[1] = (xyz[1] - ijkOut[1]) * factor;
    xyz[2] = (xyz[2] - ijkOut[2]) * factor;
    // Appears to be faster to check for found here than always calling sample.
    return BaseT::mAcc.found() ? BaseT::sample(xyz, val) : val[0][0][0];
}

template<typename AccT>
template<typename RealT, template<typename...> class Vec3T>
inline Vec3T<typename AccT::ValueType> SampleFromVoxels<AccT, 1, false>::gradient(Vec3T<RealT> xyz, const int in_maxDepth) const
{
    ValueT val[2][2][2];
    using CVT = typename CoordT::ValueType;
    CoordT ijk;
#if defined(NANOVDB_USE_INTRINSICS) && !defined( __CUDA_ARCH__) && (defined(_MSC_VER) || defined(__GNUC__) || defined(__clang__))
    int32_t ijktemp[4];
    __m128 packedval = _mm_set_ps(xyz[0], xyz[1], xyz[2], xyz[0]);
    __m128 floorval = _mm_floor_ps(packedval);
    __m128i intcvt = _mm_cvtps_epi32(floorval);
    _mm_storeu_si128(reinterpret_cast<__m128i*>(ijktemp), intcvt);
    ijk[0] = static_cast<CVT>(ijktemp[3]);
    ijk[1] = static_cast<CVT>(ijktemp[2]);
    ijk[2] = static_cast<CVT>(ijktemp[1]);
#else
    for (int i = 0; i < 3; ++i) {
        ijk[i] = static_cast<CVT>(nanovdb::Floor(xyz[i]));
    }
#endif
    nanovdb::Vec3f vIJK;
    CoordT tijk;
    int shift;
    float factor;
    BaseT::stencil(ijk, val, vIJK, tijk, shift, factor, in_maxDepth);
    // Account for varying voxelsize in the uvw offset
    xyz[0] = (xyz[0] - vIJK[0]) * factor;
    xyz[1] = (xyz[1] - vIJK[1]) * factor;
    xyz[2] = (xyz[2] - vIJK[2]) * factor;
    return BaseT::gradient(xyz, val, Traits::kMaxDepth - shift);
}

template<typename AccT>
template<typename RealT, template<typename...> class Vec3T>
bool SampleFromVoxels<AccT, 1, false>::zeroCrossing(Vec3T<RealT> xyz, const int in_maxDepth) const
{
    ValueT val[2][2][2];
    using CVT = typename CoordT::ValueType;
    CoordT ijk;
#if defined(NANOVDB_USE_INTRINSICS) && !defined( __CUDA_ARCH__) && (defined(_MSC_VER) || defined(__GNUC__) || defined(__clang__))
    int32_t ijktemp[4];
    __m128 packedval = _mm_set_ps(xyz[0], xyz[1], xyz[2], xyz[0]);
    __m128 floorval = _mm_floor_ps(packedval);
    __m128i intcvt = _mm_cvtps_epi32(floorval);
    _mm_storeu_si128(reinterpret_cast<__m128i*>(ijktemp), intcvt);
    ijk[0] = static_cast<CVT>(ijktemp[3]);
    ijk[1] = static_cast<CVT>(ijktemp[2]);
    ijk[2] = static_cast<CVT>(ijktemp[1]);
#else
    for (int i = 0; i < 3; ++i) {
        ijk[i] = static_cast<CVT>(nanovdb::Floor(xyz[i]));
    }
#endif
    nanovdb::Vec3f vIJK;
    CoordT tijk;
    int shift;
    float factor;
    BaseT::stencil(ijk, val, vIJK, tijk, shift, factor, in_maxDepth);
    return BaseT::zeroCrossing(val);
}

} // namespace nanovtt

#endif // NANOVTT_SAMPLE_FROM_VOXELS_H_HAS_BEEN_INCLUDED
