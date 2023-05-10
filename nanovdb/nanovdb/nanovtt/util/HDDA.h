// Copyright Contributors to the OpenVDB Project
// SPDX-License-Identifier: MPL-2.0

/// @file HDDA.h
///
/// @authors Ken Museth and Autodesk

#ifndef NANOVTT_HDDA_H_HAS_BEEN_INCLUDED
#define NANOVTT_HDDA_H_HAS_BEEN_INCLUDED

#define ENFORCE_FORWARD_STEPPING

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wold-style-cast"
#pragma clang diagnostic ignored "-Wnewline-eof"
#endif
#include <nanovdb/util/HDDA.h>
#include <nanovdb/nanovtt/NanoVTT.h>
#if defined(__clang__)
#pragma clang diagnostic pop
#endif

namespace nanovtt {

/// @brief A Digital Differential Analyzer specialized for OpenVDB grids
/// @note Conceptually similar to Bresenham's line algorithm applied
/// to a 3D Ray intersecting OpenVDB nodes or voxels. Log2Dim = 0
/// corresponds to a voxel and Log2Dim a tree node of size 2^Log2Dim.
///
/// @note The Ray template class is expected to have the following
/// methods: test(time), t0(), t1(), invDir(), and  operator()(time).
/// See the example Ray class above for their definition.
///
/// @note The HDDA class is copied to nanovtt to support an additional
/// update function.
template<typename RayT, typename CoordT = nanovdb::Coord>
class HDDA
{
public:
    using RealType = typename RayT::RealType;
    using RealT = RealType;
    using Vec3Type = typename RayT::Vec3Type;
    using Vec3T = Vec3Type;
    using CoordType = CoordT;

    /// @brief Default ctor
    HDDA() = default;

    /// @brief ctor from ray and dimension at which the DDA marches
    __hostdev__ HDDA(const RayT& ray, int dim) { this->init(ray, dim); }

    /// @brief Re-initializes the HDDA
    __hostdev__ void init(const RayT& ray, RealT startTime, RealT maxTime, int dim)
    {
        assert(startTime <= maxTime);
        mDim = dim;
        mT0 = startTime;
        mT1 = maxTime;
        const Vec3T &pos = ray(mT0), &dir = ray.dir(), &inv = ray.invDir();
        mVoxel = nanovdb::RoundDown<CoordT>(pos) & (~(dim - 1));
        for (int axis = 0; axis < 3; ++axis) {
            if (dir[axis] <= RealT(0) && dir[axis] >= RealT(0)) { //handles dir = +/- 0
                mNext[axis] = nanovdb::Maximum<RealT>::value(); //i.e. disabled!
                mStep[axis] = 0;
            } else if (inv[axis] > 0) {
                mStep[axis] = 1;
                mNext[axis] = mT0 + (mVoxel[axis] + dim - pos[axis]) * inv[axis];
                mDelta[axis] = inv[axis];
            } else {
                mStep[axis] = -1;
                mNext[axis] = mT0 + (mVoxel[axis] - pos[axis]) * inv[axis];
                mDelta[axis] = -inv[axis];
            }
        }
    }

    /// @brief Simular to init above except it uses the bounds of the input ray
    __hostdev__ void init(const RayT& ray, int dim) { this->init(ray, ray.t0(), ray.t1(), dim); }

    ///////////////////////////////////////////////////////////
    // NON-STANDARD API BEGIN
    ///////////////////////////////////////////////////////////

    /// @brief Updates the HDDA to march with the specified dimension at the specified pos and voxel
    __hostdev__ bool update(const RayT& ray, int dim, const Vec3T& pos, CoordT voxel)
    {
        voxel &= (~(dim - 1));
        if (voxel == mVoxel && mDim == dim) return false;

        mDim = dim;
        mVoxel = voxel;

        for (int axis = 0; axis < 3; ++axis) {
            if (mStep[axis] == 0)
                continue;
            mNext[axis] = mT0 + (mVoxel[axis] - pos[axis]) * ray.invDir()[axis];
            if (mStep[axis] > 0)
                mNext[axis] += dim * ray.invDir()[axis];
        }

        return true;
    }

    ///////////////////////////////////////////////////////////
    // NON-STANDARD API END
    ///////////////////////////////////////////////////////////

    /// @brief Updates the HDDA to march with the specified dimension
    __hostdev__ bool update(const RayT& ray, int dim)
    {
        if (mDim == dim)
            return false;
        mDim = dim;
        const Vec3T &pos = ray(mT0), &inv = ray.invDir();
        mVoxel = nanovdb::RoundDown<CoordT>(pos) & (~(dim - 1));
        for (int axis = 0; axis < 3; ++axis) {
            if (mStep[axis] == 0)
                continue;
            mNext[axis] = mT0 + (mVoxel[axis] - pos[axis]) * inv[axis];
            if (mStep[axis] > 0)
                mNext[axis] += dim * inv[axis];
        }

        return true;
    }

    __hostdev__ int dim() const { return mDim; }

    /// @brief Increment the voxel index to next intersected voxel or node
    /// and returns true if the step in time does not exceed maxTime.
    __hostdev__ bool step()
    {
        const int axis = MinIndex(mNext);
        switch (axis) {
        case 0:
            return step<0>();
        case 1:
            return step<1>();
        default:
            return step<2>();
        }
    }

    /// @brief Return the index coordinates of the next node or voxel
    /// intersected by the ray. If Log2Dim = 0 the return value is the
    /// actual signed coordinate of the voxel, else it is the origin
    /// of the corresponding VDB tree node or tile.
    /// @note Incurs no computational overhead.
    __hostdev__ const CoordT& voxel() const { return mVoxel; }

    /// @brief Return the time (parameterized along the Ray) of the
    /// first hit of a tree node of size 2^Log2Dim.
    /// @details This value is initialized to startTime or ray.t0()
    /// depending on the constructor used.
    /// @note Incurs no computational overhead.
    __hostdev__ RealType time() const { return mT0; }

    /// @brief Return the maximum time (parameterized along the Ray).
    __hostdev__ RealType maxTime() const { return mT1; }

    /// @brief Return the time (parameterized along the Ray) of the
    /// second (i.e. next) hit of a tree node of size 2^Log2Dim.
    /// @note Incurs a (small) computational overhead.
    __hostdev__ RealType next() const
    {
#if 1
        return Min(mT1, Min(mNext[0], Min(mNext[1], mNext[2])));
#else
#if 1 //def __CUDA_ARCH__
        return fminf(mT1, fminf(mNext[0], fminf(mNext[1], mNext[2])));
#else
        return std::min(mT1, std::min(mNext[0], std::min(mNext[1], mNext[2])));
#endif
#endif
    }

private:
    // helper to implement the general form
    template<int axis>
    __hostdev__ bool step()
    {
#ifdef ENFORCE_FORWARD_STEPPING
        //if (mNext[axis] <= mT0) mNext[axis] += mT0 - mNext[axis] + fmaxf(mNext[axis]*1.0e-6f, 1.0e-6f);
        //if (mNext[axis] <= mT0) mNext[axis] += mT0 - mNext[axis] + (mNext[axis] + 1.0f)*1.0e-6f;
        if (mNext[axis] <= mT0) {
            mNext[axis] += mT0 - 0.999999f * mNext[axis] + 1.0e-6f;
        }
#endif
        mT0 = mNext[axis];
        mNext[ axis] += mDim * mDelta[axis];
        mVoxel[axis] += mDim * mStep[ axis];
        return mT0 <= mT1;
    }

    int32_t mDim;
    RealT   mT0, mT1; // min and max allowed times
    CoordT  mVoxel, mStep; // current voxel location and step to next voxel location
    Vec3T   mDelta, mNext; // delta time and next time
}; // class HDDA


static constexpr float kHDDADelta = 0.0001f;
static constexpr int kHDDAMaxBisections = 10;

/////////////////////////////////////////// ZeroCrossing ////////////////////////////////////////////

/// @brief returns true if the ray intersects a zero-crossing at the voxel level of the grid in the accessor
///        The empty-space ray-marching is performed at all levels of the tree using an
///        HDDA. If an intersection is detected, then ijk is updated with the index coordinate of the closest
///        voxel after the intersection point, v contains the grid values at ijk, and t is set to the time of
///        the intersection along the ray.
///
///        The NanoVTT::ZeroCrossing implementation is a patch of NanoVDB's ZeroCrossing that works for adaptive volumes.
template<typename RayT, typename AccT>
inline __hostdev__ bool ZeroCrossing(RayT& ray, AccT& acc, nanovdb::Coord& ijk, typename AccT::ValueType& v, float& t) {
    if (!ray.clip(acc.root().bbox()) || ray.t1() > 1e20f)
        return false; // clip ray to bbox
    ray.setMinTime(ray.t0());
    ijk = nanovdb::RoundDown<nanovdb::Coord>(ray.start()); // first hit of bbox
    const auto                 v0 = acc.getValue(ijk);
    HDDA<RayT, nanovdb::Coord> hdda(ray, acc.getDim());

    // We leave out the call to is_active as it always returns true inside the
    // bbox of the level set.
    while (hdda.step()) {
        const float v1 = acc.getValue(hdda.voxel());
        ijk = hdda.voxel();
        // We need to update the dim and adjust the voxel before we retrieve the value
        // as we could have entered a finer voxel.
        hdda.update(ray, acc.getDim());
        if (ijk != hdda.voxel()) {
            v = acc.getValue(hdda.voxel());
        }
        else {
            v = v1;
        }
        // At resolution jumps there can be inconsistencies in the sampling
        // that lead to false positives for the zero crossing. The second
        // check avoids this.
        if (v * v0 < 0 && v * v1 >= 0) { // zero crossing
            ijk = hdda.voxel();
            t   = hdda.time();
            return true;
        }
    }
    return false;
}

template<typename RayT, typename SamplerT>
inline __hostdev__ bool zeroCrossingSubVoxelQuadratic(
    const float isoValue, RayT& ray, SamplerT& sampler, float& t0)
{
#ifndef __CUDACC__
    using std::min;
    using std::max;
#endif

    using Vec3T = typename RayT::Vec3Type;
    using RealT = typename RayT::RealType;
    if (!ray.clip(sampler.mAcc.root().bbox()) || ray.t1() > 1e20f)
        return false; // clip ray to bbox

    // First we step delta along the ray to make sure we end up inside the voxel we're going to step through.
    // Then we can just use dim from that voxel. We also ensure that the first voxel is always inside the bbox.
    ray.setMinTime(ray.t0() + kHDDADelta);
    if (ray.t0() > ray.t1())
        return false;

    float pVal = sampler.mAcc.getValue(nanovdb::RoundDown<nanovdb::Coord>(ray.eye() + ray.dir() * ray.t0()));
    bool inside = (pVal < isoValue);
    float pTmin = ray.t0();

    // Setup HDDA
    // The assumption here now is that we always start inside a voxel inside the bbox and therefore
    // the hdda.dim() and hdda.voxel(), hdda.time() and hdda.next() will always correspond to the current voxel.
    HDDA<RayT, nanovdb::Coord> hdda(ray, sampler.mAcc.getDim());

    do {
        Vec3T ijkf = Vec3T(static_cast<RealT>(hdda.voxel().x()), static_cast<RealT>(hdda.voxel().y()), static_cast<RealT>(hdda.voxel().z()));
        float dimf = static_cast<RealT>(hdda.dim());

        // Here we first make sure that the ray start point is clamped to be slightly inside the hdda voxel.
        // This may not be the cause due to numerical roundoff error and can lead to false positives.
        Vec3T pos = ray.eye() + ray.dir() * hdda.time();
        pos[0] = Min(Max(ijkf[0] + kHDDADelta, pos[0]), ijkf[0] + dimf);
        pos[1] = Min(Max(ijkf[1] + kHDDADelta, pos[1]), ijkf[1] + dimf);
        pos[2] = Min(Max(ijkf[2] + kHDDADelta, pos[2]), ijkf[2] + dimf);

        // Retrieve stencil for voxel.
        {
            // Note that this call modifies pos.
            Vec3T pos2 = pos;
            sampler.cache(pos2);
        }
        // Update dim in hdda.
        hdda.update(ray, sampler.mAcc.getDim(), pos, sampler.getCachedCoord());
        ijkf = Vec3T(static_cast<RealT>(hdda.voxel().x()), static_cast<RealT>(hdda.voxel().y()), static_cast<RealT>(hdda.voxel().z()));
        dimf = static_cast<RealT>(hdda.dim());

        // Clamp end point of voxel to be exactly within voxel such that we can re-use
        // the stencil for interpolation.
        t0 = hdda.next();

        pos = ray.eye() + ray.dir() * t0;
        pos[0] = Min(Max(ijkf[0], pos[0]), ijkf[0] + dimf);
        pos[1] = Min(Max(ijkf[1], pos[1]), ijkf[1] + dimf);
        pos[2] = Min(Max(ijkf[2], pos[2]), ijkf[2] + dimf);

        float pVal2 = (sampler.mAcc.found() ? sampler.sampleCachedVoxel(pos) : pVal);
        const bool inside2 = (pVal2 < isoValue);

        if (sampler.mAcc.found() && (inside != inside2 || sampler.zeroCrossing()))
        {
            // Using the secant method directly here will often not converge
            // because there will be a local minimum or maximum in the level set
            // function along the ray which can prevent features like corners
            // from being captured correctly. Instead we make a quadratic fit
            // to the cubic level set function along the ray and find the
            // closest root of this analytically.

            // Fit a quadratic with re-scaling of the equation
            if (pVal <= 0 && pVal >= 0) {
                t0 = pTmin;
                return true;
            }
            else {
                // Evaluate at midpoint between a and b.
                pos = ray.eye() + ray.dir() * (pTmin + t0) * 0.5f;
                pos[0] = Min(Max(ijkf[0], pos[0]), ijkf[0] + dimf);
                pos[1] = Min(Max(ijkf[1], pos[1]), ijkf[1] + dimf);
                pos[2] = Min(Max(ijkf[2], pos[2]), ijkf[2] + dimf);
                float pValHalf = sampler.sampleCachedVoxel(pos);

                // normalize function values, so our quadratic is interpolating
                // (1,fHalf,f1)
                float r0 = 1.f / pVal;
                float fHalf = r0 * pValHalf, f1 = r0 * pVal2;

                float t = 4 * fHalf - f1 - 1; // equal to f'(0)+2
                float d = t * t - 4 * f1; // note: absolutely has to be >=0 if f1<=0
                                    // (i.e. we saw a sign change)
                if (d >= 0 && (t <= 0 || f1 <= 0)) {
                    d = sqrt(d) - t;
                    d = 2.f / (2 + d);
                    t0 = d * (t0 - pTmin) + pTmin;
                    return true;
                }
            } // if (pVal <= 0 && pVal >= 0) .. else

        } // else if (sampler.zeroCrossing())

        pVal = pVal2;
        pTmin = t0;
    } while (hdda.step());

    return false;
} // zeroCrossingSubVoxelQuadratic


template<typename RayT, typename SamplerT>
inline __hostdev__ bool zeroCrossingSubVoxel(
    const float isoValue, RayT& ray, SamplerT& sampler, float& t0)
{
#ifndef __CUDACC__
    using std::min;
    using std::max;
#endif

    using Vec3T = typename RayT::Vec3Type;
    using RealT = typename RayT::RealType;
    if (!ray.clip(sampler.mAcc.root().bbox()) || ray.t1() > 1e20f)
        return false; // clip ray to bbox

    // First we step delta along the ray to make sure we end up inside the voxel we're going to step through.
    // Then we can just use dim from that voxel. We also ensure that the first voxel is always inside the bbox.
    ray.setMinTime(ray.t0() + kHDDADelta);
    if (ray.t0() > ray.t1())
        return false;

    float pVal = sampler.mAcc.getValue(nanovdb::RoundDown<nanovdb::Coord>(ray.eye() + ray.dir() * ray.t0()));
    bool inside = (pVal < isoValue);
    float pTmin = ray.t0();

    // Setup HDDA
    // The assumption here now is that we always start inside a voxel inside the bbox and therefore
    // the hdda.dim() and hdda.voxel(), hdda.time() and hdda.next() will always correspond to the current voxel.
    HDDA<RayT, nanovdb::Coord> hdda(ray, sampler.mAcc.getDim());

    do {
        Vec3T ijkf = Vec3T(static_cast<RealT>(hdda.voxel().x()), static_cast<RealT>(hdda.voxel().y()), static_cast<RealT>(hdda.voxel().z()));
        float dimf = static_cast<RealT>(hdda.dim());

        // Here we first make sure that the ray start point is clamped to be slightly inside the hdda voxel.
        // This may not be the cause due to numerical roundoff error and can lead to false positives.
        Vec3T pos = ray.eye() + ray.dir() * hdda.time();
        pos[0] = Min(Max(ijkf[0] + kHDDADelta, pos[0]), ijkf[0] + dimf);
        pos[1] = Min(Max(ijkf[1] + kHDDADelta, pos[1]), ijkf[1] + dimf);
        pos[2] = Min(Max(ijkf[2] + kHDDADelta, pos[2]), ijkf[2] + dimf);

        // Retrieve stencil for voxel.
        {
            // Note that this call modifies pos.
            Vec3T pos2 = pos;
            sampler.cache(pos2);
        }
        // Update dim in hdda.
        hdda.update(ray, sampler.mAcc.getDim(), pos, sampler.getCachedCoord());
        ijkf = Vec3T(static_cast<RealT>(hdda.voxel().x()), static_cast<RealT>(hdda.voxel().y()), static_cast<RealT>(hdda.voxel().z()));
        dimf = static_cast<RealT>(hdda.dim());

        // Clamp end point of voxel to be exactly within voxel such that we can re-use
        // the stencil for interpolation.
        t0 = hdda.next();

        pos = ray.eye() + ray.dir() * t0;
        pos[0] = Min(Max(ijkf[0], pos[0]), ijkf[0] + dimf);
        pos[1] = Min(Max(ijkf[1], pos[1]), ijkf[1] + dimf);
        pos[2] = Min(Max(ijkf[2], pos[2]), ijkf[2] + dimf);

        float pVal2 = (sampler.mAcc.found() ? sampler.sampleCachedVoxel(pos) : pVal);
        const bool inside2 = (pVal2 < isoValue);

        if (inside != inside2)
        {
            // If different sign at endpoints then perform bisection for a fixed number of iterations,
            // and end with linear interpolation to find surface position.
            // To improve the convergence we take the level set phi values into account in the bisection.
            float a = pTmin;
            float b = t0;
            float phiA = static_cast<float>(fabs(pVal));
            float phiB = static_cast<float>(fabs(pVal2));
            int step2 = 0;

            // Bisection will always converge in a monotonic manner, but (depending on the endpoints)
            // the value at the new points (and hence the final interpolated point)
            // may be further from the zero crossing than the end-point values.
            while (step2 < kHDDAMaxBisections
                && (b - a) > kHDDADelta) {
                // Here we utilize the level set values to accelerate the bisection.
                // Traditional bisection:
                //const float m = a + (b - a) * 0.5f;
                // Bisection using level set values:
                assert(phiA + phiB > 0);
                const float m = a + phiA / (phiA + phiB) * (b - a);
                pos = ray.eye() + ray.dir() * m;
                // Because we're in the interior of the interval, we can most likely leave out these checks,
                // but it doesn't appear to make a difference for performance.
                pos[0] = Min(Max(ijkf[0], pos[0]), ijkf[0] + dimf);
                pos[1] = Min(Max(ijkf[1], pos[1]), ijkf[1] + dimf);
                pos[2] = Min(Max(ijkf[2], pos[2]), ijkf[2] + dimf);

                // Use accelerated sampler call.
                // Also ensures sampler does not move outside of voxel.
                pVal2 = sampler.sampleCachedVoxel(pos);
                const bool sb = (pVal2 < isoValue);
                if (inside == sb) {
                    a = m;
                    phiA = static_cast<float>(fabs(pVal2));
                }
                else {
                    b = m;
                    phiB = static_cast<float>(fabs(pVal2));
                }
                ++step2;
            } // while (step2 < kHDDAMaxBisections && (b - a) > kHDDADelta)

            // Compute final intersection by linear interpolation
            t0 = a + phiA / (phiA + phiB) * (b - a);
            return true;
        } // else if (inside != inside2)
        else if (sampler.mAcc.found() && sampler.zeroCrossing())
        {
            // Using the secant method directly here will often not converge
            // because there will be a local minimum or maximum in the level set
            // function along the ray which can prevent features like corners
            // from being captured correctly. Instead we make a quadratic fit
            // to the cubic level set function along the ray and find the
            // closest root of this analytically.

            float a = pTmin;
            float b = t0;
            float phiA = pVal;
            float phiB = pVal2;

            // Evaluate at midpoint between a and b.
            float c = (a + b) * 0.5f;
            pos = ray.eye() + ray.dir() * c;
            pos[0] = Min(Max(ijkf[0], pos[0]), ijkf[0] + dimf);
            pos[1] = Min(Max(ijkf[1], pos[1]), ijkf[1] + dimf);
            pos[2] = Min(Max(ijkf[2], pos[2]), ijkf[2] + dimf);
            float phiC = sampler.sampleCachedVoxel(pos);

            // Find the coefficients of the quadratic assuming [a;b] == [0;1]
            // (we adjust for that later).
            float pa = 2 * phiA - 4 * phiC + 2 * phiB;
            float pb = -3 * phiA + 4 * phiC - phiB;
            float pc = phiA;

            float d = pb * pb - 4 * pa * pc;

            // We only treat the case where d>0 and we have two roots.
            // We need to check if the roots are in the interval (a;b) and if
            // so, pick the one closest to a.
            if (d > 0 && pa > 0) {
                // We have an intersection
                d = nanovdb::Sqrt(d);
                float r1 = (-pb - d) / (2 * pa);
                float r2 = (-pb + d) / (2 * pa);
                r1 = ((r1 < 0) ? 2 : r1);
                r2 = ((r2 < 0) ? 2 : r2);
                d = Min(r1, r2);
                if (d >= 0.f && d <= 1.f) {
                    t0 = d * (b - a) + a;
                    return true;
                }
            }
        } // else if (sampler.zeroCrossing())

        pVal = pVal2;
        pTmin = t0;
    } while (hdda.step());

    return false;
} // zeroCrossingSubVoxel

} // namespace nanovtt

#endif // NANOVTT_HDDA_HAS_BEEN_INCLUDED
