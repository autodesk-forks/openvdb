// Copyright Contributors to the OpenVDB Project
// SPDX-License-Identifier: MPL-2.0

/*!
    \file   PNanoVTT.h

    \authors Andrew Reidmeyer and Autodesk

    \brief  PNanoVTT.h contains modified code from PNanoVDB.h and is a portable C99/GLSL/HLSL/OpenCL port of NanoVTT.h.

    \details Limitations
             \li Quantization is currently not supported by PNanoVTT, only NanoVTT.
             \li NANOVTT_SUPPORT_COALESCED_TILES cannot currently be used in combination with sampling for PNanoVTT.
             \li Currently only extended tile data (5x5x5) is supported by the PNanoVTT sampling.

*/

#ifndef NANOVTT_PNANOVTT_H_HAS_BEEN_INCLUDED
#define NANOVTT_PNANOVTT_H_HAS_BEEN_INCLUDED

#ifdef PNANOVDB_C
// Disable warnings for PNanoVDB.h
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-parameter"
#pragma clang diagnostic ignored "-Wold-style-cast"
#pragma clang diagnostic ignored "-Wcast-align"
#pragma clang diagnostic ignored "-Wshadow"
#pragma clang diagnostic ignored "-Wimplicit-float-conversion"
#pragma clang diagnostic ignored "-Wfloat-conversion"
#endif
#include <nanovdb/PNanoVDB.h>
#if defined(__clang__)
#pragma clang diagnostic pop
#endif
#endif // PNANOVDB_C

// Strides in bytes
#define PNANOVTT_TRAITS_FLOAT_EXTENDED_TILE_DATA_STRIDE_X 4
#define PNANOVTT_TRAITS_FLOAT_EXTENDED_TILE_DATA_STRIDE_Y 20
#define PNANOVTT_TRAITS_FLOAT_EXTENDED_TILE_DATA_STRIDE_Z 100

#if defined(PNANOVDB_C)
#include <assert.h>
#endif

// Disable warnings for PNanoVTT.h : required to conform to C99
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wold-style-cast"
#endif
// ------------------------------------------------ Basic Texture Types and Functions -----------------------------------------------------------

#if defined(PNANOVDB_C)

struct pnanovtt_texture_t
{
    const uint32_t* texture;
    const uint64_t dim;
};
PNANOVDB_STRUCT_TYPEDEF(pnanovtt_texture_t)

PNANOVDB_FORCE_INLINE float pnanovtt_read_texture_float(pnanovtt_texture_t tex, PNANOVDB_IN(pnanovdb_coord_t) coord)
{
    const uint64_t index = tex.dim * (tex.dim * PNANOVDB_DEREF(coord).z + PNANOVDB_DEREF(coord).y) + PNANOVDB_DEREF(coord).x;
    return ((const float*)(tex.texture))[index];
}

#elif defined(PNANOVDB_HLSL)
#define pnanovtt_texture_t Texture3D
PNANOVDB_FORCE_INLINE float pnanovtt_read_texture_float(pnanovtt_texture_t tex, PNANOVDB_IN(pnanovdb_coord_t) coord)
{
    return tex.Load(int4(PNANOVDB_DEREF(coord).x, PNANOVDB_DEREF(coord).y, PNANOVDB_DEREF(coord).z, 0));
}
#elif defined(PNANOVDB_GLSL)
#define pnanovtt_texture_t sampler3D
PNANOVDB_FORCE_INLINE float pnanovtt_read_texture_float(pnanovtt_texture_t tex, PNANOVDB_IN(pnanovdb_coord_t) coord)
{
    return texelFetch(tex, PNANOVDB_DEREF(coord), 0).x;
}
#elif defined(PNANOVDB_CL)
#define pnanovtt_texture_t /*__read_only*/ image3d_t
PNANOVDB_FORCE_INLINE float pnanovtt_read_texture_float(pnanovtt_texture_t tex, PNANOVDB_IN(pnanovdb_coord_t) coord)
{
    return read_imagef(tex, (int4)(PNANOVDB_DEREF(coord).x, PNANOVDB_DEREF(coord).y, PNANOVDB_DEREF(coord).z, 0)).x;
}
#endif

// ------------------------------------------------ Basic Functions -----------------------------------------------------------

PNANOVDB_FORCE_INLINE pnanovdb_int32_t pnanovtt_max_int32(pnanovdb_int32_t a, pnanovdb_int32_t b) { return a > b ? a : b; }
PNANOVDB_FORCE_INLINE pnanovdb_int32_t pnanovtt_min_int32(pnanovdb_int32_t a, pnanovdb_int32_t b) { return a < b ? a : b; }

#if defined(PNANOVDB_C)
PNANOVDB_FORCE_INLINE float pnanovdb_double_as_float(double v) { return (float)v; }
PNANOVDB_FORCE_INLINE float pnanovtt_abs(float f) { return abs(f); }
#elif defined(PNANOVDB_HLSL)
PNANOVDB_FORCE_INLINE float pnanovdb_double_as_float(double v) { return float(v); }
PNANOVDB_FORCE_INLINE float pnanovtt_abs(float f) { return abs(f); }
#elif defined(PNANOVDB_GLSL)
PNANOVDB_FORCE_INLINE float pnanovdb_double_as_float(double v) { return float(v); }
PNANOVDB_FORCE_INLINE float pnanovtt_abs(float f) { return abs(f); }
#elif defined(PNANOVDB_CL)
PNANOVDB_FORCE_INLINE float pnanovdb_double_as_float(double v) { return (float)v; }
PNANOVDB_FORCE_INLINE float pnanovtt_abs(float f) { return fabs(f); }
#endif

PNANOVDB_FORCE_INLINE pnanovdb_vec3_t pnanovtt_read_vec3(pnanovdb_buf_t buf, pnanovdb_address_t address)
{
    pnanovdb_vec3_t _vec3;
    _vec3.x = pnanovdb_uint32_as_float(
        pnanovdb_buf_read_uint32(buf, address.byte_offset));
    _vec3.y = pnanovdb_uint32_as_float(
        pnanovdb_buf_read_uint32(buf, address.byte_offset + 4));
    _vec3.z = pnanovdb_uint32_as_float(
        pnanovdb_buf_read_uint32(buf, address.byte_offset + 8));
    return _vec3;
}

// ------------------------------------------------ Basic Types -----------------------------------------------------------

#define PNANOVTT_TILEINDEX pnanovdb_uint32_t
#define PNANOVTT_TILEINDEX_SIZE 4
#define PNANOVTT_INVALID_TILEINDEX 0xFFFFFFFF

// basic types, type conversion
#if defined(PNANOVDB_C)
PNANOVDB_FORCE_INLINE PNANOVTT_TILEINDEX pnanovtt_min_tileindex(PNANOVTT_TILEINDEX a, PNANOVTT_TILEINDEX b) { return a < b ? a : b; }
#elif defined(PNANOVDB_HLSL)
PNANOVTT_TILEINDEX pnanovtt_min_tileindex(PNANOVTT_TILEINDEX a, PNANOVTT_TILEINDEX b) { return min(a, b); }
#elif defined(PNANOVDB_GLSL)
PNANOVTT_TILEINDEX pnanovtt_min_tileindex(PNANOVTT_TILEINDEX a, PNANOVTT_TILEINDEX b) { return min(a, b); }
#elif defined(PNANOVDB_CL)
PNANOVTT_TILEINDEX pnanovtt_min_tileindex(PNANOVTT_TILEINDEX a, PNANOVTT_TILEINDEX b) { return min(a, b); }
#endif

// ------------------------------------------------ Address Type -----------------------------------------------------------

#if defined(PNANOVDB_ADDRESS_32)
PNANOVDB_FORCE_INLINE pnanovdb_address_t pnanovdb_address_offset_neg64(pnanovdb_address_t address, pnanovdb_uint64_t byte_offset)
{
    pnanovdb_address_t ret = address;
    // lose high bits
    ret.byte_offset -= pnanovdb_uint64_low(byte_offset);
    return ret;
}

PNANOVDB_FORCE_INLINE pnanovdb_address_t pnanovdb_address_product(pnanovdb_address_t address, pnanovdb_uint32_t multiplier)
{
    pnanovdb_address_t ret = address;
    ret.byte_offset *= multiplier;
    return ret;
}

#elif defined(PNANOVDB_ADDRESS_64)
PNANOVDB_FORCE_INLINE pnanovdb_address_t pnanovdb_address_offset_neg64(pnanovdb_address_t address, pnanovdb_uint64_t byte_offset)
{
    pnanovdb_address_t ret = address;
    ret.byte_offset -= byte_offset;
    return ret;
}

PNANOVDB_FORCE_INLINE pnanovdb_address_t pnanovdb_address_product(pnanovdb_address_t address, pnanovdb_uint64_t multiplier)
{
    pnanovdb_address_t ret = address;
    ret.byte_offset *= multiplier;
    return ret;
}
#endif


// ------------------------------------------------ Buffer -----------------------------------------------------------

#if defined(PNANOVDB_BUF_C)
#include <stdint.h>
typedef uint32_t pnanovtt_grid_type_t;
#define PNANOVTT_GRID_TYPE_GET(grid_typeIn, nameIn) pnanovtt_grid_type_constants[grid_typeIn].nameIn
#elif defined(PNANOVDB_BUF_HLSL)
#define pnanovtt_grid_type_t uint
#define PNANOVTT_GRID_TYPE_GET(grid_typeIn, nameIn) pnanovtt_grid_type_constants[grid_typeIn].nameIn
#elif defined(PNANOVDB_BUF_GLSL)
#define pnanovtt_grid_type_t uint
#define PNANOVTT_GRID_TYPE_GET(grid_typeIn, nameIn) pnanovtt_grid_type_constants[grid_typeIn].nameIn
#elif defined(PNANOVDB_BUF_CL)
#define pnanovtt_grid_type_t uint
#define PNANOVTT_GRID_TYPE_GET(grid_typeIn, nameIn) pnanovtt_grid_type_constants[grid_typeIn].nameIn
#endif


// ------------------------------------------------ Core Structures -----------------------------------------------------------

#define PNANOVTT_MAGIC_NUMBER 0x4e616e6f56545430UL // "NanoVTT0" in hex - little endian (uint64_t)

#define PNANOVTT_MAJOR_VERSION_NUMBER 1 // reflects changes to the ABI and hence also the file format
#define PNANOVTT_MINOR_VERSION_NUMBER 0 //  reflects changes to the API but not ABI
#define PNANOVTT_PATCH_VERSION_NUMBER 0 //  reflects changes that does not affect the ABI or API

#define PNANOVTT_TRAITS_MAX_DEPTH 28
#define PNANOVTT_TRAITS_TILE_DATA_WIDTH 4
#define PNANOVTT_TRAITS_TILE_DATA_SIZE 64
#define PNANOVTT_TRAITS_EXTENDED_TILE_DATA_WIDTH 5
#define PNANOVTT_TRAITS_EXTENDED_TILE_DATA_SIZE 125
#define PNANOVTT_TRAITS_TILE_DOMAIN_WIDTH 1073741824

#ifdef PNANOVDB_ADDRESS_32
#define PNANOVTT_TRAITS_TILE_DATA_STRIDE_SHIFT_HIGH 16
#define PNANOVTT_TRAITS_TILE_DATA_SHIFT_SHIFT_HIGH 24
#define PNANOVTT_TRAITS_TEXTURE_COORD_SHIFT 16
#if defined(PNANOVDB_C)
#define PNANOVTT_TRAITS_TILE_DATA_OFFSET_MASK_LOW 0xFFFFFFFFu
#define PNANOVTT_TRAITS_TILE_DATA_OFFSET_MASK_HIGH ((1u << 16u) - 1u)
#define PNANOVTT_TRAITS_TILE_DATA_STRIDE_MASK_HIGH \
    (((~((1u << 16u) - 1u)) << 8u) >> 8u)
#define PNANOVTT_TRAITS_TEXTURE_COORD_MASK ((1u << 16u) - 1u)
#else
#define PNANOVTT_TRAITS_TILE_DATA_OFFSET_MASK_LOW 0xFFFFFFFFu
#define PNANOVTT_TRAITS_TILE_DATA_OFFSET_MASK_HIGH ((1u << 16) - 1u)
#define PNANOVTT_TRAITS_TILE_DATA_STRIDE_MASK_HIGH \
    (((~((1 << 16) - 1)) << 8) >> 8)
#define PNANOVTT_TRAITS_TEXTURE_COORD_MASK ((1u << 16) - 1u)
#endif
#else
#define PNANOVTT_TRAITS_TILE_DATA_STRIDE_SHIFT 48
#define PNANOVTT_TRAITS_TILE_DATA_SHIFT_SHIFT 56
#define PNANOVTT_TRAITS_TILE_DATA_SHIFT_SHIFT_HIGH 24
#define PNANOVTT_TRAITS_TILE_DATA_OFFSET_MASK ((1ull << 48) - 1ul)
#define PNANOVTT_TRAITS_TILE_DATA_SHIFT_MASK \
    kTileDataShiftMask((~((1ull << 48) - 1ul)) << 8ul)
#define PNANOVTT_TRAITS_TILE_DATA_STRIDE_MASK \
    (((~((1ull << 48) - 1ul)) << 8ul) >> 8ul)
#define PNANOVTT_TRAITS_TEXTURE_COORD_SHIFT 16
#define PNANOVTT_TRAITS_TEXTURE_COORD_MASK ((1ull << 16) - 1ul)
#endif

#define PNANOVTT_TABLE_BITS 6
#define PNANOVTT_TABLE_BITS2 12
#define PNANOVTT_TABLE_WIDTH 64
#define PNANOVTT_TABLE_SIZE 262144

struct pnanovtt_traversal_structure_t
{
    pnanovdb_coord_t m_tijk;
    pnanovdb_int32_t m_depth;
};
PNANOVDB_STRUCT_TYPEDEF(pnanovtt_traversal_structure_t)

struct pnanovtt_container_t
{
    pnanovdb_uint64_t m_buffer_size;
    pnanovdb_uint64_t m_container_size;
    PNANOVTT_TILEINDEX m_tile_count_accum[PNANOVTT_TRAITS_MAX_DEPTH + 2];
#if defined(PNANOVDB_C)
    PNANOVTT_TILEINDEX m_table[PNANOVTT_TABLE_SIZE];
#endif
    pnanovdb_coord_t m_table_min;
    pnanovdb_int32_t m_root_depth;
    pnanovdb_uint64_t m_parent_offset;
    //pnanovdb_uint32_t m_tiles[];
};
PNANOVDB_STRUCT_TYPEDEF(pnanovtt_container_t)

struct pnanovtt_container_handle_t { pnanovdb_address_t address; };
PNANOVDB_STRUCT_TYPEDEF(pnanovtt_container_handle_t)

#define PNANOVTT_CONTAINER_OFF_BUFFER_SIZE 0
#define PNANOVTT_CONTAINER_OFF_CONTAINER_SIZE 8
#define PNANOVTT_CONTAINER_OFF_TILE_COUNT_ACCUM 16
#define PNANOVTT_CONTAINER_OFF_TABLE 136
#define PNANOVTT_CONTAINER_OFF_TABLE_MIN 1048712
#define PNANOVTT_CONTAINER_OFF_ROOT_DEPTH 1048724
#define PNANOVTT_CONTAINER_OFF_PARENT_OFFSET 1048728
#define PNANOVTT_CONTAINER_OFF_TILES 1048736

PNANOVDB_FORCE_INLINE pnanovdb_uint64_t pnanovtt_container_get_buffer_size(pnanovdb_buf_t buf, pnanovtt_container_handle_t p) {
    return pnanovdb_read_uint64(buf, pnanovdb_address_offset(p.address, PNANOVTT_CONTAINER_OFF_BUFFER_SIZE));
}
PNANOVDB_FORCE_INLINE pnanovdb_uint64_t pnanovtt_container_get_container_size(pnanovdb_buf_t buf, pnanovtt_container_handle_t p) {
    return pnanovdb_read_uint64(buf, pnanovdb_address_offset(p.address, PNANOVTT_CONTAINER_OFF_CONTAINER_SIZE));
}
PNANOVDB_FORCE_INLINE pnanovdb_uint32_t pnanovtt_container_get_tile_count_accum(pnanovdb_buf_t buf, pnanovtt_container_handle_t p, pnanovdb_uint32_t i) {
    return pnanovdb_read_uint32(buf, pnanovdb_address_offset(p.address, PNANOVTT_CONTAINER_OFF_TILE_COUNT_ACCUM + i * PNANOVTT_TILEINDEX_SIZE));
}
PNANOVDB_FORCE_INLINE pnanovdb_uint32_t pnanovtt_container_get_table_entry(pnanovdb_buf_t buf, pnanovtt_container_handle_t p, PNANOVDB_IN(pnanovdb_coord_t) ijk) {
    return pnanovdb_read_uint32(buf, pnanovdb_address_offset(p.address, PNANOVTT_CONTAINER_OFF_TABLE + ((PNANOVDB_DEREF(ijk).z << PNANOVTT_TABLE_BITS2) + (PNANOVDB_DEREF(ijk).y << PNANOVTT_TABLE_BITS) + PNANOVDB_DEREF(ijk).x) * PNANOVTT_TILEINDEX_SIZE));
}
PNANOVDB_FORCE_INLINE pnanovdb_coord_t pnanovtt_container_get_table_min(pnanovdb_buf_t buf, pnanovtt_container_handle_t p) {
    return pnanovdb_read_coord(buf, pnanovdb_address_offset(p.address, PNANOVTT_CONTAINER_OFF_TABLE_MIN));
}
PNANOVDB_FORCE_INLINE pnanovdb_int32_t pnanovtt_container_get_root_depth(pnanovdb_buf_t buf, pnanovtt_container_handle_t p) {
    return pnanovdb_read_int32(buf, pnanovdb_address_offset(p.address, PNANOVTT_CONTAINER_OFF_ROOT_DEPTH));
}
PNANOVDB_FORCE_INLINE pnanovdb_uint64_t pnanovtt_container_get_parent_offset(pnanovdb_buf_t buf, pnanovtt_container_handle_t p) {
    return pnanovdb_read_uint64(buf, pnanovdb_address_offset(p.address, PNANOVTT_CONTAINER_OFF_PARENT_OFFSET));
}
PNANOVDB_FORCE_INLINE pnanovdb_uint32_t pnanovtt_container_get_active_tile_count(pnanovdb_buf_t buf, pnanovtt_container_handle_t p) {
    return pnanovdb_read_uint32(buf, pnanovdb_address_offset(p.address, PNANOVTT_CONTAINER_OFF_TILE_COUNT_ACCUM + (PNANOVTT_TRAITS_MAX_DEPTH + 1) * PNANOVTT_TILEINDEX_SIZE));
}


#ifdef NANOVTT_SUPPORT_COALESCED_TILES

PNANOVDB_FORCE_INLINE
pnanovdb_int32_t pnanovtt_get_tile_shift(pnanovdb_buf_t buf, pnanovtt_container_handle_t p, PNANOVDB_IN(PNANOVTT_TILEINDEX)
                                             in_index) {
    return 1 + (pnanovdb_uint64_high(pnanovdb_read_uint64(
               buf,
#if defined(PNANOVDB_ADDRESS_32)
               pnanovdb_address_offset(
                   p.address,
                   PNANOVTT_CONTAINER_OFF_TILES + PNANOVDB_DEREF(in_index) * 8
#else
               pnanovdb_address_offset64(
                   pnanovdb_uint32_as_uint64_low(PNANOVTT_CONTAINER_OFF_TILES) +
                   pnanovdb_uint32_as_uint64_low(PNANOVDB_DEREF(in_index)) *
                       pnanovdb_uint32_as_uint64_low(8)
#endif
                   ))) >>
           PNANOVTT_TRAITS_TILE_DATA_SHIFT_SHIFT_HIGH);
}

PNANOVDB_FORCE_INLINE
pnanovdb_int32_t pnanovtt_container_get_parent_and_shift
(pnanovdb_buf_t buf, pnanovtt_container_handle_t p, PNANOVDB_INOUT(PNANOVTT_TILEINDEX) io_index)
{
    PNANOVDB_DEREF(io_index) = pnanovdb_read_uint32(
        buf,
#if defined(PNANOVDB_ADDRESS_32)
        pnanovdb_address_offset(
            p.address,
            pnanovdb_uint64_low(pnanovtt_container_get_parent_offset(buf, p)) +
                PNANOVTT_CONTAINER_OFF_TILES +
                PNANOVDB_DEREF(io_index) * PNANOVTT_TILEINDEX_SIZE));
#else
        pnanovdb_address_offset64(
            p.address,
            pnanovtt_container_get_parent_offset(buf, p) +
                pnanovdb_uint32_as_uint64_low(PNANOVTT_CONTAINER_OFF_TILES) +
                PNANOVDB_DEREF(io_index) *
                    pnanovdb_uint32_as_uint64_low(PNANOVTT_TILEINDEX_SIZE)));
#endif
    return pnanovtt_get_tile_shift(buf, p, io_index);
}

PNANOVDB_FORCE_INLINE
pnanovdb_int32_t pnanovtt_container_get_child_index_and_shift(
        pnanovdb_buf_t              buf,
        pnanovtt_container_handle_t p,
        PNANOVTT_TILEINDEX in_index,
        PNANOVDB_IN(pnanovdb_coord_t) in_diff,
        pnanovdb_int32_t in_vw,
        PNANOVDB_INOUT(pnanovdb_coord_t)
        io_ijk,
        PNANOVDB_INOUT(PNANOVTT_TILEINDEX)
        out_index)
{
    PNANOVDB_DEREF(io_ijk).x = PNANOVDB_DEREF(in_diff).x >> in_vw;
    PNANOVDB_DEREF(io_ijk).y = PNANOVDB_DEREF(in_diff).y >> in_vw;
    PNANOVDB_DEREF(io_ijk).z = PNANOVDB_DEREF(in_diff).z >> in_vw;

#if defined(PNANOVDB_ADDRESS_32)
        pnanovdb_uint64_t i1 = pnanovdb_read_uint64(
            buf, pnanovdb_address_offset(
                p.address,
                PNANOVTT_CONTAINER_OFF_TILES +
                pnanovtt_min_tileindex(
                    in_index,
                    pnanovtt_container_get_tile_count_accum(
                        buf, p, PNANOVTT_TRAITS_MAX_DEPTH)) *
                8));

    PNANOVDB_DEREF(out_index) = pnanovdb_read_uint32(
        buf, pnanovdb_address_offset(
                 p.address, PNANOVTT_CONTAINER_OFF_TILES +
                                (pnanovdb_uint64_low(i1) &
                                 PNANOVTT_TRAITS_TILE_DATA_OFFSET_MASK_LOW) +
                                (PNANOVDB_DEREF(io_ijk).x +
                                 (PNANOVDB_DEREF(io_ijk).y << 1) +
                                 (PNANOVDB_DEREF(io_ijk).z << 2)) *
                                    PNANOVTT_TILEINDEX_SIZE));

    return 1 + (pnanovdb_uint64_high(i1) >> PNANOVTT_TRAITS_TILE_DATA_SHIFT_SHIFT_HIGH);
#else
    pnanovdb_uint64_t i1 = pnanovdb_read_uint64(
        buf,
        pnanovdb_address_offset64(
            p.address,
            pnanovdb_uint32_as_uint64_low(PNANOVTT_CONTAINER_OFF_TILES) +
            pnanovdb_uint32_as_uint64_low(
                pnanovtt_min_tileindex(
                    in_index,
                    pnanovtt_container_get_tile_count_accum(
                        buf, p, PNANOVTT_TRAITS_MAX_DEPTH))) *
            pnanovdb_uint32_as_uint64_low(8)));

    PNANOVDB_DEREF(out_index) = pnanovdb_read_uint32(
        buf, pnanovdb_address_offset64(
            p.address,
            (i1 & PNANOVTT_TRAITS_TILE_DATA_OFFSET_MASK)
            + (PNANOVDB_DEREF(io_ijk).x + (PNANOVDB_DEREF(io_ijk).y << 1) +
                (PNANOVDB_DEREF(io_ijk).z << 2)) *
            PNANOVTT_TILEINDEX_SIZE));

    return 1 + (pnanovdb_uint64(i1) >> PNANOVTT_TRAITS_TILE_DATA_SHIFT_SHIFT);
#endif
}

#else // NANOVTT_SUPPORT_COALESCED_TILES

PNANOVDB_FORCE_INLINE
void pnanovtt_container_get_parent
(pnanovdb_buf_t buf, pnanovtt_container_handle_t p, PNANOVDB_INOUT(PNANOVTT_TILEINDEX) io_index)
{
    PNANOVDB_DEREF(io_index) = pnanovdb_read_uint32(
        buf,
#if defined(PNANOVDB_ADDRESS_32)
        pnanovdb_address_offset(
            p.address,
            pnanovdb_uint64_low(pnanovtt_container_get_parent_offset(buf, p)) +
            PNANOVTT_CONTAINER_OFF_TILES +
            PNANOVDB_DEREF(io_index) * PNANOVTT_TILEINDEX_SIZE));
#else
            pnanovdb_address_offset64(
                p.address,
                pnanovtt_container_get_parent_offset(buf, p) +
                pnanovdb_uint32_as_uint64_low(PNANOVTT_CONTAINER_OFF_TILES) +
                PNANOVDB_DEREF(io_index) *
                pnanovdb_uint32_as_uint64_low(PNANOVTT_TILEINDEX_SIZE)));
#endif
}

PNANOVDB_FORCE_INLINE
void pnanovtt_container_get_child_index(
        pnanovdb_buf_t              buf,
        pnanovtt_container_handle_t p,
        PNANOVTT_TILEINDEX in_index,
        PNANOVDB_IN(pnanovdb_coord_t) in_diff,
        pnanovdb_int32_t in_vw,
        PNANOVDB_INOUT(pnanovdb_coord_t)
        io_ijk,
        PNANOVDB_INOUT(PNANOVTT_TILEINDEX)
        out_index)
{
    PNANOVDB_DEREF(io_ijk).x = PNANOVDB_DEREF(in_diff).x >> in_vw;
    PNANOVDB_DEREF(io_ijk).y = PNANOVDB_DEREF(in_diff).y >> in_vw;
    PNANOVDB_DEREF(io_ijk).z = PNANOVDB_DEREF(in_diff).z >> in_vw;

    PNANOVDB_DEREF(out_index) = pnanovdb_read_uint32(
        buf,
#if defined(PNANOVDB_ADDRESS_32)
        pnanovdb_address_offset(
            p.address,
            PNANOVTT_CONTAINER_OFF_TILES +
            pnanovdb_uint64_low(pnanovdb_read_uint64(
                buf, pnanovdb_address_offset(
                    p.address,
                    PNANOVTT_CONTAINER_OFF_TILES +
                    pnanovtt_min_tileindex(
                        in_index,
                        pnanovtt_container_get_tile_count_accum(
                            buf, p, PNANOVTT_TRAITS_MAX_DEPTH)) *
                    8)))
#else
            pnanovdb_address_offset64(
                p.address,
                PNANOVTT_CONTAINER_OFF_TILES +
                pnanovdb_read_uint64(
                buf,
                pnanovdb_address_offset64(
                    p.address,
                    pnanovdb_uint32_as_uint64(PNANOVTT_CONTAINER_OFF_TILES,
                        0) +
                    pnanovdb_uint32_as_uint64(
                        pnanovtt_min_tileindex(
                            in_index,
                            pnanovtt_container_get_tile_count_accum(
                                buf, p, PNANOVTT_TRAITS_MAX_DEPTH)),
                        0) *
                    pnanovdb_uint32_as_uint64(8, 0)))
#endif
            + (PNANOVDB_DEREF(io_ijk).x + (PNANOVDB_DEREF(io_ijk).y << 1) +
                (PNANOVDB_DEREF(io_ijk).z << 2)) *
            PNANOVTT_TILEINDEX_SIZE));
}

#endif // NANOVTT_SUPPORT_COALESCED_TILES

PNANOVDB_FORCE_INLINE PNANOVTT_TILEINDEX
pnanovtt_container_get_root(pnanovdb_buf_t              buf,
    pnanovtt_container_handle_t p,
    PNANOVDB_IN(pnanovdb_coord_t) in_key,
    PNANOVDB_INOUT(pnanovdb_coord_t) io_tijk)
{
    // Clamp in_key to tile space at the root depth.
    pnanovdb_coord_t table_min = pnanovtt_container_get_table_min(buf, p);
    pnanovdb_int32_t shift = PNANOVTT_TRAITS_MAX_DEPTH -
        pnanovtt_container_get_root_depth(buf, p) + 2;
    pnanovdb_coord_t   ijk;
    ijk.x = (PNANOVDB_DEREF(in_key).x - table_min.x) >> shift;
    ijk.y = (PNANOVDB_DEREF(in_key).y - table_min.y) >> shift;
    ijk.z = (PNANOVDB_DEREF(in_key).z - table_min.z) >> shift;
    PNANOVTT_TILEINDEX hash = ijk.x | ijk.y | ijk.z;
    PNANOVDB_DEREF(io_tijk).x = table_min.x + (ijk.x << shift);
    PNANOVDB_DEREF(io_tijk).y = table_min.y + (ijk.y << shift);
    PNANOVDB_DEREF(io_tijk).z = table_min.z + (ijk.z << shift);
    return (hash < PNANOVTT_TABLE_WIDTH
        ? pnanovtt_container_get_table_entry(buf, p, PNANOVDB_REF(ijk))
        : PNANOVTT_INVALID_TILEINDEX);
}

struct pnanovtt_tree_t
{
    pnanovdb_coord_t bbox_min;
    pnanovdb_coord_t bbox_max;
    pnanovdb_uint64_t m_container_offset;
    pnanovdb_uint64_t m_tree_size;
    //ValueType m_default_value;

};
PNANOVDB_STRUCT_TYPEDEF(pnanovtt_tree_t)
struct pnanovtt_tree_handle_t { pnanovdb_address_t address; };
PNANOVDB_STRUCT_TYPEDEF(pnanovtt_tree_handle_t)

#define PNANOVTT_TREE_OFF_BBOX_MIN 0
#define PNANOVTT_TREE_OFF_BBOX_MAX 12
#define PNANOVTT_TREE_OFF_CONTAINER_OFFSET 24
#define PNANOVTT_TREE_OFF_TREE_SIZE 32
#define PNANOVTT_TREE_OFF_DEFAULT_VALUE 40

PNANOVDB_FORCE_INLINE pnanovdb_coord_t pnanovtt_tree_get_bbox_min(pnanovdb_buf_t buf, pnanovtt_tree_handle_t p) {
    return pnanovdb_read_coord(buf, pnanovdb_address_offset(p.address, PNANOVTT_TREE_OFF_BBOX_MIN));
}
PNANOVDB_FORCE_INLINE pnanovdb_coord_t pnanovtt_tree_get_bbox_max(pnanovdb_buf_t buf, pnanovtt_tree_handle_t p) {
    return pnanovdb_read_coord(buf, pnanovdb_address_offset(p.address, PNANOVTT_TREE_OFF_BBOX_MAX));
}
PNANOVDB_FORCE_INLINE pnanovdb_uint64_t pnanovtt_tree_get_container_offset(pnanovdb_buf_t buf, pnanovtt_tree_handle_t p) {
    return pnanovdb_read_uint64(buf, pnanovdb_address_offset(p.address, PNANOVTT_TREE_OFF_CONTAINER_OFFSET));
}
PNANOVDB_FORCE_INLINE pnanovdb_uint64_t pnanovtt_tree_get_tree_size(pnanovdb_buf_t buf, pnanovtt_tree_handle_t p) {
    return pnanovdb_read_uint64(buf, pnanovdb_address_offset(p.address, PNANOVTT_TREE_OFF_TREE_SIZE));
}

struct pnanovtt_grid_type_constants_t
{
    pnanovdb_uint32_t tree_off_data;
    pnanovdb_uint32_t value_stride_bytes;
};
PNANOVDB_STRUCT_TYPEDEF(pnanovtt_grid_type_constants_t)

#if defined(PNANOVDB_GLSL_OSX)
// OSX doesn't support assignment as above
pnanovtt_grid_type_constants_t  pnanovtt_grid_type_constants[PNANOVDB_GRID_TYPE_END];

PNANOVDB_FORCE_INLINE void pnanovtt_grid_type_constants_init() {

    //    {64, 0},
    pnanovtt_grid_type_constants[PNANOVDB_GRID_TYPE_UNKNOWN].tree_off_data = 64;
    pnanovtt_grid_type_constants[PNANOVDB_GRID_TYPE_UNKNOWN].value_stride_bytes = 0;

    //    {64, 4},
    pnanovtt_grid_type_constants[PNANOVDB_GRID_TYPE_FLOAT].tree_off_data = 64;
    pnanovtt_grid_type_constants[PNANOVDB_GRID_TYPE_FLOAT].value_stride_bytes = 4;

    //    {64, 8},
    pnanovtt_grid_type_constants[PNANOVDB_GRID_TYPE_DOUBLE].tree_off_data = 64;
    pnanovtt_grid_type_constants[PNANOVDB_GRID_TYPE_DOUBLE].value_stride_bytes = 8;

    //    {64, 2},
    pnanovtt_grid_type_constants[PNANOVDB_GRID_TYPE_INT16].tree_off_data = 64;
    pnanovtt_grid_type_constants[PNANOVDB_GRID_TYPE_INT16].value_stride_bytes = 2;

    //    {64, 4},
    pnanovtt_grid_type_constants[PNANOVDB_GRID_TYPE_INT32].tree_off_data = 64;
    pnanovtt_grid_type_constants[PNANOVDB_GRID_TYPE_INT32].value_stride_bytes = 4;

    //    {64, 8},
    pnanovtt_grid_type_constants[PNANOVDB_GRID_TYPE_INT64].tree_off_data = 64;
    pnanovtt_grid_type_constants[PNANOVDB_GRID_TYPE_INT64].value_stride_bytes = 8;

    //    {64, 12},
    pnanovtt_grid_type_constants[PNANOVDB_GRID_TYPE_VEC3F].tree_off_data = 64;
    pnanovtt_grid_type_constants[PNANOVDB_GRID_TYPE_VEC3F].value_stride_bytes = 12;

    //    {64, 24},
    pnanovtt_grid_type_constants[PNANOVDB_GRID_TYPE_VEC3D].tree_off_data = 64;
    pnanovtt_grid_type_constants[PNANOVDB_GRID_TYPE_VEC3D].value_stride_bytes = 24;

    //    {64, 0},
    pnanovtt_grid_type_constants[PNANOVDB_GRID_TYPE_MASK].tree_off_data = 64;
    pnanovtt_grid_type_constants[PNANOVDB_GRID_TYPE_MASK].value_stride_bytes = 0;

    //    {64, 2},
    pnanovtt_grid_type_constants[PNANOVDB_GRID_TYPE_HALF].tree_off_data = 64;
    pnanovtt_grid_type_constants[PNANOVDB_GRID_TYPE_HALF].value_stride_bytes = 2;

    //    {64, 4},
    pnanovtt_grid_type_constants[PNANOVDB_GRID_TYPE_UINT32].tree_off_data = 64;
    pnanovtt_grid_type_constants[PNANOVDB_GRID_TYPE_UINT32].value_stride_bytes = 4;

    //    {64, 1},
    pnanovtt_grid_type_constants[PNANOVDB_GRID_TYPE_BOOLEAN].tree_off_data = 64;
    pnanovtt_grid_type_constants[PNANOVDB_GRID_TYPE_BOOLEAN].value_stride_bytes = 1;

    //    {64, 4},
    pnanovtt_grid_type_constants[PNANOVDB_GRID_TYPE_RGBA8].tree_off_data = 64;
    pnanovtt_grid_type_constants[PNANOVDB_GRID_TYPE_RGBA8].value_stride_bytes = 4;

    //    {64, 0},
    pnanovtt_grid_type_constants[PNANOVDB_GRID_TYPE_FP4].tree_off_data = 64;
    pnanovtt_grid_type_constants[PNANOVDB_GRID_TYPE_FP4].value_stride_bytes = 0;

    //    {64, 0},
    pnanovtt_grid_type_constants[PNANOVDB_GRID_TYPE_FP8].tree_off_data = 64;
    pnanovtt_grid_type_constants[PNANOVDB_GRID_TYPE_FP8].value_stride_bytes = 0;

    //    {64, 0},
    pnanovtt_grid_type_constants[PNANOVDB_GRID_TYPE_FP16].tree_off_data = 64;
    pnanovtt_grid_type_constants[PNANOVDB_GRID_TYPE_FP16].value_stride_bytes = 0;

    //    {64, 0},
    pnanovtt_grid_type_constants[PNANOVDB_GRID_TYPE_FPN].tree_off_data = 64;
    pnanovtt_grid_type_constants[PNANOVDB_GRID_TYPE_FPN].value_stride_bytes = 0;

    //    {64, 16},
    pnanovtt_grid_type_constants[PNANOVDB_GRID_TYPE_VEC4F].tree_off_data = 64;
    pnanovtt_grid_type_constants[PNANOVDB_GRID_TYPE_VEC4F].value_stride_bytes = 16;

    //    {96, 32}
    pnanovtt_grid_type_constants[PNANOVDB_GRID_TYPE_VEC4D].tree_off_data = 96;
    pnanovtt_grid_type_constants[PNANOVDB_GRID_TYPE_VEC4D].value_stride_bytes = 32;

}
#else
PNANOVDB_STATIC_CONST pnanovtt_grid_type_constants_t
pnanovtt_grid_type_constants[PNANOVDB_GRID_TYPE_END] =
{
    {64, 0},
    {64, 4},
    {64, 8},
    {64, 2},
    {64, 4},
    {64, 8},
    {64, 12},
    {64, 24},
    {64, 0},
    {64, 2},
    {64, 4},
    {64, 1},
    {64, 4},
    {64, 0},
    {64, 0},
    {64, 0},
    {64, 0},
    {64, 16},
    {96, 32}
};
PNANOVDB_FORCE_INLINE void pnanovtt_grid_type_constants_init();
PNANOVDB_FORCE_INLINE void pnanovtt_grid_type_constants_init() {}
#endif

PNANOVDB_FORCE_INLINE pnanovdb_address_t pnanovtt_tree_get_tile_data_address(pnanovdb_buf_t buf, pnanovtt_tree_handle_t p, pnanovdb_grid_type_t grid_type) {
    pnanovtt_container_handle_t container_handle;
    container_handle.address = pnanovdb_address_null();
    return pnanovdb_address_offset(
        p.address, PNANOVTT_GRID_TYPE_GET(grid_type, tree_off_data) +
        8 * (1 + pnanovtt_container_get_active_tile_count(
            buf, container_handle)));
}

PNANOVDB_FORCE_INLINE void pnanovtt_tree_get_tile_texture_coord(
    pnanovdb_grid_type_t grid_type,
    pnanovdb_buf_t         buf,
    pnanovtt_tree_handle_t p,
    PNANOVTT_TILEINDEX in_tile,
    PNANOVDB_INOUT(pnanovdb_coord_t) io_coord)
{
#ifdef PNANOVDB_ADDRESS_32
    pnanovdb_uint64_t v = pnanovdb_read_uint64(
        buf, pnanovdb_address_offset64(
            p.address, PNANOVTT_GRID_TYPE_GET(grid_type, tree_off_data) +
            pnanovdb_uint32_as_uint64(in_tile * 8, 0)));

    pnanovdb_uint32_t offset_low = pnanovdb_uint64_low(v) & PNANOVTT_TRAITS_TILE_DATA_OFFSET_MASK_LOW;
    pnanovdb_uint32_t offset_high = pnanovdb_uint64_high(v) & PNANOVTT_TRAITS_TILE_DATA_OFFSET_MASK_HIGH;

    PNANOVDB_DEREF(io_coord).x = pnanovdb_uint32_as_int32(offset_low & PNANOVTT_TRAITS_TEXTURE_COORD_MASK);
    offset_low = offset_low >> PNANOVTT_TRAITS_TEXTURE_COORD_SHIFT;
    PNANOVDB_DEREF(io_coord).y = pnanovdb_uint32_as_int32(offset_low & PNANOVTT_TRAITS_TEXTURE_COORD_MASK);
    PNANOVDB_DEREF(io_coord).z = pnanovdb_uint32_as_int32(offset_high & PNANOVTT_TRAITS_TEXTURE_COORD_MASK);
#else
    pnanovdb_uint64_t v = pnanovdb_read_uint64(
        buf, pnanovdb_address_offset64(
            p.address, PNANOVTT_GRID_TYPE_GET(grid_type, tree_off_data) +
            in_tile * 8));

    pnanovdb_uint64_t offset = v & PNANOVTT_TRAITS_TILE_DATA_OFFSET_MASK;

    PNANOVDB_DEREF(io_coord).x = pnanovdb_uint32_as_int32(pnanovdb_uint64_low(offset & PNANOVTT_TRAITS_TEXTURE_COORD_MASK));
    offset = offset >> PNANOVTT_TRAITS_TEXTURE_COORD_SHIFT;
    PNANOVDB_DEREF(io_coord).y = pnanovdb_uint32_as_int32(pnanovdb_uint64_low(offset & PNANOVTT_TRAITS_TEXTURE_COORD_MASK));
    offset = offset >> PNANOVTT_TRAITS_TEXTURE_COORD_SHIFT;
    PNANOVDB_DEREF(io_coord).z = pnanovdb_uint32_as_int32(pnanovdb_uint64_low(offset & PNANOVTT_TRAITS_TEXTURE_COORD_MASK));
#endif
}

PNANOVDB_FORCE_INLINE pnanovdb_address_t pnanovtt_tree_get_value_address_at_i(
    pnanovdb_grid_type_t grid_type,
    pnanovdb_buf_t         buf,
    pnanovtt_tree_handle_t p,
    PNANOVTT_TILEINDEX in_tile,
    pnanovdb_int32_t in_i)
{
#ifdef PNANOVDB_ADDRESS_32
    pnanovdb_uint64_t v = pnanovdb_read_uint64(
        buf, pnanovdb_address_offset64(
            p.address, PNANOVTT_GRID_TYPE_GET(grid_type, tree_off_data) +
            pnanovdb_uint32_as_uint64(in_tile * 8, 0)));

    pnanovdb_uint32_t offset_low = pnanovdb_uint64_low(v) & PNANOVTT_TRAITS_TILE_DATA_OFFSET_MASK_LOW;
    pnanovdb_uint32_t offset_high = pnanovdb_uint64_high(v) & PNANOVTT_TRAITS_TILE_DATA_OFFSET_MASK_HIGH;
    pnanovdb_uint64_t offset = pnanovdb_uint32_as_uint64(offset_low, offset_high);

    return pnanovdb_address_offset64(
        p.address,
        PNANOVTT_GRID_TYPE_GET(grid_type, tree_off_data) + offset +
        pnanovdb_uint32_as_uint64(pnanovdb_int32_as_uint32(in_i) *
            PNANOVTT_GRID_TYPE_GET(grid_type, value_stride_bytes), 0));
#else
    pnanovdb_uint64_t v = pnanovdb_read_uint64(
        buf, pnanovdb_address_offset64(
            p.address, PNANOVTT_GRID_TYPE_GET(grid_type, tree_off_data) +
            in_tile * 8));

    pnanovdb_uint64_t offset = v & PNANOVTT_TRAITS_TILE_DATA_OFFSET_MASK;

    return pnanovdb_address_offset64(
        p.address,
        PNANOVTT_GRID_TYPE_GET(grid_type, tree_off_data) + offset +
        in_i * PNANOVTT_GRID_TYPE_GET(grid_type, value_stride_bytes));
#endif
}

PNANOVDB_FORCE_INLINE pnanovdb_address_t pnanovtt_tree_get_value_address_at_ijk(
    pnanovdb_grid_type_t grid_type,
    pnanovdb_buf_t         buf,
    pnanovtt_tree_handle_t p,
    PNANOVTT_TILEINDEX in_tile,
    PNANOVDB_IN(pnanovdb_coord_t) in_ijk)
{
#ifdef PNANOVDB_ADDRESS_32
    pnanovdb_uint64_t v = pnanovdb_read_uint64(
        buf, pnanovdb_address_offset64(
            p.address, PNANOVTT_GRID_TYPE_GET(grid_type, tree_off_data) +
            pnanovdb_uint32_as_uint64(in_tile * 8, 0)));

    pnanovdb_uint32_t offset_low = pnanovdb_uint64_low(v) & PNANOVTT_TRAITS_TILE_DATA_OFFSET_MASK_LOW;
    pnanovdb_uint32_t offset_high = pnanovdb_uint64_high(v) & PNANOVTT_TRAITS_TILE_DATA_OFFSET_MASK_HIGH;
    pnanovdb_uint64_t offset = pnanovdb_uint32_as_uint64(offset_low, offset_high);

    pnanovdb_int32_t tile_data_stride =
        (pnanovdb_uint32_as_int32(pnanovdb_uint64_high(v)) &
            PNANOVTT_TRAITS_TILE_DATA_STRIDE_MASK_HIGH) >>
        PNANOVTT_TRAITS_TILE_DATA_STRIDE_SHIFT_HIGH;

    return pnanovdb_address_offset64(
        p.address,
        PNANOVTT_GRID_TYPE_GET(grid_type, tree_off_data) + offset +
        pnanovdb_uint32_as_uint64((PNANOVDB_DEREF(in_ijk).x +
            PNANOVDB_DEREF(in_ijk).y * tile_data_stride +
            PNANOVDB_DEREF(in_ijk).z * tile_data_stride * tile_data_stride) *
            PNANOVTT_GRID_TYPE_GET(grid_type, value_stride_bytes), 0));
#else
    pnanovdb_uint64_t v = pnanovdb_read_uint64(
        buf, pnanovdb_address_offset64(
            p.address, PNANOVTT_GRID_TYPE_GET(grid_type, tree_off_data) +
            in_tile * 8));

    pnanovdb_uint64_t offset = v & PNANOVTT_TRAITS_TILE_DATA_OFFSET_MASK;

    pnanovdb_int32_t  tile_data_stride =
        (v & PNANOVTT_TRAITS_TILE_DATA_STRIDE_MASK) >>
        PNANOVTT_TRAITS_TILE_DATA_STRIDE_SHIFT;

    return pnanovdb_address_offset64(
        p.address,
        PNANOVTT_GRID_TYPE_GET(grid_type, tree_off_data) + offset +
        (PNANOVDB_DEREF(in_ijk).x +
            PNANOVDB_DEREF(in_ijk).y * tile_data_stride +
            PNANOVDB_DEREF(in_ijk).z * tile_data_stride * tile_data_stride) *
        PNANOVTT_GRID_TYPE_GET(grid_type, value_stride_bytes));
#endif
}

PNANOVDB_FORCE_INLINE pnanovdb_address_t pnanovtt_tree_get_value_address(
    pnanovdb_grid_type_t grid_type,
    pnanovdb_buf_t         buf,
    pnanovtt_tree_handle_t p,
    PNANOVTT_TILEINDEX in_tile,
    PNANOVDB_IN(pnanovdb_coord_t) in_ijk,
    PNANOVDB_INOUT(pnanovdb_int32_t) io_shift,
    PNANOVDB_INOUT(pnanovdb_int32_t) io_stride,
    PNANOVDB_INOUT(pnanovdb_coord_t) io_ijk)
{
#ifdef PNANOVDB_ADDRESS_32
    pnanovdb_uint64_t v = pnanovdb_read_uint64(
        buf, pnanovdb_address_offset64(
            p.address, PNANOVTT_GRID_TYPE_GET(grid_type, tree_off_data) +
            pnanovdb_uint32_as_uint64(in_tile * 8, 0)));

    pnanovdb_uint32_t offset_low = pnanovdb_uint64_low(v) & PNANOVTT_TRAITS_TILE_DATA_OFFSET_MASK_LOW;
    pnanovdb_uint32_t offset_high = pnanovdb_uint64_high(v) & PNANOVTT_TRAITS_TILE_DATA_OFFSET_MASK_HIGH;
    pnanovdb_uint64_t offset = pnanovdb_uint32_as_uint64(offset_low, offset_high);

    pnanovdb_int32_t tile_data_stride =
        (pnanovdb_uint32_as_int32(pnanovdb_uint64_high(v)) &
            PNANOVTT_TRAITS_TILE_DATA_STRIDE_MASK_HIGH) >>
        PNANOVTT_TRAITS_TILE_DATA_STRIDE_SHIFT_HIGH;

    pnanovdb_int32_t tile_data_shift = pnanovdb_uint32_as_int32(pnanovdb_uint64_high(v)) >> PNANOVTT_TRAITS_TILE_DATA_SHIFT_SHIFT_HIGH;

    tile_data_shift = PNANOVDB_DEREF(io_shift) - tile_data_shift;
    PNANOVDB_DEREF(io_ijk).x = (PNANOVDB_DEREF(in_ijk).x >> tile_data_shift);
    PNANOVDB_DEREF(io_ijk).y = (PNANOVDB_DEREF(in_ijk).y >> tile_data_shift);
    PNANOVDB_DEREF(io_ijk).z = (PNANOVDB_DEREF(in_ijk).z >> tile_data_shift);
    PNANOVDB_DEREF(io_shift) = tile_data_shift;
    PNANOVDB_DEREF(io_stride) = tile_data_stride;

    return pnanovdb_address_offset64(
        p.address,
        PNANOVTT_GRID_TYPE_GET(grid_type, tree_off_data) + offset +
        pnanovdb_uint32_as_uint64((PNANOVDB_DEREF(io_ijk).x +
            PNANOVDB_DEREF(io_ijk).y * tile_data_stride +
            PNANOVDB_DEREF(io_ijk).z * tile_data_stride * tile_data_stride) *
            PNANOVTT_GRID_TYPE_GET(grid_type, value_stride_bytes), 0));
#else
    pnanovdb_uint64_t v = pnanovdb_read_uint64(
        buf, pnanovdb_address_offset64(
            p.address, PNANOVTT_GRID_TYPE_GET(grid_type, tree_off_data) +
            in_tile * 8));

    pnanovdb_uint64_t offset = v & PNANOVTT_TRAITS_TILE_DATA_OFFSET_MASK;

    pnanovdb_int32_t  tile_data_stride =
        (v & PNANOVTT_TRAITS_TILE_DATA_STRIDE_MASK) >>
        PNANOVTT_TRAITS_TILE_DATA_STRIDE_SHIFT;

    pnanovdb_uint32_t tile_data_shift = v >> PNANOVTT_TRAITS_TILE_DATA_SHIFT_SHIFT;

    tile_data_shift = PNANOVDB_DEREF(io_shift) - tile_data_shift;
    PNANOVDB_DEREF(io_ijk).x = (PNANOVDB_DEREF(in_ijk).x >> tile_data_shift);
    PNANOVDB_DEREF(io_ijk).y = (PNANOVDB_DEREF(in_ijk).y >> tile_data_shift);
    PNANOVDB_DEREF(io_ijk).z = (PNANOVDB_DEREF(in_ijk).z >> tile_data_shift);
    PNANOVDB_DEREF(io_shift) = tile_data_shift;
    PNANOVDB_DEREF(io_stride) = tile_data_stride;

    return pnanovdb_address_offset64(
        p.address,
        PNANOVTT_GRID_TYPE_GET(grid_type, tree_off_data) + offset +
        (PNANOVDB_DEREF(io_ijk).x +
            PNANOVDB_DEREF(io_ijk).y * tile_data_stride +
            PNANOVDB_DEREF(io_ijk).z * tile_data_stride * tile_data_stride) *
        PNANOVTT_GRID_TYPE_GET(grid_type, value_stride_bytes));
#endif
}


PNANOVDB_FORCE_INLINE void pnanovtt_tree_get_value_texture_coord(
    pnanovdb_grid_type_t grid_type,
    pnanovdb_buf_t         buf,
    pnanovtt_tree_handle_t p,
    PNANOVTT_TILEINDEX in_tile,
    PNANOVDB_IN(pnanovdb_coord_t) in_ijk,
    PNANOVDB_INOUT(pnanovdb_int32_t) io_shift,
    PNANOVDB_INOUT(pnanovdb_coord_t) io_ijk,
    PNANOVDB_INOUT(pnanovdb_coord_t) io_coord)
{
#ifdef PNANOVDB_ADDRESS_32
    pnanovdb_uint64_t v = pnanovdb_read_uint64(
        buf, pnanovdb_address_offset64(
            p.address, PNANOVTT_GRID_TYPE_GET(grid_type, tree_off_data) +
            pnanovdb_uint32_as_uint64(in_tile * 8, 0)));

    pnanovdb_uint32_t offset_low = pnanovdb_uint64_low(v) & PNANOVTT_TRAITS_TILE_DATA_OFFSET_MASK_LOW;
    pnanovdb_uint32_t offset_high = pnanovdb_uint64_high(v) & PNANOVTT_TRAITS_TILE_DATA_OFFSET_MASK_HIGH;
    PNANOVDB_DEREF(io_shift) = PNANOVDB_DEREF(io_shift) - (pnanovdb_uint32_as_int32(pnanovdb_uint64_high(v)) >> PNANOVTT_TRAITS_TILE_DATA_SHIFT_SHIFT_HIGH);

    PNANOVDB_DEREF(io_ijk).x = (PNANOVDB_DEREF(in_ijk).x >> PNANOVDB_DEREF(io_shift));
    PNANOVDB_DEREF(io_ijk).y = (PNANOVDB_DEREF(in_ijk).y >> PNANOVDB_DEREF(io_shift));
    PNANOVDB_DEREF(io_ijk).z = (PNANOVDB_DEREF(in_ijk).z >> PNANOVDB_DEREF(io_shift));

    PNANOVDB_DEREF(io_coord).x = pnanovdb_uint32_as_int32(offset_low & PNANOVTT_TRAITS_TEXTURE_COORD_MASK) + PNANOVDB_DEREF(io_ijk).x;
    offset_low = offset_low >> PNANOVTT_TRAITS_TEXTURE_COORD_SHIFT;
    PNANOVDB_DEREF(io_coord).y = pnanovdb_uint32_as_int32(offset_low & PNANOVTT_TRAITS_TEXTURE_COORD_MASK) + PNANOVDB_DEREF(io_ijk).y;
    PNANOVDB_DEREF(io_coord).z = pnanovdb_uint32_as_int32(offset_high & PNANOVTT_TRAITS_TEXTURE_COORD_MASK) + PNANOVDB_DEREF(io_ijk).z;
#else
    pnanovdb_uint64_t v = pnanovdb_read_uint64(
        buf, pnanovdb_address_offset64(
            p.address, PNANOVTT_GRID_TYPE_GET(grid_type, tree_off_data) +
            in_tile * 8));

    pnanovdb_uint64_t offset = v & PNANOVTT_TRAITS_TILE_DATA_OFFSET_MASK;
    PNANOVDB_DEREF(io_shift) = PNANOVDB_DEREF(io_shift) - (v >> PNANOVTT_TRAITS_TILE_DATA_SHIFT_SHIFT);

    PNANOVDB_DEREF(io_ijk).x = (PNANOVDB_DEREF(in_ijk).x >> PNANOVDB_DEREF(io_shift));
    PNANOVDB_DEREF(io_ijk).y = (PNANOVDB_DEREF(in_ijk).y >> PNANOVDB_DEREF(io_shift));
    PNANOVDB_DEREF(io_ijk).z = (PNANOVDB_DEREF(in_ijk).z >> PNANOVDB_DEREF(io_shift));

    PNANOVDB_DEREF(io_coord).x = pnanovdb_uint32_as_int32(pnanovdb_uint64_low(offset & PNANOVTT_TRAITS_TEXTURE_COORD_MASK)) + PNANOVDB_DEREF(io_ijk).x;
    offset = offset >> PNANOVTT_TRAITS_TEXTURE_COORD_SHIFT;
    PNANOVDB_DEREF(io_coord).y = pnanovdb_uint32_as_int32(pnanovdb_uint64_low(offset & PNANOVTT_TRAITS_TEXTURE_COORD_MASK)) + PNANOVDB_DEREF(io_ijk).y;
    offset = offset >> PNANOVTT_TRAITS_TEXTURE_COORD_SHIFT;
    PNANOVDB_DEREF(io_coord).z = pnanovdb_uint32_as_int32(pnanovdb_uint64_low(offset & PNANOVTT_TRAITS_TEXTURE_COORD_MASK)) + PNANOVDB_DEREF(io_ijk).z;
#endif
}


PNANOVDB_FORCE_INLINE pnanovdb_address_t pnanovtt_tree_get_value_address2(
    pnanovdb_grid_type_t grid_type,
    pnanovdb_buf_t         buf,
    pnanovtt_tree_handle_t p,
    PNANOVTT_TILEINDEX in_tile,
    PNANOVDB_IN(pnanovdb_coord_t) in_ijk)
{
#ifdef PNANOVDB_ADDRESS_32
    pnanovdb_uint64_t v = pnanovdb_read_uint64(
        buf, pnanovdb_address_offset(
            p.address, PNANOVTT_GRID_TYPE_GET(grid_type, tree_off_data) +
            in_tile * 8));

    pnanovdb_uint32_t offset_low = pnanovdb_uint64_low(v) & PNANOVTT_TRAITS_TILE_DATA_OFFSET_MASK_LOW;
    pnanovdb_uint32_t offset_high = pnanovdb_uint64_high(v) & PNANOVTT_TRAITS_TILE_DATA_OFFSET_MASK_HIGH;
    pnanovdb_uint64_t offset = pnanovdb_uint32_as_uint64(offset_low, offset_high);

    pnanovdb_uint32_t tile_data_stride =
        (pnanovdb_uint64_high(v) &
            pnanovdb_int32_as_uint32(PNANOVTT_TRAITS_TILE_DATA_STRIDE_MASK_HIGH)) >>
        PNANOVTT_TRAITS_TILE_DATA_STRIDE_SHIFT_HIGH;

    return pnanovdb_address_offset64(
        p.address,
        PNANOVTT_GRID_TYPE_GET(grid_type, tree_off_data) + offset +
        pnanovdb_uint32_as_uint64((PNANOVDB_DEREF(in_ijk).x +
            PNANOVDB_DEREF(in_ijk).y * tile_data_stride +
            PNANOVDB_DEREF(in_ijk).z * tile_data_stride * tile_data_stride) *
            PNANOVTT_GRID_TYPE_GET(grid_type, value_stride_bytes), 0));
#else
    pnanovdb_uint64_t v = pnanovdb_read_uint64(
        buf, pnanovdb_address_offset64(
            p.address,
            pnanovdb_uint32_as_uint64(
                PNANOVTT_GRID_TYPE_GET(grid_type, tree_off_data), 0) +
            pnanovdb_uint32_as_uint64(in_tile, 0) * 8));

    pnanovdb_uint64_t offset = v & PNANOVTT_TRAITS_TILE_DATA_OFFSET_MASK;

    pnanovdb_int32_t  tile_data_stride =
        (v & PNANOVTT_TRAITS_TILE_DATA_STRIDE_MASK) >>
        PNANOVTT_TRAITS_TILE_DATA_STRIDE_SHIFT;

    return pnanovdb_address_offset64(
        p.address,
        PNANOVTT_GRID_TYPE_GET(grid_type, tree_off_data) + offset +
        (PNANOVDB_DEREF(in_ijk).x +
            PNANOVDB_DEREF(in_ijk).y * tile_data_stride +
            PNANOVDB_DEREF(in_ijk).z * tile_data_stride * tile_data_stride) *
        PNANOVTT_GRID_TYPE_GET(grid_type, value_stride_bytes));
#endif
}

// ------------------------------------------------ ReadAccessor -----------------------------------------------------------

struct pnanovtt_readaccessor_t
{
    pnanovdb_int32_t m_depth;
    pnanovdb_int32_t m_root_depth;
    pnanovdb_int32_t m_vw;
    pnanovdb_int32_t m_tile_data_stride;
    PNANOVTT_TILEINDEX m_tile;
    pnanovdb_coord_t m_ijk;
    pnanovdb_coord_t m_tijk;
    float m_voxel_size;
    pnanovdb_bool_t m_found;
    pnanovdb_bool_t m_is_level_set;
    pnanovtt_tree_handle_t m_tree;
    pnanovtt_container_handle_t m_container;
};
PNANOVDB_STRUCT_TYPEDEF(pnanovtt_readaccessor_t)

PNANOVDB_FORCE_INLINE void pnanovtt_readaccessor_init_v2(pnanovdb_buf_t buf, pnanovdb_buf_t grid_buf, PNANOVDB_INOUT(pnanovtt_readaccessor_t) acc, pnanovdb_grid_handle_t grid, pnanovtt_tree_handle_t tree)
{
    PNANOVDB_DEREF(acc).m_tree.address = tree.address;
    // Note: the container is always stored first in the buffer
    PNANOVDB_DEREF(acc).m_container.address = pnanovdb_address_null();
#if !defined(PNANOVDB_GLSL_OSX) // Pass voxel size by argument since packDouble2x32 doesn't work on Apple arm64
    PNANOVDB_DEREF(acc).m_voxel_size = pnanovdb_double_as_float(pnanovdb_grid_get_voxel_size(grid_buf, grid, 0));
#endif
    pnanovdb_uint32_t grid_class = pnanovdb_grid_get_grid_class(grid_buf, grid);
    PNANOVDB_DEREF(acc).m_is_level_set = (grid_class == PNANOVDB_GRID_CLASS_LEVEL_SET);
    PNANOVDB_DEREF(acc).m_root_depth = PNANOVDB_DEREF(acc).m_depth = pnanovtt_container_get_root_depth(buf, PNANOVDB_DEREF(acc).m_container);
    PNANOVDB_DEREF(acc).m_vw = PNANOVTT_TRAITS_MAX_DEPTH - PNANOVDB_DEREF(acc).m_depth;
    PNANOVDB_DEREF(acc).m_tijk.x = PNANOVTT_TRAITS_TILE_DOMAIN_WIDTH;
    PNANOVDB_DEREF(acc).m_tijk.y = PNANOVTT_TRAITS_TILE_DOMAIN_WIDTH;
    PNANOVDB_DEREF(acc).m_tijk.z = PNANOVTT_TRAITS_TILE_DOMAIN_WIDTH;
    PNANOVDB_DEREF(acc).m_found = PNANOVDB_FALSE;
    PNANOVDB_DEREF(acc).m_tile = 0;
    PNANOVDB_DEREF(acc).m_tile_data_stride = PNANOVTT_TRAITS_TILE_DATA_WIDTH;
    PNANOVDB_DEREF(acc).m_ijk.x = 0;
    PNANOVDB_DEREF(acc).m_ijk.y = 0;
    PNANOVDB_DEREF(acc).m_ijk.z = 0;
}

PNANOVDB_FORCE_INLINE void pnanovtt_readaccessor_init(pnanovdb_buf_t buf, PNANOVDB_INOUT(pnanovtt_readaccessor_t) acc, pnanovdb_grid_handle_t grid, pnanovtt_tree_handle_t tree)
{
    pnanovtt_readaccessor_init_v2(buf, buf, acc, grid, tree);
}

PNANOVDB_FORCE_INLINE pnanovdb_address_t pnanovtt_readaccessor_get_value_address_at_depth_v2(pnanovdb_grid_type_t grid_type, pnanovdb_buf_t buf, pnanovdb_buf_t grid_buf, PNANOVDB_INOUT(pnanovtt_readaccessor_t) acc, PNANOVDB_IN(pnanovdb_coord_t) in_ijk, pnanovdb_int32_t in_max_depth)
{
    // This version traverses bottom up from the last tile
    PNANOVDB_DEREF(acc).m_found = PNANOVDB_TRUE;

    pnanovdb_coord_t diff;
    diff.x = PNANOVDB_DEREF(in_ijk).x - PNANOVDB_DEREF(acc).m_tijk.x;
    diff.y = PNANOVDB_DEREF(in_ijk).y - PNANOVDB_DEREF(acc).m_tijk.y;
    diff.z = PNANOVDB_DEREF(in_ijk).z - PNANOVDB_DEREF(acc).m_tijk.z;

    pnanovdb_uint32_t       hash = diff.x | diff.y | diff.z;
    // PNANOVDB_DEREF(acc).m_depth is the depth of the tile (also for coalesced tiles).
    pnanovdb_uint32_t tvw = 4 << (PNANOVTT_TRAITS_MAX_DEPTH - PNANOVDB_DEREF(acc).m_depth);
    PNANOVDB_DEREF(acc).m_vw = PNANOVTT_TRAITS_MAX_DEPTH - PNANOVDB_DEREF(acc).m_depth + 2;
    pnanovdb_int32_t mask = (0xFFFFFFFF << (PNANOVTT_TRAITS_MAX_DEPTH - PNANOVDB_DEREF(acc).m_depth + 2));
    while ((hash >= tvw || PNANOVDB_DEREF(acc).m_depth > in_max_depth) &&
        PNANOVDB_DEREF(acc).m_depth > PNANOVDB_DEREF(acc).m_root_depth) {
#ifdef NANOVTT_SUPPORT_COALESCED_TILES
        pnanovdb_int32_t t_shift = pnanovtt_container_get_parent_and_shift(buf, PNANOVDB_DEREF(acc).m_container, PNANOVDB_REF(PNANOVDB_DEREF(acc).m_tile));
        PNANOVDB_DEREF(acc).m_depth -= t_shift;
        PNANOVDB_DEREF(acc).m_vw += t_shift;
        mask = pnanovdb_uint32_as_int32(pnanovdb_int32_as_uint32(mask) << pnanovdb_int32_as_uint32(t_shift));
        tvw = (tvw << t_shift);
#else
        pnanovtt_container_get_parent(buf, PNANOVDB_DEREF(acc).m_container, PNANOVDB_REF(PNANOVDB_DEREF(acc).m_tile));
        --PNANOVDB_DEREF(acc).m_depth;
        ++PNANOVDB_DEREF(acc).m_vw;
        mask *= 2;
        tvw = (tvw << 1);
#endif
        PNANOVDB_DEREF(acc).m_tijk.x = PNANOVDB_DEREF(acc).m_tijk.x & (mask);
        PNANOVDB_DEREF(acc).m_tijk.y = PNANOVDB_DEREF(acc).m_tijk.y & (mask);
        PNANOVDB_DEREF(acc).m_tijk.z = PNANOVDB_DEREF(acc).m_tijk.z & (mask);
        diff.x = PNANOVDB_DEREF(in_ijk).x - PNANOVDB_DEREF(acc).m_tijk.x;
        diff.y = PNANOVDB_DEREF(in_ijk).y - PNANOVDB_DEREF(acc).m_tijk.y;
        diff.z = PNANOVDB_DEREF(in_ijk).z - PNANOVDB_DEREF(acc).m_tijk.z;
        hash = (diff.x) |
            (diff.y) |
            (diff.z);
    }

    // At this point we've found a containing tile or we're at the root depth.
    if ((hash >= (tvw) || PNANOVDB_DEREF(acc).m_tile == PNANOVTT_INVALID_TILEINDEX) &&
        PNANOVDB_DEREF(acc).m_depth == PNANOVDB_DEREF(acc).m_root_depth) {
        // The in_ijk is not inside the current root
        PNANOVDB_DEREF(acc).m_tile = pnanovtt_container_get_root(buf, PNANOVDB_DEREF(acc).m_container, in_ijk, PNANOVDB_REF(PNANOVDB_DEREF(acc).m_tijk));
        if (PNANOVDB_DEREF(acc).m_tile == PNANOVTT_INVALID_TILEINDEX) {
            PNANOVDB_DEREF(acc).m_vw = PNANOVTT_TRAITS_MAX_DEPTH - PNANOVDB_DEREF(acc).m_root_depth;
            PNANOVDB_DEREF(acc).m_ijk.x = 0;
            PNANOVDB_DEREF(acc).m_ijk.y = 0;
            PNANOVDB_DEREF(acc).m_ijk.z = 0;
            PNANOVDB_DEREF(acc).m_tijk.x = PNANOVTT_TRAITS_TILE_DOMAIN_WIDTH;
            PNANOVDB_DEREF(acc).m_tijk.y = PNANOVTT_TRAITS_TILE_DOMAIN_WIDTH;
            PNANOVDB_DEREF(acc).m_tijk.z = PNANOVTT_TRAITS_TILE_DOMAIN_WIDTH;
            PNANOVDB_DEREF(acc).m_found = PNANOVDB_FALSE;
            return pnanovdb_address_offset(PNANOVDB_DEREF(acc).m_tree.address, PNANOVTT_TREE_OFF_DEFAULT_VALUE);
        }
        else {
            diff.x = PNANOVDB_DEREF(in_ijk).x - PNANOVDB_DEREF(acc).m_tijk.x;
            diff.y = PNANOVDB_DEREF(in_ijk).y - PNANOVDB_DEREF(acc).m_tijk.y;
            diff.z = PNANOVDB_DEREF(in_ijk).z - PNANOVDB_DEREF(acc).m_tijk.z;
        }
    }

    // The performance of this if-statement is ok because there's typically no branch-divergence.
    if (PNANOVDB_DEREF(acc).m_depth < in_max_depth) {
        PNANOVDB_DEREF(acc).m_vw = PNANOVTT_TRAITS_MAX_DEPTH - PNANOVDB_DEREF(acc).m_depth + 1;
        PNANOVTT_TILEINDEX child_tile = PNANOVTT_INVALID_TILEINDEX;
#ifdef NANOVTT_SUPPORT_COALESCED_TILES
        pnanovdb_int32_t t_shift = pnanovtt_container_get_child_index_and_shift(
#else
        pnanovtt_container_get_child_index(
#endif
            buf,
            PNANOVDB_DEREF(acc).m_container,
            PNANOVDB_DEREF(acc).m_tile,
            PNANOVDB_REF(diff),
            PNANOVDB_DEREF(acc).m_vw,
            PNANOVDB_REF(PNANOVDB_DEREF(acc).m_ijk),
            PNANOVDB_REF(child_tile));
        pnanovdb_bool_t done = (child_tile == PNANOVTT_INVALID_TILEINDEX);
        while (!done) {
            PNANOVDB_DEREF(acc).m_tile = child_tile;
#ifdef NANOVTT_SUPPORT_COALESCED_TILES
            PNANOVDB_DEREF(acc).m_depth += t_shift;
#else
            ++PNANOVDB_DEREF(acc).m_depth;
#endif
            PNANOVDB_DEREF(acc).m_tijk.x += (PNANOVDB_DEREF(acc).m_ijk.x << PNANOVDB_DEREF(acc).m_vw);
            PNANOVDB_DEREF(acc).m_tijk.y += (PNANOVDB_DEREF(acc).m_ijk.y << PNANOVDB_DEREF(acc).m_vw);
            PNANOVDB_DEREF(acc).m_tijk.z += (PNANOVDB_DEREF(acc).m_ijk.z << PNANOVDB_DEREF(acc).m_vw);
            diff.x = PNANOVDB_DEREF(in_ijk).x - PNANOVDB_DEREF(acc).m_tijk.x;
            diff.y = PNANOVDB_DEREF(in_ijk).y - PNANOVDB_DEREF(acc).m_tijk.y;
            diff.z = PNANOVDB_DEREF(in_ijk).z - PNANOVDB_DEREF(acc).m_tijk.z;
#ifdef NANOVTT_SUPPORT_COALESCED_TILES
            PNANOVDB_DEREF(acc).m_vw -= t_shift;
            t_shift = pnanovtt_container_get_child_index_and_shift(
#else
            --PNANOVDB_DEREF(acc).m_vw;
            pnanovtt_container_get_child_index(
#endif
                buf,
                PNANOVDB_DEREF(acc).m_container,
                PNANOVDB_DEREF(acc).m_tile,
                PNANOVDB_REF(diff),
                PNANOVDB_DEREF(acc).m_vw,
                PNANOVDB_REF(PNANOVDB_DEREF(acc).m_ijk),
                PNANOVDB_REF(child_tile));
            done = (child_tile == PNANOVTT_INVALID_TILEINDEX ||
                PNANOVDB_DEREF(acc).m_depth >= in_max_depth);
        }
    }
    PNANOVDB_DEREF(acc).m_vw = PNANOVTT_TRAITS_MAX_DEPTH - PNANOVDB_DEREF(acc).m_depth;

    return pnanovtt_tree_get_value_address(grid_type,
        grid_buf, PNANOVDB_DEREF(acc).m_tree, PNANOVDB_DEREF(acc).m_tile,
        PNANOVDB_REF(diff), PNANOVDB_REF(PNANOVDB_DEREF(acc).m_vw),
        PNANOVDB_REF(PNANOVDB_DEREF(acc).m_tile_data_stride),
        PNANOVDB_REF(PNANOVDB_DEREF(acc).m_ijk));
}

PNANOVDB_FORCE_INLINE void pnanovtt_readaccessor_get_value_texture_coord_at_depth(pnanovdb_grid_type_t grid_type, pnanovdb_buf_t buf, PNANOVDB_INOUT(pnanovtt_readaccessor_t) acc, PNANOVDB_IN(pnanovdb_coord_t) in_ijk, pnanovdb_int32_t in_max_depth, PNANOVDB_INOUT(pnanovdb_coord_t) io_coord)
{
    // This version traverses bottom up from the last tile

    pnanovdb_coord_t diff;
    diff.x = PNANOVDB_DEREF(in_ijk).x - PNANOVDB_DEREF(acc).m_tijk.x;
    diff.y = PNANOVDB_DEREF(in_ijk).y - PNANOVDB_DEREF(acc).m_tijk.y;
    diff.z = PNANOVDB_DEREF(in_ijk).z - PNANOVDB_DEREF(acc).m_tijk.z;

    // PNANOVDB_DEREF(acc).m_depth is the depth of the tile (also for coalesced tiles).
    pnanovdb_uint32_t tvw = 4 << (PNANOVTT_TRAITS_MAX_DEPTH - PNANOVDB_DEREF(acc).m_depth);
    PNANOVDB_DEREF(acc).m_vw = PNANOVTT_TRAITS_MAX_DEPTH - PNANOVDB_DEREF(acc).m_depth + 2;
    pnanovdb_int32_t mask = (0xFFFFFFFF << (PNANOVTT_TRAITS_MAX_DEPTH - PNANOVDB_DEREF(acc).m_depth + 2));
    while ((pnanovdb_int32_as_uint32(diff.x | diff.y | diff.z) >= tvw || PNANOVDB_DEREF(acc).m_depth > in_max_depth) &&
        PNANOVDB_DEREF(acc).m_depth > PNANOVDB_DEREF(acc).m_root_depth) {
#ifdef NANOVTT_SUPPORT_COALESCED_TILES
        pnanovdb_int32_t t_shift = pnanovtt_container_get_parent_and_shift(buf, PNANOVDB_DEREF(acc).m_container, PNANOVDB_REF(PNANOVDB_DEREF(acc).m_tile));
        PNANOVDB_DEREF(acc).m_depth -= t_shift;
        PNANOVDB_DEREF(acc).m_vw += t_shift;
        mask = pnanovdb_uint32_as_int32(pnanovdb_int32_as_uint32(mask) << pnanovdb_int32_as_uint32(t_shift));
        tvw = (tvw << t_shift);
#else
        pnanovtt_container_get_parent(buf, PNANOVDB_DEREF(acc).m_container, PNANOVDB_REF(PNANOVDB_DEREF(acc).m_tile));
        --PNANOVDB_DEREF(acc).m_depth;
        ++PNANOVDB_DEREF(acc).m_vw;
        mask *= 2;
        tvw = (tvw << 1);
#endif
        PNANOVDB_DEREF(acc).m_tijk.x = PNANOVDB_DEREF(acc).m_tijk.x & (mask);
        PNANOVDB_DEREF(acc).m_tijk.y = PNANOVDB_DEREF(acc).m_tijk.y & (mask);
        PNANOVDB_DEREF(acc).m_tijk.z = PNANOVDB_DEREF(acc).m_tijk.z & (mask);
        diff.x = PNANOVDB_DEREF(in_ijk).x - PNANOVDB_DEREF(acc).m_tijk.x;
        diff.y = PNANOVDB_DEREF(in_ijk).y - PNANOVDB_DEREF(acc).m_tijk.y;
        diff.z = PNANOVDB_DEREF(in_ijk).z - PNANOVDB_DEREF(acc).m_tijk.z;
    }

    // At this point we've found a containing tile or we're at the root depth.
    if ((pnanovdb_int32_as_uint32(diff.x | diff.y | diff.z) >= (tvw) || PNANOVDB_DEREF(acc).m_tile == PNANOVTT_INVALID_TILEINDEX) &&
        PNANOVDB_DEREF(acc).m_depth == PNANOVDB_DEREF(acc).m_root_depth) {
        // The in_ijk is not inside the current root
        PNANOVDB_DEREF(acc).m_tile = pnanovtt_container_get_root(buf, PNANOVDB_DEREF(acc).m_container, in_ijk, PNANOVDB_REF(PNANOVDB_DEREF(acc).m_tijk));
        if (PNANOVDB_DEREF(acc).m_tile == PNANOVTT_INVALID_TILEINDEX) {
            PNANOVDB_DEREF(acc).m_vw = PNANOVTT_TRAITS_MAX_DEPTH - PNANOVDB_DEREF(acc).m_root_depth;
            PNANOVDB_DEREF(acc).m_found = PNANOVDB_FALSE;
            // Here we use the convention that the first tile data in the texture
            // is a tile data of default values.
            PNANOVDB_DEREF(io_coord).x = 0;
            PNANOVDB_DEREF(io_coord).y = 0;
            PNANOVDB_DEREF(io_coord).z = 0;
            return;
        }
        else {
            diff.x = PNANOVDB_DEREF(in_ijk).x - PNANOVDB_DEREF(acc).m_tijk.x;
            diff.y = PNANOVDB_DEREF(in_ijk).y - PNANOVDB_DEREF(acc).m_tijk.y;
            diff.z = PNANOVDB_DEREF(in_ijk).z - PNANOVDB_DEREF(acc).m_tijk.z;
        }
    }
    PNANOVDB_DEREF(acc).m_found = PNANOVDB_TRUE;

    // The performance of this if-statement is ok because there's typically no branch-divergence.
    if (PNANOVDB_DEREF(acc).m_depth < in_max_depth) {
        PNANOVDB_DEREF(acc).m_vw = PNANOVTT_TRAITS_MAX_DEPTH - PNANOVDB_DEREF(acc).m_depth + 1;
        PNANOVTT_TILEINDEX child_tile = PNANOVTT_INVALID_TILEINDEX;
#ifdef NANOVTT_SUPPORT_COALESCED_TILES
        pnanovdb_int32_t t_shift = pnanovtt_container_get_child_index_and_shift(
#else
        pnanovtt_container_get_child_index(
#endif
            buf,
            PNANOVDB_DEREF(acc).m_container,
            PNANOVDB_DEREF(acc).m_tile,
            PNANOVDB_REF(diff),
            PNANOVDB_DEREF(acc).m_vw,
            PNANOVDB_REF(PNANOVDB_DEREF(acc).m_ijk),
            PNANOVDB_REF(child_tile));
        pnanovdb_bool_t done = (child_tile == PNANOVTT_INVALID_TILEINDEX);
        while (!done) {
            PNANOVDB_DEREF(acc).m_tile = child_tile;
#ifdef NANOVTT_SUPPORT_COALESCED_TILES
            PNANOVDB_DEREF(acc).m_depth += t_shift;
#else
            ++PNANOVDB_DEREF(acc).m_depth;
#endif
            PNANOVDB_DEREF(acc).m_tijk.x += (PNANOVDB_DEREF(acc).m_ijk.x << PNANOVDB_DEREF(acc).m_vw);
            PNANOVDB_DEREF(acc).m_tijk.y += (PNANOVDB_DEREF(acc).m_ijk.y << PNANOVDB_DEREF(acc).m_vw);
            PNANOVDB_DEREF(acc).m_tijk.z += (PNANOVDB_DEREF(acc).m_ijk.z << PNANOVDB_DEREF(acc).m_vw);
            diff.x = PNANOVDB_DEREF(in_ijk).x - PNANOVDB_DEREF(acc).m_tijk.x;
            diff.y = PNANOVDB_DEREF(in_ijk).y - PNANOVDB_DEREF(acc).m_tijk.y;
            diff.z = PNANOVDB_DEREF(in_ijk).z - PNANOVDB_DEREF(acc).m_tijk.z;
#ifdef NANOVTT_SUPPORT_COALESCED_TILES
            PNANOVDB_DEREF(acc).m_vw -= t_shift;
            t_shift = pnanovtt_container_get_child_index_and_shift(
#else
            --PNANOVDB_DEREF(acc).m_vw;
            pnanovtt_container_get_child_index(
#endif
                buf,
                PNANOVDB_DEREF(acc).m_container,
                PNANOVDB_DEREF(acc).m_tile,
                PNANOVDB_REF(diff),
                PNANOVDB_DEREF(acc).m_vw,
                PNANOVDB_REF(PNANOVDB_DEREF(acc).m_ijk),
                PNANOVDB_REF(child_tile));
            done = (child_tile == PNANOVTT_INVALID_TILEINDEX ||
                PNANOVDB_DEREF(acc).m_depth >= in_max_depth);
        }
    }
    PNANOVDB_DEREF(acc).m_vw = PNANOVTT_TRAITS_MAX_DEPTH - PNANOVDB_DEREF(acc).m_depth;

    pnanovtt_tree_get_value_texture_coord(grid_type,
        buf, PNANOVDB_DEREF(acc).m_tree, PNANOVDB_DEREF(acc).m_tile,
        PNANOVDB_REF(diff), PNANOVDB_REF(PNANOVDB_DEREF(acc).m_vw),
        PNANOVDB_REF(PNANOVDB_DEREF(acc).m_ijk), io_coord);
}

PNANOVDB_FORCE_INLINE void pnanovtt_readaccessor_get_value_texture_coord_at_depth_from_root(pnanovdb_grid_type_t grid_type, pnanovdb_buf_t buf, PNANOVDB_INOUT(pnanovtt_readaccessor_t) acc, PNANOVDB_IN(pnanovdb_coord_t) in_ijk, pnanovdb_int32_t in_max_depth, PNANOVDB_INOUT(pnanovdb_coord_t) io_coord)
{
    PNANOVDB_DEREF(acc).m_tile = pnanovtt_container_get_root(
        buf, PNANOVDB_DEREF(acc).m_container, in_ijk,
        PNANOVDB_REF(PNANOVDB_DEREF(acc).m_tijk));
    PNANOVDB_DEREF(acc).m_depth = PNANOVDB_DEREF(acc).m_root_depth;
    if (PNANOVDB_DEREF(acc).m_tile == PNANOVTT_INVALID_TILEINDEX) {
        PNANOVDB_DEREF(acc).m_vw =
            PNANOVTT_TRAITS_MAX_DEPTH - PNANOVDB_DEREF(acc).m_root_depth;
        PNANOVDB_DEREF(acc).m_found = PNANOVDB_FALSE;
        // Here we use the convention that the first tile data in the texture
        // is a tile data of default values.
        PNANOVDB_DEREF(io_coord).x = 0;
        PNANOVDB_DEREF(io_coord).y = 0;
        PNANOVDB_DEREF(io_coord).z = 0;
        return;
    }

    PNANOVDB_DEREF(acc).m_found = PNANOVDB_TRUE;
    pnanovdb_coord_t diff;
    diff.x = PNANOVDB_DEREF(in_ijk).x - PNANOVDB_DEREF(acc).m_tijk.x;
    diff.y = PNANOVDB_DEREF(in_ijk).y - PNANOVDB_DEREF(acc).m_tijk.y;
    diff.z = PNANOVDB_DEREF(in_ijk).z - PNANOVDB_DEREF(acc).m_tijk.z;
    PNANOVDB_DEREF(acc).m_vw =
        PNANOVTT_TRAITS_MAX_DEPTH - PNANOVDB_DEREF(acc).m_depth + 1;
    PNANOVTT_TILEINDEX child_tile = PNANOVTT_INVALID_TILEINDEX;
#ifdef NANOVTT_SUPPORT_COALESCED_TILES
    pnanovdb_int32_t t_shift = pnanovtt_container_get_child_index_and_shift(
#else
    pnanovtt_container_get_child_index(
#endif
        buf, PNANOVDB_DEREF(acc).m_container, PNANOVDB_DEREF(acc).m_tile,
        PNANOVDB_REF(diff), PNANOVDB_DEREF(acc).m_vw,
        PNANOVDB_REF(PNANOVDB_DEREF(acc).m_ijk), PNANOVDB_REF(child_tile));
    while (child_tile != PNANOVTT_INVALID_TILEINDEX &&
        PNANOVDB_DEREF(acc).m_depth < in_max_depth)
    {
        PNANOVDB_DEREF(acc).m_tile = child_tile;
#ifdef NANOVTT_SUPPORT_COALESCED_TILES
        PNANOVDB_DEREF(acc).m_depth += t_shift;
#else
        ++PNANOVDB_DEREF(acc).m_depth;
#endif
        PNANOVDB_DEREF(acc).m_tijk.x +=
            (PNANOVDB_DEREF(acc).m_ijk.x << PNANOVDB_DEREF(acc).m_vw);
        PNANOVDB_DEREF(acc).m_tijk.y +=
            (PNANOVDB_DEREF(acc).m_ijk.y << PNANOVDB_DEREF(acc).m_vw);
        PNANOVDB_DEREF(acc).m_tijk.z +=
            (PNANOVDB_DEREF(acc).m_ijk.z << PNANOVDB_DEREF(acc).m_vw);
        diff.x = PNANOVDB_DEREF(in_ijk).x - PNANOVDB_DEREF(acc).m_tijk.x;
        diff.y = PNANOVDB_DEREF(in_ijk).y - PNANOVDB_DEREF(acc).m_tijk.y;
        diff.z = PNANOVDB_DEREF(in_ijk).z - PNANOVDB_DEREF(acc).m_tijk.z;
#ifdef NANOVTT_SUPPORT_COALESCED_TILES
        PNANOVDB_DEREF(acc).m_vw -= t_shift;
        t_shift = pnanovtt_container_get_child_index_and_shift(
#else
        --PNANOVDB_DEREF(acc).m_vw;
        pnanovtt_container_get_child_index(
#endif
            buf, PNANOVDB_DEREF(acc).m_container, PNANOVDB_DEREF(acc).m_tile,
            PNANOVDB_REF(diff), PNANOVDB_DEREF(acc).m_vw,
            PNANOVDB_REF(PNANOVDB_DEREF(acc).m_ijk), PNANOVDB_REF(child_tile));
    }
    PNANOVDB_DEREF(acc).m_vw = PNANOVTT_TRAITS_MAX_DEPTH - PNANOVDB_DEREF(acc).m_depth;

    pnanovtt_tree_get_value_texture_coord(grid_type,
        buf, PNANOVDB_DEREF(acc).m_tree, PNANOVDB_DEREF(acc).m_tile,
        PNANOVDB_REF(diff), PNANOVDB_REF(PNANOVDB_DEREF(acc).m_vw),
        PNANOVDB_REF(PNANOVDB_DEREF(acc).m_ijk), io_coord);
}

PNANOVDB_FORCE_INLINE void pnanovtt_readaccessor_get_value_texture_coord_at_coarser_depth(pnanovdb_grid_type_t grid_type, pnanovdb_buf_t buf, PNANOVDB_INOUT(pnanovtt_readaccessor_t) acc, PNANOVDB_IN(pnanovdb_coord_t) in_ijk, pnanovdb_int32_t in_max_depth, PNANOVDB_INOUT(pnanovdb_coord_t) io_coord)
{
    // This version traverses bottom up from the last tile

    pnanovdb_coord_t diff;
    diff.x = PNANOVDB_DEREF(in_ijk).x - PNANOVDB_DEREF(acc).m_tijk.x;
    diff.y = PNANOVDB_DEREF(in_ijk).y - PNANOVDB_DEREF(acc).m_tijk.y;
    diff.z = PNANOVDB_DEREF(in_ijk).z - PNANOVDB_DEREF(acc).m_tijk.z;

    // PNANOVDB_DEREF(acc).m_depth is the depth of the tile (also for coalesced tiles).
    pnanovdb_uint32_t tvw = 4 << (PNANOVTT_TRAITS_MAX_DEPTH - PNANOVDB_DEREF(acc).m_depth);
    PNANOVDB_DEREF(acc).m_vw = PNANOVTT_TRAITS_MAX_DEPTH - PNANOVDB_DEREF(acc).m_depth + 2;
    pnanovdb_int32_t mask = (0xFFFFFFFF << (PNANOVTT_TRAITS_MAX_DEPTH - PNANOVDB_DEREF(acc).m_depth + 2));
    while ((pnanovdb_int32_as_uint32(diff.x | diff.y | diff.z) >= tvw || PNANOVDB_DEREF(acc).m_depth > in_max_depth) &&
        PNANOVDB_DEREF(acc).m_depth > PNANOVDB_DEREF(acc).m_root_depth) {
#ifdef NANOVTT_SUPPORT_COALESCED_TILES
        pnanovdb_int32_t t_shift = pnanovtt_container_get_parent_and_shift(buf, PNANOVDB_DEREF(acc).m_container, PNANOVDB_REF(PNANOVDB_DEREF(acc).m_tile));
        PNANOVDB_DEREF(acc).m_depth -= t_shift;
        PNANOVDB_DEREF(acc).m_vw += t_shift;
        mask = pnanovdb_uint32_as_int32(pnanovdb_int32_as_uint32(mask) << pnanovdb_int32_as_uint32(t_shift));
        tvw = (tvw << t_shift);
#else
        pnanovtt_container_get_parent(buf, PNANOVDB_DEREF(acc).m_container, PNANOVDB_REF(PNANOVDB_DEREF(acc).m_tile));
        --PNANOVDB_DEREF(acc).m_depth;
        ++PNANOVDB_DEREF(acc).m_vw;
        mask *= 2;
        tvw = (tvw << 1);
#endif
        PNANOVDB_DEREF(acc).m_tijk.x = PNANOVDB_DEREF(acc).m_tijk.x & (mask);
        PNANOVDB_DEREF(acc).m_tijk.y = PNANOVDB_DEREF(acc).m_tijk.y & (mask);
        PNANOVDB_DEREF(acc).m_tijk.z = PNANOVDB_DEREF(acc).m_tijk.z & (mask);
        diff.x = PNANOVDB_DEREF(in_ijk).x - PNANOVDB_DEREF(acc).m_tijk.x;
        diff.y = PNANOVDB_DEREF(in_ijk).y - PNANOVDB_DEREF(acc).m_tijk.y;
        diff.z = PNANOVDB_DEREF(in_ijk).z - PNANOVDB_DEREF(acc).m_tijk.z;
    }

    // At this point we've found a containing tile (guaranteed per the assumptions).
    PNANOVDB_DEREF(acc).m_found = PNANOVDB_TRUE;
    PNANOVDB_DEREF(acc).m_vw = PNANOVTT_TRAITS_MAX_DEPTH - PNANOVDB_DEREF(acc).m_depth;

    pnanovtt_tree_get_value_texture_coord(grid_type,
        buf, PNANOVDB_DEREF(acc).m_tree, PNANOVDB_DEREF(acc).m_tile,
        PNANOVDB_REF(diff), PNANOVDB_REF(PNANOVDB_DEREF(acc).m_vw),
        PNANOVDB_REF(PNANOVDB_DEREF(acc).m_ijk), io_coord);
}

PNANOVDB_FORCE_INLINE pnanovdb_address_t pnanovtt_readaccessor_get_value_address_at_depth(pnanovdb_grid_type_t grid_type, pnanovdb_buf_t buf, PNANOVDB_INOUT(pnanovtt_readaccessor_t) acc, PNANOVDB_IN(pnanovdb_coord_t) in_ijk, pnanovdb_int32_t in_max_depth)
{
    return pnanovtt_readaccessor_get_value_address_at_depth_v2(grid_type, buf, buf, acc, in_ijk, in_max_depth);
}

PNANOVDB_FORCE_INLINE pnanovdb_address_t pnanovtt_readaccessor_get_value_address(pnanovdb_grid_type_t grid_type, pnanovdb_buf_t buf, PNANOVDB_INOUT(pnanovtt_readaccessor_t) acc, PNANOVDB_IN(pnanovdb_coord_t) in_ijk)
{
    return pnanovtt_readaccessor_get_value_address_at_depth(grid_type, buf, acc, in_ijk, PNANOVTT_TRAITS_MAX_DEPTH);
}

PNANOVDB_FORCE_INLINE pnanovdb_address_t pnanovtt_readaccessor_get_offset_extended_value_address(pnanovdb_grid_type_t grid_type, pnanovdb_buf_t buf, PNANOVDB_INOUT(pnanovtt_readaccessor_t) acc, PNANOVDB_IN(pnanovdb_coord_t) in_off)
{
    pnanovdb_coord_t ijk;
    ijk.x = PNANOVDB_DEREF(acc).m_ijk.x + PNANOVDB_DEREF(in_off).x;
    ijk.y = PNANOVDB_DEREF(acc).m_ijk.y + PNANOVDB_DEREF(in_off).y;
    ijk.z = PNANOVDB_DEREF(acc).m_ijk.z + PNANOVDB_DEREF(in_off).z;

    return pnanovtt_tree_get_value_address2(grid_type,
        buf, PNANOVDB_DEREF(acc).m_tree, PNANOVDB_DEREF(acc).m_tile,
        PNANOVDB_REF(ijk));
}

// ------------------------------------------------ ReadAccessor GetDim -----------------------------------------------------------

PNANOVDB_FORCE_INLINE pnanovdb_uint32_t pnanovtt_readaccessor_get_dim_at_depth(pnanovdb_grid_type_t grid_type, pnanovdb_buf_t buf, PNANOVDB_INOUT(pnanovtt_readaccessor_t) acc, PNANOVDB_IN(pnanovdb_coord_t) ijk, pnanovdb_int32_t depth)
{
    pnanovtt_readaccessor_get_value_address_at_depth(grid_type, buf, acc, ijk, depth);
    return 1 << PNANOVDB_DEREF(acc).m_vw;
}

PNANOVDB_FORCE_INLINE pnanovdb_uint32_t pnanovtt_readaccessor_get_dim(pnanovdb_grid_type_t grid_type, pnanovdb_buf_t buf, PNANOVDB_INOUT(pnanovtt_readaccessor_t) acc, PNANOVDB_IN(pnanovdb_coord_t) ijk)
{
    return pnanovtt_readaccessor_get_dim_at_depth(grid_type, buf, acc, ijk, PNANOVTT_TRAITS_MAX_DEPTH);
}

PNANOVDB_FORCE_INLINE pnanovdb_uint32_t pnanovtt_readaccessor_get_dim2(PNANOVDB_INOUT(pnanovtt_readaccessor_t) acc)
{
    return 1 << PNANOVDB_DEREF(acc).m_vw;
}

// ------------------------------------------------ ReadAccessor IsActive -----------------------------------------------------------

PNANOVDB_FORCE_INLINE pnanovdb_bool_t pnanovtt_readaccessor_is_texture_active_at_depth(pnanovdb_grid_type_t grid_type, pnanovdb_buf_t buf, pnanovtt_texture_t tex, PNANOVDB_INOUT(pnanovtt_readaccessor_t) acc, PNANOVDB_IN(pnanovdb_coord_t) ijk, pnanovdb_int32_t in_max_depth)
{
    if (PNANOVDB_DEREF(acc).m_is_level_set) {
        return PNANOVDB_TRUE;
    }
    else {
        pnanovdb_coord_t coord;
        pnanovtt_readaccessor_get_value_texture_coord_at_depth(grid_type, buf, acc, ijk, in_max_depth, PNANOVDB_REF(coord));
        return pnanovtt_read_texture_float(tex, PNANOVDB_REF(coord)) > 0;
    }
}

PNANOVDB_FORCE_INLINE pnanovdb_bool_t pnanovtt_readaccessor_is_active_at_depth(pnanovdb_grid_type_t grid_type, pnanovdb_buf_t buf, PNANOVDB_INOUT(pnanovtt_readaccessor_t) acc, PNANOVDB_IN(pnanovdb_coord_t) ijk, pnanovdb_int32_t in_max_depth)
{
    // When using the HDDA + nanovtt::ZeroCrossing method it's faster to always return true if m_isLevelSet is true.
    if (PNANOVDB_DEREF(acc).m_is_level_set) {
        return PNANOVDB_TRUE;
    }
    else {
        pnanovdb_address_t address = pnanovtt_readaccessor_get_value_address_at_depth(grid_type, buf, acc, ijk, in_max_depth);
        return pnanovdb_read_float(buf, address) > 0;
    }
}

PNANOVDB_FORCE_INLINE pnanovdb_bool_t pnanovtt_readaccessor_is_active(pnanovdb_grid_type_t grid_type, pnanovdb_buf_t buf, PNANOVDB_INOUT(pnanovtt_readaccessor_t) acc, PNANOVDB_IN(pnanovdb_coord_t) ijk)
{
    return pnanovtt_readaccessor_is_active_at_depth(grid_type, buf, acc, ijk, PNANOVTT_TRAITS_MAX_DEPTH);
}

// ------------------------------------------------ Float Sampling -----------------------------------------------------------

struct pnanovtt_stencil_t
{
    float m_v[8];
};
PNANOVDB_STRUCT_TYPEDEF(pnanovtt_stencil_t)

PNANOVDB_FORCE_INLINE pnanovdb_address_t pnanovtt_remap_stencil_base_address(pnanovdb_address_t base_address,
    PNANOVTT_TILEINDEX num_tiles,
    pnanovtt_tree_handle_t from_tree,
    pnanovtt_tree_handle_t to_tree,
    pnanovtt_grid_type_t from_grid_type,
    pnanovtt_grid_type_t to_grid_type,
    pnanovdb_uint32_t from_type_size,
    pnanovdb_uint32_t to_type_size)
{
#ifdef PNANOVDB_C
    // Assumptions:
    assert(from_type_size <= to_type_size);
#endif

    return pnanovdb_address_offset(
    pnanovdb_address_product(
        pnanovdb_address_offset_neg(
            base_address,
            (num_tiles + 1) * 8 +
            PNANOVTT_GRID_TYPE_GET(
                from_grid_type,
                tree_off_data) +
            from_tree.address
            .byte_offset),
        to_type_size / from_type_size),
    (num_tiles + 1) * 8 +
    PNANOVTT_GRID_TYPE_GET(
        to_grid_type,
        tree_off_data) +
    to_tree.address
    .byte_offset);
}

// Assumes that the grid point exists (ie that acc.m_found is true).
PNANOVDB_FORCE_INLINE pnanovdb_address_t pnanovtt_stencil_base_address_at_depth(pnanovdb_buf_t buf,
    pnanovdb_buf_t grid_buf,
    pnanovdb_grid_type_t grid_type,
    PNANOVDB_INOUT(pnanovtt_readaccessor_t) acc,
    PNANOVDB_IN(pnanovdb_vec3_t) in_xyz,
    PNANOVDB_INOUT(pnanovdb_vec3_t) io_frac,
    pnanovdb_int32_t in_depth)
{
    pnanovdb_coord_t ijk;
    ijk.x = pnanovdb_float_to_int32(pnanovdb_floor(PNANOVDB_DEREF(in_xyz).x));
    ijk.y = pnanovdb_float_to_int32(pnanovdb_floor(PNANOVDB_DEREF(in_xyz).y));
    ijk.z = pnanovdb_float_to_int32(pnanovdb_floor(PNANOVDB_DEREF(in_xyz).z));
    pnanovdb_address_t address = pnanovtt_readaccessor_get_value_address_at_depth_v2(grid_type, buf, grid_buf, acc, PNANOVDB_REF(ijk), in_depth);
    int              shift = PNANOVDB_DEREF(acc).m_vw;
    float factor = 2.f / (1 << (shift + 1));
    pnanovdb_vec3_t ijkf;
    ijkf.x = (PNANOVDB_DEREF(acc).m_tijk.x + (PNANOVDB_DEREF(acc).m_ijk.x << shift));
    ijkf.y = (PNANOVDB_DEREF(acc).m_tijk.y + (PNANOVDB_DEREF(acc).m_ijk.y << shift));
    ijkf.z = (PNANOVDB_DEREF(acc).m_tijk.z + (PNANOVDB_DEREF(acc).m_ijk.z << shift));
    PNANOVDB_DEREF(io_frac).x = (PNANOVDB_DEREF(in_xyz).x - ijkf.x) * factor;
    PNANOVDB_DEREF(io_frac).y = (PNANOVDB_DEREF(in_xyz).y - ijkf.y) * factor;
    PNANOVDB_DEREF(io_frac).z = (PNANOVDB_DEREF(in_xyz).z - ijkf.z) * factor;
    return address;
}

// Assumes that the grid point exists (ie that acc.m_found is true).
PNANOVDB_FORCE_INLINE void pnanovtt_stencil_base_texture_coord_at_depth(pnanovdb_buf_t buf,
    pnanovdb_grid_type_t grid_type,
    PNANOVDB_INOUT(pnanovtt_readaccessor_t) acc,
    PNANOVDB_IN(pnanovdb_vec3_t) in_xyz,
    PNANOVDB_INOUT(pnanovdb_vec3_t) io_frac,
    PNANOVDB_INOUT(pnanovdb_coord_t) io_coord,
    pnanovdb_int32_t in_depth)
{
    pnanovdb_coord_t ijk;
    ijk.x = pnanovdb_float_to_int32(pnanovdb_floor(PNANOVDB_DEREF(in_xyz).x));
    ijk.y = pnanovdb_float_to_int32(pnanovdb_floor(PNANOVDB_DEREF(in_xyz).y));
    ijk.z = pnanovdb_float_to_int32(pnanovdb_floor(PNANOVDB_DEREF(in_xyz).z));
    pnanovtt_readaccessor_get_value_texture_coord_at_depth(grid_type, buf, acc, PNANOVDB_REF(ijk), in_depth, io_coord);
    float factor = 2.f / (1 << (PNANOVDB_DEREF(acc).m_vw + 1));
    PNANOVDB_DEREF(io_frac).x = (PNANOVDB_DEREF(in_xyz).x - (PNANOVDB_DEREF(acc).m_tijk.x + (PNANOVDB_DEREF(acc).m_ijk.x << PNANOVDB_DEREF(acc).m_vw))) * factor;
    PNANOVDB_DEREF(io_frac).y = (PNANOVDB_DEREF(in_xyz).y - (PNANOVDB_DEREF(acc).m_tijk.y + (PNANOVDB_DEREF(acc).m_ijk.y << PNANOVDB_DEREF(acc).m_vw))) * factor;
    PNANOVDB_DEREF(io_frac).z = (PNANOVDB_DEREF(in_xyz).z - (PNANOVDB_DEREF(acc).m_tijk.z + (PNANOVDB_DEREF(acc).m_ijk.z << PNANOVDB_DEREF(acc).m_vw))) * factor;
}

PNANOVDB_FORCE_INLINE void pnanovtt_stencil_base_texture_coord_at_depth_from_root(pnanovdb_buf_t buf,
    pnanovdb_grid_type_t grid_type,
    PNANOVDB_INOUT(pnanovtt_readaccessor_t) acc,
    PNANOVDB_IN(pnanovdb_vec3_t) in_xyz,
    PNANOVDB_INOUT(pnanovdb_vec3_t) io_frac,
    PNANOVDB_INOUT(pnanovdb_coord_t) io_coord,
    pnanovdb_int32_t in_depth)
{
    pnanovdb_coord_t ijk;
    ijk.x = pnanovdb_float_to_int32(pnanovdb_floor(PNANOVDB_DEREF(in_xyz).x));
    ijk.y = pnanovdb_float_to_int32(pnanovdb_floor(PNANOVDB_DEREF(in_xyz).y));
    ijk.z = pnanovdb_float_to_int32(pnanovdb_floor(PNANOVDB_DEREF(in_xyz).z));
    pnanovtt_readaccessor_get_value_texture_coord_at_depth_from_root(grid_type, buf, acc, PNANOVDB_REF(ijk), in_depth, io_coord);
    float factor = 2.f / (1 << (PNANOVDB_DEREF(acc).m_vw + 1));
    PNANOVDB_DEREF(io_frac).x = (PNANOVDB_DEREF(in_xyz).x - (PNANOVDB_DEREF(acc).m_tijk.x + (PNANOVDB_DEREF(acc).m_ijk.x << PNANOVDB_DEREF(acc).m_vw))) * factor;
    PNANOVDB_DEREF(io_frac).y = (PNANOVDB_DEREF(in_xyz).y - (PNANOVDB_DEREF(acc).m_tijk.y + (PNANOVDB_DEREF(acc).m_ijk.y << PNANOVDB_DEREF(acc).m_vw))) * factor;
    PNANOVDB_DEREF(io_frac).z = (PNANOVDB_DEREF(in_xyz).z - (PNANOVDB_DEREF(acc).m_tijk.z + (PNANOVDB_DEREF(acc).m_ijk.z << PNANOVDB_DEREF(acc).m_vw))) * factor;
}

// Assumes that in_depth is coarser than or equal to the current depth of the input accessor acc,
// and that the accessor points to a valid tile in the tree.
PNANOVDB_FORCE_INLINE void pnanovtt_stencil_base_texture_coord_at_coarser_depth(pnanovdb_buf_t buf,
    pnanovdb_grid_type_t grid_type,
    PNANOVDB_INOUT(pnanovtt_readaccessor_t) acc,
    PNANOVDB_IN(pnanovdb_vec3_t) in_xyz,
    PNANOVDB_INOUT(pnanovdb_vec3_t) io_frac,
    PNANOVDB_INOUT(pnanovdb_coord_t) io_coord,
    pnanovdb_int32_t in_depth)
{
#ifdef DEBUG
    assert(in_depth <= PNANOVDB_DEREF(acc).m_depth);
#endif
    pnanovdb_coord_t ijk;
    ijk.x = pnanovdb_float_to_int32(pnanovdb_floor(PNANOVDB_DEREF(in_xyz).x));
    ijk.y = pnanovdb_float_to_int32(pnanovdb_floor(PNANOVDB_DEREF(in_xyz).y));
    ijk.z = pnanovdb_float_to_int32(pnanovdb_floor(PNANOVDB_DEREF(in_xyz).z));
    pnanovtt_readaccessor_get_value_texture_coord_at_coarser_depth(grid_type, buf, acc, PNANOVDB_REF(ijk), in_depth, io_coord);
    float factor = 2.f / (1 << (PNANOVDB_DEREF(acc).m_vw + 1));
    PNANOVDB_DEREF(io_frac).x = (PNANOVDB_DEREF(in_xyz).x - (PNANOVDB_DEREF(acc).m_tijk.x + (PNANOVDB_DEREF(acc).m_ijk.x << PNANOVDB_DEREF(acc).m_vw))) * factor;
    PNANOVDB_DEREF(io_frac).y = (PNANOVDB_DEREF(in_xyz).y - (PNANOVDB_DEREF(acc).m_tijk.y + (PNANOVDB_DEREF(acc).m_ijk.y << PNANOVDB_DEREF(acc).m_vw))) * factor;
    PNANOVDB_DEREF(io_frac).z = (PNANOVDB_DEREF(in_xyz).z - (PNANOVDB_DEREF(acc).m_tijk.z + (PNANOVDB_DEREF(acc).m_ijk.z << PNANOVDB_DEREF(acc).m_vw))) * factor;
}

PNANOVDB_FORCE_INLINE void pnanovtt_stencil_from_base_address_float(pnanovdb_buf_t grid_buf,
    pnanovdb_address_t address,
    PNANOVDB_INOUT(pnanovtt_stencil_t) io_v)
{
    // In this case all values can be looked up from the extended tile itself
    PNANOVDB_DEREF(io_v).m_v[0] = pnanovdb_read_float(grid_buf, address); // i, j, k
    PNANOVDB_DEREF(io_v).m_v[1] = pnanovdb_read_float(grid_buf, pnanovdb_address_offset(address, PNANOVTT_TRAITS_FLOAT_EXTENDED_TILE_DATA_STRIDE_Z)); // i, j, k + 1
    PNANOVDB_DEREF(io_v).m_v[2] = pnanovdb_read_float(grid_buf, pnanovdb_address_offset(address, PNANOVTT_TRAITS_FLOAT_EXTENDED_TILE_DATA_STRIDE_Y)); // i, j+1, k
    PNANOVDB_DEREF(io_v).m_v[3] = pnanovdb_read_float(grid_buf, pnanovdb_address_offset(address, PNANOVTT_TRAITS_FLOAT_EXTENDED_TILE_DATA_STRIDE_Y + PNANOVTT_TRAITS_FLOAT_EXTENDED_TILE_DATA_STRIDE_Z)); // i, j+1, k + 1
    PNANOVDB_DEREF(io_v).m_v[4] = pnanovdb_read_float(grid_buf, pnanovdb_address_offset(address, PNANOVTT_TRAITS_FLOAT_EXTENDED_TILE_DATA_STRIDE_X)); // i+1, j, k
    PNANOVDB_DEREF(io_v).m_v[5] = pnanovdb_read_float(grid_buf, pnanovdb_address_offset(address, PNANOVTT_TRAITS_FLOAT_EXTENDED_TILE_DATA_STRIDE_X + PNANOVTT_TRAITS_FLOAT_EXTENDED_TILE_DATA_STRIDE_Z)); // i+1, j, k + 1
    PNANOVDB_DEREF(io_v).m_v[6] = pnanovdb_read_float(grid_buf, pnanovdb_address_offset(address, PNANOVTT_TRAITS_FLOAT_EXTENDED_TILE_DATA_STRIDE_X + PNANOVTT_TRAITS_FLOAT_EXTENDED_TILE_DATA_STRIDE_Y)); // i+1, j+1, k
    PNANOVDB_DEREF(io_v).m_v[7] = pnanovdb_read_float(grid_buf, pnanovdb_address_offset(address, PNANOVTT_TRAITS_FLOAT_EXTENDED_TILE_DATA_STRIDE_X + PNANOVTT_TRAITS_FLOAT_EXTENDED_TILE_DATA_STRIDE_Y + PNANOVTT_TRAITS_FLOAT_EXTENDED_TILE_DATA_STRIDE_Z)); // i+1, j+1, k + 1
}

PNANOVDB_FORCE_INLINE void pnanovtt_stencil_float(pnanovdb_buf_t buf,
    pnanovdb_buf_t grid_buf,
    PNANOVDB_INOUT(pnanovtt_readaccessor_t) acc,
    PNANOVDB_IN(pnanovdb_coord_t) in_ijk,
    PNANOVDB_INOUT(pnanovtt_stencil_t) io_v,
    PNANOVDB_INOUT(pnanovdb_vec3_t) io_vijk,
    PNANOVDB_INOUT(float) io_factor,
    pnanovdb_int32_t in_depth) {
    pnanovdb_address_t address = pnanovtt_readaccessor_get_value_address_at_depth_v2(PNANOVDB_GRID_TYPE_FLOAT, buf, grid_buf, acc, in_ijk, in_depth);
    PNANOVDB_DEREF(io_v).m_v[0] = pnanovdb_read_float(grid_buf, address); // i, j, k

    if (!PNANOVDB_DEREF(acc).m_found)
        return;

    // Account for varying voxelsize in the uvw offset
    pnanovdb_coord_t tijk = PNANOVDB_DEREF(acc).m_tijk;
    int              shift = PNANOVDB_DEREF(acc).m_vw;
    PNANOVDB_DEREF(io_factor) = 2.f / (1 << (shift + 1));
    PNANOVDB_DEREF(io_vijk).x = (tijk.x + (PNANOVDB_DEREF(acc).m_ijk.x << shift));
    PNANOVDB_DEREF(io_vijk).y = (tijk.y + (PNANOVDB_DEREF(acc).m_ijk.y << shift));
    PNANOVDB_DEREF(io_vijk).z = (tijk.z + (PNANOVDB_DEREF(acc).m_ijk.z << shift));

    // In this case all values can be looked up from the extended tile itself
    PNANOVDB_DEREF(io_v).m_v[1] = pnanovdb_read_float(grid_buf, pnanovdb_address_offset(address, PNANOVTT_TRAITS_FLOAT_EXTENDED_TILE_DATA_STRIDE_Z)); // i, j, k + 1
    PNANOVDB_DEREF(io_v).m_v[2] = pnanovdb_read_float(grid_buf, pnanovdb_address_offset(address, PNANOVTT_TRAITS_FLOAT_EXTENDED_TILE_DATA_STRIDE_Y)); // i, j+1, k
    PNANOVDB_DEREF(io_v).m_v[3] = pnanovdb_read_float(grid_buf, pnanovdb_address_offset(address, PNANOVTT_TRAITS_FLOAT_EXTENDED_TILE_DATA_STRIDE_Y + PNANOVTT_TRAITS_FLOAT_EXTENDED_TILE_DATA_STRIDE_Z)); // i, j+1, k + 1
    PNANOVDB_DEREF(io_v).m_v[4] = pnanovdb_read_float(grid_buf, pnanovdb_address_offset(address, PNANOVTT_TRAITS_FLOAT_EXTENDED_TILE_DATA_STRIDE_X)); // i+1, j, k
    PNANOVDB_DEREF(io_v).m_v[5] = pnanovdb_read_float(grid_buf, pnanovdb_address_offset(address, PNANOVTT_TRAITS_FLOAT_EXTENDED_TILE_DATA_STRIDE_X + PNANOVTT_TRAITS_FLOAT_EXTENDED_TILE_DATA_STRIDE_Z)); // i+1, j, k + 1
    PNANOVDB_DEREF(io_v).m_v[6] = pnanovdb_read_float(grid_buf, pnanovdb_address_offset(address, PNANOVTT_TRAITS_FLOAT_EXTENDED_TILE_DATA_STRIDE_X + PNANOVTT_TRAITS_FLOAT_EXTENDED_TILE_DATA_STRIDE_Y)); // i+1, j+1, k
    PNANOVDB_DEREF(io_v).m_v[7] = pnanovdb_read_float(grid_buf, pnanovdb_address_offset(address, PNANOVTT_TRAITS_FLOAT_EXTENDED_TILE_DATA_STRIDE_X + PNANOVTT_TRAITS_FLOAT_EXTENDED_TILE_DATA_STRIDE_Y + PNANOVTT_TRAITS_FLOAT_EXTENDED_TILE_DATA_STRIDE_Z)); // i+1, j+1, k + 1
}

PNANOVDB_FORCE_INLINE float pnanovtt_lerp_float(float a, float b, float w)
{
    return a + w * (b - a);
}

PNANOVDB_FORCE_INLINE float pnanovtt_trilerp_float(PNANOVDB_IN(pnanovdb_vec3_t) in_xyz, PNANOVDB_IN(pnanovtt_stencil_t) in_v)
{
    return pnanovtt_lerp_float(
        pnanovtt_lerp_float(
            pnanovtt_lerp_float(PNANOVDB_DEREF(in_v).m_v[0], PNANOVDB_DEREF(in_v).m_v[1], PNANOVDB_DEREF(in_xyz).z),
            pnanovtt_lerp_float(PNANOVDB_DEREF(in_v).m_v[2], PNANOVDB_DEREF(in_v).m_v[3], PNANOVDB_DEREF(in_xyz).z),
            PNANOVDB_DEREF(in_xyz).y),
        pnanovtt_lerp_float(
            pnanovtt_lerp_float(PNANOVDB_DEREF(in_v).m_v[4], PNANOVDB_DEREF(in_v).m_v[5], PNANOVDB_DEREF(in_xyz).z),
            pnanovtt_lerp_float(PNANOVDB_DEREF(in_v).m_v[6], PNANOVDB_DEREF(in_v).m_v[7], PNANOVDB_DEREF(in_xyz).z),
            PNANOVDB_DEREF(in_xyz).y),
        PNANOVDB_DEREF(in_xyz).x);
}

PNANOVDB_FORCE_INLINE float pnanovtt_sample_float_at_depth_v2(pnanovdb_buf_t buf, pnanovdb_buf_t grid_buf, PNANOVDB_INOUT(pnanovtt_readaccessor_t) acc, pnanovdb_vec3_t in_xyz, pnanovdb_int32_t in_depth)
{
    pnanovdb_vec3_t ijk_out;
    float factor;
    pnanovtt_stencil_t v;
    pnanovdb_coord_t ijk;
    ijk.x = pnanovdb_float_to_int32(pnanovdb_floor(in_xyz.x));
    ijk.y = pnanovdb_float_to_int32(pnanovdb_floor(in_xyz.y));
    ijk.z = pnanovdb_float_to_int32(pnanovdb_floor(in_xyz.z));
    pnanovtt_stencil_float(buf, grid_buf, acc, PNANOVDB_REF(ijk), PNANOVDB_REF(v), PNANOVDB_REF(ijk_out), PNANOVDB_REF(factor), in_depth);
    // Account for varying voxelsize in the uvw offset
    in_xyz.x = (in_xyz.x - ijk_out.x) * factor;
    in_xyz.y = (in_xyz.y - ijk_out.y) * factor;
    in_xyz.z = (in_xyz.z - ijk_out.z) * factor;
    // Appears to be faster to check for found here than always calling sample.
    return PNANOVDB_DEREF(acc).m_found
        ? pnanovtt_trilerp_float(PNANOVDB_REF(in_xyz), PNANOVDB_REF(v))
        : v.m_v[0];
}

PNANOVDB_FORCE_INLINE float pnanovtt_sample_float_v2(pnanovdb_buf_t buf, pnanovdb_buf_t grid_buf, PNANOVDB_INOUT(pnanovtt_readaccessor_t) acc, pnanovdb_vec3_t in_xyz)
{
    return pnanovtt_sample_float_at_depth_v2(buf, grid_buf, acc, in_xyz, PNANOVTT_TRAITS_MAX_DEPTH);
}

PNANOVDB_FORCE_INLINE float pnanovtt_sample_float_at_depth(pnanovdb_buf_t buf, PNANOVDB_INOUT(pnanovtt_readaccessor_t) acc, pnanovdb_vec3_t in_xyz, pnanovdb_int32_t in_depth)
{
    return pnanovtt_sample_float_at_depth_v2(buf, buf, acc, in_xyz, in_depth);
}

PNANOVDB_FORCE_INLINE float pnanovtt_sample_float(pnanovdb_buf_t buf, PNANOVDB_INOUT(pnanovtt_readaccessor_t) acc, pnanovdb_vec3_t in_xyz)
{
    return pnanovtt_sample_float_at_depth_v2(buf, buf, acc, in_xyz, PNANOVTT_TRAITS_MAX_DEPTH);
}

PNANOVDB_FORCE_INLINE float pnanovtt_sample_float_stencil(PNANOVDB_IN(pnanovdb_vec3_t) in_xyz, PNANOVDB_IN(pnanovdb_vec3_t) in_ijk, const float in_factor, PNANOVDB_IN(pnanovtt_stencil_t) in_stencil)
{
    pnanovdb_vec3_t frac;
    frac.x = (PNANOVDB_DEREF(in_xyz).x - PNANOVDB_DEREF(in_ijk).x) * in_factor;
    frac.y = (PNANOVDB_DEREF(in_xyz).y - PNANOVDB_DEREF(in_ijk).y) * in_factor;
    frac.z = (PNANOVDB_DEREF(in_xyz).z - PNANOVDB_DEREF(in_ijk).z) * in_factor;
    return pnanovtt_trilerp_float(PNANOVDB_REF(frac), in_stencil);
}

// ------------------------------------------------ Vec3 Sampling -----------------------------------------------------------

struct pnanovtt_stencil_vec3_t
{
    pnanovdb_vec3_t m_v[8];
};
PNANOVDB_STRUCT_TYPEDEF(pnanovtt_stencil_vec3_t)

PNANOVDB_FORCE_INLINE void pnanovtt_stencil_from_base_address_vec3(pnanovdb_buf_t grid_buf,
    pnanovdb_address_t address,
    PNANOVDB_INOUT(pnanovtt_stencil_vec3_t) io_v)
{
    // In this case all values can be looked up from the extended tile itself
    PNANOVDB_DEREF(io_v).m_v[0] = pnanovtt_read_vec3(grid_buf, address); // i, j, k
    PNANOVDB_DEREF(io_v).m_v[1] = pnanovtt_read_vec3(grid_buf, pnanovdb_address_offset(address, 3 * (PNANOVTT_TRAITS_FLOAT_EXTENDED_TILE_DATA_STRIDE_Z))); // i, j, k + 1
    PNANOVDB_DEREF(io_v).m_v[2] = pnanovtt_read_vec3(grid_buf, pnanovdb_address_offset(address, 3 * (PNANOVTT_TRAITS_FLOAT_EXTENDED_TILE_DATA_STRIDE_Y))); // i, j+1, k
    PNANOVDB_DEREF(io_v).m_v[3] = pnanovtt_read_vec3(grid_buf, pnanovdb_address_offset(address, 3 * (PNANOVTT_TRAITS_FLOAT_EXTENDED_TILE_DATA_STRIDE_Y + PNANOVTT_TRAITS_FLOAT_EXTENDED_TILE_DATA_STRIDE_Z))); // i, j+1, k + 1
    PNANOVDB_DEREF(io_v).m_v[4] = pnanovtt_read_vec3(grid_buf, pnanovdb_address_offset(address, 3 * (PNANOVTT_TRAITS_FLOAT_EXTENDED_TILE_DATA_STRIDE_X))); // i+1, j, k
    PNANOVDB_DEREF(io_v).m_v[5] = pnanovtt_read_vec3(grid_buf, pnanovdb_address_offset(address, 3 * (PNANOVTT_TRAITS_FLOAT_EXTENDED_TILE_DATA_STRIDE_X + PNANOVTT_TRAITS_FLOAT_EXTENDED_TILE_DATA_STRIDE_Z))); // i+1, j, k + 1
    PNANOVDB_DEREF(io_v).m_v[6] = pnanovtt_read_vec3(grid_buf, pnanovdb_address_offset(address, 3 * (PNANOVTT_TRAITS_FLOAT_EXTENDED_TILE_DATA_STRIDE_X + PNANOVTT_TRAITS_FLOAT_EXTENDED_TILE_DATA_STRIDE_Y))); // i+1, j+1, k
    PNANOVDB_DEREF(io_v).m_v[7] = pnanovtt_read_vec3(grid_buf, pnanovdb_address_offset(address, 3 * (PNANOVTT_TRAITS_FLOAT_EXTENDED_TILE_DATA_STRIDE_X + PNANOVTT_TRAITS_FLOAT_EXTENDED_TILE_DATA_STRIDE_Y + PNANOVTT_TRAITS_FLOAT_EXTENDED_TILE_DATA_STRIDE_Z))); // i+1, j+1, k + 1
}

PNANOVDB_FORCE_INLINE void pnanovtt_stencil_vec3(pnanovdb_buf_t buf,
    pnanovdb_buf_t grid_buf,
    PNANOVDB_INOUT(pnanovtt_readaccessor_t) acc,
    PNANOVDB_IN(pnanovdb_coord_t) in_ijk,
    PNANOVDB_INOUT(pnanovtt_stencil_vec3_t) io_v,
    PNANOVDB_INOUT(pnanovdb_vec3_t) io_vijk,
    PNANOVDB_INOUT(float) io_factor,
    pnanovdb_int32_t in_depth) {
    pnanovdb_address_t address = pnanovtt_readaccessor_get_value_address_at_depth_v2(PNANOVDB_GRID_TYPE_VEC3F, buf, grid_buf, acc, in_ijk, in_depth);
    PNANOVDB_DEREF(io_v).m_v[0] = pnanovtt_read_vec3(grid_buf, address); // i, j, k

    if (!PNANOVDB_DEREF(acc).m_found)
        return;

    // Account for varying voxelsize in the uvw offset
    pnanovdb_coord_t tijk = PNANOVDB_DEREF(acc).m_tijk;
    int              shift = PNANOVDB_DEREF(acc).m_vw;
    PNANOVDB_DEREF(io_factor) = 2.f / (1 << (shift + 1));
    PNANOVDB_DEREF(io_vijk).x = (tijk.x + (PNANOVDB_DEREF(acc).m_ijk.x << shift));
    PNANOVDB_DEREF(io_vijk).y = (tijk.y + (PNANOVDB_DEREF(acc).m_ijk.y << shift));
    PNANOVDB_DEREF(io_vijk).z = (tijk.z + (PNANOVDB_DEREF(acc).m_ijk.z << shift));

    // In this case all values can be looked up from the extended tile itself
    PNANOVDB_DEREF(io_v).m_v[1] = pnanovtt_read_vec3(grid_buf, pnanovdb_address_offset(address, 3 * (PNANOVTT_TRAITS_FLOAT_EXTENDED_TILE_DATA_STRIDE_Z))); // i, j, k + 1
    PNANOVDB_DEREF(io_v).m_v[2] = pnanovtt_read_vec3(grid_buf, pnanovdb_address_offset(address, 3 * (PNANOVTT_TRAITS_FLOAT_EXTENDED_TILE_DATA_STRIDE_Y))); // i, j+1, k
    PNANOVDB_DEREF(io_v).m_v[3] = pnanovtt_read_vec3(grid_buf, pnanovdb_address_offset(address, 3 * (PNANOVTT_TRAITS_FLOAT_EXTENDED_TILE_DATA_STRIDE_Y + PNANOVTT_TRAITS_FLOAT_EXTENDED_TILE_DATA_STRIDE_Z))); // i, j+1, k + 1
    PNANOVDB_DEREF(io_v).m_v[4] = pnanovtt_read_vec3(grid_buf, pnanovdb_address_offset(address, 3 * (PNANOVTT_TRAITS_FLOAT_EXTENDED_TILE_DATA_STRIDE_X))); // i+1, j, k
    PNANOVDB_DEREF(io_v).m_v[5] = pnanovtt_read_vec3(grid_buf, pnanovdb_address_offset(address, 3 * (PNANOVTT_TRAITS_FLOAT_EXTENDED_TILE_DATA_STRIDE_X + PNANOVTT_TRAITS_FLOAT_EXTENDED_TILE_DATA_STRIDE_Z))); // i+1, j, k + 1
    PNANOVDB_DEREF(io_v).m_v[6] = pnanovtt_read_vec3(grid_buf, pnanovdb_address_offset(address, 3 * (PNANOVTT_TRAITS_FLOAT_EXTENDED_TILE_DATA_STRIDE_X + PNANOVTT_TRAITS_FLOAT_EXTENDED_TILE_DATA_STRIDE_Y))); // i+1, j+1, k
    PNANOVDB_DEREF(io_v).m_v[7] = pnanovtt_read_vec3(grid_buf, pnanovdb_address_offset(address, 3 * (PNANOVTT_TRAITS_FLOAT_EXTENDED_TILE_DATA_STRIDE_X + PNANOVTT_TRAITS_FLOAT_EXTENDED_TILE_DATA_STRIDE_Y + PNANOVTT_TRAITS_FLOAT_EXTENDED_TILE_DATA_STRIDE_Z))); // i+1, j+1, k + 1
}

PNANOVDB_FORCE_INLINE pnanovdb_vec3_t pnanovtt_lerp_vec3(pnanovdb_vec3_t a, pnanovdb_vec3_t b, float w)
{
    pnanovdb_vec3_t res;
    res.x = a.x + w * (b.x - a.x);
    res.y = a.y + w * (b.y - a.y);
    res.z = a.z + w * (b.z - a.z);
    return res;
}

PNANOVDB_FORCE_INLINE pnanovdb_vec3_t pnanovtt_trilerp_vec3(PNANOVDB_IN(pnanovdb_vec3_t) in_xyz, PNANOVDB_IN(pnanovtt_stencil_vec3_t) in_v)
{
    return pnanovtt_lerp_vec3(
        pnanovtt_lerp_vec3(
            pnanovtt_lerp_vec3(PNANOVDB_DEREF(in_v).m_v[0], PNANOVDB_DEREF(in_v).m_v[1], PNANOVDB_DEREF(in_xyz).z),
            pnanovtt_lerp_vec3(PNANOVDB_DEREF(in_v).m_v[2], PNANOVDB_DEREF(in_v).m_v[3], PNANOVDB_DEREF(in_xyz).z),
            PNANOVDB_DEREF(in_xyz).y),
        pnanovtt_lerp_vec3(
            pnanovtt_lerp_vec3(PNANOVDB_DEREF(in_v).m_v[4], PNANOVDB_DEREF(in_v).m_v[5], PNANOVDB_DEREF(in_xyz).z),
            pnanovtt_lerp_vec3(PNANOVDB_DEREF(in_v).m_v[6], PNANOVDB_DEREF(in_v).m_v[7], PNANOVDB_DEREF(in_xyz).z),
            PNANOVDB_DEREF(in_xyz).y),
        PNANOVDB_DEREF(in_xyz).x);
}

PNANOVDB_FORCE_INLINE pnanovdb_vec3_t pnanovtt_sample_vec3_at_depth_v2(pnanovdb_buf_t buf, pnanovdb_buf_t grid_buf, PNANOVDB_INOUT(pnanovtt_readaccessor_t) acc, pnanovdb_vec3_t in_xyz, pnanovdb_int32_t in_depth)
{
    pnanovdb_vec3_t ijk_out;
    float factor;
    pnanovtt_stencil_vec3_t v;
    pnanovdb_coord_t ijk;
    ijk.x = pnanovdb_float_to_int32(pnanovdb_floor(in_xyz.x));
    ijk.y = pnanovdb_float_to_int32(pnanovdb_floor(in_xyz.y));
    ijk.z = pnanovdb_float_to_int32(pnanovdb_floor(in_xyz.z));
    pnanovtt_stencil_vec3(buf, grid_buf, acc, PNANOVDB_REF(ijk), PNANOVDB_REF(v), PNANOVDB_REF(ijk_out), PNANOVDB_REF(factor), in_depth);
    // Account for varying voxelsize in the uvw offset
    in_xyz.x = (in_xyz.x - ijk_out.x) * factor;
    in_xyz.y = (in_xyz.y - ijk_out.y) * factor;
    in_xyz.z = (in_xyz.z - ijk_out.z) * factor;
    // Appears to be faster to check for found here than always calling sample.
    return PNANOVDB_DEREF(acc).m_found
        ? pnanovtt_trilerp_vec3(PNANOVDB_REF(in_xyz), PNANOVDB_REF(v))
        : v.m_v[0];
}

PNANOVDB_FORCE_INLINE pnanovdb_vec3_t pnanovtt_sample_vec3_v2(pnanovdb_buf_t buf, pnanovdb_buf_t grid_buf, PNANOVDB_INOUT(pnanovtt_readaccessor_t) acc, pnanovdb_vec3_t in_xyz)
{
    return pnanovtt_sample_vec3_at_depth_v2(buf, grid_buf, acc, in_xyz, PNANOVTT_TRAITS_MAX_DEPTH);
}

PNANOVDB_FORCE_INLINE pnanovdb_vec3_t pnanovtt_sample_vec3_at_depth(pnanovdb_buf_t buf, PNANOVDB_INOUT(pnanovtt_readaccessor_t) acc, pnanovdb_vec3_t in_xyz, pnanovdb_int32_t in_depth)
{
    return pnanovtt_sample_vec3_at_depth_v2(buf, buf, acc, in_xyz, in_depth);
}

PNANOVDB_FORCE_INLINE pnanovdb_vec3_t pnanovtt_sample_vec3(pnanovdb_buf_t buf, PNANOVDB_INOUT(pnanovtt_readaccessor_t) acc, pnanovdb_vec3_t in_xyz)
{
    return pnanovtt_sample_vec3_at_depth_v2(buf, buf, acc, in_xyz, PNANOVTT_TRAITS_MAX_DEPTH);
}

// ------------------------------------------------ Cached Float Sampling -----------------------------------------------------------


struct pnanovtt_cached_sampler_t
{
    pnanovdb_vec3_t m_ijk;
    float m_factor;
    pnanovtt_readaccessor_t m_acc;
    pnanovtt_stencil_t m_stencil;
};
PNANOVDB_STRUCT_TYPEDEF(pnanovtt_cached_sampler_t)


PNANOVDB_FORCE_INLINE void pnanovtt_cached_sampler_init(pnanovdb_buf_t buf, PNANOVDB_INOUT(pnanovtt_cached_sampler_t) vtt_sampler, pnanovdb_grid_handle_t grid, pnanovtt_tree_handle_t tree)
{
    PNANOVDB_DEREF(vtt_sampler).m_ijk.x = 0;
    PNANOVDB_DEREF(vtt_sampler).m_ijk.y = 0;
    PNANOVDB_DEREF(vtt_sampler).m_ijk.z = 0;
    PNANOVDB_DEREF(vtt_sampler).m_factor = 0;
    pnanovtt_readaccessor_init(
        buf, PNANOVDB_REF(PNANOVDB_DEREF(vtt_sampler).m_acc), grid, tree);
    for (int i = 0; i < 8; i++)
        PNANOVDB_DEREF(vtt_sampler).m_stencil.m_v[i] = 0;

}

PNANOVDB_FORCE_INLINE float pnanovtt_cached_sample_float_at_depth_v2(pnanovdb_buf_t buf, pnanovdb_buf_t grid_buf, PNANOVDB_INOUT(pnanovtt_cached_sampler_t) vtt_sampler, pnanovdb_vec3_t in_xyz, pnanovdb_int32_t in_depth)
{
    pnanovdb_coord_t ijk;
    ijk.x = pnanovdb_float_to_int32(pnanovdb_floor(in_xyz.x));
    ijk.y = pnanovdb_float_to_int32(pnanovdb_floor(in_xyz.y));
    ijk.z = pnanovdb_float_to_int32(pnanovdb_floor(in_xyz.z));
    pnanovtt_stencil_float(buf, grid_buf, PNANOVDB_REF(PNANOVDB_DEREF(vtt_sampler).m_acc), PNANOVDB_REF(ijk), PNANOVDB_REF(PNANOVDB_DEREF(vtt_sampler).m_stencil), PNANOVDB_REF(PNANOVDB_DEREF(vtt_sampler).m_ijk), PNANOVDB_REF(PNANOVDB_DEREF(vtt_sampler).m_factor), in_depth);
    // Account for varying voxelsize in the uvw offset
    in_xyz.x = (in_xyz.x - PNANOVDB_DEREF(vtt_sampler).m_ijk.x) * PNANOVDB_DEREF(vtt_sampler).m_factor;
    in_xyz.y = (in_xyz.y - PNANOVDB_DEREF(vtt_sampler).m_ijk.y) * PNANOVDB_DEREF(vtt_sampler).m_factor;
    in_xyz.z = (in_xyz.z - PNANOVDB_DEREF(vtt_sampler).m_ijk.z) * PNANOVDB_DEREF(vtt_sampler).m_factor;
    // Appears to be faster to check for found here than always calling sample.
    return PNANOVDB_DEREF(vtt_sampler).m_acc.m_found
        ? pnanovtt_trilerp_float(PNANOVDB_REF(in_xyz), PNANOVDB_REF(PNANOVDB_DEREF(vtt_sampler).m_stencil))
        : PNANOVDB_DEREF(vtt_sampler).m_stencil.m_v[0];
}

PNANOVDB_FORCE_INLINE float pnanovtt_cached_sample_float_v2(pnanovdb_buf_t buf, pnanovdb_buf_t grid_buf, PNANOVDB_INOUT(pnanovtt_cached_sampler_t) vtt_sampler, pnanovdb_vec3_t in_xyz)
{
    return pnanovtt_cached_sample_float_at_depth_v2(buf, grid_buf, vtt_sampler, in_xyz, PNANOVTT_TRAITS_MAX_DEPTH);
}

PNANOVDB_FORCE_INLINE float pnanovtt_cached_sample_float_at_depth(pnanovdb_buf_t buf, PNANOVDB_INOUT(pnanovtt_cached_sampler_t) vtt_sampler, pnanovdb_vec3_t in_xyz, pnanovdb_int32_t in_depth)
{
    return pnanovtt_cached_sample_float_at_depth_v2(buf, buf, vtt_sampler, in_xyz, in_depth);
}

PNANOVDB_FORCE_INLINE float pnanovtt_cached_sample_float(pnanovdb_buf_t buf, PNANOVDB_INOUT(pnanovtt_cached_sampler_t) vtt_sampler, pnanovdb_vec3_t in_xyz)
{
    return pnanovtt_cached_sample_float_at_depth_v2(buf, buf, vtt_sampler, in_xyz, PNANOVTT_TRAITS_MAX_DEPTH);
}

PNANOVDB_FORCE_INLINE float pnanovtt_cached_sample_cached_voxel_float(PNANOVDB_IN(pnanovtt_cached_sampler_t) vtt_sampler, pnanovdb_vec3_t in_xyz)
{
    // Account for varying voxelsize in the uvw offset
    in_xyz.x = (in_xyz.x - PNANOVDB_DEREF(vtt_sampler).m_ijk.x) * PNANOVDB_DEREF(vtt_sampler).m_factor;
    in_xyz.y = (in_xyz.y - PNANOVDB_DEREF(vtt_sampler).m_ijk.y) * PNANOVDB_DEREF(vtt_sampler).m_factor;
    in_xyz.z = (in_xyz.z - PNANOVDB_DEREF(vtt_sampler).m_ijk.z) * PNANOVDB_DEREF(vtt_sampler).m_factor;
    return pnanovtt_trilerp_float(PNANOVDB_REF(in_xyz), PNANOVDB_REF(PNANOVDB_DEREF(vtt_sampler).m_stencil));
}

PNANOVDB_FORCE_INLINE pnanovdb_bool_t pnanovtt_cached_sample_zero_crossing_float(PNANOVDB_IN(pnanovtt_cached_sampler_t) vtt_sampler)
{
    pnanovdb_bool_t less = PNANOVDB_DEREF(vtt_sampler).m_stencil.m_v[0] < 0;
    return (less != (PNANOVDB_DEREF(vtt_sampler).m_stencil.m_v[1] < 0)) ||
        (less != (PNANOVDB_DEREF(vtt_sampler).m_stencil.m_v[3] < 0)) ||
        (less != (PNANOVDB_DEREF(vtt_sampler).m_stencil.m_v[2] < 0)) ||
        (less != (PNANOVDB_DEREF(vtt_sampler).m_stencil.m_v[4] < 0)) ||
        (less != (PNANOVDB_DEREF(vtt_sampler).m_stencil.m_v[5] < 0)) ||
        (less != (PNANOVDB_DEREF(vtt_sampler).m_stencil.m_v[7] < 0)) ||
        (less != (PNANOVDB_DEREF(vtt_sampler).m_stencil.m_v[6] < 0));
}

// ------------------------------------------------ HDDA (extensions to PNanoVDB's HDDA) -----------------------------------------------------------

#ifdef PNANOVDB_HDDA

#define PNANOVTT_HDDA_DELTA 0.0001f
#define PNANOVTT_HDDA_MAX_BISECTIONS 10

PNANOVDB_FORCE_INLINE pnanovdb_bool_t pnanovtt_hdda_zero_crossing(
    pnanovdb_grid_type_t grid_type,
    pnanovdb_buf_t buf,
    PNANOVDB_INOUT(pnanovtt_readaccessor_t) acc,
    PNANOVDB_IN(pnanovdb_vec3_t) origin, float tmin,
    PNANOVDB_IN(pnanovdb_vec3_t) direction, float tmax,
    PNANOVDB_INOUT(float) thit,
    PNANOVDB_INOUT(float) v,
    PNANOVDB_INOUT(pnanovdb_coord_t) ijk,
    pnanovdb_int32_t max_depth
)
{
    pnanovdb_coord_t bbox_min = pnanovtt_tree_get_bbox_min(buf, PNANOVDB_DEREF(acc).m_tree);
    pnanovdb_coord_t bbox_max = pnanovtt_tree_get_bbox_max(buf, PNANOVDB_DEREF(acc).m_tree);
    pnanovdb_vec3_t bbox_minf = pnanovdb_coord_to_vec3(bbox_min);
    pnanovdb_vec3_t bbox_maxf = pnanovdb_coord_to_vec3(pnanovdb_coord_add(bbox_max, pnanovdb_coord_uniform(1)));

    pnanovdb_bool_t hit = pnanovdb_hdda_ray_clip(PNANOVDB_REF(bbox_minf), PNANOVDB_REF(bbox_maxf), origin, PNANOVDB_REF(tmin), direction, PNANOVDB_REF(tmax));
    if (!hit || tmax > 1.0e20f)
    {
        return PNANOVDB_FALSE;
    }

    // First we step delta backwards along the ray to make sure we end up outside the first voxel in the bbox.
    tmin -= PNANOVTT_HDDA_DELTA;
    pnanovdb_vec3_t pos = pnanovdb_hdda_ray_start(origin, tmin, direction);
    PNANOVDB_DEREF(ijk) = pnanovdb_hdda_pos_to_ijk(PNANOVDB_REF(pos));

    pnanovdb_address_t address = pnanovtt_readaccessor_get_value_address_at_depth(grid_type, buf, acc, ijk, max_depth);
    float v0 = pnanovdb_read_float(buf, address);
    pnanovdb_int32_t dim = pnanovdb_uint32_as_int32(pnanovtt_readaccessor_get_dim2(acc));
    pnanovdb_hdda_t hdda;
    pnanovdb_hdda_init(PNANOVDB_REF(hdda), origin, tmin, direction, tmax, dim);

    // We leave out the call to is_active as it always returns true inside the bbox of the level set.
    do {
        PNANOVDB_DEREF(ijk) = hdda.voxel;
        address = pnanovtt_readaccessor_get_value_address_at_depth(
            grid_type, buf, acc, ijk, max_depth);
        float v1 = pnanovdb_read_float(buf, address);
        // Get dim without tree access.
        dim = pnanovdb_uint32_as_int32(pnanovtt_readaccessor_get_dim2(acc));
        // We need to update the dim and adjust the voxel before we retrieve the value
        // as we could have entered a finer voxel.
        pnanovdb_hdda_update(PNANOVDB_REF(hdda), origin, direction, dim);
        if (hdda.voxel.x != PNANOVDB_DEREF(ijk).x ||
            hdda.voxel.y != PNANOVDB_DEREF(ijk).y ||
            hdda.voxel.z != PNANOVDB_DEREF(ijk).z) {
            PNANOVDB_DEREF(ijk) = hdda.voxel;
            address = pnanovtt_readaccessor_get_value_address_at_depth(
                grid_type, buf, acc, ijk, max_depth);
            PNANOVDB_DEREF(v) = pnanovdb_read_float(buf, address);
        }
        else {
            PNANOVDB_DEREF(v) = v1;
        }
        // The second check avoids false positives for the zero crossing at resolution jumps.
        if (PNANOVDB_DEREF(v) * v0 <= 0.f && PNANOVDB_DEREF(v) * v1 >= 0.f) {
            PNANOVDB_DEREF(thit) = hdda.tmin;
            return PNANOVDB_TRUE;
        }
    } while (pnanovdb_hdda_step(PNANOVDB_REF(hdda)));

    return PNANOVDB_FALSE;
}

PNANOVDB_FORCE_INLINE float pnanovtt_hdda_next(PNANOVDB_IN(pnanovdb_hdda_t) hdda) {
    return pnanovdb_min(
        PNANOVDB_DEREF(hdda).tmax,
        pnanovdb_min(PNANOVDB_DEREF(hdda).next.x,
            pnanovdb_min(PNANOVDB_DEREF(hdda).next.y,
                PNANOVDB_DEREF(hdda).next.z)));
}

PNANOVDB_FORCE_INLINE pnanovdb_bool_t
pnanovtt_hdda_update(PNANOVDB_INOUT(pnanovdb_hdda_t) hdda,
    PNANOVDB_IN(pnanovdb_vec3_t) direction,
    int dim,
    PNANOVDB_IN(pnanovdb_vec3_t) pos,
    pnanovdb_coord_t voxel)
{
    voxel.x = voxel.x & (~(dim - 1));
    voxel.y = voxel.y & (~(dim - 1));
    voxel.z = voxel.z & (~(dim - 1));

    if (voxel.x == PNANOVDB_DEREF(hdda).voxel.x &&
        voxel.y == PNANOVDB_DEREF(hdda).voxel.y &&
        voxel.z == PNANOVDB_DEREF(hdda).voxel.z &&
        PNANOVDB_DEREF(hdda).dim == dim)
        return PNANOVDB_FALSE;

    PNANOVDB_DEREF(hdda).voxel.x = voxel.x;
    PNANOVDB_DEREF(hdda).voxel.y = voxel.y;
    PNANOVDB_DEREF(hdda).voxel.z = voxel.z;
    PNANOVDB_DEREF(hdda).dim = dim;

    pnanovdb_vec3_t dir_inv = pnanovdb_vec3_div(pnanovdb_vec3_uniform(1.f), PNANOVDB_DEREF(direction));

    if (PNANOVDB_DEREF(hdda).step.x != 0)
    {
        PNANOVDB_DEREF(hdda).next.x = PNANOVDB_DEREF(hdda).tmin + (PNANOVDB_DEREF(hdda).voxel.x - PNANOVDB_DEREF(pos).x) * dir_inv.x;
        if (PNANOVDB_DEREF(hdda).step.x > 0)
        {
            PNANOVDB_DEREF(hdda).next.x += dim * dir_inv.x;
        }
    }
    if (PNANOVDB_DEREF(hdda).step.y != 0)
    {
        PNANOVDB_DEREF(hdda).next.y = PNANOVDB_DEREF(hdda).tmin + (PNANOVDB_DEREF(hdda).voxel.y - PNANOVDB_DEREF(pos).y) * dir_inv.y;
        if (PNANOVDB_DEREF(hdda).step.y > 0)
        {
            PNANOVDB_DEREF(hdda).next.y += dim * dir_inv.y;
        }
    }
    if (PNANOVDB_DEREF(hdda).step.z != 0)
    {
        PNANOVDB_DEREF(hdda).next.z = PNANOVDB_DEREF(hdda).tmin + (PNANOVDB_DEREF(hdda).voxel.z - PNANOVDB_DEREF(pos).z) * dir_inv.z;
        if (PNANOVDB_DEREF(hdda).step.z > 0)
        {
            PNANOVDB_DEREF(hdda).next.z += dim * dir_inv.z;
        }
    }

    return PNANOVDB_TRUE;
}


// Computes the ray-levelset intersection with sub-voxel precision.
// The input ray (origin, direction) as well as the output time (thit)
// and position (poshit) are in index space.
PNANOVDB_FORCE_INLINE pnanovdb_bool_t pnanovtt_hdda_zero_crossing_sub_voxel_quadratic(
    pnanovdb_buf_t       buf,
    PNANOVDB_INOUT(pnanovtt_cached_sampler_t) vtt_sampler,
    PNANOVDB_IN(pnanovdb_vec3_t) origin,
    float tmin,
    PNANOVDB_IN(pnanovdb_vec3_t) direction,
    float tmax,
    float iso_value,
    PNANOVDB_INOUT(float) thit,
    PNANOVDB_INOUT(float) vhit,
    PNANOVDB_INOUT(pnanovdb_vec3_t) poshit,
    pnanovdb_int32_t max_depth)
{
    pnanovdb_coord_t bbox_min = pnanovtt_tree_get_bbox_min(buf, PNANOVDB_DEREF(vtt_sampler).m_acc.m_tree);
    pnanovdb_coord_t bbox_max = pnanovtt_tree_get_bbox_max(buf, PNANOVDB_DEREF(vtt_sampler).m_acc.m_tree);

    pnanovdb_vec3_t bbox_minf = pnanovdb_coord_to_vec3(bbox_min);
    pnanovdb_vec3_t bbox_maxf = pnanovdb_coord_to_vec3(pnanovdb_coord_add(bbox_max, pnanovdb_coord_uniform(1)));

    pnanovdb_bool_t hit = pnanovdb_hdda_ray_clip(PNANOVDB_REF(bbox_minf), PNANOVDB_REF(bbox_maxf), origin, PNANOVDB_REF(tmin), direction, PNANOVDB_REF(tmax));
    if (!hit || tmax > 1.0e20f)
    {
        return PNANOVDB_FALSE;
    }

    // First we step PNANOVTT_HDDA_DELTA along the ray to make sure we end up inside the voxel we're going to step through.
    // Then we can just use dim from that voxel. We also ensure that the first voxel is always inside the bbox.
    tmin += PNANOVTT_HDDA_DELTA;

    // Setup HDDA
    // The assumption here now is that we always start inside a voxel inside the bbox and therefore
    // the hdda.dim and hdda.voxel, hdda.tmin and pnanovtt_hdda_next(PNANOVDB_REF(hdda)) will always correspond to the current voxel.
    pnanovdb_hdda_t hdda;
    pnanovdb_vec3_t pos = pnanovdb_hdda_ray_start(origin, tmin, direction);

    hdda.voxel = pnanovdb_hdda_pos_to_ijk(PNANOVDB_REF(pos));
    pnanovdb_address_t address =
        pnanovtt_readaccessor_get_value_address_at_depth(
            PNANOVDB_GRID_TYPE_FLOAT, buf,
            PNANOVDB_REF(PNANOVDB_DEREF(vtt_sampler).m_acc), PNANOVDB_REF(hdda.voxel),
            max_depth);

    pnanovdb_hdda_init(PNANOVDB_REF(hdda), origin, tmin, direction, tmax,
        pnanovdb_uint32_as_int32(pnanovtt_readaccessor_get_dim2(
            PNANOVDB_REF(PNANOVDB_DEREF(vtt_sampler).m_acc))));

    float pVal = pnanovdb_read_float(buf, address);

    pnanovdb_bool_t inside = (pVal < iso_value);
    float pTmin = tmin;

    do {

        pnanovdb_vec3_t ijkf = pnanovdb_coord_to_vec3(hdda.voxel);
        float dimf = pnanovdb_int32_to_float(hdda.dim);

        // Here we first make sure that the ray start point is clamped to be slightly inside the hdda voxel.
        // This may not be the cause due to numerical roundoff error and can lead to false positives.
        pos = pnanovdb_hdda_ray_start(origin, hdda.tmin, direction);
        pos.x = pnanovdb_min(pnanovdb_max(ijkf.x + PNANOVTT_HDDA_DELTA, pos.x), ijkf.x + dimf);
        pos.y = pnanovdb_min(pnanovdb_max(ijkf.y + PNANOVTT_HDDA_DELTA, pos.y), ijkf.y + dimf);
        pos.z = pnanovdb_min(pnanovdb_max(ijkf.z + PNANOVTT_HDDA_DELTA, pos.z), ijkf.z + dimf);

        // Retrieve stencil for voxel
        pnanovdb_coord_t ijk = pnanovdb_hdda_pos_to_ijk(PNANOVDB_REF(pos));
        pnanovtt_stencil_float(
            buf, buf, PNANOVDB_REF(PNANOVDB_DEREF(vtt_sampler).m_acc),
            PNANOVDB_REF(ijk),
            PNANOVDB_REF(PNANOVDB_DEREF(vtt_sampler).m_stencil),
            PNANOVDB_REF(PNANOVDB_DEREF(vtt_sampler).m_ijk),
            PNANOVDB_REF(PNANOVDB_DEREF(vtt_sampler).m_factor), max_depth);

        // Update dim in hdda.
        pnanovdb_int32_t dim = pnanovdb_uint32_as_int32(pnanovtt_readaccessor_get_dim2(PNANOVDB_REF(PNANOVDB_DEREF(vtt_sampler).m_acc)));
        pnanovtt_hdda_update(PNANOVDB_REF(hdda), direction, dim,
                             PNANOVDB_REF(pos),
                             ijk);
        dimf = pnanovdb_int32_to_float(hdda.dim);
        ijkf = pnanovdb_coord_to_vec3(hdda.voxel);

        // Clamp end point of voxel to be exactly within voxel such that we can re-use
        // the stencil for interpolation.
        tmin = pnanovtt_hdda_next(PNANOVDB_REF(hdda));

        pos = pnanovdb_hdda_ray_start(origin, tmin, direction);
        pos.x = pnanovdb_min(pnanovdb_max(ijkf.x, pos.x), ijkf.x + dimf);
        pos.y = pnanovdb_min(pnanovdb_max(ijkf.y, pos.y), ijkf.y + dimf);
        pos.z = pnanovdb_min(pnanovdb_max(ijkf.z, pos.z), ijkf.z + dimf);

        float pVal2 = ((!PNANOVDB_DEREF(vtt_sampler).m_acc.m_found) ? PNANOVDB_DEREF(vtt_sampler).m_stencil.m_v[0] :
            pnanovtt_cached_sample_cached_voxel_float(vtt_sampler, pos));

        pnanovdb_bool_t inside2 = (pVal2 < iso_value);
        if (PNANOVDB_DEREF(vtt_sampler).m_acc.m_found && (inside != inside2 || pnanovtt_cached_sample_zero_crossing_float(vtt_sampler)))
        {
            // Using the secant method directly here will often not converge
            // because there will be a local minimum or maximum in the level set
            // function along the ray which can prevent features like corners
            // from being captured correctly. Instead we make a quadratic fit
            // to the cubic level set function along the ray and find the
            // closest root of this analytically.

            pnanovdb_bool_t found = PNANOVDB_FALSE;
            float d = 0.f;

            // Fit a quadratic with re-scaling of the equation
            if (pVal <= 0 && pVal >= 0) {
                found = PNANOVDB_TRUE;
            } else {
                // Evaluate at midpoint between a and b.
                pos = pnanovdb_hdda_ray_start(origin, (pTmin + tmin) * 0.5f, direction);
                pos.x = pnanovdb_min(pnanovdb_max(ijkf.x, pos.x), ijkf.x + dimf);
                pos.y = pnanovdb_min(pnanovdb_max(ijkf.y, pos.y), ijkf.y + dimf);
                pos.z = pnanovdb_min(pnanovdb_max(ijkf.z, pos.z), ijkf.z + dimf);
                float pValHalf = pnanovtt_cached_sample_cached_voxel_float(vtt_sampler, pos);

                // normalize function values, so our quadratic is interpolating
                // (1,fHalf,f1)
                float r0    = 1.f / pVal;
                float fHalf = r0 * pValHalf, f1 = r0 * pVal2;

                float t = 4 * fHalf - f1 - 1; // equal to f'(0)+2
                d = t * t - 4 * f1; // note: absolutely has to be >=0 if f1<=0
                                    // (i.e. we saw a sign change)
                if (d >= 0 && (t <= 0 || f1 <= 0)) {
                    d     = sqrt(d) - t;
                    d     = 2.f / (2 + d);
                    found = PNANOVDB_TRUE;
                }
            } // if (pVal <= 0 && pVal >= 0) .. else

            if (found)
            {
                PNANOVDB_DEREF(thit) = d * (tmin - pTmin) + pTmin;
                pos = pnanovdb_hdda_ray_start(origin, PNANOVDB_DEREF(thit), direction);
                pos.x = pnanovdb_min(pnanovdb_max(ijkf.x, pos.x), ijkf.x + dimf);
                pos.y = pnanovdb_min(pnanovdb_max(ijkf.y, pos.y), ijkf.y + dimf);
                pos.z = pnanovdb_min(pnanovdb_max(ijkf.z, pos.z), ijkf.z + dimf);
                PNANOVDB_DEREF(vhit) = pnanovtt_cached_sample_cached_voxel_float(vtt_sampler, pos);
                PNANOVDB_DEREF(poshit) = pos;
                return PNANOVDB_TRUE;
            }
        } // else if (pnanovtt_cached_sample_zero_crossing_float(vtt_sampler))

        pVal = pVal2;
        pTmin = tmin;
    } while (pnanovdb_hdda_step(PNANOVDB_REF(hdda)));

    return PNANOVDB_FALSE;
}


// Computes the ray-levelset intersection with sub-voxel precision.
// The input ray (origin, direction) as well as the output time (thit)
// and position (poshit) are in index space.
PNANOVDB_FORCE_INLINE pnanovdb_bool_t pnanovtt_hdda_zero_crossing_sub_voxel(
    pnanovdb_buf_t       buf,
    PNANOVDB_INOUT(pnanovtt_cached_sampler_t) vtt_sampler,
    PNANOVDB_IN(pnanovdb_vec3_t) origin,
    float tmin,
    PNANOVDB_IN(pnanovdb_vec3_t) direction,
    float tmax,
    float iso_value,
    PNANOVDB_INOUT(float) thit,
    PNANOVDB_INOUT(float) vhit,
    PNANOVDB_INOUT(pnanovdb_vec3_t) poshit,
    pnanovdb_int32_t max_depth)
{
    pnanovdb_coord_t bbox_min = pnanovtt_tree_get_bbox_min(buf, PNANOVDB_DEREF(vtt_sampler).m_acc.m_tree);
    pnanovdb_coord_t bbox_max = pnanovtt_tree_get_bbox_max(buf, PNANOVDB_DEREF(vtt_sampler).m_acc.m_tree);

    pnanovdb_vec3_t bbox_minf = pnanovdb_coord_to_vec3(bbox_min);
    pnanovdb_vec3_t bbox_maxf = pnanovdb_coord_to_vec3(pnanovdb_coord_add(bbox_max, pnanovdb_coord_uniform(1)));

    pnanovdb_bool_t hit = pnanovdb_hdda_ray_clip(PNANOVDB_REF(bbox_minf), PNANOVDB_REF(bbox_maxf), origin, PNANOVDB_REF(tmin), direction, PNANOVDB_REF(tmax));
    if (!hit || tmax > 1.0e20f)
    {
        return PNANOVDB_FALSE;
    }

    // First we step PNANOVTT_HDDA_DELTA along the ray to make sure we end up inside the voxel we're going to step through.
    // Then we can just use dim from that voxel. We also ensure that the first voxel is always inside the bbox.
    tmin += PNANOVTT_HDDA_DELTA;

    // Setup HDDA
    // The assumption here now is that we always start inside a voxel inside the bbox and therefore
    // the hdda.dim and hdda.voxel, hdda.tmin and pnanovtt_hdda_next(PNANOVDB_REF(hdda)) will always correspond to the current voxel.
    pnanovdb_hdda_t hdda;
    pnanovdb_vec3_t pos = pnanovdb_hdda_ray_start(origin, tmin, direction);

    hdda.voxel = pnanovdb_hdda_pos_to_ijk(PNANOVDB_REF(pos));
    pnanovdb_address_t address =
        pnanovtt_readaccessor_get_value_address_at_depth(
            PNANOVDB_GRID_TYPE_FLOAT, buf,
            PNANOVDB_REF(PNANOVDB_DEREF(vtt_sampler).m_acc), PNANOVDB_REF(hdda.voxel),
            max_depth);

    pnanovdb_hdda_init(PNANOVDB_REF(hdda), origin, tmin, direction, tmax,
                       pnanovdb_uint32_as_int32(pnanovtt_readaccessor_get_dim2(
                           PNANOVDB_REF(PNANOVDB_DEREF(vtt_sampler).m_acc))));

    float pVal = pnanovdb_read_float(buf, address);

    pnanovdb_bool_t inside = (pVal < iso_value);
    float pTmin = tmin;

    do {

        pnanovdb_vec3_t ijkf = pnanovdb_coord_to_vec3(hdda.voxel);
        float dimf = pnanovdb_int32_to_float(hdda.dim);

        // Here we first make sure that the ray start point is clamped to be slightly inside the hdda voxel.
        // This may not be the cause due to numerical roundoff error and can lead to false positives.
        pos = pnanovdb_hdda_ray_start(origin, hdda.tmin, direction);
        pos.x = pnanovdb_min(pnanovdb_max(ijkf.x + PNANOVTT_HDDA_DELTA, pos.x), ijkf.x + dimf);
        pos.y = pnanovdb_min(pnanovdb_max(ijkf.y + PNANOVTT_HDDA_DELTA, pos.y), ijkf.y + dimf);
        pos.z = pnanovdb_min(pnanovdb_max(ijkf.z + PNANOVTT_HDDA_DELTA, pos.z), ijkf.z + dimf);

        // Retrieve stencil for voxel
        pnanovdb_coord_t ijk = pnanovdb_hdda_pos_to_ijk(PNANOVDB_REF(pos));
        pnanovtt_stencil_float(
            buf, buf, PNANOVDB_REF(PNANOVDB_DEREF(vtt_sampler).m_acc),
            PNANOVDB_REF(ijk),
            PNANOVDB_REF(PNANOVDB_DEREF(vtt_sampler).m_stencil),
            PNANOVDB_REF(PNANOVDB_DEREF(vtt_sampler).m_ijk),
            PNANOVDB_REF(PNANOVDB_DEREF(vtt_sampler).m_factor), max_depth);

        // Update dim in hdda.
        pnanovdb_int32_t dim = pnanovdb_uint32_as_int32(pnanovtt_readaccessor_get_dim2(PNANOVDB_REF(PNANOVDB_DEREF(vtt_sampler).m_acc)));
        pnanovtt_hdda_update(PNANOVDB_REF(hdda), direction, dim,
            PNANOVDB_REF(pos),
            ijk);
        dimf = pnanovdb_int32_to_float(hdda.dim);
        ijkf = pnanovdb_coord_to_vec3(hdda.voxel);

        // Clamp end point of voxel to be exactly within voxel such that we can re-use
        // the stencil for interpolation.
        tmin = pnanovtt_hdda_next(PNANOVDB_REF(hdda));

        pos = pnanovdb_hdda_ray_start(origin, tmin, direction);
        pos.x = pnanovdb_min(pnanovdb_max(ijkf.x, pos.x), ijkf.x + dimf);
        pos.y = pnanovdb_min(pnanovdb_max(ijkf.y, pos.y), ijkf.y + dimf);
        pos.z = pnanovdb_min(pnanovdb_max(ijkf.z, pos.z), ijkf.z + dimf);

        float pVal2 = ((!PNANOVDB_DEREF(vtt_sampler).m_acc.m_found) ? PNANOVDB_DEREF(vtt_sampler).m_stencil.m_v[0] :
            pnanovtt_cached_sample_cached_voxel_float(vtt_sampler, pos));

        pnanovdb_bool_t inside2 = (pVal2 < iso_value);
        if (inside != inside2)
        {
            // If different sign at endpoints then perform bisection for a fixed number of iterations,
            // and end with linear interpolation to find surface position.
            // To improve the convergence we take the level set phi values into account in the bisection.
            float a = pTmin;
            float b = tmin;
            float phiA = pnanovtt_abs(pVal);
            float phiB = pnanovtt_abs(pVal2);
            pnanovdb_int32_t step2 = 0;

            // Bisection will always converge in a monotonic manner, but (depending on the endpoints)
            // the value at the new points (and hence the final interpolated point)
            // may be further from the zero crossing than the end-point values.
            while (step2 < PNANOVTT_HDDA_MAX_BISECTIONS
                && (b - a) > PNANOVTT_HDDA_DELTA) {
                // Here we utilize the level set values to accelerate the bisection.
                // Traditional bisection:
                //const float m = a + (b - a) * 0.5f;
                // Bisection using level set values:
                float m = a + phiA / (phiA + phiB) * (b - a);
                pos = pnanovdb_hdda_ray_start(origin, m, direction);
                // Because we're in the interior of the interval, we can most likely leave out these checks,
                // but it doesn't appear to make a difference for performance.
                pos.x = pnanovdb_min(pnanovdb_max(ijkf.x, pos.x), ijkf.x + dimf);
                pos.y = pnanovdb_min(pnanovdb_max(ijkf.y, pos.y), ijkf.y + dimf);
                pos.z = pnanovdb_min(pnanovdb_max(ijkf.z, pos.z), ijkf.z + dimf);

                // Use accelerated sampler call.
                // Also ensures sampler does not move outside of voxel.
                pVal2 = pnanovtt_cached_sample_cached_voxel_float(vtt_sampler, pos);
                pnanovdb_bool_t sb = (pVal2 < iso_value);
                if (inside == sb) {
                    a = m;
                    phiA = pnanovtt_abs(pVal2);
                }
                else {
                    b = m;
                    phiB = pnanovtt_abs(pVal2);
                }
                ++step2;
            } // while (...)

            // Compute final intersection by linear interpolation
            // We know that phiA and phiB have different sign.
            PNANOVDB_DEREF(thit) = a + phiA / (phiA + phiB) * (b - a);
            pos = pnanovdb_hdda_ray_start(origin, PNANOVDB_DEREF(thit), direction);
            pos.x = pnanovdb_min(pnanovdb_max(ijkf.x, pos.x), ijkf.x + dimf);
            pos.y = pnanovdb_min(pnanovdb_max(ijkf.y, pos.y), ijkf.y + dimf);
            pos.z = pnanovdb_min(pnanovdb_max(ijkf.z, pos.z), ijkf.z + dimf);
            PNANOVDB_DEREF(vhit) = pnanovtt_cached_sample_cached_voxel_float(vtt_sampler, pos);
            PNANOVDB_DEREF(poshit) = pos;
            return PNANOVDB_TRUE;
        } // else if (inside != inside2)
        else if (PNANOVDB_DEREF(vtt_sampler).m_acc.m_found && pnanovtt_cached_sample_zero_crossing_float(vtt_sampler))
        {
            // Using the secant method directly here will often not converge
            // because there will be a local minimum or maximum in the level set
            // function along the ray which can prevent features like corners
            // from being captured correctly. Instead we make a quadratic fit
            // to the cubic level set function along the ray and find the
            // closest root of this analytically.

            float a = pTmin;
            float b = tmin;
            float phiA = pVal;
            float phiB = pVal2;

            // Evaluate at midpoint between a and b.
            float c = (a + b) * 0.5f;
            pos = pnanovdb_hdda_ray_start(origin, c, direction);
            pos.x = pnanovdb_min(pnanovdb_max(ijkf.x, pos.x), ijkf.x + dimf);
            pos.y = pnanovdb_min(pnanovdb_max(ijkf.y, pos.y), ijkf.y + dimf);
            pos.z = pnanovdb_min(pnanovdb_max(ijkf.z, pos.z), ijkf.z + dimf);
            float phiC = pnanovtt_cached_sample_cached_voxel_float(vtt_sampler, pos);

            // Find the coefficients of the quadratic assuming [a;b] == [0;1] (we adjust for that later).
            float pa = 2 * phiA - 4 * phiC + 2 * phiB;
            float pb = -3 * phiA + 4 * phiC - phiB;
            float pc = phiA;

            float d = pb * pb - 4 * pa * pc;

            // We only treat the case where d>0 and we have two roots.
            // We need to check if the roots are in the interval (a;b) and if so,
            // pick the one closest to a.
            if (d > 0 && pa > 0) {
                // We have an intersection
                d = sqrt(d);
                float r1 = (-pb - d) / (2 * pa);
                float r2 = (-pb + d) / (2 * pa);
                r1 = ((r1 < 0) ? 2 : r1);
                r2 = ((r2 < 0) ? 2 : r2);
                d = pnanovdb_min(r1, r2);
                if (d >= 0.f && d <= 1.f) {
                    PNANOVDB_DEREF(thit) = d * (b - a) + a;
                    pos = pnanovdb_hdda_ray_start(origin, PNANOVDB_DEREF(thit), direction);
                    pos.x = pnanovdb_min(pnanovdb_max(ijkf.x, pos.x), ijkf.x + dimf);
                    pos.y = pnanovdb_min(pnanovdb_max(ijkf.y, pos.y), ijkf.y + dimf);
                    pos.z = pnanovdb_min(pnanovdb_max(ijkf.z, pos.z), ijkf.z + dimf);
                    PNANOVDB_DEREF(vhit) = pnanovtt_cached_sample_cached_voxel_float(vtt_sampler, pos);
                    PNANOVDB_DEREF(poshit) = pos;
                    return PNANOVDB_TRUE;
                }
            }
        } // else if (pnanovtt_cached_sample_zero_crossing_float(vtt_sampler))

        pVal = pVal2;
        pTmin = tmin;
    } while (pnanovdb_hdda_step(PNANOVDB_REF(hdda)));

    return PNANOVDB_FALSE;
}


PNANOVDB_FORCE_INLINE pnanovdb_bool_t
first_active_texture_at_depth(pnanovdb_grid_type_t grid_type,
    pnanovdb_buf_t       buf,
    pnanovtt_texture_t   tex,
    PNANOVDB_INOUT(pnanovtt_readaccessor_t) acc,
    PNANOVDB_IN(pnanovdb_vec3_t) origin,
    float tmin,
    PNANOVDB_IN(pnanovdb_vec3_t) direction,
    PNANOVDB_INOUT(float) tmax,
    PNANOVDB_IN(pnanovdb_vec3_t) bbox_minf,
    PNANOVDB_IN(pnanovdb_vec3_t) bbox_maxf,
    PNANOVDB_INOUT(float) thit,
    pnanovdb_int32_t in_max_depth) {

    pnanovdb_bool_t hit = pnanovdb_hdda_ray_clip(
        bbox_minf, bbox_maxf, origin,
        PNANOVDB_REF(tmin), direction, tmax);
    if (!hit || PNANOVDB_DEREF(tmax) > 1.0e20f)
    {
        return PNANOVDB_FALSE;
    }

    // First we step delta along the ray to make sure we end up inside the voxel we're going to step through.
    // Then we can just use dim from that voxel. We also ensure that the first voxel is always inside the bbox.
    tmin += PNANOVTT_HDDA_DELTA;
    pnanovdb_vec3_t pos = pnanovdb_hdda_ray_start(origin, tmin, direction);
    pnanovdb_coord_t ijk = pnanovdb_hdda_pos_to_ijk(PNANOVDB_REF(pos));
    pnanovdb_int32_t dim = pnanovdb_uint32_as_int32(
        pnanovtt_readaccessor_get_dim_at_depth(grid_type, buf, acc, PNANOVDB_REF(ijk), in_max_depth));
    pnanovdb_hdda_t hdda;
    pnanovdb_hdda_init(PNANOVDB_REF(hdda), origin, tmin, direction, PNANOVDB_DEREF(tmax), dim);

    do {
        if (pnanovtt_readaccessor_is_texture_active_at_depth(grid_type, buf, tex, acc,
            PNANOVDB_REF(hdda.voxel), in_max_depth)) {
            // Here we step backwards by dim along the ray to ensure that we
            // step outside the active voxel such that the entire active voxel
            // is ray-marched.
            PNANOVDB_DEREF(thit) = hdda.tmin;
            return PNANOVDB_TRUE;
        }
        // Here we must update the hdda if dim changes
        // Get dim without tree access
        dim = pnanovdb_uint32_as_int32(pnanovtt_readaccessor_get_dim2(acc));
        pnanovdb_hdda_update(PNANOVDB_REF(hdda), origin, direction, dim);
    } while (pnanovdb_hdda_step(PNANOVDB_REF(hdda)));
    return PNANOVDB_FALSE;
}

PNANOVDB_FORCE_INLINE pnanovdb_bool_t
first_active_at_depth(pnanovdb_grid_type_t grid_type,
    pnanovdb_buf_t       buf,
    PNANOVDB_INOUT(pnanovtt_readaccessor_t) acc,
    PNANOVDB_IN(pnanovdb_vec3_t) origin,
    float tmin,
    PNANOVDB_IN(pnanovdb_vec3_t) direction,
    PNANOVDB_INOUT(float) tmax,
    PNANOVDB_IN(pnanovdb_vec3_t) bbox_minf,
    PNANOVDB_IN(pnanovdb_vec3_t) bbox_maxf,
    PNANOVDB_INOUT(float) thit,
    pnanovdb_int32_t in_max_depth) {

    pnanovdb_bool_t hit = pnanovdb_hdda_ray_clip(
        bbox_minf, bbox_maxf, origin,
        PNANOVDB_REF(tmin), direction, tmax);
    if (!hit || PNANOVDB_DEREF(tmax) > 1.0e20f)
    {
        return PNANOVDB_FALSE;
    }

    // First we step delta along the ray to make sure we end up inside the voxel we're going to step through.
    // Then we can just use dim from that voxel. We also ensure that the first voxel is always inside the bbox.
    tmin += PNANOVTT_HDDA_DELTA;
    pnanovdb_vec3_t pos = pnanovdb_hdda_ray_start(origin, tmin, direction);
    pnanovdb_coord_t ijk = pnanovdb_hdda_pos_to_ijk(PNANOVDB_REF(pos));
    pnanovdb_int32_t dim = pnanovdb_uint32_as_int32(
        pnanovtt_readaccessor_get_dim_at_depth(grid_type, buf, acc, PNANOVDB_REF(ijk), in_max_depth));
    pnanovdb_hdda_t hdda;
    pnanovdb_hdda_init(PNANOVDB_REF(hdda), origin, tmin, direction, PNANOVDB_DEREF(tmax), dim);

    do {
        if (pnanovtt_readaccessor_is_active_at_depth(grid_type, buf, acc,
            PNANOVDB_REF(hdda.voxel), in_max_depth)) {
            // Here we step backwards by dim along the ray to ensure that we
            // step outside the active voxel such that the entire active voxel
            // is ray-marched.
            PNANOVDB_DEREF(thit) = hdda.tmin;
            return PNANOVDB_TRUE;
        }
        // Here we must update the hdda if dim changes
        // Get dim without tree access
        dim = pnanovdb_uint32_as_int32(pnanovtt_readaccessor_get_dim2(acc));
        pnanovdb_hdda_update(PNANOVDB_REF(hdda), origin, direction, dim);
    } while (pnanovdb_hdda_step(PNANOVDB_REF(hdda)));
    return PNANOVDB_FALSE;
}


PNANOVDB_FORCE_INLINE pnanovdb_bool_t
first_active(pnanovdb_grid_type_t grid_type,
    pnanovdb_buf_t       buf,
    PNANOVDB_INOUT(pnanovtt_readaccessor_t) acc,
    PNANOVDB_IN(pnanovdb_vec3_t) origin,
    float tmin,
    PNANOVDB_IN(pnanovdb_vec3_t) direction,
    PNANOVDB_INOUT(float) tmax,
    PNANOVDB_IN(pnanovdb_vec3_t) bbox_minf,
    PNANOVDB_IN(pnanovdb_vec3_t) bbox_maxf,
    PNANOVDB_INOUT(float) thit)
{
    return first_active_at_depth(grid_type,
        buf,
        acc,
        origin,
        tmin,
        direction,
        tmax,
        bbox_minf,
        bbox_maxf,
        thit,
        PNANOVTT_TRAITS_MAX_DEPTH);
}

#endif // #ifdef PNANOVDB_HDDA

// ------------------------------------------------ Extensions to e.g. compile OpenCL code with cpp -----------------------------------------------------------

#ifdef PNANOVDB_CPP

typedef struct pnanovdb_bvec3_t
{
    pnanovdb_bvec3_t(const pnanovdb_bool_t in_x = false, const pnanovdb_bool_t in_y = false, const pnanovdb_bool_t in_z = false)
        : x(in_x)
        , y(in_y)
        , z(in_z)
    {}

    pnanovdb_bool_t x, y, z;
} pnanovdb_bvec3_t;

typedef struct pnanovdb_vec2_t
{
    float x, y;
} pnanovdb_vec2_t;

typedef struct pnanovdb_vec4_t
{
    pnanovdb_vec3_t xyz;
    float w;
} pnanovdb_vec4_t;


PNANOVDB_FORCE_INLINE pnanovdb_vec3_t pnanovdb_vec3_init(const float in_x, const float in_y, const float in_z)
{
    pnanovdb_vec3_t v;
    v.x = in_x;
    v.y = in_y;
    v.z = in_z;
    return v;
}

PNANOVDB_FORCE_INLINE pnanovdb_bvec3_t operator!=(const pnanovdb_vec3_t& in1, const pnanovdb_vec3_t& in2)
{
    return pnanovdb_bvec3_t((in1.x < in2.x || in1.x > in2.x), (in1.y < in2.y || in1.y > in2.y), (in1.z < in2.z || in1.z > in2.z));
}

PNANOVDB_FORCE_INLINE pnanovdb_vec3_t operator-(const pnanovdb_vec3_t& in1, const pnanovdb_vec3_t& in2)
{
    pnanovdb_vec3_t v;
    v.x = in1.x - in2.x;
    v.y = in1.y - in2.y;
    v.z = in1.z - in2.z;
    return v;
}

PNANOVDB_FORCE_INLINE pnanovdb_vec3_t operator-(const pnanovdb_vec3_t& in1)
{
    pnanovdb_vec3_t v;
    v.x = -in1.x;
    v.y = -in1.y;
    v.z = -in1.z;
    return v;
}

PNANOVDB_FORCE_INLINE pnanovdb_vec3_t operator+(const pnanovdb_vec3_t& in1, const pnanovdb_vec3_t& in2)
{
    pnanovdb_vec3_t v;
    v.x = in1.x + in2.x;
    v.y = in1.y + in2.y;
    v.z = in1.z + in2.z;
    return v;
}

PNANOVDB_FORCE_INLINE pnanovdb_vec3_t& operator+=(pnanovdb_vec3_t& in1, const pnanovdb_vec3_t& in2)
{
    in1.x += in2.x;
    in1.y += in2.y;
    in1.z += in2.z;
    return in1;
}

PNANOVDB_FORCE_INLINE pnanovdb_vec3_t& operator-=(pnanovdb_vec3_t& in1, const pnanovdb_vec3_t& in2)
{
    in1.x -= in2.x;
    in1.y -= in2.y;
    in1.z -= in2.z;
    return in1;
}

PNANOVDB_FORCE_INLINE pnanovdb_vec3_t operator*(const pnanovdb_vec3_t& in1, const pnanovdb_vec3_t& in2)
{
    pnanovdb_vec3_t v;
    v.x = in1.x * in2.x;
    v.y = in1.y * in2.y;
    v.z = in1.z * in2.z;
    return v;
}

PNANOVDB_FORCE_INLINE pnanovdb_vec3_t operator*(const pnanovdb_vec3_t& in1, const float in2)
{
    pnanovdb_vec3_t v;
    v.x = in1.x * in2;
    v.y = in1.y * in2;
    v.z = in1.z * in2;
    return v;
}

PNANOVDB_FORCE_INLINE pnanovdb_vec3_t operator*(const float in1, const pnanovdb_vec3_t& in2)
{
    pnanovdb_vec3_t v;
    v.x = in2.x * in1;
    v.y = in2.y * in1;
    v.z = in2.z * in1;
    return v;
}

PNANOVDB_FORCE_INLINE pnanovdb_vec3_t operator/(const pnanovdb_vec3_t& in1, const pnanovdb_vec3_t& in2)
{
    pnanovdb_vec3_t v;
    if (in2.x <= 0.f && in2.x >= 0.f)
        v.x = std::numeric_limits<float>::infinity();
    else
        v.x = in1.x / in2.x;
    if (in2.y <= 0.f && in2.y >= 0.f)
        v.y = std::numeric_limits<float>::infinity();
    else
        v.y = in1.y / in2.y;
    if (in2.z <= 0.f && in2.z >= 0.f)
        v.z = std::numeric_limits<float>::infinity();
    else
        v.z = in1.z / in2.z;
    return v;
}

PNANOVDB_FORCE_INLINE pnanovdb_vec3_t min(const pnanovdb_vec3_t& in1, const pnanovdb_vec3_t& in2)
{
    pnanovdb_vec3_t v;
    v.x = std::min(in1.x, in2.x);
    v.y = std::min(in1.y, in2.y);
    v.z = std::min(in1.z, in2.z);
    return v;
}

PNANOVDB_FORCE_INLINE pnanovdb_vec3_t max(const pnanovdb_vec3_t& in1, const pnanovdb_vec3_t& in2)
{
    pnanovdb_vec3_t v;
    v.x = std::max(in1.x, in2.x);
    v.y = std::max(in1.y, in2.y);
    v.z = std::max(in1.z, in2.z);
    return v;
}

PNANOVDB_FORCE_INLINE float dot(const pnanovdb_vec3_t& in1, const pnanovdb_vec3_t& in2)
{
    return in1.x * in2.x + in1.y * in2.y + in1.z * in2.z;
}

PNANOVDB_FORCE_INLINE pnanovdb_bool_t any(const pnanovdb_bvec3_t& in)
{
    return in.x || in.y || in.z;
}

PNANOVDB_FORCE_INLINE float length(const pnanovdb_vec3_t& in)
{
    return sqrtf(in.x * in.x + in.y * in.y + in.z * in.z);
}

PNANOVDB_FORCE_INLINE pnanovdb_vec3_t normalize(const pnanovdb_vec3_t& in)
{
    const float len = length(in);
    pnanovdb_vec3_t v;
    v.x = in.x / len;
    v.y = in.y / len;
    v.z = in.z / len;
    return v;
}

PNANOVDB_FORCE_INLINE float clamp(const float in, const float in_min, const float in_max)
{
    return std::max(std::min(in, in_max), in_min);
}

PNANOVDB_FORCE_INLINE float fract(const float in, float* out_floor)
{
    *out_floor = floor(in);
    return std::min(in - *out_floor, 1.f - std::numeric_limits<float>::epsilon());
}

PNANOVDB_FORCE_INLINE float saturate(const float in)
{
    return clamp(in, 0.f, 1.f);
}

PNANOVDB_FORCE_INLINE float smoothstep(const float in_x, const float in_edge0, const float in_edge1)
{
    float t;
    t = clamp((in_x - in_edge0) / (in_edge1 - in_edge0), 0.0f, 1.0f);
    return t * t * (3.0f - 2.0f * t);
}

PNANOVDB_FORCE_INLINE float RayBoxIntersectionMax(pnanovdb_vec3_t f3BBMin,
    pnanovdb_vec3_t f3BBMax,
    pnanovdb_vec3_t f3RayOrigin,
    pnanovdb_vec3_t f3RayInvDir) {
    pnanovdb_vec3_t f3OriginMin = (f3BBMin - f3RayOrigin) * f3RayInvDir;
    pnanovdb_vec3_t f3OriginMax = (f3BBMax - f3RayOrigin) * f3RayInvDir;
    pnanovdb_vec3_t f3Max = max(f3OriginMax, f3OriginMin);
    float           fEnd = std::min(f3Max.x, std::min(f3Max.y, f3Max.z));
    return fEnd;
}

#endif // PNANOVDB_CPP


#ifdef PNANOVDB_CL

typedef float2 pnanovdb_vec2_t;
typedef float4 pnanovdb_vec4_t;
typedef int3 pnanovdb_bvec3_t;

PNANOVDB_FORCE_INLINE float saturate(const float in)
{
    return clamp(in, 0.f, 1.f);
}

PNANOVDB_FORCE_INLINE float RayBoxIntersectionMax(pnanovdb_vec3_t f3BBMin,
    pnanovdb_vec3_t f3BBMax,
    pnanovdb_vec3_t f3RayOrigin,
    pnanovdb_vec3_t f3RayInvDir) {
    pnanovdb_vec3_t f3OriginMin = (f3BBMin - f3RayOrigin) * f3RayInvDir;
    pnanovdb_vec3_t f3OriginMax = (f3BBMax - f3RayOrigin) * f3RayInvDir;
    pnanovdb_vec3_t f3Max = max(f3OriginMax, f3OriginMin);
    float           fEnd = min(f3Max.x, min(f3Max.y, f3Max.z));
    return fEnd;
}

# endif // PNANOVDB_CL

#if defined(__clang__)
#pragma clang diagnostic pop
#endif

#endif // end of NANOVTT_PNANOVTT_H_HAS_BEEN_INCLUDED
