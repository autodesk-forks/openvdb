// Copyright Contributors to the OpenVDB Project
// SPDX-License-Identifier: MPL-2.0

/*!
    \file GridHandle.h

    \authors Ken Museth and Autodesk
*/

#ifndef NANOVTT_GRID_HANDLE_H_HAS_BEEN_INCLUDED
#define NANOVTT_GRID_HANDLE_H_HAS_BEEN_INCLUDED

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
#pragma clang diagnostic ignored "-Wweak-vtables"
#endif
#include <nanovdb/util/HostBuffer.h>
#include <nanovdb/util/GridHandle.h>
#if defined(__clang__)
#pragma clang diagnostic pop
#endif
#include <nanovdb/nanovtt/NanoVTT.h>

namespace nanovtt {

// --------------------------> GridHandle <------------------------------------

/// @brief This class serves to manage a raw memory buffer of a NanoVDB Grid.
///
/// @note  It is important to note that is class does NOT depend on OpenVDB.
template<typename BufferT = nanovdb::HostBuffer>
class GridHandle
{
protected:
    BufferT mBuffer;

    template<typename ValueT>
    const NanoVTTGrid<ValueT>* getGrid(uint32_t n = 0) const;

    const uint8_t* getGridPtr(uint32_t n = 0) const;

    template<typename ValueT, typename U = BufferT>
    typename std::enable_if<nanovdb::BufferTraits<U>::hasDeviceDual, const NanoVTTGrid<ValueT>*>::type
    getDeviceGrid(uint32_t n = 0) const;

    template<typename U = BufferT>
    typename std::enable_if<nanovdb::BufferTraits<U>::hasDeviceDual, const uint8_t*>::type
        getDeviceGridPtr(uint32_t n = 0) const;

    template <typename T>
    static T* no_const(const T* ptr) { return const_cast<T*>(ptr); }

public:

    /// @brief Move constructor from a buffer
    GridHandle(BufferT&& buffer) { mBuffer = std::move(buffer); }
    /// @brief Empty ctor
    GridHandle() = default;
    /// @brief Disallow copy-construction
    GridHandle(const GridHandle&) = delete;
    /// @brief Disallow copy assignment operation
    GridHandle& operator=(const GridHandle&) = delete;
    /// @brief Move copy assignment operation
    GridHandle& operator=(GridHandle&& other) noexcept
    {
        mBuffer = std::move(other.mBuffer);
        return *this;
    }
    /// @brief Move copy-constructor
    GridHandle(GridHandle&& other) noexcept { mBuffer = std::move(other.mBuffer); }
    /// @brief Default destructor
    ~GridHandle() { reset(); }
    /// @brief clear the buffer
    void reset() { mBuffer.clear(); }

    /// @brief Return a reference to the buffer
    BufferT&       buffer() { return mBuffer; }

    /// @brief Return a const reference to the buffer
    const BufferT& buffer() const { return mBuffer; }

    /// @brief Returns a non-const pointer to the data.
    ///
    /// @warning Note that the return pointer can be NULL if the GridHandle was not initialized
    uint8_t* data() { return mBuffer.data(); }

    /// @brief Returns a const pointer to the data.
    ///
    /// @warning Note that the return pointer can be NULL if the GridHandle was not initialized
    const uint8_t* data() const { return mBuffer.data(); }

    /// @brief Returns a non-const pointer to the data.
    ///
    /// @warning Note that the return pointer can be NULL if the GridHandle was not initialized
    uint8_t* deviceData() { return mBuffer.deviceData(); }

    /// @brief Returns a const pointer to the data.
    ///
    /// @warning Note that the return pointer can be NULL if the GridHandle was not initialized
    const uint8_t* deviceData() const { return mBuffer.deviceData(); }

    /// @brief Returns the size in bytes of the raw memory buffer managed by this GridHandle's allocator.
    uint64_t size() const { return mBuffer.size(); }

    /// @brief Returns a const pointer to the container encoded in this GridHandle.
    ///
    /// @warning Note that the return pointer can be NULL if the GridHandle was not initialized
    const Container* container() const { return reinterpret_cast<const Container*>(mBuffer.data()); }

    /// @brief Returns a const pointer to the device container encoded in this GridHandle.
    ///
    /// @warning Note that the return pointer can be NULL if the GridHandle was not initialized
    const Container* deviceContainer() const { return reinterpret_cast<const Container*>(deviceData()); }

    /// @brief Returns a const pointer to the @a n'th NanoVDB grid encoded in this GridHandle.
    ///
    /// @warning Note that the return pointer can be NULL if the GridHandle was not initialized, @a n is invalid
    ///
    ///          or if the template parameter does not match the specified grid!
    template<typename ValueT>
    const NanoVTTGrid<ValueT>* grid(uint32_t n = 0) const { return this->template getGrid<ValueT>(n); }

    /// @brief Returns a pointer to the @a n'th  NanoVDB grid encoded in this GridHandle.
    ///
    /// @warning Note that the return pointer can be NULL if the GridHandle was not initialized, @a n is invalid
    ///          or if the template parameter does not match the specified grid!
    template<typename ValueT>
    NanoVTTGrid<ValueT>* grid(uint32_t n = 0) { return no_const(this->template getGrid<ValueT>(n)); }

    /// @brief Returns a const pointer to the @a n'th NanoVDB grid encoded in this GridHandle.
    ///
    /// @warning Note that the return pointer can be NULL if the GridHandle was not initialized or @a n is invalid!
    const uint8_t* gridPtr(uint32_t n = 0) const { return this->getGridPtr(n); }

    /// @brief Returns a pointer to the @a n'th  NanoVDB grid encoded in this GridHandle.
    ///
    /// @warning Note that the return pointer can be NULL if the GridHandle was not initialized or @a n is invalid!
    uint8_t* gridPtr(uint32_t n = 0) { return no_const(this->getGridPtr(n)); }

    /// @brief Return a const pointer to the @a n'th grid encoded in this GridHandle on the device, e.g. GPU
    ///
    /// @warning Note that the return pointer can be NULL if the GridHandle was not initialized, @a n is invalid
    ///          or if the template parameter does not match the specified grid!
    template<typename ValueT, typename U = BufferT>
    typename std::enable_if<nanovdb::BufferTraits<U>::hasDeviceDual, const NanoVTTGrid<ValueT>*>::type
    deviceGrid(uint32_t n = 0) const { return this->template getDeviceGrid<ValueT>(n); }

    /// @brief Return a const pointer to the @a n'th grid encoded in this GridHandle on the device, e.g. GPU
    ///
    /// @warning Note that the return pointer can be NULL if the GridHandle was not initialized, @a n is invalid
    ///          or if the template parameter does not match the specified grid!
    template<typename ValueT, typename U = BufferT>
    typename std::enable_if<nanovdb::BufferTraits<U>::hasDeviceDual, NanoVTTGrid<ValueT>*>::type
    deviceGrid(uint32_t n = 0) { return no_const(this->template getDeviceGrid<ValueT>(n)); }

    /// @brief Return a const pointer to the @a n'th grid encoded in this GridHandle on the device, e.g. GPU
    ///
    /// @warning Note that the return pointer can be NULL if the GridHandle was not initialized or @a n is invalid!
    template<typename U = BufferT>
    typename std::enable_if<nanovdb::BufferTraits<U>::hasDeviceDual, const uint8_t*>::type
    deviceGridPtr(uint32_t n = 0) const { return this->getDeviceGridPtr(n); }

    /// @brief Return a const pointer to the @a n'th grid encoded in this GridHandle on the device, e.g. GPU
    ///
    /// @warning Note that the return pointer can be NULL if the GridHandle was not initialized or @a n is invalid!
    template<typename U = BufferT>
    typename std::enable_if<nanovdb::BufferTraits<U>::hasDeviceDual, uint8_t*>::type
    deviceGridPtr(uint32_t n = 0) { return no_const(this->getDeviceGridPtr(n)); }

    /// @brief Upload the grid to the device, e.g. from CPU to GPU
    ///
    /// @note This method is only available if the buffer supports devices
    template<typename U = BufferT>
    typename std::enable_if<nanovdb::BufferTraits<U>::hasDeviceDual, void>::type
    deviceUpload(void* stream = nullptr, bool sync = true) { mBuffer.deviceUpload(stream, sync); }

    /// @brief Download the grid to from the device, e.g. from GPU to CPU
    ///
    /// @note This method is only available if the buffer supports devices
    template<typename U = BufferT>
    typename std::enable_if<nanovdb::BufferTraits<U>::hasDeviceDual, void>::type
    deviceDownload(void* stream = nullptr, bool sync = true) { mBuffer.deviceDownload(stream, sync); }

    /// @brief Return true if this handle is empty, i.e. has no allocated memory
    bool empty() const { return size() == 0; }

    /// @brief Return true if this handle contains a grid
    operator bool() const { return !this->empty(); }

    /// @brief Returns a const pointer to the n'th grid meta data (see definition above).
    ///
    /// @warning Note that the return pointer can be NULL if the GridHandle was
    /// not initialized or less than n grids are stored in this GridHandle.
    const GridMetaData* gridMetaData(uint32_t n = 0) const
    {
        const uint8_t* data = mBuffer.data();
        if (data == nullptr)
            return nullptr;
        const Container* container = reinterpret_cast<const Container*>(data);
        data += container->memUsage();
        const GridMetaData* grid = reinterpret_cast<const GridMetaData*>(data);
        if (n >= grid->gridCount()) {// un-initialized or index is out of range
            return nullptr;
        }
        while (n != grid->gridIndex()) {
            data += grid->gridSize();
            grid = reinterpret_cast<const GridMetaData*>(data);
        }
        return grid;
    }

    /// @brief Returns the GridType of the n'th grid handled by this instance,
    /// and GridType::End if empty or less than n grids are stored in this
    /// GridHandle.
    nanovdb::GridType gridType(uint32_t n = 0) const
    {
        const uint8_t* data = mBuffer.data();
        if (data == nullptr)
            return nanovdb::GridType::End;
        const Container* container = reinterpret_cast<const Container*>(data);
        data += container->memUsage();
        const GridMetaData* grid = reinterpret_cast<const GridMetaData*>(data);
        if (n >= grid->gridCount()) {// un-initialized or index is out of range
            return nanovdb::GridType::End;
        }
        while (n != grid->gridIndex()) {
            data += grid->gridSize();
            grid = reinterpret_cast<const GridMetaData*>(data);
        }
        return grid->gridType();
    }

    /// @brief Return the number of grids contained in this buffer
    uint32_t gridCount() const
    {
        auto const* const ptr = this->gridMetaData();
        return ptr ? ptr->gridCount() : 0;
    }

}; // GridHandle

// --------------------------> Implementation of private methods in GridHandle <------------------------------------

template<typename BufferT>
template<typename ValueT>
inline const NanoVTTGrid<ValueT>* GridHandle<BufferT>::getGrid(uint32_t index) const
{
    using GridT = const NanoVTTGrid<ValueT>;
    const uint8_t* data = mBuffer.data();
    if (data == nullptr)
        return nullptr;
    const Container *container = reinterpret_cast<const Container*>(data);
    data += container->memUsage();
    const GridT *grid = reinterpret_cast<const GridT*>(data);
    if (index >= grid->gridCount()) {// un-initialized or index is out of range
        return nullptr;
    }
    while(index != grid->gridIndex()) {
        data += grid->gridSize();
        grid  = reinterpret_cast<const GridT*>(data);
    }
    return grid->gridType() == nanovdb::mapToGridType<ValueT>() ? grid : nullptr;
}

template<typename BufferT>
template<typename ValueT, typename U>
inline typename std::enable_if<nanovdb::BufferTraits<U>::hasDeviceDual, const NanoVTTGrid<ValueT>*>::type
GridHandle<BufferT>::getDeviceGrid(uint32_t index) const
{
    using GridT = const NanoVTTGrid<ValueT>;
    const uint8_t* data = mBuffer.data();
    const uint8_t* dev = mBuffer.deviceData();
    if (data == nullptr || dev == nullptr)
        return nullptr;
    const Container *container = reinterpret_cast<const Container*>(data);
    data += container->memUsage();
    GridT *grid = reinterpret_cast<GridT*>(data);
    if (index >= grid->gridCount()) {// un-initialized or index is out of range
        return nullptr;
    }
    dev += container->memUsage();
    while(index != grid->gridIndex()) {
        data += grid->gridSize();
        dev  += grid->gridSize();
        grid  = reinterpret_cast<GridT*>(data);
    }
    return grid->gridType() == nanovdb::mapToGridType<ValueT>() ? reinterpret_cast<GridT*>(dev) : nullptr;
}

template<typename BufferT>
inline const uint8_t* GridHandle<BufferT>::getGridPtr(uint32_t index) const
{
    const uint8_t* data = mBuffer.data();
    if (data == nullptr)
        return nullptr;
    const Container* container = reinterpret_cast<const Container*>(data);
    data += container->memUsage();
    const GridMetaData* grid = reinterpret_cast<const GridMetaData*>(data);
    if (index >= grid->gridCount()) {// un-initialized or index is out of range
        return nullptr;
    }
    while (index != grid->gridIndex()) {
        data += grid->gridSize();
        grid = reinterpret_cast<const GridMetaData*>(data);
    }
    return data;
}

template<typename BufferT>
template<typename U>
inline typename std::enable_if<nanovdb::BufferTraits<U>::hasDeviceDual, const uint8_t*>::type
GridHandle<BufferT>::getDeviceGridPtr(uint32_t index) const
{
    const uint8_t* data = mBuffer.data();
    const uint8_t* dev = mBuffer.deviceData();
    if (data == nullptr || dev == nullptr)
        return nullptr;
    const Container* container = reinterpret_cast<const Container*>(data);
    data += container->memUsage();
    dev += container->memUsage();
    const GridMetaData* grid = reinterpret_cast<const GridMetaData*>(data);
    if (index >= grid->gridCount()) {// un-initialized or index is out of range
        return nullptr;
    }
    while (index != grid->gridIndex()) {
        data += grid->gridSize();
        dev += grid->gridSize();
        grid = reinterpret_cast<const GridMetaData*>(data);
    }
    return dev;
}

} // namespace nanovtt

#endif // NANOVTT_GRID_HANDLE_H_HAS_BEEN_INCLUDED
