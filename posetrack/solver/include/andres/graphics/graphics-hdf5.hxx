#pragma once
#ifndef ANDRES_GRAPHICS_GRAPHICS_HDF5_HXX
#define ANDRES_GRAPHICS_GRAPHICS_HDF5_HXX

#include "hdf5.h"

#include "graphics.hxx"

namespace andres {
namespace graphics {
namespace hdf5 {

enum FileAccessMode {READ_ONLY, READ_WRITE};
enum HDF5Version {HDF5_VERSION_DEFAULT, HDF5_VERSION_LATEST};

template<class T> struct HDF5Type { hid_t type() const { throw std::runtime_error("conversion of this type to HDF5 type not implemented."); } };
template<> struct HDF5Type<char> { hid_t type() const { return H5T_NATIVE_CHAR; } };
template<> struct HDF5Type<unsigned char> { hid_t type() const { return H5T_NATIVE_UCHAR; } };
template<> struct HDF5Type<short> { hid_t type() const { return H5T_NATIVE_SHORT; } };
template<> struct HDF5Type<unsigned short> { hid_t type() const { return H5T_NATIVE_USHORT; } };
template<> struct HDF5Type<int> { hid_t type() const { return H5T_NATIVE_INT; } };
template<> struct HDF5Type<unsigned int> { hid_t type() const { return H5T_NATIVE_UINT; } };
template<> struct HDF5Type<long> { hid_t type() const { return H5T_NATIVE_LONG; } };
template<> struct HDF5Type<unsigned long> { hid_t type() const { return H5T_NATIVE_ULONG; } };
template<> struct HDF5Type<long long> { hid_t type() const { return H5T_NATIVE_LLONG; } };
template<> struct HDF5Type<unsigned long long> { hid_t type() const { return H5T_NATIVE_ULLONG; } };
template<> struct HDF5Type<float> { hid_t type() const { return H5T_NATIVE_FLOAT; } };
template<> struct HDF5Type<double> { hid_t type() const { return H5T_NATIVE_DOUBLE; } };
template<> struct HDF5Type<long double> { hid_t type() const { return H5T_NATIVE_LDOUBLE; } };

/// Create an HDF5 file.
///
/// \param filename Name of the file.
/// \param hdf5version HDF5 version tag.
///
/// \returns HDF5 handle
///
/// \sa openFile(), closeFile()
///
inline hid_t
createFile(
    const std::string& filename,
    HDF5Version hdf5version = HDF5_VERSION_DEFAULT
) {
    hid_t version = H5P_DEFAULT;
    if(hdf5version == HDF5_VERSION_LATEST) {
        version = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_libver_bounds(version, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST);
    }

    hid_t fileHandle = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, version);
    if(fileHandle < 0) {
        throw std::runtime_error("Could not create HDF5 file: " + filename);
    }

    return fileHandle;
}

/// Open an HDF5 file.
///
/// \param filename Name of the file.
/// \param fileAccessMode File access mode.
/// \param hdf5version HDF5 version tag.
///
/// \returns HDF5 handle
///
/// \sa closeFile(), createFile()
///
inline hid_t
openFile(
    const std::string& filename,
    FileAccessMode fileAccessMode = READ_ONLY,
    HDF5Version hdf5version = HDF5_VERSION_DEFAULT
) {
    hid_t access = H5F_ACC_RDONLY;
    if(fileAccessMode == READ_WRITE) {
        access = H5F_ACC_RDWR;
    }

    hid_t version = H5P_DEFAULT;
    if(hdf5version == HDF5_VERSION_LATEST) {
        version = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_libver_bounds(version, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST);
    }

    hid_t fileHandle = H5Fopen(filename.c_str(), access, version);
    if(fileHandle < 0) {
        throw std::runtime_error("Could not open HDF5 file: " + filename);
    }

    return fileHandle;
}

/// Close an HDF5 file
///
/// \param handle Handle to the HDF5 file.
///
/// \sa openFile(), createFile()
///
inline void closeFile(
    const hid_t& handle
) {
    H5Fclose(handle);
}

/// Create an HDF5 group.
///
/// \param parentHandle HDF5 handle on the parent group or file.
/// \param groupName Name of the group.
/// \returns HDF5 handle on the created group
///
/// \sa openGroup(), closeGroup()
///
inline hid_t
createGroup(
    const hid_t& parentHandle,
    const std::string& groupName
) {
    hid_t groupHandle = H5Gcreate(parentHandle, groupName.c_str(),
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if(groupHandle < 0) {
        throw std::runtime_error("Could not create HDF5 group.");
    }
    return groupHandle;
}

/// Open an HDF5 group.
///
/// \param parentHandle HDF5 handle on the parent group or file.
/// \param groupName Name of the group.
/// \returns HDF5 handle on the opened group.
///
/// \sa createGroup(), closeGroup()
///
inline hid_t
openGroup(
    const hid_t& parentHandle,
    const std::string& groupName
) {
    hid_t groupHandle = H5Gopen(parentHandle, groupName.c_str(), H5P_DEFAULT);
    if(groupHandle < 0) {
        throw std::runtime_error("Could not open HDF5 group.");
    }
    return groupHandle;
}

/// Close an HDF5 group.
///
/// \param handle HDF5 handle on group to close.
///
/// \sa openGroup(), createGroup()
///
inline void
closeGroup(
    const hid_t& handle
) {
    H5Gclose(handle);
}

/// Save a vector to an HDF5 dataset.
///
template<class T>
inline void
save(
    const hid_t parentHandle,
    const std::string datasetName,
    const std::vector<T>& data
) {
    hsize_t shape[] = {data.size()};
    hid_t dataspace = H5Screate_simple(1, shape, NULL);
    if(dataspace < 0) {
        throw std::runtime_error("could not create HDF5 dataspace.");
    }
    HDF5Type<T> typeMemory;
    hid_t dataset = H5Dcreate(parentHandle, datasetName.c_str(), typeMemory.type(), dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if(dataset < 0) {
        H5Sclose(dataspace);
        throw std::runtime_error("could not create HDF5 dataset.");
    }
    hid_t status = H5Dwrite(dataset, typeMemory.type(), H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());
    if(status < 0) {
        H5Dclose(dataset);
        H5Sclose(dataspace);
        throw std::runtime_error("could not write to HDF5 dataset.");
    }
}

/// Load a vector from an HDF5 dataset.
///
template<class T>
inline void
load(
    const hid_t parentHandle,
    const std::string datasetName,
    std::vector<T>& data
) {
    // open dataset and get types
    hid_t dataset = H5Dopen(parentHandle, datasetName.c_str(), H5P_DEFAULT);
    if(dataset < 0) {
        throw std::runtime_error("could not open HDF5 dataset.");
    }
    hid_t typeFile = H5Dget_type(dataset);

    // get dimension and shape
    hid_t filespace = H5Dget_space(dataset);
    int dimension = H5Sget_simple_extent_ndims(filespace);
    if(dimension != 1) {
        throw std::runtime_error("HDF5 dataset is not one-dimensional.");
    }
    hsize_t size = 0;
    herr_t status = H5Sget_simple_extent_dims(filespace, &size, NULL);
    if(status < 0) {
        H5Dclose(dataset);
        H5Tclose(typeFile);
        H5Sclose(filespace);
        throw std::runtime_error("could not get shape of HDF5 dataset.");
    }

    // read
    data.resize(size);
    status = H5Dread(dataset, HDF5Type<T>().type(), H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());

    // close dataset and types
    H5Dclose(dataset);
    H5Tclose(typeFile);
    H5Sclose(filespace);
    if(status < 0) {
        throw std::runtime_error("could not read from HDF5 dataset 'points'.");
    }
}

template<class T, class S>
class HDF5Type<andres::graphics::PointProperty<T, S> > {
private:
    typedef andres::graphics::PointProperty<T, S> PointPropertyType;

public:
    HDF5Type()
        : type_(H5Tcreate(H5T_COMPOUND, sizeof(PointPropertyType)))
        { H5Tinsert(type_, "visibility", HOFFSET(PointPropertyType, visibility_), HDF5Type<andres::graphics::Bit>().type());
          H5Tinsert(type_, "color", HOFFSET(PointPropertyType, color_), HDF5Type<andres::graphics::Color>().type());
          H5Tinsert(type_, "alpha", HOFFSET(PointPropertyType, alpha_), HDF5Type<andres::graphics::Color>().type()); }
    ~HDF5Type()
        { H5Tclose(type_); }
    hid_t type() const
        { return type_; }

private:
    hid_t type_;
};

template<class T, class S>
class HDF5Type<andres::graphics::Point<T, S> > {
private:
    typedef T value_type;
    typedef S size_type;
    typedef andres::graphics::Point<T, S> PointType;

public:
    HDF5Type()
        : type_(H5Tcreate(H5T_COMPOUND, sizeof(PointType)))
        { H5Tinsert(type_, "vector", HOFFSET(PointType, r_), HDF5Type<value_type>().type());
          H5Tinsert(type_, "property-index", HOFFSET(PointType, propertyIndex_), HDF5Type<size_type>().type()); }
    ~HDF5Type()
        { H5Tclose(type_); }
    hid_t type() const
        { return type_; }

private:
    hid_t type_;
};

template<class T, class S>
class HDF5Type<andres::graphics::LineProperty<T, S> > {
private:
    typedef andres::graphics::LineProperty<T, S> LinePropertyType;

public:
    HDF5Type()
        : type_(H5Tcreate(H5T_COMPOUND, sizeof(LinePropertyType)))
        { H5Tinsert(type_, "visibility", HOFFSET(LinePropertyType, visibility_), HDF5Type<andres::graphics::Bit>().type());
          H5Tinsert(type_, "color", HOFFSET(LinePropertyType, color_), HDF5Type<andres::graphics::Color>().type());
          H5Tinsert(type_, "alpha", HOFFSET(LinePropertyType, alpha_), HDF5Type<andres::graphics::Color>().type()); }
    ~HDF5Type()
        { H5Tclose(type_); }
    hid_t type() const
        { return type_; }

private:
    hid_t type_;
};


template<class T, class S>
class HDF5Type<andres::graphics::Line<T, S> > {
private:
    typedef T value_type;
    typedef S size_type;
    typedef andres::graphics::Line<T, S> LineType;

public:
    HDF5Type()
        : type_(H5Tcreate(H5T_COMPOUND, sizeof(LineType)))
        { H5Tinsert(type_, "point-indices", HOFFSET(LineType, pointIndices_), HDF5Type<size_type>().type());
          H5Tinsert(type_, "property-index", HOFFSET(LineType, propertyIndex_), HDF5Type<size_type>().type()); }
    ~HDF5Type()
        { H5Tclose(type_); }
    hid_t type() const
        { return type_; }

private:
    hid_t type_;
};

template<class T, class S>
void
save(
    const hid_t parentHandle,
    const andres::graphics::Graphics<T, S>& graphics
) {
    save(parentHandle, "point-properties", graphics.pointProperties());
    save(parentHandle, "points", graphics.points());
    save(parentHandle, "line-properties", graphics.lineProperties());
    save(parentHandle, "lines", graphics.lines());
}

template<class T, class S>
void
load(
    const hid_t parentHandle,
    andres::graphics::Graphics<T, S>& graphics
) {
    typedef andres::graphics::Graphics<T, S> GraphicsType;
    typedef typename GraphicsType::PointType PointType;
    typedef typename GraphicsType::PointsVector PointsVector;
    typedef typename GraphicsType::PointPropertiesVector PointPropertiesVector;
    typedef typename GraphicsType::LineType LineType;
    typedef typename GraphicsType::LinesVector LinesVector;
    typedef typename GraphicsType::LinePropertiesVector LinePropertiesVector;

    PointPropertiesVector pointProperties;
    PointsVector points;
    LinePropertiesVector lineProperties;
    LinesVector lines;

    load(parentHandle, "point-properties", pointProperties);
    load(parentHandle, "points", points);
    load(parentHandle, "line-properties", lineProperties);
    load(parentHandle, "lines", lines);

    graphics.assign(pointProperties, points, lineProperties, lines);
}

} // namespace hdf5
} // graphics
} // namespace andres

#endif // #ifndef ANDRES_GRAPHICS_GRAPHICS_HDF5_HXX
