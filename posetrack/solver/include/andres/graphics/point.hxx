#pragma once
#ifndef ANDRES_GRAPHICS_POINT_HXX
#define ANDRES_GRAPHICS_POINT_HXX

#include <cstddef>

#include "types.hxx"

namespace andres {
namespace graphics {

namespace hdf5 {
    template<class T> class HDF5Type;
}

template<class T = float, class S = std::size_t>
class Point {
public:
    typedef T value_type;
    typedef S size_type;

    Point()
        : Point(0, 0, 0)
        {}
    Point(const value_type x, const value_type y, const value_type z, const size_type propertyIndex = 0)
        : propertyIndex_(propertyIndex)
        { r_[0] = x; r_[1] = y; r_[2] = z; }
    value_type& operator[](const size_type index)
        { assert(index < 3); return r_[index]; }
    size_type& propertyIndex()
        { return propertyIndex_; }
    value_type operator[](const size_type index) const
        { assert(index < 3); return r_[index]; }
    size_type propertyIndex() const
        { return propertyIndex_; }
    bool operator==(const Point<T, S>& other) const
        { return r_[0] == other.r_[0]
            && r_[1] == other.r_[1]
            && r_[2] == other.r_[2]
            && propertyIndex_ == other.propertyIndex_;
        }

private:
    value_type r_[3];
    size_type propertyIndex_;

friend class andres::graphics::hdf5::HDF5Type<Point<T, S> >;
};

template<class T = float, class S = std::size_t>
class PointProperty {
public:
    typedef T value_type;
    typedef S size_type;

    PointProperty(const Bit visibility = true)
        : PointProperty(visibility, 0, 0, 0, 255)
        {}
    PointProperty(const Bit visibility, const Color r, const Color g, const Color b, const Color alpha)
        :   visibility_(visibility)
        { color_[0] = r; color_[1] = g; color_[2] = b; alpha_ = alpha; }
    Bit& visibility()
        { return visibility_; }
    Color& color(const size_type j)
        { assert(j < 3); return color_[j]; }
    Color& alpha()
        { return alpha_; }

    Bit visibility() const
        { return visibility_; }
    Color color(const size_type j) const
        { assert(j < 3); return color_[j]; }
    Color alpha() const
        { return alpha_; }
    bool operator==(const PointProperty<T, S>& other) const
        { return visibility_ == other.visibility_
            && color_[0] == other.color_[0]
            && color_[1] == other.color_[1]
            && color_[2] == other.color_[2]
            && alpha_ == other.alpha_;
        }

private:
    Bit visibility_;
    Color color_[3];
    Color alpha_;

friend class andres::graphics::hdf5::HDF5Type<PointProperty<T, S> >;
};

} // namespace graphics
} // namespace andres

#endif // #ifndef ANDRES_GRAPHICS_POINT_HXX
