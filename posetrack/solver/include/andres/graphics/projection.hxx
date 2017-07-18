#pragma once
#ifndef ANDRES_GRAPHICS_PROJECTION_HXX
#define ANDRES_GRAPHICS_PROJECTION_HXX

#include <cstddef>

namespace andres {
namespace graphics {

template<class T = float, class S = std::size_t>
struct OrthogonalProjection {
    typedef T value_type;
    typedef S size_type;

    OrthogonalProjection()
        : OrthogonalProjection(1.0, 0.0, 0.5, 0.0, 1.0, 0.5)
        {}
    OrthogonalProjection(
        const value_type p00, const value_type p01, const value_type p02,
        const value_type p10, const value_type p11, const value_type p12
    )
        {
            parameters_[0][0] = p00; parameters_[0][1] = p01; parameters_[0][2] = p02;
            parameters_[1][0] = p10; parameters_[1][1] = p11; parameters_[1][2] = p12;
        }
    void operator()(const value_type x, const value_type y, const value_type z, value_type& r, value_type& s) const
        {
            r = parameters_[0][0] * x + parameters_[0][1] * y + parameters_[0][2] * z;
            s = parameters_[1][0] * x + parameters_[1][1] * y + parameters_[1][2] * z;
        }
    value_type parameters_[3][3];
};

} // namespace graphics
} // namespace andres

#endif // #ifndef ANDRES_GRAPHICS_PROJECTION_HXX
