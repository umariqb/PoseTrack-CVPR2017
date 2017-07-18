#pragma once
#ifndef ANDRES_GRAPHICS_SVG_HXX
#define ANDRES_GRAPHICS_SVG_HXX

#include <iostream>

#include "graphics.hxx"

namespace andres {
namespace graphics {

template<class Projection>
void
saveSVG(
    const Graphics<typename Projection::value_type, typename Projection::size_type>& g,
    const Projection& projection,
    std::ostream& output = std::cout
) {
    typedef typename Projection::value_type value_type;
    typedef typename Projection::size_type size_type;
    typedef Graphics<value_type, size_type> GraphicsType;
    typedef typename GraphicsType::PointType PointType;
    typedef typename GraphicsType::LineType LineType;
    typedef typename GraphicsType::PointPropertyType PointPropertyType;
    typedef typename GraphicsType::LinePropertyType LinePropertyType;

    // print header
    output << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>"
        << std::endl
        << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" "
        << "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">"
        << std::endl
        << "<svg version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\" "
        << "xmlns:xlink=\"http://www.w3.org/1999/xlink\""
//        << " width=" << ???
//        << " height=" << ???
        << ">"
        << std::endl;

    // plot points
    for(size_type j = 0; j < g.numberOfPoints(); ++j) {
        const PointType& point = g.point(j);
        const PointPropertyType& pointProperty = g.pointProperty(point.propertyIndex());
        if(pointProperty.visibility()) {
            value_type x = 0;
            value_type y = 0;
            projection(point[0], point[1], point[2], x, y);
            output << "<circle cx=\"" << x
                << "\" cy=\"" << y
                << "\" r=\"1pt\" style=\""
                << "stroke:rgb(" << static_cast<unsigned int>(pointProperty.color(0))
                << ", " << static_cast<unsigned int>(pointProperty.color(1))
                << ", " << static_cast<unsigned int>(pointProperty.color(2)) << ");"
                << " fill:rgb(" << static_cast<unsigned int>(pointProperty.color(0))
                << ", " << static_cast<unsigned int>(pointProperty.color(1))
                << ", " << static_cast<unsigned int>(pointProperty.color(2)) << ");"
                << " stroke-width:1pt;"
                << "\"/>"
                << std::endl;
        }
    }

    // plot lines
    for(size_type j = 0; j < g.numberOfLines(); ++j) {
        const LineType& line = g.line(j);
        const LinePropertyType& lineProperty = g.lineProperty(line.propertyIndex());
        if(lineProperty.visibility()) {
            const PointType& point0 = g.point(line.pointIndex(0));
            const PointType& point1 = g.point(line.pointIndex(1));
            value_type x1 = 0;
            value_type y1 = 0;
            projection(point0[0], point0[1], point0[2], x1, y1);
            value_type x2 = 0;
            value_type y2 = 0;
            projection(point1[0], point1[1], point1[2], x2, y2);
            output << "<line x1=\"" << x1 << "\""
                << " y1=\"" << y1 << "\""
                << " x2=\"" << x2 << "\""
                << " y2=\"" << y2 << "\""
                << " style=\"stroke:rgb(" << static_cast<unsigned int>(lineProperty.color(0))
                << ", " << static_cast<unsigned int>(lineProperty.color(1))
                << ", " << static_cast<unsigned int>(lineProperty.color(2))
                << "); stroke-width:1pt;" << "\"/>"
                << std::endl;
        }
    }

    // print footer
    output << "</svg>\n";
}

} // namespace graphics
} // namespace andres

#endif // #ifndef ANDRES_GRAPHICS_SVG_HXX
