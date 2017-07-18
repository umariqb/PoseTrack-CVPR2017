#pragma once
#ifndef POSE_TRACK_DETECTION_HXX
#define POSE_TRACK_DETECTION_HXX

#include <cstddef>
namespace pt {

template<class S = std::size_t>
struct Detection {
    typedef S size_type;

    Detection()
        :  clusterIndex_(),
           status_()
        {}
    size_type clusterIndex_;
    size_type status_;
};


} // namespace pt

#endif
