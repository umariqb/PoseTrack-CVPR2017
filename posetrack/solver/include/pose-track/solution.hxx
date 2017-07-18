#pragma once
#ifndef POSE_TRACK_SOLUTION_HXX
#define POSE_TRACK_SOLUTION_HXX

#include <cstddef>
#include <vector>

#include "detection.hxx"

namespace pt {

template<class S = std::size_t>
using Solution = std::vector<Detection<S> >;

} // namespace pt

#endif
