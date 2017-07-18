#pragma once
#ifndef POSE_TRACK_PROBLEM_HXX
#define POSE_TRACK_PROBLEM_HXX

#include <cassert>
#include <cstddef>

#include <andres/marray.hxx>

#include "join.hxx"

namespace pt {

template<class T = double, class S = std::size_t>
class Problem {
public:
    typedef T value_type;
    typedef S size_type;
    typedef JoinData<size_type> JoinDataType;
    typedef JoinIndex<size_type> JoinIndexType;
    typedef JoinMap<size_type> JoinMapType;

    Problem(const size_type = 0);
    void assign(const size_type);
    void setPartClassProbability(const size_type, const size_type, const value_type);
    void setJoinProbability(const size_type, const size_type, const size_type, const size_type, const value_type, const bool);
    JoinMapType& joinMapSpatial();
    JoinMapType& joinMapTemporal();

    size_type numberOfDetections() const;
    value_type getPartClassProbability(const size_type) const;
    const JoinMapType& joinMapSpatial() const;
    const JoinMapType& joinMapTemporal() const;


private:

    size_type numberOfPartClasses_;
    size_type numberOfDetections_;
    std::vector<value_type> partClassProbabilities_;
    JoinMapType joinMapSpatial_;
    JoinMapType joinMapTemporal_;
};

template<class T, class S>
Problem<T, S>::Problem(
	const size_type numberOfDetections
)
:   numberOfDetections_(numberOfDetections)
{
}

template<class T, class S>
void
Problem<T, S>::assign(
	const size_type numberOfDetections
) {
	numberOfDetections_ = numberOfDetections;
	partClassProbabilities_.clear();
	partClassProbabilities_.resize(numberOfDetections, 0);
    joinMapSpatial_.clear();
    joinMapTemporal_.clear();
}


template<class T, class S>
inline void
Problem<T, S>::setPartClassProbability(
    const size_type detection,
    const size_type partClass,
    const value_type probability
) {
    assert(detection < numberOfDetections());
    assert(probability >= 0 && probability <= 1);

    partClassProbabilities_[detection] = probability;
}


template<class T, class S>
inline void
Problem<T, S>::setJoinProbability(
    const size_type detection0,
    const size_type detection1,
    const size_type partClass0,
    const size_type partClass1,
    const value_type probability,
    const bool isSpatialJoin
) {
    JoinIndexType joinIndex(detection0, detection1, partClass0, partClass1);
    JoinDataType joinData(probability);

    if(isSpatialJoin)
    	joinMapSpatial_[joinIndex] = joinData;
    else
    	joinMapTemporal_[joinIndex] = joinData;
}

template<class T, class S>
typename Problem<T, S>::JoinMapType&
Problem<T, S>::joinMapSpatial() {
    return joinMapSpatial_;
}

template<class T, class S>
typename Problem<T, S>::JoinMapType&
Problem<T, S>::joinMapTemporal() {
    return joinMapTemporal_;
}

template<class T, class S>
inline typename Problem<T, S>::size_type
Problem<T, S>::numberOfDetections() const {
    return numberOfDetections_;
}

template<class T, class S>
inline typename Problem<T, S>::value_type
Problem<T, S>::getPartClassProbability(
    const size_type detection
) const {
	assert(detection < partClassProbabilities_.size());
    return partClassProbabilities_[detection];
}

template<class T, class S>
const typename Problem<T, S>::JoinMapType&
Problem<T, S>::joinMapSpatial() const {
    return joinMapSpatial_;
}

template<class T, class S>
const typename Problem<T, S>::JoinMapType&
Problem<T, S>::joinMapTemporal() const {
    return joinMapTemporal_;
}

} // namespace pt

#endif
