#pragma once
#ifndef POSE_TRACK_JOIN_HXX
#define POSE_TRACK_JOIN_HXX

#include <cassert>
#include <cstddef>
#include <map>

namespace pt {

template<class S = std::size_t>
struct JoinData {
    typedef S size_type;

    JoinData(const double probability = 0.0, const size_type variableIndex = 0)
        :   probability_(probability),
            variableIndex_(variableIndex)
        {
            assert(probability_ >= 0.0 && probability_ <= 1.0);
        }
    void assign(const double probability = 0.0, const size_type variableIndex = 0)
        {
            setProbability(probability);
            setVariableIndex(variableIndex);
        }
    void setProbability(const double probability = 0.0)
        {
            assert(probability >= 0.0 && probability <= 1.0);
            probability_ = probability;
        }
    void setVariableIndex(const size_type variableIndex = size_type())
        {
            variableIndex_ = variableIndex;
        }
    size_type getVariableIndex() const
        {
            return variableIndex_;
        }
    double getProbability() const
        {
            return probability_;
        }

private:
    double probability_;
    size_type variableIndex_;
};

template<class S = std::size_t>
struct JoinIndex {
    typedef S size_type;

    JoinIndex()
        {
            assign(0, 0, 0, 0);
        }
    JoinIndex(const size_type d0, const size_type d1, const size_type c0, const size_type c1)
        {
            assign(d0, d1, c0, c1);
        }
    void assign(const size_type d0, const size_type d1, const size_type c0, const size_type c1)
        {
            if(d0 <= d1) {
                detections_[0] = d0;
                detections_[1] = d1;
                partClasses_[0] = c0;
                partClasses_[1] = c1;
            }
            else {
                detections_[0] = d1;
                detections_[1] = d0;
                partClasses_[0] = c1;
                partClasses_[1] = c0;
            }
        }
    size_type getDetection(const size_type j) const
        {
            assert(j < 2);
            return detections_[j];
        }
    size_type getPartClass(const size_type j) const
        {
            assert(j < 2);
            return partClasses_[j];
        }
    bool operator<(const JoinIndex& other) const
        {
            // lexicographic order:
            if(detections_[0] < other.detections_[0]) {
                return true;
            }
            else if(detections_[0] > other.detections_[0]){
                return false;
            }
            else {
                if(detections_[1] < other.detections_[1]) {
                    return true;
                }
                else if(detections_[1] > other.detections_[1]){
                    return false;
                }
                else {
                    if(partClasses_[0] < other.partClasses_[0]) {
                        return true;
                    }
                    else if(partClasses_[0] > other.partClasses_[0]){
                        return false;
                    }
                    else {
                        if(partClasses_[1] < other.partClasses_[1]) {
                            return true;
                        }
                        else { // >=
                            return false;
                        }
                    }
                }
            }
        }

private:
    size_type detections_[2];
    size_type partClasses_[2];
};


template<class S = std::size_t>
using JoinMap = std::map<JoinIndex<S>, JoinData<S> >;

} // namespace pt

#endif
