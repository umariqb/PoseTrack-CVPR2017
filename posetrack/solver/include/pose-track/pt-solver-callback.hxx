#pragma once
#ifndef POSE_TRACK_SOLVER_CALLBACK_HXX
#define POSE_TRACK_SOLVER_CALLBACK_HXX

#include <cstddef>
#include <limits>
#include <iostream>
#include <vector>
#include <string>

#include <andres/graph/graph.hxx>
#include <andres/graph/components.hxx>
#include "andres/graph/multicut/kernighan-lin.hxx"
#include <andres/timer.hxx>

#include "problem.hxx"
#include "solution.hxx"
#include "andres/graphics/graphics-hdf5.hxx"


namespace pt {

template<class S = size_t>
struct FeasibleSolutionCallbackEmpty {
    typedef S size_type;
    typedef Solution<size_type> SolutionType;

    void operator()(const SolutionType& solution) const
        {}
};

template<
    class ILP_SOLVER,
    class FEASBILE_SOLUTION_CALLBACK = FeasibleSolutionCallbackEmpty<>
>
class PoseTracker {
public:
    typedef ILP_SOLVER IlpSolver;
    typedef FEASBILE_SOLUTION_CALLBACK FeasibleSolutionCallback;
    typedef size_t size_type;
    typedef Problem<double, size_type> ProblemType;
    typedef Solution<size_type> SolutionType;
    typedef andres::Timer<double> TimerType;

    PoseTracker(ProblemType&, SolutionType&, const SolutionType &, FeasibleSolutionCallback = FeasibleSolutionCallback(), const bool = false, const int = 36000, const double = 1.0 / 255.0);

private:
    typedef typename ProblemType::JoinIndexType JoinIndexType;
    typedef typename ProblemType::JoinDataType JoinDataType;
    typedef typename ProblemType::JoinMapType JoinMapType;
    typedef typename andres::graph::Graph<> GraphType;
    typedef andres::graph::ComponentsBySearch<GraphType> ComponentsType;

    struct SubgraphMask {
        typedef PoseTracker<IlpSolver, FeasibleSolutionCallback> PoseTrackerType;

        SubgraphMask(const PoseTrackerType& poseTracker)
            : poseTracker_(poseTracker) {}
        bool vertex(const size_t v) const
            { return true; }
        bool edge(const size_t e) const
            {
                size_t const v0 = poseTracker_.detectionGraph_.vertexOfEdge(e, 0);
                size_t const v1 = poseTracker_.detectionGraph_.vertexOfEdge(e, 1);
                return poseTracker_.ilpSolver_.label(poseTracker_.y(v0, v1)) > .5;
            }
        const PoseTrackerType& poseTracker_;
    };

    class Callback
    :   public IlpSolver::Callback
    {
    public:        
        typedef typename IlpSolver::Callback IlpSolverCallback;
        typedef PoseTracker<IlpSolver, FeasibleSolutionCallback> PoseTrackerType;
        typedef typename PoseTrackerType::size_type size_type;
        typedef typename PoseTrackerType::ComponentsType ComponentsType;

        struct SubgraphMask {
            typedef Callback CallbackType;

            SubgraphMask(CallbackType& callback)
                : callback_(callback) {}
            bool vertex(const size_t v) const
                { return true; }
            bool edge(const size_t e) const
                {
                    const size_type v0 = callback_.poseTracker_.detectionGraph_.vertexOfEdge(e, 0);
                    const size_type v1 = callback_.poseTracker_.detectionGraph_.vertexOfEdge(e, 1);
                    return callback_.label(callback_.poseTracker_.y(v0, v1)) > .5;
                }
            Callback& callback_;
        };

        Callback(PoseTrackerType&);
        void separateAndAddLazyConstraints();

    private:
        size_type separateAndAddViolated3CycleConstraints();

        PoseTrackerType& poseTracker_;
    };

    // ilp
    void setObjectiveFunction();
    size_type addAllCouplingConstraints();
    size_type addAll3CycleConstraints();
    size_type addAllSpatioTemporalConstraints();
    size_type addPartialSolutionConstraints(const SolutionType &partial);



    void plotSolution();

    // variable indexing
    size_type y(const size_type, const size_type) const;
    size_type getIndexOffset() const;

    size_type numberOfVariablesX_ { 0 };
    size_type numberOfVariablesY_ { 0 };
    size_type numberOfVariablesT_ { 0 };
    size_type numberOfFixedVariables_ { 0 };
    std::vector<size_type> temporalEdgeIndexes;
    std::vector<size_type> temporalEdgePartClass;
    ProblemType& problem_;
    double epsilonProbability_;
    GraphType detectionGraph_;
    TimerType timer_;
    FeasibleSolutionCallback feasibleSolutionCallback_;
    IlpSolver ilpSolver_;
};

template<class ILP_SOLVER, class FEASBILE_SOLUTION_CALLBACK>
PoseTracker<ILP_SOLVER, FEASBILE_SOLUTION_CALLBACK>::PoseTracker(
    ProblemType& problem,
    SolutionType& solution,
    const SolutionType& partialSolution,
    FeasibleSolutionCallback feasibleSolutionCallback,
    const bool withSingleClusterConstraints,
    const int timeLimit,
    const double epsilonProbability

)
:   
	temporalEdgeIndexes(),
    problem_(problem),
    epsilonProbability_(epsilonProbability),
    detectionGraph_(problem.numberOfDetections()),
    feasibleSolutionCallback_(feasibleSolutionCallback)
{
    ilpSolver_.setVerbosity(true);
    ilpSolver_.setRelativeGap(0.01);
    ilpSolver_.setMemoryLimit(32);
    ilpSolver_.setTimeLimit(timeLimit);

    timer_.start();

    setObjectiveFunction();
    if(partialSolution.size() != 0)
	addPartialSolutionConstraints(partialSolution);
    addAll3CycleConstraints();
	addAllCouplingConstraints();
	addAllSpatioTemporalConstraints();


    Callback callback(*this);
    ilpSolver_.setCallback(callback);

    ilpSolver_.optimize();
    timer_.stop();

    // save solution
	ComponentsType components;
	components.build(detectionGraph_, SubgraphMask(*this));
	solution.resize(problem_.numberOfDetections());
	for (size_t d = 0; d < problem.numberOfDetections(); ++d)
	{
		solution[d].clusterIndex_ = components.labels_[d];

		if(ilpSolver_.label(d) > .5)
			solution[d].status_ = 1;
	}
}

template<class ILP_SOLVER, class FEASBILE_SOLUTION_CALLBACK>
void
PoseTracker<ILP_SOLVER, FEASBILE_SOLUTION_CALLBACK>::plotSolution() {


}


template<class ILP_SOLVER, class FEASBILE_SOLUTION_CALLBACK>
void
PoseTracker<ILP_SOLVER, FEASBILE_SOLUTION_CALLBACK>::setObjectiveFunction() {
    std::cout << "setting up objective function:" << std::endl;

    // set coefficients of variables
    numberOfVariablesX_ = problem_.numberOfDetections();
    numberOfVariablesY_ = problem_.joinMapSpatial().size();
    numberOfVariablesT_ = problem_.joinMapTemporal().size();
    std::vector<double> coefficients;
    coefficients.reserve(numberOfVariablesX_ + numberOfVariablesY_ + numberOfVariablesT_);
    coefficients.resize(numberOfVariablesX_,0);
    // set coefficients of variables x
    for (size_t d = 0; d < problem_.numberOfDetections(); ++d){
        double const p = problem_.getPartClassProbability(d);

		if (p > epsilonProbability_)
			coefficients[d] = std::log( (1.0 - p) / p );
	}

    // set coefficients of variables y
	auto its=problem_.joinMapSpatial().begin();
    for (; its!=problem_.joinMapSpatial().end(); ++its)
    {
    	JoinIndexType joinIndex = its->first;
        size_type d0 = joinIndex.getDetection(0);
        size_type d1 = joinIndex.getDetection(1);

    	double p = its->second.getProbability();

		if(p == 0) { p = 1e-15; }
		if(p == 1) { p = 1-1e-15; }

		size_type edgeIndex = getIndexOffset() + detectionGraph_.insertEdge(d0, d1);
		assert(edgeIndex == coefficients.size());
		its->second.setVariableIndex(coefficients.size());
		value_type cost = std::log( (1.0 - p) / p );
		coefficients.push_back(cost);
    }

    // set coefficients of variables t
	auto itt=problem_.joinMapTemporal().begin();
    for (; itt!=problem_.joinMapTemporal().end(); ++itt)
    {
    	JoinIndexType joinIndex = itt->first;
        size_type d0 = joinIndex.getDetection(0);
        size_type d1 = joinIndex.getDetection(1);
        size_type partClass = joinIndex.getPartClass(0);

    	double p = itt->second.getProbability();

    	// within frame temporal edges are already added during addition spatial edges
    	// no need to add them again.
    	if(detectionGraph_.findEdge(d0, d1).first){
    		continue;
    	}

		size_type edgeIndex = getIndexOffset() + detectionGraph_.insertEdge(d0, d1);
		assert(edgeIndex == coefficients.size());
		itt->second.setVariableIndex(coefficients.size());

		if(p == 0) { p = 1e-15; }
		if(p == 1) { p = 1-1e-15; }
		value_type cost = std::log( (1.0 - p) / p );
		coefficients.push_back(cost);
		temporalEdgeIndexes.push_back(edgeIndex);
		temporalEdgePartClass.push_back(partClass);
    }

    ilpSolver_.initModel(coefficients.size(), coefficients.data());
    std::cout << "   " << numberOfVariablesX_ << " variables x" << std::endl;
    std::cout << "   " << numberOfVariablesY_ << " variables y" << std::endl;
    std::cout << "   " << numberOfVariablesT_ << " variables t" << std::endl;
}

template<class ILP_SOLVER, class FEASBILE_SOLUTION_CALLBACK>
typename PoseTracker<ILP_SOLVER, FEASBILE_SOLUTION_CALLBACK>::size_type
PoseTracker<ILP_SOLVER, FEASBILE_SOLUTION_CALLBACK>::addPartialSolutionConstraints(const SolutionType &partial) {

    std::cout << "solution and problem size: " << partial.size() << " " << problem_.numberOfDetections() << std::endl;

    const size_type max_val = std::numeric_limits<size_type>::max();

    size_type n = 0;

    for(size_type d = 0; d < problem_.numberOfDetections(); ++d) { //problem_.numberOfDetections()
        auto det = partial[d];
        auto part_status = det.status_;
        auto cluster = det.clusterIndex_;
        bool suppressed = (part_status == 0 && cluster != max_val);

        //std::cout << "detection " << d << " suppressed " << suppressed << std::endl;
        const double bound = suppressed ? 0.0 : 1.0;

        if(part_status == 1 || suppressed)
		{
        	const size_type vi[] = {d};
			const double coef[] = {1.0};
			ilpSolver_.addConstraint(vi, vi + 1, coef, bound, bound); // x_{d, c} = 1
			//std::cout << "setting constraint " << d << " c " << c << " " << bound << std::endl;
			n++;
        }
    }
    numberOfFixedVariables_ = n;

    std::cout << "adding partial solution constraints: " << std::flush;

    std::cout << n << " ";
    n = 0;

    for(size_type d0 = 0; d0 < problem_.numberOfDetections()-1; ++d0) {
        auto cluster0 = partial[d0].clusterIndex_;
        auto part_status0 = partial[d0].status_;
        if (cluster0 == max_val || part_status0 == max_val) {
            //std::cout << "detection " << d0 << " is unlabeled" << std::endl;
            continue;
        }
        for(size_type d1 = d0+1; d1 < problem_.numberOfDetections(); ++d1) {
            auto cluster1 = partial[d1].clusterIndex_;
            auto part_status1 = partial[d1].status_;
            if (cluster1 == max_val || part_status1 == max_val) {
                continue;
            }
            //std::cout << "yconstr " << d0 << " " << d1 << " " << std::endl;
            std::pair<bool, size_type> e = detectionGraph_.findEdge(d0,d1);
            if(e.first){
            	const size_type vi[] = {getIndexOffset() + e.second};
            	const double c[] = {1.0};
            	const double bound = (cluster0 == cluster1) ? 1.0 : 0.0;
            	ilpSolver_.addConstraint(vi, vi + 1, c, bound, bound); // y_{d0, d1} = 1 or 0
            	//std::cout<<d0<<"\t"<<d1<<"\t"<<bound<<std::endl;
            	n++;
            }
        }
    }
    std::cout << n << std::endl;
    return n;
}


template<class ILP_SOLVER, class FEASBILE_SOLUTION_CALLBACK>
typename PoseTracker<ILP_SOLVER, FEASBILE_SOLUTION_CALLBACK>::size_type
PoseTracker<ILP_SOLVER, FEASBILE_SOLUTION_CALLBACK>::y(
    const size_type d0,
    const size_type d1
) const {
    assert(d0 < problem_.numberOfDetections());
    assert(d1 < problem_.numberOfDetections());
    std::pair<bool, size_type> edge = detectionGraph_.findEdge(d0, d1);
    assert(edge.first);
    return  getIndexOffset() + edge.second;
}

template<class ILP_SOLVER, class FEASBILE_SOLUTION_CALLBACK>
typename PoseTracker<ILP_SOLVER, FEASBILE_SOLUTION_CALLBACK>::size_type
PoseTracker<ILP_SOLVER, FEASBILE_SOLUTION_CALLBACK>::getIndexOffset() const {
    return  numberOfVariablesX_;
}

template<class ILP_SOLVER, class FEASBILE_SOLUTION_CALLBACK>
typename PoseTracker<ILP_SOLVER, FEASBILE_SOLUTION_CALLBACK>::size_type
PoseTracker<ILP_SOLVER, FEASBILE_SOLUTION_CALLBACK>::addAllCouplingConstraints()
{
    std::cout << "adding all coupling constraints: " << std::flush;

    size_t n = 0;

    size_t vi[] = { 0, 0 };
    double const coefficients[] = {1.0, -1.0};

    double const lowerBound = 0.0;
    double const upperBound = std::numeric_limits<double>::infinity();

    for (size_t d0 = 0; d0 < problem_.numberOfDetections(); ++d0)
        for (size_t d1 = d0 + 1; d1 < problem_.numberOfDetections(); ++d1)
        {
        	std::pair<bool, size_type> e = detectionGraph_.findEdge(d0, d1);
            if(e.first){
            	vi[1] = getIndexOffset() + e.second;

            	vi[0] = d0; // d0
            	ilpSolver_.addConstraint(vi, vi+2, coefficients, lowerBound, upperBound); // 0 <= (x_{d0, c} - y_{d0, d1}
            	++n;

            	vi[0] = d1; // d1
            	ilpSolver_.addConstraint(vi, vi+2, coefficients, lowerBound, upperBound); // 0 <= x_{d1, c} - y_{d0, d1}
            	++n;
            }
            else{
            }
        }

    std::cout << n << std::endl;

    return n;
}


template<class ILP_SOLVER, class FEASBILE_SOLUTION_CALLBACK>
typename PoseTracker<ILP_SOLVER, FEASBILE_SOLUTION_CALLBACK>::size_type
PoseTracker<ILP_SOLVER, FEASBILE_SOLUTION_CALLBACK>::addAll3CycleConstraints()
{
    std::cout << "adding all 3 cycle constraints: " << std::flush;

    size_t n = 0;

    double const coefficients[] = {-1.0, 1.0, 1.0};
    double const lowerBound = -std::numeric_limits<double>::infinity();
    double const upperBound = 1.0;

    size_t vi[] = { 0, 0, 0 };
    for (size_t edge = 0; edge < detectionGraph_.numberOfEdges(); ++edge)
    {
        size_t const v0 = detectionGraph_.vertexOfEdge(edge, 0);
        size_t const v1 = detectionGraph_.vertexOfEdge(edge, 1);

        if(v0 < numberOfFixedVariables_ && v1 < numberOfFixedVariables_) {continue;}

        vi[0] = y(v0, v1);

		for (size_t v2 = 0; v2 < detectionGraph_.numberOfVertices(); ++v2)
		{

			std::pair<bool, size_type> e1 = detectionGraph_.findEdge(v0,v2);
			std::pair<bool, size_type> e2 = detectionGraph_.findEdge(v1,v2);

			if(e1.first && e2.first)
			{
				vi[1] = getIndexOffset() + e1.second;
				vi[2] = getIndexOffset() + e2.second;

				ilpSolver_.addConstraint(vi, vi + 3, coefficients, lowerBound, upperBound);
				++n;
            }
		}
    }
    std::cout << n << std::endl;
    return n;
}


template<class ILP_SOLVER, class FEASBILE_SOLUTION_CALLBACK>
typename PoseTracker<ILP_SOLVER, FEASBILE_SOLUTION_CALLBACK>::size_type
PoseTracker<ILP_SOLVER, FEASBILE_SOLUTION_CALLBACK>::addAllSpatioTemporalConstraints()
{
    std::cout << "adding all spatio temporal constraints: " << std::flush;

    size_t n = 0;

    double const c[] = { 1.0, 1.0, 1.0, -1.0 };
    double const lowerBound = -std::numeric_limits<double>::infinity();
    double const upperBound = 2.0;

    size_t vi[] = { 0, 0, 0, 0};
    for (size_t e0 = 0; e0 < temporalEdgeIndexes.size()-1; ++e0)
    {
    	size_t edgeIndex0 = temporalEdgeIndexes[e0];
    	size_t partClass0 = temporalEdgePartClass[e0];
    	size_t edge0 = edgeIndex0 - getIndexOffset();
        size_t const v0 = detectionGraph_.vertexOfEdge(edge0, 0);
        size_t const v1 = detectionGraph_.vertexOfEdge(edge0, 1);

        for (size_t e1 = e0+1; e1 < temporalEdgeIndexes.size(); ++e1)
        {
        	size_t edgeIndex1 = temporalEdgeIndexes[e1];
        	size_t partClass1 = temporalEdgePartClass[e1];

        	if(partClass0 == partClass1){ continue; }

        	size_t edge1 = edgeIndex1 - getIndexOffset();
            size_t const v2 = detectionGraph_.vertexOfEdge(edge1, 0);
            size_t const v3 = detectionGraph_.vertexOfEdge(edge1, 1);

            std::pair<bool, size_type> y0 = detectionGraph_.findEdge(v0,v2);
            if(!y0.first) { continue; }
            std::pair<bool, size_type> y1 = detectionGraph_.findEdge(v1,v3);
            if(!y1.first) { continue; }

            vi[0] = edgeIndex0;
            vi[1] = edgeIndex1;
            vi[2] = getIndexOffset() + y0.second;
            vi[3] = getIndexOffset() + y1.second;
            ilpSolver_.addConstraint(vi, vi + 4, c, lowerBound, upperBound);
			++n;

            vi[2] = getIndexOffset() + y1.second;
            vi[3] = getIndexOffset() + y0.second;
            ilpSolver_.addConstraint(vi, vi + 4, c, lowerBound, upperBound);
			++n;

		/*	//std::cout<<temporalEdgePartClass[e0]<<"\t"<<temporalEdgePartClass[e1]<<std::endl;
			vi[0] = getIndexOffset() + y0.second;
			vi[1] = getIndexOffset() + y1.second;
			ilpSolver_.addConstraint(vi, vi + 2, coefficients, lowerBound, 0.);
			++n;

			vi[0] = edgeIndex0;
			ilpSolver_.addConstraint(vi, vi + 2, coefficients, lowerBound, 0.);
			++n;

			vi[0] = edgeIndex1;
			ilpSolver_.addConstraint(vi, vi + 2, coefficients, lowerBound, 0.);
			++n;*/
        }
    }
    std::cout << n << std::endl;
    return n;
}




// Callback

template<class ILP_SOLVER, class FEASBILE_SOLUTION_CALLBACK>
inline
PoseTracker<ILP_SOLVER, FEASBILE_SOLUTION_CALLBACK>::Callback::Callback(
	PoseTrackerType& poseTracker
)
:   IlpSolverCallback(poseTracker.ilpSolver_),
    poseTracker_(poseTracker)
{}

template<class ILP_SOLVER, class FEASBILE_SOLUTION_CALLBACK>
inline void
PoseTracker<ILP_SOLVER, FEASBILE_SOLUTION_CALLBACK>::Callback::separateAndAddLazyConstraints()
{

	poseTracker_.timer_.stop();
	const double elapsedSeconds = poseTracker_.timer_.elapsedSeconds();
	poseTracker_.timer_.start();

	std::cout << elapsedSeconds
		<< '\t' << IlpSolverCallback::objectiveBound()
		<< '\t' << IlpSolverCallback::objectiveBest()
		<< std::flush<<std::endl;

	//size_t const nCycle = separateAndAddViolated3CycleConstraints();
	//std::cout << '\t' << nCycle << std::flush <<std::endl;

  /*  if (nCycle == 0)
    {
        // save intermediate feasible solution
        ComponentsType components;
        components.build(poseTracker_.detectionGraph_, SubgraphMask(*this));

        typename PoseTrackerType::SolutionType solution(poseTracker_.problem_.numberOfDetections());
        for (size_t d = 0; d < poseTracker_.problem_.numberOfDetections(); ++d)
        {
            solution[d].clusterIndex_ = components.labels_[d];
        }
        poseTracker_.feasibleSolutionCallback_(solution);
    }
    std::cout << std::endl;*/
}


template<class ILP_SOLVER, class FEASBILE_SOLUTION_CALLBACK>
typename PoseTracker<ILP_SOLVER, FEASBILE_SOLUTION_CALLBACK>::Callback::size_type
PoseTracker<ILP_SOLVER, FEASBILE_SOLUTION_CALLBACK>::Callback::separateAndAddViolated3CycleConstraints() {
    size_t n = 0;

    double const coefficients[] = {-1.0, 1.0, 1.0};
    double const lowerBound = -std::numeric_limits<double>::infinity();
    double const upperBound = 1.0;

    size_t vi[] = { 0, 0, 0 };
    for (size_t edge = 0; edge < poseTracker_.detectionGraph_.numberOfEdges(); ++edge)
    {
        size_t const v0 = poseTracker_.detectionGraph_.vertexOfEdge(edge, 0);
        size_t const v1 = poseTracker_.detectionGraph_.vertexOfEdge(edge, 1);

        if(v0 < poseTracker_.numberOfFixedVariables_ && v1 < poseTracker_.numberOfFixedVariables_)
        	continue;

        vi[0] = poseTracker_.y(v0, v1);

        if (this->label(vi[0]) < .5) // cut
            for (size_t v2 = 0; v2 < poseTracker_.detectionGraph_.numberOfVertices(); ++v2)
            {
            	std::pair<bool, size_type> e1 = poseTracker_.detectionGraph_.findEdge(v2,v0);
            	std::pair<bool, size_type> e2 = poseTracker_.detectionGraph_.findEdge(v2,v1);

            	if(e1.first && e2.first)
            	{
            		vi[1] = poseTracker_.getIndexOffset() + e1.second;
            		vi[2] = poseTracker_.getIndexOffset() + e2.second;

            		if (this->label(vi[1]) > .5 && this->label(vi[2]) > .5) // join, join
            		{
						this->addLazyConstraint(vi, vi + 3, coefficients, lowerBound, upperBound);
						++n;
            		}
            	}
            }
    }
    return n;
}

} // namespace pt

#endif
