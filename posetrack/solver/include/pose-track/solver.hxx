#pragma once
#ifndef POSE_TRACK_SOLVER_HXX
#define POSE_TRACK_SOLVER_HXX

#include <cstddef>
#include <limits>
#include <iostream>
#include <vector>
#include <string>

#include <andres/graph/complete-graph.hxx>
#include <andres/graph/components.hxx>

#include "problem.hxx"
#include "solution.hxx"

namespace pt {

template<class ILP_SOLVER>
class Solver {
public:
    typedef ILP_SOLVER IlpSolver;
    typedef std::size_t size_type;
    typedef Problem<double, size_type> ProblemType;
    typedef typename ProblemType::JoinIndexType JoinIndexType;
    typedef typename ProblemType::JoinDataType JoinDataType;
    typedef typename ProblemType::JoinMapType JoinMapType;
    typedef Solution<size_type> SolutionType;

    Solver(ProblemType&, SolutionType&, const std::string&, const size_t timeLimit = 36000, const double = 1.0 / 255.0);

private:
    typedef typename andres::graph::CompleteGraph<> CompleteGraphType;
    typedef andres::graph::ComponentsBySearch<CompleteGraphType> ComponentsType;

    struct SubgraphMask {
        typedef Solver<IlpSolver> SolverType;

        SubgraphMask(const SolverType& solver)
            : solver_(solver) {}
        bool vertex(const std::size_t v) const
            { return true; }
        bool edge(const std::size_t e) const
            {
                const size_type v0 = solver_.detectionGraph_.vertexOfEdge(e, 0);
                const size_type v1 = solver_.detectionGraph_.vertexOfEdge(e, 1);
                return solver_.ilpSolver_.label(solver_.y(v0, v1)) == 1;
            }
        const SolverType& solver_;
    };

    // ilp
    void setObjectiveFunction();
    size_type addAllImpossiblePartClassConstraints();
    size_type separateAndAddViolatedUniquenessConstraints();
    size_type separateAndAddViolatedSingleClusterConstraints();
    size_type separateAndAddViolatedCouplingConstraints();
    size_type separateAndAddViolated3CycleConstraints();
    size_type addAllViolatedLinearizationConstraints();
    size_type separateAndAddViolatedMustSelectClassConstraints();

    // variable indexing
    size_type x(const size_type, const size_type) const;
    size_type y(const size_type, const size_type) const;

    size_type numberOfVariablesX_;
    size_type numberOfVariablesY_;
    size_type numberOfVariablesZ_;
    ProblemType& problem_;
    double epsilonProbability_;
    CompleteGraphType detectionGraph_;
    IlpSolver ilpSolver_;
    std::string solverIndex_;

friend struct SubgraphMask;
};

template<class ILP_SOLVER>
Solver<ILP_SOLVER>::Solver(
    ProblemType& problem,
    SolutionType& solution,
    const std::string& solverIndex,
    const size_t timeLimit,
    const double epsilonProbability
)
:   numberOfVariablesX_(),
    numberOfVariablesY_(),
    numberOfVariablesZ_(),
    problem_(problem),
    epsilonProbability_(epsilonProbability),
    detectionGraph_(problem.numberOfDetections()),
    ilpSolver_(),
    solverIndex_(solverIndex)
{
    ilpSolver_.setVerbosity(true);
    ilpSolver_.setRelativeGap(.0);//
    //ilpSolver_.setNumberOfThreads(2);
    ilpSolver_.setNumberOfThreads(1);
    //ilpSolver_.setTimeLimit(timeLimit);

    setObjectiveFunction();
    addAllImpossiblePartClassConstraints();
    addAllViolatedLinearizationConstraints();

    for(;;) { // cutting plane loop
        ilpSolver_.optimize();

        const size_type nUniqueness = separateAndAddViolatedUniquenessConstraints();
        const size_type nCoupling = separateAndAddViolatedCouplingConstraints();
        const size_type nCycle = separateAndAddViolated3CycleConstraints();
        //const size_type nMustSelectClass = separateAndAddViolatedMustSelectClassConstraints();

        // single cluster contraints
        if(solverIndex.compare("s") == 0){
            const size_type nSingleCluster = separateAndAddViolatedSingleClusterConstraints();
            if(nUniqueness + nCoupling + nCycle + nSingleCluster == 0) {
                break;
            }
        }

        // multiple clusters cases
        if(nUniqueness + nCoupling + nCycle == 0) {
            break;
        }
    }

    // save solution
    ComponentsType components;
    components.build(detectionGraph_, SubgraphMask(*this));
    solution.resize(problem_.numberOfDetections());
    for(size_type d = 0; d < problem.numberOfDetections(); ++d) {
       solution[d].clusterIndex_ = components.labels_[d];
       solution[d].partClass_ = std::numeric_limits<size_type>::max(); // suppressed
       for(size_type c = 0; c < problem.numberOfPartClasses(); ++c) {
            if(ilpSolver_.label(x(d, c)) == 1) {
                solution[d].partClass_ = c;
                break;
            }
       }
    }
}

template<class ILP_SOLVER>
void
Solver<ILP_SOLVER>::setObjectiveFunction() {
    std::cout << "setting up objective function:" << std::endl;

    numberOfVariablesX_ = problem_.numberOfDetections() * problem_.numberOfPartClasses();
    numberOfVariablesY_ = detectionGraph_.numberOfEdges();

    std::vector<double> coefficients(numberOfVariablesX_ + numberOfVariablesY_);

    // set coefficients of variables x
    for(size_type d = 0; d < problem_.numberOfDetections(); ++d)
    for(size_type c = 0; c < problem_.numberOfPartClasses(); ++c) {
        const double p = problem_.getPartClassProbability(d, c);
        if(p > epsilonProbability_) {
            const size_type vi = x(d, c);
            coefficients[vi] = std::log( (1.0 - p) / p );
        }
    }

    // introduce variables z and set their coefficients:
    // 1. the following code in known to work:
    for(size_type d0 = 0; d0 < problem_.numberOfDetections(); ++d0)
    for(size_type d1 = d0 + 1; d1 < problem_.numberOfDetections(); ++d1) {
        for(size_type c0 = 0; c0 < problem_.numberOfPartClasses(); ++c0)
        for(size_type c1 = 0; c1 < problem_.numberOfPartClasses(); ++c1) {
            const JoinIndexType joinIndex(d0, d1, c0, c1);
            auto it = problem_.joinMap().find(joinIndex);
            if(it != problem_.joinMap().end()) {
                const double p = it->second.getProbability();
                if(p > epsilonProbability_) {
                    const double c = std::log( (1.0 - p) / p );
                    it->second.setVariableIndex(coefficients.size());
                    coefficients.push_back(c);
                    ++numberOfVariablesZ_;
                }
            }
        }
    }
    // 2. the following code is experimental:
    /*
    for(size_type d0 = 0; d0 < problem_.numberOfDetections(); ++d0)
    for(size_type d1 = d0 + 1; d1 < problem_.numberOfDetections(); ++d1) {
        for(size_type c0 = 0; c0 < problem_.numberOfPartClasses(); ++c0)
        for(size_type c1 = 0; c1 < problem_.numberOfPartClasses(); ++c1) {
            const JoinIndexType joinIndex(d0, d1, c0, c1);
            auto it = problem_.joinMap().find(joinIndex);
            if(it != problem_.joinMap().end()) {
                const double p = it->second.getProbability();
                if(p > epsilonProbability_
                && problem_.getPartClassProbability(d0, c0) > epsilonProbability_
                && problem_.getPartClassProbability(d1, c1) > epsilonProbability_) {
                    const double c = std::log( (1.0 - p) / p );
                    it->second.setVariableIndex(coefficients.size());
                    coefficients.push_back(c);
                    ++numberOfVariablesZ_;
                }
            }
        }
    }
    */

    ilpSolver_.initModel(coefficients.size(), coefficients.data());

    std::cout << "   " << numberOfVariablesX_ << " variables x" << std::endl
        << "   " << numberOfVariablesY_ << " variables y" << std::endl
        << "   " << numberOfVariablesZ_ << " variables z not known to be zero" << std::endl;
}

template<class ILP_SOLVER>
typename Solver<ILP_SOLVER>::size_type
Solver<ILP_SOLVER>::addAllImpossiblePartClassConstraints() {
    std::cout << "adding all impossible part class constraints: " << std::flush;

    size_type n = 0;
    for(size_type d = 0; d < problem_.numberOfDetections(); ++d)
    for(size_type c = 0; c < problem_.numberOfPartClasses(); ++c) {
        const double p = problem_.getPartClassProbability(d, c);
        if(p <= epsilonProbability_) {
            const size_type vi[] = {x(d, c)};
            const double c[] = {1.0};
            const double lowerBound = 0.0;
            const double upperBound = 0.0;
            ilpSolver_.addConstraint(vi, vi + 1, c, lowerBound, upperBound); // x_{d, c} = 0
        }
    }

    std::cout << n << std::endl;

    return n;
}

template<class ILP_SOLVER>
typename Solver<ILP_SOLVER>::size_type
Solver<ILP_SOLVER>::separateAndAddViolatedUniquenessConstraints() {
    std::cout << "separating and adding violated uniqueness constraints: " << std::flush;

    size_type n = 0;
    size_type vi[] = {0, 0};
    const double c[] = {1.0, 1.0};
    const double lowerBound = -std::numeric_limits<double>::infinity();
    const double upperBound = 1.0;

    for(size_type d = 0; d < problem_.numberOfDetections(); ++d) {
        for(size_type c0 = 0; c0 < problem_.numberOfPartClasses(); ++c0) {
            vi[0] = x(d, c0);
            if(ilpSolver_.label(vi[0]) == 1) {
                for(size_type c1 = c0 + 1; c1 < problem_.numberOfPartClasses(); ++c1) {
                    vi[1] = x(d, c1);
                    if(ilpSolver_.label(vi[1]) == 1) {
                        ilpSolver_.addConstraint(vi, vi + 2, c, lowerBound, upperBound); // x_{d, c0} + x_{d, c1} <= 1
                        ++n;
                    }
                }
            }
        }
    }

    std::cout << n << std::endl;

    return n;
}

template<class ILP_SOLVER>
typename Solver<ILP_SOLVER>::size_type
Solver<ILP_SOLVER>::separateAndAddViolatedMustSelectClassConstraints() {
    std::cout << "separating and adding violated mustSelectClass constraints: " << std::flush;

    size_type n = 0;
    std::vector<size_type> vi(problem_.numberOfPartClasses());
    std::vector<double> c(problem_.numberOfPartClasses(), 1.0);
    const double lowerBound = 1;
    const double upperBound = std::numeric_limits<double>::infinity();
    for(size_type d0 = 0; d0 < problem_.numberOfDetections(); ++d0) {
            {
                double sum = 0;
                for(size_type c = 0; c < problem_.numberOfPartClasses(); ++c) {
                    vi[c] = x(d0, c); // d0
                    sum += ilpSolver_.label(vi[c]);
                }
                if(sum == 0) {
                    ilpSolver_.addConstraint(vi.begin(), vi.end(), c.begin(), lowerBound, upperBound); // 1 <= (sum_{c \in C} x_{d0, c})
                    ++n;
                }
            }
    }

    std::cout << n << std::endl;
    return n;
}

template<class ILP_SOLVER>
typename Solver<ILP_SOLVER>::size_type
Solver<ILP_SOLVER>::separateAndAddViolatedSingleClusterConstraints() {
    std::cout << "separating and adding violated single cluster constraints: " << std::flush;

    size_type n = 0;
    size_type vi[] = {0, 0, 0};
    const double c[] = {1.0, 1.0, -1.0};
    const double lowerBound = -std::numeric_limits<double>::infinity();
    const double upperBound = 1.0;
    for(size_type d0 = 0; d0 < problem_.numberOfDetections(); ++d0)
    for(size_type d1 = d0 + 1; d1 < problem_.numberOfDetections(); ++d1) {
        vi[2] = y(d0, d1);
        if(ilpSolver_.label(vi[2]) == 0) { // y_{d0, d1} = 0
            for(size_type c0 = 0; c0 < problem_.numberOfPartClasses(); ++c0) {
                vi[0] = x(d0, c0);
                if(ilpSolver_.label(vi[0]) == 1) { // x_{d0, c0} = 1
                    for(size_type c1 = 0; c1 < problem_.numberOfPartClasses(); ++c1) {
                        vi[1] = x(d1, c1);
                        if(ilpSolver_.label(vi[1]) == 1) { // x_{d1, c0} = 1
                            ilpSolver_.addConstraint(vi, vi + 3, c, lowerBound, upperBound); // x_{d0, c0} + x_{d1, c1} - y <= 1
                            ++n;
                        }
                    }
                }
            }
        }
    }

    std::cout << n << std::endl;

    return n;
}

template<class ILP_SOLVER>
typename Solver<ILP_SOLVER>::size_type
Solver<ILP_SOLVER>::separateAndAddViolatedCouplingConstraints() {
    std::cout << "separating and adding violated coupling constraints: " << std::flush;

    size_type n = 0;
    std::vector<size_type> vi(1 + problem_.numberOfPartClasses());
    std::vector<double> c(1 + problem_.numberOfPartClasses(), 1.0);
    c.back() = -1.0;
    const double lowerBound = 0.0;
    const double upperBound = std::numeric_limits<double>::infinity();
    for(size_type d0 = 0; d0 < problem_.numberOfDetections(); ++d0)
    for(size_type d1 = d0 + 1; d1 < problem_.numberOfDetections(); ++d1) {
        vi.back() = y(d0, d1);
        if(ilpSolver_.label(vi.back()) == 1) {
            // d0
            {
                double sum = 0;
                for(size_type c = 0; c < problem_.numberOfPartClasses(); ++c) {
                    vi[c] = x(d0, c); // d0
                    sum += ilpSolver_.label(vi[c]);
                }
                if(sum == 0) {
                    ilpSolver_.addConstraint(vi.begin(), vi.end(), c.begin(), lowerBound, upperBound); // 0 <= (sum_{c \in C} x_{d0, c}) - y_{d0, d1}
                    ++n;
                }
            }
            // d1
            {
                double sum = 0;
                for(size_type c = 0; c < problem_.numberOfPartClasses(); ++c) {
                    vi[c] = x(d1, c); // d1
                    sum += ilpSolver_.label(vi[c]);
                }
                if(sum == 0) {
                    ilpSolver_.addConstraint(vi.begin(), vi.end(), c.begin(), lowerBound, upperBound); // 0 <= (sum_{c \in C} x_{d1, c}) - y_{d0, d1}
                    ++n;
                }
            }
        }
    }

    std::cout << n << std::endl;

    return n;
}

// TODO: current implementation has run-time O(n^3). improve to O(n^2) by connected component labeling!
template<class ILP_SOLVER>
typename Solver<ILP_SOLVER>::size_type
Solver<ILP_SOLVER>::separateAndAddViolated3CycleConstraints() {
    std::cout << "separating and adding violated 3-cycle constraints: " << std::flush;

    std::size_t n = 0;
    std::size_t vi[] = {0, 0, 0};
    for(std::size_t j = 0; j < problem_.numberOfDetections(); ++j) {
        for(std::size_t k = j + 1; k < problem_.numberOfDetections(); ++k) {
            vi[0] = y(j, k);
            for(std::size_t l = k + 1; l < problem_.numberOfDetections(); ++l) {
                vi[1] = y(k, l);
                vi[2] = y(j, l);
                const double lowerBound = -1.0;
                const double upperBound = 1.0;//std::numeric_limits<double>::infinity();
                if(ilpSolver_.label(vi[0]) == 1) {
                    if(ilpSolver_.label(vi[1]) == 1) {
                        if(ilpSolver_.label(vi[2]) == 0) {
                            const double coefficients[] = {1.0, 1.0, -1.0};
                            ilpSolver_.addConstraint(vi, vi + 3, coefficients, lowerBound, upperBound);
                            ++n;
                        }
                    }
                    else {
                        if(ilpSolver_.label(vi[2]) == 1) {
                            const double coefficients[] = {1.0, -1.0, 1.0};
                            ilpSolver_.addConstraint(vi, vi + 3, coefficients, lowerBound, upperBound);
                            ++n;
                        }
                    }
                }
                else {
                    if(ilpSolver_.label(vi[1]) == 1 && ilpSolver_.label(vi[2]) == 1) {
                        const double coefficients[] = {-1.0, 1.0, 1.0};
                        ilpSolver_.addConstraint(vi, vi + 3, coefficients, lowerBound, upperBound);
                        ++n;
                    }
                }
            }
        }
    }

    std::cout << n << std::endl;

    return n;
}

template<class ILP_SOLVER>
typename Solver<ILP_SOLVER>::size_type
Solver<ILP_SOLVER>::addAllViolatedLinearizationConstraints() {
    std::cout << "adding all linearization constraints: " << std::flush;

    size_t n = 0;

    size_t vi[] = { 0, 0, 0, 0 };

    double const lowerBound = -std::numeric_limits<double>::infinity();
    double const upperBound = 2.0;
    for (size_t d0 = 0; d0 < problem_.numberOfDetections(); ++d0)
        for (size_t d1 = d0 + 1; d1 < problem_.numberOfDetections(); ++d1)
        {
            vi[0] = y(d0, d1);

            { // if y_{d0, d1} = 1
                double const c[] = { 1.0, 1.0, 1.0, -1.0 };

                for (size_t c0 = 0; c0 < problem_.numberOfPartClasses(); ++c0)
                {
                    vi[1] = x(d0, c0);

                    for (size_t c1 = 0; c1 < problem_.numberOfPartClasses(); ++c1)
                    {
                        vi[2] = x(d1, c1);

                        auto it = problem_.joinMap().find(JoinIndexType(d0, d1, c0, c1));

                        if (it->second.getProbability() <= epsilonProbability_) // if z_{d0, d1, c0, c1} is not explicitly in the ILP
                        {
                            // z_{d0, d1, c0, c1} = 0. Thus:
                            ilpSolver_.addConstraint(vi, vi + 3, c, lowerBound, upperBound);
                            ++n;
                        }
                        else // if z_{d0, d1, c0, c1} is explicitly in the ILP
                        {
                            vi[3] = it->second.getVariableIndex(); // z_{d0, d1, c0, c1}

                            ilpSolver_.addConstraint(vi, vi + 4, c, lowerBound, upperBound);
                            ++n;
                        }
                    }
                }
            }
            { // if y_{d0, d1} = 0
                const double c[] = {-1.0, 1.0};

                for (size_t c0 = 0; c0 < problem_.numberOfPartClasses(); ++c0)
                    for (size_t c1 = 0; c1 < problem_.numberOfPartClasses(); ++c1)
                    {
                        auto it = problem_.joinMap().find(JoinIndexType(d0, d1, c0, c1));
                        
                        if (it->second.getProbability() > epsilonProbability_) // if z_{d0, d1, c0, c1} is explicitly in the ILP
                        {
                            vi[1] = it->second.getVariableIndex(); // z_{d0, d1, c0, c1}

                            ilpSolver_.addConstraint(vi, vi + 2, c, lowerBound, .0);
                            ++n;
                        }
                    }
            }
        }

    double const c[] = { 1.0, -1.0 };
    for (size_t d0 = 0; d0 < problem_.numberOfDetections(); ++d0)
        for (size_t d1 = d0 + 1; d1 < problem_.numberOfDetections(); ++d1)
            for (size_t c0 = 0; c0 < problem_.numberOfPartClasses(); ++c0)
                for (size_t c1 = 0; c1 < problem_.numberOfPartClasses(); ++c1)
                {
                    auto it = problem_.joinMap().find(JoinIndexType(d0, d1, c0, c1));

                    if (it->second.getProbability() <= epsilonProbability_)
                        continue;

                    vi[0] = it->second.getVariableIndex();

                    vi[1] = x(d0, c0);

                    ilpSolver_.addConstraint(vi, vi + 2, c, lowerBound, .0);
                    ++n;

                    vi[1] = x(d1, c1);

                    ilpSolver_.addConstraint(vi, vi + 2, c, lowerBound, .0);
                    ++n;                    
                }

    std::cout << n << std::endl;

    return n;
}

template<class ILP_SOLVER>
typename Solver<ILP_SOLVER>::size_type
Solver<ILP_SOLVER>::x(
    const size_type d,
    const size_type c
) const {
    assert(d < problem_.numberOfDetections());
    assert(c < problem_.numberOfPartClasses());
    return c + d * problem_.numberOfPartClasses();
}

template<class ILP_SOLVER>
typename Solver<ILP_SOLVER>::size_type
Solver<ILP_SOLVER>::y(
    const size_type d0,
    const size_type d1
) const {
    assert(d0 < problem_.numberOfDetections());
    assert(d1 < problem_.numberOfDetections());
    return numberOfVariablesX_ + detectionGraph_.findEdge(d0, d1).second;
}

} // namespace pt

#endif
