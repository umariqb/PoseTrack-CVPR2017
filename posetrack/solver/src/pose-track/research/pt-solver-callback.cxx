#include <cstddef>
#include <string>
#include <iostream>
#include <stdexcept>

#include "andres/ilp/gurobi-callback.hxx"

#include "pt-io-hdf5.hxx"

#include "pose-track/pt-solver-callback.hxx"

struct FeasibleSolutionCallback {
    typedef std::size_t size_type;
    typedef pt::Solution<size_type> SolutionType;

    FeasibleSolutionCallback(const std::string& fileNamePrefix)
        :   fileNamePrefix_(fileNamePrefix),
            counter_()
        {}
    void operator()(const SolutionType& solution)
        {
            std::stringstream fileName;
            fileName << fileNamePrefix_ << "-INTERMEDIATE-" << counter_ << ".h5";
            saveSolution(fileName.str(), solution);
            ++counter_;
        }

    std::string fileNamePrefix_;
    size_type counter_;
};

typedef andres::ilp::Gurobi<> IlpSolver;
typedef pt::PoseTracker<IlpSolver, FeasibleSolutionCallback> Solver;

int main(int argc, char** argv) {
    int timeLimit = 36000;
    if(argc < 4) {
        std::cerr << "solver <problem.h5> <solution.h5> s/m" << std::endl;
        return 1;
    }
    
    if (argc >= 5) {
        timeLimit = atoi(argv[4]);
    }

    std::string partialSolutionFileName;

    bool stagewise = false;
    if (argc >= 6) {
        partialSolutionFileName = argv[5];
        stagewise = true;
    }

    const std::string problemFileName = argv[1];
    const std::string solutionFileName = argv[2];
    const std::string solverIndex = argv[3];

    if(solverIndex.size() != 1 || (solverIndex[0] != 's' && solverIndex[0] != 'm')) {
        std::cerr << "3rd argument must be s (single clusters) or m (multiple clusters)" << std::endl;
        return 1;
    }

    const bool withSingleClusterConstraints = (solverIndex[0] == 's');

    try {
        Problem problem;
        loadProblem(problemFileName, problem);

        Solution partialSolution;
        if(stagewise) {
            loadSolution(partialSolutionFileName, partialSolution);
        }

        Solution solution;
        FeasibleSolutionCallback feasibleSolutionCallback(solutionFileName);

        Solver solver(problem, solution, partialSolution, feasibleSolutionCallback, withSingleClusterConstraints, timeLimit);


        std::cout<<"saving solution."<<std::endl;
        saveSolution(solutionFileName, solution);
    } catch(const GRBException& exc) {
        std::cout << "GRBException " << exc.getMessage() << std::endl;
    }

    return 0;
}
