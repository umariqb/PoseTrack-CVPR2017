#include "andres/marray-hdf5.hxx"

#include "pose-track/problem.hxx"
#include "pose-track/solution.hxx"

typedef double value_type;
typedef std::size_t size_type;
typedef pt::Problem<value_type, size_type> Problem;
typedef pt::Solution<size_type> Solution;

void
loadProblem(
    const std::string problemFileName,
    Problem& problem
) {
    std::cout << "loading problem: " << std::endl;

    andres::Marray<value_type> numDetectionsInfoTable;
    andres::Marray<value_type> partClassProbabilityTable;
    andres::Marray<value_type> joinProbabilityTableTemporal;
    andres::Marray<value_type> joinProbabilityTableSpatial;
    {
        hid_t problemFile = andres::hdf5::openFile(problemFileName);
        andres::hdf5::load(problemFile, "detections-info", numDetectionsInfoTable);
        andres::hdf5::load(problemFile, "part-class-probabilities", partClassProbabilityTable);
        andres::hdf5::load(problemFile, "join-probabilities-temporal", joinProbabilityTableTemporal);
        andres::hdf5::load(problemFile, "join-probabilities-spatial", joinProbabilityTableSpatial);
        andres::hdf5::closeFile(problemFile);
    }

    size_type numberOfDetections = numDetectionsInfoTable(0);

    if(partClassProbabilityTable.dimension() != 2) {
        throw std::runtime_error("partClassProbabilityTable.dimension() != 2");
    }

    if(joinProbabilityTableTemporal.dimension() != 2) {
        throw std::runtime_error("joinProbabilityTableTemporal.dimension() != 2");
    }

    if(joinProbabilityTableSpatial.dimension() != 2) {
            throw std::runtime_error("joinProbabilityTableSpatial.dimension() != 2");
    }


    if(joinProbabilityTableTemporal.shape(1) != 4) {
        throw std::runtime_error("joinProbabilityTableTemporal.shape(1) != 4\n \expecting a 2-dimensional matrix with 4 columns: detection0, detection1, partClass, pairwiseProb");
    }

    if(joinProbabilityTableSpatial.shape(1) != 5) {
        throw std::runtime_error("joinProbabilityTableSpatial.shape(1) != 5\n \expecting a 2-dimensional matrix with 5 columns: detection0, detection1, partClass1, partClass2, pairwiseProb");
    }

    problem.assign(numberOfDetections);

    // set unary probabilities
    for(size_type d0 = 0; d0 < numberOfDetections; ++d0){
        problem.setPartClassProbability(d0, partClassProbabilityTable(d0,0), partClassProbabilityTable(d0,1));
    }

    // set binaries for temporal connections
    for(size_type j = 0; j < joinProbabilityTableTemporal.shape(0); ++j) {
        problem.setJoinProbability(
			joinProbabilityTableTemporal(j, 0),
			joinProbabilityTableTemporal(j, 1),
			joinProbabilityTableTemporal(j, 2),
			joinProbabilityTableTemporal(j, 2),
			joinProbabilityTableTemporal(j, 3),
			false
        );
    }

    // set binaries for spatial connections
    for(size_type j = 0; j < joinProbabilityTableSpatial.shape(0); ++j) {
        problem.setJoinProbability(
			joinProbabilityTableSpatial(j, 0),
			joinProbabilityTableSpatial(j, 1),
			joinProbabilityTableSpatial(j, 2),
			joinProbabilityTableSpatial(j, 3),
			joinProbabilityTableSpatial(j, 4),
			true
        );
    }

    std::cout  << numberOfDetections << " total number of detections." <<
    	joinProbabilityTableTemporal.shape(0)+joinProbabilityTableSpatial.shape(0) <<
    	" total number of edges in the graph." << std::endl;

}

void
saveSolution(
    const std::string solutionFileName,
    const Solution& solution
) {
    size_type shape[] = {solution.size(), 2};
    andres::Marray<size_type> m(shape, shape + 2);
    for(size_type d = 0; d < solution.size(); ++d) {
        m(d, 0) = solution[d].status_;
        m(d, 1) = solution[d].clusterIndex_;
    }

    /*
    std::cout << "solution:" << std::endl
        << std::cout << m.transposedView().asString() << std::endl;
    */

    hid_t file = andres::hdf5::createFile(solutionFileName);
    andres::hdf5::save(file, "part-tracks", m);
    andres::hdf5::closeFile(file);
}

void
loadSolution(
    const std::string solutionFileName,
    Solution& solution
) {
    typedef pt::Detection<size_type> Detection;
    andres::Marray<size_type> m;
    {
        hid_t file = andres::hdf5::openFile(solutionFileName);
        andres::hdf5::load(file, "detection-parts-and-clusters", m);
        andres::hdf5::closeFile(file);
    }

    const size_type numberOfDetections = m.shape(0);
    solution.resize(numberOfDetections);
    for(size_type d = 0; d < numberOfDetections; ++d) {
        solution[d].status_ = m(d, 0);
        solution[d].clusterIndex_ = m(d, 1);
    }
}


