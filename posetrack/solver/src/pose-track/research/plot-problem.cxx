#include <iostream>

#include "andres/graphics/graphics-hdf5.hxx"

#include "io-hdf5.hxx"

typedef andres::graphics::Graphics<> Graphics;

int main(
    int argc,
    char** argv
) {
    if(argc != 5) {
        std::cerr << "usage: plot-problem <problem.h5> <graphics.h5> <class0> <class1>" << std::endl;
        return 1;
    }

    const std::string problemFileName = argv[1];
    const std::string graphicsFileName = argv[2];
    const size_t partClass0 = atoi(argv[3]);
    const size_t partClass1 = atoi(argv[4]);
    const double minProbability = 0.5;

    // load problem (probabilities)
    Problem problem;
    loadProblem(problemFileName, problem);

    // load problem (geometry)
    andres::Marray<double> coordinates;
    {
        hid_t file = andres::hdf5::openFile(problemFileName);
        andres::hdf5::load(file, "coordinates-vertices", coordinates);
        andres::hdf5::closeFile(file);
    }
    std::cout << coordinates.shape(0)
        << " " << coordinates.shape(1)
        << std::endl;

    // draw graphics
    Graphics graphics;
    for(size_t j = 0; j < problem.numberOfDetections(); ++j) {
        graphics.definePoint(coordinates(j, 0) / 640.0, coordinates(j, 1)  / 640.0, 0);
    }
    for(auto it = problem.joinMap().begin(); it != problem.joinMap().end(); ++it) {
        const pose::JoinIndex<>& index = it->first;
        const pose::JoinData<>& data = it->second;

        const size_t d0 = index.getDetection(0);
        const size_t c0 = index.getPartClass(0);

        const size_t d1 = index.getDetection(1);
        const size_t c1 = index.getPartClass(1);

        if(c0 == partClass0 && c1 == partClass1
        || c0 == partClass1 && c1 == partClass0) {
            const double p = data.getProbability() // pairwise
                * problem.getPartClassProbability(d0, c0) // unary
                * problem.getPartClassProbability(d1, c1); // unary
            if(p > minProbability) {
                unsigned char gray = p * 255.0;
                size_t property = graphics.defineLineProperty(true, gray, 255 - gray, 0);
                graphics.defineLine(d0, d1, property);
            }
        }
    }

    // save graphics
    {
        hid_t file = andres::graphics::hdf5::createFile(graphicsFileName);
        andres::graphics::hdf5::save(file, graphics);
        andres::graphics::hdf5::closeFile(file);
    }

    return 0;
}
