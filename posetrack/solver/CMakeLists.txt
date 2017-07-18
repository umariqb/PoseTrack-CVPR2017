cmake_minimum_required(VERSION 2.6)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)
project(tracking)
include_directories(graph/include graphics/include include)
file(GLOB headers
    include/andres/*.hxx
    include/andres/graph/*.hxx
    include/andres/graph/hdf5/*.hxx
    include/andres/graphics/*.hxx
    include/andres/ilp/*.hxx
    include/pose/*.hxx
)

file(GLOB pt-headers
    include/andres/*.hxx
    include/andres/graph/*.hxx
    include/andres/graph/hdf5/*.hxx
    include/andres/graphics/*.hxx
    include/andres/ilp/*.hxx
    include/pose-track/*.hxx
)


enable_testing()

##############################################################################
# GNU g++ specific settings
##############################################################################
if(CMAKE_COMPILER_IS_GNUCXX)
    add_definitions(-std=c++11)
endif()

##############################################################################
# OpenMP
##############################################################################
find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    message(STATUS "OpenMP found")
else()
    message("OpenMP not found")
endif()

##############################################################################
# HDF5
##############################################################################
find_package(HDF5)
if(HDF5_FOUND)
    include_directories(${HDF5_INCLUDE_DIRS})
endif()

##############################################################################
# Gurobi
##############################################################################
find_package(GUROBI)
if(GUROBI_FOUND)
    include_directories(${GUROBI_INCLUDE_DIR})
endif()

##############################################################################
# GLUT
##############################################################################
find_package(GLUT)


##############################################################################
# targets: pose
##############################################################################

# research
if(HDF5_FOUND)
    if(GUROBI_FOUND)
        add_executable(pt-solver-callback src/pose-track/research/pt-solver-callback.cxx ${pt-headers})
        target_link_libraries(pt-solver-callback ${HDF5_LIBRARIES} ${GUROBI_LIBRARIES})
    endif()
endif()
