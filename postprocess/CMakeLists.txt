cmake_minimum_required(VERSION 3.5)
set (CMAKE_CXX_STANDARD 11)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")


project(lammpsIPCpostprocess)

set(ALL_FILES
    IPC.hpp
    IPCpostprocess.cpp
    IPCpostprocess.hpp
    main.cpp
    IPCpostprocessPotential.hpp
    IPCpostprocessPotential.cpp
    dataAnalysisClasses/IPCpostprocessNeighbourAnalysis.hpp
    dataAnalysisClasses/IPCpostprocessNeighbourAnalysis.cpp
    dataAnalysisClasses/IPCpostprocessOrientationsAnalysis.hpp
    dataAnalysisClasses/IPCpostprocessOrientationsAnalysis.cpp
    dataAnalysisClasses/IPCpairCorrelation.hpp
    dataAnalysisClasses/IPCpairCorrelation.cpp
    dataAnalysisClasses/IPCpostprocessEccentricityHistogram.hpp
    dataAnalysisClasses/IPCpostprocessEccentricityHistogram.cpp
    dataAnalysisClasses/IPCpostprocessMeanSquaredDisplacement.hpp
    dataAnalysisClasses/IPCpostprocessMeanSquaredDisplacement.cpp
    dataAnalysisClasses/IPCpostprocessAxialityDeviationHistogram.hpp
    dataAnalysisClasses/IPCpostprocessAxialityDeviationHistogram.cpp
)

add_executable(${PROJECT_NAME} ${ALL_FILES})
