#include "IPCpostprocessMeanSquaredDisplacement.hpp"
#include <iostream>

IPCmsd::IPCmsd(std::string const& outputFileName, Ensemble const& system, Triad const& pboxSide) {
    meanSquaredDisplFile.open(outputFileName);
    ipcCentersPreviousPositions.resize(system.size(), {0.0, 0.0, 0.0});
    updatePreviousPositions(system);
    displacementOfEachIPCs.resize(system.size(), {0.0, 0.0, 0.0});
    timeCounter = 0;

    boxSide = pboxSide;
}

void IPCmsd::accumulateAndPrint(const Ensemble &system) {
    computeMSD(system);
    updatePreviousPositions(system);
}


void IPCmsd::computeMSD(const Ensemble &system) {
    ++timeCounter;
    //std::cout << timeCounter << std::endl;
    Triad meanSquaredDisplacement_d = {0.0, 0.0, 0.0};

    for (IPC ipc: system) {
        for (int d: DIMENSIONS) {
            double delta_xd = ipc.ipcCenter.x[d] - ipcCentersPreviousPositions[ipc.number][d];
            relativePBC(delta_xd);
            displacementOfEachIPCs[ipc.number][d] += delta_xd;
            meanSquaredDisplacement_d[d] += std::pow(displacementOfEachIPCs[ipc.number][d],2);
        }
    }
    double meanSquaredDisplacement = 0.0;
    for (int d: DIMENSIONS) {
        meanSquaredDisplacement += boxSide[d]*meanSquaredDisplacement_d[d];
    }
    meanSquaredDisplacement /= (int)system.size();

    meanSquaredDisplFile << timeCounter << "\t" << meanSquaredDisplacement << std::endl;
}

void IPCmsd::updatePreviousPositions(const Ensemble &system) {
    // set the current as the previous iteration state
    for (IPC ipc: system) {
        for (int d: DIMENSIONS) {
            ipcCentersPreviousPositions[ipc.number][d] = ipc.ipcCenter.x[d];
        }
    }
}
