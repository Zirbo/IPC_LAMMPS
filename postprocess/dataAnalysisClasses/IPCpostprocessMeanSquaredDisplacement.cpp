#include "IPCpostprocessMeanSquaredDisplacement.hpp"

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
    double meanSquaredDisplacement = 0.0;

    for (IPC ipc: system) {
        for (int d: DIMENSIONS) {
            double delta_xj = ipc.ipcCenter.x[d] - ipcCentersPreviousPositions[ipc.number][d];
            relativePBC(delta_xj);
            displacementOfEachIPCs[ipc.number][d] += boxSide[d]*delta_xj;
            meanSquaredDisplacement += std::pow(displacementOfEachIPCs[ipc.number][d],2);
        }
    }
    meanSquaredDisplacement /= (int)system.size();
    meanSquaredDisplFile << timeCounter << "\t" << meanSquaredDisplacement << "\n";
}

void IPCmsd::updatePreviousPositions(const Ensemble &system) {
    // set the current as the previous iteration state
    for (IPC ipc: system) {
        for (int d: DIMENSIONS) {
            ipcCentersPreviousPositions[ipc.number][d] = ipc.ipcCenter.x[d];
        }
    }
}
