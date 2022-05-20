#ifndef __IPCPOSTPROCESS_HEADER_INCLUDED__
#define __IPCPOSTPROCESS_HEADER_INCLUDED__

/*---------------------------------------------------------------------------------------
 * Postprocess for Inverse Patchy Colloids simulations with variable number of patches.
 *
 * asda
 *---------------------------------------------------------------------------------------*/
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <limits>

#include "IPC.hpp"
#include "IPCpostprocessPotential.hpp"

struct IPCpostprocessInitializer {
    std::string trajFilename;
    std::string inputFilename;
    std::string potDirName;
    int startStep;
    int finalStep;

    IPCpostprocessInitializer() {
        startStep = 0;
        finalStep = std::numeric_limits<int>::max();
    }
};

class IPCpostprocess {
public:
    IPCpostprocess(IPCpostprocessInitializer const& initializer);
    void run();


private:
    IPCpostprocess();

    // input trajectory
    std::ifstream trajectoryFile;
    // state point
    int nIPCs;
    Triad boxSide;
    Ensemble ipcs;
    VectorOfTriads ipcOrientationsFirstPatch;
    VectorOfTriads ipcOrientationsSecndPatch;
    std::vector<double> ipcEccentricities;
    // geometry
    double ipcRadius, patchRadius, patchEccentricity, interactionRange;

    int startStep, finalStep, currentStep;

    IPCpotential potential;

    void readFirstConfiguration();
    bool readNewConfiguration();
    void readIPCconfiguration();
    inline void relativePBC(double & x) {  x -= std::round(x);  }
    void computeOrientations();
    //void computePerfectOrientations();
    void printFinalOrientations(std::string const& outFileName);
};

#endif //__IPCPOSTPROCESS_HEADER_INCLUDED__
