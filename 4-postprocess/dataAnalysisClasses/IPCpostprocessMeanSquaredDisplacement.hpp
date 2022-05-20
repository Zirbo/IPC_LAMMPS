#ifndef __IPCPOSTPROCESS_MSD_HEADER_INCLUDED__
#define __IPCPOSTPROCESS_MSD_HEADER_INCLUDED__

#include <vector>
#include <array>
#include <string>
#include <fstream>
#include <cmath>

#include "../IPC.hpp"

class IPCmsd {
public:
    IPCmsd(std::string const& outputFileName, Ensemble const& system, Triad const& pboxSide);
    void accumulateAndPrint(Ensemble const& system);

private:
    IPCmsd();

    std::ofstream meanSquaredDisplFile;
    VectorOfTriads ipcCentersPreviousPositions;
    VectorOfTriads displacementOfEachIPCs;
    int timeCounter;
    Triad boxSide;

    void computeMSD(Ensemble const& system);
    void updatePreviousPositions(Ensemble const& system);
    inline void relativePBC(double & x) {  x -= std::round(x);  }
};

#endif //__IPCPOSTPROCESS_MSD_HEADER_INCLUDED__
