#ifndef __IPCPOSTPROCESS_AXIALITYDEVIATION_HEADER_INCLUDED__
#define __IPCPOSTPROCESS_AXIALITYDEVIATION_HEADER_INCLUDED__

#include <vector>
#include <string>
#include "../IPC.hpp"

class IPCAxialityDeviationHistogram {
public:
    IPCAxialityDeviationHistogram();
    void accumulate(VectorOfTriads const& firstPatchOrientations, VectorOfTriads const& secndPatchOrientations);
    void print(std::string const& outputFileName, int nIPCs);

private:

    std::vector<double> axialityDeviationHistogram;
    int numberOfBins;
    double maxSampling;
    int totalSamples;
};

#endif //__IPCPOSTPROCESS_AXIALITYDEVIATION_HEADER_INCLUDED__
