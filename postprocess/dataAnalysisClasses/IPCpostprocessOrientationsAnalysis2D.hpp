#ifndef __IPCPOSTPROCESS_ORIENTATIONSANALYSIS_HEADER_INCLUDED__
#define __IPCPOSTPROCESS_ORIENTATIONSANALYSIS_HEADER_INCLUDED__

#include <vector>
//#include <cmath>

#include "../IPC.hpp"
#include "../IPCpostprocessPotential.hpp"

class IPCorientationsAnalysis2D {
public:
    IPCorientationsAnalysis2D()
        : orientationHistogramSize(40) {
        orientationsHistogram.resize(orientationHistogramSize, std::vector<double>(2*orientationHistogramSize, 0.));
    }
    void accumulate(VectorOfTriads const& ipcOrientations);
    void print(std::string const& outputFileName, const int nIPCs);

private:
    std::vector<std::vector<double>> orientationsHistogram;
    int orientationHistogramSize;
    int totalSamples;

    inline void relativePBC(double & x) {  x -= std::round(x);  }
};

#endif //__IPCPOSTPROCESS_ORIENTATIONSANALYSIS_HEADER_INCLUDED__
