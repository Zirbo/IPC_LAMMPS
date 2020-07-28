#ifndef __IPCPOSTPROCESS_ORIENTATIONSANALYSIS_SPHERE_HEADER_INCLUDED__
#define __IPCPOSTPROCESS_ORIENTATIONSANALYSIS_SPHERE_HEADER_INCLUDED__

#include <vector>
#include <cmath>

#include "../IPC.hpp"
#include "../IPCpostprocessPotential.hpp"

class IPCorientationsAnalysis3D {
public:
    IPCorientationsAnalysis3D(int nBins)
        : orientationHistogramSize(nBins) {
        side = 2*orientationHistogramSize;
        square = side*side;
        orientationsHistogram.resize(square*side, 0.);
    }
    void accumulate(VectorOfTriads const& ipcOrientations);
    void print(std::string const& outputFileName);

private:
    std::vector<double> orientationsHistogram;
    int orientationHistogramSize;
    int totalSamples;
    int side;
    int square;

    inline void relativePBC(double & x) {  x -= std::round(x);  }
};

#endif //__IPCPOSTPROCESS_ORIENTATIONSANALYSIS_SPHERE_HEADER_INCLUDED__
