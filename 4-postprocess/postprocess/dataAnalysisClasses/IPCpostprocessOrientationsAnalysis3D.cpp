#include <cmath>
#include <fstream>
#include <iomanip>

#include "IPCpostprocessOrientationsAnalysis3D.hpp"

void IPCorientationsAnalysis3D::accumulate(VectorOfTriads const& ipcOrientations) {

    for (auto const& ipcOrientation: ipcOrientations) {
        int x = orientationHistogramSize + int(ipcOrientation[0]*orientationHistogramSize);
        int y = orientationHistogramSize + int(ipcOrientation[1]*orientationHistogramSize);
        int z = orientationHistogramSize + int(ipcOrientation[2]*orientationHistogramSize);
        ++orientationsHistogram[x + y*side + z*square];
    }

    ++totalSamples;
}

void IPCorientationsAnalysis3D::print(std::string const& outputFileName, int const nIPCs) {
    std::ofstream outputFile(outputFileName);
    outputFile << std::scientific << std::setprecision(6);

    const double binSize = 1./orientationHistogramSize;
    double norm = 1./(totalSamples*nIPCs);

    for (int x = 0; x < 2*orientationHistogramSize; ++x) {
        for (int y = 0; y < 2*orientationHistogramSize; ++y) {
            for (int z = 0; z < 2*orientationHistogramSize; ++z) {
                int output = orientationsHistogram[x + y*side + z*square];
                if (output == 0)
                    continue;
                outputFile << (x - orientationHistogramSize)*binSize << "\t"
                           << (y - orientationHistogramSize)*binSize << "\t"
                           << (z - orientationHistogramSize)*binSize << "\t"
                           << norm*output << "\t"
                           << std::log10(norm*output) << "\n";
            }
        }
    }
    outputFile.close();
}
