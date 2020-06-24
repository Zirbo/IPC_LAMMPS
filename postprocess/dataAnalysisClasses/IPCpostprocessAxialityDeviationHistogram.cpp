#include <fstream>
#include <iostream>
#include <cmath>

#include "IPCpostprocessAxialityDeviationHistogram.hpp"

IPCAxialityDeviationHistogram::IPCAxialityDeviationHistogram() {
    numberOfBins = 1000;
    maxSampling = -0.995;
    totalSamples = 0;
    axialityDeviationHistogram.resize(numberOfBins + 1, 0.);
}

void IPCAxialityDeviationHistogram::accumulate(VectorOfTriads const& firstPatchOrientations, VectorOfTriads const& secndPatchOrientations) {
    const double scaling = -numberOfBins/(1. + maxSampling);

    for (int i = 0; i < (int)firstPatchOrientations.size(); ++i) {
        double cosine = 0.;
        for (int d: DIMENSIONS) {
            cosine += firstPatchOrientations[i][d]*secndPatchOrientations[i][d];
        }

        int cosint = 0;
        if (cosine > maxSampling)
            cosint = numberOfBins;
        else
            cosint = int( (cosine - maxSampling)*scaling );
        axialityDeviationHistogram[cosint] += 1;
    }
    ++totalSamples;
}

void IPCAxialityDeviationHistogram::print(const std::string &outputFileName, int nIPCs) {
    const double norm = 1./(totalSamples*nIPCs);
    const double inorm = -(1. + maxSampling)/numberOfBins;
    std::ofstream outputFile(outputFileName);
    const double degrees = 180./M_PI;

    for(int i = 0; i < numberOfBins; ++i) {
        double place = (inorm*(i+1)) + maxSampling;
        outputFile << place << "\t" << std::acos(place)*degrees << "\t" << norm*axialityDeviationHistogram[i] << "\n";
    }

}
