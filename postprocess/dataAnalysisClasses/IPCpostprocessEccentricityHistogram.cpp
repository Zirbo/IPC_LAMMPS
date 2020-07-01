#include <fstream>
#include <cmath>
#include "IPCpostprocessEccentricityHistogram.hpp"

IPCeccentricityHistogram::IPCeccentricityHistogram(double p_fixedEccentricity) {
    semiamplitude = 200;
    percentage = 0.20;
    eccentricityHistogram.resize(2*semiamplitude + 1);
    fixedEccentricity = p_fixedEccentricity;
    inverseSemiamplitude = semiamplitude/(fixedEccentricity*percentage);
    totalSamples = 0;
}

void IPCeccentricityHistogram::accumulate(std::vector<double> const& eccentricities) {
    for(double ecc: eccentricities) {
        int pos = std::round( (ecc - fixedEccentricity)*inverseSemiamplitude ) + semiamplitude;
        if (pos < 0)    pos = 0;
        if (pos > 2*semiamplitude)   pos = 2*semiamplitude + 1;
        ++eccentricityHistogram[pos];
    }
    ++totalSamples;
}

void IPCeccentricityHistogram::print(const std::string &outputFileName, int nIPCs) {
    const double norm = 1./(totalSamples*nIPCs);
    const double inorm = 1./inverseSemiamplitude;
    std::ofstream outputFile(outputFileName);

    for(int i = 0; i < (int)eccentricityHistogram.size(); ++i) {
        outputFile << fixedEccentricity + ( i - semiamplitude )*inorm
                   << "\t" << eccentricityHistogram[i]*norm << "\n";
    }

}
