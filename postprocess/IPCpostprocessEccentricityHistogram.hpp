#ifndef __IPCPOSTPROCESS_ECCENTRICITYHISTOGRAM_HEADER_INCLUDED__
#define __IPCPOSTPROCESS_ECCENTRICITYHISTOGRAM_HEADER_INCLUDED__

#include <vector>
#include <string>

class IPCeccentricityHistogram {
public:
    IPCeccentricityHistogram(double p_fixedEccentricity);
    void accumulate(const std::vector<double> &eccentricities);
    void print(std::string const& outputFileName, int nIPCs);

private:
    IPCeccentricityHistogram();

    std::vector<long> eccentricityHistogram;
    double fixedEccentricity;
    int semiamplitude;
    double percentage;
    double inverseSemiamplitude;
    int totalSamples;
};

#endif //__IPCPOSTPROCESS_ECCENTRICITYHISTOGRAM_HEADER_INCLUDED__
