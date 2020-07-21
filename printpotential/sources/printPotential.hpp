#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>

class PotentialForLammps {
public:
    PotentialForLammps(std::string const& inputFileName);
    void computeSiteSitePotentials();
    void printLAMMPSpotentialsToFile(std::string const& outputDirName);
    void printPotentialsToFileForVisualization(std::string const& outputDirName);

    static void error(std::string const & errorMessage) {
    std::cerr << errorMessage;
    exit(1);
}


private:
    double e_BB, e_Bs1, e_Bs2, e_s1s1, e_s1s2, e_s2s2, e_min;
    double firstPatchEccentricity, firstPatchRadius;
    double secndPatchEccentricity, secndPatchRadius;
    double ipcRadius, interactionRange;
    double fakeHSdiameter, fakeHScoefficient, fakeHSexponent;

    double samplingStep, cutoffValue;

    std::vector<double> uHS, uBB, uBs1, uBs2, us1s2, us1s1, us2s2;
    std::vector<double> fHS, fBB, fBs1, fBs2, fs1s2, fs1s1, fs2s2;

    double computeOmega(double Ra, double Rb, double rab);
    double computeOmegaRadialDerivative(double Ra, double Rb, double rab);
    void printPotentialsToFileForVisualizationSingleOrientation(const int potentialPrintingStep);
};
