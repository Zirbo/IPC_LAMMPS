#include <cmath>
#include <fstream>
#include <iomanip>

#include "IPCpostprocessOrientationsAnalysis2D.hpp"

void IPCorientationsAnalysis2D::accumulate(VectorOfTriads const& ipcOrientations) {
    // polar -> angle with the z-axis; azimuth -> angle with the x-axis
    // accumulute in a histogram the polar and azimuth angles of each IPC
    const double binSize = M_PI/orientationHistogramSize;

    for (auto const& ipcOrientation: ipcOrientations) {
        double polarAngle = std::acos(ipcOrientation[2]);
        const double polarScaling = 1./std::sin(polarAngle);
        double azimuthAngle = M_PI + std::atan2(ipcOrientation[1]*polarScaling, ipcOrientation[0]*polarScaling);

        // we assume patch symmetry, so if cosPolar < 0 we do polar += pi/2 and azimuth += pi
/*        if (ipcOrientation[2] < 0) {
            polarAngle = M_PI - polarAngle;
            azimuthAngle += M_PI;
            if (azimuthAngle > 2*M_PI)
                azimuthAngle -= 2*M_PI;
        }*/

        const int polarBin = int(polarAngle/binSize);
        const int azimuthBin = int(azimuthAngle/binSize);
        ++orientationsHistogram[polarBin][azimuthBin];
    }

    ++totalSamples;
}

void IPCorientationsAnalysis2D::print(std::string const& outputFileName, const int nIPCs) {
    std::ofstream outputFile(outputFileName);
    outputFile << std::scientific << std::setprecision(6);

    double norm = 1./(totalSamples*nIPCs);
    const double binSize = M_PI/orientationHistogramSize;

    for (int a = 0; a < 2*orientationHistogramSize; ++a) {
        for (int p = 0; p < orientationHistogramSize; ++p) {
            outputFile << a*binSize << "\t" << p*binSize << "\t" << orientationsHistogram[p][a]*norm << "\t" << std::log10(orientationsHistogram[p][a]*norm) << "\n";
        }
    }
    outputFile.close();
}
