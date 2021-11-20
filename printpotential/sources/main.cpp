#include <unistd.h>
#include <iostream>
#include <cstring>

#include "printPotential.hpp"

static void usage() {
    std::cerr << "Usage:\n\t-j for janus instead of normal IPCs\n"
        << "\t-i input file\n\t-o output dir\n";
}

int main( int argc, char *argv[] ) {
    IpcType type = IpcType::IPC;
    std::string inputFileName;
    std::string outputDirName;

    int opt = 0;
    while ((opt = getopt(argc, argv, "ji:o:")) != -1) {
        switch(opt) {
        case 'j':
            type = IpcType::JANUS;
            break;
        case 'i':
            inputFileName = optarg;
            break;
        case 'o':
            outputDirName = optarg;
            break;
        default:
            usage();
            return EXIT_FAILURE;
        }
    }

    if (inputFileName.empty() || outputDirName.empty()) {
        usage();
        return EXIT_FAILURE;
    }

    PotentialForLammps potential(inputFileName, type);
    potential.printLAMMPSpotentialsToFile(outputDirName);
    potential.printRadialPotentialsToFile(outputDirName);
    potential.printAngularPotentialsToFile(outputDirName);
}
