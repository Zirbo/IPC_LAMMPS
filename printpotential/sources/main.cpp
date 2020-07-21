#include "printPotential.hpp"


int main( int argc, char *argv[] ) {
    // parse arguments
    if (argc != 3)
        PotentialForLammps::error("Wrong number of arguments.\n");
    const std::string inputFileName(argv[1]);
    const std::string outputDirName(argv[2]);

    PotentialForLammps potential(inputFileName);
    potential.computeSiteSitePotentials();
    potential.printLAMMPSpotentialsToFile(outputDirName);
}
