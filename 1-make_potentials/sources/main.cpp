#include <unistd.h>

#include <cstring>
#include <iostream>

#include "printPotential.hpp"

static void
usage()
{
  std::cerr << "LAMMPS potentials generator for Inverse Patchy Colloids\n"
            << "made by S.Ferrari under the supervision of "
            << "Prof. E. Bianchi of the T.U. Wien :)\n\n"
            << "Usage:\n"
            << "\t-c start from contact values\n"
            << "\t-e start from epsilons\n"
            << "\t-s symmetry, allowed values are janus, symm, asymm,\n"
            << "\t\t which stand for Janus, Symmetric and Asymmetric OSPCs\n"
            << "\t-m mapping, can be g (geometric) or e (exponential)\n"
            << "\t-i inputfile\n"
            << "\t-o output dir\n"
            << "All parameters are required, "
            << "except for -c and -e which are mutually exclusive.\n";
  exit(EXIT_FAILURE);
}

int
main(int argc, char* argv[])
{
  Symmetry symmetry = Symmetry::NONE;
  Colloid colloid = Colloid::OSPC;
  Mapping mapping = Mapping::NONE;
  std::string inputFileName;
  std::string outputDirName;
  bool startFromContactValues = false;
  bool startFromEpsilons = false;

  int opt = 0;
  while ((opt = getopt(argc, argv, "hceps:m:i:o:")) != -1) {
    switch (opt) {
      case 'h':
        usage();
        break;
      case 'c':
        startFromContactValues = true;
        break;
      case 'e':
        startFromEpsilons = true;
        break;
      case 'p':
        colloid = Colloid::IPC;
        break;
      case 's':
        if (strcmp(optarg, "janus") == 0) {
          symmetry = Symmetry::JANUS;
        } else if (strcmp(optarg, "symm") == 0) {
          symmetry = Symmetry::SYMMETRIC;
        } else if (strcmp(optarg, "asymm") == 0) {
          symmetry = Symmetry::ASYMMETRIC;
        }
        break;
      case 'm':
        if (strcmp(optarg, "g") == 0) {
          mapping = Mapping::GEOMETRIC;
        } else if (strcmp(optarg, "e") == 0) {
          mapping = Mapping::EXPONENTIAL;
        }
        break;
      case 'i':
        inputFileName = optarg;
        break;
      case 'o':
        outputDirName = optarg;
        break;
      default:
        usage();
    }
  }

  if (inputFileName.empty() || outputDirName.empty()) {
    // these parameters are mandatory
    std::cerr << "ERROR: input or output directory is empty!\n\n";
    usage();
  }
  if (startFromContactValues == startFromEpsilons) {
    // either both or none given, but we only want one!
    std::cerr << "ERROR: only one of -c or -e is allowed!\n\n";
    usage();
  }
  if (symmetry == Symmetry::NONE || mapping == Mapping::NONE) {
    std::cerr << "ERROR: symmetry or mapping missing\n\n";
    usage();
  }

  try {
    PotentialForLammps potential(
      inputFileName, symmetry, colloid, mapping, startFromContactValues);
    potential.printLAMMPSpotentialsToFile(outputDirName);
    potential.printRadialPotentialsToFile(outputDirName);
    potential.printAngularPotentialsToFile(outputDirName);
    potential.printPotentialAlongPathToFile(outputDirName);
    potential.printRecapFile(outputDirName);
  } catch (std::runtime_error& e) {
    std::cerr << "ERROR: " << e.what() << "!\n\n";
    usage();
  }
}
