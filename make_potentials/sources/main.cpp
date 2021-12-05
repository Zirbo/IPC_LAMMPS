#include <unistd.h>

#include <cstring>
#include <iostream>

#include "printPotential.hpp"

static void usage() {
  std::cerr << "Usage:\n"
            << "\t-c start from contact values\n"
            << "\t-e start from epsilons\n"
            << "\t-m model (janus, ipc, asym-ipc)\n"
            << "\t-i inputfile\n"
            << "\t-o output dir\n"
            << "All parameters are required,\n"
            << "except for -c and -e which are mutually exclusive.";
  exit(EXIT_FAILURE);
}

int main(int argc, char* argv[]) {
  IpcType type = IpcType::NONE;
  std::string inputFileName;
  std::string outputDirName;
  bool startFromContactValues = false;
  bool startFromEpsilons = false;

  int opt = 0;
  while ((opt = getopt(argc, argv, "cem:i:o:")) != -1) {
    switch (opt) {
      case 'c':
        startFromContactValues = true;
        break;
      case 'e':
        startFromEpsilons = true;
        break;
      case 'm':
        if (strcmp(optarg, "janus") == 0) {
          type = IpcType::JANUS;
        } else if (strcmp(optarg, "ipc") == 0) {
          type = IpcType::IPC;
        } else if (strcmp(optarg, "asym-ipc") == 0) {
          type = IpcType::ASYM_IPC;
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

  if (inputFileName.empty() || outputDirName.empty() || type == IpcType::NONE) {
    // these parameters are mandatory
    usage();
  }
  if (startFromContactValues == startFromEpsilons) {
    // either both or none given
    usage();
  }

  PotentialForLammps potential(inputFileName, type, startFromContactValues);
  potential.printLAMMPSpotentialsToFile(outputDirName);
  potential.printRadialPotentialsToFile(outputDirName);
  potential.printAngularPotentialsToFile(outputDirName);
}
