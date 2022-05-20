#include <iostream>
#include <string>
#include <sstream>
#include "IPCpostprocess.hpp"


int main ( int argc, char *argv[] ) {

    std::stringstream helpMessage;
    helpMessage << "Mandatory arguments:\n - trajectory file\n - inputfile\n - directory were the potentials are stored\n\n"
                << "Optional arguments:\n - time step where to start from\n - time step where to end"
                << "   (if not given, the program goes from start to finish)\n";

    IPCpostprocessInitializer initializer;
    if(argc == 4 || argc == 6) {
        initializer.trajFilename = argv[1];
        initializer.inputFilename = argv[2];
        initializer.potDirName = argv[3];
        if (argc == 6) {
            try {
                initializer.startStep = std::stoi(argv[4]);
                initializer.finalStep = std::stoi(argv[5]);
            } catch(...) {
                std::cout << "Error in parsing the start and/or final time step\n";
                exit(1);
            }
        }
    } else {
        std::cout << helpMessage.str();
        exit(1);
    }

    IPCpostprocess postProcess(initializer);
    postProcess.run();
}
