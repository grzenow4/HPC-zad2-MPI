#include <cmath>
#include <iostream>
#include <vector>

#include <mpi.h>

#include "utils.hh"

void printUsage(const std::string& progName) {
    std::cerr << "Usage: " << progName << " [-a sparse_matrix_file_a] [-b sparse_matrix_file_b] [-v] [-g g_value] [-t 2D|3D|balanced] [-l value]\n";
}

InputOptions parseInput(int argc, char *argv[]) {
    InputOptions options;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];

        if (arg == "-a") {
            if (i + 1 < argc) {
                options.fileA = argv[++i];
            } else {
                printUsage(argv[0]);
                MPI_Finalize();
                exit(EXIT_FAILURE);
            }
        } else if (arg == "-b") {
            if (i + 1 < argc) {
                options.fileB = argv[++i];
            } else {
                printUsage(argv[0]);
                MPI_Finalize();
                exit(EXIT_FAILURE);
            }
        } else if (arg == "-v") {
            options.verbose = true;
        } else if (arg == "-t") {
            if (i + 1 < argc) {
                options.type = argv[++i];
                if (options.type != "2D" && options.type != "3D") {
                    std::cerr << "Unsupported option " << options.type << "\n";
                    MPI_Finalize();
                    exit(EXIT_FAILURE);
                }
            } else {
                printUsage(argv[0]);
                MPI_Finalize();
                exit(EXIT_FAILURE);
            }
        } else if (arg == "-l") {
            if (i + 1 < argc) {
                options.layers = std::stoi(argv[++i]);
            } else {
                printUsage(argv[0]);
                MPI_Finalize();
                exit(EXIT_FAILURE);
            }
        } else if (arg == "-g") {
            if (i + 1 < argc) {
                options.gValue = std::stod(argv[++i]);
            } else {
                printUsage(argv[0]);
                MPI_Finalize();
                exit(EXIT_FAILURE);
            }
        } else {
            printUsage(argv[0]);
            MPI_Finalize();
            exit(EXIT_FAILURE);
        }
    }

    return options;
}

void printOptions(InputOptions options) {
    std::cout << "Options:\n"
              << "fileA: " << options.fileA << "\n"
              << "fileB: " << options.fileB << "\n"
              << "verbose: " << (options.verbose ? "true" : "false") << "\n"
              << "type: " << options.type << "\n"
              << "layers: " << options.layers << "\n"
              << "gValue: " << options.gValue << "\n";
}

bool isSquare(int n) {
    if (n < 0) {
        return false;
    }

    int sr = std::sqrt(n);
    return (sr * sr == n);
}
