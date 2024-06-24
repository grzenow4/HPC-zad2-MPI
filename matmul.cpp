#include <cassert>
#include <iostream>
#include <string>
#include <vector>

#include <mpi.h>

#include "matrix.hh"
#include "grid.hh"
#include "utils.hh"

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int numProcesses, myRank;
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    struct timespec spec;
    clock_gettime(CLOCK_REALTIME, &spec);
    srand(spec.tv_nsec);

    InputOptions options = parseInput(argc, argv);

    if (numProcesses % options.layers != 0) {
        std::cerr << "Layers must divide numProcesses\n";
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    if (!isSquare(numProcesses / options.layers)) {
        std::cerr << "Processes must form a square in each layer\n";
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    {
        Grid grid(numProcesses, myRank, options.layers);
        grid.readMatrices(options.fileA, options.fileB);
        if (options.type == "2D") {
            grid.SUMMA_2D(true);
        } else {
            grid.SUMMA_3D();
        }
        if (options.gValue > 0) {
            grid.printGValue(options.gValue);
        }
        if (options.verbose) {
            grid.printMatrix();
        }
    } // ~Grid() must be called before MPI_Finalize();

    MPI_Finalize();
    return 0;
}
