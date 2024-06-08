#include <cassert>
#include <string>
#include <vector>

#include <mpi.h>

#include "matrix.hh"
#include "grid.hh"
#include "utils.hh"

#include <random>

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int numProcesses, myRank;
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    if (!isSquare(numProcesses)) {
        std::cerr << "Error wrong number of processes\n";
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    struct timespec spec;
    clock_gettime(CLOCK_REALTIME, &spec);
    srand(spec.tv_nsec);

    InputOptions options = parseInput(argc, argv);

    {
        Grid grid(numProcesses, myRank);
        grid.readMatrices(options.fileA, options.fileB);
        grid.SUMMA_2D();
        grid.printMatrix();
    } // ~Grid() must be called before MPI_Finalize();

    MPI_Finalize();
    return 0;
}
