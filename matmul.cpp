#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <mpi.h>

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
    printOptions(options);

    MPI_Finalize();
    return 0;
}
