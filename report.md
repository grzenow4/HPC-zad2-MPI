# Author: Grzegorz Nowakowski

## Project structure overview
My solution consists of three (and a half) modules:
1. Utils: This module parses the program input parameters and checks if they are correct.
1. Matrix: This module contains the `Matrix` class that stores a sparse matrix in a `std::vector<MatrixElement>` and implements `add`, `multiply` and `transform` operations on matrices in such sparse representation. `MatrixElement` is a struct that represents a single matrix element `{row, col, val}`.
1. Grid: This module represents the process grid with all operations specified in the task description and in the paper. It is responsible for reading the input matrices in CSR format, distributing them among processes, performing `SUMMA` algorithms on these matrices (with the help of the previous module), and also gathering `gValue` and all other elements of the result matrix to be printed.
1. `matmul.cpp`: Imports all essential functions and uses them in the `main()` function.

## Intensity of the problem

### Communication complexity of SUMMA3D
|                  | A-Bcast                      | B-Bcast                      | AllToall      |
|------------------|:----------------------------:|:----------------------------:|:-------------:|
| Per process data | nnz(A) / p                   | nnz(B) / p                   | FLOPs / p     |
| Comm. size       | sqrt(p / l)                  | sqrt(p / l)                  | l             |
| Latency          | L * log(p / l)               | L * log(p / l)               | L * l         |
| Bandwidth        | b * nnz(A) / p               | b * nnz(B) / p               | b * FLOPs / p |
| How many times   | sqrt(p / l)                  | sqrt(p / l)                  | 1             |
| Total latency    | L * log(p / l) * sqrt(p / l) | L * log(p / l) * sqrt(p / l) | L * l         |
| Total bandwidth  | b * nnz(A) / sqrt(p * l)     | b * nnz(A) / sqrt(p * l)     | b * FLOPs / p |

### Computational complexity of SUMMA3D
|                | LocalMultiply             | Merge Layer            | Merge Fiber        |
|----------------|:-------------------------:|:----------------------:|:------------------:|
| Complexity     | FLOPs / (p * sqrt(p / l)) | FLOPs / p * log(p / l) | FLOPs / p * log(l) |
| How many times | sqrt(p / l)               | 1                      | 1                  |
| Total          | FLOPs / p                 | FLOPs / p * log(p / l) | FLOPs / p * log(l) |

## Description of my solution

### Reading the input
Since the first process may read the whole matrix, it indeed does. After that, it counts the elements that will be sent to every other process and sends that number to the appropriate processes with `MPI_Send`. Then, it asynchronously sends all elements to each process that they belong to using `MPI_Isend`. In the end, it waits for each send to complete.

Other processes first use `MPI_Recv` to receive how many elements they are going to expect, and then receive that many elements with `MPI_Recv`.

### 2D-SUMMA
This algorithm works as described in the paper, using two similar helper functions, `broadcastRow` and `broadcastCol`, so I will focus on the first one. Since the logic in this algorithm is the same as in the paper, I will discuss only the communication model I used.

At first, we split processes in `MPI_COMM_WORLD` to `layers * k` comms, where `k^2` is the number of processes in a layer. Since processes from the same row in the same layer are together in `comm`, the root process (`stage` rank in comm) broadcasts the information about its matrix size with `MPI_Ibcast`. After that, using `MPI_Ibcast` once again, the root process sends its submatrix to all other processes in this row.

Using non-blocking broadcast is the best way to perform these steps of the algorithm since there is very low latency in comparison to making `k - 1` sends and receives. Also, there is no option for a program to deadlock.

### 3D-SUMMA
Firstly, `2D-SUMMA` is invoked, so we skip this step as it is described above.

There are two main parts of this algotihm after `2D-SUMMA` completes:
- `colSplit` splits the result matrix from `2D-SUMMA` into `l` parts, but it does not include any MPI code, so it is a bit boring.
- `allToAll`:
    1. Splits the processes by their ranks inside layers to form fibers.
    1. Using `MPI_Alltoall`, processes share with each other in a fiber how many elements the `i-th` process is going to send.
    1. After all processes have recvCounts, they compute displacements arrays and with the `MPI_Alltoallv` function, they send and receive their submatrices to and from each other. After this operation, every process has all `l` submatrices from its fiber that belong to it.
    1. Having these submatrices, each process locally adds all of them, creating submatrix `C`, which is a part of the whole result matrix.

### Gathering results

#### GValue
Each process computes its local `gValue`. Then, using `MPI_Reduce`, we sum up all values to obtain the global result, which is printed by the process 0.

#### Result matrix
Using `MPI_Gather`, process 0 gathers the number of elements in every result matrix from each process. Then, it computes the displacements array to prepare for the `MPI_Gatherv` function, which gathers all elements into one array. `MPI_Gatherv` is a tricky workaround for using `numProcesses - 1` receives. Even if non-blocking sends/recvs are used, there would be higher latency due to the increased number of communications.

## Performance tests
On matrix with 1,000 rows and columns and 100,000 non-zero elements it takes approximately:
- 0.16s to read the input and distribute elements among processes with 24x24 grid,
- 0.02s to perform 2D-SUMMA algorithm and gather gValue,
- 0.27s to read the input and distribute elements among processes with 12x12x10 grid,
- 0.04s to perform 3D-SUMMA algorithm and gather gValue.

On matrix with 20,000 rows and columns and 1,000,000 non-zero elements it takes approximately:
- 1s to read the input and distribute elements among processes with 24x24 grid,
- 0.55s to perform 2D-SUMMA algorithm and gather gValue,
- 1.12s to read the input and distribute elements among processes with 12x12x10 grid,
- 0.12s to perform 3D-SUMMA algorithm and gather gValue.

On matrix with 100,000 rows and columns and 10,000,000 non-zero elements it takes approximately:
- 12s to read the input and distribute elements among processes with 24x24 grid,
- 8.84s to perform 2D-SUMMA algorithm and gather gValue,
- 12s to read the input and distribute elements among processes with 12x12x10 grid,
- 4.7s to perform 3D-SUMMA algorithm and gather gValue.
