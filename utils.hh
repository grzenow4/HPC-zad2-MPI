#pragma once

#include <string>

struct InputOptions {
    std::string fileA;
    std::string fileB;
    bool verbose = false;
    std::string type;
    int layers = 1;
    double gValue = -1;
};

void printUsage(const std::string& progName);

InputOptions parseInput(int argc, char *argv[]);

void printOptions(InputOptions options);

bool isSquare(int n);
