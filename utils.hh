#pragma once

#include <string>

struct InputOptions {
    std::string fileA;
    std::string fileB;
    bool verbose;
    std::string type;
    int layers;
    int gValue;
};

void printUsage(const std::string& progName);

InputOptions parseInput(int argc, char *argv[]);

void printOptions(InputOptions options);

bool isSquare(int n);