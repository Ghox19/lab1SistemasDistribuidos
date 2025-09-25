#include "Network.h"
#include <cstdlib>  // rand()
#include <ctime>    // time()
#include <iostream>

Network::Network(int size, double diffCoeff, double dampCoeff, double dt)
    : networkSize(size), diffusionCoeff(diffCoeff), dampingCoeff(dampCoeff), timestep(dt)
{
    for (int i = 0; i < size; ++i)
        nodes.emplace_back(i);
}

void Network::initializeRandomNetwork() {
    std::srand(std::time(nullptr));

    for (int i = 0; i < networkSize; ++i) {
        nodes[i].clearNeighbors();
        int numNeighbors = 1 + std::rand() % (networkSize - 1);
        for (int n = 0; n < numNeighbors; ++n) {
            int neighborId;
            do {
                neighborId = std::rand() % networkSize;
            } while (neighborId == i);
            nodes[i].addNeighbor(neighborId);
        }
    }
}

void Network::initializeRegularNetwork(int dimensions) {
    if (dimensions != 1) {
        std::cerr << "Solo implementado para dimensión 1 (línea)" << std::endl;
        return;
    }

    for (int i = 0; i < networkSize; ++i) {
        if (i > 0)
            nodes[i].addNeighbor(i - 1);
        if (i < networkSize - 1)
            nodes[i].addNeighbor(i + 1);
    }
}

