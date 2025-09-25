#ifndef NETWORK_H
#define NETWORK_H

#include <vector>
#include "Node.h"

class Network {
private:
    std::vector<Node> nodes;
    int networkSize;
    double diffusionCoeff;
    double dampingCoeff;
    double timestep;

public:
    Network(int size, double diffCoeff, double dampCoeff, double dt);

    void initializeRandomNetwork();
    void initializeRegularNetwork(int dimensions);

    // Getters para los parámetros físicos y de simulación
    double getDiffusionCoeff() const { return diffusionCoeff; }
    double getDampingCoeff() const { return dampingCoeff; }
    double getTimestep() const { return timestep; }

    std::vector<Node>& getNodes() { return nodes; }
    const std::vector<Node>& getNodes() const { return nodes; }
    int getSize() const { return networkSize; }
};

#endif // NETWORK_H
