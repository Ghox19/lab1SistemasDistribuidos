#ifndef NETWORK_H
#define NETWORK_H

#include <vector>      // Para std::vector
#include "Node.h"      // Para la clase Node

class Network {
private:
    std::vector<Node> nodes;
    int networkSize;
    double diffusionCoeff;
    double dampingCoeff;
public:
    Network(int size, double diffCoeff, double dampCoeff)
        : networkSize(size), diffusionCoeff(diffCoeff), dampingCoeff(dampCoeff)
    {
        for (int i = 0; i < size; ++i)
            nodes.emplace_back(i);
    }

    void initializeRandomNetwork();
    void initializeRegularNetwork(int dimensions);

    std::vector<Node>& getNodes() { return nodes; }
    const std::vector<Node>& getNodes() const { return nodes; }
    int getSize() const { return networkSize; }
};

#endif // NETWORK_H