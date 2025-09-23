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

    void propagateWaves();                     // método base (secuencial)
    void propagateWaves(int scheduleType);    // versión paralela (por implementar)
    void propagateWaves(int scheduleType, int chunkSize);
    void propagateWavesCollapse();

    std::vector<Node>& getNodes() { return nodes; }
    const std::vector<Node>& getNodes() const { return nodes; }
    int getSize() const { return networkSize; }
};

#endif // NETWORK_H
