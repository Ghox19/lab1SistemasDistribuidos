#ifndef NETWORK_H
#define NETWORK_H

#include <vector>
#include "Node.h"
#include <omp.h>

class Network {
private:
    std::vector<Node> nodes;
    int networkSize;
    double diffusionCoeff;
    double dampingCoeff;
    double noiseCoeff;
    double timestep;
    int getCenter();

public:
    Network(int size, double diffCoeff, double dampCoeff, double noiseCoeff, double dt);
    void initializeRandomNetwork();
    void initializeRegularNetwork(int dimensions);
    double getDiffusionCoeff() const { return diffusionCoeff; }
    double getDampingCoeff() const { return dampingCoeff; }
    double getNoiseCoeff() const { return noiseCoeff; }
    double getTimestep() const { return timestep; }

    std::vector<Node>& getNodes() { return nodes; }
    const std::vector<Node>& getNodes() const { return nodes; }
    int getSize() const { return networkSize; }
    void updateNoiseCoeff(double newNoise) { noiseCoeff = newNoise; }


    void propagateWaves();
    void propagateWaves(int scheduleType);
    void propagateWaves(int scheduleType, int chunkSize);
    void propagateWavesCollapse();
};

#endif // NETWORK_H
