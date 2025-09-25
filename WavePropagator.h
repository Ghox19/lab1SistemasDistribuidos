#ifndef WAVEPROPAGATOR_H
#define WAVEPROPAGATOR_H

#include "Network.h"

class WavePropagator {
private:
    Network& network;

public:
    WavePropagator(Network& net) : network(net) {}

    // Función base secuencial de propagación
    void propagateWaves();

    // Versiones paralelas con scheduling (por implementar)
    void propagateWaves(int scheduleType);
    void propagateWaves(int scheduleType, int chunkSize);
    void propagateWavesCollapse();
};

#endif // WAVEPROPAGATOR_H
