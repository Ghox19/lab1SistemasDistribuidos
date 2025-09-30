#ifndef WAVEPROPAGATOR_H
#define WAVEPROPAGATOR_H

#include "Network.h"
#include <omp.h>

class WavePropagator {
private:
    Network& network;

public:
    WavePropagator(Network& net) : network(net) {}

    // 1. Métodos de propagación básicos
    void propagateWaves();
    void propagateWaves(int scheduleType);
    void propagateWaves(int scheduleType, int chunkSize);
    void propagateWavesCollapse();

    // 2. Integración temporal
    void integrateEuler();
    void integrateEuler(int syncType);
    void integrateEuler(int syncType, bool useBarrier);

    // 3. Cálculo de energía
    double calculateEnergy();
    double calculateEnergy(int method);
    double calculateEnergy(int method, bool usePrivate);

    // 4. Procesamiento de nodos
    void processNodes();
    void processNodes(int taskType);
    void processNodes(int taskType, bool useSingle);

    // 5. Métodos específicos para cláusulas únicas
    void simulatePhasesBarrier();
    void parallelInitializationSingle();
    double calculateMetricsFirstprivate();
    void calculateFinalStateLastprivate();
};

#endif // WAVEPROPAGATOR_H
