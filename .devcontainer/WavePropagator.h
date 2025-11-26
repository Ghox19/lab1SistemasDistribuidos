#ifndef WAVEPROPAGATOR_H
#define WAVEPROPAGATOR_H

#include "Network.h"
#include <omp.h>

class WavePropagator {
private:
    Network& network;

public:
    WavePropagator(Network& net) : network(net) {}

    // 1. Integración temporal
    void integrateEulerCore(std::vector<double>& newAmplitudes, double& totalEnergy, int syncType) ;
    void integrateEuler();
    void integrateEuler(int syncType);
    void integrateEuler(int syncType, bool useBarrier);

    // 2. Cálculo de energía
    double calculateEnergy();
    double calculateEnergy(int method);
    double calculateEnergy(int method, bool usePrivate);

    // 3. Procesamiento de nodos
    void processNodes();
    void processNodes(int taskType);
    void processNodes(int taskType, bool useSingle);

    // 4. Métodos específicos para cláusulas únicas
    void parallelInitializationSingle();
    double calculateMetricsFirstprivate();
};

#endif // WAVEPROPAGATOR_H
