#ifndef METRICS_CALCULATOR_H
#define METRICS_CALCULATOR_H

#include <Network.h>

class MetricsCalculator {
private:
    Network& network;

public:
    MetricsCalculator(Network& net) : network(net) {};

    double calculateTotalEnergy();
    double calculatePropagationSpeed();              // Nuevo método para limpiar vecinos
    double calculateAverageAmplitude();

};

#endif // NODE_H
