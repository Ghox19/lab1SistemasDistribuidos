#include "MetricsCalculator.h"
#include <omp.h>

// Cálculo de energía total del sistema
double MetricsCalculator::calculateTotalEnergy() {
    double sum = 0.0;
    std::vector<Node> nodes = network.getNodes();
    int n = network.getSize();

    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        Node node = nodes[i];
        double amplitude = node.getAmplitude(); 
        sum += amplitude * amplitude;
    }

    return sum;
};


// Cálculo de amplitud promedio del sistema
double MetricsCalculator::calculateAverageAmplitude() {
    double sum = 0.0;
    std::vector<Node> nodes = network.getNodes();
    int n = network.getSize();

    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        Node node = nodes[i];
        sum += node.getAmplitude(); 
    }

    return sum / n;
}