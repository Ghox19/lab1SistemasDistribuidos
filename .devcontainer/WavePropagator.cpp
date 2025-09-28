#include "WavePropagator.h"
#include <omp.h>
#include <vector>

void WavePropagator::propagateWaves() {
    std::vector<double> newAmplitudes(network.getSize(), 0.0);
    double D = network.getDiffusionCoeff();
    double gamma = network.getDampingCoeff();
    double dt = network.getTimestep();

    // Para cada nodo calcular la ecuación
    for (int i = 0; i < network.getSize(); ++i) {
        double sum_neighbors = 0.0;
        double Ai = network.getNodes()[i].getAmplitude();

        // Sumar diferencias con vecinos
        for (int neighborId : network.getNodes()[i].getNeighbors()) {
            sum_neighbors += network.getNodes()[neighborId].getAmplitude() - Ai;
        }

        double diffusion = D * sum_neighbors;
        double damping = -gamma * Ai;
        double source = 0.0; // No hay fuente externa en esta versión
        double delta = diffusion + damping + source;

        newAmplitudes[i] = Ai + delta * dt;
    }

    // Actualizar amplitudes después del cálculo para todos los nodos
    for (int i = 0; i < network.getSize(); ++i) {
        network.getNodes()[i].updateAmplitude(newAmplitudes[i]);
    }
}

void WavePropagator::propagateWaves(int scheduleType) {
    // Por implementar
}

void WavePropagator::propagateWaves(int scheduleType, int chunkSize) {
    std::vector<double> newAmplitudes(network.getSize(), 0.0);
    double D = network.getDiffusionCoeff();
    double gamma = network.getDampingCoeff();
    double dt = network.getTimestep();

    omp_sched_t scheduleKind;
    switch (scheduleType) {
        case 0: scheduleKind = omp_sched_static; break;
        case 1: scheduleKind = omp_sched_dynamic; break;
        case 2: scheduleKind = omp_sched_guided; break;
        case 3: scheduleKind = omp_sched_auto; break;
        default: scheduleKind = omp_sched_static; break;
    }
    omp_set_schedule(scheduleKind, chunkSize);

    #pragma omp parallel for schedule(runtime)
    for (int i = 0; i < network.getSize(); ++i) {
        double sum_neighbors = 0.0;
        double Ai = network.getNodes()[i].getAmplitude();

        for (int neighborId : network.getNodes()[i].getNeighbors()) {
            sum_neighbors += network.getNodes()[neighborId].getAmplitude() - Ai;
        }

        double diffusion = D * sum_neighbors;
        double damping = -gamma * Ai;
        double source = 0.0; // Sin fuente externa en esta versión
        double delta = diffusion + damping + source;

        newAmplitudes[i] = Ai + delta * dt;
    }

    // Actualizar amplitudes después del cálculo para todos los nodos (no paralelizado para evitar race condition)
    for (int i = 0; i < network.getSize(); ++i) {
        network.getNodes()[i].updateAmplitude(newAmplitudes[i]);
    }
}

void WavePropagator::propagateWavesCollapse() {
    // Por implementar
}
