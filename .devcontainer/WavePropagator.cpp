#include "WavePropagator.h"
#include <vector>
#include <omp.h>
#include "Network.h"


void WavePropagator::integrateEulerCore(std::vector<double>& newAmplitudes, double& totalEnergy, int syncType) {
    if (syncType == 0) { 
        #pragma omp parallel for
        for (int i = 0; i < network.getSize(); ++i) {
            double Ai = network.getNodes()[i].getAmplitude();
            newAmplitudes[i] = Ai * 0.99;
            #pragma omp atomic
            totalEnergy += Ai * Ai;
        }
    } else if (syncType == 1) { 
        #pragma omp parallel for
        for (int i = 0; i < network.getSize(); ++i) {
            double Ai = network.getNodes()[i].getAmplitude();
            newAmplitudes[i] = Ai * 0.99;
            #pragma omp critical
            {
                totalEnergy += Ai * Ai;
            }
        }
    } else if (syncType == 2) { 
        #pragma omp parallel
        {
            #pragma omp for nowait
            for (int i = 0; i < network.getSize(); ++i) {
                double Ai = network.getNodes()[i].getAmplitude();
                newAmplitudes[i] = Ai * 0.99;
            }
        }
    }
}

void WavePropagator::integrateEuler() {
    network.propagateWaves(); 
}

void WavePropagator::integrateEuler(int syncType) {
    std::vector<double> newAmplitudes(network.getSize(), 0.0);
    double totalEnergy = 0.0;

    integrateEulerCore(newAmplitudes, totalEnergy, syncType);

    for (int i = 0; i < network.getSize(); ++i) {
        network.getNodes()[i].updateAmplitude(newAmplitudes[i]);
    }
}

void WavePropagator::integrateEuler(int syncType, bool useBarrier) {
    std::vector<double> newAmplitudes(network.getSize(), 0.0);
    double totalEnergy = 0.0;

    #pragma omp parallel
    {
        integrateEulerCore(newAmplitudes, totalEnergy, syncType);

        if (useBarrier) {
            #pragma omp barrier
        }
    }

    for (int i = 0; i < network.getSize(); ++i) {
        network.getNodes()[i].updateAmplitude(newAmplitudes[i]);
    }
}


// 3. CÁLCULO DE ENERGÍA
double WavePropagator::calculateEnergy() {
    double totalEnergy = 0.0;
    
    for (int i = 0; i < network.getSize(); ++i) {
        double amp = network.getNodes()[i].getAmplitude();
        totalEnergy += amp * amp;
    }
    
    return totalEnergy;
}

double WavePropagator::calculateEnergy(int method) {
    double totalEnergy = 0.0;

    if (method == 0) { // reduce
        #pragma omp parallel for reduction(+:totalEnergy)
        for (int i = 0; i < network.getSize(); ++i) {
            double amp = network.getNodes()[i].getAmplitude();
            totalEnergy += amp * amp;
        }
    } else if (method == 1) { // atomic
        #pragma omp parallel for
        for (int i = 0; i < network.getSize(); ++i) {
            double amp = network.getNodes()[i].getAmplitude();
            double localEnergy = amp * amp;
            
            #pragma omp atomic
            totalEnergy += localEnergy;
        }
    }

    return totalEnergy;
}

double WavePropagator::calculateEnergy(int method, bool usePrivate) {
    double totalEnergy = 0.0;
    double privateEnergy;

    if (usePrivate) {
        #pragma omp parallel private(privateEnergy)
        {
            privateEnergy = 0.0;
        
            #pragma omp for
            for (int i = 0; i < network.getSize(); ++i) {
                double amp = network.getNodes()[i].getAmplitude();
                privateEnergy += amp * amp;
            }
            
            #pragma omp atomic
            totalEnergy += privateEnergy;
        }
    } else {
        return calculateEnergy(method);
    }

    return totalEnergy;
}

// 4. PROCESAMIENTO DE NODOS
void WavePropagator::processNodes() {
    for (int i = 0; i < network.getSize(); ++i) {
        double amp = network.getNodes()[i].getAmplitude();
        network.getNodes()[i].updateAmplitude(amp * 0.99);
    }
}

void WavePropagator::processNodes(int taskType) {
    if (taskType == 0) { // task
        #pragma omp parallel
        {
            #pragma omp single
            {
                for (int i = 0; i < network.getSize(); ++i) {
                    #pragma omp task
                    {
                        double amp = network.getNodes()[i].getAmplitude();
                        network.getNodes()[i].updateAmplitude(amp * 0.99);
                    }
                }
            }
        }
    } else if (taskType == 1) { // parallel for
        #pragma omp parallel for
        for (int i = 0; i < network.getSize(); ++i) {
            double amp = network.getNodes()[i].getAmplitude();
            network.getNodes()[i].updateAmplitude(amp * 0.99);
        }
    }
}

void WavePropagator::processNodes(int taskType, bool useSingle) {
    if (useSingle) {
        #pragma omp parallel
        {
            #pragma omp single
            {
                processNodes(taskType);
            }
        }
    } else {
        processNodes(taskType);
    }
}

// 5. MÉTODOS ESPECÍFICOS PARA CLÁUSULAS ÚNICAS


void WavePropagator::parallelInitializationSingle() {
    #pragma omp parallel
    {
        #pragma omp single
        {
            // Inicialización que solo debe hacer un thread
            for (int i = 0; i < network.getSize(); ++i) {
                network.getNodes()[i].updateAmplitude(1.0);
            }
        }
    }
}

double WavePropagator::calculateMetricsFirstprivate() {
    double initialValue = 1.0;
    double totalMetric = 0.0;

    #pragma omp parallel for firstprivate(initialValue) reduction(+:totalMetric)
    for (int i = 0; i < network.getSize(); ++i) {
        double amp = network.getNodes()[i].getAmplitude();
        double metric = amp * initialValue;
        totalMetric += metric;
        initialValue *= 0.99; // Modificación local
    }

    return totalMetric;
}


