#include "WavePropagator.h"
#include <vector>
#include <omp.h>
#include "Network.h"

// INTEGRACIÓN TEMPORAL
// Funcion base para integración temporal de Euler
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

// Integración temporal de Euler basica ya implementada en Network
void WavePropagator::integrateEuler() {
    network.propagateWaves(); 
}

// Sobrecarga de integración con tipos de sincronización atomic / critical / nowait
void WavePropagator::integrateEuler(int syncType) {
    std::vector<double> newAmplitudes(network.getSize(), 0.0);
    double totalEnergy = 0.0;

    integrateEulerCore(newAmplitudes, totalEnergy, syncType);

    for (int i = 0; i < network.getSize(); ++i) {
        network.getNodes()[i].updateAmplitude(newAmplitudes[i]);
    }
}

// Sobrecarga de integración con opción de barrera
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


// CALCULO DE ENERGÍA
// Cálculo básico de energía de manera secuencial
double WavePropagator::calculateEnergy() {
    double totalEnergy = 0.0;
    
    for (int i = 0; i < network.getSize(); ++i) {
        double amp = network.getNodes()[i].getAmplitude();
        totalEnergy += amp * amp;
    }
    
    return totalEnergy;
}

// Sobrecarga de cálculo de energía con métodos paralelos, siendo estos atomic y critical
double WavePropagator::calculateEnergy(int method) {
    double totalEnergy = 0.0;

    if (method == 0) { 
        #pragma omp parallel for reduction(+:totalEnergy)
        for (int i = 0; i < network.getSize(); ++i) {
            double amp = network.getNodes()[i].getAmplitude();
            totalEnergy += amp * amp;
        }
    } else if (method == 1) {
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

// Sobrecarga de cálculo de energía con opción de variable privada utilizando atomic
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

// PROCESAMIENTO DE NODOS
// Procesamiento básico de nodos de manera secuencial
void WavePropagator::processNodes() {
    for (int i = 0; i < network.getSize(); ++i) {
        double amp = network.getNodes()[i].getAmplitude();
        network.getNodes()[i].updateAmplitude(amp * 0.99); // Ejemplo simple de amortiguacion
    }
}

// Sobrecarga de procesamiento de nodos con tipos de tareas, utilizando tasking y parallel for
void WavePropagator::processNodes(int taskType) {
    if (taskType == 0) { 
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
    } else if (taskType == 1) { 
        #pragma omp parallel for
        for (int i = 0; i < network.getSize(); ++i) {
            double amp = network.getNodes()[i].getAmplitude();
            network.getNodes()[i].updateAmplitude(amp * 0.99);
        }
    }
}

// Sobrecarga de procesamiento de nodos con opción de cláusula single
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

// MÉTODOS ESPECÍFICOS PARA CLÁUSULAS ÚNICAS
// Simulación de fases con barrera
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

// Sobrecarga de cálculo de métricas con cláusula firstprivate
double WavePropagator::calculateMetricsFirstprivate() {
    double initialValue = 1.0;
    double totalMetric = 0.0;

    #pragma omp parallel for firstprivate(initialValue) reduction(+:totalMetric)
    for (int i = 0; i < network.getSize(); ++i) {
        double amp = network.getNodes()[i].getAmplitude();
        double metric = amp * initialValue;
        totalMetric += metric;
        initialValue *= 0.99; 
    }

    return totalMetric;
}


