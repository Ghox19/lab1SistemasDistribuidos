#include "WavePropagator.h"
#include <vector>
#include <omp.h>

void WavePropagator::propagateWaves() {
    std::vector<double> newAmplitudes(network.getSize(), 0.0);
    double D = network.getDiffusionCoeff();
    double gamma = network.getDampingCoeff();
    double dt = network.getTimestep();

    for (int i = 0; i < network.getSize(); ++i) {
        double sum_neighbors = 0.0;
        double Ai = network.getNodes()[i].getAmplitude();

        for (int neighborId : network.getNodes()[i].getNeighbors()) {
            sum_neighbors += network.getNodes()[neighborId].getAmplitude() - Ai;
        }

        double diffusion = D * sum_neighbors;
        double damping = -gamma * Ai;
        double delta = diffusion + damping;

        newAmplitudes[i] = Ai + delta * dt;
    }

    for (int i = 0; i < network.getSize(); ++i) {
        network.getNodes()[i].updateAmplitude(newAmplitudes[i]);
    }
}

void WavePropagator::propagateWaves(int scheduleType) {
    // Crear vector temporal para almacenar las nuevas amplitudes calculadas
    // Esto evita condiciones de carrera (race conditions) al separar lectura y escritura
    std::vector<double> newAmplitudes(network.getSize(), 0.0);
    
    // Obtener los parámetros físicos de la simulación desde la red
    double D = network.getDiffusionCoeff();      // Coeficiente de difusión
    double gamma = network.getDampingCoeff();    // Coeficiente de amortiguación  
    double dt = network.getTimestep();           // Paso temporal

    // BRANCH 1: Scheduling estático
    if (scheduleType == 0) {
        // static: Divide las iteraciones en chunks iguales entre threads
        // Cada thread recibe un bloque contiguo de iteraciones
        // Mejor para cargas de trabajo uniformes
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < network.getSize(); ++i) {
            // Variables locales a cada thread (automáticamente privadas)
            double sum_neighbors = 0.0;
            double Ai = network.getNodes()[i].getAmplitude();  // Amplitud actual del nodo i

            // Calcular suma de diferencias con vecinos
            // Esta es la parte de difusión: la onda se propaga hacia vecinos
            for (int neighborId : network.getNodes()[i].getNeighbors()) {
                // (A_vecino - A_i): diferencia que impulsa la difusión
                sum_neighbors += network.getNodes()[neighborId].getAmplitude() - Ai;
            }

            // Aplicar ecuación física de propagación de ondas:
            // A_i(t+dt) = A_i(t) + dt * [D * difusión - γ * amortiguación]
            double diffusion = D * sum_neighbors;    // Término de difusión
            double damping = -gamma * Ai;           // Término de amortiguación (negativo)
            
            // Método de Euler explícito para integración temporal
            newAmplitudes[i] = Ai + (diffusion + damping) * dt;
        }
        
    // BRANCH 2: Scheduling dinámico    
    } else if (scheduleType == 1) {
        // dynamic: Asigna iteraciones a threads conforme terminan trabajo previo
        // Mejor balanceo para cargas de trabajo irregulares
        // Más overhead pero mejor distribución si nodos tienen diferente complejidad
        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < network.getSize(); ++i) {
            // Misma lógica física que en static
            double sum_neighbors = 0.0;
            double Ai = network.getNodes()[i].getAmplitude();

            for (int neighborId : network.getNodes()[i].getNeighbors()) {
                sum_neighbors += network.getNodes()[neighborId].getAmplitude() - Ai;
            }

            double diffusion = D * sum_neighbors;
            double damping = -gamma * Ai;
            newAmplitudes[i] = Ai + (diffusion + damping) * dt;
        }
        
    // BRANCH 3: Scheduling guiado
    } else if (scheduleType == 2) {
        // guided: Empieza con chunks grandes y los reduce exponencialmente
        // Combina ventajas de static (menos overhead) con dynamic (balanceo)
        // Chunks iniciales grandes → eficiencia, chunks finales pequeños → balanceo
        #pragma omp parallel for schedule(guided)
        for (int i = 0; i < network.getSize(); ++i) {
            // Misma lógica física que en las otras ramas
            double sum_neighbors = 0.0;
            double Ai = network.getNodes()[i].getAmplitude();

            for (int neighborId : network.getNodes()[i].getNeighbors()) {
                sum_neighbors += network.getNodes()[neighborId].getAmplitude() - Ai;
            }

            double diffusion = D * sum_neighbors;
            double damping = -gamma * Ai;
            newAmplitudes[i] = Ai + (diffusion + damping) * dt;
        }
    }

    // FASE DE ACTUALIZACIÓN (ejecutada por un solo thread, secuencialmente)
    // Una vez calculadas todas las nuevas amplitudes, actualizamos los nodos
    // Esto se hace fuera de la región paralela para evitar condiciones de carrera
    for (int i = 0; i < network.getSize(); ++i) {
        network.getNodes()[i].updateAmplitude(newAmplitudes[i]);
    }
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
        newAmplitudes[i] = Ai + (diffusion + damping) * dt;
    }

    for (int i = 0; i < network.getSize(); ++i) {
        network.getNodes()[i].updateAmplitude(newAmplitudes[i]);
    }
}

void WavePropagator::propagateWavesCollapse() {
    std::vector<std::vector<double>> grid(100, std::vector<double>(100, 0.0));
    
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < 100; ++i) {
        for (int j = 0; j < 100; ++j) {
            // Simulación en grid 2D con collapse
            grid[i][j] = (i + j) * 0.1;
        }
    }
}

// 2. INTEGRACIÓN TEMPORAL

void WavePropagator::integrateEuler() {
    propagateWaves(); // Método básico
}

void WavePropagator::integrateEuler(int syncType) {
    std::vector<double> newAmplitudes(network.getSize(), 0.0);
    double totalEnergy = 0.0;

    if (syncType == 0) { // atomic
        #pragma omp parallel for
        for (int i = 0; i < network.getSize(); ++i) {
            double Ai = network.getNodes()[i].getAmplitude();
            newAmplitudes[i] = Ai * 0.99; // Simple decay
            
            #pragma omp atomic
            totalEnergy += Ai * Ai;
        }
    } else if (syncType == 1) { // critical
        #pragma omp parallel for
        for (int i = 0; i < network.getSize(); ++i) {
            double Ai = network.getNodes()[i].getAmplitude();
            newAmplitudes[i] = Ai * 0.99;
            
            #pragma omp critical
            {
                totalEnergy += Ai * Ai;
            }
        }
    } else if (syncType == 2) { // nowait
        #pragma omp parallel
        {
            #pragma omp for nowait
            for (int i = 0; i < network.getSize(); ++i) {
                double Ai = network.getNodes()[i].getAmplitude();
                newAmplitudes[i] = Ai * 0.99;
            }
        }
    }

    for (int i = 0; i < network.getSize(); ++i) {
        network.getNodes()[i].updateAmplitude(newAmplitudes[i]);
    }
}

void WavePropagator::integrateEuler(int syncType, bool useBarrier) {
    integrateEuler(syncType);
    
    if (useBarrier) {
        #pragma omp parallel
        {
            #pragma omp barrier
            // Sincronización explícita
        }
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

    if (usePrivate) {
        #pragma omp parallel
        {
            double privateEnergy = 0.0;
            
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

void WavePropagator::simulatePhasesBarrier() {
    #pragma omp parallel
    {
        // Fase 1: Calcular nuevas amplitudes
        #pragma omp for
        for (int i = 0; i < network.getSize(); ++i) {
            double amp = network.getNodes()[i].getAmplitude();
            // Procesamiento fase 1
        }

        #pragma omp barrier

        // Fase 2: Actualizar amplitudes
        #pragma omp for
        for (int i = 0; i < network.getSize(); ++i) {
            double amp = network.getNodes()[i].getAmplitude();
            network.getNodes()[i].updateAmplitude(amp * 0.99);
        }
    }
}

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

void WavePropagator::calculateFinalStateLastprivate() {
    double finalValue = 0.0;

    #pragma omp parallel for lastprivate(finalValue)
    for (int i = 0; i < network.getSize(); ++i) {
        double amp = network.getNodes()[i].getAmplitude();
        finalValue = amp; // El último thread preserva este valor
    }

    // finalValue contiene el valor de la última iteración
}
