#include "Network.h"
#include <cstdlib>  // rand()
#include <ctime>    // time()
#include <iostream>

Network::Network(int size, double diffCoeff, double dampCoeff, double dt)
    : networkSize(size), diffusionCoeff(diffCoeff), dampingCoeff(dampCoeff), timestep(dt)
{
    for (int i = 0; i < size; ++i)
        nodes.emplace_back(i);
}

void Network::initializeRandomNetwork() {
    std::srand(std::time(nullptr));

    for (int i = 0; i < networkSize; ++i) {
        nodes[i].clearNeighbors();
        int numNeighbors = 1 + std::rand() % (networkSize - 1);
        for (int n = 0; n < numNeighbors; ++n) {
            int neighborId;
            do {
                neighborId = std::rand() % networkSize;
            } while (neighborId == i);
            nodes[i].addNeighbor(neighborId);
        }
    }
}

void Network::initializeRegularNetwork(int dimensions) {
    if (dimensions != 1) {
        std::cerr << "Solo implementado para dimensión 1 (línea)" << std::endl;
        return;
    }

    for (int i = 0; i < networkSize; ++i) {
        if (i > 0)
            nodes[i].addNeighbor(i - 1);
        if (i < networkSize - 1)
            nodes[i].addNeighbor(i + 1);
    }
}


void Network::propagateWaves() {
    std::vector<double> newAmplitudes(getSize(), 0.0);
    double D = getDiffusionCoeff();
    double gamma = getDampingCoeff();
    double dt = getTimestep();

    for (int i = 0; i < getSize(); ++i) {
        double sum_neighbors = 0.0;
        double Ai = getNodes()[i].getAmplitude();

        for (int neighborId : getNodes()[i].getNeighbors()) {
            sum_neighbors += getNodes()[neighborId].getAmplitude() - Ai;
        }

        double diffusion = D * sum_neighbors;
        double damping = -gamma * Ai;
        double delta = diffusion + damping;

        newAmplitudes[i] = Ai + delta * dt;
    }

    for (int i = 0; i < getSize(); ++i) {
        getNodes()[i].updateAmplitude(newAmplitudes[i]);
    }
}

void Network::propagateWaves(int scheduleType) {
    // Crear vector temporal para almacenar las nuevas amplitudes calculadas
    // Esto evita condiciones de carrera (race conditions) al separar lectura y escritura
    std::vector<double> newAmplitudes(getSize(), 0.0);
    
    // Obtener los parámetros físicos de la simulación desde la red
    double D = getDiffusionCoeff();      // Coeficiente de difusión
    double gamma = getDampingCoeff();    // Coeficiente de amortiguación  
    double dt = getTimestep();           // Paso temporal

    // BRANCH 1: Scheduling estático
    if (scheduleType == 0) {
        // static: Divide las iteraciones en chunks iguales entre threads
        // Cada thread recibe un bloque contiguo de iteraciones
        // Mejor para cargas de trabajo uniformes
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < getSize(); ++i) {
            // Variables locales a cada thread (automáticamente privadas)
            double sum_neighbors = 0.0;
            double Ai = getNodes()[i].getAmplitude();  // Amplitud actual del nodo i

            // Calcular suma de diferencias con vecinos
            // Esta es la parte de difusión: la onda se propaga hacia vecinos
            for (int neighborId : getNodes()[i].getNeighbors()) {
                // (A_vecino - A_i): diferencia que impulsa la difusión
                sum_neighbors += getNodes()[neighborId].getAmplitude() - Ai;
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
        for (int i = 0; i < getSize(); ++i) {
            // Misma lógica física que en static
            double sum_neighbors = 0.0;
            double Ai = getNodes()[i].getAmplitude();

            for (int neighborId : getNodes()[i].getNeighbors()) {
                sum_neighbors += getNodes()[neighborId].getAmplitude() - Ai;
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
        for (int i = 0; i < getSize(); ++i) {
            // Misma lógica física que en las otras ramas
            double sum_neighbors = 0.0;
            double Ai = getNodes()[i].getAmplitude();

            for (int neighborId : getNodes()[i].getNeighbors()) {
                sum_neighbors += getNodes()[neighborId].getAmplitude() - Ai;
            }

            double diffusion = D * sum_neighbors;
            double damping = -gamma * Ai;
            newAmplitudes[i] = Ai + (diffusion + damping) * dt;
        }
    }

    // FASE DE ACTUALIZACIÓN (ejecutada por un solo thread, secuencialmente)
    // Una vez calculadas todas las nuevas amplitudes, actualizamos los nodos
    // Esto se hace fuera de la región paralela para evitar condiciones de carrera
    for (int i = 0; i < getSize(); ++i) {
        getNodes()[i].updateAmplitude(newAmplitudes[i]);
    }
}


void Network::propagateWaves(int scheduleType, int chunkSize) {
    std::vector<double> newAmplitudes(getSize(), 0.0);
    double D = getDiffusionCoeff();
    double gamma = getDampingCoeff();
    double dt = getTimestep();

    omp_sched_t scheduleKind;
    switch (scheduleType) {
        case 0: scheduleKind = omp_sched_static; break;
        case 1: scheduleKind = omp_sched_dynamic; break;
        case 2: scheduleKind = omp_sched_guided; break;
        default: scheduleKind = omp_sched_static; break;
    }
    omp_set_schedule(scheduleKind, chunkSize);

    #pragma omp parallel for schedule(runtime)
    for (int i = 0; i < getSize(); ++i) {
        double sum_neighbors = 0.0;
        double Ai = getNodes()[i].getAmplitude();

        for (int neighborId : getNodes()[i].getNeighbors()) {
            sum_neighbors += getNodes()[neighborId].getAmplitude() - Ai;
        }

        double diffusion = D * sum_neighbors;
        double damping = -gamma * Ai;
        newAmplitudes[i] = Ai + (diffusion + damping) * dt;
    }

    for (int i = 0; i < getSize(); ++i) {
        getNodes()[i].updateAmplitude(newAmplitudes[i]);
    }
}

void Network::propagateWavesCollapse() {
    std::vector<std::vector<double>> grid(100, std::vector<double>(100, 0.0));
    
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < 100; ++i) {
        for (int j = 0; j < 100; ++j) {
            // Simulación en grid 2D con collapse
            grid[i][j] = (i + j) * 0.1;
        }
    }
}