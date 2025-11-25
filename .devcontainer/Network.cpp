#include "Network.h"
#include <cstdlib>  // rand()
#include <ctime>    // time()
#include <cmath>
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
    if (dimensions > 2 || dimensions < 1) {
        std::cerr << "Solo implementado para una o dos dimensiones" << std::endl;
        return;
    }
    // Caso unidimensional
    if (dimensions == 1) {
        for (int i = 0; i < networkSize; ++i) {
            if (i > 0)
                nodes[i].addNeighbor(i - 1);
            if (i < networkSize - 1)
                nodes[i].addNeighbor(i + 1);
        }
        return;
    }

    // Caso bidimensional
    int rows = static_cast<int>(std::sqrt(this->getSize()));
    while (rows > 1 && this->getSize() % rows != 0) {
        --rows;
    }
    int cols = this->getSize() / rows;


    for (int r = 0; r < rows; ++r) {
        for (int c = 0; c < cols; ++c) {
            int idx = r * cols + c;
            // vecino superior
            if (r > 0) {
                int up = (r - 1) * cols + c;
                nodes[idx].addNeighbor(up);
            }
            // vecino inferior
            if (r < rows - 1) {
                int down = (r + 1) * cols + c;
                nodes[idx].addNeighbor(down);
            }
            // vecino izquierdo
            if (c > 0) {
                int left = r * cols + (c - 1);
                nodes[idx].addNeighbor(left);
            }
            // vecino derecho
            if (c < cols - 1) {
                int right = r * cols + (c + 1);
                nodes[idx].addNeighbor(right);
            }
        }
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
    std::vector<double> newAmplitudes(getSize(), 0.0);
    
    // Obtener los parámetros físicos de la simulación desde la red
    double D = getDiffusionCoeff();      // Coeficiente de difusión
    double gamma = getDampingCoeff();    // Coeficiente de amortiguación  
    double dt = getTimestep();           // Paso temporal

    // BRANCH 1: Scheduling estático
    if (scheduleType == 0) {
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < getSize(); ++i) {
            double sum_neighbors = 0.0;
            double Ai = getNodes()[i].getAmplitude();  // Amplitud actual del nodo i
            for (int neighborId : getNodes()[i].getNeighbors()) {
                sum_neighbors += getNodes()[neighborId].getAmplitude() - Ai;
            }

            // A_i(t+dt) = A_i(t) + dt * [D * difusión - γ * amortiguación]
            double diffusion = D * sum_neighbors;    
            double damping = -gamma * Ai;          
            
            // Método de Euler explícito para integración temporal
            newAmplitudes[i] = Ai + (diffusion + damping) * dt;
        }
        
    // BRANCH 2: Scheduling dinámico    
    } else if (scheduleType == 1) {
        #pragma omp parallel for schedule(dynamic)
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
        
    // BRANCH 3: Scheduling guiado
    } else if (scheduleType == 2) {
        #pragma omp parallel for schedule(guided)
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
    int side = static_cast<int>(std::sqrt(getSize()));
    if (side * side != getSize()) {
        std::cerr << "Error: La red debe ser cuadrada para usar collapse en 2D." << std::endl;
        return;
    }

    std::vector<double> newAmplitudes(getSize(), 0.0);
    double D = getDiffusionCoeff();
    double gamma = getDampingCoeff();
    double dt = getTimestep();

    #pragma omp parallel for collapse(2)
    for (int i = 0; i < side; ++i) {
        for (int j = 0; j < side; ++j) {
            int idx = i * side + j;
            double sum_neighbors = 0.0;
            double Ai = getNodes()[idx].getAmplitude();

            if (i > 0)        sum_neighbors += getNodes()[(i - 1) * side + j].getAmplitude() - Ai;   
            if (i < side - 1) sum_neighbors += getNodes()[(i + 1) * side + j].getAmplitude() - Ai; 
            if (j > 0)        sum_neighbors += getNodes()[i * side + (j - 1)].getAmplitude() - Ai;  
            if (j < side - 1) sum_neighbors += getNodes()[i * side + (j + 1)].getAmplitude() - Ai;   

            double diffusion = D * sum_neighbors;
            double damping = -gamma * Ai;
            newAmplitudes[idx] = Ai + (diffusion + damping) * dt;
        }
    }

    for (int i = 0; i < getSize(); ++i) {
        getNodes()[i].updateAmplitude(newAmplitudes[i]);
    }
}
