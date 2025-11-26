#include "Network.h"
#include <cstdlib>  // rand()
#include <ctime>    // time()
#include <cmath>
#include <iostream>

// Constructor de la red
Network::Network(int size, double diffCoeff, double dampCoeff,double noiseCoeff, double dt)
    : networkSize(size), diffusionCoeff(diffCoeff), dampingCoeff(dampCoeff), noiseCoeff(noiseCoeff), timestep(dt)
{
    for (int i = 0; i < size; ++i)
        nodes.emplace_back(i);
}
    

// Inicialización de la red aleatoria
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

// Inicialización de la red regular (1D o 2D)
void Network::initializeRegularNetwork(int dimensions) {
    if (dimensions > 2 || dimensions < 1) {
        std::cerr << "Solo implementado para una o dos dimensiones" << std::endl;
        return;
    }
    // Caso unidimensional
    if (dimensions == 1) {
        for (int i = 0; i < networkSize; ++i) {
            nodes[i].updateAmplitude(0.0);
            if (i > 0)
                nodes[i].addNeighbor(i - 1);
            if (i < networkSize - 1)
                nodes[i].addNeighbor(i + 1);
        }
        for (int i = 0; i < networkSize; ++i) nodes[i].updateAmplitude(0.0);
        return;
    }

    // Caso bidimensional
    int rows = static_cast<int>(std::sqrt(this->getSize()));
    while (rows > 1 && this->getSize() % rows != 0) {
        --rows;
    }
    int cols = this->getSize() / rows;

    // Se definen todos los ve1cinos en una malla 2D
    for (int r = 0; r < rows; ++r) {
        for (int c = 0; c < cols; ++c) {
            int idx = r * cols + c;
            nodes[idx].updateAmplitude(0.0);
            if (r > 0) {
                int up = (r - 1) * cols + c;
                nodes[idx].addNeighbor(up);
            }
            if (r < rows - 1) {
                int down = (r + 1) * cols + c;
                nodes[idx].addNeighbor(down);
            }
            if (c > 0) {
                int left = r * cols + (c - 1);
                nodes[idx].addNeighbor(left);
            }
            if (c < cols - 1) {
                int right = r * cols + (c + 1);
                nodes[idx].addNeighbor(right);
            }
        }
    }
    for (int i = 0; i < networkSize; ++i) nodes[i].updateAmplitude(0.0);
}


#include <cmath>

// Comun para todas: cálculo del centro
double rowsize_d = std::sqrt(static_cast<double>(getSize()));
int rowsize = static_cast<int>(std::round(rowsize_d));
int center = static_cast<int>(std::round((rowsize / 2.0) * rowsize + (rowsize / 2.0)));

// Versión secuencial
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

        if (i == center) newAmplitudes[i] += noiseCoeff;
    }
    for (int i = 0; i < getSize(); ++i)
        getNodes()[i].updateAmplitude(newAmplitudes[i]);
}

// Paralelo con scheduleType
void Network::propagateWaves(int scheduleType) {
    std::vector<double> newAmplitudes(getSize(), 0.0);
    double D = getDiffusionCoeff();
    double gamma = getDampingCoeff();
    double dt = getTimestep();

    double rowsize_d = std::sqrt(static_cast<double>(getSize()));
    int rowsize = static_cast<int>(std::round(rowsize_d));
    int center = static_cast<int>(std::round((rowsize / 2.0) * rowsize + (rowsize / 2.0)));

    if (scheduleType == 0) {
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < getSize(); ++i) {
            double sum_neighbors = 0.0, Ai = getNodes()[i].getAmplitude();
            for (int neighborId : getNodes()[i].getNeighbors())
                sum_neighbors += getNodes()[neighborId].getAmplitude() - Ai;
            double diffusion = D * sum_neighbors, damping = -gamma * Ai, delta = diffusion + damping;
            newAmplitudes[i] = Ai + delta * dt;
        }
    } else if (scheduleType == 1) {
        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < getSize(); ++i) { /* igual que arriba */ }
    } else if (scheduleType == 2) {
        #pragma omp parallel for schedule(guided)
        for (int i = 0; i < getSize(); ++i) { /* igual que arriba */ }
    }
    // Solo el centro: suma la fuente periódica (fuera del for paralelizado)
    newAmplitudes[center] += noiseCoeff;
    for (int i = 0; i < getSize(); ++i)
        getNodes()[i].updateAmplitude(newAmplitudes[i]);
}

// Paralelo por chunks
void Network::propagateWaves(int scheduleType, int chunkSize) {
    std::vector<double> newAmplitudes(getSize(), 0.0);
    double D = getDiffusionCoeff();
    double gamma = getDampingCoeff();
    double dt = getTimestep();

    double rowsize_d = std::sqrt(static_cast<double>(getSize()));
    int rowsize = static_cast<int>(std::round(rowsize_d));
    int center = static_cast<int>(std::round((rowsize / 2.0) * rowsize + (rowsize / 2.0)));

    omp_sched_t scheduleKind = omp_sched_static;
    if (scheduleType == 1) scheduleKind = omp_sched_dynamic;
    else if (scheduleType == 2) scheduleKind = omp_sched_guided;
    omp_set_schedule(scheduleKind, chunkSize);

    #pragma omp parallel for schedule(runtime)
    for (int i = 0; i < getSize(); ++i) {
        double sum_neighbors = 0.0, Ai = getNodes()[i].getAmplitude();
        for (int neighborId : getNodes()[i].getNeighbors())
            sum_neighbors += getNodes()[neighborId].getAmplitude() - Ai;
        double diffusion = D * sum_neighbors, damping = -gamma * Ai, delta = diffusion + damping;
        newAmplitudes[i] = Ai + delta * dt;
    }
    newAmplitudes[center] += noiseCoeff;
    for (int i = 0; i < getSize(); ++i)
        getNodes()[i].updateAmplitude(newAmplitudes[i]);
}

// Paralelo 2D con collapse
void Network::propagateWavesCollapse() {
    int side = static_cast<int>(std::round(std::sqrt(getSize())));
    if (side * side != getSize()) {
        std::cerr << "Error: La red debe ser cuadrada para usar collapse en 2D." << std::endl;
        return;
    }
    int center = static_cast<int>(std::round((side / 2.0) * side + (side / 2.0)));

    std::vector<double> newAmplitudes(getSize(), 0.0);
    double D = getDiffusionCoeff();
    double gamma = getDampingCoeff();
    double dt = getTimestep();

    #pragma omp parallel for collapse(2)
    for (int i = 0; i < side; ++i) {
        for (int j = 0; j < side; ++j) {
            int idx = i * side + j;
            double sum_neighbors = 0.0, Ai = getNodes()[idx].getAmplitude();
            if (i > 0)        sum_neighbors += getNodes()[(i - 1) * side + j].getAmplitude() - Ai;   
            if (i < side - 1) sum_neighbors += getNodes()[(i + 1) * side + j].getAmplitude() - Ai; 
            if (j > 0)        sum_neighbors += getNodes()[i * side + (j - 1)].getAmplitude() - Ai;  
            if (j < side - 1) sum_neighbors += getNodes()[i * side + (j + 1)].getAmplitude() - Ai;   
            double diffusion = D * sum_neighbors, damping = -gamma * Ai, delta = diffusion + damping;
            newAmplitudes[idx] = Ai + delta * dt;
        }
    }
    newAmplitudes[center] += noiseCoeff;
    for (int i = 0; i < getSize(); ++i)
        getNodes()[i].updateAmplitude(newAmplitudes[i]);
}
