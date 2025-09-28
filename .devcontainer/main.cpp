#include <iostream>
#include "Network.h"
#include "Node.h"
#include "WavePropagator.h"

int main() {
    double D = 0.1;
    double gamma = 0.01;
    double dt = 0.01;

    // Crear red de 3 nodos: 0, 1, 2
    Network network(3, D, gamma, dt);

    // Limpiar vecinos y conectar: nodo 1 con 0 y 2
    network.getNodes()[0].clearNeighbors();
    network.getNodes()[1].clearNeighbors();
    network.getNodes()[2].clearNeighbors();

    network.getNodes()[1].addNeighbor(0);
    network.getNodes()[1].addNeighbor(2);

    // Asignar amplitudes: nodo 1=1.0, vecino 0=0.8, vecino 2=1.2
    network.getNodes()[0].updateAmplitude(0.8);
    network.getNodes()[1].updateAmplitude(1.0);
    network.getNodes()[2].updateAmplitude(1.2);

    // Mostrar amplitudes iniciales
    std::cout << "Amplitudes iniciales:" << std::endl;
    for (int i = 0; i < 3; ++i) {
        std::cout << "Nodo " << i << ": " << network.getNodes()[i].getAmplitude() << std::endl;
    }

    // Invocar propagación paralela
    WavePropagator propagator(network);
    int scheduleType = 0;   // 0: static, 1: dynamic, 2: guided, 3: auto
    int chunkSize = 1;
    propagator.propagateWaves(scheduleType, chunkSize);

    // Mostrar amplitudes tras la propagación
    std::cout << "\nAmplitudes después de 1 paso paralelo:" << std::endl;
    for (int i = 0; i < 3; ++i) {
        std::cout << "Nodo " << i << ": " << network.getNodes()[i].getAmplitude() << std::endl;
    }

    // Validación manual
    double neighbors_diff = (0.8 - 1.0) + (1.2 - 1.0);
    double diffusion = D * neighbors_diff;
    double damping = -gamma * 1.0;
    double total_delta = dt * (diffusion + damping);
    double expected = 1.0 + total_delta;
    std::cout << "\nCálculo manual nodo 1: " << std::endl;
    std::cout << "Dif. vecinos = " << neighbors_diff << std::endl;
    std::cout << "Difusión = " << diffusion << std::endl;
    std::cout << "Amortiguación = " << damping << std::endl;
    std::cout << "Delta total = " << total_delta << std::endl;
    std::cout << "Nuevo valor esperado = " << expected << std::endl;

    return 0;
}
