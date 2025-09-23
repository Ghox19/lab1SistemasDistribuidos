#include <iostream>
#include "Network.h"
#include "Node.h"

// Utiliza los valores del ejemplo: D=0.1, gamma=0.01, dt=0.01
int main() {
    double D = 0.1;       // coeficiente de difusión
    double gamma = 0.01;  // coeficiente de amortiguación
    double dt = 0.01;     // paso temporal

    // Crear una red de 3 nodos: 0, 1, 2
    Network network(3, D, gamma, dt);

    // Limpiar vecinos y asignar conectividad: nodo 1 con vecinos 0 y 2
    network.getNodes()[0].clearNeighbors();
    network.getNodes()[1].clearNeighbors();
    network.getNodes()[2].clearNeighbors();

    network.getNodes()[1].addNeighbor(0);
    network.getNodes()[1].addNeighbor(2);

    // Asignar amplitudes como dice el ejemplo:
    // Nodo 1: 1.0, vecino 0: 0.8, vecino 2: 1.2
    network.getNodes()[0].updateAmplitude(0.8);
    network.getNodes()[1].updateAmplitude(1.0);
    network.getNodes()[2].updateAmplitude(1.2);

    // Mostrar amplitudes iniciales
    std::cout << "Amplitudes iniciales:" << std::endl;
    for (int i = 0; i < 3; ++i) {
        std::cout << "Nodo " << i << ": " << network.getNodes()[i].getAmplitude() << std::endl;
    }

    // Ejecutar un solo paso de propagación
    network.propagateWaves();

    // Mostrar amplitudes tras la propagación
    std::cout << "\nAmplitudes después de 1 paso:" << std::endl;
    for (int i = 0; i < 3; ++i) {
        std::cout << "Nodo " << i << ": " << network.getNodes()[i].getAmplitude() << std::endl;
    }

    // Opcional: mostrar el cálculo detallado para el nodo central
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
