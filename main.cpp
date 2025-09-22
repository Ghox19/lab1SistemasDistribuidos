#include <iostream>
#include <vector>
#include "Node.h"
#include "Network.h"

int main() {
    // Parámetros de la red
    int networkSize = 10;
    double diffusionCoeff = 0.1;  // No usado en esta prueba simple
    double dampingCoeff = 0.01;   // No usado en esta prueba simple

    // Crear la red
    Network network(networkSize, diffusionCoeff, dampingCoeff);

    // Conectar nodos en línea (topología simple 1D)
    for (int i = 0; i < networkSize; ++i) {
        if (i > 0)
            network.getNodes()[i].addNeighbor(i - 1);
        if (i < networkSize - 1)
            network.getNodes()[i].addNeighbor(i + 1);
    }

    // Asignar amplitudes iniciales
    for (int i = 0; i < networkSize; ++i) {
        double initialAmplitude = (i == networkSize / 2) ? 1.0 : 0.0;
        network.getNodes()[i].updateAmplitude(initialAmplitude);
    }

    // Mostrar información de cada nodo
    for (const Node& node : network.getNodes()) {
        std::cout << "Nodo " << node.getDegree() << " con vecinos: ";
        for (int neighbor : node.getNeighbors()) {
            std::cout << neighbor << " ";
        }
        std::cout << "\nAmplitud actual: " << node.getAmplitude() 
                  << ", amplitud previa: " << node.getPreviousAmplitude() << std::endl;
    }

    return 0;
}
