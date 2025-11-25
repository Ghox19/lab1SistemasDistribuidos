#include <iostream>
#include <fstream>
#include <cmath>
#include "Network.h"
#include "WavePropagator.h"

int main() {
    const int networkWidth = 100;
    const int networkHeight = 100;
    const int networkSize = networkWidth * networkHeight;
    const int totalSteps = 1000;
    const double dt = 0.01;
    const double D = 0.1;
    const double gamma = 0.01;
    const double amplitude0 = 1.0;
    const double sourceMagnitude = 0.1;
    const double omega = 2.0 * 3.14159265359;

    Network network(networkSize, D, gamma, dt);

    // Vecinos 2D
    for (int i = 0; i < networkSize; ++i) {
        network.getNodes()[i].clearNeighbors();
        int x = i % networkWidth, y = i / networkWidth;
        if (y > 0) network.getNodes()[i].addNeighbor(i - networkWidth);
        if (y < networkHeight - 1) network.getNodes()[i].addNeighbor(i + networkWidth);
        if (x > 0) network.getNodes()[i].addNeighbor(i - 1);
        if (x < networkWidth - 1) network.getNodes()[i].addNeighbor(i + 1);
    }
    // Amplitudes iniciales
    for (int i = 0; i < networkSize; ++i) network.getNodes()[i].updateAmplitude(0.0);
    int center = (networkHeight / 2) * networkWidth + (networkWidth / 2);
    network.getNodes()[center].updateAmplitude(amplitude0);


    std::ofstream out("salida2d.dat");
    out << "# Simulacion 2D " << networkWidth << "x" << networkHeight << " pasos=" << totalSteps << " dt=" << dt << std::endl;
    out << "# D=" << D << " gamma=" << gamma << " amp0=" << amplitude0 << std::endl;
    out << "# Fuente: " << sourceMagnitude << " sin(" << omega << "t)" << std::endl;
    out << "# Formato: step time x y amp" << std::endl;

    for (int step = 0; step < totalSteps; ++step) {
        double time = step * dt;
        double source = sourceMagnitude * std::sin(omega * time);
        // Si quieres usar la fuente, agrega el término aquí:
        // network.getNodes()[center].updateAmplitude(network.getNodes()[center].getAmplitude() + dt * source);

        network.propagateWaves(1);

        // Guardar todos los pasos sin filtro
        for (int y = 0; y < networkHeight; ++y) {
            for (int x = 0; x < networkWidth; ++x) {
                int idx = y * networkWidth + x;
                out << step << " " << time << " " << x << " " << y << " " << network.getNodes()[idx].getAmplitude() << "\n";
            }
        }
        out << std::endl; // línea vacía entre snapshots
    }
    out.close();

    std::cout << "Simulación completa y archivo .dat generado con TODOS los pasos." << std::endl;
    return 0;
}
