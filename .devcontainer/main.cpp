#include <iostream>
#include <fstream>
#include <cmath>
#include "Network.h"

int main(int argc, char* argv[]) {
    const int networkWidth = 100;
    const int networkHeight = 100;
    const int networkSize = networkWidth * networkHeight;
    const int totalSteps = 1000;
    const double dt = 0.01;
    const double D = 0.1;
    const double gamma = 0.01;
    const double amplitude0 = 1.0;
    const double sourceMagnitude = 0.1;
    const double frequency = 1.0; 
    const double omega = 2.0 * M_PI * frequency; 

    Network network(networkSize, D, gamma, omega, dt);

    network.initializeRegularNetwork(2);

    if (argc > 1) {
        std::string mode = argv[1];
        if (mode == std::string("-benchmark")) {
            std::cout << "Ejecutando benchmark...\n";
            // Código para modo benchmark
        } else if (mode == std::string("-analysis")) {
            std::cout << "Ejecutando análisis...\n";
            // Código para modo análisis
        }
    } else {
        std::cout << "Ejecutando simulación normal...\n";

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

            // Ahora prueba la función con collapse(2)
            network.propagateWavesCollapse();

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

        std::cout << "Simulación completa y archivo .dat generado con TODOS los pasos usando collapse(2)." << std::endl;

    }
    return 0;
    
}
