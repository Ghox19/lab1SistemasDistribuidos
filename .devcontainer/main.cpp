#include <iostream>
#include <fstream>
#include <cmath>
#include "Network.h"
#include "Benchmark.h"

int run_benchmarks() {
    Benchmark benchmark;
    benchmark.executeBenchmark();
    return 0;
}

int run_analysis() {
    Benchmark benchmark;
    benchmark.executeScalingAnalysis();
    return 0;
}

int run_default() {
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


    // Amplitudes iniciales
    for (int i = 0; i < networkSize; ++i) network.getNodes()[i].updateAmplitude(0.0);
    int center = (networkHeight / 2) * networkWidth + (networkWidth / 2);
    network.getNodes()[center].updateAmplitude(amplitude0);

    std::ofstream out("wave_evolution.dat");
    out << "# Simulacion 2D " << networkWidth << "x" << networkHeight << " pasos=" << totalSteps << " dt=" << dt << std::endl;
    out << "# D=" << D << " gamma=" << gamma << " amp0=" << amplitude0 << std::endl;
    out << "# Fuente: " << sourceMagnitude << " sin(" << omega << "t)" << std::endl;
    out << "# Formato: step time x y amp" << std::endl;

    for (int step = 0; step < totalSteps; ++step) {
        double time = step * dt;
        
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

    std::cout << "Simulación completa y archivo .dat generado con TODOS los pasos." << std::endl;

    return 0;
}

int main(int argc, char* argv[]) {

     if (argc < 2) {
        std::cout << "Ejecutando simulación normal...\n";
        run_default();
        return 0;
    }

    std::string mode = argv[1];

    if (mode == "-benchmark") {
        std::cout << "Ejecutando benchmark...\n";
        run_benchmarks();
        return 0;
    } else if (mode == "-analysis") {
        std::cout << "Ejecutando análisis...\n";
        run_analysis();
        return 0;
    } else {
        std::cerr << "Opcion desconocida: " << mode << "\n";
        std::cerr << "Use -benchmark o -analysis\n";
        return 1;
    }

    return 0;
    
}
