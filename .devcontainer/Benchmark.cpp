#include "Benchmark.h"
#include "Network.h"
#include "WavePropagator.h"
#include <cmath>
#include <vector>
#include <numeric>


Benchmark::~Benchmark() {
    if (resultsFile.is_open())
        resultsFile.close();
}

void Benchmark::openResultsFile(const std::string& filename) {
    resultsFile.open(filename, std::ios::out | std::ios::trunc);
}

void Benchmark::closeResultsFile() {
    if (resultsFile.is_open())
        resultsFile.close();
}



performanceMeasurement Benchmark::measureTime(std::vector<double> times) {
     // Calculate mean
    double sum = std::accumulate(times.begin(), times.end(), 0.0);
    double mean = sum / times.size();
    
    // Calculate standard deviation
    double variance_sum = 0.0;
    for (double value : times) {
        double diff = value - mean;
        variance_sum += diff * diff;
    }
    
    double variance = variance_sum / times.size();
    double stddev = std::sqrt(variance);

    return {mean, stddev};
}


performanceMeasurement Benchmark::measureSpeedup(std::vector<double> serialTimes, std::vector<double> parallelTimes) {
    performanceMeasurement serialMeasure = measureTime(serialTimes);
    performanceMeasurement parallelMeasure = measureTime(parallelTimes);

    double speedup = serialMeasure.time / parallelMeasure.time;

    double squareSerialErr = std::pow(serialMeasure.error / serialMeasure.time, 2);
    double squareParallelErr = std::pow(parallelMeasure.error / parallelMeasure.time, 2);

    double speedupErr = speedup * std::sqrt(squareSerialErr + squareParallelErr);

    return {speedup, speedupErr};
}

performanceMeasurement Benchmark::measureEfficiency(performanceMeasurement speedup, int currentThreads) {
    double efficiency = speedup.time / currentThreads;
    double error = speedup.error / currentThreads;

    return {efficiency, error};
}

amdahl Benchmark::calculateAmdahl(double initialTime, performanceMeasurement serialTime, int currentThreads) {
    double serialFraction = serialTime.time / (serialTime.time + initialTime);
    double rightPartDenominator = (1 - serialFraction) / currentThreads;
    double denominator = serialFraction + rightPartDenominator;
    double amdahl = 1 / denominator;

    return {amdahl, serialFraction};
}

performanceAnalysis Benchmark::analyzePerformance(std::vector<double> serialTimes, std::vector<double> parallelTimes, double initialTime, int currentThreads) {
    performanceMeasurement serialTime = measureTime(serialTimes);
    performanceMeasurement parallelTime = measureTime(parallelTimes);
    performanceMeasurement speedup = measureSpeedup(serialTimes, parallelTimes);
    performanceMeasurement efficiency = measureEfficiency(speedup, currentThreads);
    // todo
    amdahl amdahl = calculateAmdahl(initialTime, serialTime, currentThreads);
    return {serialTime, parallelTime, speedup, efficiency, amdahl};
}

void Benchmark::writeSchedulesPerformanceHeader() {
    resultsFile << "##################################################\n"
                << "#           Comparación de Schedules             #\n"
                << "##################################################\n"
                << "\n"
                << "N° de threads, tiempo serial, tiempo paralelo, speedup, eficiencia, amdahl, tipo de schedule, tamaño de chunk\n";
    return;
};

void Benchmark::writeSynchronizationPerformanceHeader() {
    resultsFile << "\n"
                << "##################################################\n"
                << "#         Comparación de Sincronización          #\n"
                << "##################################################\n"
                << "\n"
                << "N° de threads, tiempo serial, tiempo paralelo, speedup, eficiencia, amdahl, tipo de sincronizacion, useBarrier\n";
    return;
};

void Benchmark::writeDataPerformanceHeader() {
    resultsFile << "\n"
                << "##################################################\n"
                << "#               Manejo de datos                  #\n"
                << "##################################################\n"
                << "\n"
                << "N° de threads, tiempo serial, tiempo paralelo, speedup, eficiencia, amdahl, tipo de metodo, usePrivate\n";
    return;
};

void Benchmark::writeAdvancedSynchronizationPerformanceHeader() {
    resultsFile << "\n"
                << "##################################################\n"
                << "#            Sincronización avanzada             #\n"
                << "##################################################\n"
                << "\n"
                << "N° de threads, tiempo serial, tiempo paralelo, speedup, eficiencia, amdahl, tipo de tarea, useSingle\n";
    return;
};

void Benchmark::writeScalabilityHeader() {
    resultsFile << "\n"
                << "##################################################\n"
                << "#         Escalabilidad y Fraccion serial        #\n"
                << "##################################################\n"
                << "\n"
                << "N° de threads, tiempo serial, tiempo paralelo, speedup, eficiencia, fracción serial, amdahl\n";
    return;
};

void Benchmark::writeScalabilityResults(performanceAnalysis performance, int currentThreads) {
    resultsFile << currentThreads << ", "
                << performance.serialTime.time << " ± " << performance.serialTime.error << ", "
                << performance.parallelTime.time << " ± " << performance.parallelTime.error << ", "
                << performance.speedup.time << " ± " << performance.speedup.error << ", "
                << performance.efficiency.time << " ± " << performance.efficiency.error << ", "
                << performance.amdahl.serialFraction << ", "
                << performance.amdahl.amdahl << "\n";
    return;
};


void Benchmark::writeSchedulesPerformance(performanceAnalysis performance, int scheduleType, int chunkSize) {
    resultsFile << threads << ", "
                << performance.serialTime.time << " ± " << performance.serialTime.error << ", "
                << performance.parallelTime.time << " ± " << performance.parallelTime.error << ", "
                << performance.speedup.time << " ± " << performance.speedup.error << ", "
                << performance.efficiency.time << " ± " << performance.efficiency.error << ", "
                << performance.amdahl.amdahl << ", "
                << scheduleType << ", "
                << chunkSize << "\n";
    return;
};

void Benchmark::writeSynchronizationPerformance(performanceAnalysis performance, int syncType, bool useBarrier) {
    resultsFile << threads << ", "
                << performance.serialTime.time << " ± " << performance.serialTime.error << ", "
                << performance.parallelTime.time << " ± " << performance.parallelTime.error << ", "
                << performance.speedup.time << " ± " << performance.speedup.error << ", "
                << performance.efficiency.time << " ± " << performance.efficiency.error << ", "
                << performance.amdahl.amdahl << ", "
                << syncType << ", "
                << useBarrier << "\n";
    return;
};

void Benchmark::writeDataPerformance(performanceAnalysis performance, int methodType, bool usePrivate) {
    resultsFile << threads << ", "
                << performance.serialTime.time << " ± " << performance.serialTime.error << ", "
                << performance.parallelTime.time << " ± " << performance.parallelTime.error << ", "
                << performance.speedup.time << " ± " << performance.speedup.error << ", "
                << performance.efficiency.time << " ± " << performance.efficiency.error << ", "
                << performance.amdahl.amdahl << ", "
                << methodType << ", "
                << usePrivate << "\n";
    return;
};
            
void Benchmark::writeAdvancedSynchronizationPerformance(performanceAnalysis performance, int taskType, bool useSingle) {
    resultsFile << threads << ", "
                << performance.serialTime.time << " ± " << performance.serialTime.error << ", "
                << performance.parallelTime.time << " ± " << performance.parallelTime.error << ", "
                << performance.speedup.time << " ± " << performance.speedup.error << ", "
                << performance.efficiency.time << " ± " << performance.efficiency.error << ", "
                << performance.amdahl.amdahl << ", "
                << taskType << ", "
                << useSingle << "\n";
    return;
};

void Benchmark::writeEnergyConservationHeader() {
    resultsFile << "\n"
                << "##################################################\n"
                << "#             Conservación de Energía            #\n"
                << "##################################################\n"
                << "\n"
                << "Step, tiempo, energía total\n";
    return;
};

void Benchmark::writeEnergyConservation(double step, double time, double energy) {
    resultsFile << step << ", "
                << time << ", "
                << energy << "\n";
    return;
};

void Benchmark::executeBenchmark() {
    openResultsFile("benchmark_results.dat");
    writeSchedulesPerformanceHeader();
    compareSchedules();
    writeSynchronizationPerformanceHeader();
    compareSynchronization();
    writeDataPerformanceHeader();
    compareSynchronization();
    writeAdvancedSynchronizationPerformanceHeader();
    compareAdvancedSynchronization();
    closeResultsFile();
};

void Benchmark::executeScalingAnalysis() {
    openResultsFile("scaling_analysis.dat");
    writeScalabilityHeader();
    analyzeScalability();
    closeResultsFile();
};

void Benchmark::analyzeScalability() {
    double initialT0 = omp_get_wtime();
    Network network(networkSize, D, gamma, dt);
    network.initializeRegularNetwork(2);
    network.getNodes()[0].updateAmplitude(amplitude0);

    double initialT1 = omp_get_wtime();
    double initialTime = initialT1 - initialT0;

    std::vector<int> threadOptions = {1, 2, 3, 4};
    

    // Tiempos paralelos
    for (int currentThreads : threadOptions) {
        // Tiempos seriales
        omp_set_num_threads(1);
        std::vector<double> serialTimes;
        for (int j = 0; j < repeats; j++) {
            double t0 = omp_get_wtime();
            network.propagateWaves(0, 16);
            double t1 = omp_get_wtime();
            double time = t1 - t0;
            serialTimes.push_back(time);
        }

        omp_set_num_threads(currentThreads);
        std::vector<double> parallelTimes;
        for (int j = 0; j < repeats; j++) {
            double t0 = omp_get_wtime();
            network.propagateWaves(0, 16);
            double t1 = omp_get_wtime();
            double time = t1 - t0;
            parallelTimes.push_back(time);
        }

        performanceAnalysis performance = analyzePerformance(serialTimes, parallelTimes, initialTime, currentThreads);
        writeScalabilityResults(performance, currentThreads);

    }
};

            
void Benchmark::compareSchedules() {
    double initialT0 = omp_get_wtime();
    Network network(networkSize, D, gamma, dt);
    network.initializeRegularNetwork(2);
    network.getNodes()[0].updateAmplitude(amplitude0);

    double initialT1 = omp_get_wtime();
    double initialTime = initialT1 - initialT0;


    // Tiempos seriales
    for (int scheduleType : scheduleTypes) {
        for (int chunkSize: chunkSizes) {
            omp_set_num_threads(1);
            std::vector<double> serialTimes;
            for (int j = 0; j < repeats; j++) {
                double t0 = omp_get_wtime();
                network.propagateWaves(scheduleType, chunkSize);
                double t1 = omp_get_wtime();
                double time = t1 - t0;
                serialTimes.push_back(time);
            }
            omp_set_num_threads(threads);
            std::vector<double> parallelTimes;
            for (int j = 0; j < repeats; j++) {
                double t0 = omp_get_wtime();
                network.propagateWaves(scheduleType, chunkSize);
                double t1 = omp_get_wtime();
                double time = t1 - t0;
                parallelTimes.push_back(time);
            }

            performanceAnalysis performance =  analyzePerformance(serialTimes, parallelTimes, initialTime, threads);
            writeSchedulesPerformance(performance, scheduleType, chunkSize);
        }
    }
};

void Benchmark::compareSynchronization() {

    double initialT0 = omp_get_wtime();
    Network network(networkSize, D, gamma, dt);
    network.initializeRegularNetwork(2);
    network.getNodes()[0].updateAmplitude(amplitude0);
    WavePropagator propagator(network);

    double initialT1 = omp_get_wtime();
    double initialTime = initialT1 - initialT0;

    for (int syncType : syncTypes) {
        for (bool useBarrier: booleans) {
            // Tiempos seriales
            omp_set_num_threads(1);
            std::vector<double> serialTimes;
            for (int j = 0; j < repeats; j++) {
                double t0 = omp_get_wtime();
                propagator.integrateEuler(syncType, useBarrier);
                double t1 = omp_get_wtime();
                double time = t1 - t0;
                serialTimes.push_back(time);
            }
            // Tiempos paralelos
            omp_set_num_threads(threads);
            std::vector<double> parallelTimes;
            for (int j = 0; j < repeats; j++) {
                double t0 = omp_get_wtime();
                propagator.integrateEuler(syncType, useBarrier);
                double t1 = omp_get_wtime();
                double time = t1 - t0;
                parallelTimes.push_back(time);
            }

            performanceAnalysis performance = analyzePerformance(serialTimes, parallelTimes, initialTime, threads);
            writeSynchronizationPerformance(performance, syncType, useBarrier);
        }
    }
};



void Benchmark::compareData() {
    double initialT0 = omp_get_wtime();
    Network network(networkSize, D, gamma, dt);
    network.initializeRegularNetwork(2);
    network.getNodes()[0].updateAmplitude(amplitude0);
    WavePropagator propagator(network);
    
    double initialT1 = omp_get_wtime();
    double initialTime = initialT1 - initialT0;
    
    for (int dataManagementType : dataManagementTypes) {
        for (bool usePrivate: booleans) {
            // Tiempos seriales
            omp_set_num_threads(1);
            std::vector<double> serialTimes;
            for (int j = 0; j < repeats; j++) {
                double t0 = omp_get_wtime();
                propagator.calculateEnergy(0, 1);
                double t1 = omp_get_wtime();
                double time = t1 - t0;
                serialTimes.push_back(time);
            }

            // Tiempos paralelos
            omp_set_num_threads(threads);
            std::vector<double> parallelTimes;
            for (int j = 0; j < repeats; j++) {
                double t0 = omp_get_wtime();
                propagator.calculateEnergy(0, 1);
                double t1 = omp_get_wtime();
                double time = t1 - t0;
                parallelTimes.push_back(time);
            }
            performanceAnalysis performance = analyzePerformance(serialTimes, parallelTimes, initialTime, threads);
            writeDataPerformance(performance, dataManagementType, usePrivate);


        }
    }
};






void Benchmark::compareAdvancedSynchronization() {
    double initialT0 = omp_get_wtime();
    Network network(networkSize, D, gamma, dt);
    network.initializeRegularNetwork(2);
    network.getNodes()[0].updateAmplitude(amplitude0);
    WavePropagator propagator(network);

    double initialT1 = omp_get_wtime();
    double initialTime = initialT1 - initialT0;

    for (int taskType : taskTypes) {
        for (bool useSingle: booleans) {
            // Tiempos seriales
            omp_set_num_threads(1);
            std::vector<double> serialTimes;
            for (int j = 0; j < repeats; j++) {
                double t0 = omp_get_wtime();
                propagator.processNodes(taskType, useSingle);
                double t1 = omp_get_wtime();
                double time = t1 - t0;
                serialTimes.push_back(time);
            }
            // Tiempos paralelos
            omp_set_num_threads(threads);
            std::vector<double> parallelTimes;
            for (int j = 0; j < repeats; j++) {
                double t0 = omp_get_wtime();
                propagator.processNodes(taskType, useSingle);
                double t1 = omp_get_wtime();
                double time = t1 - t0;
                parallelTimes.push_back(time);
            }

            performanceAnalysis performance = analyzePerformance(serialTimes, parallelTimes, initialTime, threads);
            writeAdvancedSynchronizationPerformance(performance, taskType, useSingle);
        }
    }
    
};

void Benchmark::measureEnergyConservation() {
    openResultsFile("energy_conservation.dat");
    writeEnergyConservationHeader();
    int totalSteps = 1000;
    const double sourceMagnitude = 0.1;
    const double omega = 2.0 * 3.14159265359;

    Network network(networkSize, D, gamma, dt);
    network.initializeRegularNetwork(2);
    network.getNodes()[0].updateAmplitude(amplitude0);
    WavePropagator propagator(network);

    for (int step = 0; step < totalSteps; ++step) {
        double time = step * dt;
        // double source = sourceMagnitude * std::sin(omega * time);
        // network.getNodes()[0].updateAmplitude(network.getNodes()[0].getAmplitude() + dt * source);

        network.propagateWaves();
        double totalEnergy = propagator.calculateEnergy();
        writeEnergyConservation(step, time, totalEnergy);
    };

    closeResultsFile();

    
};