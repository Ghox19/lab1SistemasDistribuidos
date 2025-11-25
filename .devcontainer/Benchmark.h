#ifndef BENCHMARK_H
#define BENCHMARK_H

#include <vector>
#include <fstream>
#include "Node.h"

struct performanceMeasurement {
    double time;
    double error;
};

struct amdahl {
    double amdahl;
    double serialFraction;
};

struct performanceAnalysis {
    performanceMeasurement serialTime;
    performanceMeasurement parallelTime;
    performanceMeasurement speedup;
    performanceMeasurement efficiency;
    amdahl amdahl;
};

class Benchmark {
private:
    std::ofstream resultsFile;
    int repeats = 10;
    int threads = 4;
    int networkSize = 100000;
    double dt = 0.01;
    double D = 0.1;
    double gamma = 0.01;
    double amplitude0 = 1.0;
    std::vector<int> chunkSizes = {1, 4, 16, 64};
    std::vector<int> scheduleTypes = {0, 1, 2};
    std::vector<int> syncTypes = {0, 1, 2};
    std::vector<int> dataManagementTypes = {0, 1};
    std::vector<int> taskTypes = {0, 1};
    std::vector<bool> booleans = {false, true};
    performanceMeasurement measureTime(std::vector<double> times);
    performanceMeasurement measureSpeedup(std::vector<double> serialTimes, std::vector<double> parallelTimes);
    performanceMeasurement measureEfficiency(performanceMeasurement speedup, int currentThreads);
    amdahl calculateAmdahl(double initialTime, performanceMeasurement serialTime, int currentThreads);
    performanceAnalysis analyzePerformance(std::vector<double> serialTimes, std::vector<double> parallelTimes, double initialTime, int currentThreads);
    void writeSchedulesPerformanceHeader();
    void writeSynchronizationPerformanceHeader();
    void writeDataPerformanceHeader();
    void writeAdvancedSynchronizationPerformanceHeader();
    void writeSchedulesPerformance(performanceAnalysis performance, int scheduleType, int chunkSize);
    void writeSynchronizationPerformance(performanceAnalysis performance, int syncType, bool useBarrier);
    void writeDataPerformance(performanceAnalysis performance, int methodType, bool usePrivate);
    void writeAdvancedSynchronizationPerformance(performanceAnalysis performance, int taskType, bool useSingle);
    void compareSchedules();
    void compareSynchronization();
    void compareData();
    void compareAdvancedSynchronization();
    void analyzeScalability();
    void writeScalabilityHeader();
    void writeScalabilityResults(performanceAnalysis performance, int currentThreads);
    void openResultsFile(const std::string& filename);
    void closeResultsFile();
    void writeEnergyConservationHeader();
    void writeEnergyConservation(double step, double time, double energy);
    
    
public:
    Benchmark() = default;
    ~Benchmark();
    void executeScalingAnalysis();
    void executeBenchmark();
    void measureEnergyConservation();

};

#endif
