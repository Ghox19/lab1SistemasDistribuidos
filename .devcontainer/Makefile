CXX = g++
CXXFLAGS = -Wall -Wextra -O3 -fopenmp -std=c++17
LDFLAGS = -fopenmp

TARGET = main
SOURCES = main.cpp Node.cpp Network.cpp WavePropagator.cpp MetricsCalculator.cpp Benchmark.cpp
HEADERS = Node.h Network.h WavePropagator.h MetricsCalculator.h Benchmark.h

$(TARGET): $(SOURCES) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SOURCES) $(LDFLAGS)

clean:
	rm -f $(TARGET) *.o *.dat *.png

benchmark: $(TARGET)
	./$(TARGET) -benchmark

analysis: $(TARGET)
	./$(TARGET) -analysis

.PHONY: clean benchmark analysis