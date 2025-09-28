#ifndef NODE_H
#define NODE_H

#include <vector>

class Node {
private:
    int id;
    double amplitude;
    double previousAmplitude;
    std::vector<int> neighbors;

public:
    Node(int nodeId);

    void addNeighbor(int neighborId);
    void clearNeighbors();              // Nuevo método para limpiar vecinos
    void updateAmplitude(double newAmplitude);

    double getAmplitude() const;
    double getPreviousAmplitude() const;

    const std::vector<int>& getNeighbors() const;  // Para lectura
    std::vector<int>& getNeighbors();              // Para modificación

    int getDegree() const;
};

#endif // NODE_H
