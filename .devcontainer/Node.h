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

    // Métodos para gestionar vecinos y amplitud
    void addNeighbor(int neighborId);
    void clearNeighbors();           
    void updateAmplitude(double newAmplitude);

    // Métodos para obtener amplitudes
    double getAmplitude() const;
    double getPreviousAmplitude() const;

    // Métodos para obtener vecinos
    const std::vector<int>& getNeighbors() const;
    std::vector<int>& getNeighbors();              

    int getDegree() const;
};

#endif // NODE_H
