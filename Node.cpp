#include "Node.h"

Node::Node(int nodeId) : id(nodeId), amplitude(0.0), previousAmplitude(0.0) {}

void Node::addNeighbor(int neighborId) {
    neighbors.push_back(neighborId);
}

void Node::clearNeighbors() {
    neighbors.clear();
}

void Node::updateAmplitude(double newAmplitude) {
    previousAmplitude = amplitude;
    amplitude = newAmplitude;
}

double Node::getAmplitude() const {
    return amplitude;
}

double Node::getPreviousAmplitude() const {
    return previousAmplitude;
}

const std::vector<int>& Node::getNeighbors() const {
    return neighbors;
}

std::vector<int>& Node::getNeighbors() {
    return neighbors;
}

int Node::getDegree() const {
    return neighbors.size();
}
