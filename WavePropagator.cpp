#include "WavePropagator.h"
#include <vector>

void WavePropagator::propagateWaves() {
    /*
    Esta función implementa la ecuación diferencial en forma discreta para
    actualizar la amplitud A_i de cada nodo i en la red:

    A_i(t + Δt) = A_i(t) + Δt * [ D * Σ_{j ∈ vecinos(i)} (A_j(t) - A_i(t))
                               - γ * A_i(t)
                               + S_i(t) ]

    donde:
    - D es el coeficiente de difusión,
    - γ es el coeficiente de amortiguación,
    - S_i(t) es la fuente externa (aquí asumida 0),
    - Σ_{j ∈ vecinos(i)} es la suma sobre los nodos vecinos de i.

    El procedimiento es:
    1. Para cada nodo i, calcular la suma de diferencias de amplitud con sus vecinos.
    2. Multiplicar la suma por el coeficiente de difusión D.
    3. Calcular el término de amortiguación, proporcional a la amplitud actual.
    4. Asumir fuente 0 para esta versión.
    5. Calcular el cambio total multiplicado por el paso temporal Δt.
    6. Calcular la nueva amplitud como la suma del valor actual y el cambio calculado.
    7. Actualizar todas las amplitudes al finalizar el ciclo para evitar dependencias.

    Este enfoque asegura que la actualización se base en valores del tiempo actual (t),
    manteniendo la estabilidad y evitando la dependencia de valores ya modificados
    en el paso t+Δt.

    Implicación física:
    - La difusión hace que la amplitud tienda a igualarse con vecinos.
    - La amortiguación representa pérdida energética, disminuyendo amplitudes.
    - La ausencia de fuente implica que no hay entrada de nueva energía.
    */

    std::vector<double> newAmplitudes(network.getSize(), 0.0);
    double D = network.getDiffusionCoeff();
    double gamma = network.getDampingCoeff();
    double dt = network.getTimestep();

    // Para cada nodo calcular la ecuación
    for (int i = 0; i < network.getSize(); ++i) {
        double sum_neighbors = 0.0;
        double Ai = network.getNodes()[i].getAmplitude();

        // Sumar diferencias con vecinos
        for (int neighborId : network.getNodes()[i].getNeighbors()) {
            sum_neighbors += network.getNodes()[neighborId].getAmplitude() - Ai;
        }

        double diffusion = D * sum_neighbors;
        double damping = -gamma * Ai;
        double source = 0.0; // No hay fuente externa en esta versión
        double delta = diffusion + damping + source;

        newAmplitudes[i] = Ai + delta * dt;
    }

    // Actualizar amplitudes después del cálculo para todos los nodos
    for (int i = 0; i < network.getSize(); ++i) {
        network.getNodes()[i].updateAmplitude(newAmplitudes[i]);
    }
}

void WavePropagator::propagateWaves(int scheduleType) {
    // Por implementar
}

void WavePropagator::propagateWaves(int scheduleType, int chunkSize) {
    // Por implementar
}

void WavePropagator::propagateWavesCollapse() {
    // Por implementar
}
