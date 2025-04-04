#include "DrainageNetwork.h"
#include "BayesianOptimizer.h"
#include <iostream>
int main() {
    const int NUM_NODES = 10;
    const int NUM_PIPES = 15;
    const int NUM_BASINS = 2;
    DrainageNetwork network(NUM_NODES, NUM_PIPES, NUM_BASINS);
    size_t dim = 5 * NUM_PIPES + NUM_BASINS + 2;
    std::vector<double> lb(dim);
    std::vector<double> ub(dim);
    size_t idx = 0;
    for (int i = 0; i < NUM_PIPES; ++i) {
        lb[idx] = 0.15; ub[idx++] = 2.0;
        lb[idx] = 0.001; ub[idx++] = 0.05;
        lb[idx] = 0.01; ub[idx++] = 0.015;
        lb[idx] = 0.1; ub[idx++] = 1.0;
        lb[idx] = 0.01; ub[idx++] = 0.5;
    }
    for (int i = 0; i < NUM_BASINS; ++i) {
        lb[idx] = 100.0; ub[idx++] = 10000.0;
    }
    lb[idx] = 1.0; ub[idx++] = 50.0;
    lb[idx] = 0.3; ub[idx++] = 0.9;
    BayesianOptimizer optimizer(dim, lb, ub);
    auto objective = [&](const std::vector<double>& params) {
        network.setParameters(params);
        double rainfallIntensity = 50.0;
        double performance = network.evaluatePerformance(rainfallIntensity);
        double cost = network.calculateCost();
        return performance - cost / 1e8;
    };
    optimizer.optimize(objective, 100, network);
    std::cout << "Optimal performance-cost score: " << optimizer.getBestValue() << std::endl;
    std::cout << "Total cost: $" << network.calculateCost() << std::endl;
    std::cout << "Flow capacity: " << network.calculateFlowCapacity() << " m³/s" << std::endl;
    std::cout << "Performance ratio: " << network.evaluatePerformance(50.0) << std::endl;
    auto bestParams = optimizer.getBestParameters();
    idx = 0;
    std::cout << "\nOptimal pipe parameters:\n";
    for (size_t i = 0; i < NUM_PIPES; ++i) {
        std::cout << "Pipe " << i << ":\n";
        std::cout << "  Diameter: " << bestParams[idx++] << "m\n";
        std::cout << "  Slope: " << bestParams[idx++] << " (" << (bestParams[idx-1] * 100) << "%)\n";
        std::cout << "  Roughness: " << bestParams[idx++] << "\n";
        std::cout << "  Runoff Coeff: " << bestParams[idx++] << "\n";
        std::cout << "  Inlet Capacity: " << bestParams[idx++] << " m³/s\n";
    }
    std::cout << "\nOptimal basin capacities:\n";
    for (size_t i = 0; i < NUM_BASINS; ++i) {
        std::cout << "Basin " << i << ": " << bestParams[idx++] << " m³\n";
    }
    std::cout << "Infiltration Rate: " << bestParams[idx++] << " mm/h\n";
    std::cout << "Sediment Efficiency: " << bestParams[idx++] << "\n";
    return 0;
}
