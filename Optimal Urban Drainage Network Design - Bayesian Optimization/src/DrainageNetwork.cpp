#include "DrainageNetwork.h"
#include <cmath>
#include <numeric>
DrainageNetwork::DrainageNetwork(int numNodes, int numPipes, int numBasins)
    : rainfallDuration(1.0), rainfallFrequency(10.0) {
    nodes.resize(numNodes);
    pipes.resize(numPipes);
    basins.resize(numBasins);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1000.0);
    for (auto& node : nodes) {
        node.elevation = dis(gen) / 10.0;
        node.demand = dis(gen) / 10000.0;
        node.x = dis(gen);
        node.y = dis(gen);
    }
    for (int i = 0; i < numPipes; ++i) {
        pipes[i].diameter = 0.3;
        pipes[i].length = dis(gen) / 10.0;
        pipes[i].slope = 0.01;
        pipes[i].roughness = 0.013;
        pipes[i].upstreamNode = i % numNodes;
        pipes[i].downstreamNode = (i + 1) % numNodes;
        pipes[i].runoffCoeff = 0.5;
        pipes[i].inletCapacity = 0.05;
        pipes[i].catchmentArea = dis(gen) * 100.0;
        pipes[i].imperviousness = 0.7;
    }
    for (auto& basin : basins) {
        basin.capacity = 1000.0;
        basin.outletRate = 0.1;
    }
}
double DrainageNetwork::calculatePipeFlow(const Pipe& pipe, double depth) const {
    double R = depth / 2.0;
    double A = M_PI * depth * depth / 4.0;
    return (1.0 / pipe.roughness) * A * std::pow(R, 2.0/3.0) * std::sqrt(pipe.slope);
}
double DrainageNetwork::calculateRunoff(const Pipe& pipe, double rainfallIntensity) const {
    double intensityMperS = rainfallIntensity / 3600000.0;
    double infiltrationLoss = infiltrationRate / 3600000.0;
    double effectiveIntensity = std::max(0.0, intensityMperS - infiltrationLoss);
    return pipe.runoffCoeff * effectiveIntensity * pipe.catchmentArea * pipe.imperviousness;
}
double DrainageNetwork::calculateVelocity(const Pipe& pipe, double flow) const {
    double A = M_PI * pipe.diameter * pipe.diameter / 4.0;
    return flow / A;
}
double DrainageNetwork::calculateCost() const {
    double totalCost = 0.0;
    for (const auto& pipe : pipes) {
        double pipeCost = (150.0 * pipe.diameter + 50.0 * pipe.slope) * pipe.length;
        double inletCost = 1000.0 * pipe.inletCapacity;
        totalCost += pipeCost + inletCost;
    }
    for (const auto& basin : basins) {
        totalCost += 50.0 * basin.capacity;
    }
    return totalCost;
}
double DrainageNetwork::calculateFlowCapacity() const {
    double minCapacity = std::numeric_limits<double>::max();
    for (const auto& pipe : pipes) {
        double capacity = calculatePipeFlow(pipe, pipe.diameter);
        minCapacity = std::min(minCapacity, std::min(capacity, pipe.inletCapacity));
    }
    return minCapacity;
}
double DrainageNetwork::evaluatePerformance(double rainfallIntensity) const {
    double capacity = calculateFlowCapacity();
    double totalRunoff = 0.0;
    for (const auto& pipe : pipes) {
        totalRunoff += calculateRunoff(pipe, rainfallIntensity);
    }
    double basinStorageEffect = 0.0;
    for (const auto& basin : basins) {
        double stored = std::min(basin.capacity / (rainfallDuration * 3600.0), totalRunoff);
        basinStorageEffect += stored;
    }
    double effectiveRunoff = totalRunoff - basinStorageEffect;
    double performance = capacity / (effectiveRunoff + 0.001);
    double velocityPenalty = 0.0;
    for (const auto& pipe : pipes) {
        double velocity = calculateVelocity(pipe, calculateRunoff(pipe, rainfallIntensity));
        if (velocity < 0.6) velocityPenalty += 0.1;
        if (velocity > 3.0) velocityPenalty += 0.1;
    }
    double qualityScore = sedimentEfficiency * 0.5;
    return performance - velocityPenalty + qualityScore;
}
void DrainageNetwork::setParameters(const std::vector<double>& params) {
    size_t nPipes = pipes.size();
    size_t nBasins = basins.size();
    size_t idx = 0;
    for (size_t i = 0; i < nPipes && idx < params.size(); ++i) {
        pipes[i].diameter = std::max(0.15, std::min(2.0, params[idx++]));
        pipes[i].slope = std::max(0.001, std::min(0.05, params[idx++]));
        pipes[i].roughness = std::max(0.01, std::min(0.015, params[idx++]));
        pipes[i].runoffCoeff = std::max(0.1, std::min(1.0, params[idx++]));
        pipes[i].inletCapacity = std::max(0.01, std::min(0.5, params[idx++]));
    }
    for (size_t i = 0; i < nBasins && idx < params.size(); ++i) {
        basins[i].capacity = std::max(100.0, std::min(10000.0, params[idx++]));
    }
    if (idx < params.size()) infiltrationRate = std::max(1.0, std::min(50.0, params[idx++]));
    if (idx < params.size()) sedimentEfficiency = std::max(0.3, std::min(0.9, params[idx++]));
}
