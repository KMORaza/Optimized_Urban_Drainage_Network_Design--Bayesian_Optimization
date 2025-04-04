#ifndef DRAINAGE_NETWORK_H
#define DRAINAGE_NETWORK_H
#include <vector>
#include <random>
class DrainageNetwork {
public:
    struct Pipe {
        double diameter;
        double length;
        double slope;
        double roughness;
        int upstreamNode;
        int downstreamNode;
        double runoffCoeff;
        double inletCapacity;
        double catchmentArea;
        double imperviousness;
    };
    struct Node {
        double elevation;
        double demand;
        double x, y;
    };
    struct DetentionBasin {
        double capacity;
        double outletRate;
    };
    DrainageNetwork(int numNodes, int numPipes, int numBasins);
    double calculateCost() const;
    double calculateFlowCapacity() const;
    double evaluatePerformance(double rainfallIntensity) const;
    void setParameters(const std::vector<double>& params);
    const std::vector<Pipe>& getPipes() const {
        return pipes;
        }
    const std::vector<Node>& getNodes() const {
        return nodes;
        }
    const std::vector<DetentionBasin>& getBasins() const {
        return basins;
        }
private:
    std::vector<Pipe> pipes;
    std::vector<Node> nodes;
    std::vector<DetentionBasin> basins;
    double rainfallDuration;
    double rainfallFrequency;
    double infiltrationRate;
    double sedimentEfficiency;
    double calculatePipeFlow(const Pipe& pipe, double depth) const;
    double calculateRunoff(const Pipe& pipe, double rainfallIntensity) const;
    double calculateVelocity(const Pipe& pipe, double flow) const;
};
#endif
