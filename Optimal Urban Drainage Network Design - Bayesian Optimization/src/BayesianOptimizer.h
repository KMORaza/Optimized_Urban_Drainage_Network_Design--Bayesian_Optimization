#ifndef BAYESIAN_OPTIMIZER_H
#define BAYESIAN_OPTIMIZER_H
#include <vector>
#include <functional>
#include "DrainageNetwork.h"
class BayesianOptimizer {
public:
    BayesianOptimizer(size_t dim, const std::vector<double>& lb,
                     const std::vector<double>& ub);

    void optimize(const std::function<double(const std::vector<double>&)>& objective,
                 int maxIterations, DrainageNetwork& network);

    std::vector<double> getBestParameters() const { return bestParams; }
    double getBestValue() const { return bestValue; }
private:
    struct GPPoint {
        std::vector<double> x;
        double y;
    };
    size_t dimension;
    std::vector<double> lowerBounds;
    std::vector<double> upperBounds;
    std::vector<GPPoint> observations;
    std::vector<double> bestParams;
    double bestValue;
    double kernel(const std::vector<double>& x1, const std::vector<double>& x2) const;
    std::vector<double> sampleNextPoint() const;
    double expectedImprovement(const std::vector<double>& x) const;
    void updateGaussianProcess();
};
#endif
