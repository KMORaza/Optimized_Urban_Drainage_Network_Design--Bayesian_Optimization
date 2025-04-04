#include "BayesianOptimizer.h"
#include <random>
#include <Eigen/Dense>
using namespace Eigen;
BayesianOptimizer::BayesianOptimizer(size_t dim, const std::vector<double>& lb, const std::vector<double>& ub)
    : dimension(dim), lowerBounds(lb), upperBounds(ub), bestValue(-std::numeric_limits<double>::max()) {
    bestParams.resize(dim, 0.0);
}
double BayesianOptimizer::kernel(const std::vector<double>& x1, const std::vector<double>& x2) const {
    double sum = 0.0;
    for (size_t i = 0; i < dimension; ++i) {
        double diff = x1[i] - x2[i];
        sum += diff * diff;
    }
    return std::exp(-sum / (2.0 * dimension));
}
void BayesianOptimizer::updateGaussianProcess() {
    size_t n = observations.size();
    MatrixXd K(n, n);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            K(i, j) = kernel(observations[i].x, observations[j].x);
        }
    }
    K += MatrixXd::Identity(n, n) * 1e-8;
}
double BayesianOptimizer::expectedImprovement(const std::vector<double>& x) const {
    double mean = 0.0;
    double variance = 1.0;
    for (const auto& obs : observations) {
        mean += kernel(x, obs.x) * obs.y;
    }
    mean /= observations.size();
    double z = (mean - bestValue) / std::sqrt(variance);
    return (mean - bestValue) * 0.5 * (1.0 + std::erf(z / std::sqrt(2.0))) +
           std::sqrt(variance) * std::exp(-z * z / 2.0) / std::sqrt(2.0 * M_PI);
}
std::vector<double> BayesianOptimizer::sampleNextPoint() const {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::vector<double> bestX(dimension);
    double bestEI = -std::numeric_limits<double>::max();
    for (int i = 0; i < 100; ++i) {
        std::vector<double> x(dimension);
        for (size_t j = 0; j < dimension; ++j) {
            std::uniform_real_distribution<> dis(lowerBounds[j], upperBounds[j]);
            x[j] = dis(gen);
        }
        double ei = expectedImprovement(x);
        if (ei > bestEI) {
            bestEI = ei;
            bestX = x;
        }
    }
    return bestX;
}
void BayesianOptimizer::optimize(const std::function<double(const std::vector<double>&)>& objective,
                                int maxIterations, DrainageNetwork& network) {
    std::random_device rd;
    std::mt19937 gen(rd());
    for (int i = 0; i < 5; ++i) {
        std::vector<double> x(dimension);
        for (size_t j = 0; j < dimension; ++j) {
            std::uniform_real_distribution<> dis(lowerBounds[j], upperBounds[j]);
            x[j] = dis(gen);
        }
        network.setParameters(x);
        double y = objective(x);
        observations.push_back({x, y});
        if (y > bestValue) {
            bestValue = y;
            bestParams = x;
        }
    }
    for (int iter = 0; iter < maxIterations - 5; ++iter) {
        updateGaussianProcess();
        std::vector<double> nextPoint = sampleNextPoint();
        network.setParameters(nextPoint);
        double y = objective(nextPoint);
        observations.push_back({nextPoint, y});
        if (y > bestValue) {
            bestValue = y;
            bestParams = nextPoint;
        }
    }
}
