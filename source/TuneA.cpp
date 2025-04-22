// TuneMethodA.cpp
// Implements Method A: two-parameter (global + shape) correction for tube lengths
// based on the structure of Tune.cpp

#include "Yin.h"
#include "sfx.h"
#include "Birl.h"
#include <juce_audio_devices/juce_audio_devices.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <limits>

extern float MSEwithPenalty(const float* desiredFreqs);
extern void printFreqs();

// Method A tuning: two-parameter correction (global + shape)
void tuneMethodA(const float* desiredFreqs, int gridSize) {
    const int N = NUM_OF_TONEHOLES + 1;

    // 1) Read theoretical tube lengths
    std::vector<double> L_phys(N);
    for (int j = 0; j < N; ++j)
        L_phys[j] = birl::SFXPhysicalModelGetTubeLength(j);

    // 2) Global direction v1 = L_phys / ||L_phys||
    double norm_L = 0.0;
    for (double L : L_phys) norm_L += L * L;
    norm_L = std::sqrt(norm_L);
    std::vector<double> v1(N, 0.0);
    if (norm_L > 0.0) {
        for (int j = 0; j < N; ++j)
            v1[j] = L_phys[j] / norm_L;
    }

    // 3) Local shape v2 via one SPSA step
    double delta = 1e-4;
    std::vector<double> v_rand(N);
    for (int j = 0; j < N; ++j)
        v_rand[j] = (std::rand() % 2 ? 1.0 : -1.0);

    // Perturb +delta
    for (int j = 0; j < N; ++j)
        birl::SFXPhysicalModelSetTubeLength(j, L_phys[j] + delta * v_rand[j]);
    double lossPlus = MSEwithPenalty(desiredFreqs);
    // Perturb -delta
    for (int j = 0; j < N; ++j)
        birl::SFXPhysicalModelSetTubeLength(j, L_phys[j] - delta * v_rand[j]);
    double lossMinus = MSEwithPenalty(desiredFreqs);
    // Restore originals
    for (int j = 0; j < N; ++j)
        birl::SFXPhysicalModelSetTubeLength(j, L_phys[j]);

    std::vector<double> grad(N);
    for (int j = 0; j < N; ++j)
        grad[j] = (lossPlus - lossMinus) / (2.0 * delta * v_rand[j]);
    double norm_grad = 0.0;
    for (double g : grad) norm_grad += g * g;
    norm_grad = std::sqrt(norm_grad);
    std::vector<double> v2(N, 0.0);
    if (norm_grad > 0.0) {
        for (int j = 0; j < N; ++j)
            v2[j] = grad[j] / norm_grad;
    }

    // 4) Grid search over alpha, beta ≤ 0
    double bestLoss = std::numeric_limits<double>::infinity();
    double bestAlpha = 0.0, bestBeta = 0.0;
    double alphaMin = -norm_L, alphaMax = 0.0;
    double betaMin  = -norm_L, betaMax  = 0.0;

    for (int ai = 0; ai < gridSize; ++ai) {
        double alpha = alphaMin + (alphaMax - alphaMin) * ai / (gridSize - 1);
        for (int bi = 0; bi < gridSize; ++bi) {
            double beta = betaMin + (betaMax - betaMin) * bi / (gridSize - 1);
            // Apply correction
            for (int j = 0; j < N; ++j) {
                double newL = L_phys[j] + alpha * v1[j] + beta * v2[j];
                birl::SFXPhysicalModelSetTubeLength(j, std::max(newL, 0.0));
            }
            double loss = MSEwithPenalty(desiredFreqs);
            if (loss < bestLoss) {
                bestLoss = loss;
                bestAlpha = alpha;
                bestBeta = beta;
            }
        }
        printf("alpha: %.6f, best loss: %.6f\n", alpha, bestLoss);
    }

    // 5) Apply best found correction
    for (int j = 0; j < N; ++j) {
        double newL = L_phys[j] + bestAlpha * v1[j] + bestBeta * v2[j];
        birl::SFXPhysicalModelSetTubeLength(j, std::max(newL, 0.0));
    }

    printf("Method A applied: alpha=%.6f, beta=%.6f, final loss=%.6f\n", bestAlpha, bestBeta, bestLoss);
    printFreqs();
}

// Example usage:

