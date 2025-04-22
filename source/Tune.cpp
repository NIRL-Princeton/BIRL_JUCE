//
// Created by Milan Sastry on 3/19/25.
//


#include "Yin.h"
#include "sfx.h"
#include "Birl.h"
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>


Yin yinz = Yin(44100.0f, 1024, 0.1);
double tubeLengths[NUM_OF_TONEHOLES+1];

void resetFingers(){
    for (int i = 0; i < NUM_OF_TONEHOLES; i++)
    {
        birl::fingers[i] = 0.0f;
    }
}

// runs the physical model for the shortest amount of time that will produce a buffer that can be used to find frequency for a tonehole
float simulate(juce::AudioBuffer<float>& buffer, int iterations, int numHolesClosed, int numAveraged){
    resetFingers();
    std::vector<float> freqs(iterations);

    birl::breathArray[0] = 1.0f;
    birl::breathArray[1] = 1.0f;
    birl::SFXPhysicalModelSetBreathPressure(0.75f);

    for (int i = 0; i < numHolesClosed; i++)
    {
        birl::fingers[i] = 1.0f;
    }

    float outputSample[2] = { 0.0f, 0.0f };
    buffer.clear();
    for(int i = 0; i < iterations; i++){
        // buffer.clear (0, 0, buffer.getNumSamples());
        // buffer.clear (1, 0, buffer.getNumSamples());
        float* leftChannel = buffer.getWritePointer(0);
        float* rightChannel = buffer.getWritePointer(1);

        for(int j = 0; j < buffer.getNumSamples(); j++){
            outputSample[0] = leftChannel[j];
            outputSample[1] = rightChannel[j];

            birl::SFXPhysicalModelPMTick(outputSample);

            buffer.setSample(0, j, outputSample[0]);
            buffer.setSample(1, j, outputSample[1]);
        }
        float pitch = yinz.getPitch (buffer);
        if (pitch > 0.0f) freqs[i] = pitch;
        // if (pitch > 0.0f) DBG("Detected Pitch: " << pitch << " Hz");
    }
    float sum = 0.0f;
    for (int i = 0; i < numAveraged; i++)
    {
        sum+=freqs[iterations - numAveraged + i];
    }
    for (int i = 0; i < numHolesClosed; i++)
    {
        birl::fingers[i] = 0.0f;
    }
    birl::SFXPhysicalModelSetBreathPressure(0.75f);
    birl::breathArray[0] = 0.0f;
    birl::breathArray[1] = 0.0f;

    return sum / static_cast<float>(numAveraged);
}
float* getFreqs(){
    juce::AudioBuffer<float> buffer(2, 1024);
    static float freqs[NUM_OF_TONEHOLES+1];
    for (int i = 0; i < 10; i++){
        float pitch = simulate(buffer, 6, i,3);
        freqs[i] = LEAF_frequencyToMidi (pitch);
        //printf ("hole %d: %f \n",i,pitch);
    }
    return freqs;
}

// mean squared error
    float MSEwithPenalty(const float* desiredFreqs){
        float* freqs = getFreqs();
        float errors[10];
        float sum = 0.0f;
        float maxErr = 0.0f;

        for (int i = 0; i < 10; i++){
            errors[i] = fabs(freqs[i] - desiredFreqs[i]);
            sum += errors[i];
            if (errors[i] > maxErr) maxErr = errors[i];
        }

        float mean = sum / 10.0f;

        float variance = 0.0f;
        for (int i = 0; i < 10; i++){
            float diff = errors[i] - mean;
            variance += diff * diff;
        }
        float stddev = sqrt(variance / 10.0f);

        // Weighted loss: MSE-like base + spread penalty + max outlier penalty
        float mse = 0.0f;
        for (int i = 0; i < 10; i++){
            mse += errors[i] * errors[i];
        }
        mse /= 10.0f;

        float loss = 0.6f * mse + 0.25f * stddev + 0.15f * maxErr;
        return loss;
    }

void printFreqs()
{
    printf("pitches in cents\n");
    const float* freqs = getFreqs();
    for (int i = 0; i < 10; i++)
    {
        printf("hole %d: %f\n", i, freqs[i]);
    }
    printf("\n");
}

void printDiffs(const float* desiredFreqs, float* freqs)
{
    printf("differentials in cents\n");
    for (int i = 0; i < 10; i++)
    {
        float freq = freqs[i];
        float desired = desiredFreqs[i];
        float diff = freq - desired;
        printf("hole %d: %f\n", i, diff);
    }
    printf("\n");
}

void spsaGradientDescent(const float* desiredFreqs,
                         float epsilon,
                         float learningRate,
                         int numIterations,
                         float tolerance,
                         float momentum)
{
    // Initialize previous error and velocity for each tube.
    float prevError = std::numeric_limits<float>::max();
    double velocity[10] = {0.0};

    // Main optimization loop.
    for (int n = 0; n < numIterations; n++)
    {
        // float currentError = MSEwithPenalty(desiredFreqs);
        // if (std::abs(prevError - currentError) < tolerance)
        // {
        //     printf("Converged at iteration %d with loss: %f\n", n, currentError);
        //     break;
        // }

        // Get current tube lengths.
        double initialLengths[10];
        for (int i = 0; i < 10; i++)
        {
            initialLengths[i] = birl::SFXPhysicalModelGetTubeLength(i);
        }

        // Generate a random perturbation vector delta with entries ±1.
        double delta[10];
        for (int i = 0; i < 10; i++)
        {
            delta[i] = (rand() % 2 == 0) ? 1.0 : -1.0;
        }

        // Compute loss for positive perturbation.
        for (int i = 0; i < 10; i++)
        {
            double perturbedLengths = initialLengths[i] + epsilon * delta[i];
            birl::SFXPhysicalModelSetTubeLength(i, perturbedLengths);
        }
        float lossPlus = MSEwithPenalty(desiredFreqs);

        // Compute loss for negative perturbation.
        for (int i = 0; i < 10; i++)
        {
            double perturbedLengths = initialLengths[i] - epsilon * delta[i];
            birl::SFXPhysicalModelSetTubeLength(i, perturbedLengths);
        }
        float lossMinus = MSEwithPenalty(desiredFreqs);

        // restore initial lengths.
        for (int i = 0; i < 10; i++)
        {
            birl::SFXPhysicalModelSetTubeLength(i, initialLengths[i]);
        }

        // Compute SPSA gradient estimate for each tube.
        float gradient[10];
        for (int i = 0; i < 10; i++)
        {
            gradient[i] = (lossPlus - lossMinus) / (2.0f * epsilon * delta[i]);
        }

        // Update momentum (velocity) and then update tube lengths.
        for (int i = 0; i < 10; i++)
        {
            // Update velocity: combine previous velocity with the current gradient.
            velocity[i] = momentum * velocity[i] + gradient[i];

            // Update length using the momentum term.
            double newLength = initialLengths[i] - learningRate * velocity[i];
            if(newLength < 0.0)
                newLength = 0.0;

            tubeLengths[i] = newLength;
            birl::SFXPhysicalModelSetTubeLength(i, newLength);
        }

        float currentError = 0.5f*(lossPlus + lossMinus);
        if (fabs(prevError - currentError) < tolerance) break;
        prevError = currentError;

        if (n % 10 == 0)
        {
            printf("Iteration %d, Error: %f\n", n, currentError);
            for (int i = 0; i < 10; i++)
            {
                double len = birl::SFXPhysicalModelGetTubeLength(i);
                printf("Tube %d length: %f, gradient: %f\n", i, len, gradient[i]);
            }
        }
    }

    for (int i = 0; i < 10; i++)
    {
        printf("Final tube length %d: %f\n", i, tubeLengths[i]);
    }
    printFreqs();
}


// Multi-sample SPSA: averages over numSamples perturbation directions
void multiSampleSPSA(const float* desiredFreqs,
                                 float epsilon,
                                 float learningRate,
                                 float decayFactor,
                                 int numIterations,
                                 float tolerance,
                                 int numSamples,
                                 float momentum){
    float prevError = std::numeric_limits<float>::max();

    // Initialize velocity vector for momentum.
    double velocity[10] = {0.0};

    for (int iter = 0; iter < numIterations; iter++)
    {
        float currentError = MSEwithPenalty(desiredFreqs);

        if (std::abs(prevError - currentError) < tolerance)
        {
            printf("Converged at iteration %d with loss: %f\n", iter, currentError);
            break;
        }

        double initialLengths[10];
        for (int i = 0; i < 10; i++)
        {
            initialLengths[i] = birl::SFXPhysicalModelGetTubeLength(i);
        }

        float gradientSum[10] = {0.0f};
        int sampleNumber = numSamples;

        for (int sample = 0; sample < sampleNumber; sample++)
        {
            double delta[10];
            for (int i = 0; i < 10; i++)
            {
                delta[i] = (rand() % 2 == 0) ? 1.0 : -1.0;
            }

            for (int i = 0; i < 10; i++)
            {
                double newLength = initialLengths[i] + epsilon * delta[i];
                birl::SFXPhysicalModelSetTubeLength(i, newLength);
            }
            float lossPlus = MSEwithPenalty(desiredFreqs);

            for (int i = 0; i < 10; i++)
            {
                double newLength = initialLengths[i] - epsilon * delta[i];
                birl::SFXPhysicalModelSetTubeLength(i, newLength);
            }
            float lossMinus = MSEwithPenalty(desiredFreqs);

            for (int i = 0; i < 10; i++)
            {
                birl::SFXPhysicalModelSetTubeLength(i, initialLengths[i]);
            }

            for (int i = 0; i < 10; i++)
            {
                float sampleGradient = (lossPlus - lossMinus) / (2.0f * epsilon * delta[i]);
                gradientSum[i] += sampleGradient;
            }
        }

        // Average the gradient over the samples.
        float avgGradient[10];
        for (int i = 0; i < 10; i++)
        {
            avgGradient[i] = gradientSum[i] / numSamples;
        }

        for (int i = 0; i < 10; i++)
        {
            velocity[i] = momentum * velocity[i] + avgGradient[i];

            double newLength = initialLengths[i] - learningRate * velocity[i];
            if (newLength < 0.0)
                newLength = 0.0;

            tubeLengths[i] = newLength;
            birl::SFXPhysicalModelSetTubeLength(i, newLength);
        }

        prevError = currentError;

        if (iter % 10 == 0)
        {
            printf("Iteration %d, Error: %f, Learning Rate: %f\n", iter, currentError,learningRate);
            for (int i = 0; i < 10; i++)
            {
                double len = birl::SFXPhysicalModelGetTubeLength(i);
                printf("Tube %d length: %f, avg gradient: %f", i, len, avgGradient[i]);
            }
        }
    }

    for (int i = 0; i < 10; i++)
    {
        printf("Final tube length %d: %f\n", i, tubeLengths[i]);
    }
    printFreqs();
}

