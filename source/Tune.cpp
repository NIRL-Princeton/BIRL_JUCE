//
// Created by Milan Sastry on 3/19/25.
//


#include "Yin.h"
#include "sfx.h"
#include "Birl.h"
#include <juce_audio_devices/juce_audio_devices.h>
#include "leaf.h"
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>


// LEAF leaf;
Yin yinz = Yin(44100.0f, 1024, 0.1);
double tubeLengths[NUM_OF_TONEHOLES+1];
//int sampleRate = 48000;

// initializes all the stuff to simulate the model, like the delay lines, toneholes, etc.
void initialize( LEAF leaf){
    // LEAF_init(&leaf,sampleRate, birl::medium_memory, 500000, []() {return (float)rand() / RAND_MAX; });
    // tMempool_init(&birl::smallPool, birl::small_memory, 80328, &leaf);
    // tMempool_init(&birl::largePool, birl::large_memory, 33554432, &leaf);
    // birl::initGlobalSFXObjects(leaf);
    birl::SFXPhysicalModelPMAlloc(leaf);
    // for (int i = 0; i < NUM_OF_TONEHOLES; i++)
    // {
    //     tubeLengths[i] = arr[i];
    //     printf ("length of tube %d: %f\n",i,arr[i]);
    // }

}


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
    float sum =0.0f;
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
    //initialize(48000.0f,leaf);
    // printf ("Running simulation\n");
    static float freqs[NUM_OF_TONEHOLES+1];
    for (int i = 0; i < 10; i++){
        float pitch = simulate(buffer, 6, i,3);
        freqs[i] = pitch;
        //printf ("hole %d: %f \n",i,pitch);
    }
    return freqs;
}

// mean squared error
    float lossFunction(const float* desiredFreqs){
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
    printf("frequencies in Hz\n");
    const float* freqs = getFreqs();
    for (int i = 0; i < 10; i++)
    {
        printf("hole %d: %f\n", i, freqs[i]);
    }
    printf("\n");
}

void printDiffs(const float* desiredFreqs, float* freqs)
{
    printf("differentials in Hz\n");
    for (int i = 0; i < 10; i++)
    {
        float freq = freqs[i];
        float desired = desiredFreqs[i];
        float diff = freq - desired;
        printf("hole %d: %f\n", i, diff);
    }
    printf("\n");
}

float* calculateGradient(const float* desiredFreqs, float epsilon, float currentError){
    float baseError = currentError;
    static float gradient[10];
    static double temp[10];

    for (int i = 0; i < 10; i++)
    {
        double perturbedLengths[10];
        for (int j = 0; j < 10; j++)
        {
            double currentLength = birl::SFXPhysicalModelGetTubeLength (j);
            perturbedLengths[j] = currentLength;
        }
        perturbedLengths[i]+=epsilon;
        for (int j = 0; j < 10; j++)
        {
            temp[j] = birl::SFXPhysicalModelGetTubeLength (j);
            birl::SFXPhysicalModelSetTubeLength(j, perturbedLengths[j]);
        }
        float perturbedError = lossFunction(desiredFreqs);
        for (int j = 0; j < 10; j++)
        {
            birl::SFXPhysicalModelSetTubeLength(j, temp[j]);
        }
        gradient[i] = (perturbedError - baseError) / epsilon;

    }
    return gradient;
}

void gradientDescent(const float* desiredFreqs, float epsilon, float learningRate, int numIterations, float tolerance,LEAF leaf){

    //double* tubeLengths = tubeLengths_;
    initialize (leaf);
    float prevError = std::numeric_limits<float>::max();
   for (int n = 0; n < numIterations; n++){
        float currentError = lossFunction(desiredFreqs);
        if (std::abs(prevError - currentError) < tolerance)
        {
            printf ("converged\n");
            break;
        }
        float* gradient = calculateGradient(desiredFreqs, epsilon, currentError);
       //float* gradient = calculateCentralGradient(desiredFreqs, epsilon);
        for (int i = 0; i < 10; i++)
        {
            double newLength = birl::SFXPhysicalModelGetTubeLength(i) - learningRate * gradient[i];
            birl::SFXPhysicalModelSetTubeLength (i, newLength);

        }
        prevError = currentError;
        if (n % 10 == 0)
        {
            printf ("Error: %f Iteration %d\n", currentError,n);
            for (int i = 0; i < 10; i++)
            {
                double length = birl::SFXPhysicalModelGetTubeLength(i);
                printf ("lengths %d: %f \n",i,length);
            }
            for (int i =0; i < 10; i++)
            {
                printf("gradient[%d]: %f \n",i,gradient[i]);
            }
            printf ("\n");
        }


    }
    for (int i = 0; i < 10; i++)
    {
        printf ("lengths %d: %f \n",i,tubeLengths[i]);
    }
}


void spsaGradientDescent(const float* desiredFreqs,
                         float epsilon,
                         float learningRate,
                         int numIterations,
                         float tolerance,
                         LEAF leaf)
{
    // Initialize the physical model and tube lengths
    //initialize(leaf);
    float prevError = std::numeric_limits<float>::max();

    // Main optimization loop
    for (int n = 0; n < numIterations; n++)
    {
        // Evaluate the current loss
        float currentError = lossFunction(desiredFreqs);
        if (std::abs(prevError - currentError) < tolerance)
        {
            printf("Converged at iteration %d with loss: %f\n", n, currentError);
            break;
        }

        // Generate a random perturbation vector with entries ±1
        double delta[10];
        for (int i = 0; i < 10; i++)
        {
            delta[i] = (rand() % 2 == 0) ? 1.0 : -1.0;
        }

        // Save the original tube lengths
        double originalLengths[10];
        for (int i = 0; i < 10; i++)
        {
            originalLengths[i] = birl::SFXPhysicalModelGetTubeLength(i);
        }

        // Compute loss at the positive perturbation: tubeLengths + epsilon * delta
        for (int i = 0; i < 10; i++)
        {
            double newLength = originalLengths[i] + epsilon * delta[i];
            birl::SFXPhysicalModelSetTubeLength(i, newLength);
        }
        float lossPlus = lossFunction(desiredFreqs);

        // Compute loss at the negative perturbation: tubeLengths - epsilon * delta
        for (int i = 0; i < 10; i++)
        {
            double newLength = originalLengths[i] - epsilon * delta[i];
            birl::SFXPhysicalModelSetTubeLength(i, newLength);
        }
        float lossMinus = lossFunction(desiredFreqs);

        // Compute the SPSA gradient estimate for each tube length:
        float gradient[10];
        for (int i = 0; i < 10; i++)
        {
            gradient[i] = (lossPlus - lossMinus) / (2.0f * epsilon * delta[i]);
        }

        // Restore original tube lengths before applying the update
        for (int i = 0; i < 10; i++)
        {
            birl::SFXPhysicalModelSetTubeLength(i, originalLengths[i]);
        }

        for (int i = 0; i < 10; i++)
        {
            double newLength = originalLengths[i] - learningRate * gradient[i];
            if(newLength < 0.0)
                newLength = 0.0;
            tubeLengths[i] = newLength;
            birl::SFXPhysicalModelSetTubeLength(i, newLength);
        }

        prevError = currentError;

        if (true)
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
}

void spsaGradientDescentMom(const float* desiredFreqs,
                         float epsilon,
                         float learningRate,
                         int numIterations,
                         float tolerance,
                         float momentum,      // e.g., 0.9
                         LEAF leaf)
{
    // Initialize previous error and velocity for each tube.
    float prevError = std::numeric_limits<float>::max();
    double velocity[10] = {0.0};

    // Main optimization loop.
    for (int n = 0; n < numIterations; n++)
    {
        float currentError = lossFunction(desiredFreqs);
        if (std::abs(prevError - currentError) < tolerance)
        {
            printf("Converged at iteration %d with loss: %f\n", n, currentError);
            break;
        }

        // Get current tube lengths.
        double params[10];
        for (int i = 0; i < 10; i++)
        {
            params[i] = birl::SFXPhysicalModelGetTubeLength(i);
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
            double newParam = params[i] + epsilon * delta[i];
            birl::SFXPhysicalModelSetTubeLength(i, newParam);
        }
        float lossPlus = lossFunction(desiredFreqs);

        // Compute loss for negative perturbation.
        for (int i = 0; i < 10; i++)
        {
            double newParam = params[i] - epsilon * delta[i];
            birl::SFXPhysicalModelSetTubeLength(i, newParam);
        }
        float lossMinus = lossFunction(desiredFreqs);

        // Restore original parameters.
        for (int i = 0; i < 10; i++)
        {
            birl::SFXPhysicalModelSetTubeLength(i, params[i]);
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

            // Update parameter using the momentum term.
            double newParam = params[i] - learningRate * velocity[i];
            if(newParam < 0.0)
                newParam = 0.0;  // enforce non-negative tube length

            // Update global tube length if applicable and set the new tube length.
            tubeLengths[i] = newParam;
            birl::SFXPhysicalModelSetTubeLength(i, newParam);
        }

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


// Multi-sample SPSA with momentum: averages over numSamples perturbation directions,
// then uses momentum to update parameters.
void multiSampleSPSAWithMomentum(const float* desiredFreqs,
                                 float epsilon,
                                 float initialLearningRate,
                                 float decayFactor,
                                 int numIterations,
                                 float tolerance,
                                 int numSamples,
                                 float momentum,  // e.g., 0.9
                                 LEAF leaf)
{
    float prevError = std::numeric_limits<float>::max();

    // Initialize velocity vector for momentum.
    double velocity[10] = {0.0};

    for (int iter = 0; iter < numIterations; iter++)
    {
        float learningRate = initialLearningRate;



        float currentError = lossFunction(desiredFreqs);
        if (std::abs(prevError - currentError) < tolerance)
        {
            printf("Converged at iteration %d with loss: %f\n", iter, currentError);
            break;
        }

        // Get current parameters.
        double params[10];
        for (int i = 0; i < 10; i++)
        {
            params[i] = birl::SFXPhysicalModelGetTubeLength(i);
        }

        // Initialize gradient accumulator.
        float gradientSum[10] = {0.0f};
        int sampleNumber = numSamples;

        // if (iter > numIterations/8)
        // {
        //     sampleNumber = numSamples+2;
        // }
        // if (iter > numIterations/4)
        // {
        //     sampleNumber = numSamples+3;
        // }
        // if (iter > numIterations/2)
        // {
        //     sampleNumber = numSamples+4;
        // }


        // Loop over multiple SPSA samples.
        for (int sample = 0; sample < sampleNumber; sample++)
        {
            // Generate a random perturbation vector (delta) with entries ±1.
            double delta[10];
            for (int i = 0; i < 10; i++)
            {
                delta[i] = (rand() % 2 == 0) ? 1.0 : -1.0;
            }

            // Evaluate loss for positive perturbation.
            for (int i = 0; i < 10; i++)
            {
                double newParam = params[i] + epsilon * delta[i];
                birl::SFXPhysicalModelSetTubeLength(i, newParam);
            }
            float lossPlus = lossFunction(desiredFreqs);

            // Evaluate loss for negative perturbation.
            for (int i = 0; i < 10; i++)
            {
                double newParam = params[i] - epsilon * delta[i];
                birl::SFXPhysicalModelSetTubeLength(i, newParam);
            }
            float lossMinus = lossFunction(desiredFreqs);

            // Restore original parameters.
            for (int i = 0; i < 10; i++)
            {
                birl::SFXPhysicalModelSetTubeLength(i, params[i]);
            }

            // Compute and accumulate the gradient estimate for this sample.
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

        // Update parameters using momentum.
        for (int i = 0; i < 10; i++)
        {
            // Update momentum: velocity = momentum * velocity + avgGradient.
            velocity[i] = momentum * velocity[i] + avgGradient[i];

            // Update parameter using the momentum term.
            double newParam = params[i] - learningRate * velocity[i];
            if (newParam < 0.0)
                newParam = 0.0;  // enforce non-negative tube length

            tubeLengths[i] = newParam; // update global tube length if applicable.
            birl::SFXPhysicalModelSetTubeLength(i, newParam);
        }

        prevError = currentError;

        if (iter % 10 == 0)
        {
            printf("Iteration %d, Error: %f, Learning Rate: %f\n", iter, currentError,learningRate);
            for (int i = 0; i < 10; i++)
            {
                double len = birl::SFXPhysicalModelGetTubeLength(i);
                printf("Tube %d length: %f, avg gradient: %f, velocity: %f\n", i, len, avgGradient[i], velocity[i]);
            }
        }
    }

    for (int i = 0; i < 10; i++)
    {
        printf("Final tube length %d: %f\n", i, tubeLengths[i]);
    }
    printFreqs();
}


// hole 0: 509.600769
// hole 1: 449.589996
// hole 2: 404.052002
// hole 3: 374.100647
// hole 4: 334.878479
// hole 5: 303.636108
// hole 6: 273.128235
// hole 7: 256.884705
// hole 8: 232.439407
// hole 9: 212.101227
//


// gradientDescent (targetFrequencies, 1e-4, 5e-7,1000, 1e-6, leaf);
// Error: 60.541729 Iteration 999
// lengths 0: 16.106490
// lengths 1: 2.738496
// lengths 2: 2.584404
// lengths 3: 1.116656
// lengths 4: 2.984964
// lengths 5: 2.760160
// lengths 6: 3.789892
// lengths 7: 1.727938
// lengths 8: 3.810752
// lengths 9: 5.030512


//spsaGradientDescent (targetFrequencies, 1e-4, 8e-7, numIterations, tolerance, leaf);
// Final tube length 0: 16.139052
// Final tube length 1: 2.759116
// Final tube length 2: 2.599197
// Final tube length 3: 1.115447
// Final tube length 4: 2.982172
// Final tube length 5: 2.772683
// Final tube length 6: 3.783491
// Final tube length 7: 1.716655
// Final tube length 8: 3.812374
// Final tube length 9: 5.025982

//multiSampleSPSAWithMomentum (targetFrequencies, 0.0002, 0.00001,0.99,2000, 1e-6, 2,0.95, leaf);
// halving learning rate after 500 iterations

// Iteration 1210, Error: 8.924214, Learning Rate: 0.000002
// Tube 0 length: 16.368984, avg gradient: -91.303589, velocity: -289.427818
// Tube 1 length: 3.076642, avg gradient: 38.453346, velocity: 60.042284
// Tube 2 length: 2.446632, avg gradient: 5.279783, velocity: 13.585908
// Tube 3 length: 1.222204, avg gradient: 61.603783, velocity: 26.969868
// Tube 4 length: 2.490758, avg gradient: 99.776985, velocity: 69.110650
// Tube 5 length: 3.059859, avg gradient: -34.979591, velocity: 253.465919
// Tube 6 length: 3.976828, avg gradient: 91.303589, velocity: 296.947733
// Tube 7 length: 1.886947, avg gradient: -91.303589, velocity: 328.601855
// Tube 8 length: 3.835278, avg gradient: -8.753538, velocity: -51.598539
// Tube 9 length: 5.380689, avg gradient: -0.280142, velocity: 167.906314
//
//
//
// hole 0: 640.687866
// hole 1: 554.488770
// hole 2: 498.921417
// hole 3: 458.898712
// hole 4: 416.222198
// hole 5: 376.047272
// hole 6: 336.898346
// hole 7: 316.362885
// hole 8: 286.682098
// hole 9: 261.571289


//Iteration 5050, Error: 1.525883, 200hz fundamental
// Tube 0 length: 21.956817, gradient: -4.733205
// Tube 1 length: 3.737189, gradient: 4.733205
// Tube 2 length: 3.480582, gradient: 4.733205
// Tube 3 length: 1.709943, gradient: -4.733205
// Tube 4 length: 3.984575, gradient: 4.733205
// Tube 5 length: 3.974303, gradient: 4.733205
// Tube 6 length: 5.280642, gradient: -4.733205
// Tube 7 length: 2.557904, gradient: 4.733205
// Tube 8 length: 5.225859, gradient: -4.733205
// Tube 9 length: 6.601649, gradient: -4.733205

// local optimization, doesn't really work
void bfgsGradientDescent(const float* desiredFreqs, float epsilon, int maxIterations, float tolerance)
{
    const int n = 10;
// Get current tube lengths into vector x
double x[n];
for (int i = 0; i < n; i++) {
    x[i] = birl::SFXPhysicalModelGetTubeLength(i);
}

// Initialize the inverse Hessian approximation H as identity
double H[n][n];
for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
        H[i][j] = (i == j) ? 1.0 : 0.0;
    }
}

// Evaluate the initial loss and gradient
double f_val = lossFunction(desiredFreqs);
double g[n];
{
    float* grad = calculateGradient(desiredFreqs, epsilon, f_val);
    for (int i = 0; i < n; i++)
        g[i] = grad[i];
}

// Main BFGS iterations
int iter;
for (iter = 0; iter < maxIterations; iter++) {
    // Compute the norm of the gradient
    double norm_g = 0.0;
    for (int i = 0; i < n; i++)
        norm_g += g[i] * g[i];
    norm_g = sqrt(norm_g);
    if (norm_g < tolerance) {
        printf("BFGS converged at iteration %d with norm(g) = %f\n", iter, norm_g);
        break;
    }

    // Compute search direction p = -H * g
    double p[n];
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            sum += H[i][j] * g[j];
        }
        p[i] = -sum;
    }

    // Backtracking line search parameters
    double alpha = 1.0;
    double c = 1e-4;
    // Save current x as x_old
    double x_old[n];
    for (int i = 0; i < n; i++)
        x_old[i] = x[i];

    // Compute dot product g^T * p
    double dot_gp = 0.0;
    for (int i = 0; i < n; i++)
        dot_gp += g[i] * p[i];

    // Lambda to set tube lengths from an array
    auto set_x = [&](double arr[]) {
        for (int i = 0; i < n; i++) {
            // Ensure nonnegative tube lengths.
            double newVal = (arr[i] < 0.0) ? 0.0 : arr[i];
            birl::SFXPhysicalModelSetTubeLength(i, newVal);
        }
    };

    // Evaluate the loss at the current x
    set_x(x_old);
    double f_old = lossFunction(desiredFreqs);

    // Perform backtracking line search
    double x_new[n];
    while (true) {
        for (int i = 0; i < n; i++) {
            x_new[i] = x_old[i] + alpha * p[i];
            if (x_new[i] < 0.0)
                x_new[i] = 0.0;
        }
        set_x(x_new);
        double f_new = lossFunction(desiredFreqs);
        if (f_new <= f_old + c * alpha * dot_gp)
            break;
        alpha *= 0.5;
        if (alpha < 1e-8)
            break;
    }

    // Update x: x_new becomes the new x.
    for (int i = 0; i < n; i++)
        x[i] = x_new[i];
    set_x(x);

    // Compute new gradient g_new at updated x
    double g_new[n];
    {
        float* grad_new = calculateGradient(desiredFreqs, epsilon, lossFunction(desiredFreqs));
        for (int i = 0; i < n; i++)
            g_new[i] = grad_new[i];
    }

    // Compute differences: s = x_new - x_old, y = g_new - g
    double s[n], y[n];
    for (int i = 0; i < n; i++) {
        s[i] = x[i] - x_old[i];
        y[i] = g_new[i] - g[i];
    }

    // Compute rho = 1 / (y^T * s)
    double ys = 0.0;
    for (int i = 0; i < n; i++)
        ys += y[i] * s[i];
    double rho = (ys != 0.0) ? 1.0 / ys : 1e+10;

    // Update the inverse Hessian approximation H using the BFGS formula:
    // H_new = (I - rho * s * y^T) * H * (I - rho * y * s^T) + rho * s * s^T
    double delta[n][n];
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            delta[i][j] = ((i == j) ? 1.0 : 0.0) - rho * s[i] * y[j];
        }
    }
    double M[n][n];
    // Compute M = delta * H
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            double sum = 0.0;
            for (int k = 0; k < n; k++){
                sum += delta[i][k] * H[k][j];
            }
            M[i][j] = sum;
        }
    }
    double term1[n][n];
    // Compute term1 = M * delta^T (note: delta^T[i][j] = delta[j][i])
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            double sum = 0.0;
            for (int k = 0; k < n; k++){
                sum += M[i][k] * delta[j][k];
            }
            term1[i][j] = sum;
        }
    }
    double term2[n][n];
    // Compute term2 = rho * s * s^T
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            term2[i][j] = rho * s[i] * s[j];
        }
    }
    double H_new[n][n];
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            H_new[i][j] = term1[i][j] + term2[i][j];
        }
    }
    // Copy H_new into H for the next iteration.
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            H[i][j] = H_new[i][j];
        }
    }
    // Update the gradient vector for next iteration.
    for (int i = 0; i < n; i++){
        g[i] = g_new[i];
    }

    printf("BFGS Iteration %d, Loss: %f, Norm(g): %f, alpha: %f\n", iter, lossFunction(desiredFreqs), norm_g, alpha);
    for (int i = 0; i < n; i++){
        printf("Tube %d length: %f\n", i, x[i]);
    }
}

printf("Final tube lengths from BFGS:\n");
for (int i = 0; i < n; i++){
     printf("Tube %d length: %f\n", i, x[i]);
}

}

const int PARAM_DIM = 10;
std::default_random_engine generator;

// === Sample Candidates ===
std::vector<std::vector<double>> samplePopulation(
    const std::vector<double>& mean,
    const std::vector<std::vector<double>>& cov,
    double sigma,
    int lambda)
{
    std::vector<std::vector<double>> population(lambda, std::vector<double>(PARAM_DIM));
    std::normal_distribution<double> dist(0.0, 1.0);

    for (int i = 0; i < lambda; i++) {
        for (int j = 0; j < PARAM_DIM; j++) {
            double z = dist(generator);
            population[i][j] = mean[j] + sigma * std::sqrt(cov[j][j]) * z;
        }
    }
    return population;
}

// SMA-ES genetic algorithm
std::vector<std::pair<std::vector<double>, double>> evaluateFitness(
    const std::vector<std::vector<double>>& population,
    const float* desiredFreqs)
{
    std::vector<std::pair<std::vector<double>, double>> evaluated;
    for (auto& individual : population) {
        for (int i = 0; i < PARAM_DIM; i++) {
            birl::SFXPhysicalModelSetTubeLength(i, individual[i]);
        }
        float loss = lossFunction(desiredFreqs);
        evaluated.push_back({ individual, loss });
    }
    std::sort(evaluated.begin(), evaluated.end(), [](auto& a, auto& b) {
        return a.second < b.second;
    });
    return evaluated;
}

std::vector<double> updateMean(
    const std::vector<std::pair<std::vector<double>, double>>& sortedPop,
    int mu)
{
    std::vector<double> newMean(PARAM_DIM, 0.0);
    double weight = 1.0 / mu;

    for (int i = 0; i < mu; i++) {
        for (int j = 0; j < PARAM_DIM; j++) {
            newMean[j] += weight * sortedPop[i].first[j];
        }
    }
    return newMean;
}

void cmaesTuning(const float* desiredFreqs)
{
    int lambda = 20;
    int mu = 10;
    double sigma = 0.5;
    int iterations = 100;
    double tubeLengths[NUM_OF_TONEHOLES+1];

    std::vector<double> mean(PARAM_DIM);
    for (int i = 0; i < PARAM_DIM; i++) {
        mean[i] = birl::SFXPhysicalModelGetTubeLength(i);
    }

    std::vector<std::vector<double>> cov(PARAM_DIM, std::vector<double>(PARAM_DIM, 0.0));
    for (int i = 0; i < PARAM_DIM; i++) cov[i][i] = 1.0;

    for (int gen = 0; gen < iterations; gen++) {
        auto population = samplePopulation(mean, cov, sigma, lambda);
        auto evaluated = evaluateFitness(population, desiredFreqs);
        mean = updateMean(evaluated, mu);

        for (int i = 0; i < PARAM_DIM; i++) {
            birl::SFXPhysicalModelSetTubeLength(i, mean[i]);
        }

        printf("Generation %d, Best Loss: %f\n", gen, evaluated[0].second);
        for (int i = 0; i < PARAM_DIM; i++) {
            printf("Tube %d: %f\n", i, mean[i]);
        }
    }

    for (int i = 0; i < PARAM_DIM; i++) {
        tubeLengths[i] = mean[i];
    }
}