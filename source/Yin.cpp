//
// Created by Milan Sastry on 3/15/25.
//

#include "Yin.h"
#include <juce_audio_devices/juce_audio_devices.h>

Yin::Yin(const float rate, const int size, const double thresh)
    : sampleRate(rate), bufferSize(size), threshold(thresh)
{
    yinBuffer.resize(bufferSize/2, 0.0f);
}

float Yin::getPitch(const juce::AudioBuffer<float>& buffer)
{
    float pitch;

    difference(buffer);
    normalizeDifference();
    int tauEstimate = absoluteThreshold();

    if (tauEstimate != -1)
    {
        float betterTau = parabolicInterpolation (tauEstimate);
        pitch = sampleRate / betterTau;
    }
    else pitch = -1.0f;
    return pitch;
}

void Yin::difference (const juce::AudioBuffer<float>& buffer)
{
    int tau;
    for (tau = 0; tau < bufferSize/2; tau++) {
        yinBuffer[tau] = 0;
    }
    for (tau = 1; tau < bufferSize/2; tau++) {
        for (int index = 0; index < bufferSize/2; index++) {
            const float delta = buffer.getSample (channel, index) - buffer.getSample (channel, index + tau);
            yinBuffer[tau] += delta * delta;
        }
    }
}

void Yin::normalizeDifference()
{
    yinBuffer[0] = 1.0f;
    float sum = 0.0f;
    for (int tau = 1; tau < bufferSize/2; tau++)
    {
        sum += yinBuffer[tau];
        yinBuffer[tau] *= yinBuffer[tau] / sum;
    }
}

int Yin::absoluteThreshold()
{
    int tau;
    for (tau = 2; tau < bufferSize/2; tau++)
    {
        if (yinBuffer[tau] < threshold)
        {
            while (tau + 1 < bufferSize/2 && yinBuffer[tau + 1] < yinBuffer[tau]) tau++;
            break;
        }
    }
    if (tau == bufferSize/2 || yinBuffer[tau] >= threshold) tau = -1;
    return tau;
}

float Yin::parabolicInterpolation(int tauEstimate)
{
    float betterTau;
    int x0;
    int x2;

    if (tauEstimate < 1) x0 = tauEstimate;
    else x0 = tauEstimate - 1;

    if (tauEstimate + 1 < bufferSize/2) x2 = tauEstimate + 1;
    else x2 = tauEstimate;

    if (x0 == tauEstimate)
    {
        if (yinBuffer[tauEstimate] <= yinBuffer[x2]) betterTau = tauEstimate;
        else betterTau =  x2;
    }
    else if (x2 == tauEstimate)
    {
        if (yinBuffer[tauEstimate] <= yinBuffer[x0]) betterTau = tauEstimate;
        else betterTau = x0;
    }
    else
    {
        const float s0 = yinBuffer[x0];
        const float s1 = yinBuffer[tauEstimate];
        const float s2 = yinBuffer[x2];
        betterTau = tauEstimate + (s2 - s0) / (2 * (2 * s1 - s2 - s0));
    }
    return betterTau;
}








