//
// Created by Milan Sastry on 3/15/25.
//

#ifndef YIN_H
#define YIN_H

#include <juce_audio_devices/juce_audio_devices.h>
#include <vector>
constexpr int channel = 0; //use the left channel

class Yin {
    public:
        Yin(float rate, int size, double thresh);
        float getPitch(const juce::AudioBuffer<float>& buffer);

    private:
        float sampleRate;
        int bufferSize;
        double threshold;
        std::vector<float> yinBuffer;

        void difference(const juce::AudioBuffer<float>& buffer);
        void normalizeDifference();
        int absoluteThreshold();
        float parabolicInterpolation(int tauEstimate);
};

#endif //YIN_H
