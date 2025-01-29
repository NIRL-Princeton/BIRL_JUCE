#ifndef PROCESSOR_H_
#define PROCESSOR_H_


#include "Filters.h"
#include "Tube.h"
#include <juce_audio_processors/juce_audio_processors.h>
#include <juce_osc/juce_osc.h>
#include <juce_audio_devices/juce_audio_devices.h>
#include "melatonin_audio_sparklines/melatonin_audio_sparklines.h"
#include "leaf-envelopes.h"
#include "leaf-oscillators.h"
#include "PhysicalModel.h"
#include "leaf.h"
using namespace juce;
#define MEM_SIZE 50000
#define OPENHOLE_FREQ 349.0

class BirlOneHoleAudioProcessor  : public AudioProcessor, public PhysicalModel
{
public:
    BirlOneHoleAudioProcessor();
    ~BirlOneHoleAudioProcessor();

    double gain;
    //    double frequency;
    //    double filterfrequency;
    AudioProcessorValueTreeState tree;
    int keys_[NUM_OF_KEYS];
    bool buttons_[NUM_OF_BUTTONS];
    int maxKeyArg[NUM_OF_KEYS];
    int minKeyArg[NUM_OF_KEYS];

    double keysGradient [NUM_OF_KEYS];


    void prepareToPlay (double sampleRate, int samplesPerBlock) override;
    void releaseResources() override;
#ifndef JucePlugin_PreferredChannelConfigurations
    bool isBusesLayoutSupported (const BusesLayout& layouts) const override;
#endif
    void processBlock (AudioBuffer<float>&, MidiBuffer&) override;

    AudioProcessorEditor* createEditor() override;
    bool hasEditor() const override;

    const String getName() const override;
    bool acceptsMidi() const override;
    bool producesMidi() const override;
    bool isMidiEffect() const override;
    double getTailLengthSeconds() const override;

    int getNumPrograms() override;
    int getCurrentProgram() override;
    void setCurrentProgram (int index) override;
    const String getProgramName (int index) override;
    void changeProgramName (int index, const String& newName) override;

    //==============================================================================
    void getStateInformation (MemoryBlock& destData) override;
    void setStateInformation (const void* data, int sizeInBytes) override;
    //==============================================================================
    void synthesize();
    void parameterSmooth();
    void scaleValues();
    void changeTuning();
    void buttonListen();
    //==============================================================================
    tSVF filter1{};
    tSVF filter2{};
    float modAmount{10};

private:
    PhysicalModel birl;
    /* Envelopes */
    tADSR breathPosEnv;

    tExpSmooth XY1_XEnv;
    tExpSmooth XY1_YEnv;
    tExpSmooth XY2_XEnv;
    tExpSmooth XY2_YEnv;
    tExpSmooth keysEnv[11];


    //    tSawtooth saw1;
    //    tSawtooth saw2;
    //    tSawtooth saw3;

    tCycle sin{};
    tCycle mod{};

    tExpSmooth smoother{};

    // detune array should take MIDI
    float detune[3] = {0.999, 1.0, 1.001};



    char myMemory[MEM_SIZE];
    float audioBufferTick{};

    //freq modulation
    float modTick;

    float smoothGain;
    LEAF leaf;
    //==============================================================================
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (BirlOneHoleAudioProcessor)
};


#endif
