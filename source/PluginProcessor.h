#ifndef PROCESSOR_H_
#define PROCESSOR_H_

#include "sfx.h"
#include "Filters.h"
#include "Tube.h"
#include <juce_audio_processors/juce_audio_processors.h>
#include <juce_osc/juce_osc.h>
#include <juce_audio_devices/juce_audio_devices.h>
#include "leaf.h"
//==============================================================================
/**
*/
class BirlAudioProcessor  : public juce::AudioProcessor, private juce::OSCReceiver::Listener<juce::OSCReceiver::MessageLoopCallback>, public juce::AsyncUpdater, public juce::MidiInputCallback, public juce::MidiKeyboardStateListener
{
public:
    //==============================================================================
    BirlAudioProcessor();
    ~BirlAudioProcessor() override;
    //==============================================================================
    juce::MidiKeyboardState keyboardState;
    void handleIncomingMidiMessage(juce::MidiInput* source, const juce::MidiMessage& message) override;
    void handleNoteOn(juce::MidiKeyboardState*, int midiChannel, int midiNoteNumber, float velocity) override;
    void handleNoteOff(juce::MidiKeyboardState*, int midiChannel, int midiNoteNumber, float velocity) override;
    void setMidiInput();
    void oscMessageReceived(const juce::OSCMessage&) override;
    void oscBundleReceived(const juce::OSCBundle&) override;
    int getArgValue(const juce::OSCArgument&);
    float getToneholeValue(int index);
    float getBreathArrayValue(int index);
    //==============================================================================
    void handleAsyncUpdate() override;
    //==============================================================================
    void prepareToPlay (double sampleRate, int samplesPerBlock) override;
    void releaseResources() override;
   #ifndef JucePlugin_PreferredChannelConfigurations
    bool isBusesLayoutSupported (const BusesLayout& layouts) const override;
   #endif
    void processBlock (juce::AudioBuffer<float>&, juce::MidiBuffer&) override;
    //==============================================================================
    juce::AudioProcessorEditor* createEditor() override;
    bool hasEditor() const override;
    //==============================================================================
    const juce::String getName() const override;
    bool acceptsMidi() const override;
    bool producesMidi() const override;
    bool isMidiEffect() const override;
    double getTailLengthSeconds() const override;

    //==============================================================================
    int getNumPrograms() override;
    int getCurrentProgram() override;
    void setCurrentProgram (int index) override;
    const juce::String getProgramName (int index) override;
    void changeProgramName (int index, const juce::String& newName) override;

    //==============================================================================
    void getStateInformation (juce::MemoryBlock& destData) override;
    void setStateInformation (const void* data, int sizeInBytes) override;
    
    int controlNumber;
    int prevControlNumber;
    bool loading;
    float interpVal = 0.0f;

    
    enum Breath {EXHALE, INHALE, NOSIGNAL};
    int breathing = NOSIGNAL;
    
    float currentBreathPos;
    float maxBreathPos = 0.001;
    float currentBreathNeg;

    float maxBreathNeg = 0.001;

    float breath = 0.001;
    
    float audio[2];

    juce::AudioProcessorValueTreeState parameters;

private:
    

    
    /* for osc */
    bool isAddingFromMidiInput = false;
    juce::OSCReceiver rcvr;
    int port;
    juce::TextButton connectButton{"Connect"};
    bool isExhaling = false;
    LEAF leaf;
    //==============================================================================
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (BirlAudioProcessor)
};

#endif
