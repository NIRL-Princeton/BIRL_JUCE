#ifndef EDITOR_H
#define EDITOR_H

#include "PluginProcessor.h"
#include "BirlOSCReceiver.h"

class BirlOneHoleAudioProcessorEditor  : public AudioProcessorEditor

{
public:
    BirlOneHoleAudioProcessorEditor (BirlOneHoleAudioProcessor&);
    ~BirlOneHoleAudioProcessorEditor();
    void paint (Graphics&) override;
    void resized() override;

private:
    BirlOneHoleAudioProcessor& processor;
    BirlOSCReceiver slidersAndKeys;
    ComboBox tuningMenu_;


    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (BirlOneHoleAudioProcessorEditor)
};

#endif


