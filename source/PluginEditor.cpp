#include "PluginEditor.h"
#include "PluginProcessor.h"

BirlOneHoleAudioProcessorEditor::BirlOneHoleAudioProcessorEditor (BirlOneHoleAudioProcessor& p)
    : AudioProcessorEditor (&p), processor(p), slidersAndKeys(p)
{
    setSize (500, 300);

    slidersAndKeys.setBounds(0, 0, 500, 300);
    addAndMakeVisible(slidersAndKeys);

    //    animate.setBounds(180, 180, 80, 80);
    //    addAndMakeVisible(animate);

    tuningMenu_.addItem("Just Intonation", 1);
    tuningMenu_.addItem("Equal Tempered", 2);
    tuningMenu_.addItem("Meantone", 3);
    tuningMenu_.addItem("Highland Bagpipe", 4);
    tuningMenu_.addItem("Add Custom", 5);
    tuningMenu_.setJustificationType(Justification::centred);
    tuningMenu_.setBounds(100, 300, 300, 30);
    //    addAndMakeVisible(tuningMenu_);

    std::unique_ptr<AudioProcessorValueTreeState::ComboBoxAttachment> tuningSelection
        (new AudioProcessorValueTreeState::ComboBoxAttachment(processor.tree, "tuning",
            tuningMenu_));
}

BirlOneHoleAudioProcessorEditor::~BirlOneHoleAudioProcessorEditor()
{
}

void BirlOneHoleAudioProcessorEditor::paint (Graphics& g)
{
    g.fillAll(Colours::palegreen);
    g.setColour(Colours::white);
    g.setFont(18.0f);

}

void BirlOneHoleAudioProcessorEditor::resized()
{
}
