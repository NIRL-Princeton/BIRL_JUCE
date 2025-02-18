#include "PluginEditor.h"
#include "BinaryData.h"
#include <iostream>
//==============================================================================
BirlAudioProcessorEditor::BirlAudioProcessorEditor (BirlAudioProcessor& p)
    : AudioProcessorEditor (&p), processor (p), control_parameter (processor.parameters, "control", control_combo_box), synth_parameter (processor.parameters, "synth", synth_combo_box)
{
    // if the synth is physical model do this
    panel = Drawable::createFromImageData (BinaryData::panelpm, BinaryData::panelpmSize);
    panel->setBounds (Rectangle<int> (0, 0, 700, 700));
    addAndMakeVisible (*panel);

    // birl::initParamNames();

    /* envelope knobs for physical model */

    for (int i = 0; i < 21; i++)
    {
        knobs.add (new Slider());
        knobs[i]->setLookAndFeel (&knob);
        knobs[i]->setSliderStyle (Slider::RotaryVerticalDrag);
        knobs[i]->setTextBoxStyle (Slider::TextBoxBelow, false, 0, 0);
        knobs[i]->setTextBoxIsEditable (true);
        knobs[i]->setRange (0, 1);
        addAndMakeVisible (knobs[i]);
    }

    knobs[0]->setBounds (160, 190, 40, 40); // dc blocker 1
    knobs[1]->setBounds (215, 190, 40, 40); // dc blocker 2
    knobs[2]->setBounds (282, 190, 40, 40); // shaper drive
    knobs[3]->setBounds (337, 190, 40, 40); // shaper mix

    knobs[4]->setBounds (160, 290, 40, 40); // pf cutoff 1
    knobs[5]->setBounds (215, 290, 40, 40); // pf Q 1
    knobs[6]->setBounds (282, 290, 40, 40); // pf cutoff 2
    knobs[7]->setBounds (337, 290, 40, 40); // pf Q 2

    knobs[8]->setBounds (160, 395, 40, 40); // lp cutoff 1
    knobs[9]->setBounds (215, 395, 40, 40); // lp Q 1
    knobs[10]->setBounds (282, 395, 40, 40); // lp cutoff 2
    knobs[11]->setBounds (337, 395, 40, 40); // lp Q 2

    knobs[12]->setBounds (153, 495, 40, 40); // bq
    knobs[13]->setBounds (200, 495, 40, 40); // bq2
    knobs[14]->setBounds (247, 495, 40, 40); // bq3
    knobs[15]->setBounds (294, 495, 40, 40); // bq4
    knobs[16]->setBounds (341, 495, 40, 40); // bq5

    knobs[17]->setBounds (409, 191, 40, 40); // noise cutoff
    knobs[18]->setBounds (409, 280, 40, 40); // noise Q
    knobs[19]->setBounds (409, 351, 40, 40); // noise gain

    knobs[20]->setBounds (419, 470, 60, 60); // master gain

    sliderAttachments.add (new juce::AudioProcessorValueTreeState::SliderAttachment (p.parameters, "dcblocker1", *knobs[0]));
    sliderAttachments.add (new juce::AudioProcessorValueTreeState::SliderAttachment (p.parameters, "dcblocker2", *knobs[1]));
    sliderAttachments.add (new juce::AudioProcessorValueTreeState::SliderAttachment (p.parameters, "shaper_drive", *knobs[2]));
    sliderAttachments.add (new juce::AudioProcessorValueTreeState::SliderAttachment (p.parameters, "shaper_mix", *knobs[3]));

    sliderAttachments.add (new juce::AudioProcessorValueTreeState::SliderAttachment (p.parameters, "pf1_cutoff", *knobs[4]));
    sliderAttachments.add (new juce::AudioProcessorValueTreeState::SliderAttachment (p.parameters, "pf1_q", *knobs[5]));
    sliderAttachments.add (new juce::AudioProcessorValueTreeState::SliderAttachment (p.parameters, "pf2_cutoff", *knobs[6]));
    sliderAttachments.add (new juce::AudioProcessorValueTreeState::SliderAttachment (p.parameters, "pf2_q", *knobs[7]));

    sliderAttachments.add (new juce::AudioProcessorValueTreeState::SliderAttachment (p.parameters, "lp1_cutoff", *knobs[8]));
    sliderAttachments.add (new juce::AudioProcessorValueTreeState::SliderAttachment (p.parameters, "lp2_q", *knobs[9]));
    sliderAttachments.add (new juce::AudioProcessorValueTreeState::SliderAttachment (p.parameters, "lp2_cutoff", *knobs[10]));
    sliderAttachments.add (new juce::AudioProcessorValueTreeState::SliderAttachment (p.parameters, "lp2_q", *knobs[11]));

    sliderAttachments.add (new juce::AudioProcessorValueTreeState::SliderAttachment (p.parameters, "biquad_coeff1", *knobs[12]));
    sliderAttachments.add (new juce::AudioProcessorValueTreeState::SliderAttachment (p.parameters, "biquad_coeff2", *knobs[13]));
    sliderAttachments.add (new juce::AudioProcessorValueTreeState::SliderAttachment (p.parameters, "biquad_coeff3", *knobs[14]));
    sliderAttachments.add (new juce::AudioProcessorValueTreeState::SliderAttachment (p.parameters, "biquad_coeff4", *knobs[15]));
    sliderAttachments.add (new juce::AudioProcessorValueTreeState::SliderAttachment (p.parameters, "biquad_coeff5", *knobs[16]));

    sliderAttachments.add (new juce::AudioProcessorValueTreeState::SliderAttachment (p.parameters, "noise_bp_cutoff", *knobs[17]));
    sliderAttachments.add (new juce::AudioProcessorValueTreeState::SliderAttachment (p.parameters, "noise_bp_q", *knobs[18]));
    sliderAttachments.add (new juce::AudioProcessorValueTreeState::SliderAttachment (p.parameters, "noise_gain", *knobs[19]));

    sliderAttachments.add (new juce::AudioProcessorValueTreeState::SliderAttachment (p.parameters, "gain", *knobs[20]));

    //    gainAttachment = std::make_unique<AudioProcessorValueTreeState::SliderAttachment>(processor.parameters, "gain", knobs[20]);

    diagram.setBounds (0, 184, 140, 385);
    addAndMakeVisible (diagram);

    breathSensors.setBounds (0, 0, 300, 300);
    addAndMakeVisible (breathSensors);

    /* if the synth is ss do this
      * panel = Drawable::createFromImageData(BinaryData::panelsynth_svg, BinaryData::panelsynth_svgSize);
      * panel->setBounds(Rectangle<int>(0, 0, 700, 700));
      * addAndMakeVisible(*panel);
      */
    //
    //
    //    control_parameter = std::make_unique<AudioProcessorValueTreeState::ComboBoxAttachment>(processor.parameters, "control", control_combo_box);
    //    synth_parameter = std::make_unique<AudioProcessorValueTreeState::ComboBoxAttachment>(processor.parameters, "synth", synth_combo_box);
    //
    //
    //    else if (i < 8) {
    //        knobs[i]->setBounds(Rectangle<int>(250+50*i, 255, 40, 40));
    //    }
    //    else if (i < 10) {
    //        knobs[i]->setBounds(450+50*(i-8), 400, 40, 40);
    //    }
    //    else if (i < 12) {
    //        knobs[i]->setBounds(450+50*(i-10), 450, 40, 40);
    //    }
    //    else if (i < 15) {
    //        knobs[i]->setBounds(500+50*(i-12), 500, 40, 40);
    //    }
    //    else {
    //        knobs[i]->setBounds(500, 500, 60, 60);
    //    }

    /* oscillator buttons */

    //    Path path;
    //    path.addEllipse(0, 0, 30, 30);
    //
    //    for (int i = 0; i < 8; i++) {
    //        buttons.add(new ShapeButton("", Colours::darkgreen, Colours::lightgreen, Colours::lightgreen));
    //        buttons[i]->setShape(path, true, true, true);
    //        addAndMakeVisible(buttons[i]);
    //    }
    //
    //    buttons[0]->setBounds(280+20, 200, 22, 22);
    //    buttons[1]->setBounds(310+20, 200, 22, 22);
    //    buttons[2]->setBounds(340+20, 200, 22, 22);
    //    buttons[3]->setBounds(370+20, 200, 22, 22);
    //    buttons[4]->setBounds(280+20, 236, 22, 22);
    //    buttons[5]->setBounds(310+20, 236, 22, 22);
    //    buttons[6]->setBounds(340+20, 236, 22, 22);
    //    buttons[7]->setBounds(370+20, 236, 22, 22);

    setSize (700, 700);

    //    addAndMakeVisible(control_menu);
    //    control_menu.setBounds(Rectangle<int>(155, 130, 130, 20));
    //    control_menu.addItem("Physical Model", 1);
    //    control_menu.addItem("Rule-Based", 2);
    //    control_menu.addItem("Neural Net", 3);
    //    control_menu.setColour(ComboBox::ColourIds::backgroundColourId, Colours::transparentWhite);
    //    control_menu.setColour(ComboBox::ColourIds::textColourId, Colours::white);
    //    control_menu.setColour(ComboBox::ColourIds::buttonColourId, Colours::black);
    //    control_menu.setColour(ComboBox::ColourIds::arrowColourId, Colours::white);
    //    control_menu.setColour(ComboBox::ColourIds::outlineColourId, Colours::transparentWhite);
    //    control_menu.setJustificationType(Justification::centred);
    //    control_menu.setSelectedId(1, dontSendNotification);
    //    control_menu.onChange = [this] { controlChanged(); };

    //    addAndMakeVisible(*control_combo_box);
    //    control_combo_box->setBounds(Rectangle<int>(155, 130, 130, 20));
    //    control_combo_box->addItem("Physical Model", 1);
    //    control_combo_box->addItem("Rule-Based", 2);
    //    control_combo_box->addItem("Neural Net", 3);
    //    control_combo_box->setColour(ComboBox::ColourIds::backgroundColourId, Colours::transparentWhite);
    //    control_combo_box->setColour(ComboBox::ColourIds::textColourId, Colours::white);
    //    control_combo_box->setColour(ComboBox::ColourIds::buttonColourId, Colours::black);
    //    control_combo_box->setColour(ComboBox::ColourIds::arrowColourId, Colours::white);
    //    control_combo_box->setColour(ComboBox::ColourIds::outlineColourId, Colours::transparentWhite);
    //    control_combo_box->setJustificationType(Justification::centred);
    //    control_combo_box->setSelectedId(1, dontSendNotification);
    //    control_combo_box->onChange = [this] { controlChanged(); };

    //    addAndMakeVisible(synth_menu);
    //    synth_menu.setBounds(Rectangle<int>(310, 130, 150, 20));
    //    synth_menu.addItem("Physical Model", 1);
    //    synth_menu.addItem("Subtractive Synthesis", 2);
    //    synth_menu.setColour(ComboBox::ColourIds::backgroundColourId, Colours::transparentWhite);
    //    synth_menu.setColour(ComboBox::ColourIds::textColourId, Colours::white);
    //    synth_menu.setColour(ComboBox::ColourIds::buttonColourId, Colours::black);
    //    synth_menu.setColour(ComboBox::ColourIds::arrowColourId, Colours::white);
    //    synth_menu.setColour(ComboBox::ColourIds::outlineColourId, Colours::transparentWhite);
    //    synth_menu.setJustificationType(Justification::centred);
    //    synth_menu.setSelectedId(1, dontSendNotification);
    //    synth_menu.onChange = [this] { synthChanged(); };

    //    addAndMakeVisible(*synth_combo_box);
    //    synth_combo_box->setBounds(Rectangle<int>(310, 130, 150, 20));
    //    synth_combo_box->addItem("Physical Model", 1);
    //    synth_combo_box->addItem("Subtractive Synthesis", 2);
    //    synth_combo_box->setColour(ComboBox::ColourIds::backgroundColourId, Colours::transparentWhite);
    //    synth_combo_box->setColour(ComboBox::ColourIds::textColourId, Colours::white);
    //    synth_combo_box->setColour(ComboBox::ColourIds::buttonColourId, Colours::black);
    //    synth_combo_box->setColour(ComboBox::ColourIds::arrowColourId, Colours::white);
    //    synth_combo_box->setColour(ComboBox::ColourIds::outlineColourId, Colours::transparentWhite);
    //    synth_combo_box->setJustificationType(Justification::centred);
    //    synth_combo_box->setSelectedId(1, dontSendNotification);
    //    synth_combo_box->onChange = [this] { controlChanged(); };

    startTimer (100);
}

BirlAudioProcessorEditor::~BirlAudioProcessorEditor()
{
    for (int i = 0; i < NUM_KNOBS; i++)
    {
        knobs[i]->setLookAndFeel (nullptr);
    }
}

//==============================================================================
void BirlAudioProcessorEditor::paint (juce::Graphics& g)
{
    g.setGradientFill (ColourGradient (Colour (30, 30, 30), juce::Point<float> (0, 0), Colour (10, 10, 10), juce::Point<float> (0, 700), false));
    g.fillRect (0, 0, 700, 700);

    // posbreath sensor
    g.setColour (Colours::darkgreen);
    g.fillEllipse (87, 40, 22, 22);
    // if you're inhaling!
    //        g.setGradientFill(ColourGradient(Colours::lightgreen, juce::Point<float>(98, 51), Colours::transparentWhite, juce::Point<float>(87, 40), true));
    //        g.fillEllipse(83, 36, 30, 30);

    //negbreath sensor
    g.setColour (Colours::darkred);
    g.fillEllipse (87, 76, 22, 22);
    //        g.setGradientFill(ColourGradient(Colours::pink, juce::Point<float>(98, 87), Colours::transparentWhite, juce::Point<float>(87, 76), true));
    //        g.fillEllipse(83, 72, 30, 30);
}

void BirlAudioProcessorEditor::resized()
{
    // This is generally where you'll want to lay out the positions of any
    // subcomponents in your editor..
}
void BirlAudioProcessorEditor::controlChanged()
{
    if (control_combo_box.getSelectedIdAsValue() == 4 || control_combo_box.getSelectedIdAsValue() == 5)
    {
        // do stuff?
    }
}

void BirlAudioProcessorEditor::synthChanged() {}

void BirlAudioProcessorEditor::changeToneholeGraphicSize (Graphics& g, int toneholeIndex)
{
    g.setColour (Colours::white);
    //    float x =
    g.fillEllipse (0, 0, 30, 30);
}

void BirlAudioProcessorEditor::timerCallback()
{
    for (int i = 0; i < 9; i++)
    {
        newDiameter[i] = 2 * (int) sqrt (16 * 16 * processor.getToneholeValue (i));
        diagram.setDiameter (i, newDiameter[i]);
        diagram.setX (i, diagram.getX_ (i) + 0.5 * (diagram.getD_ (i) - newDiameter[i]));
        diagram.setY (i, diagram.getY_ (i) + 0.5 * (diagram.getD_ (i) - newDiameter[i]));
    }
    for (int i = 0; i < 2; i++)
    {
        newBreathDiameter[i] = 2 * (int) sqrt (15 * 15 * processor.getBreathArrayValue (i));
        breathSensors.setDiameter (i, newBreathDiameter[i]);
        breathSensors.setX (i, breathSensors.getX_ (i) + 0.5 * (breathSensors.getD_ (i) - newBreathDiameter[i]));
        breathSensors.setY (i, breathSensors.getY_ (i) + 0.5 * (breathSensors.getD_ (i) - newBreathDiameter[i]));
    }
    repaint (diagram.getBounds());
    repaint (breathSensors.getBounds());
}
