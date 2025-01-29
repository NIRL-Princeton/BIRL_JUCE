#ifndef EDITOR_H
#define EDITOR_H

#include "PluginProcessor.h"
#include <juce_gui_basics/juce_gui_basics.h>
#include <juce_audio_processors/juce_audio_processors.h>
using namespace juce;

class OtherLookAndFeel : public LookAndFeel_V4
{
    void drawRotarySlider(Graphics& g, int x, int y, int width, int height, float sliderPos, const float rotaryStartAngle, const float rotaryEndAngle, Slider& slider) override
    {
        auto radius = jmin(width / 2, height / 2) - 4.0f;
        auto centreX = x + width * 0.5f;
        auto centreY = y + height * 0.5f;
        auto rx = centreX - radius;
        auto ry = centreY - radius;
        auto rw = radius * 2.0f;
        auto angle = rotaryStartAngle + sliderPos * (rotaryEndAngle - rotaryStartAngle);
        Array<Colour> colours{ Colour(117, 117, 117), Colour(136, 136, 137), Colour(140, 139, 139), Colour(168, 167, 167), Colour(166, 163, 163), Colour(208, 208, 208), Colour(207, 207, 207) };
        //auto colorOne = Colour(117, 117, 117);

        Path p;
        auto pointerLength = radius * 0.33f;
        auto pointerThickness = 2.0f;
        p.addRectangle(-pointerThickness * 0.5f, -radius, pointerThickness, pointerLength);
        p.applyTransform(AffineTransform::rotation(angle).translated(centreX, centreY));
        g.setColour(Colours::dimgrey);
        g.fillEllipse(rx - 2, ry - 2, rw + 4, rw + 6);
        g.setColour(Colours::white);
        g.fillEllipse(rx, ry, rw, rw);
        g.setColour(Colours::black);
        g.fillPath(p);

    }
    
    Label* createSliderTextBox (Slider& slider) override
    {
        auto* l = LookAndFeel_V2::createSliderTextBox (slider);
        l->setColour (Label::textColourId, Colours::white);
        return l;
    }

};

class BirlInteractiveDiagram : public Component {
public:
    BirlInteractiveDiagram() {}
    void paint(Graphics& g) override {
        g.setColour(Colours::white);
        for (int i = 0; i < 9; i++) {
            g.fillEllipse(x[i], y[i], d[i], d[i]);
        }
    }
    void setDiameter (int index, int diameter) {
        d[index] = diameter;
    }
    void setX (int index, int xcoor) {
        x[index] = xcoor;
    }
    void setY (int index, int ycoor) {
        y[index] = ycoor;
    }
    int getD_(int index) { return d_[index]; }
    int getX_(int index) { return x_[index]; }
    int getY_(int index) { return y_[index]; }
private:
    const int x_[9] = {82, 82, 82, 94, 84, 84, 84, 66, 66};
    const int y_[9] = {0, 40, 79, 116, 174, 215, 255, 289, 329};
    const int d_[9] = {32, 32, 32, 32, 32, 32, 32, 32, 32};
    
    int x[9] = {82, 82, 82, 94, 84, 84, 84, 66, 66};
    int y[9] = {0, 40, 79, 116, 174, 215, 255, 289, 329};
    int d[9] = {32, 32, 32, 32, 32, 32, 32, 32, 32};
    
};

class BirlInteractiveBreathSensor : public Component {
public:
    BirlInteractiveBreathSensor() {}
    void paint(Graphics& g) override {
        g.setGradientFill(ColourGradient(Colours::lightgreen, juce::Point<float>(x[0] + d[0] * 0.5, y[0] + d[0] * 0.5), Colours::transparentWhite, juce::Point<float>(x[0] + lightRadius, y[0] + lightRadius), true));
        g.fillEllipse(x[0], y[0], d[0], d[0]);
        g.setGradientFill(ColourGradient(Colours::pink, juce::Point<float>(x[1] + d[1] * 0.5, y[1] + d[1] * 0.5), Colours::transparentWhite, juce::Point<float>(x[1] + lightRadius, y[1] + lightRadius), true));
        g.fillEllipse(x[1], y[1], d[1], d[1]);
    }
    void setDiameter (int index, int diameter) {
        d[index] = diameter;
    }
    void setX (int index, int xcoor) {
        x[index] = xcoor;
    }
    void setY (int index, int ycoor) {
        y[index] = ycoor;
    }
    int getD_(int index) { return d_[index]; }
    int getX_(int index) { return x_[index]; }
    int getY_(int index) { return y_[index]; }

private:
    const int x_[2] = {83, 83};
    const int y_[2] = {36, 72};
    const int d_[2] = {30, 30};
    const int lightRadius = 2;
    
    int x[2] = {83, 83};
    int y[2] = {36, 72};
    int d[2] = {30, 30};
};


class BirlAudioProcessorEditor  : public AudioProcessorEditor, public Timer
{
public:
    BirlAudioProcessorEditor (BirlAudioProcessor&);
    ~BirlAudioProcessorEditor();

    //==============================================================================
    void paint (Graphics&) override;

    void resized() override;
//    void sliderValueChanged(Slider* slider) override;
    void controlChanged();
//    void buttonClicked(Button* button) override;
//    void buttonStateChanged(Button* button) override;
    void timerCallback() override;
    
    void synthChanged();
    void changeToneholeGraphicSize(Graphics&, int);

private:
    #define NUM_KNOBS 4 // test

    BirlAudioProcessor& processor;
    BirlInteractiveDiagram diagram;
    BirlInteractiveBreathSensor breathSensors;
    
    // ui stuff
    std::unique_ptr<Drawable> logo;
    std::unique_ptr<Drawable> panel;
    Image baseline;
    
    int newDiameter[9];
    int newBreathDiameter[2];
    
    Array<Slider*> knobs;
    Array<ShapeButton*> buttons;
//    ComboBox control_menu;
//    ComboBox synth_menu;
    ComboBox control_combo_box;
    ComboBox synth_combo_box;
    ShapeButton* train_button;
    
    OtherLookAndFeel knob;
    
//    std::unique_ptr<ComboBoxParameterAttachment> control_parameter;
//    std::unique_ptr<ComboBoxParameterAttachment> synth_parameter;
    juce::AudioProcessorValueTreeState::ComboBoxAttachment control_parameter;
    juce::AudioProcessorValueTreeState::ComboBoxAttachment synth_parameter;
    
    
    
    
public:
    // attachments
    
//    std::unique_ptr<juce::AudioProcessorValueTreeState::SliderAttachment> gainAttachment;


    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (BirlAudioProcessorEditor)
};


#endif
