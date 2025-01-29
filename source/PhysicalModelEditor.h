#ifndef PHYSICAL_MODEL_EDITOR_H
#define PHYSICAL_MODEL_EDITOR_H

#include <juce_gui_basics/juce_gui_basics.h>
#include <math.h>
#include "PhysicalModel.h"


class PhysicalModelEditor : public Component
{
public:
    PhysicalModelEditor()
    {
        addAndMakeVisible(connectButton);
        connectButton.setBounds(400, 20, 80, 30);

        setSize (320, 200);
    }
    ~PhysicalModelEditor(){}

    void paint (Graphics& g) override
    {
        g.fillAll (Colours::black);
    }

    void resized() override
    {
    }

protected:
    Slider keySliders[NUM_OF_KEYS];
    TextButton buttonButtons[NUM_OF_KEYS];
    TextButton connectButton {"Connect"};

    
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (PhysicalModelEditor)
};
#endif
