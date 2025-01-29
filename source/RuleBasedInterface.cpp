/*
  ==============================================================================

    RuleBasedInterface.cpp
    Created: 21 Jun 2020 8:57:08pm
    Author:  Emily Liushen

  ==============================================================================
*/

#include "RuleBasedInterface.h"

RuleBasedInterface::RuleBasedInterface(BirlOscillatorAudioProcessor& p) : processor (p)
{
    populateKeySliders();
    populateButtonButtons();
    
    breathpos.setSliderStyle(Slider::LinearVertical);
    breathpos.setBounds(400, 30, 20, 200);
    breathpos.setRange(0.0, 1.0);
    breathpos.setValue(0.0);
    addAndMakeVisible(breathpos);
    
    maxbreathposlabel.setBounds(400, 10, 100, 20);
    maxbreathposlabel.setText("MAX: ", dontSendNotification);
    addAndMakeVisible(maxbreathposlabel);
    
    breathneg.setSliderStyle(Slider::LinearVertical);
    breathneg.setBounds(450, 30, 20, 200);
    breathneg.setRange(0.0, 1.0);
    breathneg.setValue(0.0);
    addAndMakeVisible(breathneg);

    maxbreathneglabel.setBounds(450, 10, 100, 20);
    maxbreathneglabel.setText("MAX: ", dontSendNotification);
    addAndMakeVisible(maxbreathneglabel);
    
    oscReceiver.addListener(this);
    connectButton.setBounds(400, 300, 60, 20);
    addAndMakeVisible(connectButton);
    connectButton.onClick = [this] { connectButtonClicked(); };
    
    monday.setBounds(200, 100, 200, 200);
    monday.setRange(20.0, 3000.0);
    monday.setSliderStyle(Slider::LinearVertical);
    monday.addListener(this);
    addAndMakeVisible(monday);
}
RuleBasedInterface::~RuleBasedInterface()
{
}
//==============================================================================
/**
 */
    void RuleBasedInterface::paint (Graphics& g) {
        
    }
    void RuleBasedInterface::resized() {
    
    }
    void RuleBasedInterface::populateKeySliders()
    {
        for(int i = 0; i < 9; ++i){
            keySliders[i].setSliderStyle(Slider::LinearVertical);
            keySliders[i].setRange(0.0, 1.0);
            keySliders[i].setBounds(100+30*i, 30, 20, 100);
            keySliders[i].setNumDecimalPlacesToDisplay(0);
            addAndMakeVisible(keySliders[i]);
        }
    }
           
    void RuleBasedInterface::populateButtonButtons()
    {
        for(int i = 0; i < 9; ++i) {
            buttonButtons[i].setBounds(30, 5+30*i, 30, 20);
            buttonButtons[i].setButtonText("OFF");
            buttonButtons[i].changeWidthToFitText();
            addAndMakeVisible(buttonButtons[i]);
        }
    }
    
//==============================================================================
/** OSC LISTENER FUNCTIONS
        Receive OSC arguments from the Birl in the categories of keys, buttons, breathpos, breathneg, abs.
        Calls functions corresponding to the categories.
 */
    void RuleBasedInterface::oscMessageReceived (const OSCMessage& message)
    {
        // received a complete list of key input
        if (message.getAddressPattern().toString()=="/birl/keys")
        {
            // every time zero-init the key index to begin new grouping of input
            int keyIndex{0};
            // parse through all the arguments
            for (auto* arg = message.begin(); arg != message.end(); ++arg)
            {
                // what to do for each argument while parsing through all 9 keys
                if (keyIndex < 9)
                { // which it will be the first 9 times this is run
                    // adjust the maximum key argument
                    if (getArgValue(*arg) > maxKeyArg[keyIndex])
                    {
                        maxKeyArg[keyIndex] = getArgValue(*arg);
                    }
                    // clip
                    keys[keyIndex] = (float) getArgValue(*arg)
                                        / (maxKeyArg[keyIndex] * CLIP_RATIO);
                    if (keys[keyIndex] > 1.0){
                        keys[keyIndex] = 1.0;
                    }
                    // set the slider corresponding to the key index
                    // to the arg value.
                    keySliders[keyIndex].setValue(keys[keyIndex],
                                                  dontSendNotification);
                    adjustMidi();
                    ++keyIndex;
                }
            }
        }
           
        // received a complete list of button input
        if (message.getAddressPattern().toString()=="/birl/buttons")
        {
            int buttonIndex{0};
            for (auto* arg = message.begin(); arg != message.end(); ++arg)
            {
                // what to do for each arg while parsing through all 9 buttons
                if (buttonIndex < 9) {
                    // which it will be the first 9 times this is run
                    // set the slider corresponding to the key index = arg value
                    octaveButtonPressed();
                    transpose();
                    if (getArgValue(*arg) == 1) {
                        buttonButtons[buttonIndex].setButtonText("ON");
                        buttons[buttonIndex] = 1;
                    } else {
                        buttonButtons[buttonIndex].setButtonText("OFF");
                        buttons[buttonIndex] = 0;
                    }
                    ++buttonIndex;
                }
            }
        }
           
        //received breathpos value
        if(message.getAddressPattern().toString()=="/birl/breathpos")
        {
            for (auto* arg = message.begin(); arg != message.end(); ++arg)
            {
                // set currentbreathpos to the /birl/breathpos input
                currentbreathpos = getArgValue(*arg);
                // determine whether the user is inhaling or exhaling
                if (currentbreathpos > 0)
                {
                    isExhaling = true;
                }

                // adjust the maxbreathpos if necessary
                if(currentbreathpos > maxbreathpos)
                {
                    setMaxBreathPosValue(currentbreathpos);
                }
                // set the value of the slider to the fraction equal to
                // (currentbreathvalue - absbreath)/(maxbreathvalue-absbreath)
                breathpos.setValue(currentbreathpos/maxbreathpos,
                                   dontSendNotification);
                if(isExhaling){
                    processor.gain = (currentbreathpos/maxbreathpos);
//                    processor.filterfrequency = LEAF_clip(30.0, processor.gain
//                                    * 15 * processor.frequency, 17000.0);
                }
            }
        }
        //received breathneg value
        if(message.getAddressPattern().toString()=="/birl/breathneg")
        {
            for (auto* arg = message.begin(); arg != message.end(); ++arg)
            {
                currentbreathneg = getArgValue(*arg);
                if (currentbreathneg > 0)
                {
                    isExhaling = false;
                }
                if(currentbreathneg > maxbreathneg)
                {
                    setMaxBreathNegValue(currentbreathneg);
                }
                breathneg.setValue(currentbreathneg/maxbreathneg,
                                   dontSendNotification);
                if (!isExhaling)
                {
                    processor.gain = (currentbreathneg/maxbreathneg);
//                    processor.filterfrequency = LEAF_clip(30.0, processor.gain
//                                     * 15 * processor.frequency, 17000.0);
                }
            }
        }
        if(message.getAddressPattern().toString()=="/birl/abs")
        {
            for(auto* arg = message.begin(); arg != message.end(); ++arg)
            {
                absbreath = getArgValue(*arg);
            }
        }
        triggerAsyncUpdate();
    }
           
    void RuleBasedInterface::oscBundleReceived (const OSCBundle& bundle)
    {
        for (auto* element = bundle.begin(); element != bundle.end(); ++element)
          {
            if (element->isMessage())
                oscMessageReceived (element->getMessage());
            else if (element->isBundle())
                oscBundleReceived(element->getBundle());
        }
        triggerAsyncUpdate();
    }

    int RuleBasedInterface::getArgValue(const OSCArgument& arg)
    {
        if (arg.isInt32()) {
            return arg.getInt32();
        }
        else
        {
            handleArgTypeError();
            return 0;
        }
    }
    void RuleBasedInterface::setMaxBreathPosValue(double newMaxPos)
    {
        maxbreathpos = newMaxPos;
        updateMaxBreathLabel(&maxbreathposlabel, maxbreathpos);
    }
    void RuleBasedInterface::setMaxBreathNegValue(double newMaxNeg)
    {
        maxbreathneg = newMaxNeg;
        updateMaxBreathLabel(&maxbreathneglabel, maxbreathneg);
    }
    void RuleBasedInterface::updateMaxBreathLabel(Label* label, double newMax)
    {
        label->setText("MAX: " + (String)(newMax), dontSendNotification);
        addAndMakeVisible(label);
    }
    void RuleBasedInterface::handleAsyncUpdate() {}
       
//==============================================================================
/** Functions for the Connect/Disconnect button.
    Checks whether Birl is properly connected and updates button text.
 */
    void RuleBasedInterface::connect() {
        int portToConnect = BIRL_PORT;
        if (oscReceiver.connect(portToConnect)) {
            currentPortNumber = portToConnect;
            connectButton.setButtonText("Disconnect");
        }
    }
    void RuleBasedInterface::disconnect()
    {
        if (oscReceiver.disconnect()) {
            currentPortNumber = -1;
            connectButton.setButtonText("Connect");
        }
    }
    bool RuleBasedInterface::isConnected() const
    {
        return currentPortNumber == BIRL_PORT;
    }
    void RuleBasedInterface::connectButtonClicked()
    {
        if (!isConnected())
            connect();
        else
            disconnect();
    }
    void RuleBasedInterface::handleArgTypeError()
    {
        AlertWindow::showMessageBoxAsync(AlertWindow::WarningIcon,
            "Argument Type Error",
            "Error: Birl input not integer values",
            "OK");
    }
    
//==============================================================================
/** Listens to the fingerings and calculates a double to adjust the current MIDI value.
 */
    void RuleBasedInterface::adjustMidi()
    {
        double midiAdjustment = 12.0;
        if (keys[0] > 0.0)
        {
            midiAdjustment -= keys[0];
            if (keys[0] == 1.0)
            {
                midiAdjustment -= 2 * keys[1];
                if (keys[1] == 1.0)
                {
                    midiAdjustment -= 2 * keys[2];
                }
                
            }
            if (keys[1] == 1.0 && keys[2] == 1.0)
            {
                if (keys[4] > 0.0)
                {
                    midiAdjustment -= 2 * keys[4];
                    if (keys[4] == 1.0)
                    {
                        midiAdjustment -= keys[5];
                        if (keys[5] == 1.0)
                        {
                            midiAdjustment -= 2 * keys[6];
                            if (keys[6] == 1.0)
                            {
                                midiAdjustment -= keys[7];
                                midiAdjustment -= 2 * keys[8];
                            }
                        }
                    }
                }
                if (keys[4] == 0.0 && keys[5] > 0.0)
                {
                    midiAdjustment -= keys[5];
                }

            }
        }
        else if (keys[1] > 0.0)
        {
            midiAdjustment -= keys[1];
        }
        midiAdjustment += keys[3];

        calculateFrequency(midiAdjustment + currentMidi + octaveTransposition);
    }
    
//==============================================================================
/** Takes the result of the calculations in adjustMidi() and sets the processor frequency.
    Uses LEAF's MIDI-to-Frequency function.
 
    The  basic algorithm:
    1. Separate MIDI floating point input into its integer ("floor") and non-integer ("decimals") components.
    2. Calculate the pitch class (0-11) from the integer component.
    3. Multiply floor variable by the proper scale based on the current base pitch and scale choice).
    4. Interpolate the additional offset by multiplying the non-integer component by the difference
        between the MIDI floor note and the note below it that is partially covered.
    5. Convert the MIDI floating point value into a frequency.
    6. Set processor gain to that frequency.
 */

    void RuleBasedInterface::calculateFrequency(double midiAdjustment)
    {
        
        // split midiAdjustment into its floor and decimal part
        double floor;
        double decimals = modf(midiAdjustment, &floor);
//        DBG("floor = " + (String)floor + ", decimals = " + (String)decimals);

        
        // calculate pitch class (0-11) from intpart
        int pitchclass = (int)floor % 12 - 5;
//        DBG("pitchclass = " + (String)pitchclass);

        // declare our notes
        double note{};
        double noteAnotherFinger{};
        
        /* multiple the int part by the proper scale
         * (based on current base pitch and scale choice)
         * if no scale choice, default to equal-tempered
         */
        
        int scaleType = getScaleType();
        if (scaleType == 0)
        {
            note = normalscale[pitchclass];
            noteAnotherFinger = normalscale[pitchclass-1];
        }
        if (scaleType == 1)
        {
            note = otonalscale[pitchclass];
            noteAnotherFinger = otonalscale[pitchclass-1];
        }
        if (scaleType == 2)
        {
            note = utonalscale[pitchclass];
            noteAnotherFinger = utonalscale[pitchclass-1];
                
        }
        if (scaleType == 4)
        {
            note = meantonescale[pitchclass];
            noteAnotherFinger = meantonescale[pitchclass-1];
        }
        else
        {
            note = equaltemperedscale[pitchclass];
            noteAnotherFinger = equaltemperedscale[pitchclass-1];
        }

        // set processor frequency
        processor.frequency = LEAF_midiToFrequency(floor
                                    + decimals * (note - noteAnotherFinger));
        
//        DBG("midiAdjustment = " + (String)midiAdjustment
//            + ", note = " + (String)note
//            + ", noteAnotherFinger = " + (String)noteAnotherFinger);
//        DBG("processor.frequency = " + (String) processor.frequency);
    }
    
//==============================================================================
/** Functions for the Birl buttons.
 */

    void RuleBasedInterface::octaveButtonPressed()
    {
        if (buttons[0] == 1) {
            octaveTransposition = 12.0;
        } else if (buttons[1] == 1){
            octaveTransposition = -12.0;
        } else {
            octaveTransposition = 0.0;
        }
    }
    
    void RuleBasedInterface::transpose()
    {
        if (buttons[2] == 0)
        {
            upIsLocked = false;
        }
        if (buttons[2] == 1 && !upIsLocked){
            ++currentMidi;
            upIsLocked = true;
//            DBG("transposition: " + (String)currentMidi);
        }
        if (buttons[3] == 0)
        {
            downIsLocked = false;
        }
        if (buttons[3] == 1 && !downIsLocked){
            --currentMidi;
            downIsLocked = true;
//            DBG("transposition: " + (String)currentMidi);
        }
    }

    int RuleBasedInterface::getScaleType()
    {
        int tempScaleButtons[5]{0, 0, 0, 0, 0};

        for (int tempScaleButtonsIndex = 0; tempScaleButtonsIndex < 5;
             ++tempScaleButtonsIndex) {
            if (buttons[tempScaleButtonsIndex + 4] == 1)
            {
                tempScaleType = tempScaleButtonsIndex;
                tempScaleButtons[tempScaleButtonsIndex] = 1;
            }
            else {
                tempScaleButtons[tempScaleButtonsIndex] = 0;
            }
        }
//        DBG("Scale Type: " + (String)(tempScaleType));
        return tempScaleType;
    }
void RuleBasedInterface::sliderValueChanged(Slider* slider)
{
    if (slider == &monday)
    {
        processor.modAmount = slider->getValue();

    }
}



