#include "PluginEditor.h"

#include "ui.h"


namespace birl
{
    uint16_t XY_values[NUM_XY_CONTROLS];

    tExpSmooth xy[NUM_XY_CONTROLS];
    float smoothedXY[NUM_XY_CONTROLS];
    
    uint8_t knobPage = 0;
    uint8_t numPages[ControlNil];
    
    uint8_t buttonValues[NUM_BUTTONS];
    uint8_t buttonValuesPrev[NUM_BUTTONS];
    uint8_t cleanButtonValues[NUM_BUTTONS];
    uint32_t buttonHysteresis[NUM_BUTTONS];
    uint32_t buttonHysteresisThreshold = 5;
    uint32_t buttonCounters[NUM_BUTTONS]; // How long a button has been in its current state
    uint32_t buttonHoldThreshold = 200;
    uint32_t buttonHoldMax = 200;
    
    int8_t writeKnobFlag = -1;
    int8_t writeButtonFlag = -1;
    int8_t writeActionFlag = -1;
    
#define NUM_CHARACTERS_PER_CONTROL_NAME 10
    const char* controlNames[ControlNil];
    const char* controlNamesDetails[ControlNil];
    const char* shortControlNames[ControlNil];
    
    const char* knobParamNames[ControlNil][NUM_SYNTH_KNOB_VALUES];

    int8_t currentParamIndex = -1;
    uint8_t orderedParams[8];
    
    uint8_t buttonActionsSFX[NUM_BUTTONS+1][ActionNil];
    uint8_t buttonActionsUI[NUM_BUTTONS+1][ActionNil];
    float displayValues[NUM_SYNTH_KNOB_VALUES];
    int8_t cvAddParam[ControlNil];
    char* (*buttonActionFunctions[ControlNil])(BirlButton, ButtonAction);
    
    BirlControlType currentControl = (BirlControlType)0;
    BirlControlType previousControl = ControlNil;
    uint8_t loadingControl = 0;
    
    void initParamNames()
    {
        for (int i = 0; i < ControlNil; i++)
        {
            cvAddParam[i] = -1;
        }
        
        controlNames[PhysicalModelPM] = "PHYSICAL_MODEL";
        shortControlNames[PhysicalModelPM] = "PMPM";
        controlNamesDetails[PhysicalModelPM] = "FULL_BIRL";
        numPages[PhysicalModelPM] = 1;
        knobParamNames[PhysicalModelPM][0] = "GAIN";
        knobParamNames[PhysicalModelPM][1] = "FUNDAMENTAL";
        knobParamNames[PhysicalModelPM][2] = "NUM_TONEHOLES";
        knobParamNames[PhysicalModelPM][3] = "DCBLOCKER1";
        knobParamNames[PhysicalModelPM][4] = "DCBLOCKER2";
        knobParamNames[PhysicalModelPM][5] = "BIQUAD_COEFF1";
        knobParamNames[PhysicalModelPM][6] = "BIQUAD_COEFF2";
        knobParamNames[PhysicalModelPM][7] = "BIQUAD_COEFF3";
        knobParamNames[PhysicalModelPM][8] = "BIQUAD_COEFF4";
        knobParamNames[PhysicalModelPM][9] = "BIQUAD_COEFF5";
        knobParamNames[PhysicalModelPM][10] = "PF1_CUTOFF";
        knobParamNames[PhysicalModelPM][11] = "PF1_Q";
        knobParamNames[PhysicalModelPM][12] = "PF2_CUTOFF";
        knobParamNames[PhysicalModelPM][13] = "PF2_Q";
        knobParamNames[PhysicalModelPM][14] = "LP1_CUTOFF";
        knobParamNames[PhysicalModelPM][15] = "LP1_Q";
        knobParamNames[PhysicalModelPM][16] = "LP2_CUTOFF";
        knobParamNames[PhysicalModelPM][17] = "LP2_Q";
        knobParamNames[PhysicalModelPM][18] = "SHAPER_DRIVE";
        knobParamNames[PhysicalModelPM][19] = "SHAPER_MIX";
        knobParamNames[PhysicalModelPM][20] = "NOISE_GAIN";
        knobParamNames[PhysicalModelPM][21] = "NOISE_BP_CUTOFF";
        knobParamNames[PhysicalModelPM][22] = "NOISE_BP_Q";

        controlNames[RuleBasedSynth] = "RULE_BASED_SYNTH";
        shortControlNames[RuleBasedSynth] = "RBS";
        controlNamesDetails[RuleBasedSynth] = "TOY_BIRL";
        knobParamNames[RuleBasedSynth][0] = "GAIN";
        knobParamNames[RuleBasedSynth][1] = "ENV1_ATTACK";
        knobParamNames[RuleBasedSynth][2] = "ENV1_DECAY";
        knobParamNames[RuleBasedSynth][3] = "ENV1_SUSTAIN";
        knobParamNames[RuleBasedSynth][4] = "ENV1_RELEASE";
        knobParamNames[RuleBasedSynth][5] = "ENV2_ATTACK";
        knobParamNames[RuleBasedSynth][6] = "ENV2_DECAY";
        knobParamNames[RuleBasedSynth][7] = "ENV2_SUSTAIN";
        knobParamNames[RuleBasedSynth][8] = "ENV2_RELEASE";
        knobParamNames[RuleBasedSynth][9] = "SVF1_CUTOFF";
        knobParamNames[RuleBasedSynth][10] = "SVF1_Q";
        knobParamNames[RuleBasedSynth][11] = "SVF2_CUTOFF";
        knobParamNames[RuleBasedSynth][12] = "SVF2_Q";
        knobParamNames[RuleBasedSynth][13] = "NOISE_GAIN";
        knobParamNames[RuleBasedSynth][14] = "NOISE_BP_CUTOFF";
        knobParamNames[RuleBasedSynth][15] = "NOISE_BP_Q";

    }
    
    void buttonCheck(void)
    {
        for (int i = 0; i < NUM_BUTTONS; i++) {
            if (buttonValues[i] != buttonValuesPrev[i])
            {
                buttonHysteresis[i]++;
            }
            if (cleanButtonValues[i] == 1) {
                buttonActionsSFX[i][ActionHoldContinuous] = true;
                buttonActionsUI[i][ActionHoldContinuous] = true;
                writeButtonFlag = i;
                writeActionFlag = ActionHoldContinuous;
            }
            if (buttonHysteresis[i] < buttonHysteresisThreshold) {
                if (buttonCounters[i] < buttonHoldMax) buttonCounters[i]++;
                if ((buttonCounters[i] >= buttonHoldThreshold) && (cleanButtonValues[i] == 1)){
                    buttonActionsSFX[i][ActionHoldInstant] = true;
                    buttonActionsUI[i][ActionHoldInstant] = true;
                    writeButtonFlag = i;
                    writeActionFlag = ActionHoldInstant;
                }
            }
            else
            {
                cleanButtonValues[i] = buttonValues[i];
                buttonHysteresis[i] = 0;
                buttonCounters[i] = 0;
                
                if (cleanButtonValues[i] == 1)
                {
                    buttonActionsSFX[i][ActionPress] = true;
                    buttonActionsUI[i][ActionPress] = true;
                    writeButtonFlag = i;
                    writeActionFlag = ActionPress;
                }
                else if (cleanButtonValues[i] == 0)
                {
                    buttonActionsSFX[i][ActionRelease] = true;
                    buttonActionsUI[i][ActionRelease] = true;
                    writeButtonFlag = i;
                    writeActionFlag = ActionRelease;
                }
                buttonValuesPrev[i] = buttonValues[i];
            }
        }
        
        if (buttonActionsUI[ButtonPrevControl][ActionPress] == 1)
        {
            previousControl = currentControl;

            if (currentControl <= 0) currentControl = (BirlControlType)((int)ControlNil - 1);
            else currentControl = (BirlControlType)((int)currentControl - 1);

            loadingControl = 1;
            clearButtonActions();
        }
        
        if (buttonActionsUI[ButtonNextControl][ActionPress] == 1)
        {
            previousControl = currentControl;
            
            if (currentControl >= ControlNil - 1) currentControl = (BirlControlType)0;
            else currentControl = (BirlControlType)((int)currentControl + 1);
            
            loadingControl = 1;
            clearButtonActions();
        }
        
        if (buttonActionsUI[ButtonA][ActionPress] == 1)
        {
            if (currentTuning == 0)
            {
               currentTuning = NUM_TUNINGS - 1;
            }
            else
            {
               currentTuning = (currentTuning - 1);
            }
            changeTuning();
            buttonActionsUI[ButtonA][ActionPress] = 0;
        }
        if (buttonActionsUI[ButtonB][ActionPress] == 1)
        {
            currentTuning = (currentTuning + 1) % NUM_TUNINGS;
            changeTuning();
            buttonActionsUI[ButtonB][ActionPress] = 0;
        }
    }
    
    void clearButtonActions()
    {
        for (int b = 0; b < ButtonNil; b++)
        {
            for (int a = 0; a < ActionNil; a++)
            {
                buttonActionsUI[b][a] = 0;
                buttonActionsSFX[b][a] = 0;
                writeButtonFlag = -1;
                writeActionFlag = -1;
            }
        }
    }
    void changeTuning() {
        for (int i = 0; i < 12; i++)
        {
            centsDeviation[i] = tuningPresets[currentTuning][i];
        }
        if (currentTuning == 0)
        {
            // stuff
        }
        else
        {
            // stuff
        }
        if (currentControl == PhysicalModelPM)
        {
            // tune accordingly?
        }
    }
    
//    void incrementPage()
//    {
//        knobPage = (knobPage + 1) % numPages[currentControl];
//        setKnobValues(controlKnobValues[currentControl] + (knobPage * KNOB_PAGE_SIZE));
//    }
//    void decrementPage()
//    {
//        if (knobPage == 0) knobPage = numPages[currentControl] - 1;
//        else knobPage--;
//        setKnobValues(controlKnobValues[currentControl] + (knobPage * KNOB_PAGE_SIZE));
//    }
    void resetKnobValues()
    {
        // idk what to do here
    }
    void setKnobValues(float* values)
    {
        for (int i = 0; i < KNOB_PAGE_SIZE; i++) {
            int knob = i;
            if (knob + (knobPage * KNOB_PAGE_SIZE) == cvAddParam[currentControl]) {
                knob = KNOB_PAGE_SIZE;
            }
//            knobActive[knob] = 0;
        }
    }
    void setKnobValue(int knob, float value)
    {
        if (knob + (knobPage * KNOB_PAGE_SIZE) == cvAddParam[currentControl]) {
            knob = KNOB_PAGE_SIZE;
        }
//            knobActive[knob] = 0;

    }
    void deactivateKnob(int knob)
    {
        // knobActive[knob] = 0;
    }
    void deactivateAllKnobs()
    {
        /* for (int i = 0; i < NUM_ADC_CHANNELS; i++)
            {
                knobActive[i] = 0;
                }
         */
    }
    
    }


