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

    uint8_t buttonActionsSFX[NUM_BUTTONS + 1][ActionNil];
    uint8_t buttonActionsUI[NUM_BUTTONS + 1][ActionNil];
    float displayValues[NUM_SYNTH_KNOB_VALUES];
    int8_t cvAddParam[ControlNil];
    char* (*buttonActionFunctions[ControlNil]) (BirlButton, ButtonAction);

    BirlControlType currentControl = (BirlControlType) 0;
    BirlControlType previousControl = ControlNil;
    uint8_t loadingControl = 0;

    void initParamNames()
    {
        for (int i = 0; i < ControlNil; i++)
        {
            cvAddParam[i] = -1;
        }

        controlNames[PhysicalModelPM] = "physical_model";
        shortControlNames[PhysicalModelPM] = "pmpm";
        controlNamesDetails[PhysicalModelPM] = "full_birl";
        numPages[PhysicalModelPM] = 1;
        knobParamNames[PhysicalModelPM][0] = "gain";
        knobParamNames[PhysicalModelPM][1] = "fundamental";
        knobParamNames[PhysicalModelPM][2] = "num_toneholes";
        knobParamNames[PhysicalModelPM][3] = "dcblocker1";
        knobParamNames[PhysicalModelPM][4] = "dcblocker2";
        knobParamNames[PhysicalModelPM][5] = "biquad_coeff1";
        knobParamNames[PhysicalModelPM][6] = "biquad_coeff2";
        knobParamNames[PhysicalModelPM][7] = "biquad_coeff3";
        knobParamNames[PhysicalModelPM][8] = "biquad_coeff4";
        knobParamNames[PhysicalModelPM][9] = "biquad_coeff5";
        knobParamNames[PhysicalModelPM][10] = "pf1_cutoff";
        knobParamNames[PhysicalModelPM][11] = "pf1_q";
        knobParamNames[PhysicalModelPM][12] = "pf2_cutoff";
        knobParamNames[PhysicalModelPM][13] = "pf2_q";
        knobParamNames[PhysicalModelPM][14] = "lp1_cutoff";
        knobParamNames[PhysicalModelPM][15] = "lp1_q";
        knobParamNames[PhysicalModelPM][16] = "lp2_cutoff";
        knobParamNames[PhysicalModelPM][17] = "lp2_q";
        knobParamNames[PhysicalModelPM][18] = "shaper_drive";
        knobParamNames[PhysicalModelPM][19] = "shaper_mix";
        knobParamNames[PhysicalModelPM][20] = "noise_gain";
        knobParamNames[PhysicalModelPM][21] = "noise_bp_cutoff";
        knobParamNames[PhysicalModelPM][22] = "noise_bp_q";

        controlNames[RuleBasedSynth] = "rule_based_synth";
        shortControlNames[RuleBasedSynth] = "rbs";
        controlNamesDetails[RuleBasedSynth] = "toy_birl";
        knobParamNames[RuleBasedSynth][0] = "gain";
        knobParamNames[RuleBasedSynth][1] = "env1_attack";
        knobParamNames[RuleBasedSynth][2] = "env1_decay";
        knobParamNames[RuleBasedSynth][3] = "env1_sustain";
        knobParamNames[RuleBasedSynth][4] = "env1_release";
        knobParamNames[RuleBasedSynth][5] = "env2_attack";
        knobParamNames[RuleBasedSynth][6] = "env2_decay";
        knobParamNames[RuleBasedSynth][7] = "env2_sustain";
        knobParamNames[RuleBasedSynth][8] = "env2_release";
        knobParamNames[RuleBasedSynth][9] = "svf1_cutoff";
        knobParamNames[RuleBasedSynth][10] = "svf1_q";
        knobParamNames[RuleBasedSynth][11] = "svf2_cutoff";
        knobParamNames[RuleBasedSynth][12] = "svf2_q";
        knobParamNames[RuleBasedSynth][13] = "noise_gain";
        knobParamNames[RuleBasedSynth][14] = "noise_bp_cutoff";
        knobParamNames[RuleBasedSynth][15] = "noise_bp_q";
    }

    void buttonCheck (void)
    {
        for (int i = 0; i < NUM_BUTTONS; i++)
        {
            if (buttonValues[i] != buttonValuesPrev[i])
            {
                buttonHysteresis[i]++;
            }
            if (cleanButtonValues[i] == 1)
            {
                buttonActionsSFX[i][ActionHoldContinuous] = true;
                buttonActionsUI[i][ActionHoldContinuous] = true;
                writeButtonFlag = i;
                writeActionFlag = ActionHoldContinuous;
            }
            if (buttonHysteresis[i] < buttonHysteresisThreshold)
            {
                if (buttonCounters[i] < buttonHoldMax)
                    buttonCounters[i]++;
                if ((buttonCounters[i] >= buttonHoldThreshold) && (cleanButtonValues[i] == 1))
                {
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

            if (currentControl <= 0)
                currentControl = (BirlControlType) ((int) ControlNil - 1);
            else
                currentControl = (BirlControlType) ((int) currentControl - 1);

            loadingControl = 1;
            clearButtonActions();
        }

        if (buttonActionsUI[ButtonNextControl][ActionPress] == 1)
        {
            previousControl = currentControl;

            if (currentControl >= ControlNil - 1)
                currentControl = (BirlControlType) 0;
            else
                currentControl = (BirlControlType) ((int) currentControl + 1);

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
    void changeTuning()
    {
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
    //     void setKnobValues(float* values)
    //     {
    //         for (int i = 0; i < KNOB_PAGE_SIZE; i++) {
    //             int knob = i;
    //             if (knob + (knobPage * KNOB_PAGE_SIZE) == cvAddParam[currentControl]) {
    //                 knob = KNOB_PAGE_SIZE;
    //             }
    // //            knobActive[knob] = 0;
    //         }
    //     }
    //     void setKnobValue(int knob, float value)
    //     {
    //         if (knob + (knobPage * KNOB_PAGE_SIZE) == cvAddParam[currentControl]) {
    //             knob = KNOB_PAGE_SIZE;
    //         }
    // //            knobActive[knob] = 0;
    //
    //     }
    void deactivateKnob (int knob)
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
