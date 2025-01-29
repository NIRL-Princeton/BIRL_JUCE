#ifndef UI_H_
#define UI_H_

#include "tunings.h"
#include "sfx.h"

#include <stdint.h>
namespace birl
{
    #define NUM_CONTROL_VALUES 5
    #define NUM_SYNTH_KNOB_VALUES 31
    #define NUM_XY_CONTROLS 4
    #define NUM_BUTTONS 9
    #define KNOB_PAGE_SIZE 5
    
    typedef enum _BirlControlType
    {
        PhysicalModelPM = 0,
        RuleBasedPM,
        RuleBasedSynth,
        NeuralNetPM,
        NeuralNetSynth,
        ControlNil
    } BirlControlType;

    typedef enum _BirlButton
    {
        ButtonOctaveUp = 0,
        ButtonOctaveDown,
        ButtonPrevControl,
        ButtonNextControl,
        ButtonA,
        ButtonB,
        ButtonC,
        ButtonD,
        ButtonE,
        ButtonNil
    } BirlButton;
    
    typedef enum _ButtonAction
    {
        ActionPress = 0,
        ActionRelease,
        ActionHoldInstant,
        ActionHoldContinuous,
        ActionNil
    } ButtonAction;
    
    extern tExpSmooth xy[NUM_XY_CONTROLS];
    extern uint16_t XY_values[NUM_XY_CONTROLS];
    extern float smoothedXY[NUM_XY_CONTROLS];
    extern uint8_t buttonValues[NUM_BUTTONS];
    
    extern int8_t writeKnobFlag;
    extern int8_t writeButtonFlag;
    extern int8_t writeActionFlag;

    
    extern BirlControlType currentControl;
    extern BirlControlType previousControl;

    extern uint8_t loadingControl;
    extern uint8_t loadingSynth;

    
    // Display values
    extern const char* controlNames[ControlNil];
    extern const char* controlNamesDetails[ControlNil];
    extern const char* shortControlNames[ControlNil];

    extern const char* knobParamNames[ControlNil][NUM_SYNTH_KNOB_VALUES];
    extern float displayValues[NUM_SYNTH_KNOB_VALUES];
    extern int8_t cvAddParam[ControlNil];
    extern uint8_t knobPage;
    extern uint8_t buttonActionsUI[NUM_BUTTONS+1][ActionNil];
    extern uint8_t buttonActionsSFX[NUM_BUTTONS+1][ActionNil];
    extern char* (*buttonActionFunctions[ControlNil])(BirlButton, ButtonAction);

    void initParamNames(void);
    void buttonCheck(void);
    void xyCheck(void);
    void clearButtonActions(void);
    void changeTuning(void);
    void writeCurrentPresetToFlash(void);
    void incrementPage(void);
    void decrementPage(void);
    void resetKnobValues(void);
    void setKnobValues(float* values);
    void setKnobValue(int knob, float value);
    void deactivateKnob(int knob);
    void deactivateAllKnobs(void);
    
    const char* UIConnectButton(BirlButton button, ButtonAction action);

    }


#endif

