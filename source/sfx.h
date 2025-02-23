#ifndef SFX_H
#define SFX_H


//#include "PluginEditor.h"
#include "ui.h"
#include "Filters.h"
#include "Tube.h"
#include "Tuning.h"
#include "tunings.h"
#include "juce_audio_basics/juce_audio_basics.h"
#include "leaf-physical.h"


#define SMALL_MEM_SIZE 80328
#define MED_MEM_SIZE 519000
#define LARGE_MEM_SIZE 33554432 // size of SDRAM IC




namespace birl {

extern char small_memory[SMALL_MEM_SIZE];
extern char medium_memory[MED_MEM_SIZE];
extern char large_memory[LARGE_MEM_SIZE];


extern tMempool smallPool;
extern tMempool largePool;

/* physical model */
extern Tube tubes[MAX_TONEHOLES + 1];
extern tPoleZero toneHoles[MAX_TONEHOLES];
//extern DCFilter *dcblocker1;
//extern DCFilter *dcblocker2;
extern tHighpass dcblocker1;
extern tHighpass dcblocker2;
extern tBiQuad biquad;
extern tSVF pf1;
extern tSVF pf2;
extern tSVF lp1;
extern tSVF lp2;
extern tSVF noiseBP;
extern tStereoRotation tubeRot[NUM_OF_TONEHOLES+1];
 extern   float pedal1;
extern    float pedal2;

extern float controlKnobValues[ControlNil][NUM_SYNTH_KNOB_VALUES];
extern uint8_t knobActive[NUM_SYNTH_KNOB_VALUES];
extern float fingers[NUM_OF_TONEHOLES];
extern float maxToneholeArg[NUM_OF_TONEHOLES];
extern bool buttons[NUM_OF_BUTTONS];

extern float breathArray[2];

//extern PlayMode samplerMode;
extern float sampleLength;

extern uint32_t freeze;

void initGlobalSFXObjects(LEAF &leaf);

/* physical model physical model */

void SFXPhysicalModelPMAlloc(LEAF &leaf);
void SFXPhysicalModelPMFrame(juce::AudioBuffer<float>& buffer);
void SFXPhysicalModelPMTick(float* input);
void SFXPhysicalModelPMFree(void);
void SFXPhysicalModelSetToneholeRadius(int index, float radius);
void SFXPhysicalModelSetTonehole(int index, float newValue);
void SFXPhysicalModelSetBreathPressure(float input);
void SFXPhysicalModelCalcTHCoeffs();
void SFXPhysicalModelTune(float fundamental);
void SFXPhysicalModelRetune(float fundamental);
float SFXPhysicalModelInterpolateLinear(float a, float b, float alpha);


/* slide birl */
void SFXRuleBasedPMAlloc();
void SFXRuleBasedPMFrame();
void SFXRuleBasedPMTick(float* input);
void SFXRuleBasedPMFree(void);

/* toy birl */
void SFXRuleBasedSynthAlloc();
void SFXRuleBasedSynthFrame();
void SFXRuleBasedSynthTick(float* input);
void SFXRuleBasedSynthFree(void);

/* neural net pm */
void SFXNeuralNetPMAlloc();
void SFXNeuralNetPMFrame();
void SFXNeuralNetPMTick(float* input);
void SFXNeuralNetPMFree(void);

/* neural net synth */
void SFXNeuralNetSynthAlloc();
void SFXNeuralNetSynthFrame();
void SFXNeuralNetSynthTick(float* input);
void SFXNeuralNetSynthFree(void);

/* MIDI functions */
void noteOn(int key, int velocity);
void noteOff(int key, int velocity);
void pitchBend(int data);
void sustainOn(void);
void sustainOff(void);
void toggleBypass(void);
void toggleSustain(void);


void calculateFreq(int voice);
float calculateTunedMidiNote(float tempNote); // for detuning?

void calculateNoteArray(void);
float nearestNote(float period);
float nearestNoteWithHysteresis(float note, float hysteresis);

void clearNotes(void);
void ctrlInput(int ctrl, int value);

}

#endif
