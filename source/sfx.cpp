#include "sfx.h"
#include "leaf.h"
namespace birl
{

char small_memory[SMALL_MEM_SIZE];
char medium_memory[MED_MEM_SIZE];
char large_memory[LARGE_MEM_SIZE];


tMempool smallPool;
tMempool largePool;

float defaultControlKnobValues[ControlNil][NUM_SYNTH_KNOB_VALUES];
float controlKnobValues[ControlNil][NUM_SYNTH_KNOB_VALUES];
//uint8_t knobActive[NUM_ADC_CHANNELS];
float prevDisplayValues[NUM_SYNTH_KNOB_VALUES];


#define EXP_BUFFER_SIZE 128
float expBuffer[EXP_BUFFER_SIZE];
float expBufferSizeMinusOne = EXP_BUFFER_SIZE - 1;

#define DECAY_EXP_BUFFER_SIZE 512
float decayExpBuffer[DECAY_EXP_BUFFER_SIZE];
float decayExpBufferSizeMinusOne = DECAY_EXP_BUFFER_SIZE - 1;

float toneholes[NUM_OF_TONEHOLES];
float maxToneholeArg[NUM_OF_TONEHOLES];
bool buttons[NUM_OF_BUTTONS];

uint8_t dummyMemory[128];

float randomCall()
{
    return 0.0f;
}

void initGlobalSFXObjects(LEAF &leaf)
{
    LEAF_init(&leaf, 44100.0f, medium_memory, MED_MEM_SIZE, [](){return (float)rand()/RAND_MAX;});
    tMempool_init (&smallPool, small_memory, SMALL_MEM_SIZE, &leaf);
    tMempool_init (&largePool, large_memory, LARGE_MEM_SIZE, &leaf);

    LEAF_generate_exp(expBuffer, 1000.0f, -1.0f, 0.0f, -0.0008f, EXP_BUFFER_SIZE);
    // exponential buffer rising from 0 to 1
    LEAF_generate_exp(decayExpBuffer, 0.001f, 0.0f, 1.0f, -0.0008f, DECAY_EXP_BUFFER_SIZE);
    // exponential decay buffer falling from 1 to 0
    
    defaultControlKnobValues[PhysicalModelPM][0] = 0.2f;        // gain
    defaultControlKnobValues[PhysicalModelPM][1] = 100.0f;      // fundamental
    defaultControlKnobValues[PhysicalModelPM][2] = 9.0f;        // num_toneholes
    defaultControlKnobValues[PhysicalModelPM][3] = 0.995f;      // dcblocker1
    defaultControlKnobValues[PhysicalModelPM][4] = 0.995f;      // dcblocker2
    defaultControlKnobValues[PhysicalModelPM][5] = 0.169301f;   // biquad_coeff1
    defaultControlKnobValues[PhysicalModelPM][6] = 0.338601f;   // biquad_coeff2
    defaultControlKnobValues[PhysicalModelPM][7] = 0.169301f;   // biquad_coeff3
    defaultControlKnobValues[PhysicalModelPM][8] = -0.482013f;  // biquad_coeff4
    defaultControlKnobValues[PhysicalModelPM][9] = 0.186622f;   // biquad_coeff5
    defaultControlKnobValues[PhysicalModelPM][10] = 2000.0f;    // pf1_cutoff
    defaultControlKnobValues[PhysicalModelPM][11] = 0.5f;       // pf1_q
    defaultControlKnobValues[PhysicalModelPM][12] = 1000.0f;    // pf2_cutoff
    defaultControlKnobValues[PhysicalModelPM][13] = 1.0f;       // pf2_q
    defaultControlKnobValues[PhysicalModelPM][14] = 5000.0f;    // lp1_cutoff
    defaultControlKnobValues[PhysicalModelPM][15] = 0.5f;       // lp1_q
    defaultControlKnobValues[PhysicalModelPM][16] = 5000.0f;    // lp2_cutoff
    defaultControlKnobValues[PhysicalModelPM][17] = 0.5f;       // lp2_q
    defaultControlKnobValues[PhysicalModelPM][18] = 0.0f;       // shaper drive
    defaultControlKnobValues[PhysicalModelPM][19] = 0.0f;       // shaper mix
    defaultControlKnobValues[PhysicalModelPM][20] = 0.2f;       // noise gain
    defaultControlKnobValues[PhysicalModelPM][21] = 16000.0f;   // noise_bp_cutoff
    defaultControlKnobValues[PhysicalModelPM][22] = 1.0f;       // noise_bp_q
    
    defaultControlKnobValues[RuleBasedSynth][0] = 1.0f;
    defaultControlKnobValues[RuleBasedSynth][1] = 1.0f;
    defaultControlKnobValues[RuleBasedSynth][2] = 0.5f;
    defaultControlKnobValues[RuleBasedSynth][3] = 0.0f;
    defaultControlKnobValues[RuleBasedSynth][4] = 0.0f;
    defaultControlKnobValues[RuleBasedSynth][5] = 0.0f;
    defaultControlKnobValues[RuleBasedSynth][6] = 1.0f;
    defaultControlKnobValues[RuleBasedSynth][7] = 0.0f;
    defaultControlKnobValues[RuleBasedSynth][8] = 0.0f;
    defaultControlKnobValues[RuleBasedSynth][9] = 22000.0f;
    defaultControlKnobValues[RuleBasedSynth][10] = 0.7f;
    defaultControlKnobValues[RuleBasedSynth][11] = 1.0f;
    defaultControlKnobValues[RuleBasedSynth][12] = 0.0f;
    defaultControlKnobValues[RuleBasedSynth][13] = 0.0f;
    defaultControlKnobValues[RuleBasedSynth][14] = 1.0f;
    defaultControlKnobValues[RuleBasedSynth][15] = 0.9f;
    
    defaultControlKnobValues[RuleBasedSynth][16] = 1.0f;
    defaultControlKnobValues[RuleBasedSynth][17] = 2.0f;
    defaultControlKnobValues[RuleBasedSynth][18] = 0.5f;
    defaultControlKnobValues[RuleBasedSynth][19] = 0.0f;
    defaultControlKnobValues[RuleBasedSynth][20] = 0.0f;
    defaultControlKnobValues[RuleBasedSynth][21] = 0.0f;
    defaultControlKnobValues[RuleBasedSynth][22] = 0.0f;
    defaultControlKnobValues[RuleBasedSynth][23] = 0.0f;
    defaultControlKnobValues[RuleBasedSynth][24] = 22000.0f;
    defaultControlKnobValues[RuleBasedSynth][25] = 0.7f;
    defaultControlKnobValues[RuleBasedSynth][26] = 0.0f;
    defaultControlKnobValues[RuleBasedSynth][27] = 0.0f;
    defaultControlKnobValues[RuleBasedSynth][28] = 0.0f;
    defaultControlKnobValues[RuleBasedSynth][29] = 0.0f;
    defaultControlKnobValues[RuleBasedSynth][30] = 0.9f;
    
    for (int c = 0; c < ControlNil; c++)
    {
        for (int v = 0; v < sizeof(controlKnobValues[c]); v++) {
            controlKnobValues[c][v] = defaultControlKnobValues[c][v];
        }
    }
}

/* physical model internal poly */
Tube *tubes[MAX_TONEHOLES - 1];
tPoleZero toneHoles[MAX_TONEHOLES];
//DCFilter *dcblocker1;
//DCFilter *dcblocker2;
tHighpass dcblocker1;
tHighpass dcblocker2;
tBiQuad biquad;
tSVF pf1;
tSVF pf2;
tSVF lp1;
tSVF lp2;
tSVF noiseBP;

float pam_[MAX_TONEHOLES];
float pbp_[MAX_TONEHOLES];
float pthp_[MAX_TONEHOLES];
float scatter_[MAX_TONEHOLES];
float thCoeff_[MAX_TONEHOLES];
float rth_[MAX_TONEHOLES];
float tubeLengths_[MAX_TONEHOLES - 1];

float rb;
float outputGain;
float noiseGain;
float breathPressure;
float prevBreathPressure;
float breathArray[2];

float fundamental;
int tubeIndex;
int numToneholes;

float count;
float min;
float max;
float mDrive;
float shaperMix;

void SFXPhysicalModelPMAlloc(LEAF &leaf)
{
    for (int i = 0; i < MAX_TONEHOLES; i++) {
        tPoleZero_initToPool(&toneHoles[i], &smallPool);
    }
    
    outputGain = defaultControlKnobValues[PhysicalModelPM][0];
    noiseGain = defaultControlKnobValues[PhysicalModelPM][20];
    
    fundamental = defaultControlKnobValues[PhysicalModelPM][1];
    tubeIndex = 0;
    numToneholes = defaultControlKnobValues[PhysicalModelPM][2];
    
    breathPressure = 0.0;
    prevBreathPressure = 0.0;
    count = 0;
    min = 0.0;
    max = 0.0;
    mDrive = 0.0;
    shaperMix = 0.0;
    
    for (int i = 0; i < MAX_TONEHOLES - 1; i++) {
//        tubes[i] = initTube(3); // IDK???
        tLinearDelay_init(&tubes[i]->upper, 100, 512, &leaf);
        tLinearDelay_init(&tubes[i]->lower, 100, 512, &leaf);
    }
//    dcblocker1 = initDCFilter(defaultControlKnobValues[PhysicalModelPM][3]);
//    dcblocker2 = initDCFilter(defaultControlKnobValues[PhysicalModelPM][4]);
    tHighpass_initToPool(&dcblocker1, 13.0, &smallPool);
    tHighpass_initToPool(&dcblocker2, 13.0, &smallPool);
    // choice of 13 Hz for cutoff due to translation from DC Blocker to LEAF.
    
    tBiQuad_initToPool(&biquad, &smallPool);
    tBiQuad_setCoefficients(biquad, defaultControlKnobValues[PhysicalModelPM][5], defaultControlKnobValues[PhysicalModelPM][6], defaultControlKnobValues[PhysicalModelPM][7], defaultControlKnobValues[PhysicalModelPM][8], defaultControlKnobValues[PhysicalModelPM][9]);
    
    tSVF_initToPool(&pf1,     SVFTypePeak,      defaultControlKnobValues[PhysicalModelPM][10], defaultControlKnobValues[PhysicalModelPM][11], &smallPool);
    tSVF_initToPool(&pf2,     SVFTypePeak,      defaultControlKnobValues[PhysicalModelPM][12], defaultControlKnobValues[PhysicalModelPM][13], &smallPool);
    tSVF_initToPool(&lp1,     SVFTypeLowpass,   defaultControlKnobValues[PhysicalModelPM][14], defaultControlKnobValues[PhysicalModelPM][15], &smallPool);
    tSVF_initToPool(&lp2,     SVFTypeLowpass,   defaultControlKnobValues[PhysicalModelPM][16], defaultControlKnobValues[PhysicalModelPM][17], &smallPool);
    tSVF_initToPool(&noiseBP, SVFTypeBandpass,  defaultControlKnobValues[PhysicalModelPM][21], defaultControlKnobValues[PhysicalModelPM][22], &smallPool);
}

void SFXPhysicalModelSetToneholeRadius(int index, float radius) {
    if (radius < MIN_TONEHOLE_RADIUS || radius > MAX_TONEHOLE_RADIUS) {
        juce::AlertWindow::showMessageBoxAsync(juce::AlertWindow::WarningIcon,
            "Radius out of range",
            "Error: Radius is out of range",
            "OK");
        return;
    }
    rth_[index] = radius;
    scatter_[index] = -pow(radius,2) / ( pow(radius,2) + 2*pow(rb,2) );

    // Calculate toneHole coefficients.
    double te = radius;
    thCoeff_[index] = (te*2*(SRATE*OVERSAMPLE) - C_m) / (te*2*(SRATE*OVERSAMPLE) + C_m);
}
void SFXPhysicalModelSetTonehole(int index, float newValue) {
    double new_coeff;
    if (newValue <= 0.0)
        new_coeff = 0.9995;
    else if (newValue >= 1.0)
        new_coeff = thCoeff_[index];
    else
        new_coeff = newValue * (thCoeff_[index] - 0.9995) + 0.9995;
    tPoleZero_setA1(toneHoles[index], -new_coeff);
    tPoleZero_setB0(toneHoles[index], new_coeff);
}
void SFXPhysicalModelSetBreathPressure(float input) {
    breathPressure = input;
}
void SFXPhysicalModelCalcTHCoeffs() {
    // Calculate the initial tonehole three-port scattering coefficient.
    scatter_[0] = -pow(rth_[0],2) / ( pow(rth_[0],2) + 2*pow(rb,2) );
    
    // Calculate tonehole coefficients and set for initially open.
    float te = 1.4 * rth_[0]; // effective length of the open hole
    thCoeff_[0] = (te*2*(SRATE*OVERSAMPLE) - C_m) / (te*2*(SRATE*OVERSAMPLE) + C_m);
//    thCoeff_[0] = (rth_[0]*2*(SRATE*OVERSAMPLE) - C_m) / (rth_[0]*2*(SRATE*OVERSAMPLE) + C_m);

    // Initialize toneholes.
    tPoleZero_setA1(toneHoles[0], -thCoeff_[0]);
    tPoleZero_setB0(toneHoles[0], thCoeff_[0]);
    tPoleZero_setB1(toneHoles[0], -1.0);
}
void SFXPhysicalModelTune(float fundamental) {
    float effectiveLength = calcLS(fundamental);
    int cutLength = calcLC(effectiveLength);
    double boreDiameter = calcd1(cutLength, effectiveLength);
    
    tubeLengths_[0] = calclH(0, boreDiameter, calcLSh(0, fundamental));
//    C.O. : FINAL TUBE USES FRACTIONAL DELAY LENGTH.
//    tubeLengths_[1] = calclH(1, boreDiameter, calcLSh(1, fundamental));
    tubeLengths_[1] = calcLSh(1, fundamental);

    
    if (tubeLengths_[0] == 0 || tubeLengths_[1] == 0) {
        juce::AlertWindow::showMessageBoxAsync(juce::AlertWindow::WarningIcon,
            "Int Cast Error",
            "Error: Cut length cast down to 0. Use a different tuning or try oversampling.",
            "OK");
        return;
    }
    // Initialize tonehole.
    tPoleZero_initToPool(&toneHoles[0], &smallPool);

    // Initialize tube A.
//    tubes[0] = initTube(tubeLengths_[0]);
    
    tLinearDelay_setDelay(tubes[0]->upper, tubeLengths_[0]);
    tLinearDelay_setDelay(tubes[0]->lower, tubeLengths_[0]);
    
    // Initialize tube B.
//    tLinearDelay_init (&endTube_[0], tubeLengths_[1], (int)(tubeLengths_[1]+1));
//    tLinearDelay_init (&endTube_[1], tubeLengths_[1], (int)(tubeLengths_[1]+1));
    
    // Main bore radius, in centimeters.
    rb = boreDiameter / 2.0;
    
    // Calculate the toneHole radii.
    rth_[0] = 0.5 * calcdH(0, boreDiameter, tubeLengths_[0], fundamental);
    
    // Calculate the tonehole coefficients.
    SFXPhysicalModelCalcTHCoeffs();
}
//void SFXPhysicalModelRetune(float fundamental) {
//    double effectiveLength = calcLS(fundamental);
//    int cutLength = calcLC(effectiveLength);
//    double boreDiameter = calcd1(cutLength, effectiveLength);
//
//    double old_upper = *(tubes[0]->upper->data);
//    double old_lower = *(tubes[0]->lower->data);
//
//
//    tubeLengths_[0] = calclH(0, boreDiameter, calcLSh(0, fundamental));
////    tubeLengths_[0] = calcLSh(0, f);
//
//    freeTube(tubes[0]);
////    freeFracTube(ftubes_[0]);
//
//    tubes[0] = initTube(tubeLengths_[0]);
////    ftubes_[0] = initFracTube(tubeLengths_[0]);
//
//
//    inputDelayLine(tubes[0]->upper, old_upper);
//    inputDelayLine(tubes[0]->lower, old_lower);
////    tLinearDelay_tickIn(&(ftubes_[0]->upper), old_upper);
////    tLinearDelay_tickIn(&(ftubes_[0]->lower), old_lower);
//    rb = boreDiameter / 2.0;
//}
float SFXPhysicalModelInterpolateLinear(float a, float b, float alpha) {
    return (alpha * a) + ((1.0-alpha) * b);
}

void SFXPhysicalModelPMFrame() {

}


void SFXPhysicalModelPMTick(float* input) {
    float sample = 0.0f;
    float bellReflected;

//    double pap;
//    double pbm;
//    double pthm;
//    double scatter;
//
    float breath = breathPressure;
    float noise = (double) rand() / (double) RAND_MAX;
    
//    noise = noiseGain * (inputSVFBand(noiseBP, noise));
    noise = noiseGain * tSVF_tick(noiseBP, noise);
    breath += breath * noise;
    
    float pressureDiff = tLinearDelay_tickOut(tubes[0]->lower) - breath;
//    float pressureDiff = accessDelayLine(&tubes[0]->lower) - breath;
//    double pressureDiff = tLinearDelay_tickOut(&(ftubes_[0]->lower)) - breath;
    float reedLookup = pressureDiff * reedTable (pressureDiff);
    
    breath = LEAF_clip(-1.0, breath + reedLookup, 1.0);
//    breath = interpolateLinear(shaper(breath, m_drive), breath, shaperMix);

    breath = tSVF_tick(pf1, breath);
    breath = tSVF_tick(lp1, breath);
    breath = tHighpass_tick(dcblocker1, breath);
//    breath = inputSVFPeak(pf_, breath);
//    breath = inputSVFLP(lp_, breath);
//    breath = inputDCFilter(dcBlocker_, breath);
    

    
    /* bell filters */
    
    float bell = tLinearDelay_tickOut(tubes[0]->upper);
//    float bell = accessDelayLine(tubes[0]->upper);
//    double bell = tLinearDelay_tickOut(&(ftubes_[0]->upper));
    
    // Reflection = Inversion + gain reduction + lowpass filtering.
    bell = tSVF_tick(pf2, bell);
    bell = tSVF_tick(lp2, bell);
    bell = tHighpass_tick(dcblocker2, bell);
//    bell = inputSVFLP(lp2_, bell);
//    bell = inputDCFilter(dcBlocker2, bell);
    
    sample = bell;
    bellReflected = bell * -0.995;

    tLinearDelay_tickIn(tubes[0]->upper, breath);
    tLinearDelay_tickIn(tubes[0]->lower, bellReflected);

//    inputDelayLine(tubes[0]->upper, breath);
//    inputDelayLine(tubes[0]->lower, bellReflected);
    
//    tLinearDelay_tickIn(&(ftubes_[0]->upper), breath);
//    tLinearDelay_tickIn(&(ftubes_[0]->lower), bellReflected);
//    tanhClip(outsamp);
    
    sample = tanhf(sample);
    input[0] = sample;
    input[1] = sample;
}

void SFXPhysicalModelPMFree(void) {
    for (int i = 0; i < MAX_TONEHOLES; i++) {
        tPoleZero_free(&toneHoles[i]);
    }
    for (int i = 0; i < MAX_TONEHOLES - 1; i++)
        freeTube(tubes[i]);
    tHighpass_free(&dcblocker1);
    tHighpass_free(&dcblocker2);
    tBiQuad_free(&biquad);
    tSVF_free(&pf1);
    tSVF_free(&pf2);
    tSVF_free(&lp1);
    tSVF_free(&lp2);
    tSVF_free(&noiseBP);
}

/* slide birl */
void SFXRuleBasedPMAlloc() {
    
}
void SFXRuleBasedPMFrame() {}
void SFXRuleBasedPMTick(float* input) {}
void SFXRuleBasedPMFree(void) {}

/* toy birl */
#define NUM_OF_OSCILLATORS 2
#define NUM_OF_FILTERS 2
#define NUM_OF_SMOOTHERS 3

tCycle sine[NUM_OF_OSCILLATORS];
tSawtooth saw[NUM_OF_OSCILLATORS];
tSquare square[NUM_OF_OSCILLATORS];
tTriangle triangle[NUM_OF_OSCILLATORS];
tEnvelope envelopes[NUM_OF_OSCILLATORS];

//tSVF filters[NUM_OF_FILTERS];
tSVF lowpass[NUM_OF_FILTERS];
tSVF bandpass[NUM_OF_FILTERS];
tSVF highpass[NUM_OF_FILTERS];
tSVF notch[NUM_OF_FILTERS];
tSVF peak[NUM_OF_FILTERS];

tExpSmooth smoothers[NUM_OF_SMOOTHERS];
tSVF noise;

float sample[NUM_OF_OSCILLATORS];
float frequency[NUM_OF_OSCILLATORS];
int currentOctave = 0;
int prevOctave = 0;


float openholeMIDI = 53.0f;
float octaveTransposition = 0.0f;

void SFXRuleBasedSynthAlloc() {

    for (int i = 0; i < NUM_OF_OSCILLATORS; i++) {
        tCycle_initToPool(&sine[i], &smallPool);
        tSawtooth_initToPool(&saw[i], &smallPool);
        tSquare_initToPool(&square[i], &smallPool);
        tTriangle_initToPool(&triangle[i], &smallPool);
        tEnvelope_initToPool(&envelopes[i], 0.0f, 0.0f, 0, &smallPool); // loop?
    }
    
    for (int i = 0; i < NUM_OF_FILTERS; i++) {
        tSVF_initToPool(&lowpass[i], SVFTypeLowpass, 22000.0, 0.7f, &smallPool); // default
        tSVF_initToPool(&bandpass[i], SVFTypeBandpass, 3000.0, 0.7f, &smallPool); // default
        tSVF_initToPool(&highpass[i], SVFTypeHighpass, 0.0, 0.7f, &smallPool); // default
        tSVF_initToPool(&notch[i], SVFTypeNotch, 110.0, 0.7f, &smallPool); // default
        tSVF_initToPool(&peak[i], SVFTypePeak, 110.0, 0.7f, &smallPool); // default
    }
    
    for (int i = 0; i < NUM_OF_SMOOTHERS; i++) {
        tExpSmooth_initToPool(&smoothers[i], 0.0f, 0.02f, &smallPool); // 20 ms default
    }

    tSVF_initToPool(&noise, SVFTypeBandpass, 10000.0f, 0.3f, &smallPool);
}
void SFXRuleBasedSynthFrame() {

    if (buttons[ButtonOctaveUp] == 1) {
        octaveTransposition = 12.0;
    } else if (buttons[ButtonOctaveDown] == 1){
        octaveTransposition = -12.0;
    } else {
        octaveTransposition = 0.0;
    }

//    if (buttonActionsSFX[ButtonOctaveUp][ActionHoldContinuous] == 1) {
//        ++currentOctave;
//        buttonActionsSFX[ButtonOctaveUp][ActionPress] = 0;
//    }
//    if (buttonActionsSFX[ButtonOctaveDown][ActionPress] == 1) {
//        --currentOctave;
//        buttonActionsSFX[ButtonOctaveDown][ActionPress] = 0;
//    }
//    displayValues[0] = controlKnobValues[RuleBasedSynth][0]; // osc on

    float midiAdjustment = 12.0;
    if (toneholes[0] > 0.0)
    {
      midiAdjustment -= toneholes[0];
      if (toneholes[0] == 1.0)
      {
          midiAdjustment -= 2 * toneholes[1];
          if (toneholes[1] == 1.0)
          {
              midiAdjustment -= 2 * toneholes[2];
          }
          
      }
      if (toneholes[1] == 1.0 && toneholes[2] == 1.0)
      {
          if (toneholes[4] > 0.0)
          {
              midiAdjustment -= 2 * toneholes[4];
              if (toneholes[4] == 1.0)
              {
                  midiAdjustment -= toneholes[5];
                  if (toneholes[5] == 1.0)
                  {
                      midiAdjustment -= 2 * toneholes[6];
                      if (toneholes[6] == 1.0)
                      {
                          midiAdjustment -= toneholes[7];
                          midiAdjustment -= 2 * toneholes[8];
                      }
                  }
              }
          }
          if (toneholes[4] == 0.0 && toneholes[5] > 0.0)
          {
              midiAdjustment -= toneholes[5];
          }
      }
    }
    else if (toneholes[1] > 0.0)
    {
      midiAdjustment -= toneholes[1];
    }
    midiAdjustment += toneholes[3];

    float floor;
    float decimals = modf(midiAdjustment + openholeMIDI + octaveTransposition, &floor);

    
    // calculate pitch class (0-11) from intpart
    int pitchclass = (int) floor % 12 - 5;
//        DBG("pitchclass = " + (String)pitchclass);

    // declare our notes
    float note;
    float noteAnotherFinger;
    
    /* multiple the int part by the proper scale
     * (based on current base pitch and scale choice)
     * if no scale choice, default to equal-tempered
     */
    
    if (buttons[ButtonPrevControl] == 1) {
          // transpose down
      }
    if (buttons[ButtonNextControl] == 1) {
        // transpose up
    }
    
    if (buttons[ButtonA] == 1) {
        ++currentTuning;
        buttons[ButtonA] = 0;
    }
    
    if (buttons[ButtonB]== 1) {
        --currentTuning;
        buttons[ButtonB] = 0;
    }
    
    if (buttons[ButtonC] == 1) {
        // left free
    }
    
    if (buttons[ButtonD] == 1) {
        // left free
    }
    
    if (buttons[ButtonE] == 1) {
        // left free
    }
    
    note = tuningPresets[currentTuning][pitchclass];
    noteAnotherFinger = tuningPresets[currentTuning][pitchclass-1];
    
    float pitchinMIDI = floor + decimals * (note - noteAnotherFinger);
    float pitchinFreq = LEAF_midiToFrequency(pitchinMIDI);
    
    for (int i = 0; i < NUM_OF_OSCILLATORS; i++) {
        frequency[i] = pitchinFreq;
    }
    

    // frequency = LEAF_midiToFrequency(floor + decimals * (note - noteAnotherFinger));
    
//        DBG("midiAdjustment = " + (String)midiAdjustment
//            + ", note = " + (String)note
//            + ", noteAnotherFinger = " + (String)noteAnotherFinger);
//        DBG("processor.frequency = " + (String) processor.frequency);
    
}


void SFXRuleBasedSynthTick(float* input) {
    

    if (controlKnobValues[2][0] == 1.0) { // if the first oscillator is on

        // manage octave transposition for osc1
        float osc1_freq_midi = LEAF_frequencyToMidi(frequency[0]);
        float osc1_octave = controlKnobValues[2][3]; // this is gonna be a value from -3.0 to 3.0, default is 0
        osc1_freq_midi = osc1_freq_midi + 12.0 * osc1_octave;
        frequency[0] = LEAF_midiToFrequency(osc1_freq_midi);
        
        // semi/detune for osc1 * needs tweaking!!!!
        float osc1_semi = controlKnobValues[2][4];
        float osc1_detuned = controlKnobValues[2][5]; // value from -3.0 to 3.0, default is 0
        osc1_freq_midi = LEAF_frequencyToMidi(frequency[0]) + osc1_semi + osc1_detuned * (float)(1.0/3.0f);
        frequency[0] = LEAF_midiToFrequency(osc1_freq_midi);
        
        // manage octave transposition for osc2
        float osc2_freq_midi = LEAF_frequencyToMidi(frequency[1]);
        float osc2_octave = controlKnobValues[2][19];
        osc2_freq_midi = osc2_freq_midi + 12.0 * osc2_octave;
        frequency[1] = LEAF_midiToFrequency(osc2_freq_midi);
        
        // semi/detune for osc2 * needs tweaking !!!
        float osc2_semi = controlKnobValues[2][20];
        float osc2_detuned = controlKnobValues[2][21]; // value from -3.0 to 3.0, default is 0
        osc2_freq_midi = LEAF_frequencyToMidi(frequency[1]) + osc2_semi + osc2_detuned * (float)(1.0/3.0f);
        frequency[1] = LEAF_midiToFrequency(osc1_freq_midi);
        
    
        switch ((int) controlKnobValues[2][1]) { // osc 1 waveform
            case 0: // sine
                tCycle_setFreq(sine[0], frequency[0]);
                sample[0] = tCycle_tick(sine[0]);
                break;
            case 1: // saw
                tSawtooth_setFreq(saw[0], frequency[0]);
                sample[0] = tSawtooth_tick(saw[0]);
                break;
            case 2: // square
                tSquare_setFreq(square[0], frequency[0]);
                sample[0] = tSquare_tick(square[0]);
                break;
            case 3: // triangle
                tTriangle_setFreq(triangle[0], frequency[0]);
                sample[0] = tTriangle_tick(triangle[0]);
                break;
            default:
                break;
        }
    } // end osc1
    
    if (controlKnobValues[2][16] == 1) {
        switch ((int) controlKnobValues[2][17]) { // osc 2 waveform
            case 0: // sine
                tCycle_setFreq(sine[1], frequency[1]);
                sample[1] = tCycle_tick(sine[1]);
                break;
            case 1: // saw
                tSawtooth_setFreq(saw[1], frequency[1]);
                sample[1] = tSawtooth_tick(saw[1]);
                break;
            case 2: // square
                tSquare_setFreq(square[1], frequency[1]);
                sample[1] = tSquare_tick(square[1]);
                break;
            case 3: // triangle
                tTriangle_setFreq(triangle[1], frequency[1]);
                sample[1] = tTriangle_tick(triangle[1]);
                break;
            default:
                break;
        }
    } // end osc2
    

    /* FILTERS */
    
    float f1_to_f2 = controlKnobValues[2][8];
    float f1_freq = controlKnobValues[2][9];
    float f1_q = controlKnobValues[2][10];
    
    if (controlKnobValues[2][6] == 1.0) { // if filter 1 is on
        switch ((int) controlKnobValues[2][7]) {
            case 0: // lowpass
                tSVF_setFreqAndQ(lowpass[0], f1_freq, f1_q);
                sample[0] = tSVF_tick(lowpass[0], sample[0]);
                break;
            case 1: // bandpass
                tSVF_setFreqAndQ(bandpass[0], f1_freq, f1_q);
                sample[0] = tSVF_tick(bandpass[0], sample[0]);
                break;
            case 2: // highpass
                tSVF_setFreqAndQ(highpass[0], f1_freq, f1_q);
                sample[0] = tSVF_tick(highpass[0], sample[0]);
                break;
            case 3: // notch
                tSVF_setFreqAndQ(notch[0], f1_freq, f1_q);
                sample[0] = tSVF_tick(notch[0], sample[0]);
                break;
            case 4: // peak
                tSVF_setFreqAndQ(peak[0], f1_freq, f1_q);
                sample[0] = tSVF_tick(peak[0], sample[0]);
                break;
            default:
                break;
        }
    } // end filter 1
    
    float f2_freq = controlKnobValues[2][24];
    float f2_q = controlKnobValues[2][25];
    
    if (controlKnobValues[2][22] == 1.0) { // if filter 2 is on
        switch ((int) controlKnobValues[2][23]) {
            case 0: // lowpass
                tSVF_setFreqAndQ(lowpass[1], f2_freq, f2_q);
                sample[1] = tSVF_tick(lowpass[1], sample[1]);
                break;
            case 1: // bandpass
                tSVF_setFreqAndQ(bandpass[1], f2_freq, f2_q);
                sample[1] = tSVF_tick(bandpass[1], sample[1]);
                break;
            case 2: // highpass
                tSVF_setFreqAndQ(highpass[1], f2_freq, f2_q);
                sample[1] = tSVF_tick(highpass[1], sample[1]);
                break;
            case 3: // notch
                tSVF_setFreqAndQ(notch[1], f2_freq, f2_q);
                sample[1] = tSVF_tick(notch[1], sample[1]);
                break;
            case 4: // peak
                tSVF_setFreqAndQ(peak[1], f2_freq, f2_q);
                sample[1] = tSVF_tick(peak[1], sample[1]);
                break;
            default:
                break;
        }
    } // end filter 2
    

    /* GAIN */
    
    float master_gain = controlKnobValues[0][0];
    
    float osc1_gain = controlKnobValues[2][2]; 
    float osc2_gain = controlKnobValues[2][18];

    tExpSmooth_setDest(smoothers[0], osc1_gain);
    sample[0] *= tExpSmooth_tick(smoothers[0]) * 0.5f;
    tExpSmooth_setDest(smoothers[1], osc2_gain);
    sample[1] *= tExpSmooth_tick(smoothers[1]) * 0.5f;
    
    float outsample = sample[0] + sample[1];
    
//    tExpSmooth_setDest(&smoother1, controlKnobValues[0][0]);
//    float smootherGain = tExpSmooth_tick(&smoother1);
//    float filterFrequency = LEAF_clip(30.0, smootherGain * 15 * 400, 17000.0);

    /* NOISE */
    float synth_noise_gain = controlKnobValues[0][20];
    float synth_noise_cutoff = controlKnobValues[0][21];
    float synth_noise_q = controlKnobValues[0][22];
    tSVF_setFreqAndQ(noise, synth_noise_cutoff, synth_noise_q);
    outsample += tSVF_tick(noise, outsample) * synth_noise_gain;
    outsample *= master_gain;
    
    outsample = tanhf(outsample);
    input[0] = outsample;
    input[1] = outsample;
}
void SFXRuleBasedSynthFree(void) {
    

    for (int i = 0; i < NUM_OF_OSCILLATORS; i++) {
        tCycle_free(&sine[i]);
        tSawtooth_free(&saw[i]);
        tSquare_free(&square[i]);
        tTriangle_free(&triangle[i]);
        tEnvelope_free(&envelopes[i]); // loop?
    }
    
    for (int i = 0; i < NUM_OF_FILTERS; i++) {
        tSVF_free(&lowpass[i]); // default
        tSVF_free(&bandpass[i]);
        tSVF_free(&highpass[i]);
        tSVF_free(&notch[i]);
        tSVF_free(&peak[i]);
    }
    
    for (int i = 0; i < NUM_OF_SMOOTHERS; i++) {
        tExpSmooth_free(&smoothers[i]); // 20 ms default
    }
    tSVF_free(&noise);
}



/* neural net pm */
void SFXNeuralNetPMAlloc() {}
void SFXNeuralNetPMFrame() {}
void SFXNeuralNetPMTick(float* input) {}
void SFXNeuralNetPMFree(void) {}

/* neural net synth */
void SFXNeuralNetSynthAlloc() {}
void SFXNeuralNetSynthFrame() {}
void SFXNeuralNetSynthTick(float* input) {}
void SFXNeuralNetSynthFree(void) {}

/* MIDI functions */
void noteOn(int key, int velocity) {}
void noteOff(int key, int velocity) {}
void pitchBend(int data) {}
void sustainOn(void) {}
void sustainOff(void) {}
void toggleBypass(void) {}
void toggleSustain(void) {}

void calculateFreq(int voice) {}

//float calculateTunedMidiNote(float tempNote) {} // for detuning?

//void calculateNoteArray(void) {}
//float nearestNote(float period) {}
//float nearestNoteWithHysteresis(float note, float hysteresis) {}

void clearNotes(void) {}
void ctrlInput(int ctrl, int value) {}

    
    
}
