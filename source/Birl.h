#ifndef BIRL_H_
#define BIRL_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const double CLIP_RATIO = 0.7;
const int BIRL_PORT = 1234;
#define NUM_OF_TONEHOLES 9
#define NUM_OF_BUTTONS 9
#define MAX_TONEHOLES 9
const int FRONT_TUBES = 1;
//const int MAX_TONEHOLES = NUM_OF_KEYS;
//const int MAX_TUBES = NUM_OF_KEYS + FRONT_TUBES;
const int SAMPLE_INDEX = 0;
const int MAX_TUBE_LENGTH = 100;
const double MIN_TONEHOLE_RADIUS = 0.0001;
const double MAX_TONEHOLE_RADIUS = 0.004;
const double RB_TWEAK_FACTOR = 0.0001;
const double RTH_TWEAK_FACTOR = 0.005;

static double originalRth_[MAX_TONEHOLES];
static double rth_[MAX_TONEHOLES];
static double tubeLengths_[MAX_TONEHOLES+1];


const int OVERSAMPLE = 1;

const double C_cm = 34500.0;
const double C_m = 345.0;
static double sRate_ = 441000; // sample rate is actually 48000 i think?
const double reedTableOffset = 0.7;
const double reedTableSlope = -0.3;
/* const double MIN_D1 = 0.05; */

// c major with Eb and Bb
// const float targetFrequencies[10] = {640.99, 553.89, 495.93, 467.40, 415.63, 375.86, 334.92, 316.51, 284.78, 262.00};
// const float targetFres[10] = {490.929557, 427.466852, 382.329735, 359.846659, 320.182398, 289.295447, 257.733345, 243.083612, 218.831734,200.0000};
#endif

