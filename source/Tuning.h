#ifndef TUNING_H_
#define TUNING_H_

#include "Birl.h"
#include <juce_core/juce_core.h>
#include <juce_gui_basics/juce_gui_basics.h>
// ONEHOLE. TWO NOTES.
const int NUM_NOTES = 10;

// CONSTS.
const double MIN_D1 = 0.5;
const double MIN_DH = 0.03;
const double DH_FIRST_GUESS = 0.03;

// = {2.519840, 2.244920, 2.000000, 1.887750, 1.681790, 1.498310, 1.334830, 1.259920, 1.122460, 1.000000, 0.9439};
const double tuning[] = {10.0/4.0, 18.0/8.0, 2.0/1.0, 15.0/8.0, 5.0/3.0, 3.0/2.0, 4.0/3.0, 5.0/4.0, 9.0/8.0, 1.0, 15.0/16.0};

static inline double convertTocm(double samps) {
    return (samps * (C_cm / (SRATE*OVERSAMPLE)));
}

static double convertToSamples(double cm) {
    return (cm * ((SRATE*OVERSAMPLE) / C_cm));
}


static double calcg(int thNum) {
    // CHANGED THIS FROM NUM_NOTES - 1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (thNum < 0 || thNum > NUM_NOTES) {
        printf("calcg: thNum out of bounds: %d\n", thNum);
        return 0.0;
    }
    return (tuning[thNum] / tuning[thNum+1]) - 1.0;
}

static double calcd1(int LC, double LS) {
    return (LS - (double) LC)/(0.3 * convertToSamples(2.54));
    /* return (LS - (double) LC)/0.3; */
}

// In samples.
static int calcLC(double LS) {
    double LC = (int) LS;
    double d1 = calcd1(LC, LS);
    while (d1 < MIN_D1) {
        printf("calcLC: d1 = %f for this value of LC so we're shortening LC!!!!!!\n", d1);
        LC -= 1.0;
        d1 = calcd1(LC, LS);
    }
    return LC;
}

// In samples.
static double calcLS(double Fc) {
    return (SRATE*OVERSAMPLE)/(4.0 * Fc);
}

static int calclL(double d1, int thNum, double LS) {
    double dH = 1.0 * OVERSAMPLE;
    double g = calcg(thNum);
    double LSh = (1.0/tuning[thNum]) * LS;
    double LBh = dH * ((d1*d1)/(dH*dH)) - 0.45*d1;
    double z = 0.5 * g * sqrt(1 + 4*(LBh/(g*LSh))) - 0.5*g;
    return (int) (LSh - (z*LSh));
}

// In samples iff LSh is in samples.
static double calcLBh(int thNum, double LSh, int lLint) {
    double g = calcg(thNum);
    double lL = (double) lLint;

    double gLSh = g * LSh;
    double nmrtr1 = (LSh + 0.5*gLSh - lL) * (LSh + 0.5*gLSh - lL);
    nmrtr1 /= gLSh;
    return nmrtr1 - (gLSh / 4.0);
}

static double calcdH(int thNum, double d1, double LS, int lL) {
    double LSh = (1.0/tuning[thNum]) * LS;
    double LBh = calcLBh(thNum, LSh, lL);
    return (d1*d1) / (LBh + 0.45 * d1);
}

static double checkTuning(double d1, double dH, double LSh, double lL, double g) {
    double LBh = dH * ((d1*d1)/(dH*dH)) - 0.45*d1;
    double z = 0.5 * g * sqrt(1 + 4*(LBh/(g*LSh))) - 0.5*g;
    return (SRATE*OVERSAMPLE)/(4 * (lL + (z*LSh)));
}
/*
// Calculates effective length of Birl in centimeters, given fundamental frequency of tube in Hertz.
static inline double calcLS (double Fc) {
    return C_cm / (4.0 * Fc);
}

// Calculates the bore diameter, in centimeters.
static inline double calcd1 (int LC, double LS) {
     return (LS - (double) LC) / 0.3;
}

//is this right???? -JS
// Gives an int length of Birl near LS in centimeters.
static inline int calcLC(double LS) {
    int LC = (int) LS;
    double d1 = calcd1(LC, LS);
    while (d1 < MIN_D1) {
        LC -= 1;
        d1 = calcd1(LC, LS);
    }
    return LC;
}



static inline double calcg (int index) {
    if (index < 0 || index > NUM_NOTES) {
        juce::AlertWindow::showMessageBoxAsync(juce::AlertWindow::WarningIcon,
            "Range Error",
            "Error: Index is out of range",
            "OK");
        return 0.0;
    }
    if (index == 1)
        return tuning[index];
    else
        return tuning[index] / tuning[index+1];
}

static int calclL(double d1, int thNum, double LS) {
    double dH = 1.0 * OVERSAMPLE;
    double g = calcg(thNum);
    double LSh = (1.0/tuning[thNum]) * LS;
    double LBh = dH * ((d1*d1)/(dH*dH)) - 0.45*d1;
    double z = 0.5 * g * sqrt(1 + 4*(LBh/(g*LSh))) - 0.5*g;
    return (int) (LSh - (z*LSh));
}

// TONE HOLES.

// Calculates the effective length of the tube with a cut length of lH.

static inline double calcLSh (int index, double Fc) {
    double frequency = tuning[index] * Fc;
    double LSh = C_cm / (4.0 * frequency);
    return LSh;
}

// Calculates the tonehole diameter.
static inline double calcdH(int index, double d1, int lH, double Fc) {
    // Calculate LSh.
    double LSh = calcLSh(index, Fc);
    // Calculate z.
    double z = (double) lH / LSh;
    // Calculate g.
    double g = calcg(index);
    // Calculate LBh.
    double LBh;
    double LS = calcLS(Fc);
    double x2 = g * LS / 4.0; // Try LS as well.
    double x1 = pow(((z + 0.5*g) / (0.5*g)), 2) - 1;
    LBh = x1 * x2;
    // Solve for dH.
    double dH;
    double chunk = (LBh + 0.45 * d1) / (d1*d1);
    double radical = sqrt(1.0 + 4.0 * lH * chunk);
    dH = (1 + radical) / (2 * chunk);
    return dH;
}

// The cut length of the tube, in centimeters.
static inline int calclH (int index, double d1, double LSh) {
    int lH = (int) LSh;
    double dH = calcdH(index, d1, lH, LSh);
    while (dH < MIN_DH) {
        printf("calclH: dH = %f for this value of LC so we're shortening LC\n", dH);
        lH -= 1;
        dH = calcdH(index, d1, lH, LSh);
    }
    return lH;
}


*/
#endif


