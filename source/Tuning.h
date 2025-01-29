#ifndef TUNING_H_
#define TUNING_H_

#include "Birl.h"
#include <juce_core/juce_core.h>
#include <juce_gui_basics/juce_gui_basics.h>
// ONEHOLE. TWO NOTES.
const int NUM_NOTES = 2;

// CONSTS.
const double MIN_D1 = 0.5;
const double MIN_DH = 0.03;
const double DH_FIRST_GUESS = 0.03;

static double tuning[] = {1.0, 1.0};

// Calculates effective length of Birl in centimeters, given fundamental frequency of tube in Hertz.
static inline double calcLS (double Fc) {
    return C_cm / (4.0 * Fc);
}

// Calculates the bore diameter, in centimeters.
static inline double calcd1 (int LC, double LS) {
     return (LS - (double) LC) / 0.3;
}

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

#endif


