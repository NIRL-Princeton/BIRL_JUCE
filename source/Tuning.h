#ifndef TUNING_H_
#define TUNING_H_

#include "Birl.h"
#include <cmath>
#include <juce_core/juce_core.h>
#include <juce_gui_basics/juce_gui_basics.h>

//this is the right file
// = {2.519840, 2.244920, 2.000000, 1.887750, 1.681790, 1.498310, 1.334830, 1.259920, 1.122460, 1.000000, 0.9439};
const double tuning[] = {10.0/4.0, 18.0/8.0, 2.0/1.0, 15.0/8.0, 5.0/3.0, 3.0/2.0, 4.0/3.0, 5.0/4.0, 9.0/8.0, 1.0, 15.0/16.0};
// const double tuning[] = {2.0/1.0, 15.0/8.0, 5.0/3.0, 3.0/2.0, 4.0/3.0, 5.0/4.0, 9.0/8.0, 1.0, 15.0/16.0};

// const double tuning[] = {1.8877, 1.6818, 1.4982, 1.3348, 1.2599, 1.224, 1.0, 0.9439, 0.9449};
// const double tuning[] = {1.0594, 1.0594, 1.0594, 1.0594, 1.0594, 1.0594, 1.0594, 1.0594, 1.0594, 1.0594, 1.0594};
// const double tuning[] = {2.0, 740.0/392.0, 698.46/392.0, 695.26/392.0, 622.25/392.0, 587.33/392.0, 554.37/392.0, 523.25/392.0, 466.16/392.0};
// const double tuning[] = {20.0/4.0, 36.0/8.0, 4.0/1.0, 30.0/8.0, 10.0/3.0, 6.0/2.0, 8.0/3.0, 10.0/4.0, 18.0/8.0, 2.0, 30.0/16.0};


const double BORE_DIAMETER = 0.2;
/* We're keeping the tonehole diameter proportional to the bore diameter in line with Forster p. 237. */
// const double TONEHOLE_DIAMETER = (15.3 / 19.0) * BORE_DIAMETER;

/* Julius O. Smith's website: Tonehole Filter Design. These are in cm, converted from mm. */
// const double BORE_DIAMETER = 1.89;
const double TONEHOLE_DIAMETER = 0.4765;
const double MINIMUM_CHIMNEY = 0.34;
const double TONEHOLE_RADIUS_CURVATURE = 0.05;
const double FICTITIOUS_DUCT = (MINIMUM_CHIMNEY + TONEHOLE_DIAMETER) * pow ((BORE_DIAMETER / TONEHOLE_DIAMETER), 2.0) - 0.45 * BORE_DIAMETER;



static inline double convertTocm(double samps) {
    return (samps * (C_cm / (SRATE*OVERSAMPLE)));
}

static double convertToSamples(double cm) {
    return (cm * ((SRATE*OVERSAMPLE) / C_cm));
}


/* -----------------------------------------------------------------------------------/
/* CALCULATIONS FOR RATIO BETWEEN SUCCESSIVE TONE-HOLES.
/* -----------------------------------------------------------------------------------*/

// static double calc_g() {
//     return pow (2.0, 1.0/12.0); //equal temperament
// }

static double calc_g(int tonehole_number)
{
    return tuning[tonehole_number] - 1.0;
    // return pow (2.0, (1.0/12.0)) - 1;

}

/* -----------------------------------------------------------------------------------/
/* CALCULATIONS FOR TUBE LENGTHS.
/* -----------------------------------------------------------------------------------*/
/**  Calculates the length of an imaginary closed-open tube.
 * See Equation (1): acoustic length = c / (4* fundamental_frequency)
 */

static double calc_acoustic_length(double fundamental_frequency) {
    return C_cm / (4.0 * fundamental_frequency); // returns length in centimeters
}

/* In actual wind instruments, however, air pressure waves
 * actually travel a small distance beyond the open end
 * of the tube. Therefore, the fundamental wavelength
 * will actually be slightly longer and the fundamental
 * frequency slightly lower than what would be expected from
 * the calculations using calc_theoretical_length.
 *
 * This implies a discrepancy between the "acoustic length"
 * of the tube and the "cut length" of the tube, which refers
 * to the actual length. We call this difference the
 * "end correction" of the tube and label it delta(l_T).
 *
 * The intuition is that for a closed-closed tube,
 * delta(l_T) = 0, since the pressure waves do reflect
 * exactly at the ends of the tube.
 *
 * A closed-open tube of cut length l_T will produce a
 * fundamental frequency equal to that produced by a
 * closed-closed tube of the corresponding effective length L_S.
 *
 * We refer to this fictitious tube as the substitution tube.
 *
 * L_S = l_T + delta(l)
 *
 * Translation: The length of the fictitious closed-closed
 * "substitution tube" is equal to the actual length of the
 * closed-open tube + an end correction delta(l).
 */

/* Calculates the "cut length" of a closed-open tube that
 * corresponds to the acoustic length of a closed-closed
 * tube with the same frequency.
 *
 * See Equation (2): end correction = 0.3 * bore diameter
 */

static float calc_cut_length(double fundamental_frequency, double bore_diameter)
{
    return calc_acoustic_length(fundamental_frequency) - 0.3 * bore_diameter;
}

/* "This phenomenon applies to the length of the tube at each tonehole as well.
 * We refer to the length from the mouthpiece end of the tube to the center
 * of any given tonehole [i] as l_H, which is equal to the effective length of
 * the tube when all toneholes before tonehole [i] are closed.
 * This action is equivalent to shortening the cut length of the tube
 * to a value denoted by l_H, which has a corresponding effective
 * length L_S(h)."
 *
 * See Equation (3): delta(l_H) = z * L_S(h)
 */

/* Calculates the effective length at tonehole i, L_S(h).
 * First, the desired pitch output must be chosen from the tuning array.
 * L_S(h) can be determined using Equation (1) based upon the desire frequency.
 */
static double acoustic_length_at_tonehole (int tonehole_number, double desired_frequency)
{
	return C_m / (4 * desired_frequency) * 100.0; //cm
}

/* Tests different values of tonehole diameter in order to find the cut length.*/

static double cut_length_at_tonehole (int tonehole_number, double desired_frequency, double g)
{
    /* Now let's cut a little of the tube. We'll call this delta(l_H). */
    double cut_guess = -1.0; // start guessing here

    /* Based on Equation (3), we'll estimate a z value. */
    double z_guess = cut_guess / acoustic_length_at_tonehole (tonehole_number, desired_frequency); //cm

    /* We want to find a cut_length that results in what z should be. */
    double z_should_be = 0.5 * g * sqrt(1 + (4 * FICTITIOUS_DUCT) / (g * acoustic_length_at_tonehole (tonehole_number, desired_frequency))) - 0.5 * g;
    /* Now we need to loop different values of cut_guess until z_guess = z based upon Equation (4). */
    while (z_should_be > z_guess)
    {
        cut_guess += 0.00001;

        z_guess = cut_guess / acoustic_length_at_tonehole (tonehole_number, desired_frequency);
    }
    /* Now with the new value of z, use Equation (3) to find the proper cut length. */
    printf("z_guess = %f, z_should_be = %f\n", z_guess, z_should_be);
    double cut_length_for_frequency = acoustic_length_at_tonehole (tonehole_number, desired_frequency) - cut_guess;
    return cut_length_for_frequency; //cm
}



/* -----------------------------------------------------------------------------------/
/* BELOW IS THE PREVIOUS VERSION OF CODE THAT TESTS DIFFERENT VALUES OF D1.
/* -----------------------------------------------------------------------------------/
 * It tests different values of d1 and recuts the tube based
 * upon successively shorter values of d1 longer than
 * the MIN_D1 in Birl.h.
 */
/* "Then, a reasonable starting value must be chosen for either dH or lH
 * which can be used to solve (3).
 * We will guess a value for the diameter of each tonehole, dH."
 */
static int recut_from_d1(double effective_length) {
    int cut_length = (int) effective_length;
    double d1 = (effective_length - (double) cut_length) / 0.3;
    while (d1 < MIN_D1) {
        printf("calcLC: d1 = %f for this value of LC so we're shortening LC!!!!!!\n", d1);
        cut_length -= 1.0;
        d1 = (effective_length - (double) cut_length) / 0.3;
    }
    return cut_length;
}

// In samples if LSh is in samples.
// static double calcLBh(int tonehole_number, double LSh, int lLint) {
//     double g = calcg(tonehole_number);
//     double lL = (double) lLint;
//
//     double gLSh = g * LSh;
//     double nmrtr1 = (LSh + 0.5*gLSh - lL) * (LSh + 0.5*gLSh - lL);
//     nmrtr1 /= gLSh;
//     return nmrtr1 - (gLSh / 4.0);
// }
//
// static double calcdH(int thNum, double d1, double LS, int lL) {
//     double LSh = (1.0/tuning[thNum]) * LS;
//     double LBh = calcLBh(thNum, LSh, lL);
//     return (d1*d1) / (LBh + 0.45 * d1);
// }
//
// static double checkTuning(double d1, double dH, double LSh, double lL, double g) {
//     double LBh = dH * ((d1*d1)/(dH*dH)) - 0.45*d1;
//     double z = 0.5 * g * sqrt(1 + 4*(LBh/(g*LSh))) - 0.5*g;
//     return (SRATE*OVERSAMPLE)/(4 * (lL + (z*LSh)));
// }
/*


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


