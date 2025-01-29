#include "PhysicalModel.h"

PhysicalModel::PhysicalModel()
{

    numToneHoles_ = NUM_OF_KEYS;
    numTubes_ = NUM_OF_KEYS+1;
    tubeIndex_ = 0;
    ::setTuning(EQUAL_TEMPERED);

//    tune(239.9);
    tubes_[0] = initTube(25);
    
    
    dcBlocker_ = initDCFilter(0.995);
    dcBlocker2_ = initDCFilter(0.995);
    biquad_ = initBiquad();
    biquadSetCoeffs(biquad_, 0.169301, 0.338601, 0.169301, -0.482013, 0.186622);
    pf_ = initSVF(2000.0, 0.5);
    lp_ = initSVF(5000.0, 0.5);
    pf2_ = initSVF(1000.0, 1.0);
    lp2_ = initSVF(5000.0, 0.5);
    noiseBP_ = initSVF(16000.0, 1.0);
        
    outputGain_ = 1.0;
    noiseGain_ = 0.2;
    
    breathPressure_ = 0.0;
    prevBreathPressure_ = 0.0;
    count = 0;
    min = 0.0;
    max = 0.0;
    m_drive_ = 0.0;
    shaperMix_ = 0.0;
    srand((unsigned)(NULL));

}
PhysicalModel::~PhysicalModel(){}


void PhysicalModel::setToneHoleIndex(int index) {
    if (index < 0 || index >= numToneHoles_) {
        printf("index out of range: %d\n", index);
        return;
    }
    toneHoleIndex_ = index;
}
void PhysicalModel::setToneHole(double newValue) {
    double new_coeff;
    int index = toneHoleIndex_;
    double thCoeff = thCoeff_[index];
    
    if ( newValue <= 0.0 )
        new_coeff = 0.9995;
    else if ( newValue >= 1.0 )
        new_coeff = thCoeff;
    else
        new_coeff = ( newValue * (thCoeff - 0.9995) ) + 0.9995;
    setA1PoleZero( toneHoles_[index], new_coeff );
    setB0PoleZero( toneHoles_[index], new_coeff );
}

void PhysicalModel::setToneHoleRadius(double newRadius) {
     if (newRadius < MIN_TONEHOLE_RADIUS || newRadius > MAX_TONEHOLE_RADIUS) {
         printf("Radius is too big or too small: %f\n", newRadius);
         return;
     }
     int index = toneHoleIndex_;

     rth_[index] = newRadius;
     scatter_[index] = -pow(newRadius,2) / ( pow(newRadius,2) + 2*pow(rb_,2) );

     // Calculate toneHole coefficients.
     double te = newRadius;    // effective length of the open hole
     thCoeff_[index] = (te*2*(SRATE*OVERSAMPLE) - C_m) / (te*2*(SRATE*OVERSAMPLE) + C_m);
 }
void PhysicalModel::setLengthIndex(int index) {
    if (index < 0 || index >= numTubes_) {
        printf("index out of range: %d\n", index);
        return;
    }
    tubeIndex_ = index;
}
void PhysicalModel::setLength(int newlen) {
    if (newlen < 0 || newlen >= MAX_TUBE_LENGTH) {
        printf("length out of range: %d\n", newlen);
        return;
    }
    int index = tubeIndex_;

    tubeLengths_[index] = newlen;
    tubes_[index] = initTube(newlen);
}
void PhysicalModel::setBreathPressure(double input) {
    breathPressure_ = input;
}

void PhysicalModel::setLPCutoff(double cut) {
    cut = clip(cut, 30.0, 16000.0);
    printf("lp cutoff: %f\n", cut);
    setCutoffSVF(lp_, cut);
    setCutoffSVF(lp2_, cut);
}
void PhysicalModel::setLPQ(double Q) {
    Q = clip(Q, 0.0, 1.0);
    printf("lp Q: %f\n", Q);
    setQSVF(lp_, Q);
    setQSVF(lp2_, Q);
}

void PhysicalModel::setPFCutoff(double cut) {
    cut = clip(cut, 30.0, 16000.0);
    printf("pf cutoff: %f\n", cut);
    setCutoffSVF(pf_, cut);
    setCutoffSVF(pf2_, cut);
}

void PhysicalModel::setPFQ(double Q) {
    Q = clip(Q, 0.0, 1.0);
    printf("pf Q: %f\n", Q);
    setQSVF(pf_, Q);
    setQSVF(pf2_, Q);
}

void PhysicalModel::setNoiseBPCutoff(double cut) {
    cut = clip(cut, 30.0, 16000.0);
    printf("noiseBP cutoff: %f\n", cut);
    setCutoffSVF(noiseBP_, cut);
}

void PhysicalModel::setNoiseGain(double gain) {
    gain = clip(gain, 0.0, 1.0);
    printf("noise gain: %f\n", gain);
    noiseGain_ = gain;
}

void PhysicalModel::setNoiseBPQ(double Q) {
    Q = clip(Q, 0.0, 1.0);
    printf("noiseBP Q: %f\n", Q);
    setQSVF(noiseBP_, Q);
}

void PhysicalModel::setShaperDrive(double d) {
    d = clip(d, 0.0, 1.0);
    printf("m_drive_: %f\n", d);
    m_drive_ = d;
}

void PhysicalModel::setShaperMix(double d) {
    d = clip(d, 0.0, 1.0);
    printf("shaper mix: %f\n", shaperMix_);
    shaperMix_ = d;
}

void PhysicalModel::setTuningWrapper(int t) {
    ::setTuning((tuningSystem) t);
}

void PhysicalModel::setNumToneHoles(int toneholes){
    numToneHoles_ = toneholes;
    toneHoleIndex_ = toneholes + 1;
}

void PhysicalModel::setCustomTuning(double freqs[NUM_NOTES]) {
//    double freqs[] = {f9, f8, f7, f6, f5, f4, f3, f2, f1, f0, fn1};
    populateCustomTuning(freqs);
    tune(freqs[NUM_NOTES - 1]);
}
void PhysicalModel::calcTHCoeffs() {
    // Calculate initial tone hole three-port scattering coefficients
    for (int i = 0; i < MAX_TONEHOLES; i++) {
        scatter_[i] = -pow(rth_[i],2) / ( pow(rth_[i],2) + 2*pow(rb_,2) );

        // Calculate toneHole coefficients and set for initially open.
        thCoeff_[i] = (rth_[i]*2*(SRATE*OVERSAMPLE) - C_m) / (rth_[i]*2*(SRATE*OVERSAMPLE) + C_m);

        // Initialize tone holes.
        toneHoles_[i] = initPoleZero();
        setA1PoleZero(toneHoles_[i], -thCoeff_[i]);
        setB0PoleZero(toneHoles_[i], thCoeff_[i]);
        setB1PoleZero(toneHoles_[i], -1.0);
    }
}
void PhysicalModel::tune(double Fc) {
//    double LS = calcLS(Fc);
//    int LC = calcLC(LS);
//    double d1 = calcd1(LC, LS);
//    printf("LC: %d\n", LC);
//    printf("d1: %f\n", d1);
//
//    int prevlL = 0;
//
//    for (int i = 0; i < numTubes_; i++) {
////        printf("tubelength %d: %d\n", i, calclL(d1, i, LS) - prevlL);
////        tubeLengths_[i] = calclL(d1, i, LS) - prevlL;
//        tubeLengths_[i] = calclL(d1, i, LS);
//
//        if (tubeLengths_[i] == 0) {
//            printf("ERROR: Integer delay line lengths clash!!!!! Use a different tuning or try oversampling.\n");
//            return;
//        }
//        // if (tubes_[i] != NULL) {
//        //     printf("WANNA BE FREEEEEEE\n");
//        //     freeTube(tubes_[i]);
//        // }
//
//        if (i == 0) {
//            tubes_[i] = initTube(tubeLengths_[i]);
//            // prevlL += tubeLengths_[i] + correction;
//            prevlL += tubeLengths_[i];
//        } else {
//            tubes_[i] = initTube(tubeLengths_[i]);
//            prevlL += tubeLengths_[i];
//        }
//        printf("th %d: lL = %d\n", i, tubeLengths_[i]);
//    }
//    originalrb_ = d1 / 2.0;    // main bore radius, in cm.
//    rb_ = originalrb_;
//    int lL = tubeLengths_[0];
//    for (int i = 0; i < numToneHoles_; i++) {
//        originalRth_[i] = calcdH(i, d1, LS, lL)/2.0;
//        rth_[i] = originalRth_[i];
//        lL += tubeLengths_[i+1];
//    }
//    calcTHCoeffs();
}

double PhysicalModel::interpolateLinear(double a, double b, double alpha) {
    return (alpha * a) + ((1.0-alpha) * b);
}
 
#if 0

double PhysicalModel::birlTick()
{
    double breathInterp[OVERSAMPLE];
    for (int i = 0; i < OVERSAMPLE; i++) {
        breathInterp[i] = interpolateLinear(breathPressure_, prevBreathPressure_, (double) (i+1) / (double) OVERSAMPLE);
    }
    prevBreathPressure_ = breathPressure_;

    double pap;
    double pbm;
    double pthm;
    double outsamp = 0.0;
    double scatter;
    double bellReflected;

    for (int t = 0; t < OVERSAMPLE; t++) {
        double breath = breathInterp[t];
        double noise = ((double)rand()/(double)RAND_MAX);
        noise = noiseGain_ * (inputSVFBand(noiseBP_, noise));
        // double noise = noiseGain_ * rand();
        breath += breath * noise;

        // Calculate the differential pressure = reflected - mouthpiece pressures
        double pressureDiff = accessDelayLine(tubes_[0]->lower) - breath;
        double reedLookup = pressureDiff * reedTable( pressureDiff );
//        breath = tanhClip(breath + reedLookup);
        breath = LEAF_clip(-1.0, breath + reedLookup, 1.0);


        // debug("%d\n", tubeLengths_[0]);
        // Helps reduce high-pitched noise.
        // breath = inputBiquad(biquad_, breath);
//        breath = interpolateLinear(shaper(breath, m_drive_), breath, shaperMix_);
        // breath = inputSVFPeak(pf_, breath);
        // breath = inputSVFLP(lp_, breath);
        // breath = inputDCFilter(dcBlocker_, breath);

        if (breath > 1.0 || breath < -1.0) {
            printf("breath going out of bounds of -1 to 1: %f\n", breath);
        }

        for (int i = 0; i < numToneHoles_; i++) {
            // Index in tubes_[] of tube positioned before toneHoles[i].
            int a = i + FRONT_TUBES - 1;
            // Index in tubes_[] of tube positioned after toneHoles[i].
            int b = i + FRONT_TUBES;

            // Three-port junction scattering.
            pap = accessDelayLine(tubes_[a]->upper);
            pbm = accessDelayLine(tubes_[b]->lower);
            pthm = z1PoleZero(toneHoles_[i]);

            scatter = scatter_[i] * (pap + pbm - (2 * pthm));
            pbp_[i] = pap + scatter;
            pam_[i] = pbm + scatter;
            pthp_[i] = pap + scatter + pbm - pthm;

            // Sample output at tubes_[SAMPLE_INDEX].
            if (a == SAMPLE_INDEX) {
                outsamp += pap + pam_[i];
            }

            // Bell reflection at last tube.
            if (i == numToneHoles_ - 1) {
                double bell = accessDelayLine(tubes_[b]->upper);
                // double bell2 = shaper(bell1, m_drive_);
                // double bell4 = inputSVFPeak(pf2_, bell1);
                // bell = inputSVFLP(lp2_, bell);
                // // Reflection = Inversion + gain reduction + lowpass filtering.
                // bell = inputSVFLP(lp2_, bell);
                // bell = inputDCFilter(dcBlocker2_, bell);
                bellReflected = tSVF_tick(&inLoopLP, (bell * -0.995));
                // bellReflected = filter_.tick(bell * -0.995);
            }
        }

        // Perform all inputs at the end so that we're not altering
        // state prior to calculations.
        for (int i = 0; i < numToneHoles_; i++) {
            // Index in tubes_[] of tube positioned before toneHoles[i].
            int a = i + FRONT_TUBES - 1;
            // Index in tubes_[] of tube positioned after toneHoles[i].
            int b = i + FRONT_TUBES;

            inputPoleZero(toneHoles_[i], pthp_[i]);
            // toneHoles_[i]->tick(pthp_[i]);
            inputDelayLine(tubes_[a]->lower, pam_[i]);
            inputDelayLine(tubes_[b]->upper, pbp_[i]);
        }

        inputDelayLine(tubes_[0]->upper, breath);
        inputDelayLine(tubes_[numTubes_-1]->lower, bellReflected);

    }

    count++;

    // debug("%.5f -0-> ", breath);
    // for (int i = 0; i < numToneHoles_; i++) {
    //     debug("\t%.5f\t-%d->", pbp_[i], i+1);
    // }
    // debug("\n");
    // for (int i = 0; i < numToneHoles_; i++) {
    //     debug("%.5f\t<-%d-\t", pam_[i], i);
    // }
    // debug("%.5f\n", bellReflected);
    // debug("range: %f - %f\n", min, max);
    // debug("\n");


    // Clipping, account for oversampling, and gain.
    outsamp /= (double) OVERSAMPLE;
    outsamp = tanhClip(outsamp);
    outsamp *= outputGain_;
    return outsamp;
}

#endif

double PhysicalModel::birlTick() {
    double outsamp = 0.0;
    double bellReflected;
    
    double breath = breathPressure_;
    double noise = (double) rand() / (double) RAND_MAX;
    noise = noiseGain_ * (inputSVFBand(noiseBP_, noise));
    breath += breath * noise;
    
    // REED
    double pressureDiff = accessDelayLine(tubes_[0]->lower) - breath;
    double reedLookup = pressureDiff * reedTable (pressureDiff);
    breath = LEAF_clip(-1.0, breath + reedLookup, 1.0);
    
    double bell = accessDelayLine(tubes_[0]->upper);
    outsamp = bell;
    bellReflected = bell * -0.995;
    
    inputDelayLine(tubes_[0]->upper, breath);
    inputDelayLine(tubes_[0]->lower, bellReflected);
    tanhClip(outsamp);
    
    return outsamp;
}


