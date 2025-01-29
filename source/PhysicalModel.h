#ifndef PHYSICAL_MODEL_H
#define PHYSICAL_MODEL_H

#include "Birl.h"
#include "Filters.h"
#include "Tube.h"
#include "Tuning.h"
#include "leaf-filters.h"


class PhysicalModel
{
public:
    PhysicalModel();
    ~PhysicalModel();
    
    void setNumToneHoles(int);
    void setToneHoleIndex(int);
    void setToneHole(double);
    void setToneHoleRadius(double);
    void setLengthIndex(int);
    void setLength(int);
    void setBreathPressure(double);
    void setLPCutoff(double);
    void setLPQ(double);
    void setPFCutoff(double);
    void setPFQ(double);
    void setNoiseBPCutoff(double);
    void setNoiseGain(double);
    void setNoiseBPQ(double);
    void setShaperDrive(double);
    void setShaperMix(double);
    void setTuningWrapper(int);
    void setCustomTuning(double freqs[NUM_NOTES]);
    void calcTHCoeffs();
    void tune(double);
    double interpolateLinear(double a, double b, double alpha);
    tSVF inLoopLP; // in loop low pass
    double birlTick();


    
protected:
    int numTubes_;
    int tubeIndex_;
    int numToneHoles_;
    int toneHoleIndex_;
    
    double pam_[MAX_TONEHOLES];
    double pbp_[MAX_TONEHOLES];
    double pthp_[MAX_TONEHOLES];
    double scatter_[MAX_TONEHOLES];
    double thCoeff_[MAX_TONEHOLES];
    
    Tube *tubes_ [MAX_TUBES];
    PoleZero *toneHoles_[MAX_TONEHOLES];
    
    
    DCFilter *dcBlocker_;
    DCFilter *dcBlocker2_;
    Biquad *biquad_;
    Biquad *biquad2_;
    SVF *pf_;
    SVF *lp_;
    SVF *pf2_;
    SVF *lp2_;
    SVF *noiseBP_;
        
    double rb_;
    double originalrb_;
    double rhGain_;
    double outputGain_;
    double noiseGain_;
    double m_drive_;
    double shaperMix_;
    
    double breathPressure_{0};
    double prevBreathPressure_{0};
    
    int count;
    double min, max;
};

#endif
