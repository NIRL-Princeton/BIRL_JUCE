#ifndef TUBE_H_
#define TUBE_H_

#include "sfx.h"
//extern LEAF leaf;
//
//LEAF leaf;
typedef struct _DelayLine {
    double *data;
    int length;
    double *pointer;
    double *end;
} DelayLine;

static inline DelayLine *initDelayLine(int len) {
    if (len < 0) {
        printf("Tried to initialize delay line of length 0\n");
        return NULL;
    }
    DelayLine *dl = (DelayLine *) calloc(len, sizeof(DelayLine));
    dl->length = len;
    dl->data = (double *) calloc(len, len * sizeof(double));
    dl->pointer = dl->data;
    dl->end = dl->data + len - 1;
    return dl;
}

static inline void freeDelayLine(DelayLine *dl) {
    if (dl && dl->data)
        free(dl->data);
    dl->data = 0;
    free(dl);
}

// Don't call this more than once in tick()!
static inline void inputDelayLine(DelayLine *dl, double insamp) {
    double *ptr = dl->pointer;
    *ptr = insamp;
//    ptr++;
//    if (ptr > dl->end)
//        ptr = dl->data;
    dl->pointer = ptr;
}

static inline double accessDelayLine(DelayLine *dl) {
    return *(dl->pointer);
}


typedef struct _Tube {
//    DelayLine *upper, *lower;
    tLinearDelay upper, lower;
} Tube;


static inline Tube *initTube(int len, LEAF &leaf) {
    Tube *tube = (Tube *) calloc(1, sizeof(Tube));
//    tube->upper = initDelayLine(len);
//    tube->lower = initDelayLine(len);
    tLinearDelay_init(&tube->upper, len, len+1, &leaf);
    tLinearDelay_init(&tube->lower, len, len+1, &leaf);
    return tube;
}

static inline void freeTube(Tube *tube) {
    // freeDelayLine(tube->upper);
    // freeDelayLine(tube->lower);
    tLinearDelay_free(&tube->upper);
    tLinearDelay_free(&tube->lower);
    free(tube);
}

typedef struct _FracTube {
    tLinearDelay upper, lower;
} FracTube;

static inline FracTube *initFracTube(float len, LEAF &leaf) {
    FracTube *fractube = (FracTube *) calloc (1, sizeof(FracTube));
    tLinearDelay_init(&(fractube->upper), len, int(len+1), &leaf);
    tLinearDelay_init(&(fractube->lower), len, int(len+1), &leaf);
    return fractube;
}

static inline void freeFracTube(FracTube *fractube) {
    tLinearDelay_free(&(fractube->upper));
    tLinearDelay_free(&(fractube->lower));
    free(fractube);
}

#endif
