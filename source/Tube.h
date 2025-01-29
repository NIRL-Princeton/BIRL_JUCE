#ifndef BIRL_TUBE_H
#define BIRL_TUBE_H

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
    DelayLine *dl = (DelayLine *)calloc(len, sizeof(DelayLine));
    dl->length = len;
    dl->data = (double *)calloc(len, len * sizeof(double));
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
    ptr++;
    if (ptr > dl->end)
        ptr = dl->data;
    dl->pointer = ptr;
}

static inline double accessDelayLine(DelayLine *dl) {
    return *(dl->pointer);
}

typedef struct _Tube {
    DelayLine *upper, *lower;
} Tube;

static inline Tube *initTube(int len) {
    Tube *tube = (Tube *) calloc(1, sizeof(Tube));
    tube->upper = initDelayLine(len);
    tube->lower = initDelayLine(len);
    return tube;
}

static inline void freeTube(Tube *tube) {
    freeDelayLine(tube->upper);
    freeDelayLine(tube->lower);
    free(tube);
}

#endif
