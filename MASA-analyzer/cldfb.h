/*--------------------------------------------------------------------------------*
 * MASA analyzer                                                                  *
 * ----------------------------------                                             *
 * (C) 2024 Nokia Technologies Ltd. See LICENSE.md for license.                  *
 *                                                                                *
 *--------------------------------------------------------------------------------*/

#ifndef MASAANALYZER_CLDFB_H
#define MASAANALYZER_CLDFB_H

#include "pocketfft/pocketfft.h"

#define CLDFB_BAND_COUNT 60
#define CLDFB_SLOT_COUNT 16

enum CldfbDelayMode {
    CLDFB_5_MILLISECONDS = 0,
    CLDFB_1_25_MILLISECONDS = 1
};

enum CldfbDirection {
    CLDFB_FORWARD = 0,
    CLDFB_INVERSE = 1
};

typedef struct {
    cfft_plan fftSetup;
    float scale;
    float *prototypeFilter;
    float windowModulatorReal[120];
    float windowModulatorImag[120];
    float FDmodulatorReal[60];
    float FDmodulatorImag[60];
    float tdBuffer[10][60];
    float realBuf120samples[120];
    float complexBuf60SamplesReal[60];
    float complexBuf60SamplesImag[60];
    int reorder[60];
    int bufIdx;
} cldfb_struct;

void cldfbInit(cldfb_struct *h, enum CldfbDelayMode mode, enum CldfbDirection direction);

void cldfbDestroy(cldfb_struct *h);

void cldfbForward(cldfb_struct *h, float *inTD, float *outFDReal, float *outFDImag);

void cldfbInverse(cldfb_struct *h, float *inFDReal, float *inFDImag, float *outTD);

#endif //MASAANALYZER_CLDFB_H
