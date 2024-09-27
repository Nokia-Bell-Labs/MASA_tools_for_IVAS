/*--------------------------------------------------------------------------------*
 * MASA analyzer                                                                  *
 * ----------------------------------                                             *
 * (C) 2024 Nokia Technologies Ltd. See LICENSE.md for license.                  *
 *                                                                                *
 *--------------------------------------------------------------------------------*/

#include <math.h>
#include "masaCommonHeader.h"
#include "cldfb.h"
#include "cldfb_prototype_filters.h"
#include "pocketfft/pocketfft.h"


void cldfbInit(cldfb_struct *h, enum CldfbDelayMode mode, enum CldfbDirection direction) {
    if (direction == CLDFB_FORWARD) {
        // Forward
        h->fftSetup = make_cfft_plan(60);
        // The scale and phase corrections 
        h->scale = 1.0f / 0.75f;
        for (int s = 0; s < 60; s++) {
            float phase = (float) MASA_PI * (1.0f + (1.0f - 1.0f / 60.0f) * ((float) s + 0.25f));
            if (mode == CLDFB_1_25_MILLISECONDS) {
                phase -= 0.5f * (float) MASA_PI;
            }
            h->FDmodulatorReal[s] = cosf(phase);
            h->FDmodulatorImag[s] = sinf(phase);
        }
    } else {
        // Inverse
        h->fftSetup = make_cfft_plan(60);
        // The scale and phase corrections 
        h->scale = 0.75f;
        for (int s = 0; s < 60; s++) {
            float phase = (float) MASA_PI * (1.0f + (1.0f + 1.0f / 60.0f) * ((float) s + 0.25f));
            h->FDmodulatorReal[s] = cosf(phase);
            h->FDmodulatorImag[s] = sinf(phase);
        }
    }

    // Read proto filter
    if (mode == CLDFB_5_MILLISECONDS) {
        h->prototypeFilter = CLDFB_PROTOTYPE_FILTER_5_00_MILLISECONDS;
    } else {
        h->prototypeFilter = CLDFB_PROTOTYPE_FILTER_1_25_MILLISECONDS;
    }

    // Modulator for 0.5 bin frequency shift
    for (int s = 0; s < 120; s++) {
        float phase = -(float) MASA_PI * ((float) s / 120.0f + 0.5f);
        h->windowModulatorReal[s] = cosf(phase);
        h->windowModulatorImag[s] = sinf(phase);
    }

    // Formulate reorder indexing. Relates to how FFT bins are wrapped
    for (int b = 0; b < 60; b++) {
        int idx;
        if (b % 2 == 0) {
            idx = b / 2;
        } else {
            idx = 59 - (b - 1) / 2;
        }
        h->reorder[b] = idx;
    }

    // Flush time domain buffer that is has 10 frames * 60 samples history
    for (int frame = 0; frame < 10; frame++) {
        for (int s = 0; s < 60; s++) {
            h->tdBuffer[frame][s] = 0.0f;
        }
    }
    h->bufIdx = 0;
}

// Note: complex-to-complex fft
void fft(cfft_plan fftSetup, float *inCpxReal, float *inCpxImag, float *outCpxReal, float *outCpxImag) {
    double temp[120];
    for (int b = 0; b < 60; b++) {
        temp[b*2] = inCpxReal[b];
        temp[b*2+1] = inCpxImag[b];
    }
    cfft_forward(fftSetup, temp, 1);
    for (int b = 0; b < 60; b++) {
        outCpxReal[b] = (float) temp[b*2];
        outCpxImag[b] = (float) temp[b*2+1];
    }
}

void cldfbDestroy(cldfb_struct *h)
{
    if (h != NULL)
    {
        destroy_cfft_plan(h->fftSetup);
    }
}


// Note: complex-to-complex ifft
void ifft(cfft_plan fftSetup, float *inCpxReal, float *inCpxImag, float *outCpxReal, float *outCpxImag) {
    double temp[120];
    for (int b = 0; b < 60; b++) {
        temp[b*2] = inCpxReal[b];
        temp[b*2+1] = inCpxImag[b];
    }
    cfft_backward(fftSetup, temp, 1/60.0f);
    for (int b = 0; b < 60; b++) {
        outCpxReal[b] = (float) temp[b*2];
        outCpxImag[b] = (float) temp[b*2+1];
    }
}

void cldfbForward(cldfb_struct *h, float *inTD, float *outFDReal, float *outFDImag) {
    // Flush buffers
    for (int s = 0; s < 120; s++) {
        h->realBuf120samples[s] = 0.0f;
    }
    for (int s = 0; s < 60; s++) {
        h->complexBuf60SamplesReal[s] = 0.0f;
        h->complexBuf60SamplesImag[s] = 0.0f;
    }

    // Collect input
    for (int s = 0; s < 60; s++) {
        h->tdBuffer[h->bufIdx][s] = inTD[s];
    }

    // Increment circular buffer index
    h->bufIdx++;
    h->bufIdx %= 10;
    int idx = h->bufIdx;

    // Apply prototype filter that is 600 length, and wrap to 120-length buffer
    int s600 = 0;
    for (int frame = 0; frame < 10; frame++) {
        for (int s = 0; s < 60; s++) {
            h->realBuf120samples[s600 % 120] += h->prototypeFilter[s600] * h->tdBuffer[idx][s];
            s600++;
        }
        idx++;
        idx %= 10;
    }

    // Apply half-band complex modulation and wrap to 60-length buffer
    for (int s = 0; s < 120; s++) {
        h->complexBuf60SamplesReal[s % 60] += h->windowModulatorReal[s] * h->realBuf120samples[s];
        h->complexBuf60SamplesImag[s % 60] += h->windowModulatorImag[s] * h->realBuf120samples[s];
    }

    // Convert to frequency bins
    float tmpCpxReal[60];
    float tmpCpxImag[60];
    float temp;
    fft(h->fftSetup, h->complexBuf60SamplesReal, h->complexBuf60SamplesImag, tmpCpxReal, tmpCpxImag);

    // Apply phase corrections
    for (int b = 0; b < 60; b++) {
        temp = tmpCpxReal[b] * h->FDmodulatorReal[b] - tmpCpxImag[b] * h->FDmodulatorImag[b];
        tmpCpxImag[b] = tmpCpxReal[b] * h->FDmodulatorImag[b] + tmpCpxImag[b] * h->FDmodulatorReal[b];
        tmpCpxReal[b] = temp;
        if (b >= 30) {
            tmpCpxReal[b] = -tmpCpxReal[b];
        }
    }

    // Bin reordering, relates to how FFT is used, and scaling
    for (int b = 0; b < 60; b++) {
        outFDReal[b] = h->scale * tmpCpxReal[h->reorder[b]];
        outFDImag[b] = h->scale * tmpCpxImag[h->reorder[b]];
    }
}


void cldfbInverse(cldfb_struct *h, float *inFDReal, float *inFDImag, float *outTD) {
    // Bin reordering, relates to how FFT is used
    float tmpCpxReal[60];
    float tmpCpxImag[60];
    for (int b = 0; b < 60; b++) {
        tmpCpxReal[h->reorder[b]] = inFDReal[b];
        tmpCpxImag[h->reorder[b]] = inFDImag[b];
    }

    // Apply phase corrections
    for (int b = 0; b < 60; b++) {
        if (b >= 30) {
            tmpCpxReal[b] = -tmpCpxReal[b];
        }
        float temp = tmpCpxReal[b] * h->FDmodulatorReal[b] - tmpCpxImag[b] * h->FDmodulatorImag[b];
        tmpCpxImag[b] = tmpCpxReal[b] * h->FDmodulatorImag[b] + tmpCpxImag[b] * h->FDmodulatorReal[b];
        tmpCpxReal[b] = temp;
    }
    // Convert to frequency bins
    ifft(h->fftSetup, tmpCpxReal, tmpCpxImag, h->complexBuf60SamplesReal, h->complexBuf60SamplesImag);

    // Flush the last frame, already output in the overlap-add processing
    int idx = (h->bufIdx - 1 + 10) % 10;
    for (int s = 0; s < 60; s++) {
        h->tdBuffer[idx][s] = 0.0f;
    }
    idx = h->bufIdx;

    // Apply prototype filter and inverse complex modulation
    int s600 = 0;
    for (int frame = 0; frame < 10; frame++) {
        for (int s = 0; s < 60; s++) {
            float temp = h->complexBuf60SamplesReal[s] * h->windowModulatorReal[s600 % 120];
            temp += h->complexBuf60SamplesImag[s] * h->windowModulatorImag[s600 % 120];
            temp *= h->prototypeFilter[599 - s600];
            h->tdBuffer[idx][s] = temp;
            s600++;
        }
        idx++;
        idx %= 10;
    }

    // Scale output
    for (int s = 0; s < 60; s++) {
        outTD[s] = h->scale * h->tdBuffer[h->bufIdx][s];
    }
    h->bufIdx++;
    h->bufIdx %= 10;
}
