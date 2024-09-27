/*--------------------------------------------------------------------------------*
 * MASA analyzer                                                                  *
 * ----------------------------------                                             *
 * (C) 2024 Nokia Technologies Ltd. See LICENSE.md for license.                  *
 *                                                                                *
 *--------------------------------------------------------------------------------*/

#include "masaCommonHeader.h"

const short unsigned bandBottomLimits[MASA_BAND_COUNT] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9,  10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 25, 30, 40};
const short unsigned bandTopLimits[MASA_BAND_COUNT]    = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 25, 30, 40, 60};

const float eigenMobileEqB[3] = {0.942999149551865f, -1.30698914427719f,  0.503968729438150f};
const float eigenMobileEqA[3] = {1.0f,               -1.279632424997809f, 0.477592250072517f};

const uint8_t defaultFormatDescriptor[8] = {0x49, 0x56, 0x41, 0x53, 0x4D, 0x41, 0x53, 0x41}; // "IVASMASA"
