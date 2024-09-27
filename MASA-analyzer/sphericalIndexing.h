/*--------------------------------------------------------------------------------*
 * MASA analyzer                                                                  *
 * ----------------------------------                                             *
 * (C) 2024 Nokia Technologies Ltd. See LICENSE.md for license.                  *
 *                                                                                *
 *--------------------------------------------------------------------------------*/

#ifndef MASAANALYZER_SPHERICALINDEXING_H
#define MASAANALYZER_SPHERICALINDEXING_H

#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include "masaCommonHeader.h"

#define SPH_CB_SIZE_16BIT 122

typedef struct
{
    float theta_cb[SPH_CB_SIZE_16BIT];
    int16_t no_theta;
    int32_t no_phi[SPH_CB_SIZE_16BIT];  /* here no_phi contains the cummulated values of points on circles */
} SphericalGridData;

void generateSphericalGrid(SphericalGridData *data);

void indexDirection(
    uint16_t *sphIndex,               /* o  : output index for direction */
    float theta,                      /* i  : input elevation to be quantized */
    float phi,                        /* i  : input azimuth to be  quantized */
    const SphericalGridData *gridData /* i  : generated grid data */
);

#endif //MASAANALYZER_SPHERICALINDEXING_H
