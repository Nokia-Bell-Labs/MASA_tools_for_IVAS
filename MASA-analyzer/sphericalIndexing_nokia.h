/*--------------------------------------------------------------------------------*
 * MASA analyzer                                                                  *
 * ----------------------------------                                             *
 * (C) 2024 Nokia Technologies Ltd.. See LICENSE.md for license.                  *
 *                                                                                *
 *--------------------------------------------------------------------------------*/

#ifndef MASAANALYZER_SPHERICALINDEXING_NOKIA_H
#define MASAANALYZER_SPHERICALINDEXING_NOKIA_H

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
} SphericalGridData_nokia;

void generateSphericalGrid_nokia(SphericalGridData_nokia *data);

void deindexDirection_nokia(
    uint16_t sphIndex,             /* i  : Spherical index */
    const SphericalGridData_nokia *data, /* i  : Prepared spherical grid */
    float *theta,                  /* o  : Elevation */
    float *phi                     /* o  : Azimuth */
);

void indexDirection_nokia(
    uint16_t *sphIndex,               /* o  : output index for direction */
    float theta,                      /* i  : input elevation to be quantized */
    float phi,                        /* i  : input azimuth to be  quantized */
    const SphericalGridData_nokia *gridData /* i  : generated grid data */
);

void quantize_theta_phi_nokia(
    int16_t no_th,             /* i  : elevation codebook size */
    const int32_t* no_phi_loc, /* i  : number of azimuth values for each elevation codeword */
    float abs_theta,           /* i  : absolute value of elevation to be quantized */
    int16_t *id_phi,           /* o  : azimuth index */
    float phi,                 /* i  : input azimuth value; to be quantized  */
    int16_t *id_theta,         /* o  : elevation index */
    float *phi_q               /* o  : rotated quantized azimuth */
);

int16_t quantize_phi_nokia(    /* o  : index azimuth */
    float phi,           /* i  : azimuth value */
    int16_t flag_delta,  /* i  : flag indicating if the azimuth codebook is translated or not */
    float *phi_hat,      /* o  : quantized azimuth */
    int16_t n            /* i  : azimuth codebook size */
);

float direction_distance_nokia( /* o  : distortion value */
    float theta,          /* i  : elevation absolute value */
    float theta_hat,      /* i  : quantized elevation value in absolute value */
    float phi,            /* i  : azimuth value */
    float phi_hat         /* i  : quantized azimuth value */
);

int16_t quantize_theta_nokia( /* o  : index of quantized value */
    float x,            /* i  : theta value to be quantized */
    int16_t no_cb,      /* i  : number of codewords */
    float *xhat         /* o  : quantized value */
);

#endif //MASAANALYZER_SPHERICALINDEXING_NOKIA_H
