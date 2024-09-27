/*--------------------------------------------------------------------------------*
 * MASA analyzer                                                                  *
 * ----------------------------------                                             *
 * (C) 2024 Nokia Technologies Ltd. See LICENSE.md for license.                  *
 *                                                                                *
 *--------------------------------------------------------------------------------*/

#ifndef MASAANALYZER_MASACOMMONHEADER_H
#define MASAANALYZER_MASACOMMONHEADER_H

#include <stdint.h>
#include "cldfb.h"

#define FRAME_SIZE 960u
#define SUBFRAME_SIZE 240u   // To obtain 5 ms time resolution
#define SUBFRAME_COUNT ((FRAME_SIZE) / (SUBFRAME_SIZE))
#define CLDFB_SLOTS_PER_SUBFRAME ((CLDFB_SLOT_COUNT) / (SUBFRAME_COUNT))
#define MASA_BAND_COUNT 24u
#define MASA_MAXIMUM_DIRECTIONS 2u
#define MAX_TRANSPORT_COUNT 2u

#define EIGENMIKE_CHANNEL_COUNT 32u
#define FOA_CHANNEL_COUNT 4u
#define HOA2_CHANNEL_COUNT 9u

// POSIX pi value
#define MASA_PI 3.14159265358979323846264338327950288f
#define MASA_EPSILON 1e-12f


typedef enum
{
    MASA_TRANSPORT_MONO,
    MASA_TRANSPORT_STEREO
} TransportType;


typedef struct
{
    uint8_t formatDescriptor[8];   // 8x 8 bits
    uint8_t numberOfDirections;  // 1 bit
    uint8_t numberOfChannels;    // 1 bit
    uint8_t sourceFormat;        // 2 bits
    uint8_t transportDefinition; // 3 bits
    uint8_t channelAngle;        // 3 bits
    uint8_t channelDistance;     // 6 bits
    uint8_t channelLayout;       // 3 bits, used only when sourceFormat == bit value 10

} MasaDescriptiveMetadata;

typedef struct
{
    uint16_t directionIndex[SUBFRAME_COUNT][MASA_BAND_COUNT];
    uint8_t directToTotalRatio[SUBFRAME_COUNT][MASA_BAND_COUNT];
    uint8_t spreadCoherence[SUBFRAME_COUNT][MASA_BAND_COUNT];
} MasaDirectionSpatialMetadata;

typedef struct
{
    uint8_t diffuseToTotalRatio[SUBFRAME_COUNT][MASA_BAND_COUNT];
    uint8_t surroundCoherence[SUBFRAME_COUNT][MASA_BAND_COUNT];
    uint8_t remainderToTotalRatio[SUBFRAME_COUNT][MASA_BAND_COUNT];
} MasaCommonSpatialMetadata;

typedef struct
{
    MasaDescriptiveMetadata descriptiveMeta;
    MasaDirectionSpatialMetadata directionalSpatialMeta[MASA_MAXIMUM_DIRECTIONS];
    MasaCommonSpatialMetadata commonSpatialMeta;

} MasaMetadata;


extern const unsigned short bandBottomLimits[MASA_BAND_COUNT];
extern const unsigned short bandTopLimits[MASA_BAND_COUNT];

extern const float eigenMobileEqB[3];
extern const float eigenMobileEqA[3];

extern const uint8_t defaultFormatDescriptor[8];

#endif //MASAANALYZER_MASACOMMONHEADER_H
