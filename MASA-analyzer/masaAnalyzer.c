/*--------------------------------------------------------------------------------*
 * MASA analyzer                                                                  *
 * ----------------------------------                                             *
 * (C) 2024 Nokia Technologies Ltd. See LICENSE.md for license.                  *
 *                                                                                *
 *--------------------------------------------------------------------------------*/

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "cldfb.h"
#include "masaCommonHeader.h"
#include "sphericalIndexing.h"

// Local constants
#define EIGEN_TO_FOA_FILTER_FILENAME  "eigen_to_foa_cldfb_domain_filters.bin"
#define EIGEN_TO_HOA2_FILTER_FILENAME "eigen_to_hoa2_cldfb_domain_filters.bin"
#define SECTOR_FILTER_FILENAME        "sector_filters.bin"
#define EIGEN_MOBILE_LEFT_TRANSPORT 13u
#define EIGEN_MOBILE_RIGHT_TRANSPORT 29u
#define EIGEN_MOBILE_CENTER_TRANSPORT 0u

static const float smoothAlphaSh = 0.9048f;


// Local structures
typedef enum
{
    ANALYZER_1DIR_ANALYSIS,  // Performs single direction analysis
    ANALYZER_2DIR_ANALYSIS   // Performs two direction analysis using sector-based processing
} AnalyzerDirectionMode;

typedef enum
{
    ANALYZER_EIGENMIKE,
    ANALYZER_EIGENMIKE_PLUS_TRANSPORT,
    ANALYZER_FOA_PLUS_TRANSPORT,
    ANALYZER_HOA2_PLUS_TRANSPORT,
    ANALYZER_FOA,
    ANALYZER_HOA2
} AnalyzerInputFormat;

typedef struct
{
    AnalyzerInputFormat inputFormat;
    short numberOfInputChannels;
    AnalyzerDirectionMode dirMode;
    short directionCount;
    TransportType transportType;
    short transportCount;

    cldfb_struct* cldfb[EIGENMIKE_CHANNEL_COUNT];
    SphericalGridData sphericalGridData;

    float eigenToFoaFilterReal[CLDFB_BAND_COUNT][EIGENMIKE_CHANNEL_COUNT][FOA_CHANNEL_COUNT];
    float eigenToFoaFilterImag[CLDFB_BAND_COUNT][EIGENMIKE_CHANNEL_COUNT][FOA_CHANNEL_COUNT];

    float eigenToHoa2FilterReal[CLDFB_BAND_COUNT][EIGENMIKE_CHANNEL_COUNT][HOA2_CHANNEL_COUNT];
    float eigenToHoa2FilterImag[CLDFB_BAND_COUNT][EIGENMIKE_CHANNEL_COUNT][HOA2_CHANNEL_COUNT];

    float sectorFilters[HOA2_CHANNEL_COUNT][FOA_CHANNEL_COUNT][2];

    FILE* inputFile;
    FILE* outputPcmFile;
    FILE* outputMetadataFile;

    int16_t inputRaw[FRAME_SIZE*(EIGENMIKE_CHANNEL_COUNT+2)];
    float inputData[FRAME_SIZE*(EIGENMIKE_CHANNEL_COUNT+2)];

    float eigenSignalReal[EIGENMIKE_CHANNEL_COUNT][CLDFB_SLOT_COUNT][CLDFB_BAND_COUNT];
    float eigenSignalImag[EIGENMIKE_CHANNEL_COUNT][CLDFB_SLOT_COUNT][CLDFB_BAND_COUNT];

    float*** foaSignalReal;
    float*** foaSignalImag;

    float hoa2SignalReal[HOA2_CHANNEL_COUNT][CLDFB_SLOT_COUNT][CLDFB_BAND_COUNT];
    float hoa2SignalImag[HOA2_CHANNEL_COUNT][CLDFB_SLOT_COUNT][CLDFB_BAND_COUNT];

    float*** sectorSignal1Real;
    float*** sectorSignal1Imag;
    float*** sectorSignal2Real;
    float*** sectorSignal2Imag;

    float subframeEnergy[MASA_MAXIMUM_DIRECTIONS][MASA_BAND_COUNT];
    float subframeShEnergy[MASA_BAND_COUNT][4];
    float subframeIntensity[MASA_MAXIMUM_DIRECTIONS][MASA_BAND_COUNT][3];
    float subframeIntensityLength[MASA_MAXIMUM_DIRECTIONS][MASA_BAND_COUNT];

    float azimuth[MASA_MAXIMUM_DIRECTIONS][SUBFRAME_COUNT][MASA_BAND_COUNT];
    float elevation[MASA_MAXIMUM_DIRECTIONS][SUBFRAME_COUNT][MASA_BAND_COUNT];
    float directToTotalRatio[MASA_MAXIMUM_DIRECTIONS][SUBFRAME_COUNT][MASA_BAND_COUNT];
    float spreadCoherence[MASA_MAXIMUM_DIRECTIONS][SUBFRAME_COUNT][MASA_BAND_COUNT];
    float diffuseToTotalRatio[SUBFRAME_COUNT][MASA_BAND_COUNT];
    float surroundCoherence[SUBFRAME_COUNT][MASA_BAND_COUNT];

    int16_t* transportSignal;
    MasaMetadata outputMetadata;

    float eigenMobileEqHistory[2][2];

    float transportSignalReal[MAX_TRANSPORT_COUNT][CLDFB_SLOT_COUNT][CLDFB_BAND_COUNT];
    float transportSignalImag[MAX_TRANSPORT_COUNT][CLDFB_SLOT_COUNT][CLDFB_BAND_COUNT];

    float transEnergy[CLDFB_BAND_COUNT];
    float refEnergy[CLDFB_BAND_COUNT];
    float generalCoherence[MASA_BAND_COUNT];

} AnalyzerData;


// Local function prototypes
static void openAnalyzer(AnalyzerData* data);
static void closeAnalyzer(AnalyzerData* data);
static void createAnalysisSignals(AnalyzerData* data);
static void performSpatialAnalysis(AnalyzerData* data);
static void createTransportSignal(AnalyzerData* data);
static void formDescriptiveMetadata(AnalyzerData* data);
static void formSpatialMetadata(AnalyzerData* data);
static void writeTransportFrame(FILE* file, const int16_t* signal, size_t signalLength);
static void writeMetadataFrame(FILE* file, MasaMetadata* meta);
static void clearAnalysisBuffers(AnalyzerData* data);

static void complexMult(float* outRe, float* outIm, float in1Re, float in1Im, float in2Re, float in2Im);
static float fclampf(float in, float lowLim, float highLim);
static size_t write_uint16_t_LE(const uint16_t* data, size_t nitems, FILE* stream);
static void printUsage();


/*
 * masaAnalyzer main function
 *
 * Arguments:
 * transport channel mode
 * direction analysis mode
 * input format
 * (optional descriptive meta for external transports)
 * input file
 * output pcm
 * output meta
 */
int main(int argc, char* argv[])
{
    // Check correct number of command line arguments
    if (argc < 7)
    {
        printUsage();
        exit(EXIT_FAILURE);
    }

    AnalyzerData analyzerData = {0};

    size_t parseCount = 1;

    // Parse transport mode
    if (strcmp(argv[parseCount], "-mono") == 0)
    {
        analyzerData.transportType = MASA_TRANSPORT_MONO;
        analyzerData.transportCount = 1;
    }
    else if (strcmp(argv[parseCount], "-stereo") == 0)
    {
        analyzerData.transportType = MASA_TRANSPORT_STEREO;
        analyzerData.transportCount = 2;
    }
    else
    {
        printUsage();
        exit(EXIT_FAILURE);
    }
    parseCount++;

    // Parse 1dir vs. 2dir analysis mode
    if (strcmp(argv[parseCount], "-1dir") == 0)
    {
        analyzerData.dirMode = ANALYZER_1DIR_ANALYSIS;
        analyzerData.directionCount = 1;
    }
    else if (strcmp(argv[parseCount], "-2dir") == 0)
    {
        analyzerData.dirMode = ANALYZER_2DIR_ANALYSIS;
        analyzerData.directionCount = 2;
    }
    else
    {
        printUsage();
        exit(EXIT_FAILURE);
    }
    parseCount++;

    // Parse input format
    if (strcmp(argv[parseCount], "-eigen") == 0)
    {
        analyzerData.inputFormat = ANALYZER_EIGENMIKE;
        analyzerData.numberOfInputChannels = 32;
    }
    else if (strcmp(argv[parseCount], "-eigenext") == 0)
    {
        analyzerData.inputFormat = ANALYZER_EIGENMIKE_PLUS_TRANSPORT;
        analyzerData.numberOfInputChannels = 32 + analyzerData.transportCount;
    }
    else if (strcmp(argv[parseCount], "-foaext") == 0)
    {
        analyzerData.inputFormat = ANALYZER_FOA_PLUS_TRANSPORT;
        analyzerData.numberOfInputChannels = 4 + analyzerData.transportCount;
    }
    else if (strcmp(argv[parseCount], "-hoa2ext") == 0)
    {
        analyzerData.inputFormat = ANALYZER_HOA2_PLUS_TRANSPORT;
        analyzerData.numberOfInputChannels = 9 + analyzerData.transportCount;
    }
    else if (strcmp(argv[parseCount], "-foa") == 0)
    {
        analyzerData.inputFormat = ANALYZER_FOA;
        analyzerData.numberOfInputChannels = 4;
    }
    else if (strcmp(argv[parseCount], "-hoa2") == 0)
    {
        analyzerData.inputFormat = ANALYZER_HOA2;
        analyzerData.numberOfInputChannels = 9;
    }
    else
    {
        printUsage();
        exit(EXIT_FAILURE);
    }
    parseCount++;

    // Parse optional descriptive metadata for external transports
    if (analyzerData.inputFormat == ANALYZER_EIGENMIKE_PLUS_TRANSPORT || analyzerData.inputFormat == ANALYZER_FOA_PLUS_TRANSPORT || analyzerData.inputFormat == ANALYZER_HOA2_PLUS_TRANSPORT)
    {
        if (strcmp(argv[parseCount], "-descmeta") == 0)
        {
            // Do another check that we have enough parameters
            if (argc != 11)
            {
                printUsage();
                exit(EXIT_FAILURE);
            }
            char* temp;
            parseCount++;
            analyzerData.outputMetadata.descriptiveMeta.transportDefinition = (uint8_t) strtol(argv[parseCount], &temp, 2);
            parseCount++;
            analyzerData.outputMetadata.descriptiveMeta.channelAngle = (uint8_t) strtol(argv[parseCount], &temp, 2);
            parseCount++;
            analyzerData.outputMetadata.descriptiveMeta.channelDistance = (uint8_t) strtol(argv[parseCount], &temp, 2);
            parseCount++;
        }
        else
        {
            // Do another check that we have enough parameters
            if (argc != 7)
            {
                printUsage();
                exit(EXIT_FAILURE);
            }
            analyzerData.outputMetadata.descriptiveMeta.transportDefinition = 0;
            analyzerData.outputMetadata.descriptiveMeta.channelAngle = 0;
            analyzerData.outputMetadata.descriptiveMeta.channelDistance = 0;
        }
    }

    // Give warning if 2dir analysis is not possible
    if ((analyzerData.inputFormat == ANALYZER_FOA_PLUS_TRANSPORT || analyzerData.inputFormat == ANALYZER_FOA) && analyzerData.dirMode == ANALYZER_2DIR_ANALYSIS )
    {
        printf("Warning!\n");
        printf("2 direction analysis requires HOA2 or Eigenmike input. Changing to 1 direction analysis.\n");
        analyzerData.dirMode = ANALYZER_1DIR_ANALYSIS;
        analyzerData.directionCount = 1;
    }

    // Simple sanity check that the input file is not a wave file.
    // This is to avoid usage errors.
    const char *extension = strrchr(argv[parseCount], '.');
    if (strcmp(extension, ".wav") == 0)
    {
        printf("Warning!\n");
        printf("Wave file extension detected. Please ensure that you are using a converted PCM file as input.\n");
    }

    analyzerData.inputFile = fopen(argv[parseCount], "rb");
    if (analyzerData.inputFile == NULL)
    {
        fprintf(stderr, "Failed to open input file for reading\n");
        exit(EXIT_FAILURE);
    }
    parseCount++;

    analyzerData.outputPcmFile = fopen(argv[parseCount], "wb");
    if (analyzerData.outputPcmFile == NULL)
    {
        fprintf(stderr, "Failed to open output pcm file for writing\n");
        exit(EXIT_FAILURE);
    }
    parseCount++;

    analyzerData.outputMetadataFile = fopen(argv[parseCount], "wb");
    if (analyzerData.outputMetadataFile == NULL)
    {
        fprintf(stderr, "Failed to open output metadata file for writing\n");
        exit(EXIT_FAILURE);
    }

    openAnalyzer(&analyzerData);

    size_t readCount = 0;

    // Main processing loop
    while ((readCount = fread(analyzerData.inputRaw, sizeof(int16_t), FRAME_SIZE * analyzerData.numberOfInputChannels, analyzerData.inputFile)) > 0)
    {
        // Check if we got a full frame. This is a rudimentary check to verify
        // we don't process something wrong. Ensure that input is correct length
        // suppress this warning.
        if (readCount < FRAME_SIZE * analyzerData.numberOfInputChannels)
        {
            printf("Warning, input file was not aligned with %d channels and 960 time samples.\n", analyzerData.numberOfInputChannels);
            printf("Truncating the last frame out.\n");
            break;
        }

        // Convert to float and scale range to [-1,1[
        for (size_t n = 0; n < readCount; ++n)
        {
            analyzerData.inputData[n] = (float) analyzerData.inputRaw[n] / (float) INT16_MAX;
        }

        createAnalysisSignals(&analyzerData);
        performSpatialAnalysis(&analyzerData);
        createTransportSignal(&analyzerData);
        formDescriptiveMetadata(&analyzerData);
        formSpatialMetadata(&analyzerData);

        writeTransportFrame(analyzerData.outputPcmFile, analyzerData.transportSignal, FRAME_SIZE * analyzerData.transportCount);
        writeMetadataFrame(analyzerData.outputMetadataFile, &(analyzerData.outputMetadata));

        clearAnalysisBuffers(&analyzerData);
    }

    closeAnalyzer(&analyzerData);

    return 0;
}


static void openAnalyzer(AnalyzerData* data)
{
    // Load Eigen to FOA filters
    FILE* file = fopen(EIGEN_TO_FOA_FILTER_FILENAME, "rb");
    if (file == NULL)
    {
        fprintf(stderr, "Failed to open Eigenmike to FOA filter file for reading\n");
        exit(EXIT_FAILURE);
    }

    float readTemp[CLDFB_BAND_COUNT*EIGENMIKE_CHANNEL_COUNT*HOA2_CHANNEL_COUNT*2]; // Reserve space for largest data
    unsigned long out = fread(&readTemp, sizeof(float), CLDFB_BAND_COUNT*EIGENMIKE_CHANNEL_COUNT*FOA_CHANNEL_COUNT*2, file);
    fclose(file);
    if (out != CLDFB_BAND_COUNT*EIGENMIKE_CHANNEL_COUNT*FOA_CHANNEL_COUNT*2)
    {
        fprintf(stderr, "Eigen to FOA filter read failed.\n");
        exit(EXIT_FAILURE);
    }

    float* readPointer = &(readTemp[0]);
    for (size_t n = 0; n < CLDFB_BAND_COUNT; n++)
    {
        for (size_t m = 0; m < EIGENMIKE_CHANNEL_COUNT; m++)
        {
            for (size_t k = 0; k < FOA_CHANNEL_COUNT; k++)
            {
                data->eigenToFoaFilterReal[n][m][k] = *readPointer++;
                data->eigenToFoaFilterImag[n][m][k] = *readPointer++;
            }
        }
    }

    // Load Eigen to HOA2 filters & sector filters if working with 2 direction analysis
    if (data->dirMode == ANALYZER_2DIR_ANALYSIS)
    {
        file = fopen(EIGEN_TO_HOA2_FILTER_FILENAME, "rb");
        if (file == NULL)
        {
            fprintf(stderr, "Failed to open Eigenmike to HOA2 filter file for reading\n");
            exit(EXIT_FAILURE);
        }

        out = fread(&readTemp, sizeof(float), CLDFB_BAND_COUNT*EIGENMIKE_CHANNEL_COUNT*HOA2_CHANNEL_COUNT*2, file);
        if (out != CLDFB_BAND_COUNT*EIGENMIKE_CHANNEL_COUNT*HOA2_CHANNEL_COUNT*2)
        {
            fprintf(stderr, "Eigen to HOA2 filter read failed.\n");
            exit(EXIT_FAILURE);
        }
        fclose(file);

        readPointer = &(readTemp[0]);
        for (size_t n = 0; n < CLDFB_BAND_COUNT; n++)
        {
            for (size_t m = 0; m < EIGENMIKE_CHANNEL_COUNT; m++)
            {
                for (size_t k = 0; k < HOA2_CHANNEL_COUNT; k++)
                {
                    data->eigenToHoa2FilterReal[n][m][k] = *readPointer++;
                    data->eigenToHoa2FilterImag[n][m][k] = *readPointer++;
                }
            }
        }

        file = fopen(SECTOR_FILTER_FILENAME, "rb");
        if (file == NULL)
        {
            fprintf(stderr, "Failed to open sector filter file for reading\n");
            exit(EXIT_FAILURE);
        }

        out = fread(&readTemp, sizeof(float), HOA2_CHANNEL_COUNT*FOA_CHANNEL_COUNT*2, file);
        if (out != HOA2_CHANNEL_COUNT*FOA_CHANNEL_COUNT*2)
        {
            fprintf(stderr, "Sector filter read failed.\n");
            exit(EXIT_FAILURE);
        }
        fclose(file);

        readPointer = &(readTemp[0]);
        for (size_t n = 0; n < HOA2_CHANNEL_COUNT; n++)
        {
            for (size_t m = 0; m < FOA_CHANNEL_COUNT; m++)
            {
                for (size_t k = 0; k < 2; k++)
                {
                    data->sectorFilters[n][m][k] = *readPointer++;
                }
            }
        }

    }

    for (size_t n = 0; n < EIGENMIKE_CHANNEL_COUNT; n++)
    {
        data->cldfb[n] = (cldfb_struct*) malloc(sizeof(cldfb_struct));
        cldfbInit(data->cldfb[n], CLDFB_1_25_MILLISECONDS, CLDFB_FORWARD);
    }

    data->eigenMobileEqHistory[0][0] = 0.0f;
    data->eigenMobileEqHistory[0][1] = 0.0f;
    data->eigenMobileEqHistory[1][0] = 0.0f;
    data->eigenMobileEqHistory[1][1] = 0.0f;

    data->transportSignal = (int16_t*) calloc(FRAME_SIZE * data->transportCount, sizeof(int16_t));

    generateSphericalGrid(&(data->sphericalGridData));

    // Reserve storage structures
    data->foaSignalReal = (float***) malloc(FOA_CHANNEL_COUNT * sizeof(float**));
    data->foaSignalImag = (float***) malloc(FOA_CHANNEL_COUNT * sizeof(float**));
    if (data->foaSignalReal != NULL && data->foaSignalImag != NULL)
    {
        for (size_t n = 0; n < FOA_CHANNEL_COUNT; n++)
        {
            data->foaSignalReal[n] = (float**) malloc(CLDFB_SLOT_COUNT * sizeof(float*));
            data->foaSignalImag[n] = (float**) malloc(CLDFB_SLOT_COUNT * sizeof(float*));
            if (data->foaSignalReal[n] != NULL && data->foaSignalImag[n] != NULL)
            {
                for (size_t m = 0; m < CLDFB_SLOT_COUNT; m++)
                {
                    data->foaSignalReal[n][m] = (float*) calloc(CLDFB_BAND_COUNT, sizeof(float));
                    data->foaSignalImag[n][m] = (float*) calloc(CLDFB_BAND_COUNT, sizeof(float));
                }
            }
        }
    }

    if (data->dirMode == ANALYZER_2DIR_ANALYSIS)
    {
        data->sectorSignal1Real = (float***) malloc(FOA_CHANNEL_COUNT * sizeof(float**));
        data->sectorSignal1Imag = (float***) malloc(FOA_CHANNEL_COUNT * sizeof(float**));
        if (data->sectorSignal1Real != NULL && data->sectorSignal1Imag != NULL)
        {
            for (size_t n = 0; n < FOA_CHANNEL_COUNT; n++)
            {
                data->sectorSignal1Real[n] = (float**) malloc(CLDFB_SLOT_COUNT * sizeof(float*));
                data->sectorSignal1Imag[n] = (float**) malloc(CLDFB_SLOT_COUNT * sizeof(float*));
                if (data->sectorSignal1Real[n] != NULL && data->sectorSignal1Imag[n] != NULL)
                {
                    for (size_t m = 0; m < CLDFB_SLOT_COUNT; m++)
                    {
                        data->sectorSignal1Real[n][m] = (float*) calloc(CLDFB_BAND_COUNT, sizeof(float));
                        data->sectorSignal1Imag[n][m] = (float*) calloc(CLDFB_BAND_COUNT, sizeof(float));
                    }
                }
            }
        }
        data->sectorSignal2Real = (float***) malloc(FOA_CHANNEL_COUNT * sizeof(float**));
        data->sectorSignal2Imag = (float***) malloc(FOA_CHANNEL_COUNT * sizeof(float**));
        if (data->sectorSignal2Real != NULL && data->sectorSignal2Imag != NULL)
        {
            for (size_t n = 0; n < FOA_CHANNEL_COUNT; n++)
            {
                data->sectorSignal2Real[n] = (float**) malloc(CLDFB_SLOT_COUNT * sizeof(float*));
                data->sectorSignal2Imag[n] = (float**) malloc(CLDFB_SLOT_COUNT * sizeof(float*));
                if (data->sectorSignal2Real[n] != NULL && data->sectorSignal2Imag[n] != NULL)
                {
                    for (size_t m = 0; m < CLDFB_SLOT_COUNT; m++)
                    {
                        data->sectorSignal2Real[n][m] = (float*) calloc(CLDFB_BAND_COUNT, sizeof(float));
                        data->sectorSignal2Imag[n][m] = (float*) calloc(CLDFB_BAND_COUNT, sizeof(float));
                    }
                }
            }
        }
    }
}


static void closeAnalyzer(AnalyzerData* data)
{
    fclose(data->inputFile);
    fclose(data->outputPcmFile);
    fclose(data->outputMetadataFile);

    for (size_t n = 0; n < EIGENMIKE_CHANNEL_COUNT; n++)
    {
        cldfbDestroy(data->cldfb[n]);
        free(data->cldfb[n]);
        data->cldfb[n] = NULL;
    }

    free(data->transportSignal);

    // Free storage structures
    if (data->foaSignalReal != NULL && data->foaSignalImag != NULL)
    {
        for (size_t n = 0; n < FOA_CHANNEL_COUNT; n++)
        {
            if (data->foaSignalReal[n] != NULL && data->foaSignalImag[n] != NULL)
            {
                for (size_t m = 0; m < CLDFB_SLOT_COUNT; m++)
                {
                    free(data->foaSignalReal[n][m]);
                    free(data->foaSignalImag[n][m]);
                }
            }
            free(data->foaSignalReal[n]);
            free(data->foaSignalImag[n]);
        }
    }
    free(data->foaSignalReal);
    free(data->foaSignalImag);

    if (data->dirMode == ANALYZER_2DIR_ANALYSIS)
    {
        if (data->sectorSignal1Real != NULL && data->sectorSignal1Imag != NULL)
        {
            for (size_t n = 0; n < FOA_CHANNEL_COUNT; n++)
            {
                if (data->sectorSignal1Real[n] != NULL && data->sectorSignal1Imag[n] != NULL)
                {
                    for (size_t m = 0; m < CLDFB_SLOT_COUNT; m++)
                    {
                        free(data->sectorSignal1Real[n][m]);
                        free(data->sectorSignal1Imag[n][m]);
                    }
                }
                free(data->sectorSignal1Real[n]);
                free(data->sectorSignal1Imag[n]);
            }
        }
        free(data->sectorSignal1Real);
        free(data->sectorSignal1Imag);

        if (data->sectorSignal2Real != NULL && data->sectorSignal2Imag != NULL)
        {
            for (size_t n = 0; n < FOA_CHANNEL_COUNT; n++)
            {
                if (data->sectorSignal2Real[n] != NULL && data->sectorSignal2Imag[n] != NULL)
                {
                    for (size_t m = 0; m < CLDFB_SLOT_COUNT; m++)
                    {
                        free(data->sectorSignal2Real[n][m]);
                        free(data->sectorSignal2Imag[n][m]);
                    }
                }
                free(data->sectorSignal2Real[n]);
                free(data->sectorSignal2Imag[n]);
            }
        }
        free(data->sectorSignal2Real);
        free(data->sectorSignal2Imag);
    }
}

static void performSpatialAnalysis(AnalyzerData* data)
{
    // Perform spatial analysis in subframes to obtain 5 ms time resolution
    for (size_t sf = 0; sf < SUBFRAME_COUNT; sf++)
    {
        // Compute SH energy
        float meanCoherence = 0.0f;
        for (size_t band = 0; band < MASA_BAND_COUNT; band++)
        {
            float tempBandShEnergy[4] = {0};
            for (size_t bin = bandBottomLimits[band]; bin < bandTopLimits[band]; bin++)
            {
                for (size_t c = 0; c < FOA_CHANNEL_COUNT; c++)
                {
                    for (size_t n = 0; n < CLDFB_SLOTS_PER_SUBFRAME; n++)
                    {
                        size_t ts = sf * CLDFB_SLOTS_PER_SUBFRAME + n;
                        tempBandShEnergy[c] += data->foaSignalReal[c][ts][bin]*data->foaSignalReal[c][ts][bin];
                        tempBandShEnergy[c] += data->foaSignalImag[c][ts][bin]*data->foaSignalImag[c][ts][bin];
                    }
                }
            }

            for (size_t c = 0; c < FOA_CHANNEL_COUNT; c++)
            {
                data->subframeShEnergy[band][c] = tempBandShEnergy[c] * (1.0f - smoothAlphaSh) + data->subframeShEnergy[band][c] * smoothAlphaSh;
            }

            data->generalCoherence[band] =
                fclampf(1.0f - ((data->subframeShEnergy[band][1] + data->subframeShEnergy[band][2] + data->subframeShEnergy[band][3]) / data->subframeShEnergy[band][0]),
                    0.0f, 1.0f);

            meanCoherence += data->generalCoherence[band];
        }
        meanCoherence /= (float)MASA_BAND_COUNT;
        data->generalCoherence[0] = meanCoherence;

        // Perform directional analysis
        for (size_t dir = 0; dir < data->directionCount; dir++)
        {
            float*** analysisSigReal = NULL;
            float*** analysisSigImag = NULL;

            if (data->dirMode == ANALYZER_1DIR_ANALYSIS)
            {
                analysisSigReal = data->foaSignalReal;
                analysisSigImag = data->foaSignalImag;
            }
            else if (data->dirMode == ANALYZER_2DIR_ANALYSIS && dir == 0)
            {
                analysisSigReal = data->sectorSignal1Real;
                analysisSigImag = data->sectorSignal1Imag;
            }
            else if (data->dirMode == ANALYZER_2DIR_ANALYSIS && dir == 1)
            {
                analysisSigReal = data->sectorSignal2Real;
                analysisSigImag = data->sectorSignal2Imag;
            }

            // Compute energy
            for (size_t band = 0; band < MASA_BAND_COUNT; band++)
            {
                float tempBandEnergy = 0;
                for (size_t bin = bandBottomLimits[band]; bin < bandTopLimits[band]; bin++)
                {
                    for (size_t c = 0; c < FOA_CHANNEL_COUNT; c++)
                    {
                        for (size_t n = 0; n < CLDFB_SLOTS_PER_SUBFRAME; n++)
                        {
                            size_t ts = sf * CLDFB_SLOTS_PER_SUBFRAME + n;
                            tempBandEnergy += analysisSigReal[c][ts][bin]*analysisSigReal[c][ts][bin];
                            tempBandEnergy += analysisSigImag[c][ts][bin]*analysisSigImag[c][ts][bin];
                        }
                    }
                }
                tempBandEnergy /= 2.0f;
                data->subframeEnergy[dir][band] = tempBandEnergy;
            }

            // Compute intensity
            for (size_t band = 0; band < MASA_BAND_COUNT; band++)
            {
                float tempBandIntensity[3] = {0};
                for (size_t bin = bandBottomLimits[band]; bin < bandTopLimits[band]; bin++)
                {
                    for (size_t n = 0; n < CLDFB_SLOTS_PER_SUBFRAME; n++)
                    {
                        size_t ts = sf * CLDFB_SLOTS_PER_SUBFRAME + n;

                        // Real[W * conj(X)]
                        tempBandIntensity[0] += analysisSigReal[0][ts][bin]*analysisSigReal[3][ts][bin];
                        tempBandIntensity[0] += analysisSigImag[0][ts][bin]*analysisSigImag[3][ts][bin];

                        // Real[W * conj(Y)]
                        tempBandIntensity[1] += analysisSigReal[0][ts][bin]*analysisSigReal[1][ts][bin];
                        tempBandIntensity[1] += analysisSigImag[0][ts][bin]*analysisSigImag[1][ts][bin];

                        // Real[W * conj(Z)]
                        tempBandIntensity[2] += analysisSigReal[0][ts][bin]*analysisSigReal[2][ts][bin];
                        tempBandIntensity[2] += analysisSigImag[0][ts][bin]*analysisSigImag[2][ts][bin];
                    }
                }
                data->subframeIntensity[dir][band][0] = tempBandIntensity[0];
                data->subframeIntensity[dir][band][1] = tempBandIntensity[1];
                data->subframeIntensity[dir][band][2] = tempBandIntensity[2];

                data->subframeIntensityLength[dir][band] = sqrtf(
                    tempBandIntensity[0] * tempBandIntensity[0]
                    + tempBandIntensity[1] * tempBandIntensity[1]
                    + tempBandIntensity[2] * tempBandIntensity[2]);
            }

            // Compute azimuth, elevation, and energy ratio
            for (size_t band = 0; band < MASA_BAND_COUNT; band++)
            {
                // Azimuth
                data->azimuth[dir][sf][band] = atan2f(data->subframeIntensity[dir][band][1],
                                                      data->subframeIntensity[dir][band][0]) / MASA_PI * 180;

                // Elevation
                data->elevation[dir][sf][band] = atan2f(data->subframeIntensity[dir][band][2],
                                                        sqrtf(data->subframeIntensity[dir][band][0]*data->subframeIntensity[dir][band][0]
                                                              + data->subframeIntensity[dir][band][1]*data->subframeIntensity[dir][band][1]))
                                                 / MASA_PI * 180;

                // Direct-to-total energy ratio
                data->directToTotalRatio[dir][sf][band] = fclampf(data->subframeIntensityLength[dir][band]
                                                                  / (data->subframeEnergy[dir][band] + MASA_EPSILON), 0.0f, 1.0f);


                data->spreadCoherence[dir][sf][band] = data->directToTotalRatio[dir][sf][band] * data->generalCoherence[band];
            }
        }

        // Formulate surround coherence
        for (size_t band = 0; band < MASA_BAND_COUNT; band++)
        {
            data->surroundCoherence[sf][band] = 1.0f;
            float meanRatio = 0.0f;
            for (size_t dir = 0; dir < data->directionCount; dir++)
            {
                meanRatio += data->directToTotalRatio[dir][sf][band] / (float) data->directionCount;
            }
            data->surroundCoherence[sf][band] = (1.0f - meanRatio) * data->generalCoherence[band];
        }

        // Scale ratios if working with two directions as above analysis assumes separate diffuse ratios instead of a single one
        if (data->dirMode == ANALYZER_2DIR_ANALYSIS)
        {
            for (size_t band = 0; band < MASA_BAND_COUNT; band++)
            {
                float energyDir[MASA_MAXIMUM_DIRECTIONS];
                float energyNonDir[MASA_MAXIMUM_DIRECTIONS];
                float totalEnergy = 0.0f;
                for (size_t dir = 0; dir < data->directionCount; dir++)
                {
                    energyDir[dir] = data->directToTotalRatio[dir][sf][band] * data->subframeEnergy[dir][band];
                    energyNonDir[dir] = (1.0f - data->directToTotalRatio[dir][sf][band]) * data->subframeEnergy[dir][band];
                    totalEnergy += energyDir[dir] + energyNonDir[dir];
                }

                if (totalEnergy > 0.0f)
                {
                    data->diffuseToTotalRatio[sf][band] = 1.0f;
                    for (size_t dir = 0; dir < data->directionCount; dir++)
                    {
                        data->directToTotalRatio[dir][sf][band] = energyDir[dir] / totalEnergy;
                        data->diffuseToTotalRatio[sf][band] -= data->directToTotalRatio[dir][sf][band];
                    }
                }
                else
                {
                    data->diffuseToTotalRatio[sf][band] = 1.0f;
                    for (size_t dir = 0; dir < data->directionCount; dir++)
                    {
                        data->directToTotalRatio[dir][sf][band] = 0.0f;
                    }
                }
            }
        }
        else if (data->dirMode == ANALYZER_1DIR_ANALYSIS)
        {
            for (size_t band = 0; band < MASA_BAND_COUNT; band++)
            {
                data->diffuseToTotalRatio[sf][band] = 1.0f - data->directToTotalRatio[0][sf][band];
            }
        }
    }
}

static void createAnalysisSignals(AnalyzerData* data)
{
    // Obtain signals for analysis in CLDFB domain
    if (data->inputFormat == ANALYZER_EIGENMIKE || data->inputFormat == ANALYZER_EIGENMIKE_PLUS_TRANSPORT)
    {
        // Go to frequency domain
        for (size_t n = 0; n < EIGENMIKE_CHANNEL_COUNT; n++)
        {
            float sigTime[FRAME_SIZE];
            for (size_t m = 0; m < FRAME_SIZE; m++)
            {
                sigTime[m] = data->inputData[m * data->numberOfInputChannels + n];
            }

            for (size_t m = 0; m < CLDFB_SLOT_COUNT; m++)
            {
                cldfbForward(data->cldfb[n], &sigTime[CLDFB_BAND_COUNT*m], data->eigenSignalReal[n][m], data->eigenSignalImag[n][m]);
            }
        }

        // Obtain FOA signals
        for (size_t n = 0; n < CLDFB_BAND_COUNT; n++)
        {
            for (size_t m = 0; m < EIGENMIKE_CHANNEL_COUNT; m++)
            {
                for (size_t k = 0; k < FOA_CHANNEL_COUNT; k++)
                {
                    for (size_t l = 0; l < CLDFB_SLOT_COUNT; l++)
                    {
                        float tempReal = 0;
                        float tempImag = 0;
                        complexMult(&tempReal, &tempImag,
                                    data->eigenSignalReal[m][l][n],
                                    data->eigenSignalImag[m][l][n],
                                    data->eigenToFoaFilterReal[n][m][k],
                                    data->eigenToFoaFilterImag[n][m][k]);

                        data->foaSignalReal[k][l][n] += tempReal;
                        data->foaSignalImag[k][l][n] += tempImag;
                    }
                }
            }
        }

        // If working with 2 direction analysis, also form HOA2 signals and the two sector signals
        if (data->dirMode == ANALYZER_2DIR_ANALYSIS)
        {
            for (size_t n = 0; n < CLDFB_BAND_COUNT; n++)
            {
                for (size_t m = 0; m < EIGENMIKE_CHANNEL_COUNT; m++)
                {
                    for (size_t k = 0; k < HOA2_CHANNEL_COUNT; k++)
                    {
                        for (size_t l = 0; l < CLDFB_SLOT_COUNT; l++)
                        {
                            float tempReal = 0;
                            float tempImag = 0;
                            complexMult(
                                &tempReal,
                                &tempImag,
                                data->eigenSignalReal[m][l][n],
                                data->eigenSignalImag[m][l][n],
                                data->eigenToHoa2FilterReal[n][m][k],
                                data->eigenToHoa2FilterImag[n][m][k]);

                            data->hoa2SignalReal[k][l][n] += tempReal;
                            data->hoa2SignalImag[k][l][n] += tempImag;
                        }
                    }
                }
            }

            for (size_t n = 0; n < CLDFB_BAND_COUNT; n++)
            {
                for (size_t m = 0; m < CLDFB_SLOT_COUNT; m++)
                {
                    for (size_t foa_c = 0; foa_c < FOA_CHANNEL_COUNT; foa_c++)
                    {
                        for (size_t hoa2_c = 0; hoa2_c < HOA2_CHANNEL_COUNT; hoa2_c++)
                        {
                            data->sectorSignal1Real[foa_c][m][n] += data->hoa2SignalReal[hoa2_c][m][n] * data->sectorFilters[hoa2_c][foa_c][0];
                            data->sectorSignal1Imag[foa_c][m][n] += data->hoa2SignalImag[hoa2_c][m][n] * data->sectorFilters[hoa2_c][foa_c][0];
                            data->sectorSignal2Real[foa_c][m][n] += data->hoa2SignalReal[hoa2_c][m][n] * data->sectorFilters[hoa2_c][foa_c][1];
                            data->sectorSignal2Imag[foa_c][m][n] += data->hoa2SignalImag[hoa2_c][m][n] * data->sectorFilters[hoa2_c][foa_c][1];
                        }
                    }
                }
            }
        }
    }
    else // We are in FOA/HOA format directly for analysis
    {
        short numberOfAnalysisChannels;

        if (data->inputFormat == ANALYZER_FOA_PLUS_TRANSPORT || data->inputFormat == ANALYZER_HOA2_PLUS_TRANSPORT)
        {
            numberOfAnalysisChannels = data->numberOfInputChannels - data->transportCount;
        }
        else
        {
            numberOfAnalysisChannels = data->numberOfInputChannels;
        }

        // Go to frequency domain
        for (size_t n = 0; n < numberOfAnalysisChannels; n++)
        {
            float sigTime[FRAME_SIZE];
            for (size_t m = 0; m < FRAME_SIZE; m++)
            {
                sigTime[m] = data->inputData[m * data->numberOfInputChannels + n];
            }

            for (size_t m = 0; m < CLDFB_SLOT_COUNT; m++)
            {
                if (data->inputFormat == ANALYZER_FOA_PLUS_TRANSPORT || data->inputFormat == ANALYZER_FOA)
                {
                    cldfbForward(data->cldfb[n], &sigTime[CLDFB_BAND_COUNT*m], data->foaSignalReal[n][m], data->foaSignalImag[n][m]);
                }
                else if (data->inputFormat == ANALYZER_HOA2_PLUS_TRANSPORT || data->inputFormat == ANALYZER_HOA2)
                {
                    cldfbForward(data->cldfb[n], &sigTime[CLDFB_BAND_COUNT*m], data->hoa2SignalReal[n][m], data->hoa2SignalImag[n][m]);
                }
            }
        }

        if (data->inputFormat == ANALYZER_HOA2_PLUS_TRANSPORT || data->inputFormat == ANALYZER_HOA2)
        {
            // Copy FOA-part of HOA2 into FOA structure
            for (size_t n = 0; n < FOA_CHANNEL_COUNT; n++)
            {
                for (size_t m = 0; m < CLDFB_SLOT_COUNT; m++)
                {
                    for (size_t k = 0; k < CLDFB_BAND_COUNT; k++)
                    {
                        data->foaSignalReal[n][m][k] = data->hoa2SignalReal[n][m][k];
                        data->foaSignalImag[n][m][k] = data->hoa2SignalImag[n][m][k];
                    }
                }
            }

            // Create sector signals
            if (data->dirMode == ANALYZER_2DIR_ANALYSIS)
            {
                for (size_t n = 0; n < CLDFB_BAND_COUNT; n++)
                {
                    for (size_t m = 0; m < CLDFB_SLOT_COUNT; m++)
                    {
                        for (size_t foa_c = 0; foa_c < FOA_CHANNEL_COUNT; foa_c++)
                        {
                            for (size_t hoa2_c = 0; hoa2_c < HOA2_CHANNEL_COUNT; hoa2_c++)
                            {
                                data->sectorSignal1Real[foa_c][m][n] += data->hoa2SignalReal[hoa2_c][m][n] * data->sectorFilters[hoa2_c][foa_c][0];
                                data->sectorSignal1Imag[foa_c][m][n] += data->hoa2SignalImag[hoa2_c][m][n] * data->sectorFilters[hoa2_c][foa_c][0];
                                data->sectorSignal2Real[foa_c][m][n] += data->hoa2SignalReal[hoa2_c][m][n] * data->sectorFilters[hoa2_c][foa_c][1];
                                data->sectorSignal2Imag[foa_c][m][n] += data->hoa2SignalImag[hoa2_c][m][n] * data->sectorFilters[hoa2_c][foa_c][1];
                            }
                        }
                    }
                }
            }
        }
    }
}

static void createTransportSignal(AnalyzerData* data)
{
    if (data->inputFormat == ANALYZER_EIGENMIKE)
    {
        if (data->transportType == MASA_TRANSPORT_STEREO)
        {
            for (size_t n = 0; n < FRAME_SIZE; n++)
            {
                // Direct-Form II filter
                float v = 0;
                float y = 0;

                // Left channel
                v = data->inputData[n * EIGENMIKE_CHANNEL_COUNT + EIGEN_MOBILE_LEFT_TRANSPORT]
                    - eigenMobileEqA[1] * data->eigenMobileEqHistory[0][0]
                    - eigenMobileEqA[2] * data->eigenMobileEqHistory[0][1];

                y = eigenMobileEqB[0] * v + eigenMobileEqB[1] * data->eigenMobileEqHistory[0][0]
                    + eigenMobileEqB[2] * data->eigenMobileEqHistory[0][1];

                data->transportSignal[n * data->transportCount] = (int16_t) roundf(y * INT16_MAX);

                data->eigenMobileEqHistory[0][1] = data->eigenMobileEqHistory[0][0];
                data->eigenMobileEqHistory[0][0] = v;

                // Right channel
                v = data->inputData[n * EIGENMIKE_CHANNEL_COUNT + EIGEN_MOBILE_RIGHT_TRANSPORT]
                    - eigenMobileEqA[1] * data->eigenMobileEqHistory[1][0]
                    - eigenMobileEqA[2] * data->eigenMobileEqHistory[1][1];

                y = eigenMobileEqB[0] * v + eigenMobileEqB[1] * data->eigenMobileEqHistory[1][0]
                    + eigenMobileEqB[2] * data->eigenMobileEqHistory[1][1];

                data->transportSignal[n * data->transportCount + 1] = (int16_t) roundf(y * INT16_MAX);

                data->eigenMobileEqHistory[1][1] = data->eigenMobileEqHistory[1][0];
                data->eigenMobileEqHistory[1][0] = v;
            }
        }
        else // mono
        {
            for (size_t n = 0; n < FRAME_SIZE; n++)
            {
                // Direct-Form II filter
                float v = 0;
                float y = 0;

                v = data->inputData[n * EIGENMIKE_CHANNEL_COUNT + EIGEN_MOBILE_CENTER_TRANSPORT]
                    - eigenMobileEqA[1] * data->eigenMobileEqHistory[0][0]
                    - eigenMobileEqA[2] * data->eigenMobileEqHistory[0][1];

                y = eigenMobileEqB[0] * v + eigenMobileEqB[1] * data->eigenMobileEqHistory[0][0]
                    + eigenMobileEqB[2] * data->eigenMobileEqHistory[0][1];

                data->transportSignal[n * data->transportCount] = (int16_t) roundf(y * INT16_MAX);

                data->eigenMobileEqHistory[0][1] = data->eigenMobileEqHistory[0][0];
                data->eigenMobileEqHistory[0][0] = v;
            }
        }
    }
    else if (data->inputFormat == ANALYZER_FOA || data->inputFormat == ANALYZER_HOA2)
    {
        int16_t leftSample;
        int16_t rightSample;

        for (size_t n = 0; n < FRAME_SIZE; n++)
        {
            leftSample = 0.5 * (data->inputRaw[n * data->numberOfInputChannels] + data->inputRaw[n * data->numberOfInputChannels + 1]);
            rightSample = 0.5 * (data->inputRaw[n * data->numberOfInputChannels] - data->inputRaw[n * data->numberOfInputChannels + 1]);

            if (data->transportType == MASA_TRANSPORT_STEREO)
            {
                data->transportSignal[n * data->transportCount] = leftSample;
                data->transportSignal[n * data->transportCount + 1] = rightSample;
            }
            else
            {
                data->transportSignal[n * data->transportCount] = data->inputRaw[n * data->numberOfInputChannels];
            }
        }

    }
    else // We are in FOA/HOA2/Eigenmike + transport mode and we can simply copy signals for transport
    {
        short firstTransportIndex;
        if (data->inputFormat == ANALYZER_FOA_PLUS_TRANSPORT)
        {
            firstTransportIndex = 4;
        }
        else if (data->inputFormat == ANALYZER_HOA2_PLUS_TRANSPORT)
        {
            firstTransportIndex = 9;
        }
        else if (data->inputFormat == ANALYZER_EIGENMIKE_PLUS_TRANSPORT)
        {
            firstTransportIndex = 32;
        }

        for (size_t n = 0; n < FRAME_SIZE; n++)
        {
            data->transportSignal[n * data->transportCount] = data->inputRaw[n * data->numberOfInputChannels + firstTransportIndex];
        }
        if (data->transportType == MASA_TRANSPORT_STEREO)
        {
            for (size_t n = 0; n < FRAME_SIZE; n++)
            {
                data->transportSignal[n * data->transportCount + 1] = data->inputRaw[n * data->numberOfInputChannels + firstTransportIndex + 1];
            }
        }
    }
}


static void formDescriptiveMetadata(AnalyzerData* data)
{
    MasaDescriptiveMetadata* meta = &(data->outputMetadata.descriptiveMeta);

    // All values set here are already in the final bit value format instead of
    // the comparable functional format.

    // Set channel layout to zero as it is not used.
    meta->channelLayout = 0;

    // Set bits to 'IVASMASA'
    for (size_t n = 0; n < 8; n++)
    {
        meta->formatDescriptor[n] = defaultFormatDescriptor[n];
    }

    switch (data->dirMode)
    {
        case ANALYZER_1DIR_ANALYSIS:
            meta->numberOfDirections = 0x0; // One direction
            break;
        case ANALYZER_2DIR_ANALYSIS:
            meta->numberOfDirections = 0x1; // Two directions
            break;
    }

    switch (data->transportType)
    {
        case MASA_TRANSPORT_MONO:
            meta->numberOfChannels = 0x0; // One channel
            break;

        case MASA_TRANSPORT_STEREO:
            meta->numberOfChannels = 0x1; // Two channels
            break;
    }

    meta->sourceFormat = 0x1; // Microphone grid
    switch (data->inputFormat)
    {
        case ANALYZER_EIGENMIKE:
            if (data->transportType == MASA_TRANSPORT_STEREO)
            {
                meta->transportDefinition = 0x0; // Bit value 000, Unknown/Other as mic on surface has frequency dependent directivity
                meta->channelAngle = 0x1;        // Bit value 001, 90°
                meta->channelDistance = 0x1E;    // 0.0823 channel distance, Eigenmike is 0.084 diameter
            }
            break;
        case ANALYZER_EIGENMIKE_PLUS_TRANSPORT:
        case ANALYZER_HOA2_PLUS_TRANSPORT:
        case ANALYZER_FOA_PLUS_TRANSPORT:
            if (data->transportType == MASA_TRANSPORT_STEREO)
            {
                // These parameters are already set correctly in command line parsing
                // meta->transportDefinition
                // meta->channelAngle
                // meta->channelDistance
            }
            break;
        case ANALYZER_FOA:
            meta->sourceFormat = 0x3; // Ambisonics
            if (data->transportType == MASA_TRANSPORT_STEREO)
            {
                meta->transportDefinition = 0x3; // Bit value 011, Cardioid
                meta->channelAngle = 0x1;        // Bit value 001, 90°
                // meta->channelDistance         // Transport channels are assumed to be coincident, therefore no 'Channel distance' field is specified
            }
            break;
        case ANALYZER_HOA2:
            meta->sourceFormat = 0x3; // Ambisonics
            if (data->transportType == MASA_TRANSPORT_STEREO)
            {
                meta->transportDefinition = 0x3; // Bit value 011, Cardioid
                meta->channelAngle = 0x1;        // Bit value 001, 90°
                // meta->channelDistance         // Transport channels are assumed to be coincident, therefore no 'Channel distance' field is specified
            }
            break;
        default:
            assert("Unsupported case for descriptive metadata.");
    }
}


static void formSpatialMetadata(AnalyzerData* data)
{
    MasaDirectionSpatialMetadata* dirMeta = data->outputMetadata.directionalSpatialMeta;
    MasaCommonSpatialMetadata* comMeta = &(data->outputMetadata.commonSpatialMeta);

    for (size_t sf = 0; sf < SUBFRAME_COUNT; sf++)
    {
        for (size_t b = 0; b < MASA_BAND_COUNT; b++)
        {
            for (size_t dir = 0; dir < data->directionCount; dir++)
            {
                indexDirection(&(dirMeta[dir].directionIndex[sf][b]), data->elevation[dir][sf][b], data->azimuth[dir][sf][b], &(data->sphericalGridData));
                dirMeta[dir].directToTotalRatio[sf][b] = (uint8_t) floorf(data->directToTotalRatio[dir][sf][b] * UINT8_MAX);
                dirMeta[dir].spreadCoherence[sf][b] = (uint8_t) floorf(data->spreadCoherence[dir][sf][b] * UINT8_MAX);
            }

            comMeta->diffuseToTotalRatio[sf][b] = (uint8_t) ceilf(data->diffuseToTotalRatio[sf][b] * UINT8_MAX);
            comMeta->surroundCoherence[sf][b] = (uint8_t) floorf(data->surroundCoherence[sf][b] * UINT8_MAX);
            comMeta->remainderToTotalRatio[sf][b] = 0; // No remainder analysis

            // Correct ratios if we have mismatch from one. This is to ensure that ratios sum to UINT8_MAX (255) exactly.
            // This requirement by the format that ratios should sum to total of one.
            if (dirMeta[0].directToTotalRatio[sf][b] + dirMeta[1].directToTotalRatio[sf][b] + comMeta->diffuseToTotalRatio[sf][b] != UINT8_MAX)
            {
                comMeta->diffuseToTotalRatio[sf][b] = UINT8_MAX;
                for (size_t dir = 0; dir < data->directionCount; dir++)
                {
                    comMeta->diffuseToTotalRatio[sf][b] -= dirMeta[dir].directToTotalRatio[sf][b];
                }
            }
        }
    }
}

static void writeTransportFrame(FILE* file, const int16_t* signal, const size_t signalLength)
{
    if (file == NULL)
    {
        fprintf(stderr, "Transport file not open.\n");
        exit(EXIT_FAILURE);
    }

    if (fwrite(signal, sizeof(int16_t), signalLength, file) < signalLength)
    {
        fprintf(stderr, "Transport file write failure.\n");
        exit(EXIT_FAILURE);
    }
}


static void writeMetadataFrame(FILE* file, MasaMetadata* meta)
{
    if (file == NULL)
    {
        fprintf(stderr, "Metadata file not open.\n");
        exit(EXIT_FAILURE);
    }

    if (fwrite(&(meta->descriptiveMeta.formatDescriptor), sizeof(uint8_t), 8, file) != 8)
    {
        fprintf(stderr, "Metadata file write failure.\n");
        exit(EXIT_FAILURE);
    }

    // Construct a single uint16_t from other descriptive meta
    uint16_t descMetaTemp = 0;
    descMetaTemp += meta->descriptiveMeta.numberOfDirections << 15u;
    descMetaTemp += meta->descriptiveMeta.numberOfChannels << 14u;
    descMetaTemp += meta->descriptiveMeta.sourceFormat << 12u;
    if (meta->descriptiveMeta.sourceFormat == 0x0 || meta->descriptiveMeta.sourceFormat == 0x1)
    {
        descMetaTemp += meta->descriptiveMeta.transportDefinition << 9u;
        descMetaTemp += meta->descriptiveMeta.channelAngle << 6u;
        descMetaTemp += meta->descriptiveMeta.channelDistance;
    }
    else if (meta->descriptiveMeta.sourceFormat == 0x2)
    {
        descMetaTemp += meta->descriptiveMeta.channelLayout << 9u;
        // 9 LSB remain at zero
    }
    else if (meta->descriptiveMeta.sourceFormat == 0x3)
    {
        descMetaTemp += meta->descriptiveMeta.transportDefinition << 9u;
        descMetaTemp += meta->descriptiveMeta.channelAngle << 6u;
        // 6 LSB remain at zero
    }

    if (write_uint16_t_LE(&(descMetaTemp), 1, file) != 1)
    {
        fprintf(stderr, "Metadata file write failure.\n");
        exit(EXIT_FAILURE);
    }

    for (size_t sf = 0; sf < SUBFRAME_COUNT; sf++)
    {
        // Direction 1
        if (write_uint16_t_LE(meta->directionalSpatialMeta[0].directionIndex[sf], MASA_BAND_COUNT, file) != MASA_BAND_COUNT)
        {
            fprintf(stderr, "Metadata file write failure.\n");
            exit(EXIT_FAILURE);
        }

        if (fwrite(meta->directionalSpatialMeta[0].directToTotalRatio[sf], sizeof(uint8_t), MASA_BAND_COUNT, file) != MASA_BAND_COUNT)
        {
            fprintf(stderr, "Metadata file write failure.\n");
            exit(EXIT_FAILURE);
        }

        if (fwrite(meta->directionalSpatialMeta[0].spreadCoherence[sf], sizeof(uint8_t), MASA_BAND_COUNT, file) != MASA_BAND_COUNT)
        {
            fprintf(stderr, "Metadata file write failure.\n");
            exit(EXIT_FAILURE);
        }

        if (meta->descriptiveMeta.numberOfDirections == 0x1u)
        {
            // Direction 2
            if (write_uint16_t_LE(meta->directionalSpatialMeta[1].directionIndex[sf], MASA_BAND_COUNT, file) != MASA_BAND_COUNT)
            {
                fprintf(stderr, "Metadata file write failure.\n");
                exit(EXIT_FAILURE);
            }

            if (fwrite(meta->directionalSpatialMeta[1].directToTotalRatio[sf], sizeof(uint8_t), MASA_BAND_COUNT, file) != MASA_BAND_COUNT)
            {
                fprintf(stderr, "Metadata file write failure.\n");
                exit(EXIT_FAILURE);
            }

            if (fwrite(meta->directionalSpatialMeta[1].spreadCoherence[sf], sizeof(uint8_t), MASA_BAND_COUNT, file) != MASA_BAND_COUNT)
            {
                fprintf(stderr, "Metadata file write failure.\n");
                exit(EXIT_FAILURE);
            }
        }

        // Common spatial meta
        if (fwrite(meta->commonSpatialMeta.diffuseToTotalRatio[sf], sizeof(uint8_t), MASA_BAND_COUNT, file) != MASA_BAND_COUNT)
        {
            fprintf(stderr, "Metadata file write failure.\n");
            exit(EXIT_FAILURE);
        }

        if (fwrite(meta->commonSpatialMeta.surroundCoherence[sf], sizeof(uint8_t), MASA_BAND_COUNT, file) != MASA_BAND_COUNT)
        {
            fprintf(stderr, "Metadata file write failure.\n");
            exit(EXIT_FAILURE);
        }

        if (fwrite(meta->commonSpatialMeta.remainderToTotalRatio[sf], sizeof(uint8_t), MASA_BAND_COUNT, file) != MASA_BAND_COUNT)
        {
            fprintf(stderr, "Metadata file write failure.\n");
            exit(EXIT_FAILURE);
        }
    }
}

static void clearAnalysisBuffers(AnalyzerData* data)
{
    for (size_t n = 0; n < CLDFB_BAND_COUNT; n++)
    {
        for (size_t k = 0; k < FOA_CHANNEL_COUNT; k++)
        {
            for (size_t l = 0; l < CLDFB_SLOT_COUNT; l++)
            {
                data->foaSignalReal[k][l][n] = 0;
                data->foaSignalImag[k][l][n] = 0;
                if (data->dirMode == ANALYZER_2DIR_ANALYSIS)
                {
                    data->sectorSignal1Real[k][l][n] = 0;
                    data->sectorSignal1Imag[k][l][n] = 0;
                    data->sectorSignal2Real[k][l][n] = 0;
                    data->sectorSignal2Imag[k][l][n] = 0;
                }
            }
        }
    }

    for (size_t n = 0; n < CLDFB_BAND_COUNT; n++)
    {
        for (size_t k = 0; k < HOA2_CHANNEL_COUNT; k++)
        {
            for (size_t l = 0; l < CLDFB_SLOT_COUNT; l++)
            {
                data->hoa2SignalReal[k][l][n] = 0;
                data->hoa2SignalImag[k][l][n] = 0;
            }
        }
    }
}

static void complexMult(float* outRe, float* outIm, const float in1Re, const float in1Im, const float in2Re, const float in2Im)
{
    *outRe = in1Re * in2Re - in1Im * in2Im;
    *outIm = in1Re * in2Im + in2Re * in1Im;
}


static float fclampf(float in, float lowLim, float highLim)
{
    return fmaxf(lowLim, fminf(highLim, in));
}


// Helper function to write unsigned 16-bit values explicitly as little endian into file using fwrite underneath
static size_t write_uint16_t_LE(const uint16_t* data, size_t nitems, FILE* stream)
{
    size_t itemsWritten = 0;
    uint8_t temp[2];
    for (size_t n = 0; n < nitems; n++)
    {
        temp[0] = (uint8_t) (data[n] & 0xFF);
        temp[1] = (uint8_t) (data[n] >> 8);
        if (fwrite(temp, sizeof(uint8_t), 2, stream) == 2)
        {
            itemsWritten++;
        }
    }

    return itemsWritten;
}

static void printUsage()
{
    printf("MASA analyzer\n\n");
    printf("masaAnalyzer transportMode directionMode inputFormat (descMeta) inputFile outputPcmFile outputMetadataFile\n");
    printf("Valid transport modes are -mono and -stereo\n");
    printf("Valid direction modes are -1dir and -2dir\n");
    printf("Valid input formats are -eigen, -eigenext, -foa, -foaext, -hoa2, and -hoa2ext\n");
    printf("For \"ext\" input formats, (descMeta) can optionally be used to give transport description for external mics. This is given as -descmeta transportDefinition channelAngle channelDistance\n");
    printf("Input file is 48 kHz LEI signed 16-bit PCM with interleaved channels.\n");
    printf("For -eigen, this is 32 channels of Eigenmike input.\n");
    printf("For -eigenext, this is 32 channels of Eigenmike input and 1 or 2 channels of transport signal.\n");
    printf("For -foa, this is 4 FOA channels (ACN order SN3D).\n");
    printf("For -foaext, this is 4 FOA channels (ACN order SN3D), and 1 or 2 channels of transport signal.\n");
    printf("For -hoa2, this is 9 HOA2 channels (ACN order SN3D).\n");
    printf("For -hoa2ext, this is 9 HOA2 channels (ACN order SN3D), and 1 or 2 channels of transport signal.\n");
    printf("Check README.md for more information.\n\n");
}
