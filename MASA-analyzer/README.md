# MASA analyzer tool

## Description

This repository provides a source code distribution of the MASA analyzer tool. This tool creates 3GPP IVAS standard 
compliant MASA format signals from specific input format signals. The standard compliant MASA format is defined in
3GPP IVAS specification TS 26.258 Annex A.

If you are using the tool in a publication, please cite [^1].

## License

Copyright (C) 2024 Nokia Technologies Ltd.

Source code and data binaries are licensed under the ["Metadata-Assisted Spatial Audio Format (MASA) Analyzer Tool Non-commercial License Version 1.0"](LICENSE.md)
with the exception of the contents of the `MASA-analyzer/pocketfft` directory which is licensed under the [BSD-3-Clause license](pocketfft/LICENSE.md).


## Supported platforms

There is no specific supported platform but being standard C implementation, this tool should be widely applicable. 
We have tested and used the software at least on the following platforms:

- Windows 10 with Visual Studio 17
- MacOS with Clang 15 (using Unix Makefiles)


## Dependencies

The MASA analyzer tool depends on an FFT-library. The FFT-library provided with this distribution is the PocketFFT
C-version from a specific commit (SHA 81d171a6d5562e3aaa2c73489b70f564c633ff81) that is provided in 
[GitHub](https://github.com/malfet/pocketfft) or in [MPCDF Gitlab](https://gitlab.mpcdf.mpg.de/mtr/pocketfft).

The FFT-library can be changed to something more suitable or faster with relatively minor modifications. It is
expected that the results of the spatial analysis may have minor differences when the underlying FFT changes as
the TF-domain signal values will be slightly different.

## How to use

### Building

Any standard C capable tool can be used to build the project. For convenience, a CMake-project has been defined for 
the project and can be used to generate platform-specific build folders. A few examples run from this folder are:

The default build system (Unix Makefiles unless changed)
```
cmake -S . -B build-folder
```

Visual Studio 17
```
cmake -S . -B build-folder -G "Visual Studio 15 2017"
```

Xcode
```
cmake -S . -B build-folder -G "Xcode"
```

You can see the available build systems with command `cmake --help`.

To build the produced build, the simplest way is with

```
cmake --build build-folder
```

although you may want to add some additional options depending on if building a release build or using parallel processing. 
Please check the CMake instructions for further options like release-builds and parallel processing builds.


### Running

First, ensure that all the `.bin` files are located in the same folder from which you are invoking the tools. This should 
preferably be the same folder where the built tool executables are. This is to ensure that the tools find the required data 
binaries. The necessary data binaries are 

- `eigen_to_foa_cldfb_domain_filters.bin`
- `eigen_to_hoa2_cldfb_domain_filters.bin`
- `sector_filters.bin`

MASA analyzer can generate MASA format with 1 or 2 transport channels and 1 or 2 directions in the metadata. These are set
with the command-line arguments. In addition, the input format needs to be defined (see below for further description of 
input format options in the next section). An example invocation of the MASA analyzer for 2 transport, 1 direction mode 
with Eigenmike-em32 input is the following.

```
masaAnalyzer -stereo -1dir -eigen input.pcm transport.pcm metadata.met
```

The output is a PCM file containing 1 or 2 transport channels (16-bit integer, channels interleaved) and associated metadata
file. The metadata file contains consecutive independent metadata frames that correspond to exactly 20 ms of audio in the
transport channel file. The input and output of the analyzer use 48 kHz sampling rate. Hence, the audio frame corresponding
to the metadata is exactly 960 samples long per channel. As IVAS supports also 16 kHz and 32 kHz sample rates, it is a valid
approach to resample the PCM file with any suitable tool to fit the desired IVAS input sample rate. The metadata does not
need to be modified for this purpose.

The full usage for the MASA analyzer tool is following:

```
masaAnalyzer transportMode directionMode inputFormat (descMeta) inputFile outputPcmFile outputMetadataFile
Valid transport modes are -mono and -stereo
Valid direction modes are -1dir and -2dir
Valid input formats are -eigen, -eigenext, -foa, -foaext, -hoa2, and -hoa2ext
For "ext" input formats, (descMeta) can optionally be used to give transport description for external mics. This is given as 
-descmeta transportDefinition channelAngle channelDistance
Input file is 48 kHz LEI signed 16-bit PCM with interleaved channels.
For -eigen, this is 32 channels of Eigenmike input.
For -eigenext, this is 32 channels of Eigenmike input and 1 or 2 channels of transport signal.
For -foa, this is 4 FOA channels (ACN order SN3D).
For -foaext, this is 4 FOA channels (ACN order SN3D), and 1 or 2 channels of transport signal.
For -hoa2, this is 9 HOA2 channels (ACN order SN3D).
For -hoa2ext, this is 9 HOA2 channels (ACN order SN3D), and 1 or 2 channels of transport signal.
```

### Input options for analyzer

The input for MASA analyzer is provided as PCM signals where each sample is represented by a 16-bit signed integer value.
As multiple channels are present, they are provided in interleaved format such that the PCM bitstream contains following
order

```
|ch 1 samp 1|ch 2 samp 1|...|ch N  samp 1|ch 1 samp 2|ch 2 samp 2|...|ch N  samp 2|...
```

MASA analyzer supports two approaches for providing input. The first approach is to give direct spatial signals for
analysis. Transport signals are generated from the spatial signals. The options for this approach are:

`-eigen`, where the input file has 32 channel [Eigenmike-em32](https://eigenmike.com/eigenmike-em32) microphone signal.
`-foa`, input is FOA (ACN order, SN3D normalization) signal with 4 channels. Only 1dir analysis is supported for this input.
`-hoa2`, input is HOA2 (ACN order, SN3D normalization) signal with 9 channels.
For both -foa and -hoa2 options, omnidirectional 0th-order component is applied as mono transport,
and two cardioids with azimuth angles +90° and -90° (i.e., pointing left and right) are generated for stereo transports.

The second approach is to provide separately the spatial signal for analysis and the transport signals to be included in
the MASA format. The options in this case are:

`-eigenext`, otherwise similar as `-eigen` but transport channels are included as extra channels at the end. Thus,
the file should have 33 (`-mono`) or 34 (`-stereo`) channels in total.
`-foaext`, input is FOA (ACN order, SN3D normalization) with transport channels included in the end. Thus,
the files should have 5 or 6 channels in total. Only 1dir analysis is supported for this input.
`-hoa2ext`, input HOA2 (ACN order, SN3D normalization) with transport channels included in the end. Thus,
the files should have 10 or 11 channels in total.

In all cases where separate transports are provided, they can be constructed in any desired way but a few assumptions
are made. First, the analysis signal and the transport signals need to be in sync. Second, with stereo transport, it
is assumed that the first channel corresponds to the left channel and the second channel corresponds to the right channel.

When any of the "ext" modes are used with stereo transports, it is possible to optionally give the transport channel
describing descriptive metadata parameters as part of the command line. These are given directly after the inputFormat
argument in the form `-descmeta transportDefinition channelAngle channelDistance` as direct bit values corresponding to
the desired metadata. For example, following full command may be used to have FOA + two stereo transports with the
transports being cardioid mics at ±90° angles and with 9.62 cm distance.

```
masaAnalyzer -stereo -1dir -foaext -descmeta 011 001 100000 input.pcm transport.pcm metadata.met
```

For Eigenmike input, these values are automatically set as 000 001 011110 which are respectively: "Unknown/Other", ±90°,
and 0.0823 (Eigenmike has diameter of 8.4 cm). The possible values are presented in TS 26.258 Annex A.

For FOA and HOA2 inputs, the descriptive metadata is set as 11 011 001 0000000, where the first two bits defines the
sourceFormat ("Ambisonics") and the rest of the bits defines the transportDefinition ("Cardioid"), channelAngle (±90°) and
channelDistance (this field is zero padded, as the transport channels are assumed to be coincident)


## How to use with IVAS

The 3GPP IVAS codec supports MASA format as an input format to the IVAS encoder and the IVAS renderer. The IVAS codec is 
algorithmically detailed in the specification TS 26.253 and the floating point reference code is available in specification 
TS 26.258. The specifications can be found directly with search engines or via website

[3GPP specifications](https://www.3gpp.org/specifications-technologies/specifications-by-series)

The floating point reference code has extensive instructions in readme.txt. However, to get started, the following set of 
commands can be used to obtain static binaural rendering from MASA format with 2 transport channels.

Using encoder & decoder
```
IVAS_enc -masa 2 input_file.met 128000 48 input_file.pcm ivas_bitstream.bit
IVAS_dec BINAURAL 48 ivas_bitstream.bit output_file_binaural.wav
```

Using renderer
```
IVAS_rend -of BINAURAL -fs 48 -if MASA2 -i input_file.pcm -im input_file.met -o output_file_binaural.wav
```

## References

[^1] Paulus, J., Laaksonen, L., Pihlajakuja, T., Laitinen, M.-V., Vilkamo, J., Vasilache, A., *Metadata-assisted spatial audio (MASA) - an overview,* in Proc. of the IEEE International Symposium on the Internet of Sounds, Erlangen, Germany, 2024.


## Version history

- 1.0
  - Initial version
