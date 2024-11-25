# `pyivasmasa_io`
`pyivasmasa_io` is a Python library for reading and writing MASA (Metadata-Assisted Spatial Audio) metadata files.
These files are supported by the 3GPP IVAS (Third Generation Partnership Project Immersive Voice and Audio Services) codec[^1].
The file format is described in Annex A of the IVAS C-reference specification[^2].

An overview of the MASA audio format in general, estimating MASA metadata parameters from audio signals, and a description of the MASA metadata file format is provided in [^3].

If you are using the code in a publication, please cite [^3].

## Dependencies to be installed by the user
- Python >3.8. The development and testing has been done with version 3.11.
- NumPy. The exact version should not be critical. The development and testing has been done with version 1.28.

[`requirements.txt`](requirements.txt) is provided and can be used to install the dependencies, e.g., with `pipenv` by
```shell
pipenv install -r requirements.txt --python=3.11
```
(or whatever Python version should be used).

Or installing the dependencies into an existing virtual environment with `pip` by
```shell
pip install -r requirements.txt
```

## Installing `pyivasmasa_io` as a library dependency
- Installation manually using `pipenv` from a cloned repository:
  ```shell
  git clone https://github.com/Nokia-Bell-Labs/MASA_tools_for_IVAS.git
  pipenv install MASA_tools_for_IVAS/MASA-IO
  ```

- Similarly, installation manually using `pip` from a cloned repository:
  ```shell
  git clone https://github.com/Nokia-Bell-Labs/MASA_tools_for_IVAS.git
  pip install MASA_tools_for_IVAS/MASA-IO
  ```
  
- Installing directly from [GitHub version](https://github.com/Nokia-Bell-Labs/MASA_tools_for_IVAS/) with `pipenv` 
  ```shell
  pipenv install "pyivasmasa_io @ git+https://github.com/Nokia-Bell-Labs/MASA_tools_for_IVAS.git#subdirectory=MASA-IO"
  ```

- Similarly, installing directly using `pip`
  ```shell
  pip install "pyivasmasa_io @ git+https://github.com/Nokia-Bell-Labs/MASA_tools_for_IVAS.git#subdirectory=MASA-IO"
  ```

- Using the [GitHub version](https://github.com/Nokia-Bell-Labs/MASA_tools_for_IVAS) as a dependency.
  In the `Pipfile` file use the dependency
  ```
  [packages]
  pyivasmasa_io = {git = "https://github.com/Nokia-Bell-Labs/MASA_tools_for_IVAS.git", subdirectory = "MASA-IO"}
  ```

## Usage
### Data structures
- `MasaMetadata` represents one MASA spatial metadata frame corresponding to 20 ms of audio. It consists of
  - `MasaDescriptiveMetadata`, which is a header field containing information, e.g., of the source material, number of transport audio signal channels, and number of directional fields in the current frame,
  - `MasaCommonSpatialMetadata`, which contains the non-directional time/frequency-selective metadata of the frame, including
    - `diffuseToTotalRatio`, which is a (4, 24) `numpy.ndarray` of `numpy.float32` in the range 0..1 describing the diffuse-to-total energy ratio,
    - `surroundCoherence`, which is a (4, 24) `numpy.ndarray` of `numpy.float32` in the range 0..1 describing the surround coherence in the frame,
    - `remainderToTotalRatio`, which is a (4, 24) `numpy.ndarray` of `numpy.float32` in the range 0..1 describing the non-direction and non-diffuse energy in the frame.
  - one or two `MasaDirectionSpatialMetadata`, which contain the directional metadata consisting of
    - `directionIndex`, which is a (4, 24) `numpy.ndarray` of `numpy.uint16` of the spherical indices defining the sound direction azimuth and elevation
    - `directToTotalRatio`, which is a (4, 24) `numpy.ndarray` of `numpy.float32` in the range 0..1 describing the direct-to-total energy ratio,
    - `spreadCoherence`, which is a (4, 24) `numpy.ndarray` of `numpy.float32` in the range 0..1 describing the spread coherence of the sound,
    - `azimuth`, which is a (4, 24) `numpy.ndarray` of `numpy.float32` in the range -180..+180 describing the direction of arrival azimuth of the sound in degrees (0 is front-center and positive values are towards left)
    - `elevation`, which is a (4, 24) `numpy.ndarray` of `numpy.float32` in the range -90..+90 describing the direction of arrival elevation of the sound in degrees (0 is horizontal plane and positive values up)

### Reading a MASA metadata file
```Python
import pyivasmasa_io

py_io_lib = pyivasmasa_io.PyIvasMasaIO()
meta_frames = py_io_lib.ivasmasa_read(masa_meta_filename=masa_meta_filename)
```
- `meta_frames` will be a `list` of `MasaMetadata` instances. Each element corresponds to one MASA metadata frame.

### Reading a range of MASA metadata file
```Python
meta_frames = py_io_lib.ivasmasa_read(masa_meta_filename=masa_meta_filename,
                                      start_frame=5,
                                      n_frames=10,
                                      skip_mode='fast')
```
- This reads `10` frames starting from the index `5` (0-based indexing), using `fast` search.
  - `fast` search assumes that each metadata frame contains the same number of directional fields as the first one in the file.
  - `slow` search parses the descriptive metadata in each frame to see the number of directional fields in the frame to be able to skip correctly. This is useful if the file contains mixed 1- and 2-directional metadata frames.

### Writing a MASA metadata file
```Python
py_io_lib.ivasmasa_write(masa_meta_filename=meta_filename, 
                         meta_frames=meta_frames, 
                         convert_from_angles=False)
```
- `meta_frames` is a `list` of `MasaMetadata` instances.
- When `convert_from_angles == True`, the direction information in the `MasaDirectionSpatialMetadata` fields `azimuth` and `elevation` are used to construct the `directionIndex` information.
- When `convert_from_angles == False`, the direction information in the `MasaDirectionSpatialMetadata` field `directionIndex` is used and the contents of the fields `azimuth` and `elevation` are ignored.

### Concatenating metadata values along time axis into matrices
```Python
import numpy as np

dir_idx1_mat = np.concatenate([one_meta_frame.dir1Meta.directionIndex for one_meta_frame in meta_frames], axis=0).transpose()
ratio1_mat = np.concatenate([one_meta_frame.dir1Meta.directToTotalRatio for one_meta_frame in meta_frames], axis=0).transpose()
spr_coh1_mat = np.concatenate([one_meta_frame.dir1Meta.spreadCoherence for one_meta_frame in meta_frames], axis=0).transpose()
azi1_mat = np.concatenate([one_meta_frame.dir1Meta.azimuth for one_meta_frame in meta_frames], axis=0).transpose()
ele1_mat = np.concatenate([one_meta_frame.dir1Meta.elevation for one_meta_frame in meta_frames], axis=0).transpose()
```
- All the resulting matrices have 24 rows (i.e., bands) and `4*len(meta_frames)` (i.e., number of sub-frames) columns.

## License
Copyright (C) 2024 Nokia Technologies Ltd.

Licensed under the ["Metadata-Assisted Spatial Audio Format (MASA) metadata file reading and writing library (MASA-IO) Non-commercial License"](LICENSE.md). 
See [`LICENSE.md`](LICENSE.md) for details.

## References
[^1] 3GPP, *"Codec for Immersive Voice and Audio Services; Detailed algorithmic description incl. RTP payload format and SDP parameter definition"*, 3rd Generation Partnership Project (3GPP)}, Technical Specification (TS) 26.253, 2024, v18.0.0. https://portal.3gpp.org/desktopmodules/Specifications/SpecificationDetails.aspx?specificationId=3319

[^2] 3GPP, *"Codec for Immersive Voice and Audio Services; C code (floating-point)"*, 3rd Generation Partnership Project (3GPP), Technical Specification (TS) 26.258, 2024, v18.1.0. https://portal.3gpp.org/desktopmodules/Specifications/SpecificationDetails.aspx?specificationId=3324

[^3] Paulus, J., Laaksonen, L., Pihlajakuja, T., Laitinen, M.-V., Vilkamo, J., Vasilache, A., *Metadata-assisted spatial audio (MASA) - an overview,* in Proc. of the IEEE International Symposium on the Internet of Sounds, Erlangen, Germany, 2024.


