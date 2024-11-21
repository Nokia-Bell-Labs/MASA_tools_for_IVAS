"""
  py_ivasmasa_io.py

  Python library for reading and writing IVASMASA metadata files.

  Copyright (C) 2024 Nokia Technologies Ltd.
  See LICENSE.md for license.
"""

from functools import cache  # spherical index -> azi+ele caching
import io
import math
import struct  # for bytes to numbers
from typing import BinaryIO, Optional, Union

import numpy as np
import numpy.typing as npt


class MasaMetaBase:
    """ A common base class for the MASA metadata classes. """
    @property
    def _fields(self):
        """ Returns a tuple of attribute names for the class.
        Compatibility feature with NamedTuple. """
        return tuple(vars(self).keys())


class MasaDescriptiveMetadata(MasaMetaBase):
    """ MASA Descriptive Metadata structure """
    def __init__(self,
                 numberOfDirections: int,
                 numberOfChannels: int,
                 formatDescriptor: Optional[bytes] = None,
                 sourceFormat: Optional[int] = 0,
                 transportDefinition: Optional[int] = 0,
                 channelAngle: Optional[int] = 0,
                 channelDistance: Optional[int] = 0,
                 channelLayout: Optional[int] = 0):
        """
        Parameters
        ----------
        numberOfDirections: int {0, 1}
            number of spatial directions in the metadata - 1
        numberOfChannels: int {0, 1}
            number of channels in the transport signal - 1
        formatDescriptor : bytearray of 8 bytes, optional
            expected content b'IVASMASA'
        sourceFormat: int, optional
            original content source format:
              0: default / other
              1: microphone grid
              2: channel-based
              3: Ambisonics
        transportDefinition: int, optional
            transport definition for sourceFormat "default" / "microphone grid" / "Ambisonics"
            0: unknown / other, default
            1: omni
            2: subcardioid
            3: cardioid
            4: supercardioid
            5: hypercardioid
            6: dipole
            7: binaural
        channelAngle: int, optional
            channel angles for directive patterns in sourceFormats "default" / "microphone grid" / "Ambisonics"
            0: unspecified, default
            1: +-90 degrees
            2: +-70 degrees, XY-stereo
            3: +-55 degrees, XY-stereo, ORF-stereo
            4: +-45 degrees, NOS-stereo, XY-stereo, Blumlein pair
            5: +-30 degrees
            6: +-0 degrees, AB-stereo, needs spacing information
            7: reserved
        channelDistance: int, optional
            channel distance for sourceFormats "default" / "microphone grid"
            0: unspecified, default
            1: 0 m / coincident
            2: <0.01 m
            3-62: channel distance: d = (100**(1/59))**(channelDistance-3) / 100
            63: >1 m
        channelLayout: int, optional
            channel layout for sourceFormat "channel-based"
            0: unknown / other, default
            1: other planar
            2: 2.0, CICP2
            3: 5.1, CICP6
            4: 5.1+2, CICP14
            5: 5.1+4, CICP16
            6: 7.1, CICP12
            7: 7.1+4, CICP19
        """
        if formatDescriptor is None:
            self.formatDescriptor = PyIvasMasaIO.DEFAULT_FORMAT_DESCRIPTOR
        else:
            self.formatDescriptor = formatDescriptor

        if numberOfDirections in (0, 1):
            self.numberOfDirections = numberOfDirections
        else:
            raise AttributeError('numberOfDirections must be 0 or 1 (0-based value)')

        if numberOfChannels in (0, 1):
            self.numberOfChannels = numberOfChannels
        else:
            raise AttributeError('numberOfChannels must be 0 or 1 (0-based value)')

        if sourceFormat is not None:
            if not (0 <= sourceFormat <= 3):
                raise AttributeError('sourceFormat must be in range 0..3')
            else:
                self.sourceFormat = sourceFormat
        else:
            self.sourceFormat = 0

        if transportDefinition is not None:
            if not (0 <= transportDefinition <= 7):
                raise AttributeError('transportDefinition must be in range 0..7')
            else:
                self.transportDefinition = transportDefinition
        else:
            self.transportDefinition = 0

        if channelAngle is not None:
            if not (0 <= channelAngle <= 6):
                raise AttributeError('channelAngle must be in range 0..6')
            else:
                self.channelAngle = channelAngle
        else:
            self.channelAngle = 0

        if channelDistance is not None:
            if not (0 <= channelDistance <= 63):
                raise AttributeError('channelDistance must be in range 0..63')
            else:
                self.channelDistance = channelDistance
        else:
            self.channelDistance = 0

        if channelLayout is not None:
            if not (0 <= channelLayout <= 7):
                raise AttributeError('channelLayout must be in range 0..7')
            else:
                self.channelLayout = channelLayout
        else:
            self.channelLayout = 0


class MasaCommonSpatialMetadata(MasaMetaBase):
    """ MASA Common Spatial Metadata structure.
    Each field is an np.array(dtype=np.float32, shape=(SUBFRAME_COUNT, MASA_BAND_COUNT)
    """
    def __init__(self,
                 diffuseToTotalRatio: Optional[np.ndarray] = None,
                 surroundCoherence: Optional[np.ndarray] = None,
                 remainderToTotalRatio: Optional[np.ndarray] = None):
        """
        Parameters
        ----------
        diffuseToTotalRatio: numpy.ndarray of shape (4, 24) of np.float32, optional
            diffuse-to-total energy ratio values in the range 0..1 for each sub-frame and band
        surroundCoherence: numpy.ndarray of shape (4, 24) of np.float32, optional
            surround coherence values in the range 0..1 for each sub-frame and band
        remainderToTotalRatio: numpy.ndarray of shape (4, 24) of np.float32, optional
            remainder-to-total energy ratio values in the range 0..1 for each sub-frame and band

        If the constructor arguments are None, fully-diffuse parameters are created.
        """
        if diffuseToTotalRatio is None:
            self.diffuseToTotalRatio = np.ones(shape=(PyIvasMasaIO.SUBFRAME_COUNT,
                                                      PyIvasMasaIO.MASA_BAND_COUNT),
                                               dtype=np.float32)
        else:
            self.diffuseToTotalRatio = diffuseToTotalRatio

        if surroundCoherence is None:
            self.surroundCoherence = np.zeros(shape=(PyIvasMasaIO.SUBFRAME_COUNT,
                                                     PyIvasMasaIO.MASA_BAND_COUNT),
                                              dtype=np.float32)
        else:
            self.surroundCoherence = surroundCoherence

        if remainderToTotalRatio is None:
            self.remainderToTotalRatio = np.zeros(shape=(PyIvasMasaIO.SUBFRAME_COUNT,
                                                         PyIvasMasaIO.MASA_BAND_COUNT),
                                                  dtype=np.float32)
        else:
            self.remainderToTotalRatio = remainderToTotalRatio


class MasaDirectionSpatialMetadata(MasaMetaBase):
    """ MASA Directional Spatial Metadata structure.

    directionIndex is np.array(dtype=np.uint16, shape=(SUBFRAME_COUNT, MASA_BAND_COUNT)
    Other fields are np.array(dtype=np.float32, shape=(SUBFRAME_COUNT, MASA_BAND_COUNT)
    """
    def __init__(self,
                 directionIndex: Optional[np.ndarray] = None,
                 directToTotalRatio: Optional[np.ndarray] = None,
                 spreadCoherence: Optional[np.ndarray] = None,
                 azimuth: Optional[np.ndarray] = None,
                 elevation: Optional[np.ndarray] = None):
        """
        Parameters
        ----------
        directionIndex: numpy.ndarray of shape (4, 24) of np.uint16, optional
            direction indices for each sub-frame and band.
            if None, the indices are computed from azimuth and elevation. NB: this will have performance impact
        directToTotalRatio: numpy.ndarray of shape (4, 24) of np.float32, optional
           direct-to-total energy ratio values in the range 0..1 for each sub-frame and band
        spreadCoherence: numpy.ndarray of shape (4, 24) of np.float32, optional
            spread coherence values in the range 0..1 for each sub-frame and band
        azimuth: numpy.ndarray of shape (4, 24) of np.float32, optional
            azimuth angle values in the range -180..+180 for each sub-frame and band.
            0 is the center, positive values to left, negative to right
        elevation: numpy.ndarray of shape (4, 24) of np.float32, optional
            elevation angle values in the range -90..90 for each sub-frame and band.
            0 is the horizontal plane, positive values up, negative values down
        """
        if directToTotalRatio is None:
            self.directToTotalRatio = np.zeros(shape=(PyIvasMasaIO.SUBFRAME_COUNT,
                                                      PyIvasMasaIO.MASA_BAND_COUNT),
                                               dtype=np.float32)
        else:
            self.directToTotalRatio = directToTotalRatio

        if spreadCoherence is None:
            self.spreadCoherence = np.zeros(shape=(PyIvasMasaIO.SUBFRAME_COUNT,
                                                   PyIvasMasaIO.MASA_BAND_COUNT),
                                            dtype=np.float32)
        else:
            self.spreadCoherence = spreadCoherence

        if azimuth is None:
            self.azimuth = np.zeros(shape=(PyIvasMasaIO.SUBFRAME_COUNT,
                                           PyIvasMasaIO.MASA_BAND_COUNT),
                                    dtype=np.float32)
        else:
            self.azimuth = azimuth

        if elevation is None:
            self.elevation = np.zeros(shape=(PyIvasMasaIO.SUBFRAME_COUNT,
                                             PyIvasMasaIO.MASA_BAND_COUNT),
                                      dtype=np.float32)
        else:
            self.elevation = elevation

        if directionIndex is None:
            # no indices given, so compute them from azimuth and elevation
            self.directionIndex = np.zeros(shape=(PyIvasMasaIO.SUBFRAME_COUNT,
                                                  PyIvasMasaIO.MASA_BAND_COUNT),
                                           dtype=np.uint16)
            sph_idx = SphericalGrid()
            for sf_idx in range(PyIvasMasaIO.SUBFRAME_COUNT):
                for band_idx in range(PyIvasMasaIO.MASA_BAND_COUNT):
                    self.directionIndex[sf_idx, band_idx] = sph_idx.index_direction(phi=self.azimuth[sf_idx, band_idx],
                                                                                    theta=self.elevation[sf_idx, band_idx])
        else:
            self.directionIndex = directionIndex


class MasaMetadata(MasaMetaBase):
    """ Full MASA Metadata frame. """
    def __init__(self,
                 descriptiveMeta: MasaDescriptiveMetadata,
                 commonSpatialMeta: Optional[MasaCommonSpatialMetadata],
                 dir1Meta: MasaDirectionSpatialMetadata,
                 dir2Meta: Optional[MasaDirectionSpatialMetadata] = None):
        """
        Parameters
        ----------
        descriptiveMeta: MasaDescriptiveMetadata
        commonSpatialMeta: MasaCommonSpatialMetadata, optional
            if None, default common metadata is generated with
            surround coherence = 0, diffuse-to-total energy ratio = 1 - sum(direct-to-total energy ratio)
        dir1Meta: MasaDirectionSpatialMetadata
        dir2Meta: MasaDirectionSpatialMetadata, optional
            if None, create default directional metadata with 0 direct-to-total energy ratio for the second direction.

        NB: The consistency of the information is not validated. E.g., the number of directions in the descriptive
        metadata is not compared against the given directional metadata.
        """
        self.descriptiveMeta = descriptiveMeta
        self.dir1Meta = dir1Meta

        if dir2Meta is not None:
            self.dir2Meta = dir2Meta
        else:
            self.dir2Meta = MasaDirectionSpatialMetadata()

        if commonSpatialMeta is None:
            # create default common spatial metadata with diffuse-to-total energy ratio computed from data
            diffuse_ratio = 1.0 - self.dir1Meta.directToTotalRatio - self.dir2Meta.directToTotalRatio
            self.commonSpatialMeta = MasaCommonSpatialMetadata(diffuseToTotalRatio=diffuse_ratio,
                                                               surroundCoherence=None,
                                                               remainderToTotalRatio=None)

        else:
            self.commonSpatialMeta = commonSpatialMeta


class SphericalGrid:
    """ Spherical grid implementation for mapping between azimuth and elevation angles and spherical grid index.

    The explicit use of single-precision datatype is for closer conformance with the reference implementation.
    """
    # constants
    SPH_CB_SIZE_16BIT = 122

    MASA_NO_INDEX = 32767
    NO_POINTS_EQUATOR = 430
    ANGLE_AT_EQUATOR = np.float32(0.012894427382667)
    ANGLE_AT_EQUATOR_DEG = np.float32(0.738796268264740)
    INV_ANGLE_AT_EQUATOR_DEG = np.float32(1.353553128183453)

    MASA_PI = np.float32(3.14159265358979323846264338327950288)

    BA = np.array((1.698533296919981e+02, 1.223976966547744e+02), dtype=np.float32)
    DELTA = np.array((8.667598146620304e+05, 1.364579281970680e+06), dtype=np.float32)
    DIV = np.array((-0.181726144905283, -0.096108553968620), dtype=np.float32)
    A4 = np.array((-11.005571053314377, -20.809802222734728), dtype=np.float32)

    def __init__(self) -> None:
        """ Constructor """
        n_plus = np.array((18, 29, 45, 61, 71, 85, 88, 92, 95, 107), dtype=np.int16)
        n_minus = np.array((15, 19, 62, 65, 83, 89, 91, 98, 106, 109), dtype=np.int16)

        self.no_phi = np.zeros((SphericalGrid.SPH_CB_SIZE_16BIT, ), dtype=np.int32)  # int32_t[SPH_CB_SIZE_16BIT]

        n1 = np.zeros((SphericalGrid.SPH_CB_SIZE_16BIT, ), dtype=np.int16)
        n1[0] = SphericalGrid.NO_POINTS_EQUATOR
        n1[1:-1] = self.fround(np.float32(422.546)
                               * np.cos(SphericalGrid.ANGLE_AT_EQUATOR
                                        * np.arange(start=1,
                                                    stop=SphericalGrid.SPH_CB_SIZE_16BIT-1,
                                                    dtype=np.float32
                                                    )
                                        )
                               ).astype(np.int16)

        n1[n_plus] += 1
        n1[n_minus] -= 1

        self.no_phi[0] = n1[0]
        for i in range(1, SphericalGrid.SPH_CB_SIZE_16BIT-1):
            self.no_phi[i] = self.no_phi[i-1] + 2 * n1[i]

        self.no_phi[-1] = self.no_phi[-2] + 2
        self.no_theta = SphericalGrid.SPH_CB_SIZE_16BIT

    @cache
    def deindex_direction(self, sph_index: np.uint16) -> tuple[np.float32, np.float32]:
        """ Mapping from spherical index to azimuth and elevation.

        :param sph_index: spherical index to be decoded
        :type sph_index: np.uint16
        :return azimuth: in degrees
        :rtype azimuth: np.float32
        :return elevation: in degrees
        :rtype elevation: np.float32
        """
        if not isinstance(sph_index, np.uint16):
            sph_index = np.uint16(sph_index)

        if sph_index >= 57300:
            ba_crt = SphericalGrid.BA[1]
            div_crt = SphericalGrid.DIV[1]
            a4_crt = SphericalGrid.A4[1]
            del_crt = SphericalGrid.DELTA[1]
        else:
            ba_crt = SphericalGrid.BA[0]
            div_crt = SphericalGrid.DIV[0]
            a4_crt = SphericalGrid.A4[0]
            del_crt = SphericalGrid.DELTA[0]

        estim = ba_crt + div_crt * np.sqrt(del_crt + a4_crt * sph_index)
        if estim > SphericalGrid.SPH_CB_SIZE_16BIT - 1:
            estim = np.float32(SphericalGrid.SPH_CB_SIZE_16BIT - 1)

        id_th = np.int16(self.roundf(estim))
        if id_th < 0:
            id_th = np.int16(0)

        if id_th == 0:
            base_low = np.int32(0)
            base_up = self.no_phi[0]
            n_crt = base_up
        else:
            base_low = self.no_phi[id_th - 1]
            base_up = self.no_phi[id_th]
            n_crt = np.int16(np.right_shift(base_up - base_low, 1))

        sign_theta = np.float32(1.0)

        if sph_index < base_low:
            id_th -= 1

            if id_th == 0:
                base_low = np.int32(0)
                base_up = self.no_phi[0]
                n_crt = base_up
            else:
                base_up = base_low
                base_low = self.no_phi[id_th - 1]
                n_crt = np.int16(np.right_shift(base_up - base_low, 1))

        elif sph_index >= base_up:
            id_th += 1
            base_low = base_up
            base_up = self.no_phi[id_th]
            n_crt = np.int16(np.right_shift(base_up - base_low, 1))

        id_phi = sph_index - base_low
        if sph_index - base_low >= n_crt:
            id_phi -= n_crt
            sign_theta = np.float32(-1.0)

        if id_th == 0:
            theta = np.float32(0.0)
            phi = np.float32(sph_index) * np.float32(360.0) / np.float32(n_crt) - np.float32(180.0)
        else:
            if id_th == self.no_theta - 1:
                phi = np.float32(-180.0)
                theta = np.float32(90.0 * sign_theta)
            else:
                theta = np.float32(id_th) * SphericalGrid.ANGLE_AT_EQUATOR_DEG * sign_theta
                if id_th % 2 == 0:
                    phi = np.float32(np.float32(id_phi) * 360.0 / np.float32(n_crt) - 180.0)
                else:
                    phi = np.float32((np.float32(id_phi) + 0.5) * 360.0 / np.float32(n_crt) - 180.0)

        return phi, theta

    def index_direction(self, phi: np.float32, theta: np.float32) -> np.uint16:
        """ Encode azimuth and elevation into spherical index

        :param phi: input azimuth to be indexed, in degrees
        :type phi: np.float32
        :param theta: input elevation to be indexed, in degrees
        :type theta: np.float32
        :return idx_sph: spherical index
        :rtype idx_sph: np.uint16
        """
        if not isinstance(phi, np.float32):
            phi = np.float32(phi)
        if not isinstance(theta, np.float32):
            theta = np.float32(theta)

        phi = np.float32(phi + 180.0)
        if theta < 0:
            abs_theta = -theta
            sign_th = -1
        else:
            abs_theta = theta
            sign_th = 1

        id_phi, id_th = self._quantize_theta_phi(abs_theta, phi)

        #/* Starting from Equator, alternating positive and negative */
        if id_th > 0:
            cum_n_pos = self.no_phi[id_th - 1]
            cum_n_neg = cum_n_pos + np.right_shift(self.no_phi[id_th] - self.no_phi[id_th - 1], 1, dtype=np.int32)
        else:
            cum_n_pos = 0
            cum_n_neg = 0

        if id_th == 0:
            idx_sph = np.uint16(id_phi)
        else:
            if id_th == self.no_theta - 1:
                id_phi = 0

            if sign_th > 0:
                idx_sph = np.uint16(cum_n_pos + id_phi)
            else:
                idx_sph = np.uint16(cum_n_neg + id_phi)

        return idx_sph

    def _quantize_theta_phi(self, abs_theta: np.float32, phi: np.float32) -> tuple[np.int16, np.int16]:
        """
        :param abs_theta: absolute value of elevation to be quantized, in degrees
        :type abs_theta: np.float32
        :param phi: input azimuth value to be quantized, in degrees
        :type: np.float32

        :return id_phi: azimuth index
        :rtype id_phi: np.int16
        :return id_theta: elevation index
        :rtype id_theta: np.int16
        """
        no_th = self.no_theta
        no_phi_loc = self.no_phi

        if not isinstance(abs_theta, np.float32):
            abs_theta = np.float32(abs_theta)
        if not isinstance(phi, np.float32):
            phi = np.float32(phi)

        id_th, theta_hat = SphericalGrid._quantize_theta(abs_theta, no_th)

        if no_th > 1:
            if id_th == 0:
                id_th1 = 1
            elif id_th == no_th - 1:
                id_th1 = no_th - 2
            else:
                id_th1 = id_th - 1
                id_th2 = id_th + 1

                diff1 = np.float32(abs_theta - id_th1 * SphericalGrid.ANGLE_AT_EQUATOR_DEG)
                diff2 = np.float32(abs_theta - id_th2 * SphericalGrid.ANGLE_AT_EQUATOR_DEG)

                if id_th2 == no_th - 1:
                    diff2 = np.float32(abs_theta - 90.0)

                if np.greater(np.abs(diff1), np.abs(diff2)):
                    id_th1 = id_th2

            theta_hat1 = np.float32(SphericalGrid.ANGLE_AT_EQUATOR_DEG * id_th1)

            if id_th1 == no_th - 1:
                theta_hat1 = np.float32(90.0)

            if no_phi_loc[id_th] > 1:
                if id_th > 0:
                    no_phi = np.uint16(np.right_shift(no_phi_loc[id_th] - no_phi_loc[id_th - 1], 1))
                else:
                    no_phi = np.uint16(no_phi_loc[id_th])

                id_ph, phi_hat = SphericalGrid._quantize_phi(phi, np.int16(id_th % 2 == 1), no_phi)

            else:
                id_ph = SphericalGrid.MASA_NO_INDEX
                phi_hat = np.float32(180.0)

            if (no_phi_loc[id_th1] > 1) and (id_ph < SphericalGrid.MASA_NO_INDEX):
                if id_th1 > 0:
                    no_phi = np.uint16(np.right_shift(no_phi_loc[id_th1] - no_phi_loc[id_th1 - 1], 1))
                else:
                    no_phi = np.uint16(no_phi_loc[id_th1])

                id_ph1, phi_hat1 = SphericalGrid._quantize_phi(phi, np.int16(id_th1 % 2 == 1), no_phi)

                # NB: this is not always behaving identically with the C-version
                d = SphericalGrid._direction_distance(abs_theta, theta_hat, phi, phi_hat)
                d1 = SphericalGrid._direction_distance(abs_theta, theta_hat1, phi, phi_hat1)

                if d1 > d:
                    id_ph = id_ph1
                    id_th = id_th1  # /* theta_hat = ANGLE_AT_EQUATOR_DEG * id_th1; */

        else:
            if id_th > 0:
                no_phi = np.uint16(np.right_shift(no_phi_loc[id_th] - no_phi_loc[id_th - 1], 1))
            else:
                no_phi = np.uint16(no_phi_loc[id_th])

            id_ph, phi_hat = SphericalGrid._quantize_phi(phi, (id_th % 2 == 1), no_phi)

        id_phi = id_ph
        id_theta = id_th

        return id_phi, id_theta

    @staticmethod
    def _quantize_theta(x: np.float32, no_cb: int) -> tuple[np.int16, np.float32]:
        """
        :param x: theta value to be quantized
        :type x: np.float32
        :param no_cb: number of codewords
        :type no_cb: np.int16
        :return idx: index of quantized value
        :rtype idx: np.int16
        :return xhat: quantized value
        :rtype xhat: np.float32
        """
        if not isinstance(x, np.float32):
            x = np.float32(x)

        imin = np.int16(x * SphericalGrid.INV_ANGLE_AT_EQUATOR_DEG + 0.5)

        if imin >= no_cb - 1:
            imin = np.int16(no_cb - 1)
            diff1 = np.float32(x - 90.0)
            diff2 = np.float32(x - SphericalGrid.ANGLE_AT_EQUATOR_DEG * (imin - 1))

            if np.greater(np.abs(diff1), np.abs(diff2)):
                imin -= 1
                xhat = np.float32(imin * SphericalGrid.ANGLE_AT_EQUATOR_DEG)

            else:
                xhat = np.float32(90.0)

        else:
            xhat = np.float32(imin * SphericalGrid.ANGLE_AT_EQUATOR_DEG)

        return imin, xhat

    @staticmethod
    def _quantize_phi(phi: np.float32, flag_delta: bool, n: np.int16) -> tuple[np.int16, np.float32]:
        """
        :param phi: azimuth value to be quantized
        :type phi: np.float32
        :param flag_delta: flag indicating if the azimuth codebook is translated or not
        :type flag_delta: bool
        :param n: azimuth codebook size
        :type n: np.int16
        :return id_phi: index azimuth
        :rtype id_phi: np.int16
        :return phi_hat: quantized azimuth
        :rtype phi_hat: np.float32
        """
        if not isinstance(phi, np.float32):
            phi = np.float32(phi)

        delta_phi = np.float32(360.0 / n)
        two = np.float32(2.0)

        if n == 1:
            phi_hat = np.float32(0.0)
            id_phi = np.int16(0)

        else:
            if flag_delta:
                dd = np.float32(delta_phi / two)
            else:
                dd = np.float32(0.0)

            id_phi = np.int16((phi - dd + delta_phi / two) / delta_phi)

            if id_phi == n:
                id_phi = 0

            if id_phi == -1:
                id_phi = n - 1

            phi_hat = np.float32(id_phi * delta_phi + dd)

        return id_phi, phi_hat

    @staticmethod
    def _direction_distance(theta: np.float32, theta_hat: np.float32, phi: np.float32, phi_hat: np.float32) -> np.float32:
        """ Distance between two directions given in azimuth and elevation.

        :param theta: elevation absolute value, in radians
        :type theta: np.float32
        :param theta_hat: quantized elevation value in absolute value, in radians
        :type theta_hat: np.float32
        :param phi: azimuth value, in radians
        :type phi: np.float32
        :param phi_hat: quantized azimuth value, in radians
        :type phi_hat: np.float32
        :return d: distortion value
        :rtype d: np.float32

        Note that this does not return exactly the same value as the C-implementation due to numerical precision.
        """
        if not isinstance(theta, np.float32):
            theta = np.float32(theta)
        if not isinstance(theta_hat, np.float32):
            theta_hat = np.float32(theta_hat)

        theta *= SphericalGrid.MASA_PI / 180
        theta_hat *= SphericalGrid.MASA_PI / 180

        d = (np.sin(theta, dtype=np.float32) * np.sin(theta_hat, dtype=np.float32)
             + np.cos(theta, dtype=np.float32) * np.cos(theta_hat, dtype=np.float32)
             * np.cos((phi - phi_hat) * SphericalGrid.MASA_PI / 180, dtype=np.float32))

        return d

    @staticmethod
    def fround(x: npt.NDArray[float]) -> npt.NDArray[float]:
        """ Rounding similar to C's "round away from zero". """
        a = np.abs(x)
        b = np.floor(a) + np.floor(2 * (a % 1))
        if isinstance(x, np.ndarray):
            return np.where(x > 0,
                            b,
                            -b)
        else:
            if x > 0:
                return b
            else:
                return -b

    @staticmethod
    def roundf(x: float) -> int:
        """ Rounding similar to C's "round away from zero". """
        a = abs(x)
        b = math.floor(a) + math.floor(2 * (a % 1))

        if x > 0:
            return b
        else:
            return -b


class PyIvasMasaIO:
    """  A class for reading IVAS-MASA metadata files. """
    DEFAULT_FORMAT_DESCRIPTOR = b'IVASMASA'  # (0x49, 0x56, 0x41, 0x53, 0x4D, 0x41, 0x53, 0x41). NB: other implementations use uint8s
    MASA_BAND_COUNT = 24
    SUBFRAME_COUNT = 4
    FRAME_LEN_SAMPLES = 960

    # sizes of 1dir and 2dir frames
    # header + desc + n_sf*(n_bands*(dir1 + dir2 + common))
    N_BYTES_SHORT = 682  # 8 + 2 + 4*(24*(2+1+1 + 0 + 1+1+1))
    N_BYTES_LONG = 1066  # 8 + 2 + 4*(24*(2+1+1 + 2+1+1 + 1+1+1))

    # data structure for one sub-frame of metadata in the cases of 1-directional metadata
    onedir_spatial_meta_dtype = np.dtype([('dir1_idx', '<u2', MASA_BAND_COUNT),  # uint16, LE *24
                                          ('dir1_ratios', '<u1', MASA_BAND_COUNT),  # uint8, LE *24
                                          ('dir1_spr_coh', '<u1', MASA_BAND_COUNT),  # uint8, LE *24
                                          ('diff_ratio', '<u1', MASA_BAND_COUNT),  # uint8, LE *24
                                          ('sur_coh', '<u1', MASA_BAND_COUNT),  # uint8, LE *24
                                          ('rem_ratio', '<u1', MASA_BAND_COUNT),  # uint8, LE *24
                                          ])

    # data structure for one sub-frame of metadata in the cases of 2-directional metadata
    twodir_spatial_meta_dtype = np.dtype([('dir1_idx', '<u2', MASA_BAND_COUNT),  # uint16, LE *24
                                          ('dir1_ratios', '<u1', MASA_BAND_COUNT),  # uint8, LE *24
                                          ('dir1_spr_coh', '<u1', MASA_BAND_COUNT),  # uint8, LE *24
                                          ('dir2_idx', '<u2', MASA_BAND_COUNT),  # uint16, LE *24
                                          ('dir2_ratios', '<u1', MASA_BAND_COUNT),  # uint8, LE *24
                                          ('dir2_spr_coh', '<u1', MASA_BAND_COUNT),  # uint8, LE *24
                                          ('diff_ratio', '<u1', MASA_BAND_COUNT),  # uint8, LE *24
                                          ('sur_coh', '<u1', MASA_BAND_COUNT),  # uint8, LE *24
                                          ('rem_ratio', '<u1', MASA_BAND_COUNT),  # uint8, LE *24
                                          ])

    def __init__(self):
        """ Constructor """
        super().__init__()

        self.sph_grid = SphericalGrid()
        self.uint8_max = np.float32(np.iinfo(np.uint8).max)

    def ivasmasa_write(self,
                       masa_meta_filename: str,
                       meta_frames: list[MasaMetadata],
                       convert_from_angles: bool = False) -> None:
        """ Write MASA metadata into a file.

        :param masa_meta_filename: path to the MASA metadata file
        :type masa_meta_filename: str
        :param meta_frames: MASA metadata frames
        :type meta_frames: list[MasaMetadata]
        :param convert_from_angles: compute direction indices from azimuth and elevation (False->use given dir_idx)
        :type convert_from_angles: bool
        """
        # construct full metadata in byte buffer
        byte_buff = b''
        for one_frame in meta_frames:
            byte_buff += PyIvasMasaIO._get_format_descriptor()
            byte_buff += PyIvasMasaIO._form_desc_meta(one_frame.descriptiveMeta)
            if one_frame.descriptiveMeta.numberOfDirections == 0:
                byte_buff += self._form_masa_meta(dir1_meta=one_frame.dir1Meta,
                                                  dir2_meta=None,
                                                  common_meta=one_frame.commonSpatialMeta,
                                                  convert_from_angles=convert_from_angles)
            else:
                byte_buff += self._form_masa_meta(dir1_meta=one_frame.dir1Meta,
                                                  dir2_meta=one_frame.dir2Meta,
                                                  common_meta=one_frame.commonSpatialMeta,
                                                  convert_from_angles=convert_from_angles)

        # write the buffer into a file
        with open(masa_meta_filename, 'wb') as f:
            f.write(byte_buff)

    @staticmethod
    def _get_format_descriptor() -> bytes:
        """  Return b'IVASMASA' format descriptor. """
        return PyIvasMasaIO.DEFAULT_FORMAT_DESCRIPTOR

    @staticmethod
    def _form_desc_meta(descriptiveMeta: MasaDescriptiveMetadata) -> bytes:
        """ Form 2-byte descriptive metadata from data properties. """
        # construct a single uint16_t from other descriptive meta
        descMetaTemp = 0
        descMetaTemp += descriptiveMeta.numberOfDirections << 15
        descMetaTemp += descriptiveMeta.numberOfChannels << 14
        descMetaTemp += descriptiveMeta.sourceFormat << 12
        if (descriptiveMeta.sourceFormat == 0x0) or (descriptiveMeta.sourceFormat == 0x1):
            descMetaTemp += descriptiveMeta.transportDefinition << 9
            descMetaTemp += descriptiveMeta.channelAngle << 6
            descMetaTemp += descriptiveMeta.channelDistance

        elif descriptiveMeta.sourceFormat == 0x2:
            descMetaTemp += descriptiveMeta.channelLayout << 9
            # 9 LSB remain at zero

        elif descriptiveMeta.sourceFormat == 0x3:
            descMetaTemp += descriptiveMeta.transportDefinition << 9
            descMetaTemp += descriptiveMeta.channelAngle << 6
            # 6 LSB remain at zero

        return struct.pack('<H', descMetaTemp)  # LE uint16_t

    @staticmethod
    def _sanitize_azi_ele(azi: np.ndarray, ele: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """ Wrap azimuth at +/-180 degrees. If elevation exceed +-90, wrap it
        and adjust the matching azimuth.

        :param azi: azimuth angle in degrees
        :type azi: np.ndarray of np.float32
        :param ele: elevation angle in degrees
        :type ele: np.ndarray of np.float32
        :return: tuple[np.ndarray[np.float32], np.ndarray[np.float32]]. sanitized angles
        """
        if not isinstance(azi, np.ndarray):
            azi = np.array([azi])

        if not isinstance(ele, np.ndarray):
            ele = np.array([ele])

        azi = azi.astype(np.float32)
        ele = ele.astype(np.float32)

        if np.any(np.abs(azi) > 180.0) or np.any(np.abs(ele) > 90.0):
            ele = np.mod(ele + 180.0, 360.0) - 180.0  # eliminate extra wraps

            flip_pos = np.nonzero(np.abs(ele) > 90.0)
            ele[flip_pos] = np.sign(ele[flip_pos]) * 180.0 - ele[flip_pos]
            azi[flip_pos] += 180.0  # when ele changed sign, need to flip azi by 180 degrees

            azi = np.mod(azi + 180.0, 360.0) - 180.0

            return azi, ele

        else:
            return azi, ele

    def _form_masa_meta(self,
                        dir1_meta: MasaDirectionSpatialMetadata,
                        dir2_meta: Union[None, MasaDirectionSpatialMetadata],
                        common_meta: MasaCommonSpatialMetadata,
                        convert_from_angles: bool = False) -> bytes:
        """ Form actual MASA metadata (one or two directional metadata, common metadata) into a byte buffer.

        :param dir1_meta: direction 1 metadata
        :type dir1_meta: MasaDirectionSpatialMetadata
        :param dir2_meta: direction 2 metadata
        :type dir2_meta: MasaDirectionSpatialMetadata or None
        :param common_meta: common metadata
        :type common_meta: MasaCommonSpatialMetadata
        :param convert_from_angles: flag if azimuth and elevation should be converted to direction index or the
           provided direction indices used instead
        :type convert_from_angles: bool
        """
        byte_buff = b''

        for sf in range(PyIvasMasaIO.SUBFRAME_COUNT):
            if convert_from_angles:
                azi_data, ele_data = PyIvasMasaIO._sanitize_azi_ele(azi=dir1_meta.azimuth[sf, :],
                                                                    ele=dir1_meta.elevation[sf, :])

                dir_idx_data = np.zeros((PyIvasMasaIO.MASA_BAND_COUNT,), dtype=np.uint16)
                for band_idx in range(PyIvasMasaIO.MASA_BAND_COUNT):
                    dir_idx_data[band_idx] = self.sph_grid.index_direction(np.float32(SphericalGrid.roundf(2.0 * azi_data[band_idx]) / 2.0),
                                                                           np.float32(SphericalGrid.roundf(2.0 * ele_data[band_idx]) / 2.0))

            else:
                dir_idx_data = dir1_meta.directionIndex[sf, :]

            tmp_data = dir_idx_data.astype('<H')  # uint16 LE
            byte_buff += tmp_data.tobytes()

            dir1_ratio = np.ceil(dir1_meta.directToTotalRatio[sf, :] * self.uint8_max).astype(np.uint8)

            byte_buff += dir1_ratio.tobytes()

            tmp_data = np.floor(dir1_meta.spreadCoherence[sf, :] * self.uint8_max).astype(np.uint8)
            byte_buff += tmp_data.tobytes()

            if dir2_meta is not None:
                if convert_from_angles:
                    azi_data, ele_data = PyIvasMasaIO._sanitize_azi_ele(azi=dir2_meta.azimuth[sf, :],
                                                                        ele=dir2_meta.elevation[sf, :])

                    dir_idx_data = np.zeros((PyIvasMasaIO.MASA_BAND_COUNT,), dtype=np.uint16)
                    for band_idx in range(PyIvasMasaIO.MASA_BAND_COUNT):
                        dir_idx_data[band_idx] = self.sph_grid.index_direction(np.float32(SphericalGrid.roundf(2.0 * azi_data[band_idx]) / 2.0),
                                                                               np.float32(SphericalGrid.roundf(2.0 * ele_data[band_idx]) / 2.0))

                else:
                    dir_idx_data = dir2_meta.directionIndex[sf, :]

                tmp_data = dir_idx_data.astype('<H')  # uint16 LE
                byte_buff += tmp_data.tobytes()

                dir2_ratio = np.ceil(dir2_meta.directToTotalRatio[sf, :] * self.uint8_max).astype(np.uint8)
                byte_buff += dir2_ratio.tobytes()

                tmp_data = np.floor(dir2_meta.spreadCoherence[sf, :] * self.uint8_max).astype(np.uint8)
                byte_buff += tmp_data.tobytes()

            else:
                dir2_ratio = np.zeros_like(dir1_ratio)

            diff_ratio = np.floor(common_meta.diffuseToTotalRatio[sf, :] * self.uint8_max).astype(np.uint8)
            tmp_data = np.floor(common_meta.surroundCoherence[sf, :] * self.uint8_max).astype(np.uint8)
            rem_ratio = np.floor(common_meta.remainderToTotalRatio[sf, :] * self.uint8_max).astype(np.uint8)

            # handle direct and diffuse ratios correctly
            to_mod = (dir1_ratio + dir2_ratio + diff_ratio + rem_ratio) != 255
            rem_zero = rem_ratio[to_mod] == 0
            rem_nonzero = np.logical_not(rem_zero)
            rem_ratio[to_mod][rem_zero] = (255
                                           - dir1_ratio[to_mod][rem_zero]
                                           - dir2_ratio[to_mod][rem_zero]
                                           - diff_ratio[to_mod][rem_zero])
            diff_ratio[to_mod][rem_nonzero] = (255
                                               - (dir1_ratio[to_mod][rem_nonzero] + dir2_ratio[to_mod][rem_nonzero]))

            byte_buff += diff_ratio.tobytes()
            byte_buff += tmp_data.tobytes()
            byte_buff += rem_ratio.tobytes()

        return byte_buff

    def ivasmasa_read(self, masa_meta_filename: str,
                      start_frame: int = 0,
                      n_frames: int = -1,
                      skip_mode: str = 'fast') -> Union[None, list[MasaMetadata]]:
        """ Read MASA metadata file or a range from it.

        :param masa_meta_filename: path to the MASA metadata file
        :type masa_meta_filename: str
        :param start_frame: start reading from this frame. default=0
        :type start_frame: int
        :param n_frames: number of frames to read (or try to read). n_frames <= 0 returns all read frames
        :type n_frames: int
        :param skip_mode: when reading start from non-beginning, the seek can be done 'fast' assuming all frames have
                          the same number of directions, or 'slow" reading the header of each frame for the skip
        :type skip_mode: str
        :return: None or a list of MasaMetadata
        """
        all_frames = []
        read_count = 0
        # pre-load the whole file into memory
        with open(masa_meta_filename, 'rb') as f_:
            f = io.BytesIO(f_.read())

        while (n_frames <= 0) or ((n_frames > 0) and (read_count < n_frames)):
            if PyIvasMasaIO._read_format_descriptor(f)[1]:  # try to read 'IVASMASA' header
                desc_ok, desc_meta = PyIvasMasaIO._read_desc_meta(f)  # read descriptive metadata

                if start_frame > 0:
                    # need to perform seek()
                    if skip_mode == 'fast':
                        # read n_dirs from the first frame and assume to be constant
                        n_dirs = desc_meta.numberOfDirections  # 0-based
                        if n_dirs == 0:
                            skip_sz = PyIvasMasaIO.N_BYTES_SHORT
                        else:
                            skip_sz = PyIvasMasaIO.N_BYTES_LONG

                        n_bytes_skip = start_frame * skip_sz
                        f.seek(n_bytes_skip, 1)  # whence=1 => relative seek

                    else:
                        # read n_dirs from each frame and skip accordingly
                        skip_count = 0
                        while skip_count < start_frame:
                            n_dirs = desc_meta.numberOfDirections  # 0-based
                            if n_dirs == 0:
                                skip_sz = PyIvasMasaIO.N_BYTES_SHORT
                            else:
                                skip_sz = PyIvasMasaIO.N_BYTES_LONG

                            n_bytes_skip = skip_sz - 2 - len(PyIvasMasaIO.DEFAULT_FORMAT_DESCRIPTOR)  # skip only the spatial data
                            f.seek(n_bytes_skip, 1)  # whence=1 => relative seek
                            skip_count += 1

                            if PyIvasMasaIO._read_format_descriptor(f)[1]:
                                desc_ok, desc_meta = self._read_desc_meta(f)
                            else:
                                break

                    start_frame = -1  # seek() has been done now

                if desc_ok:
                    dir_ok, dir1_meta, dir2_meta, common_meta = self._read_masa_meta(f, n_dirs=desc_meta.numberOfDirections)

                    if dir_ok:
                        frame_meta = MasaMetadata(descriptiveMeta=desc_meta,
                                                  commonSpatialMeta=common_meta,
                                                  dir1Meta=dir1_meta,
                                                  dir2Meta=dir2_meta)
                        all_frames.append(frame_meta)
                        read_count += 1
                    else:
                        # problems parsing directional metadata
                        break
                else:
                    # problems parsing descriptive metadata
                    break
            else:
                # problems reading format descriptor
                break

        return all_frames

    @staticmethod
    def _read_format_descriptor(f: BinaryIO) -> tuple[bool, bool]:
        """ Try to read b'IVASMASA' from the file.

        :param f: handle to open Masa metadata file
        :type f: BinaryIO
        :return: read_ok: bool if file could be read, decriptor_ok: bool if descriptor matches the expected one
        """
        if this_bytes := f.read(len(PyIvasMasaIO.DEFAULT_FORMAT_DESCRIPTOR)):
            return True, this_bytes == PyIvasMasaIO.DEFAULT_FORMAT_DESCRIPTOR
        else:
            return False, False

    @staticmethod
    def _read_desc_meta(f: BinaryIO) -> tuple[bool, Union[None, MasaDescriptiveMetadata]]:
        """ Read descriptive metadata from file.

        :param f: handle to open Masa metadata file
        :type f: BinaryIO
        :return: read_ok: bool if file could be read, None or MasaDescriptiveMetadata structure
        """
        if this_bytes := f.read(2):
            temp = struct.unpack('<H', this_bytes)[0]  # LE uint16_t

            source_format = (temp >> 12) & 0x3

            if (source_format == 0x0) or (source_format == 0x1):
                transport_definition = (temp >> 9) & 0x7
                channel_angle = (temp >> 6) & 0x7
                channel_distance = temp & 0x3F
                channel_layout = 0  # set to zero as unused

            elif source_format == 0x2:
                channel_layout = (temp >> 9) & 0x7
                transport_definition = 0  # set to zero as unused
                channel_angle = 0  # set to zero as unused
                channel_distance = 0  # set to zero as unused

            else:  # source_format == 0x3:
                transport_definition = (temp >> 9) & 0x7
                channel_angle = (temp >> 6) & 0x7
                channel_distance = 0  # set to zero as unused
                channel_layout = 0  # set to zero as unused

            descriptive_meta = MasaDescriptiveMetadata(formatDescriptor=PyIvasMasaIO.DEFAULT_FORMAT_DESCRIPTOR,
                                                       numberOfDirections=(temp >> 15) & 0x1,
                                                       numberOfChannels=(temp >> 14) & 0x1,
                                                       sourceFormat=source_format,
                                                       transportDefinition=transport_definition,
                                                       channelAngle=channel_angle,
                                                       channelDistance=channel_distance,
                                                       channelLayout=channel_layout)

            return True, descriptive_meta
        else:
            return False, None

    def _read_masa_meta(self, f: BinaryIO, n_dirs: int) \
            -> tuple[bool, Union[None, MasaDirectionSpatialMetadata],
            Union[None, MasaDirectionSpatialMetadata],
            Union[None, MasaCommonSpatialMetadata]]:
        """ Read the actual Masa metadata.

        :param f: open file handle to Masa metadata file
        :type f: BinaryIO
        :param n_dirs: number of directions to read. 0-based
        :type n_dirs: int
        :return: read_ok: bool if file could be read, dir1 metadata, dir2 metadata (with zeros if n_dirs==0), common metadata
        """
        # read all 4 sub-frames in one go
        if n_dirs == 0:  # 0-based
            read_dtype = self.onedir_spatial_meta_dtype
        else:
            read_dtype = self.twodir_spatial_meta_dtype

        # get the bytes from the buffer and check that enough data was obtained
        n_bytes_to_read = read_dtype.itemsize * PyIvasMasaIO.SUBFRAME_COUNT
        tmp_bytes = f.read(n_bytes_to_read)

        if (tmp_bytes is None) or (len(tmp_bytes) != n_bytes_to_read):
            print(f'WARNING: Not enough content available in the buffer!')
            return False, None, None, None

        tmp_data = np.frombuffer(tmp_bytes,
                                 dtype=read_dtype,
                                 count=PyIvasMasaIO.SUBFRAME_COUNT)

        if len(tmp_data) < PyIvasMasaIO.SUBFRAME_COUNT:
            # TODO: add checks within 4 sub-frames
            return False, None, None, None

        # concatenate from struct fields into matrices of size (SUBFRAME_COUNT, MASA_BAND_COUNT)
        mat_dir1_idx = np.stack([one_sf['dir1_idx'] for one_sf in tmp_data], axis=0).astype(np.uint16)
        mat_dir1_ratios = (np.stack([one_sf['dir1_ratios'] for one_sf in tmp_data], axis=0) / self.uint8_max).astype(np.float32)
        mat_dir1_spr_coh = (np.stack([one_sf['dir1_spr_coh'] for one_sf in tmp_data], axis=0) / self.uint8_max).astype(np.float32)

        if n_dirs == 0:  # 0-based
            mat_dir2_idx = np.zeros((PyIvasMasaIO.SUBFRAME_COUNT, PyIvasMasaIO.MASA_BAND_COUNT), dtype=np.uint16)
            mat_dir2_ratios = np.zeros((PyIvasMasaIO.SUBFRAME_COUNT, PyIvasMasaIO.MASA_BAND_COUNT), dtype=np.float32)
            mat_dir2_spr_coh = np.zeros((PyIvasMasaIO.SUBFRAME_COUNT, PyIvasMasaIO.MASA_BAND_COUNT), dtype=np.float32)
        else:
            mat_dir2_idx = np.stack([one_sf['dir2_idx'] for one_sf in tmp_data], axis=0).astype(np.uint16)
            mat_dir2_ratios = (np.stack([one_sf['dir2_ratios'] for one_sf in tmp_data], axis=0) / self.uint8_max).astype(np.float32)
            mat_dir2_spr_coh = (np.stack([one_sf['dir2_spr_coh'] for one_sf in tmp_data], axis=0) / self.uint8_max).astype(np.float32)

        mat_diff_ratio = (np.stack([one_sf['diff_ratio'] for one_sf in tmp_data], axis=0) / self.uint8_max).astype(np.float32)
        mat_sur_coh = (np.stack([one_sf['sur_coh'] for one_sf in tmp_data], axis=0) / self.uint8_max).astype(np.float32)
        mat_rem_ratio = (np.stack([one_sf['rem_ratio'] for one_sf in tmp_data], axis=0) / self.uint8_max).astype(np.float32)

        # implementation using the per-element spherical index decoding. this allows using value memoization
        mat_dir1_azi = np.zeros((PyIvasMasaIO.SUBFRAME_COUNT, PyIvasMasaIO.MASA_BAND_COUNT), dtype=np.float32)
        mat_dir1_ele = np.zeros((PyIvasMasaIO.SUBFRAME_COUNT, PyIvasMasaIO.MASA_BAND_COUNT), dtype=np.float32)
        mat_dir2_azi = np.zeros((PyIvasMasaIO.SUBFRAME_COUNT, PyIvasMasaIO.MASA_BAND_COUNT), dtype=np.float32)
        mat_dir2_ele = np.zeros((PyIvasMasaIO.SUBFRAME_COUNT, PyIvasMasaIO.MASA_BAND_COUNT), dtype=np.float32)

        for sf_idx in range(PyIvasMasaIO.SUBFRAME_COUNT):
            for b_idx in range(PyIvasMasaIO.MASA_BAND_COUNT):
                mat_dir1_azi[sf_idx, b_idx], mat_dir1_ele[sf_idx, b_idx] = self.sph_grid.deindex_direction(mat_dir1_idx[sf_idx, b_idx])

        if n_dirs > 0:  # having this as a separate loop instead of being inside the loop above makes this faster for 1dir
            for sf_idx in range(PyIvasMasaIO.SUBFRAME_COUNT):
                for b_idx in range(PyIvasMasaIO.MASA_BAND_COUNT):
                    mat_dir2_azi[sf_idx, b_idx], mat_dir2_ele[sf_idx, b_idx] = self.sph_grid.deindex_direction(mat_dir2_idx[sf_idx, b_idx])

        mat_dir1_azi = self.sph_grid.fround(2.0 * mat_dir1_azi) / 2.0
        mat_dir1_ele = self.sph_grid.fround(2.0 * mat_dir1_ele) / 2.0

        mat_dir2_azi = self.sph_grid.fround(2.0 * mat_dir2_azi) / 2.0
        mat_dir2_ele = self.sph_grid.fround(2.0 * mat_dir2_ele) / 2.0

        dir1_meta = MasaDirectionSpatialMetadata(directionIndex=mat_dir1_idx,
                                                 directToTotalRatio=mat_dir1_ratios,
                                                 spreadCoherence=mat_dir1_spr_coh,
                                                 azimuth=mat_dir1_azi,
                                                 elevation=mat_dir1_ele)

        dir2_meta = MasaDirectionSpatialMetadata(directionIndex=mat_dir2_idx,
                                                 directToTotalRatio=mat_dir2_ratios,
                                                 spreadCoherence=mat_dir2_spr_coh,
                                                 azimuth=mat_dir2_azi,
                                                 elevation=mat_dir2_ele)

        common_meta = MasaCommonSpatialMetadata(diffuseToTotalRatio=mat_diff_ratio,
                                                surroundCoherence=mat_sur_coh,
                                                remainderToTotalRatio=mat_rem_ratio)

        return True, dir1_meta, dir2_meta, common_meta


