#  Python library for reading and writing IVASMASA metadata files.
#
#  Copyright (C) 2024 Nokia Technologies Ltd.
#  See LICENSE.md for license.
from setuptools import setup

from __init__ import __version__

setup(
    name = 'pyivasmasa_io',
    version =  __version__,
    py_modules = ['pyivasmasa_io'],
    license = 'Metadata-Assisted Spatial Audio Format (MASA) metadata file reading and writing library (MASA-IO) Non-commercial License',
    description = 'A library for reading and writing MASA (Metadata-assisted spatial audio) metadata files that conform to 3GPP TS 26.258.',
    url = 'https://github.com/Nokia-Bell-Labs/MASA_tools_for_IVAS',
    install_requires = ['numpy'],
)