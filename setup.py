#!/usr/bin/env python

from setuptools import setup
import sys

compat_pkgs = []
if sys.version_info.major < 3:
    compat_pkgs.append("pathlib2")

setup(
    name="Nuc Processing",
    version="0.0.1",
    description="A program to generate contacts from Hi-C data",
    install_requires=['numpy', 'scipy'] + compat_pkgs,
    packages=['nuc_processing'],
    include_package_data=True,
    entry_points={'console_scripts': [
        'nuc_contact_map=nuc_processing.NucContactMap:main',
        'nuc_contact_probability=nuc_processing.NucContactProbability:main',
        'nuc_process=nuc_processing.NucProcess:main',
    ]},
)
