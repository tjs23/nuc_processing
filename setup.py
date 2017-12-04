#!/usr/bin/env python

from setuptools import setup

setup(
    name="Nuc Processing",
    version="0.0.1",
    description="A program to generate contacts from Hi-C data",
    install_requires=['numpy', 'scipy'],
    packages=['nuc_processing'],
    include_package_data=True,
)
