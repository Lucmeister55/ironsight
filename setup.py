#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import setup

# Define package info
name = "ironsight"
description = "Ferroptosis sensitivity prediction suite downstream to Dorado for Oxford Nanopore DNA sequencing datasets"
with open("README.md", "r") as fh:
    long_description = fh.read()

# Collect info in a dictionnary for setup.py
setup(
    name=name,
    description="Differential methylation calling suite for Nanopore methylation calls PycoMeth",
    version="1.0.0",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Lucmeister55/ironsight,",
    author="Luca Visser",
    author_email="Luca.Visser@UGent.be",
    license="GPL",
    python_requires=">=3.10",
    classifiers=[
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3"
    ],
    install_requires=[],
    packages=["ironsight"],
    package_dir={"ironsight": "ironsight"},
    package_data={name: ["templates/*"]},
    entry_points={"console_scripts": ["ironsight=ironsight.__main__:main"]},
)