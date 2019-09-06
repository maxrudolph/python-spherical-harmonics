# python-spherical-harmonics
Utilities for working with spherical harmonics in python.

Routines related to spherical harmonic representations of data.
The routines here are designed to work with Thorsten Becker's .ab file format
which is used with the HC mantle circulation code.

The primary utility of this software is to plot fields from spherical harmonic coefficients and to calculate spherical harmonic expansions of fields. The method used to calculate spherical harmonic expansions here is least-squares fitting. It has only been tested for spherical harmonic degrees up to and including l=15. Instabilities may arise at higher degrees. For fields that lack an exact representation as spherical harmonics, data should be defined on an equispaced mesh rather than a regular (e.g. degree-by-degree) grid. 

A (small) test suite is included in ab_pack_tests.py.

Requirements:
- Python 3
- numpy
- scipy
- matplotlib
- cartopy

(optional) Thorsten Becker's sh_syn tool (included with HC):
+ https://github.com/geodynamics/HC
or here:
+ https://github.com/thwbecker/shansyn
