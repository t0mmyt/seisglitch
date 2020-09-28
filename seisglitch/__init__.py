#!/usr/bin/env python

# Copyright 2019 John-Robert Scholz
#
# This file is part of Seisglitch.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# -*- coding: utf-8 -*-



"""
**SEISglitch** is a Python package based on ObsPy_.
Its purpose is to detect "glitches" (also referred to as "long-period disurbances") in the raw seismic data of the VBB (very broadband) 
and SP (short period) seismometers, both part of the SEIS instrument of NASA's Discovery mission "InSight" to planet Mars. 
The package allows to detect glitches, plot them in different ways, and remove them from the seismic raw data. 
**SEISglitch** corresponds to the MPS method detailed in Scholz et al. (2020, see below), all other methods
(ISAE, UCLA, IPGP) were implemented in MATLAB and their essential scripts are also delivered with this software.
**SEISglitch** comes further shipped with some useful features to ease data handling (e.g. download, decimation, time conversion).

| Find the actual Python code at:
| https://pss-gitlab.math.univ-paris-diderot.fr/data-processing-wg/seisglitch

The code is based on the following, peer-reviewed paper: 

    Scholz, J.-R., Widmer-Schnidrig, R., P. Davis, P. Lognonne, B. Pinot, R. F. Garcia, et al. 
    “Detection, Analysis and Removal of Glitches from InSight’s Seismic Data from Mars.” Earth and Space Science, submitted (2020).

If you find ``SEISglitch`` useful, used the MATLAB alternatives, corrected data as to the considerations outlined in the paper, or 
used the corrected data provided along with this package, please consider citing this and the following paper:

    Lognonné, P., W. B. Banerdt, W. T. Pike, D. Giardini, U. Christensen, R. F. Garcia, T. Kawamura, et al. “Constraints on the 
    Shallow Elastic and Anelastic Structure of Mars from InSight Seismic Data.” Nature Geoscience 13, no. 3 (March 2020): 
    213–20. https://doi.org/10.1038/s41561-020-0536-y.

.. _ObsPy: https://github.com/obspy/obspy/wiki
"""



### VERSION
__version__ = '1.0.0'
__author__  = 'John-Robert Scholz'