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
Its purpose is to detect "glitches" on the VBB (very broadband) and SP (short period)
seismometers, both part of the SEIS instrument package of NASA's InSight discovery mission to planet Mars. 
The package allows to plot the detected glitches, plot them in different ways, and remove them from the seismic raw data. 
Finally, there are some useful features implemented to ease data handling (e.g. download, decimation, time conversion).

| 
| Find the actual Python code at:
| https://pss-gitlab.math.univ-paris-diderot.fr/data-processing-wg/seisglitch
| 

The code is based on the following, peer-reviewed paper: 

    Scholz, J.-R., Widmer-Schnidrig, R., P. Davis, P. Lognonne, B. Pinot, R. F. Garcia, et al. 
    “Detection, Analysis and Removal of Glitches from InSight’s Seismic Data from Mars.” Earth and Space Science, submitted (2020).

| 

If you find ``SEISglitch`` useful, used the MATLAB alternatives, corrected data as to the considerations outlined in the paper, or 
used the corrected data provided along with this package, please consider citing. :)


.. _ObsPy: https://github.com/obspy/obspy/wiki
"""



### VERSION
__version__ = '1.0.0'
__author__  = 'John-Robert Scholz'