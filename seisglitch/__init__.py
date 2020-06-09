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
Its purpose is to detect glitches on InSight's seismometers VBB (very broadband) and SP (short period)
that are both part of the SEIS instrument package, plot the detected glitches in different ways,
and remove the glitches.

| 
| Find the actual Python code at:
| https://pss-gitlab.math.univ-paris-diderot.fr/data-processing-wg/seisglitch
| 

The code is based on the following, peer-reviewed paper. If you find ``SEISglitch`` useful, please consider citing. :)

    Scholz, J.-R., Widmer-Schnidrig, R., P. Davis, P. Lognonne, B. Pinot, R. F. Garcia, Francis Nimmo, et al. 
    “Detection, Analysis and Removal of Glitches from InSight’s Seismic Data from Mars.” Journal of Geophysical Research: 
    Planets submitted (2020).


.. _ObsPy: https://github.com/obspy/obspy/wiki
"""



### VERSION
__version__ = '0.0.3'
__author__  = 'John-Robert Scholz'