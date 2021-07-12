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
*SEISglitch* is a seismological Python package based on the ObsPy_ library.
Its purpose is to handle "glitches"  - also referred to as "long-period data disturbances" - present in the raw seismic data of the VBB (very broadband) 
and SP (short period) seismometers, both part of the *SEIS* instrument of NASA's "InSight_" Discovery mission to planet Mars. 
The package allows to detect and remove glitches from the seismic raw data, and provides plot funtionality to analyse their behaviour in detail. 
SEISglitch corresponds to the 'MPS' method detailed in Scholz et al. (2020, see below). All other methods
('ISAE', 'UCLA', 'IPGP') were implemented in MATLAB and their essential scripts are also shipped with this software.
SEISglitch further comes with some useful features to ease data handling, e.g. download, processing, and time conversion.
Potentially, major functionalities of this toolbox could also be used to treat the same 
data disturbances occurring on other seismic stations, however, this has not been tested.

| Find the Python source code at:
| https://pss-gitlab.math.univ-paris-diderot.fr/data-processing-wg/seisglitch

The package is based on the following, peer-reviewed paper: 
::

    Scholz, J.‐R., Widmer‐Schnidrig, R., Davis, P., Lognonné, P., Pinot, B., Garcia, R. F., 
    et al. (2020). Detection, analysis and removal of glitches from InSight’s seismic 
    data from Mars. Earth and Space Science, 7, e2020EA001317‐T. 
    https://doi.org/10.1029/2020EA001317

If you find ``seisglitch`` useful, used the MATLAB alternatives, corrected data as to the considerations outlined in the paper, or 
used the corrected data provided along with this package, please consider citing this and the following paper:
::

    Lognonné, P., W. B. Banerdt, W. T. Pike, D. Giardini, U. Christensen, R. F. Garcia, 
    T. Kawamura, et al. “Constraints on the Shallow Elastic and Anelastic Structure of
    Mars from InSight Seismic Data.” Nature Geoscience 13, no. 3 (March 2020), 213–20. 
    https://doi.org/10.1038/s41561-020-0536-y.

You can also help the Obspy developers by citing_ their work that the seisglitch package heavily relies on.

Note: Along with this package also come deglitched data for a selection of seismic events. 
These corrected raw data will automatically be on your machine if you follow the :ref:`installation`.
If you just want to quickly download the deglitched data, go here_.

.. _ObsPy: https://github.com/obspy/obspy/wiki
.. _InSight: https://mars.nasa.gov/insight/
.. _citing: https://github.com/obspy/obspy/wiki#acknowledging
.. _here: https://pss-gitlab.math.univ-paris-diderot.fr/data-processing-wg/seisglitch/tree/master/DEGLITCHED_DATA
"""



### VERSION
__version__ = '1.0.1'
__author__  = 'John-Robert Scholz'