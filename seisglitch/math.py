#!/usr/bin/env python

# Copyright 2019 John-Robert Scholz
#
# This file is part of Ppol.
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


#####  python modules import  #####
import numpy as np


def covariance_matrix(*observables, demean_before=False, bias=False):
    
    r"""
    Calculate covariance matrix using numpy.
    Refer to:
        http://docs.scipy.org/doc/numpy/reference/generated/numpy.cov.html
    
    Basic:
        COVA_mat = Y'\*Y/N (biased) OR Y'\*Y/(N-1) (unbiased), where Y is centered data matrix with 
        len(observables) as columns and observations i .. N as rows
    """

    
    ### put passed data of observable into matrix (check of each observable has N datapoints is later)
    matrix = np.array( observables )
    
    
    ### check if observable objects have same length
    try:
        matrix.shape[1]
    except IndexError:      # can't index number of columns because they are different for the rows --> not same length of passed data
        sys.exit( "ERROR: Lengths of observables don't match." )
    

    ### de-meaning
    if demean_before:
        # centered data matrix with columns as observables (2 or 3) and rows as observations i .. N
        mean        = matrix.mean(axis=1)   # mean along rows
        matrix      = (matrix.T - mean).T   # remove mean from rows with tiny hack
    

    ### Calculate covariance matrix : Y'*Y/N (biased) OR Y'*Y/(N-1) (unbiased) 
    cova_matrix = np.cov( matrix, bias=bias )
    
    return cova_matrix
def rotate_2D(comp_1, comp_2, angle, clockwise=True):

    r"""
    This function rotates 2 data arrays (e.g., seismological component data) 
    into the desired coordinate system using `angle` (clockwise direction by default).

    The given algorithm works for 'comp_2' being oriented 90° clockwise to 'comp_1'.

    ^  
    \| comp_1
    \|
    \|
    \|
    o------> comp_2


    :type comp_1:  np.array
    :param comp_1: List with floats representing data.
                   Note 'comp_1' is oriented 90° counter-clockwise to 'comp_2'.

    :type comp_2:  np.array
    :param comp_2: List with floats representing data.
                   Note 'comp_2' is oriented 90° clockwise to 'comp_1'.

    :type angle:  float
    :param angle: Angle about which both data traces are rotated. Can be negative.

    :type clockwise:  bool
    :param clockwise: if True  :         clockwise rotation of both data traces 
                      if False : counter-clockwise rotation of both data traces 

    :type return:  np.array, np.array
    :param return: rotated data lists
    """

    if not clockwise:
        angle = 360-angle               # convert angle to as it would be clockwise rotation


    comp_1_new =  comp_1*np.cos(angle*np.pi/180) + comp_2*np.sin(angle*np.pi/180)
    comp_2_new = -comp_1*np.sin(angle*np.pi/180) + comp_2*np.cos(angle*np.pi/180)

    return comp_1_new, comp_2_new
def unit_vector(vec):

    """ 
    Returns the unit vector of the vector (=array).  
    """

    return vec / np.linalg.norm(vec)
def circ_stat(*angles, axis = None):

    """
    The arithmetic mean / median is not always appropriate for angles.
    For this, polar angles are converted to cartesian coordinates, the mean / median
    is then calculated, and the result finally converted back to polar coordinates.

    Check for example:
    https://en.wikipedia.org/wiki/Directional_statistics
    
    ---------

    INPUT: *angles: either angles to be transformed seperated by comma
                    or one list with all angles inside. 

                    All values must be floats, not strings.


    :param mean: circular mean of angles

    :param ang_deviation: corresponds to one standard deviation
                          This estimate is used for circular mean errors in the
                          Ppol package.
                          Paper:
                            Berens, P., 2009. CircStat: a MATLAB toolbox for circular 
                            statistics, J. Stat. Softw., 31(10), 1–21.

    :param circ_STD: circular standard deviation, according to common sources.
                     This estimate is NOT used for circular mean errors in the
                     Ppol package.

    :param median: circular median of angles. If len(list) = odd, median is 
                   value in the middle. If len(list) = even, left &
                   right values are interpolated.

    :param mad: median absolute deviation, according to common sources
                This estimate is NOT used for circular median errors in the
                Ppol package.

    :param smad: mad times 1.4826, is approx. standad deviation
                 (for normally, gaussian distributed data)
                 This estimate is used for circular mean errors in the
                 Ppol package.                 
                 Paper: 
                   Rousseeuw, P.J. & Croux, C., 1993. Alternatives to the median 
                   absolute deviation, J. Am. Stat. Assoc., 88(424), 1273–1283

    :return: mean, ang_deviation, circ_STD, median, mad, smad
    """



    ### CONVERT ANGLES FROM DEGREE INTO RADIANTS
    if isinstance( angles[0], (list, tuple, np.ndarray)):
        angles = angles[0]
    angles_rad = np.array( angles ) / 180*np.pi



    ### CIRCUCLAR MEAN
    mean_rad = np.arctan2( np.mean( np.sin(angles_rad),axis ), np.mean( np.cos(angles_rad),axis ) )
    mean     = (mean_rad * 180/np.pi) % 360

    # calc length of  mean resultant vector of  circular distribution
    if np.ma.isMaskedArray(angles_rad) and angles_rad.mask.shape!=():
        N = np.sum(~angles_rad.mask,axis)
    else:
        if axis is None:
            N = angles_rad.size
        else:
            N = angles_rad.shape[axis]
    R = np.sqrt( np.sum(np.sin(angles_rad),axis)**2 + np.sum(np.cos(angles_rad),axis)**2 ) / N

    # error estimates
    circ_Var      = 1-R                                         # circular variance
    circ_STD      = np.sqrt( -2*np.log(R) ) * 180/np.pi         # circular STD: values between 0 and infinity (without term: 180/np.pi)
    ang_deviation = np.sqrt(  2*(1-R)     ) * 180/np.pi         # values between 0 and sqrt(2)  (without term: 180/np.pi) check:  Berens, P. (2009). CircStat: A MATLAB toolbox for circular statistics
    #print('circ stat: ', mean, ang_deviation, circ_STD)


    ### CIRCULAR MEDIAN
    median_rad = np.arctan2( np.median( np.sin(angles_rad),axis ), np.median( np.cos(angles_rad),axis ) )
    median     = (median_rad * 180/np.pi) % 360

    mad_rad = np.arctan2( np.median( np.absolute(np.sin(angles_rad-median_rad)) ), np.median( np.absolute( np.cos(angles_rad-median_rad) )) )
    mad = mad_rad * 180/np.pi
    smad = 1.4826 * mad



    ### RETURN
    return mean, ang_deviation, circ_STD, median, mad, smad
def fit_function1(x, a, b, c, d, e):
    
    """
    `x` data should be in radiants.
    """
    
    return a + b * np.sin(x) + c * np.cos(x) + d * np.sin(2*x) + e * np.cos(2*x)    # 
def split_into_connecting_parts(array, distance=1):

    """
    """

    i    = 0
    last = 0

    while True:
        try:
            if array[i+1]-array[i] != distance:
                yield array[last:i+1]
                last = i+1
        
            i += 1
                
        except: # last index reached
            yield array[last:]
            break
def moving_window(data, window_length_in_samples=100, step_in_samples=50, equal_end=True):

    """
    Yield data of moving windows according to given parameters.
    """

    i = 0

    while True:

        # window indices
        index_window_start = i * step_in_samples
        index_window_end   = index_window_start + window_length_in_samples


        # latest when start time of window is >= number samples, then no more windows
        i += 1
        if index_window_start >= len(data):
            break
       
        # insist last window is of length `window_length_in_samples` 
        if equal_end:
           
            if len(data[index_window_start:index_window_end]) < window_length_in_samples:
                break

        yield index_window_start, index_window_end, data[index_window_start:index_window_end]
def rotate_2D(comp_1, comp_2, angle, clockwise=False):

    """
    This function rotates 2 data traces (i.e. numpy arrays, e.g., seismological data) 
    into the desired coordinate system using 'angle' (clockwise direction by default).

    The given algorithm works for 'comp_2' being oriented 90° clockwise to 'comp_1'.

    ^  
    | comp_1
    |
    |
    |
    o------> comp_2


    :type comp_1:  np.array
    :param comp_1: List with floats representing data.
                   Note 'comp_1' is oriented 90° counter-clockwise to 'comp_2'.

    :type comp_2:  np.array
    :param comp_2: List with floats representing data.
                   Note 'comp_2' is oriented 90° clockwise to 'comp_1'.

    :type angle:  float
    :param angle: Angle about which both data traces are rotated. Can be negative.

    :type clockwise:  bool
    :param clockwise: if True  :         clockwise rotation of both data traces 
                      if False : counter-clockwise rotation of both data traces 

    :type return:  np.array, np.array
    :param return: rotated data lists
    """

    if not clockwise:
        angle = 360-angle               # convert angle to as it would be clockwise rotation


    comp_1_new =  comp_1*np.cos(angle*np.pi/180) + comp_2*np.sin(angle*np.pi/180)
    comp_2_new = -comp_1*np.sin(angle*np.pi/180) + comp_2*np.cos(angle*np.pi/180)

    return comp_1_new, comp_2_new
def normalise(data, scale_to_between=[]):

    """
    Normalise passed data (array-like).
    `scale_to_between` is list with 2 elements.

    Returns data as was if length of data is one or two.
    """

    data = np.asarray( data )
    
    if len(data)==0 or len(data)==1:
        return data

    if isinstance(scale_to_between, (int, float)):
        scale_to_between = [scale_to_between]

    if scale_to_between:
        if len(scale_to_between) == 1:
            scale_to_between = [0, scale_to_between[0]]
        scale_to_between.sort()

        scale  = np.abs(scale_to_between[-1]-scale_to_between[0]) / 2.
        drange = np.max(data)-np.min(data)
        data   = data * 2 / drange
        data   = data - np.max(data)+1             # data have y values between [-1,1]
        data  *= scale                          # data have y values between [-1/scale,1/scale] 
        data  += scale_to_between[-1]-scale     # eacht trace has y values filling range `scale_to_between`
    
    else:
        data /= np.max(np.abs(data))

    return data