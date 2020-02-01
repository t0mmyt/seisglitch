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
import os
import sys
import scipy
import inspect
import warnings
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
import numpy as np


#####  matplotlib modules import  #####
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, FancyArrowPatch


#####  obspy modules import  #####
from obspy.core.stream import Trace
from obspy.signal import rotate


#####  seisglitch utils import  #####
from seisglitch.util import quick_plot, _Arrow3D
from seisglitch.math import covariance_matrix, rotate_2D, unit_vector


# data gaps
def ppol_fit(fit_func, BAZs, orientations, sigma=None, absolute_sigma=True, p0=None):

    """
    Perfom Ppol fit via a specified fit function.

    Detailed discription ...


    Parameters
    ----------
    func : function
        should return y-values in radiants.
    BAZs : list
        jkfdslagwfsil
    orientations : list
        dsftlögjedsflh
    sigma : bool
        sklgjsag
    absolute_sigma : bool
        djskaygj
    p0 : list
        dffe

    .. _scipy.optimize.curve_fit: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html
    """

    BAZs_rad           = BAZs * np.pi/180
    Orients_rad        = orientations * np.pi/180
    num_fit_parameters = len(inspect.signature(fit_func).parameters)-1   #-1, because x-values do not count (want only number fit parameters)


    # convert absolute sigmas into radiants
    if np.array(sigma).any():
        sigma_rad = sigma * np.pi/180
    else:
        sigma_rad = np.array([])
    

    # perform fit
    try:
        popt, pcov         = scipy.optimize.curve_fit(fit_func, BAZs_rad, Orients_rad, sigma=sigma_rad, absolute_sigma=absolute_sigma, p0=p0)
        perr               = np.sqrt(np.diag( pcov ))

        if not np.isnan(popt).all() and np.isnan(perr).all():
            raise ValueError

        return popt * 180/np.pi, perr * 180/np.pi

    
    except TypeError:       # e.g. number measurements m < num_fit_parameters
        print(u'WARNING: Ppol fit: %s measurements < %s fit parameters.' % (len(BAZs), num_fit_parameters))
        print(u'                   Fit values and errors set to nan. Adjust your conditions for better fits.')
        
        popt = np.ones( num_fit_parameters ) * np.nan
        perr = np.ones( num_fit_parameters ) * np.nan
        return popt, perr


    except ValueError:      # ppol fit too bad
        print(u'WARNING: Ppol fit: all errors are nan, indicating individual measurements do NOT cover enough backazimuthal range.')
        print(u'                   Fit values and errors set to nan. Adjust your conditions for better fits.')        
        
        popt = np.ones( num_fit_parameters ) * np.nan
        perr = np.ones( num_fit_parameters ) * np.nan
        return popt, perr
class ppol():

    """
    TO DO:
        - data gaps
    """

    def __init__(self, comp_N=None, comp_E=None, comp_Z=None, stream=None, demean=True, **kwargs):
        traces       = self._assign_traces(comp_N=comp_N, comp_E=comp_E, comp_Z=comp_Z, stream=stream)
        self.trace_N = traces[0].copy()
        self.trace_E = traces[1].copy()
        self.trace_Z = traces[2].copy()

        if demean:
            self.trace_N.detrend('demean')
            self.trace_E.detrend('demean')
            self.trace_Z.detrend('demean')
            self.demeaned = True
        else:
            self.demeaned = False

        self._sanity_check()
        self.results = self.calc(**kwargs)
    def __str__(self, tag='', print_3D=True):

        """
        When calling print(...)
        """

        if print_3D:
            string =    u'\n' \
                      + u'PPOL CALCULATIONS: %s   \n' % tag \
                      + u'   %11s : %s            \n' % ('demeaned',   self.demeaned) \
                      + u'   %11s : %s            \n' % ('NPTS',       len(self.trace_N.data)) \
                      + u'   %11s : %-5.1f ± %3.1f\n' % ('BAZ_2D',     self.results[0], self.results[1]) \
                      + u'   %11s : %-5.1f ± %3.1f\n' % ('BAZ_3D',     self.results[2], self.results[3]) \
                      + u'   %11s : %-5.1f ± %3.1f\n' % ('INC_2D',     self.results[4], self.results[5]) \
                      + u'   %11s : %-5.1f ± %3.1f\n' % ('INC_3D',     self.results[6], self.results[7]) \
                      + u'   %11s : %-10.1g       \n' % ('SNR_HOR_2D', self.results[8]) \
                      + u'   %11s : %-10.1g       \n' % ('SNR_3D',     self.results[9]) \
                      + u'   %11s : %-10.1g       \n' % ('SNR_RZp_2D', self.results[10]) \
                      + u'   %11s : %-10.1g       \n' % ('SNR_RZp_3D', self.results[11]) \
                      + u'   %11s : %-10.3f       \n' % ('POL_HOR_2D', self.results[12]) \
                      + u'   %11s : %-10.3f       \n' % ('POL_3D',     self.results[13]) \
                      + u'   %11s : %-10.3f       \n' % ('POL_RZp_2D', self.results[14]) \
                      + u'   %11s : %-10.3f       \n' % ('POL_RZp_3D', self.results[15]) \
                      + u'   %11s : %s            \n' % ('EigVecs_2D', self.results[20][0]) \
                      + u'   %11s : %s            \n' % (''          , self.results[20][1]) \
                      + u'   %11s : %s            \n' % ('EigVals_2D', self.results[21]) \
                      + u'   %11s : %s            \n' % ('EigVecs_3D', self.results[22][0]) \
                      + u'   %11s : %s            \n' % (''          , self.results[22][1]) \
                      + u'   %11s : %s            \n' % (''          , self.results[22][2]) \
                      + u'   %11s : %s              ' % ('EigVals_3D', self.results[23])
        else:
            string =    u'\n' \
                      + u'PPOL CALCULATIONS: %s   \n' % tag \
                      + u'   %11s : %s            \n' % ('demeaned',   self.demeaned) \
                      + u'   %11s : %s            \n' % ('NPTS',       len(self.trace_N.data)) \
                      + u'   %11s : %-5.1f ± %3.1f\n' % ('BAZ_2D',     self.results[0], self.results[1]) \
                      + u'   %11s : %-10.1g       \n' % ('SNR_HOR_2D', self.results[8]) \
                      + u'   %11s : %-10.3f       \n' % ('POL_HOR_2D', self.results[12]) \
                      + u'   %11s : %s            \n' % ('EigVecs_2D', self.results[20][0]) \
                      + u'   %11s : %s            \n' % (''          , self.results[20][1]) \
                      + u'   %11s : %s              ' % ('EigVals_2D', self.results[21])

        return string
    def _assign_traces(self, comp_N=None, comp_E=None, comp_Z=None, stream=None):

        """
        """

        if stream is not None:
            stream_N = stream.select(component='N') or stream.select(component='Q') or stream.select(component='U') or stream.select(component='1') or stream.select(component='T')
            stream_E = stream.select(component='E') or stream.select(component='T') or stream.select(component='V') or stream.select(component='2') or stream.select(component='R')
            stream_Z = stream.select(component='Z') or stream.select(component='L') or stream.select(component='W') or stream.select(component='3')

            if not stream_N and not stream_E:
                print(u'No idea how to perform polarization analysis on components: %s' % ', '.join( [tr.id for tr in stream] ))
                sys.exit()

            else:
                trace_N = stream_N[0]
                trace_E = stream_E[0]
                if stream_Z is not None:
                    trace_Z = stream_Z[0]
                else:
                    trace_Z = Trace()


        elif comp_N is not None and comp_E is not None:

            if isinstance(comp_N, Trace):
                trace_N = comp_N
            elif isinstance(comp_N, (list, tuple, np.ndarray)):
                trace_N = Trace(data=np.asarray(comp_N), header={'channel':'N'})
            else:
                print(u'`comp_N` must be either a ObsPy `Trace` object or a `list`-like object containing your waveform data.')
                sys.exit()

            if isinstance(comp_E, Trace):
                trace_E = comp_E
            elif isinstance(comp_E, (list, tuple, np.ndarray)):
                trace_E = Trace(data=np.asarray(comp_E), header={'channel':'E'})
            else:
                print(u'`comp_E` must be either a ObsPy `Trace` object or a `list`-like object containing your waveform data.')
                sys.exit()

            if comp_Z is not None:
                if isinstance(comp_Z, Trace):
                    trace_Z = comp_Z
                elif isinstance(comp_Z, (list, tuple, np.ndarray)):
                    trace_Z = Trace(data=np.asarray(comp_Z), header={'channel':'Z'})
                else:
                    print(u'`comp_Z` must be either a ObsPy `Trace` object or a `list`-like object containing your waveform data.')
                    sys.exit()                    
            else:
                trace_Z = Trace()


        else:
            print(u'You must either specify an ObsPy `stream` object containing at least two horizontal components')
            print(u'or `comp_N` and `comp_E` (both either as ObsPy `Trace` objects or lists).')
            sys.exit()


        return trace_N , trace_E, trace_Z
    def _sanity_check(self):

        
        len_N   = self.trace_N.stats.npts
        len_E   = self.trace_E.stats.npts
        len_Z   = self.trace_Z.stats.npts

        start_N = self.trace_N.stats.starttime
        start_E = self.trace_E.stats.starttime
        start_Z = self.trace_Z.stats.starttime

        if len_Z != 0:
            if len_N!=len_E or len_N!=len_Z or len_E!=len_Z:
                print(u'ERROR: Data lengths of components are not equal. Cannot perform ppol analysis.')
                sys.exit()

        else:
            if len_N!=len_E:
                print(u'ERROR: Data lengths of components are not equal. Cannot perform ppol analysis.')
                sys.exit()


        if len_Z != 0:
            if start_N!=start_E or start_N!=start_Z or start_E!=start_Z:
                print(u'WARNING: Data do not start at the same time. Analysis is performed nevertheless.')
                print('  '+self.trace_N)
                print('  '+self.trace_E)
                print('  '+self.trace_Z)

        else:
            if start_N!=start_E:
                print(u'WARNING: Data do not start at the same time. Analysis is performed nevertheless.')
                print('  '+self.trace_N)
                print('  '+self.trace_E)
    def calc(self, bias=False, fix_angles='EQ', Xoffset_samples_for_amplitude=None, **kwargs):

        r"""
        DO INDIVIDUAL PPOL MEASUREMENTS

        Take data arrays and perform 2-D and 3-D principle component
        analysis (2-D and 3-D PCA). Components must be of equal length
        and correspond to same start and thus end times.


        Useful for seismology, for example.

        |       ^  
        |       \| data_1
        |       \|
        |       \|
        |       \|
        |       x------------> data_2
        |     data_Z 
        |  (pointing to you)


        180° ambiguity:
        
          P-waves  -->  2-D data  -->  BAZ unknown  -->  180° ambiguity    (but if one knows expected first motion)
                   -->  3-D data  -->  BAZ unknown  -->  no 180° ambiguity (because P-wave must arrive from below, i.e. INC must >0°, or because of known, expected first motion)
          R-waves  -->  3-D data  -->  BAZ unknown  -->  no 180° ambiguity (retrogradicity, however, this algoirthmus does not treat surface / Rayleigh waves).

        Note
        ----

            An unknown event BAZ is the same as an unknown station orientation
            with a known event BAZ.

        This function does not demean data before running.
        """



        ### assign needed data
        data_1 = self.trace_N.data
        data_2 = self.trace_E.data
        data_Z = self.trace_Z.data



        ### IF NO Z-COMPONENT GIVEN
        if data_Z is None:
     

            ### 2-D phi, horizontal plane
            covariance_matrix_2D_hori      = covariance_matrix(data_1, data_2, bias=bias)                                  # does be default no demeaning internally!
            eig_val_2D, eig_vec_2D         = np.linalg.eig(covariance_matrix_2D_hori)
            index_array_descending_eig_val = eig_val_2D.argsort()[::-1]
            eig_val_2D                     = eig_val_2D[index_array_descending_eig_val]                                    # E-values descending
            eig_vec_2D                     = eig_vec_2D[:,index_array_descending_eig_val]                                  # E-vectors sorted acc. to E-values
            eig_vec_2D_1                   = eig_vec_2D[:,0]
            
            # Derived
            PHI_2D                         = (np.arctan2( eig_vec_2D_1[1], eig_vec_2D_1[0] ) * 180/np.pi )%360
            PHI_err_2D                     = np.arctan( np.sqrt( eig_val_2D[1]/eig_val_2D[0] )) * 180/np.pi
            SNR_HOR_2D                     = (eig_val_2D[0] - eig_val_2D[1]) / eig_val_2D[1]                               # De Meersman et al. (2006)
            POL_HOR_2D                     = 1 - eig_val_2D[1]/eig_val_2D[0]                                               # rectilinearity in horizontal plane (1 for linearised, 0 for circular polarisation). Jurkevics (1988)
            data_R_2D, data_T_2D           = rotate.rotate_ne_rt(data_1, data_2, PHI_2D)
            
            # Others
            PHI_3D                         = np.nan
            PHI_err_3D                     = np.nan
            INC_2D                         = np.nan
            INC_err_2D                     = np.nan
            INC_3D                         = np.nan
            INC_err_3D                     = np.nan
            SNR_3D                         = np.nan
            SNR_RZp_2D                     = np.nan
            SNR_RZp_3D                     = np.nan
            POL_3D                         = np.nan
            POL_RZp_2D                     = np.nan
            POL_RZp_3D                     = np.nan
            data_R_3D                      = np.nan*data_R_2D
            data_T_3D                      = np.nan*data_T_2D
            eig_val_3D                     = np.nan*np.ones(3)                                                            # E-values descending
            eig_vec_3D                     = np.nan*np.ones(3)                                     


            ### AMBIGUITY 180°
            if fix_angles.upper()=='EQ':         # Nothing we can do solving the 180° ambiguity in 2D (this code doesn't make use of first motion information)
                pass


            elif fix_angles.upper()=='AMP':
                data_1_offest = data_1[Xoffset_samples_for_amplitude:]
                data_2_offest = data_2[Xoffset_samples_for_amplitude:]
                
                amp_1         = data_1[np.argmax(np.abs(data_1_offest))]
                amp_2         = data_2[np.argmax(np.abs(data_2_offest))]
                
                PHI_2D_OLD    = PHI_2D

                # 2-D PHI
                if abs(amp_1)>=abs(amp_2):              # `data_1` more significant than `data_2` 
                    if amp_1>=0:                        # `data_1` positive
                        if PHI_2D>=90 and PHI_2D<180 or PHI_2D>=270 and PHI_2D<360:
                            PHI_2D = PHI_2D%180
                        else:
                            PHI_2D = PHI_2D%180 + 180
                    else:                               # `data_1` negative
                        if PHI_2D>=90 and PHI_2D<180 or PHI_2D>=270 and PHI_2D<360:
                            PHI_2D = PHI_2D%180 + 180
                        else:
                            PHI_2D = PHI_2D%180
                else:                                   # `data_2` more significant than `data_1` 
                    if amp_2>=0:                        # `data_2` positive
                        PHI_2D = PHI_2D%180 + 180
                    else:                               # `data_2` negative
                        PHI_2D = PHI_2D%180             

                # correct radial and transverse data
                data_R_2D, data_T_2D = rotate.rotate_ne_rt(data_1, data_2, PHI_2D)

                # 2-D INC
                if PHI_2D_OLD != PHI_2D:
                    INC_2D *= -1
                INC_2D = INC_2D % 180


            else:
                pass



        ### IF Z-COMPONENT GIVEN
        else:


            ### 2-D phi, horizontal plane
            covariance_matrix_2D_hori      = covariance_matrix(data_1, data_2, bias=bias)                               # does be default no demeaning internally!
            eig_val_2D, eig_vec_2D         = np.linalg.eig(covariance_matrix_2D_hori)
            index_array_descending_eig_val = eig_val_2D.argsort()[::-1]
            eig_val_2D                     = eig_val_2D[index_array_descending_eig_val]                                 # E-values descending
            eig_vec_2D                     = eig_vec_2D[:,index_array_descending_eig_val]                               # E-vectors sorted acc. to E-values
            eig_vec_2D_1                   = eig_vec_2D[:,0]
            
            # Derived
            PHI_2D                         = (np.arctan2( eig_vec_2D_1[1], eig_vec_2D_1[0] ) * 180/np.pi )%360
            PHI_err_2D                     = np.arctan( np.sqrt( eig_val_2D[1]/eig_val_2D[0] )) * 180/np.pi
            SNR_HOR_2D                     = (eig_val_2D[0] - eig_val_2D[1]) / eig_val_2D[1]                            # De Meersman et al. (2006)
            POL_HOR_2D                     = 1 - eig_val_2D[1]/eig_val_2D[0]                                            # rectilinearity in horizontal plane (1 for linearised, 0 for circular polarisation). Jurkevics (1988)
            data_R_2D, data_T_2D           = rotate.rotate_ne_rt(data_1, data_2, PHI_2D)


            ### 2-D phi, radial-vertical plane
            data_R_2D, data_T_2D             = rotate.rotate_ne_rt(data_1, data_2, PHI_2D)
            covariance_matrix_2D_radZ        = covariance_matrix(data_Z, data_R_2D, bias=bias)
            eig_val_2D_radZ, eig_vec_2D_radZ = np.linalg.eig(covariance_matrix_2D_radZ)
            index_array_descending_eig_val   = np.argsort( eig_val_2D_radZ )[::-1]
            eig_val_2D_radZ                  = eig_val_2D_radZ[index_array_descending_eig_val]                          # E-values descending
            eig_vec_2D_radZ                  = eig_vec_2D_radZ[:,index_array_descending_eig_val]                        # E-vectors sorted acc. to E-values
            eig_vec_2D_radZ_1                = eig_vec_2D_radZ[:,0]

            # derived
            INC_2D                           = np.arctan( eig_vec_2D_radZ_1[1] / eig_vec_2D_radZ_1[0] ) * 180/np.pi
            INC_err_2D                       = np.arctan( np.sqrt( eig_val_2D_radZ[1]/eig_val_2D_radZ[0] )) * 180/np.pi
            SNR_RZp_2D                       = (eig_val_2D_radZ[0] - eig_val_2D_radZ[1]) / eig_val_2D_radZ[1]           # De Meersman et al. (2006)
            POL_RZp_2D                       = 1 - eig_val_2D_radZ[1]/eig_val_2D_radZ[0]                                # rectilinearity in radial-vertical plane (1 for linearised, 0 for circular polarisation). Jurkevics (1988)


            ### 3-D
            covariance_matrix_3D             = covariance_matrix(data_1, data_2, data_Z, bias=bias)                     # does be default no demeaning internally!
            eig_val_3D, eig_vec_3D           = np.linalg.eig(covariance_matrix_3D)
            index_array_descending_eig_val   = np.argsort( eig_val_3D )[::-1]
            eig_val_3D                       = eig_val_3D[index_array_descending_eig_val]                               # E-values descending
            eig_vec_3D                       = eig_vec_3D[:,index_array_descending_eig_val]                             # E-vectors sorted acc. to E-values
            eig_vec_3D_1                     = eig_vec_3D[:,0]

            # derived
            PHI_3D                           = (np.arctan2( eig_vec_3D_1[1], eig_vec_3D_1[0] ) * 180/np.pi) % 360
            PHI_err_3D                       = np.arctan( np.sqrt( eig_val_3D[2]/(eig_val_3D[1]+eig_val_3D[0]) )) * 180/np.pi
            SNR_3D                           = abs((eig_val_3D[0] - (eig_val_3D[1] + eig_val_3D[2]))) / (eig_val_3D[1] + eig_val_3D[2])  # De Meersman et al. (2006)        
            POL_3D                           = 1 - ( eig_val_3D[1]+eig_val_3D[2] )/( 2*eig_val_3D[0] )                  # rectilinearity in 3-D. (1 for linearised, 0 for circular polarisation). Jurkevics (1988)
            data_R_3D, data_T_3D             = rotate.rotate_ne_rt(data_1, data_2, PHI_3D)


            ### 3-D phi, radial & Z-data plane
            data_R_3D, data_T_3D             = rotate.rotate_ne_rt(data_1, data_2, PHI_3D)
            covariance_matrix_3D_radZ        = covariance_matrix(data_Z, data_R_3D, bias=bias)                          # does be default no demeaning internally!
            eig_val_3D_radZ, eig_vec_3D_radZ = np.linalg.eig(covariance_matrix_3D_radZ)
            index_array_descending_eig_val   = np.argsort( eig_val_3D_radZ )[::-1]
            eig_val_3D_radZ                  = eig_val_3D_radZ[index_array_descending_eig_val]                          # E-values descending
            eig_vec_3D_radZ                  = eig_vec_3D_radZ[:,index_array_descending_eig_val]                        # E-vectors sorted acc. to E-values
            eig_vec_3D_radZ_1                = eig_vec_3D_radZ[:,0]

            # derived
            INC_3D                           = np.arctan( eig_vec_3D_radZ_1[1] / eig_vec_3D_radZ_1[0] ) * 180/np.pi
            INC_err_3D                       = np.arctan( np.sqrt( eig_val_3D_radZ[1]/eig_val_3D_radZ[0] )) * 180/np.pi
            SNR_RZp_3D                       = (eig_val_3D_radZ[0] - eig_val_3D_radZ[1]) / eig_val_3D_radZ[1]           # De Meersman et al. (2006)
            POL_RZp_3D                       = 1 - eig_val_3D_radZ[1]/eig_val_3D_radZ[0]                                # rectilinearity in radial-vertical plane (1 for linearised, 0 for circular polarisation). Jurkevics (1988)


            ### AMBIGUITY 180°
            if fix_angles.upper()=='EQ':    # However, correct baz must deliver incidence angle>0, therefore can solve ambiguity
                if INC_2D < 0:
                    PHI_2D               = (PHI_2D+180)%360
                    INC_2D               = abs(INC_2D)
                    data_R_2D, data_T_2D = rotate_2D(data_R_2D, data_T_2D, 180)
                if INC_3D < 0:
                    PHI_3D               = (PHI_3D+180)%360
                    INC_3D               = abs(INC_3D)
                    data_R_3D, data_T_3D = rotate_2D(data_R_3D, data_T_3D, 180)


            elif fix_angles.upper()=='AMP':
                data_1_offest = data_1[Xoffset_samples_for_amplitude:]
                data_2_offest = data_2[Xoffset_samples_for_amplitude:]
                data_Z_offest = data_Z[Xoffset_samples_for_amplitude:]
                
                amp_1         = data_1[np.argmax(np.abs(data_1_offest))]
                amp_2         = data_2[np.argmax(np.abs(data_2_offest))]
                amp_Z         = data_Z[np.argmax(np.abs(data_Z_offest))]
                
                PHI_2D_OLD    = PHI_2D
                PHI_3D_OLD    = PHI_3D

                # 2-D PHI
                if abs(amp_1)>=abs(amp_2):              # `data_1` more significant than `data_2` 
                    if amp_1>=0:                        # `data_1` positive
                        if PHI_2D>=90 and PHI_2D<180 or PHI_2D>=270 and PHI_2D<360:
                            PHI_2D = PHI_2D%180
                        else:
                            PHI_2D = PHI_2D%180 + 180
                    else:                               # `data_1` negative
                        if PHI_2D>=90 and PHI_2D<180 or PHI_2D>=270 and PHI_2D<360:
                            PHI_2D = PHI_2D%180 + 180
                        else:
                            PHI_2D = PHI_2D%180
                else:                                   # `data_2` more significant than `data_1` 
                    if amp_2>=0:                        # `data_2` positive
                        PHI_2D = PHI_2D%180 + 180
                    else:                               # `data_2` negative
                        PHI_2D = PHI_2D%180             
    
                # 3-D PHI
                if abs(amp_1)>=abs(amp_2):              # `data_1` more significant than `data_2` 
                    if amp_1>=0:                        # `data_1` positive
                        if PHI_3D>=90 and PHI_3D<180 or PHI_3D>=270 and PHI_3D<360:
                            PHI_3D = PHI_3D%180
                        else:
                            PHI_3D = PHI_3D%180 + 180
                    else:                               # `data_1` negative
                        if PHI_3D>=90 and PHI_3D<180 or PHI_3D>=270 and PHI_3D<360:
                            PHI_3D = PHI_3D%180 + 180
                        else:
                            PHI_3D = PHI_3D%180
                else:                                   # `data_2` more significant than `data_1` 
                    if amp_2>=0:                        # `data_2` positive
                        PHI_3D = PHI_3D%180 + 180
                    else:                               # `data_2` negative
                        PHI_3D = PHI_3D%180   

                # correct radial and transverse data
                data_R_2D, data_T_2D = rotate.rotate_ne_rt(data_1, data_2, PHI_2D)
                data_R_3D, data_T_3D = rotate.rotate_ne_rt(data_1, data_2, PHI_3D)

                # 2-D INC
                if PHI_2D_OLD != PHI_2D:
                    INC_2D *= -1
                INC_2D = INC_2D % 180

                # 3-D INC
                if PHI_3D_OLD != PHI_3D:
                    INC_3D *= -1
                INC_3D = INC_3D % 180


            else:
                pass
        


        ### RESULTS
        results = PHI_2D,     PHI_err_2D, PHI_3D,     PHI_err_3D, \
                  INC_2D,     INC_err_2D, INC_3D,     INC_err_3D, \
                  SNR_HOR_2D, SNR_3D,     SNR_RZp_2D, SNR_RZp_3D, \
                  POL_HOR_2D, POL_3D,     POL_RZp_2D, POL_RZp_3D, \
                  data_R_2D,  data_T_2D,  data_R_3D,  data_T_3D,  \
                  eig_vec_2D, eig_val_2D, eig_vec_3D, eig_val_3D



        ### ASSIGN / RETURN RESULTS
        self.results = results
        return self.results
    def plot(self, title='', verticals=(), outfile=None, show=True, **kwargs):

        """
        PLOT INDIVIDUAL PPOL MEASUREMENT


        Either provide num,py lists or Obspy traces.
        Always plots 2-D angle
        Expects that all traces have same start time and same amount of samples!
        """


        ### assign needed labels and data
        label_1 = '<--  %s  -->' %  self.trace_N.stats.channel
        label_2 = '<--  %s  -->' %  self.trace_E.stats.channel
        label_Z = '<--  %s  -->' %  self.trace_Z.stats.channel
        label_R = '<--  %s  -->' % (self.trace_E.stats.channel[:-1] + 'R')
        
        data_1  = self.trace_N.data
        data_2  = self.trace_E.data
        data_Z  = self.trace_Z.data


        ### Ppol result (some needed for plotting)
        #   only measures based on BAZ_2D printed, because they are always available
        BAZ_2D     = self.results[0]
        BAZ_err_2D = self.results[1]
        INC_2D     = self.results[4]
        INC_err_2D = self.results[5]
        data_R     = self.results[16]   # comp_R_2D
        eigvecs    = self.results[22]   # eig_vec_3D


        ### PLOT WAVEFORMS
        # Figure instance
        fig = plt.figure(figsize=(12,8))
        fig.canvas.set_window_title('Ppol plot individual measurement: %s' % title)
        fig.suptitle(title, fontsize=11)
        gs  = fig.add_gridspec(2, 3)

        ax1 = fig.add_subplot(gs[0, :])
        ax1 = quick_plot(self.trace_N, self.trace_E, self.trace_Z, verts=verticals, ylabel='Amplitudes', xlabel='Time', legend_loc='upper right', axis=ax1)


        ### SMALL HACKS to make sure for small amplitudes everything works out (precisely: arrow head of angle indicators)
        factor      = 1
        factor_str  = ''

        if data_Z is not None:
            if max( [max(abs(data_Z)), max(abs(data_1)), max(abs(data_2)), max(abs(data_R))] ) <= 1e-3:
                factor      = 1e9
                factor_str  = '%.0g' % (factor**-1)
                data_Z     *= factor
                data_1     *= factor
                data_2     *= factor
                data_R     *= factor

        else:
            if max( [max(abs(data_1)), max(abs(data_2)), max(abs(data_R))] ) <= 1e-3:
                factor      = 1e9
                factor_str  = '%.0g' % (factor**-1)
                data_1     *= factor
                data_2     *= factor


        ### PLOT PPOL
        # Variables needed for (beautiful) plotting
        colours = [[0, 0, 1-0.6*i/len(data_1)] for i in range(len(data_2))]
        
        maxi    = max( list(np.abs(data_1))+list(np.abs(data_2)) ) * 1.05
        
        xx      = np.linspace(-maxi, maxi, 100)
        yy      = 1/np.tan(BAZ_2D*np.pi/180) * xx
        yy_be   = 1/np.tan((BAZ_2D-BAZ_err_2D)*np.pi/180) * xx
        yy_af   = 1/np.tan((BAZ_2D+BAZ_err_2D)*np.pi/180) * xx
        
        x       = maxi/2*np.sin( (BAZ_2D-7)*np.pi/180 )
        y       = maxi/2*np.cos( (BAZ_2D-7)*np.pi/180 )
        x2      = maxi/2*np.sin( (BAZ_2D+2)*np.pi/180 )
        y2      = maxi/2*np.cos( (BAZ_2D+2)*np.pi/180 )

        # Set-up
        ax2 = fig.add_subplot(gs[1, 0])
        ax2.grid(ls='-.', lw=0.5)
        ax2.set(xlabel=label_2, ylabel=label_1)
        ax2.set_ylim([-maxi, maxi])
        ax2.set_xlim([-maxi, maxi])
        ax2.set_aspect('equal', adjustable='box')
        ax2.text(1, 0, title+' (%s points)' % len(data_1), ha='right', color='red', transform=ax2.transAxes, fontsize=6, bbox=dict(boxstyle='round,pad=0', fc='white', ec="white", lw=0))
        ax2.text(0, 1, factor_str, ha='center', color='black', transform=ax2.transAxes, fontsize=9, bbox=dict(boxstyle='round,pad=0', fc='white', ec="white", lw=0))
        
        ax2.spines['right'].set_color('none')
        ax2.spines['top'].set_color('none')
        ax2.spines['left'].set_color('none')
        ax2.spines['bottom'].set_color('none')  
        
        # Plot commands; data & Phi angle + erorrs
        if (BAZ_2D+BAZ_err_2D)//180 != BAZ_2D//180.001 or BAZ_2D == 0:
            ax2.fill_betweenx(yy_af, xx, facecolor='red', alpha=0.075)
            ax2.fill_betweenx(yy,    xx, facecolor='red', alpha=0.075)
        else:
            ax2.fill_between(xx, yy, yy_af, facecolor='red', alpha=0.075)       # error area of PHI
        if (BAZ_2D-BAZ_err_2D)//180 != BAZ_2D//180.001 or BAZ_2D == 0:
            ax2.fill_betweenx(yy_be, xx, facecolor='red', alpha=0.075)
            ax2.fill_betweenx(yy,    xx, facecolor='red', alpha=0.075)
        else:
            ax2.fill_between(xx, yy_be, yy, facecolor='red', alpha=0.075)       # error area of PHI
        ax2.plot([-maxi, maxi], [0, 0], 'k', lw=1)                              # centered horizontal line
        ax2.plot([0, 0], [-maxi, maxi], 'k', lw=1)                              # centered vertical line
        ax2.scatter(data_2, data_1, s=7, c=colours, zorder=3)                             # data
        ax2.plot( xx, yy, 'indianred', lw=1.5, zorder=4)                                  # PHI results
                
        # Angle arc + arrow head
        ax2.add_patch( Arc([0,0], maxi,  maxi, 90, -BAZ_2D, 0, color='indianred', lw=1.5, zorder=5))
        a = FancyArrowPatch([x,y], [x2,y2], mutation_scale=20, lw=1.5, arrowstyle="-|>", color="indianred", zorder=6)
        ax2.add_artist(a)       

        # Legend
        scatter_proxy = mpl.lines.Line2D([0],[0], c="indianred", marker='>')
        ax2.legend([scatter_proxy], ['BAZ_2D=(%.1f\u00B1%.1f)\u00b0' % (BAZ_2D, BAZ_err_2D)], numpoints=1, loc='upper right', prop={'size': 8})


        ## plot axis (part) 2 and 3, if there's Z-data
        if data_Z is not None:

            # Variables needed for (beautiful) plotting
            maxi  = max( list(np.abs(data_R))+list(np.abs(data_Z)) ) * 1.05
            
            xx    = np.linspace(-maxi, maxi, 100)
            yy    = 1/np.tan( INC_2D*np.pi/180) * xx
            yy_be = 1/np.tan((INC_2D-INC_err_2D)*np.pi/180) * xx
            yy_af = 1/np.tan((INC_2D+INC_err_2D)*np.pi/180) * xx
            
            x     = maxi/2*np.sin( (INC_2D-7)*np.pi/180 )
            y     = maxi/2*np.cos( (INC_2D-7)*np.pi/180 )
            x2    = maxi/2*np.sin( (INC_2D+2)*np.pi/180 )
            y2    = maxi/2*np.cos( (INC_2D+2)*np.pi/180 )
            
            # Set-up
            ax3 = fig.add_subplot(gs[1, 1])
            ax3.grid(ls='-.', lw=0.5, zorder=-1)
            ax3.set(xlabel=label_R, ylabel=label_Z)
            ax3.set_ylim([-maxi, maxi])
            ax3.set_xlim([-maxi, maxi])
            ax3.set_aspect('equal', adjustable='box')
            ax3.text(1, 0, title+' (%s points)' % len(data_2), horizontalalignment='right', color='red', transform=ax3.transAxes, fontsize=6, bbox=dict(boxstyle='round,pad=0', fc='white', ec="white", lw=0))
            ax3.text(0, 1, factor_str, ha='center', color='black', transform=ax3.transAxes, fontsize=9, bbox=dict(boxstyle='round,pad=0', fc='white', ec="white", lw=0))
        
            ax3.spines['right'].set_color('none')
            ax3.spines['top'].set_color('none')
            ax3.spines['left'].set_color('none')
            ax3.spines['bottom'].set_color('none')  
        
            # Plot commands; data & Phi angle + erorrs
            if (INC_2D+INC_err_2D)//180 != INC_2D//180.001 or INC_2D == 0:
                ax3.fill_betweenx(yy_af, xx, facecolor='red', alpha=0.075)
                ax3.fill_betweenx(yy,    xx, facecolor='red', alpha=0.075)
            else:
                ax3.fill_between(xx, yy, yy_af, facecolor='red', alpha=0.075)           # error area of INC
            if (INC_2D-INC_err_2D)//180 != INC_2D//180.001 or INC_2D == 0:
                ax3.fill_betweenx(yy_be, xx, facecolor='red', alpha=0.075)
                ax3.fill_betweenx(yy,    xx, facecolor='red', alpha=0.075)
            else:
                ax3.fill_between(xx, yy_be, yy, facecolor='red', alpha=0.075)           # error area of INC         
            ax3.plot([-maxi, maxi], [0, 0],  'k', lw=1)
            ax3.plot( [0, 0], [-maxi, maxi], 'k', lw=1)
            ax3.scatter(data_R, data_Z, s=5, c=colours, zorder=3)
            ax3.plot( xx, yy, 'indianred', lw=1.5, zorder=4)

            # Angle arc + arrow head
            ax3.add_patch( Arc([0,0], maxi,  maxi, 90, -INC_2D, 0, color='indianred', lw=1.5, zorder=5))
            a = FancyArrowPatch([x,y], [x2,y2], mutation_scale=20, lw=1.5, arrowstyle="-|>", color="indianred", zorder=6)
            ax3.add_artist(a)
            
            # Legend
            scatter_proxy = mpl.lines.Line2D([0],[0], c="indianred", marker='>')
            ax3.legend([scatter_proxy], ['INC_2D=(%.1f\u00B1%.1f)\u00b0' % (INC_2D, INC_err_2D)], loc='upper right', prop={'size': 8})



            ## 3-D plot
            # Variables needed for (beautiful) plotting
            mean_N   = np.average(data_1)
            mean_E   = np.average(data_2)
            mean_Z   = np.average(data_Z)
            maxi     = max( list(np.abs(data_Z))+list(np.abs(data_1))+list(np.abs(data_2)) ) * 1.05
            max_dist = np.sqrt(3*(maxi*2)**2)

            # Set-up
            ax4 = fig.add_subplot(gs[1, 2], projection='3d')
            ax4.set_xlabel(label_2, labelpad=5, fontsize=10)
            ax4.set_ylabel(label_1, labelpad=5, fontsize=10)
            ax4.zaxis.set_rotate_label(False)    # disable automatic rotation
            ax4.set_zlabel(label_Z, labelpad=5, fontsize=10, rotation=90)
            ax4.set_xlim( [-maxi,maxi] )
            ax4.set_ylim( [-maxi,maxi] )
            ax4.set_zlim( [-maxi,maxi] )
            ax4.view_init(elev=15, azim=-135)
            ax4.text(0, 1, 1, factor_str, ha='center', color='black', transform=ax4.transAxes, fontsize=9, bbox=dict(boxstyle='round,pad=0', fc='white', ec="white", lw=0))


            # Plot commands
            ax4.scatter(data_2, data_1, data_Z, c=colours, s=5, depthshade=False, zorder=3)
            ax4.scatter([mean_E], [mean_N], [mean_Z], s=20, c='white', edgecolors='k', depthshade=False, zorder=4)

            # eigenvectors
            eigvecs  = [ unit_vector(eigvecs[:,0]), unit_vector(eigvecs[:,1]), unit_vector(eigvecs[:,2]) ]
            for l, eigvec in enumerate( eigvecs ):
                length = (0.7-l*0.15) * max_dist/2      # arbitrary length, not related to actual eigen values
                a = _Arrow3D(np.asarray( [mean_E, eigvec[1]*length] ),
                             np.asarray( [mean_N, eigvec[0]*length] ),
                             np.asarray( [mean_Z, eigvec[2]*length] ),
                             mutation_scale=20, lw=1.5, arrowstyle="-|>", color="indianred", zorder=5)
                ax4.add_artist(a)

            # legend
            scatter_proxy = mpl.lines.Line2D([0],[0], c="indianred", marker='>')
            ax4.legend([scatter_proxy], ['eigen vectors'], loc='upper right', prop={'size': 8})


        ## Figure save / close
        plt.tight_layout(rect=[0, 0.02, 1, 0.98])
        if outfile:
            plt.savefig(outfile)
        if show:
            plt.show()
        plt.close()
    def display(self, tag='', print_3D=True):

        """
        Print the results of `calc`.
        """

        print(self.__str__(tag=tag, print_3D=print_3D))        


#####  _ _ N A M E _ _ = = " _ _ M A I N _ _ "  #####
if __name__ == "__main__":
    pass