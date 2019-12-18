#!/usr/bin/env python
# -*- coding: utf-8 -*-

#----------------------------------------------------------------------
#   Filename:  toolbox.py
""" Functions for seismological applications, using Python & ObsPy. """
#   Author:    John - Robert Scholz
#   Date:      Jan 2014 - present
#   Email:     john.robert.scholz@gmail.com
#   License:   TBD
#---------------------------------------------------------------------



###############  python 3.7 modules import  ######################
import os
import re
import io
import csv
import sys
import copy
import glob
import time
import pywt
import dtcwt
import scipy
import random
import fnmatch
import datetime
import threading
import collections
import math as M
import numpy as np
import pickle
import shutil
try:
	import imageio
	import hdbscan
	import seaborn as sns
	from sklearn.preprocessing import QuantileTransformer, minmax_scale
except:
	pass
from collections import Counter
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)


###############  mtplotlib modules import  ######################
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.switch_backend('TKAgg')
import matplotlib.dates as mdates
import matplotlib.ticker as mticker
import matplotlib.gridspec as gridspec
from matplotlib.widgets import TextBox
from mpl_toolkits.mplot3d import Axes3D, proj3d
from matplotlib.patches import Arc, FancyArrowPatch
from matplotlib.colors import Normalize
import matplotlib.cm as cmx


###############  obspy modules import  ######################
import obspy
from obspy import read_inventory, read_events, UTCDateTime
from obspy.io.mseed import InternalMSEEDError
from obspy.core.inventory import Inventory, Network, Station, Channel, Site, Comment
from obspy.core.event import Catalog, Origin, Event, Magnitude, CreationInfo, ResourceIdentifier
from obspy.imaging.cm import obspy_sequential
from obspy.core.trace import Trace
from obspy.core.stream import Stream
from obspy.geodetics.base import degrees2kilometers, locations2degrees, gps2dist_azimuth
from obspy.signal import rotate, trigger
from obspy.signal.filter import envelope
from obspy.taup import tau
from mars_tools.insight_time import solify, UTCify


######################  code  ######################
# convenient helpers
def c():

	""" 
	Clear screen in python. 
	""" 
	
	os.system('clear')
def pwd():
	"""
	print current working directory.
	"""
	return( os.getcwd() )
def sec2hms(seconds, digits=0):

	"""
	Convert seconds given as float into hours, minutes and seconds.
	Optional 'digits' determines position after decimal point.

	If variable 'seconds' is not understood, this function returns
	zeroes.
	
	:type sting: unicode string
	:return string: formated seconds (hh:mm:ss or hh:mm:ss.sss ..)
	"""

	try:
		seconds = float(seconds)
		m, s    = divmod(round(seconds,digits), 60)
		h, m    = divmod(m, 60)

		mask    = "{0:.%sf}" % digits
		s       = mask.format(s).zfill(digits+3)

		if digits == 0:
			return u'%s:%s:%s' % (str( int(h) ).zfill(2), str( int(m) ).zfill(2), str( int( round( float(s) ))).zfill(2))
		else:
			return u'%s:%s:%s' % (str( int(h) ).zfill(2), str( int(m) ).zfill(2), s)

	except Exception:
		return sec2hms(0, digits=digits)
def get_files(this, pattern='*'):

	"""
	Return all files in dir 'this' (and sub-directories) if 'this' is dir,
	or return file (if 'this' is file). You may use 'pattern' to search only for
	certain file name patterns (eg: file endings, station names, ..).

	Script uses os.walk and thus lists all files in a directory, meaning also
	those in sub-directories. Bash wildcards are understood, python regex not. 

	:type this: str
	:param this: file or dir (bash wildcards work)

	:type pattern: str
	:param pattern: pattern to filter file(s) (bash wildcards work)
	
	:return: list of files (or empty list, if no matching files were found)
	"""
	
	### get files from 'this' to loop over
	uber  = glob.glob(this)
	files = []

	try:
		for u in uber:
			if os.path.isfile(u):
				files.append(u)

			else:
				# walk trough (sub)-directory and put ALL files into 'files'
				for root, dirs, files_walk in os.walk(u):
					for file_walk in files_walk:

						# take care on 'pattern'
						file = os.path.join( root,file_walk )
						if fnmatch.fnmatch( os.path.basename(file), pattern ) and not file in files:
							files.append( file )

	except Exception:
		files = []

	return sorted(files)
def matlab2datetime(matlab_datenum):
	"""
	Corrects year offset usually encountered when converting from
	MatLab times to python datetime objects.

	Check explanation here:
	http://sociograph.blogspot.com/2011/04/how-to-avoid-gotcha-when-converting.html

	Basically: Python uses traditional propleptic Gregorian calendar
	           MatLab uses proleptic Gregorian calendar (ISO 8601) 
	"""
	day     = datetime.datetime.fromordinal(int(matlab_datenum))
	dayfrac = datetime.timedelta(days=matlab_datenum%1) - datetime.timedelta(days = 366)
	return day + dayfrac
def correct_CRchars(file):
	"""
	CR chars: carriage return characters
	
	See: 
	https://stackoverflow.com/questions/19425857/env-python-r-no-such-file-or-directory
	"""

	with open(file, 'rb+') as f:
		content = f.read()
		f.seek(0)
		f.write(content.replace(b'\r', b''))
		f.truncate()
def moving_window(data, window_length_in_samples=100, step_in_samples=50, equal_end=True):

	"""
	Yield data of moving windows according to given parameters.
	"""


	data = np.array( data )
	
	i = 0
	while True:

		# window indices
		window_start = i * step_in_samples
		window_end   = window_start + window_length_in_samples

		# latest when start time of window is >= number samples, then no more windows
		i += 1
		if window_start >= len(data):
			break

		# insist last window is of length `window_length_in_samples` 
		if equal_end:
			if len(data[window_start:window_end]) < window_length_in_samples:
				break


		yield window_start, window_end, data[window_start:window_end]
class Logger(object):

	""" 
	Print commands are shown in shell and written to 'logfile'
	at the same time !
	Usage: 
		sys.stdout = Logger(logfile)	# to start the magic
		sys.stdout = sys.__stdout__		# reset stdout to normal
	"""

	def __init__(self, logfile):
		self.terminal = sys.stdout
		self.logfile  = logfile
		self.log      = open(self.logfile, "w")

	def write(self, message):
		self.terminal.write(message)
		self.log.write(message)  

	def flush(self):
		#this flush method is needed for python 3 compatibility.
		#this handles the flush command by doing nothing.
		#you might want to specify some extra behavior here.
		pass 

		# general helpers
class Reprinter(object):

	"""
	Shall allow to reprint the terminal, i.e.,
	to move overprint some already printed lines.
	Would require some testion, but in principle
	is finished.
	"""

	def __init__(self):
		self.text = ''

	def moveup(self, lines):
		for _ in range(lines):
			sys.stdout.write("\x1b[A")

	def reprint(self, text):
		# Clear previous text by overwritig non-spaces with spaces
		self.moveup(self.text.count("\n"))
		sys.stdout.write(re.sub(r"[^\s]", " ", self.text))

		# Print new text
		lines = min(self.text.count("\n"), text.count("\n"))
		self.moveup(lines)
		sys.stdout.write(text)
		self.text = text
class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)

# mathematical
def union(a,b, reverse=False):

	""" 
	Form sorted union set of list a and list b. 
	Spoken: a united b
	
	:return: sorted union list of passed lists a & b / empty list if passed lists invalid
	"""
	
	if isinstance(a, tuple):
		a = list(a)
	if isinstance(b, tuple):
		b = list(b)
	
	if isinstance(a, list) and isinstance(b, list):
		liste = list( set( a + b ))
	else:
		print( colour('WARNING', Fore.YELLOW) + '%sOne or both of passed lists ' \
				'are no lists (neither tuples). Return empty list.' % tab(4, 8))
		return []
	
	if isinstance(reverse, bool):
		liste.sort(reverse=reverse)
	else:
		liste.sort(reverse=False)
	
	return liste
def inter(a,b, reverse=False):
	
	""" 
	Form intersection set of list a and list b. 
	Spoken: subset both sets share
		
	:return: sorted inter list of passed lists a & b / empty list if passed list invalid
	"""
	
	if isinstance(a, tuple):
		a = list(a)
	if isinstance(b, tuple):
		b = list(b)
	
	if isinstance(a, list) and isinstance(b, list):
		liste = list( set( [element for element in a if element in b] ))
	else:
		print( colour('WARNING', Fore.YELLOW) + '%sOne or both of passed lists ' \
				'are no lists (neither tuples). Return empty list.' % tab(4, 8))
		return []

	if isinstance(reverse, bool):
		liste.sort(reverse=reverse)
	else:
		liste.sort(reverse=False)

	return liste
def diff(a,b, reverse=False):

	""" 
	Form differtial set (relative complement). 
	Spoken: a without b
	
	:return: sorted diff list of passed lists a & b / empty list if passed lists invalid
	"""
	
	if isinstance(a, tuple):
		a = list(a)
	if isinstance(b, tuple):
		b = list(b)
	
	if isinstance(a, list) and isinstance(b, list):
		liste = list( set( [element for element in a if element not in b] ))
	else:
		print( colour('WARNING', Fore.YELLOW) + '%sOne or both of passed lists ' \
				'are no lists (neither tuples). Return empty list.' % tab(4, 8))
		return []
	
	if isinstance(reverse, bool):
		liste.sort(reverse=reverse)
	else:
		liste.sort(reverse=False)	

	return liste
def median(*angles):

	"""
    Compute 'Median Absolute Deviation' (MAD) of an array along given axis,
    using numpy.

    :type median:  float
    :param median: median of passed list. If len(list) = odd, median is 
                   value in the middle. If len(list) = even, left &
                   right values are interpolated.

   
    :type mad:  float
    :param mad: median absolute deviation, according to common sources

    :type smad:  float
    :param smad: mad times 1.4826, is approx. standad deviation
                 (for normally, gaussian distributed data)

    :return: median, mad, smad
	"""

	axis=None
	if isinstance( angles[0], (list, tuple, np.ndarray)):
		angles = angles[0]

	median = np.median( angles, axis=axis )
	mad    = np.median( np.absolute(angles - median), axis=axis )
	smad   = 1.4826 * mad

	print( str(median), str(mad), str(smad) )
	return( median, mad, smad )
def circ_stat(*angles):

	"""

	The arithmetic mean is not always appropriate for angles.
	For this, polar angles are converted to cartesian corrdinates, the arithmetic mean
	is then calculated, the result converted back to polar coordinates.
	This angle is a reasonable mean of the input angles. 

	Also, the variance of this circular mean can be calculated.

	---------

	INPUT - angles: either angles to be transformed seperated by comma
					or one list with all angles inside. 

	                All values must be floats, not strings.

	To achieve correct calculation, the input angles 
	(starting from 0, ascending clockwise from north) 
	must be converted first to polar angles 
	(starting from 0, counter clockwise from east) !

	After the mean calculation, angles are transformed back to 'North'.
	Results are rounded at the 8th digit. 
	Numpy axis is by default 'None'.

	:type final_angle: float
	:type variance: float

	:return final_angle, variance: results are printed AND returned
	"""

	axis = None

 	# angles come in degree, convert to rad
	if isinstance( angles[0], (list, tuple, np.ndarray)):
		angles = angles[0]
	angles_rad = np.array( angles ) / 180*M.pi

 	# calc circular mean
	mean_angle_rad = np.arctan2( np.mean( np.sin(angles_rad),axis ), np.mean( np.cos(angles_rad),axis ) )
	mean_angle     = mean_angle_rad*180/M.pi
	if mean_angle < 0:
		mean_angle += 360	

	# calc length of  mean resultant vector of  circular distribution
	if np.ma.isMaskedArray(angles_rad) and angles_rad.mask.shape!=():
	    N = np.sum(~angles_rad.mask,axis)
	else:
	    if axis is None:
	        N = angles_rad.size
	    else:
	        N = angles_rad.shape[axis]
	R = np.sqrt( np.sum(np.sin(angles_rad),axis)**2 + np.sum(np.cos(angles_rad),axis)**2 ) / N
	
	# VAR, angular deviation and circular standard deviation in deg 
	# (e.g.: Berens, P. (2009). CircStat: A MATLAB toolbox for circular statistics, J. Stat. Software 31, no. 10, 21 pp.)
	circ_Var      = 1-R
	ang_deviation = np.sqrt( 2 - 2*R )*180/M.pi 	 		# values between 0 and sqrt(2)  (without term: 180/M.pi)
	circ_STD      = np.sqrt( -2*M.log(R, M.e) )*180/M.pi   	# values between 0 and infinity (without term: 180/M.pi)


	### return
	#print(str(round(mean_angle, 8)), str(round(circ_Var, 8)), str(round(ang_deviation, 8)), str(round(	circ_STD, 8)))
	return(round(mean_angle, 8), round(ang_deviation, 8), round(circ_Var, 8), round( circ_STD, 8))
def covariance_matrix(*observables):
	
	"""
	Calculate covariance matrix using numpy.
	Refer to:
		http://docs.scipy.org/doc/numpy/reference/generated/numpy.cov.html
	
	Basic:
		COVA_mat = Y'*Y/N (biased) OR Y'*Y/(N-1) (unbiased), where Y is centered data matrix with 
		len(observables) as columns and observations i .. N as rows
	"""

	
	### put passed data of observable into matrix (check of each observable has N datapoints is later)
	matrix = np.array( observables )
	
	
	### check if observable objects have same length
	try:
		matrix.shape[1]
	except IndexError:		# can't index number of columns because they are different for the rows --> not same length of passed data
		sys.exit( "ERROR:     Lengths of arrays of observable don't match." )
	

	### meaning
	# calculate covariance matrix : Y'*Y/N (biased) OR Y'*Y/(N-1) (unbiased) 
	# centered data matrix with columns as observables (2 or 3) and rows as observations i .. N
	mean        = matrix.mean(axis=1)	# mean along rows
	matrix      = (matrix.T - mean).T 	# remove mean from rows with tiny hack
	cova_matrix = np.cov( matrix, bias=False )
	
	return cova_matrix
def ppol_calc(comp_1, comp_2, comp_Z=[], fix_angles='EQ', Xoffset_in_samples_for_amplitude=None):

	"""
	Take data arrays and perform 2_d and 3-D principle component
	analysis (2-D and 3-D PCA).

	Useful for seismology, for example.

		^  
		| comp_1
		|
		|
		|
		‎⊙ -----> comp_2
	  comp_z


	180° ambuigity:
	---------------
	
	  P-waves  -->  2-D data  -->  BAZ unknown  -->  no 180° ambiguity (first motion known)
               -->  3-D data  -->  BAZ unknown  -->  no 180° ambiguity (first motion known (also for sources above surface), incidence)
      R-waves  -->  3-D data  -->  BAZ unknown  -->  no 180° ambiguity (retrograde)

	Note: an unknown event BAZ is the same as an unknown station orientation
	      with a known event BAZ.
	"""


	### assign input arrays to numpy 
	# (do not demean here, this is important for option `EQ_like=False`, 
	# however you have to pass demeaned data that were demeand on a larger window beforehand!
	# Otherwise amplitude signs that option `EQ_like=False` need may be worn)
	comp_1 = np.array( comp_1 )
	comp_2 = np.array( comp_2 )
	if np.array( comp_Z ).any():
		comp_Z = np.array( comp_Z )


	### 2-D, horizontal plane
	covariance_matrix_2D_hori        = covariance_matrix(comp_1, comp_2)											# does demeaning internally!
	eig_val_2D_hori, eig_vec_2D_hori = np.linalg.eig(covariance_matrix_2D_hori)
	index_array_descending_eig_val   = eig_val_2D_hori.argsort()[::-1]
	eig_val_2D_hori                  = eig_val_2D_hori[index_array_descending_eig_val]								# E-values descending
	eig_vec_2D_hori                  = eig_vec_2D_hori[:,index_array_descending_eig_val]						    # E-vectors sorted acc. to E-values
	values                           = eig_val_2D_hori
	vectors                          = eig_vec_2D_hori
	eig_vec_2D_hori_1                = eig_vec_2D_hori[:,0]

	# derived
	phi_2D                           = (np.arctan2( eig_vec_2D_hori_1[1], eig_vec_2D_hori_1[0] ) * 180/M.pi )%360
	err_phi_2D                       = np.arctan( np.sqrt( eig_val_2D_hori[1]/eig_val_2D_hori[0] )) * 180/M.pi 	    # ...
	SNR_hori                         = (eig_val_2D_hori[0] - eig_val_2D_hori[1]) / eig_val_2D_hori[1]				# De Meersman et al. (2006)
	Rect_H                           = 1 - eig_val_2D_hori[1]/eig_val_2D_hori[0]									# rectilinearity in horizontal plane (1 for linearised, 0 for circular polarisation). Jurkevics (1988)


	if not np.array( comp_Z ).any():
		if fix_angles=='EQ':
			pass


		elif fix_angles=='AMP':
			amp_1  = comp_1[Xoffset_in_samples_for_amplitude:][np.argmax(np.abs(comp_1[Xoffset_in_samples_for_amplitude:]))]
			amp_2  = comp_2[Xoffset_in_samples_for_amplitude:][np.argmax(np.abs(comp_2[Xoffset_in_samples_for_amplitude:]))]
			
			sign_1 = np.sign( amp_1 )
			sign_2 = np.sign( amp_2 )

			if abs(amp_1)>=abs(amp_2):
				if sign_1==1:
					if phi_2D>=90 and phi_2D<=270:
						phi_2D = (phi_2D+180)%360
				else:
					if phi_2D<=90 or phi_2D>=270:
						phi_2D = (phi_2D-180)%360
			else:
				if sign_2==1:
					phi_2D = phi_2D%180
				else:
					phi_2D = (phi_2D%180)+180


		else:
			pass

		comp_R, comp_T                   = rotate.rotate_ne_rt(comp_1, comp_2, phi_2D)
		phi_3D                           = None
		err_phi_3D                       = None
		INCapp_2D                        = None
		err_INCapp_2D                    = None
		INCapp_3D                        = None
		err_INCapp_3D                    = None
		SNR_2D_radZ                      = None
		Rect_3D                          = None
		Rect_RZ                          = None
		#ccon                             = None	



	elif np.array( comp_Z ).any():

		### 3-D
		covariance_matrix_3D             = covariance_matrix(comp_1, comp_2, comp_Z)
		eig_val_3D, eig_vec_3D           = np.linalg.eig(covariance_matrix_3D)
		index_array_descending_eig_val   = np.argsort( eig_val_3D )[::-1]
		eig_val_3D                       = eig_val_3D[index_array_descending_eig_val]								# E-values descending
		eig_vec_3D                       = eig_vec_3D[:,index_array_descending_eig_val]								# E-vectors sorted acc. to E-values
		values                           = eig_val_3D
		vectors                          = eig_vec_3D
		eig_vec_3D_1                     = eig_vec_3D[:,0]

		# derived
		phi_3D                           = (np.arctan2( eig_vec_3D_1[1], eig_vec_3D_1[0] ) * 180/M.pi) % 360
		err_phi_3D                       = np.arctan( np.sqrt( eig_val_3D[2]/(eig_val_3D[1]+eig_val_3D[0]) )) * 180/M.pi 		# ...
		Rect_3D                          = 1 - ( eig_val_3D[1]+eig_val_3D[2] )/( 2*eig_val_3D[0] )					# rectilinearity in 3-D. (1 for linearised, 0 for circular polarisation). Jurkevics (1988)
		#ccon                             = eig_val_3D[0] / ( eig_val_3D[1]+eig_val_3D[2] )


		### 3-D phi, radial & Z-comp plane
		comp_R, comp_T                   = rotate.rotate_ne_rt(comp_1, comp_2, phi_3D)
		covariance_matrix_2D_radZ        = covariance_matrix(comp_Z, comp_R)
		eig_val_2D_radZ, eig_vec_2D_radZ = np.linalg.eig(covariance_matrix_2D_radZ)
		index_array_descending_eig_val   = np.argsort( eig_val_2D_radZ )[::-1]
		eig_val_2D_radZ                  = eig_val_2D_radZ[index_array_descending_eig_val]							# E-values descending
		eig_vec_2D_radZ                  = eig_vec_2D_radZ[:,index_array_descending_eig_val]						# E-vectors sorted acc. to E-values
		eig_vec_2D_radZ_1                = eig_vec_2D_radZ[:,0]

		# derived
		INCapp_3D                        = np.arctan( eig_vec_2D_radZ_1[1] / eig_vec_2D_radZ_1[0] ) * 180/M.pi
		err_INCapp_3D                    = np.arctan( np.sqrt( eig_val_2D_radZ[1]/eig_val_2D_radZ[0] )) * 180/M.pi


		### 2-D phi, radial & Z-comp plane
		comp_R, comp_T                   = rotate.rotate_ne_rt(comp_1, comp_2, phi_2D)
		covariance_matrix_2D_radZ        = covariance_matrix(comp_Z, comp_R)
		eig_val_2D_radZ, eig_vec_2D_radZ = np.linalg.eig(covariance_matrix_2D_radZ)
		index_array_descending_eig_val   = np.argsort( eig_val_2D_radZ )[::-1]
		eig_val_2D_radZ                  = eig_val_2D_radZ[index_array_descending_eig_val]							# E-values descending
		eig_vec_2D_radZ                  = eig_vec_2D_radZ[:,index_array_descending_eig_val]						# E-vectors sorted acc. to E-values
		eig_vec_2D_radZ_1                = eig_vec_2D_radZ[:,0]

		# derived
		INCapp_2D                        = np.arctan( eig_vec_2D_radZ_1[1] / eig_vec_2D_radZ_1[0] ) * 180/M.pi
		err_INCapp_2D                    = np.arctan( np.sqrt( eig_val_2D_radZ[1]/eig_val_2D_radZ[0] )) * 180/M.pi
		SNR_2D_radZ                      = (eig_val_2D_radZ[0] - eig_val_2D_radZ[1]) / eig_val_2D_radZ[1]			# De Meersman et al. (2006)
		Rect_RZ                          = 1 - eig_val_2D_radZ[1]/eig_val_2D_radZ[0]								# rectilinearity in radial-vertical plane (1 for linearised, 0 for circular polarisation). Jurkevics (1988)


		### Determined baz is ambigious by +- 180 deg. But correct baz must deliver incidence angle>0, therefore can solve ambiguity
		if fix_angles=='EQ':
			if INCapp_2D < 0:
				phi_2D         = (phi_2D+180)%360
				comp_R, comp_T = rotate_2D(comp_R, comp_T, 180)
				INCapp_2D      = abs(INCapp_2D)
			if INCapp_3D < 0:
				phi_3D         = (phi_3D+180)%360
				INCapp_3D      = abs(INCapp_3D)


		elif fix_angles=='AMP':
			amp_Z  = comp_Z[Xoffset_in_samples_for_amplitude:][np.argmax(np.abs(comp_Z[Xoffset_in_samples_for_amplitude:]))]
			amp_1  = comp_1[Xoffset_in_samples_for_amplitude:][np.argmax(np.abs(comp_1[Xoffset_in_samples_for_amplitude:]))]
			amp_2  = comp_2[Xoffset_in_samples_for_amplitude:][np.argmax(np.abs(comp_2[Xoffset_in_samples_for_amplitude:]))]
			
			sign_Z = np.sign( amp_Z )
			sign_1 = np.sign( amp_1 )
			sign_2 = np.sign( amp_2 )

			if abs(amp_1)>=abs(amp_2):
				if sign_1==1:
					if phi_2D>=90 and phi_2D<=270:
						phi_2D = (phi_2D+180)%360
				else:
					if phi_2D<=90 or phi_2D>=270:
						phi_2D = (phi_2D-180)%360
			else:
				if sign_2==1:
					phi_2D = phi_2D%180
				else:
					phi_2D = (phi_2D%180)+180

			if abs(amp_1)>=abs(amp_2):
				if sign_1==1:
					if phi_3D>=90 and phi_3D<=270:
						phi_3D = (phi_3D+180)%360
				else:
					if phi_3D<=90 or phi_3D>=270:
						phi_3D = (phi_3D-180)%360
			else:
				if sign_2==1:
					phi_3D = phi_3D%180
				else:
					phi_3D = (phi_3D%180)+180

			INCapp_2D = (abs(INCapp_2D)*sign_Z)%180 # incidence angle 0..180°
			INCapp_3D = (abs(INCapp_3D)*sign_Z)%180 # incidence angle 0..180°
		

		else:
			pass
	
	return phi_2D, phi_3D, err_phi_2D, err_phi_3D, INCapp_2D, err_INCapp_2D, INCapp_3D, err_INCapp_3D, SNR_hori, SNR_2D_radZ, Rect_3D, Rect_H, Rect_RZ, comp_R, comp_T, vectors, values
def ppol_moving(stream, window=10, step_in_samples=1, output_polar='', print_info=True):

	"""
	stream must be sorted: Z, N, E
	"""


	### ENTRY OF SCRIPT
	stream.trim_common()
	stream_string = stream[0].id[:-1] + '?'


	### EXTRACT VARIABLES NEEDED FROM STREAM OBJECT
	starttime      = stream[0].stats.starttime
	endtime        = stream[0].stats.endtime
	duration       = endtime - starttime
	delta          = stream[0].stats.delta
	SPS            = stream[0].stats.sampling_rate
	samples        = int( (endtime-starttime)/delta+1)
	data_len       = len(stream[0])
	
	window_step    = step_in_samples*delta
	window_overlap = window-window_step
	window_counter = 0
	

	### LOOP (MOVING WINDOWS)
	pol_values     = []
	line_formatter = '%6.2f %6.2f %5.2f %5.2f %6.2f %5.2f %6.2f %5.2f %7.1f %7.1f %5.3f %5.3f %5.3f %23s %23s'

	while True:
	
		# PREPARE DATA
		index_start = window_counter*step_in_samples
		index_end   = window_counter*step_in_samples+int(window*SPS)
		if index_end>=data_len:
			break	
		start  = starttime + index_start/SPS
		end    = start + window
		
		#data_Z = stream.select(component='[Z]'  )[0].data[index_start:index_end]
		#data_N = stream.select(component='[N1X]')[0].data[index_start:index_end]
		#data_E = stream.select(component='[E2Y]')[0].data[index_start:index_end]
		data_Z = stream[0].data[index_start:index_end]
		data_N = stream[1].data[index_start:index_end]
		data_E = stream[2].data[index_start:index_end]		
		

		# P-POL RESULT FOR WINDOW
		window_results = ppol_calc(data_N, data_E, data_Z)[:13] + (start, end)
		pol_values.append( line_formatter % window_results )
		
		# small output into shell
		window_counter += 1
		if window_counter%100==0:
			print(u'Window %7d:  %s - %s' % (window_counter, start, end))


	# OUTPUT FILE, IF WISHED
	if output_polar:
		header  = 'Polarizations  %s  (DATA_RANGE:  %s  %s)\n' % (stream_string, starttime, endtime)
		header += 'PH2D   PH3D 2DERR 3DERR  INC2D 2DERR  INC3D 3DERR SNR_HOR  SNR_RZ POL3D  POLH POLRZ                       START                         END'
		np.savetxt(output_polar, pol_values, header=header, fmt='%s', delimiter=None)


	# OUTPUT SHELL, IF WISHED
	if print_info:
		print()
		print(u'- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')
		print(u'Data start:       %s'          % starttime )
		print(u'Data end:         %s'          % endtime )
		print(u'Data length:      %s (h:m:s)'  % sec2hms(duration) )
		print(u'Data SPS:         %sHz (%ss)'  % (SPS, delta) )
		print(u'Data samples      %s'          % samples )
		print()
		print(u'Window length:    %ss'         % window)
		print(u'Window overlap:   %ss'         % window_overlap)
		print(u'Window step:      %ss'         % window_step)
		print(u'Number runs:      %s'          % window_counter)
		print(u'- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')


	# RETURN
	return pol_values
def xcorr(data1, data2, mode='valid'):

	"""
	Cross-correlate two data arrays using 'numpy.correlate'. 
	For this, 'data2' is shifted along 'data1'. It is not of importance
	if 'data1' contains more or less datapoints than 'data2', this script
	works for both scenarios (and equally long data arrays as well ..).
	
	The correct C-requirements are set for 'data1' & 'data2', using 
	'numpy.require'.

	Note also, amplitude scaling of data arrays don't matter, i.e.,
	don't bother normalising before exectuing this script.
	  see: https://stackoverflow.com/questions/53436231/normalized-cross-correlation-in-python

	To do the cross-correlation, choose between 3 different modes:


		VALID:
			- signals overlap fully, output dimension is thus:
			  max(A, B) - min(A, B) + 1

			#    A A A A A A
			#    B B B        = 0, index of cross-corr value in cross-corr array
			#      B B B      = 1, ..
			#        B B B    = 2, ..
			#          B B B  = 3, ..


		SAME:
			- create cross-correlation array with output length: 
			  max(A, B)
			- boundary effects may still visible
			
			#         A A A A
			#     B B B B         = 0, ..
			#       B B B B       = 1, ..
			#         B B B B     = 2, ..
			#           B B B B   = 3, ..
			

		FULL:
			- convolution at each point of overlap, output dimension:
			  (A + B - 1)

			#         A A 
			#     B B B           = 0, ..
			#       B B B         = 1, ..
			#         B B B       = 2, ..
			#           B B B     = 3, ..


	See source:
	  http://docs.scipy.org/doc/numpy/reference/generated/numpy.correlate.html


	:type data1 & data2: either list, tuple or numpy.narray
	:param data1: data1
	:param data2: data2, 'slided' along data1

	:type mode: str
	:param mode: mode of cross-correlation ('valid', 'same' or 'full')
	             default: 'valid'


	:type sample_shift: int
	:param sample_shift: number of samples to shift data2 with respect to 
	                     data1 to achieve highest correlation
	             - if = 0  :  both data arrays align with their beginning
	             - if < 0  :  data2 shifted n-samples to 'left'
	             - if > 0  :  data2 shifted n-samples to 'right'


	:type cc_factor: float, -1..1
	:param cc_factor: highest value in cross-correlation array

	:type cc_factor_abs: float, 0..1
	:param cc_factor_abs: highest value in absolute cross-correlation array

	:type corr_array: numpy.narray
	:param corr_array: each value is cross-correlation of 'data1' & 'data2'
	                   for a certain offset/sample_shift/tau between them

	:return: sample_shift, cc_factor, cc_factor_abs, corr_array
	"""


    # check variables (if they are ok to be correlated via numpy.correlate')
	if not isinstance(data1, (list, tuple, np.ndarray)) or not isinstance(data2, (list, tuple, np.ndarray)):
		print(u"ERROR:        data types neither 'list' nor 'tuple'")
		print(u"              nor 'numpy.narray'")
		return -12345, False, False, False

	if not isinstance(mode, str) or not mode.lower() in ['valid', 'same', 'full']:
		print(u"WARNING:      mode = '%s' is unvalid. Setting to: 'valid'" % mode) 
		mode = 'valid'
	mode = mode.lower()

	if len(data1) == 0 or len(data2) == 0:
		print(u'ERROR:        your data are empty.')
		return -12345, False, False, False
	
	
	# set correct flags for running C-code via numpy
	data1 = np.require(np.array( data1 ), dtype=np.float32, requirements=['C'])
	data2 = np.require(np.array( data2 ), dtype=np.float32, requirements=['C'])


    # cross-correlation, data2 'slided' along data1
	corr_array = np.correlate(data1, data2, mode=mode)


	# check if result is trivial
	if not np.any(corr_array):		# if all entries 0
		print(u'WARNING:       each single cross-correlation gave 0.')
		print(u"               'sample_shift' is thus 0, to be taken")
		print(u'               with a grain of salt !')


	# since script shall work for both scenarios - len(data1) > len(data2)
	# and vice versa - this is how it does ..
	if len(data1) < len(data2):
		corr_array = corr_array[::-1]


	# 'sample_shift' calculation, using 'corr_array' (not abs(corr_array) !)
	if mode == 'valid':
		sample_shift = np.argmax(corr_array)
	elif mode == 'same':
		sample_shift = np.argmax(corr_array) - min( len(data1),len(data2) ) / 2
	elif mode == 'full':
		sample_shift = np.argmax(corr_array) - min( len(data1),len(data2) ) + 1


	# since script shall work for both scenarios - len(data1) > len(data2)
	# and vice versa - this is how it does ..
	if len(data1) < len(data2):
		sample_shift *= -1
		
	
    # normalisation
	normalisation_factor = np.sqrt( (data1 ** 2).sum() * (data2 ** 2).sum() )
	if normalisation_factor == 0:
		print(u'ERROR:        normalisation_factor = 0! Most likely, one of')
		print(u'              the passed data arrays consists of only zeros!')
		return -12345, False, False, False


    # highest normalised (absolute) correlation value 
    # -1 <= cc_factor     <= 1
    #  0 <= cc_factor_abs <= 1
	cc_factor     = np.max(     corr_array  ) / normalisation_factor
	cc_factor_abs = np.max( abs(corr_array) ) / normalisation_factor

	return (sample_shift, cc_factor, cc_factor_abs, corr_array)
def spectrum(data, method='fft', fs=1.0, nfft=None, window_welch='hann', axis=-1):

	"""
	Plot absolute power spectrum of signal, using Scipy functions.
	Choose method between `fft` and `welch`.
	(Welch method splits input into segments, calculates
	spectral density via periodigrams, and then aerages of segments).

	Returns: freqs, ps


	For Welch method, possible windows are:
	  - boxcar, triang, blackman, hamming, hann, bartlett, flattop, parzen, bohman, blackmanharris, nuttall, 
	    barthann, kaiser (needs beta), gaussian (needs standard deviation), general_gaussian (needs power, width), 
	    slepian (needs width), dpss (needs normalized half-bandwidth), chebwin (needs attenuation), 
	    exponential (needs decay scale), tukey (needs taper fraction)

	    see: https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.get_window.html#scipy.signal.get_window


	For details:
	  - https://docs.scipy.org/doc/scipy/reference/fftpack.html
	  - https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.welch.html#scipy.signal.welch
	"""

	data = np.array( data )

	if method.lower() == 'fft':
		if nfft:
			window_length = nfft
		else:
			window_length = len(data)
		freqs = scipy.fftpack.fftfreq(window_length, 1/fs)
		ps    = pow( np.abs( scipy.fftpack.fft(data, n=nfft, axis=axis) ), 2)


	elif method.lower() == 'welch':
		freqs, ps = scipy.signal.welch(data, fs=fs, nfft=nfft, window=window_welch, scaling='spectrum', nperseg=None, noverlap=None, detrend='constant', return_onesided=True)

	return freqs, ps
def angle_vectors(vec_1, vec_2):

	"""
	The np.clip() application takes care on rounding problems.
	"""

	v1_u = unit_vector(vec_1)
	v2_u = unit_vector(vec_2)

	return np.arccos( np.clip( np.dot(v1_u, v2_u), -1.0,1.0)) * 180/M.pi
def unit_vector(vec):

    """ 
    Returns the unit vector of the vector.  
    """

    return vec / np.linalg.norm(vec)
def nearest_value_in_array(array, value):

	"""
	Return the one value from an 'array' 
	that is closest to passed 'value'.
	"""

	idx = ( np.abs( np.array(array)-value) ).argmin()
	return array[idx]
def next_power_of_2(x):

	"""
	Find next highest power of 2 number of 'x'.
	"""
	
	return 1 if x == 0 else 2**(x - 1).bit_length()
def sin_gen(freq, full_periods=1, delta=0.01, phase_shift_deg=0, amp=1, plot=False):
	
	"""
	Return x, y (=numpy arrays) of a sin function of frequency `f` and
	`full_periods` long. You can choose the `delta` (distance of two 
	points in the defition area), the`phase_shift_deg` and the 
	amplitude `amp`.
	"""

	x   = np.linspace( 0, full_periods/freq, full_periods/(freq*delta)+1, endpoint=True )
	sin = amp * np.sin(2*np.pi * (x*freq + phase_shift_deg/360))
	
	if plot:
		fig = plt.figure()
		ax  = fig.add_subplot(111)
		ax.set(xlabel='t (s)', ylabel='sin(t)')
		ax.grid(ls='-.', lw=0.5)
		ax.plot(x, sin, lw=2)
		ax.text(0.98, 0.05, "T=%.3f s,  f=%.3f Hz" % (1./freq, freq), ha='right', color='k', transform=ax.transAxes, fontsize=6)
		plt.show()

	return x, sin
def fit_sin(data, guess_freq=1):

	'''
	Fit sin to the input time sequence, and return fitting parameters "amp", "omega", 
	"phase", "offset", "freq", "period" and "fitfunc"
	'''

	data         = np.array(data)
	guess_amp    = np.std(data) * 2.**0.5
	guess_offset = np.mean(data)
	guess        = np.array( [guess_amp, 2.*np.pi*guess_freq, 0., guess_offset] )

	def sinfunc(t, A, w, p, c): 
		return A * np.sin(w*t + p) + c
	popt, pcov = scipy.optimize.curve_fit(sinfunc, np.arange(len(data)), data, p0=guess, maxfev=100000)
	A, w, p, c = popt
	f          = w/(2.*np.pi)
	fitfunc    = lambda t: A * np.sin(w*t + p) + c
	
	return {"amp": A, "omega": w, "phase": p, "offset": c, "freq": f, "period": 1./f, "fitfunc": fitfunc, "maxcov": np.max(pcov), "rawres": (guess,popt,pcov)}
def snr(data, axis=0, ddof=1):

	"""
	Signal-to-noise ratio (SNR), as by Scipy.
	"""

	data = np.array(data)
	mean = data.mean(axis)
	sd_d = data.std(axis=axis, ddof=ddof)
	return np.where(sd_d==0, 0, mean/sd_d)
def wavelet_transform(data, method='swt', wavelet='haar', mode='symmetric', level=None, set_to_zero='all', threshold=100, nmmps=50, Tmax=None, freq_ny=1., smooth_denoise=False, print_parameters=False, print_Cs=False, print_xcorr=True, plot_compare=False, plot_Cs=False, title=''):

	"""
	Use Python libraries "waveletpy" (for discrete and stationary 
	wavelet transform) and "dwctw" (for dual tree complex wavelet 
	transform) to calculate wavelet transform data of passed data (=numpy array).

	Import:
	-------
	  import pywt (for waveletpy)
	  import dtcwt (for dtcwt)


	NOTES:
	------
	  - data are padded as to `mode` before any analysis
	  - 'method' : either 'dwt', 'swt', or 'dtcwt'
	  - for possible (discrete) 'wavelet' 
	    functions for method='dwt' & 'swt', see:
		   pywt.wavelist(kind='discrete')
		for method='dtcwt', the TWO 'wavelet' functions
		are hardcoded but can be changed (see in the script, by now 
		they are 'near_sym_a' & 'qshift_a'). Details:
		   https://dtcwt.readthedocs.io/en/0.11.0/reference.html#dtcwt.coeffs.biort
		   https://dtcwt.readthedocs.io/en/0.11.0/reference.html#dtcwt.coeffs.qshift
	  - 'mode' determines used padding pattern (to next power of 2 of data length).
	    details:
	       https://pywavelets.readthedocs.io/en/latest/ref/signal-extension-modes.html#ref-modes
	    This parameter is deprecated for method 'dtcwt'.
	  - if 'level' = None:
		   dwt   : decompose at maximum level (determined internally)
		   swt   : decompose at maximum level (determined internally)
		  (for those 2, if level is manually set too high, error is thrown)
		   dtcwt : decompose at maximum level of dwt (which uses data length and 'wavelet')
	  - set_to_zero='all'   : for all levels, set detailed coefficients
	                          to 0, if their absolute value < threshold
	               ='level' : for deepest level, set detailed coefficients
	              			  to 0, if their absolute value < threshold
	                          (detailed coefficients of other levels are left untouched)
	  - smooth_denoise : True or False. If True (and threshold != 0), this happens:
	                       for detailed coefficients (cD) < 0 : cD += threshold
	                       for detailed coefficients (cD) > 0 : cD -= threshold
	                     This, of course, applies only to the deepest level or all levels,
	                     depending on variable 'level'.
	"""


	now         = time.time()
	set_to_zero = set_to_zero.lower()


	###  PREPARE  DATA  FOR  WAVELET ANALYSIS  (padding to next power of 2 of data length)
	if len(data)%2 != 0:
		# so padding can be applied to both sides without problems
		data = data[:-1]


	pad_width = 0
	if len(data) < next_power_of_2( len(data) ):
		# padding same number of elements left & right, depending on `mode`
		pad_width = int( (next_power_of_2( len(data) ) - len(data)) / 2)	# pad same number of elements left & right 
		if mode == 'symmetric':
			data = np.pad(data, pad_width, mode, reflect_type='even')
		elif mode == 'reflect':
			data = np.pad(data, pad_width, mode)
		elif mode == 'smooth' or mode == 'periodic' or mode == 'periodization':
			print('Mode "%s" as in waveletpy not implemnted in numpy, should be hand crafted.' % mode)
			print('Continue with mode SYMMETRIC.')
			mode = 'symmetric'
			data = np.pad(data, pad_width, mode)
		elif mode == 'constant':
			data = np.pad(data, pad_width, 'edge')
		elif mode == 'periodic':
			data = np.pad(data, pad_width, 'wrap')
		elif mode == 'antireflect':
			data = np.pad(data, pad_width, 'reflect', reflect_type='odd')
		else:
			print('Mode "%s" not a valid option. See online for reference. Or check first comments in script.' % mode)
			print('Continue with mode SYMMETRIC.')
			mode = 'symmetric'
			data = np.pad(data, pad_width, mode)


	###  DISCRETE  WAVELET  TRANSFORM  (DWT)
	if method.lower() == 'dwt':

		coefficients   = pywt.wavedec(data, wavelet, level=level, mode=mode)		#[cA_n, cD_n, cD_n-1, ..., cD2, cD1]
		if not level:
			level      = pywt.dwt_max_level( len(data), wavelet)

		if print_Cs:
			printer_cAs = []
			cA          = np.abs(coefficients[0])
			printer_cAs.append( [min(cA), np.median(cA), np.average(cA), max(cA), len(cA)] )

			printer_cDs  = []
			for cD in coefficients[1:]:
				cD = np.abs(cD)
				printer_cDs.append( [min(cD), np.median(cD), np.average(cD), max(cD), len(cD)] )


		if set_to_zero:
			j = 0
			maximum = []
			cA      = coefficients[0]
			for cD in coefficients[1:]:


				if set_to_zero == 'ft':
					# FIXED THRESHOLD (FT) APPROACH (absolute values > threshold --> 0)
					cD[ np.abs(cD) > threshold ] = 0
					quick_plot(copy_cD, cD, np.ones(len(cD))*threshold, -np.ones(len(cD))*threshold, data_labels=('ori','new','threshold','-threshold'), verts=verts, verts2=verts2)		
				elif set_to_zero == 'ftr':
					# FIXED THRESHOLD (FT) ROOT APPROACH (if one absolute value between 2 roots > threshold --> 0 all values between 2 roots)
					zero_crossings = np.where( np.diff( np.signbit(cD) ))[0]
					zero_crossings = np.insert(zero_crossings, len(zero_crossings), len(cD)-1)
					zero_crossings = np.insert(zero_crossings, 0, 0)
					for k in range(len(zero_crossings[:-1])):
						try:
							root_data = np.abs( cD[ zero_crossings[k]:zero_crossings[k+1]+1] )
							if max( root_data ) > threshold:
								cD[ zero_crossings[k]:zero_crossings[k+1]+1 ] = 0
						except ValueError:
							continue
					quick_plot(copy_cD, cD, np.ones(len(cD))*threshold, -np.ones(len(cD))*threshold, data_labels=('ori','new','threshold','-threshold'), verts=verts, verts2=verts2)		
				elif set_to_zero == 'malm':
					# MEDIAN ABSOLUTE LOCAL MAXIMUM (MALM) APPROACH (absolute values > x*MALM --> 0)
					zero_crossings = np.where( np.diff( np.signbit(cD) ))[0]
					zero_crossings = np.insert(zero_crossings, len(zero_crossings), len(cD)-1)
					zero_crossings = np.insert(zero_crossings, 0, 0)
					alm            = []			# asbolute local maxima
					for k in range(len(zero_crossings[:-1])):
						try:
							root_data = np.abs( cD[ zero_crossings[k]:zero_crossings[k+1]+1] )
							alm.append( max(root_data) )
						except Exception:
							continue
					_threshold_                  = threshold*np.median(alm)
					cD[ np.abs(cD) > _threshold_ ] = 0
				elif set_to_zero == 'malmr':
					#maximum.append( cD[25000] ) 
					# MEDIAN ABSOLUTE LOCAL MAXIMUM (MALM) ROOT APPROACH (if one absolute value between 2 roots > x*MALM --> 0 all values between 2 roots)
					zero_crossings = np.where( np.diff( np.signbit(cD) ))[0]
					zero_crossings = np.insert(zero_crossings, len(zero_crossings), len(cD)-1)
					zero_crossings = np.insert(zero_crossings, 0, 0)
					alm            = []
					for k in range(len(zero_crossings[:-1])):
						try:
							root_data = np.abs( cD[ zero_crossings[k]:zero_crossings[k+1]+1] )
							alm.append( max(root_data) )
						except Exception:
							continue
					_threshold_ = threshold*np.median(alm)
					
					for k in range(len(zero_crossings[:-1])):
						try:
							root_data = np.abs( cD[ zero_crossings[k]:zero_crossings[k+1]+1] )
							if max( root_data ) > _threshold_:
								cD[ zero_crossings[k]:zero_crossings[k+1]+1 ] = 0
						except ValueError:
							continue				
				elif set_to_zero == 'gmr':
					# GLOBAL MAXIMUM RATIO (GMR) APPROACH (only if max abs value > x*MALM (detect if spikes at all), then absolute values > global maximum/y --> 0)
					zero_crossings = np.where( np.diff( np.signbit(cD) ))[0]
					zero_crossings = np.insert(zero_crossings, len(zero_crossings), len(cD)-1)
					zero_crossings = np.insert(zero_crossings, 0, 0)
					alm            = []
					for k in range(len(zero_crossings[:-1])):
						try:
							root_data = np.abs( cD[ zero_crossings[k]:zero_crossings[k+1]+1] )
							alm.append( max(root_data) )
						except Exception:
							continue
					proceed_threshold = 10*np.median(alm)
					max_cD            = max(np.abs(cD))
					_threshold_       = max_cD/threshold
					
					max_cD = max(np.abs(cD))
					if max_cD > proceed_threshold: 
						cD[ np.abs(cD) > _threshold_] = 0
					quick_plot(copy_cD, cD, np.ones(len(cD))*_threshold_, -np.ones(len(cD))*_threshold_, data_labels=('ori','new','threshold','-threshold'), verts=verts, verts2=verts2)		
				elif set_to_zero == 'gmrr':
					# GLOBAL MAXIMUM RATIO (GMR) ROOT APPROACH (only if max abs value > x*MALM (detect if spikes at all), then if one absolute value between 2 roots > global maximum/y --> 0 all values between 2 roots)
					zero_crossings = np.where( np.diff( np.signbit(cD) ))[0]
					zero_crossings = np.insert(zero_crossings, len(zero_crossings), len(cD)-1)
					zero_crossings = np.insert(zero_crossings, 0, 0)
					alm            = []
					for k in range(len(zero_crossings[:-1])):
						try:
							root_data = np.abs( cD[ zero_crossings[k]:zero_crossings[k+1]+1] )
							alm.append( max(root_data) )
						except Exception:
							continue
					proceed_threshold = 10*np.median(alm)
					max_cD            = max(np.abs(cD))
					_threshold_       = max_cD/threshold
					
					for k in range(len(zero_crossings[:-1])):
						try:
							root_data = np.abs( cD[ zero_crossings[k]:zero_crossings[k+1]+1] )
							if max( root_data ) > _threshold_:
								cD[ zero_crossings[k]:zero_crossings[k+1]+1 ] = 0
						except ValueError:
							continue			
				elif set_to_zero == 'mf':
					# MEDIAN FILTER (MR) APPROACH - doesn't work so well (long compute time, spikes not removed, information loss)!				
					pad_width2 = threshold
					cD_pad     = np.pad(cD, pad_width2, 'symmetric', reflect_type='even')
					cD_new     = [np.median( cD_pad[n:n+2*pad_width2+1] ) for n in range(len(cD)) ]
					for b in range(len(cD)):
						cD[b] = cD_new[b]
				elif set_to_zero == 'mane':
					# MOVING AVERAGE N ENVELOPE (MANE) APPROACH - 
					j += 1
					if j<8:
						continue
					N     = threshold
					mane  = smooth( envelope(cD), N, 'blackman')
					#print(len(mane), len(cD))
					
					for k in range(len(cD)):
						if abs(cD[k]) > mane[k]:
							cD[k] = mane[k]*np.sign(cD[k])
				elif set_to_zero == 'sin':
					## parameters
					copy_cD     = cD.copy()	# copy original coefficients, useful for plotting of individual decompositions layers only
					
					j           +=1			# level counter (just useful)
                    #threshold   = 4.5		# threshold used to calculate glitch detection condition, here threshold acts as factor (see mmalm)
                    #nmmps       = 10		# number_moving_maxima_per_side
					nagps       = 1			# number_average_glitch_per_side
					
					Tmin_j4     = 2**(level-4)   / freq_ny
					Tmin_j      = 2**(level-j)   / freq_ny
					Tmax_j      = 2**(level-j+1) / freq_ny
					Tmax_decomp = 2**level       / freq_ny


					## detrend decompostion level with only 2 and 4 samples, respectively, pad the others for better computations
					title = u'Level: %2d (len=%6d, nmmps=%2d, nagps=%2d, Tmin=%9.1fs, Tmax=%9.1fs, fmin=%9.7fHz, fmax=%9.7fHz)' % (j, len(cD), nmmps, nagps, Tmin_j, Tmax_j, 1/Tmax_j, 1/Tmin_j)
					print(title, end=' ')
					if (Tmax and Tmax<Tmin_j) or j<=4:
						print('UNCHANGED.')
						#cD[:] = scipy.signal.detrend(cD)
						#cD[:] = np.zeros(len(cD))
						continue										
					elif Tmax and (Tmax>Tmax_decomp or Tmax>=Tmin_j4):
						print()
						print(u'WARNING:  `Tmax` (%ss) is >= Tmax of decomposition (%ss).' % (Tmax, Tmax_decomp))
						print(u'Please decrease `Tmax` or increase data length, then re-try.')
						return
					else:
						print('CHANGED.')
						#cD[:]     = scipy.signal.detrend(cD)
						pad_width3 = int( (next_power_of_2(len(cD)+1) - len(cD)) / 2)	# padding decomposition level (to treat beginning & end without problems) 
						cD2        = np.pad(cD, pad_width3, 'symmetric')				# same number of elements left & right, symmetric padding					


					## determine MMALM (Moving Median of Absolute Local Maxima)
					lm  = scipy.signal.argrelextrema( abs(cD2), np.greater)[0]
					#lmMin  = scipy.signal.argrelextrema( cD2, np.less   )[0]		# NOTE to myself: algorithm works better with abs(cD) instead of getting Maxima & Minima
					#lm     = np.sort( np.concatenate(( lmMax, lmMin )))
					lm     = np.insert(lm, len(lm), len(cD2)-1)
					lm     = np.insert(lm, 0, 0)
					#print(lm)

					mmalm = []
					for p in range( len(lm) ):
						liste = np.abs( [ cD2[ lm[ h+p ] ] for h  in np.arange( -nmmps, nmmps+1 ) if h+p>=0 and h+p<=len(lm)-1])
						if j<1:
							malm = 2.0*np.median( liste )
						else:
							malm = threshold*np.median( liste )
						#print('auffuellen x mal: %s, wert = %s' % (max(1, lm[p]-lm[p-1]), malm))
						for o in range( max(1, lm[p]-lm[p-1] )):
							mmalm.append(malm)
					mmalm = np.array( mmalm )
					mmalm = np.array( [max(mmalm) if np.isnan(m) else m for m in mmalm] )


					## correct cDs
					for l in range(len(cD2)):

						# positive glitches
						if cD2[l] >= mmalm[l]:

							list_before   = []
							index_shifter = 0
							while True:
								index_shifter += 1
								index_before   = l-index_shifter
								try:
									if cD2[index_before] <= mmalm[index_before] and cD2[index_before] >= -mmalm[index_before] and cD2[index_before-1]<=cD2[index_before] and cD2[index_before]>=cD2[index_before+1]:
										list_before.append( cD2[index_before] )
								except IndexError:
									break
								if len(list_before) == nagps:
									break

							list_after    = []
							index_shifter = 0
							while True:
								index_shifter += 1
								index_after    = l+index_shifter
								try:
									if cD2[index_after] <= mmalm[index_after] and cD2[index_after] >= -mmalm[index_after] and cD2[index_after-1]<=cD2[index_after] and cD2[index_after]>=cD2[index_after+1]:
										list_after.append( cD2[index_after] )
								except IndexError:
									break
								if len(list_after) == nagps:
									break

							list_combined = list_before + list_after
							cD2[l]        = np.median( np.array( list_combined ))


						# negative glitches
						if cD2[l] <= -mmalm[l]:

							list_before   = []
							index_shifter = 0
							while True:
								index_shifter += 1
								index_before   = l-index_shifter
								try:
									if cD2[index_before] >= -mmalm[index_before] and cD2[index_before] <= mmalm[index_before] and cD2[index_before-1]>=cD2[index_before] and cD2[index_before]<=cD2[index_before+1]:
										list_before.append( cD2[index_before] )
								except IndexError:
									break
								if len(list_before) == nagps:
									break

							list_after    = []
							index_shifter = 0
							while True:
								index_shifter += 1
								index_after    = l+index_shifter
								try:
									if cD2[index_after] >= -mmalm[index_after] and cD2[index_after] <= mmalm[index_after] and cD2[index_after-1]>=cD2[index_after] and cD2[index_after]<=cD2[index_after+1]:
										list_after.append( cD2[index_after] )
								except IndexError:
									break
								if len(list_after) == nagps:
									break

							list_combined = list_before + list_after
							cD2[l]        = np.median( np.array( list_combined ))


					## Undo padding for individual decomposition layers & mmalm (Moving Median of Absolute Local Maxima)
					cD[:] = cD2[pad_width3:len(cD2)-pad_width3]
					mmalm = mmalm[pad_width3:len(mmalm)-pad_width3]


					## Plotting individual decompistion levels
					#quick_plot(copy_cD, cD, mmalm, -mmalm, data_labels=('ori','new','mmalm','-mmalm'), title=title, xlabel='Samples of decomposed function', ylabel='Decomposition coefficient', verts=[len(cD)*0.618])
				elif len(set_to_zero.split('<'))>1 or len(set_to_zero.split('>'))>1 or len(set_to_zero.split('='))>1:
					# only decomposition levels threshold checked fulfill a condition 
					j += 1
					condition = str(j)+set_to_zero
					if eval(condition):
						cD[np.abs(cD) > threshold] = 0
				else:
					print('Option "%s" not implemented to set certain decompoistion coefficients to 0. Skipped!' % set_to_zero)
					break

				if smooth_denoise and threshold != 0:
					cD[cD > 0] -= threshold
					cD[cD < 0] += threshold

		else:
			# nothing
			print('No coefficients were set to zero!')	


		if print_Cs:
			printer_cAs_new = []
			cA              = np.abs(coefficients[0])
			printer_cAs_new.append( [min(cA), np.median(cA), np.average(cA), max(cA), len(cA)] )

			printer_cDs_new = []
			for cD in coefficients[1:]:
				cD = np.abs(cD)
				printer_cDs_new.append( [min(cD), np.median(cD), np.average(cD), max(cD), len(cD)] )


		plot_param = np.array(coefficients)[0], np.array(coefficients)[1:]
		data_recon = pywt.waverec(coefficients, wavelet, mode=mode)


	###  STATIONARY  WAVELET  TRANSFORM  (SWT)
	elif method.lower() == 'swt':

		# Approximate coefficients have no influence on reconstruction (= redundant!)
		coefficients = pywt.swt(data, wavelet, level=level)[::-1]		#coefficients = [(cA1, cD1), (cA2, cD2), ..., (cAn, cDn)]
		if not level:
			# if no level set, determine max level that was used (by default)
			level    = pywt.swt_max_level( len(data) )


		if print_Cs:
			printer_cAs  = []
			for cA_cD in coefficients:
				cA = np.abs(cA_cD[0])
				printer_cAs.append( [min(cA), np.median(cA), np.average(cA), max(cA), len(cA)] )

			printer_cDs  = []
			for cA_cD in coefficients:
				cD = np.abs(cA_cD[1])
				printer_cDs.append( [min(cD), np.median(cD), np.average(cD), max(cD), len(cD)] )


		if set_to_zero:
			j = 0
			maximum = []
			for cA, cD in coefficients:
				copy_cD = np.copy(cD)

				if set_to_zero == 'ft':
					# FIXED THRESHOLD (FT) APPROACH (absolute values > threshold --> 0)
					cD[ np.abs(cD) > threshold ] = 0			
				elif set_to_zero == 'ftr':
					# FIXED THRESHOLD (FT) ROOT APPROACH (if one absolute value between 2 roots > threshold --> 0 all values between 2 roots)
					zero_crossings = np.where( np.diff( np.signbit(cD) ))[0]
					zero_crossings = np.insert(zero_crossings, len(zero_crossings), len(cD)-1)
					zero_crossings = np.insert(zero_crossings, 0, 0)
					for k in range(len(zero_crossings[:-1])):
						try:
							root_data = np.abs( cD[ zero_crossings[k]:zero_crossings[k+1]+1] )
							if max( root_data ) > threshold:
								cD[ zero_crossings[k]:zero_crossings[k+1]+1 ] = 0
						except ValueError:
							continue
				elif set_to_zero == 'malm':
					# MEDIAN ABSOLUTE LOCAL MAXIMUM (MALM) APPROACH (absolute values > x*MALM --> 0)
					zero_crossings = np.where( np.diff( np.signbit(cD) ))[0]
					zero_crossings = np.insert(zero_crossings, len(zero_crossings), len(cD)-1)
					zero_crossings = np.insert(zero_crossings, 0, 0)
					alm            = []			# asbolute local maxima
					for k in range(len(zero_crossings[:-1])):
						try:
							root_data = np.abs( cD[ zero_crossings[k]:zero_crossings[k+1]+1] )
							alm.append( max(root_data) )
						except Exception:
							continue
					_threshold_                  = threshold*np.median(alm)
					cD[ np.abs(cD) > _threshold_ ] = 0
				elif set_to_zero == 'malmr':
					maximum.append( cD[25000] ) 
					# MEDIAN ABSOLUTE LOCAL MAXIMUM (MALM) ROOT APPROACH (if one absolute value between 2 roots > x*MALM --> 0 all values between 2 roots)
					zero_crossings = np.where( np.diff( np.signbit(cD) ))[0]
					zero_crossings = np.insert(zero_crossings, len(zero_crossings), len(cD)-1)
					zero_crossings = np.insert(zero_crossings, 0, 0)
					alm            = []
					for k in range(len(zero_crossings[:-1])):
						try:
							root_data = np.abs( cD[ zero_crossings[k]:zero_crossings[k+1]+1] )
							alm.append( max(root_data) )
						except Exception:
							continue
					_threshold_ = threshold*np.median(alm)
					
					for k in range(len(zero_crossings[:-1])):
						try:
							root_data = np.abs( cD[ zero_crossings[k]:zero_crossings[k+1]+1] )
							if max( root_data ) > _threshold_:
								cD[ zero_crossings[k]:zero_crossings[k+1]+1 ] = 0
						except ValueError:
							continue				
				elif set_to_zero == 'gmr':
					# GLOBAL MAXIMUM RATIO (GMR) APPROACH (only if max abs value > x*MALM (detect if spikes at all), then absolute values > global maximum/y --> 0)
					zero_crossings = np.where( np.diff( np.signbit(cD) ))[0]
					zero_crossings = np.insert(zero_crossings, len(zero_crossings), len(cD)-1)
					zero_crossings = np.insert(zero_crossings, 0, 0)
					alm            = []
					for k in range(len(zero_crossings[:-1])):
						try:
							root_data = np.abs( cD[ zero_crossings[k]:zero_crossings[k+1]+1] )
							alm.append( max(root_data) )
						except Exception:
							continue
					proceed_threshold = 10*np.median(alm)
					max_cD            = max(np.abs(cD))
					_threshold_       = max_cD/threshold
					
					max_cD = max(np.abs(cD))
					if max_cD > proceed_threshold: 
						cD[ np.abs(cD) > _threshold_] = 0
				elif set_to_zero == 'gmrr':
					# GLOBAL MAXIMUM RATIO (GMR) ROOT APPROACH (only if max abs value > x*MALM (detect if spikes at all), then if one absolute value between 2 roots > global maximum/y --> 0 all values between 2 roots)
					zero_crossings = np.where( np.diff( np.signbit(cD) ))[0]
					zero_crossings = np.insert(zero_crossings, len(zero_crossings), len(cD)-1)
					zero_crossings = np.insert(zero_crossings, 0, 0)
					alm            = []
					for k in range(len(zero_crossings[:-1])):
						try:
							root_data = np.abs( cD[ zero_crossings[k]:zero_crossings[k+1]+1] )
							alm.append( max(root_data) )
						except Exception:
							continue
					proceed_threshold = 10*np.median(alm)
					max_cD            = max(np.abs(cD))
					_threshold_       = max_cD/threshold
					
					for k in range(len(zero_crossings[:-1])):
						try:
							root_data = np.abs( cD[ zero_crossings[k]:zero_crossings[k+1]+1] )
							if max( root_data ) > _threshold_:
								cD[ zero_crossings[k]:zero_crossings[k+1]+1 ] = 0
						except ValueError:
							continue			
				elif set_to_zero == 'mf':
					# MEDIAN FILTER (MR) APPROACH - doesn't work so well (long compute time, spikes not removed, information loss)!				
					pad_width2 = threshold
					cD_pad     = np.pad(cD, pad_width2, 'symmetric', reflect_type='even')
					cD_new     = [np.median( cD_pad[n:n+2*pad_width2+1] ) for n in range(len(cD)) ]
					for b in range(len(cD)):
						cD[b] = cD_new[b]				
				elif len(set_to_zero.split('<'))>1 or len(set_to_zero.split('>'))>1 or len(set_to_zero.split('='))>1:
					# only decomposition levels threshold checked fulfill a condition 
					j += 1
					condition = str(j)+set_to_zero
					if eval(condition):
						cD[np.abs(cD) > threshold] = 0
				else:
					print('Option "%s" not implemented to set certain decompoistion coefficients to 0. Skipped!' % set_to_zero)
					break

				if smooth_denoise and threshold != 0:
					cD[cD > 0] -= threshold
					cD[cD < 0] += threshold
				
				#quick_plot(copy_cD, cD, np.zeros(len(cD))+threshold*np.median(alm), np.zeros(len(cD))-threshold*np.median(alm))
				#sys.exit()
		else:
			# nothing
			print('No coefficients were set to zero!')			


		if print_Cs:
			printer_cAs_new  = []
			for cA_cD in coefficients:
				cA = np.abs(cA_cD[0])
				printer_cAs_new.append( [min(cA), np.median(cA), np.average(cA), max(cA), len(cA)] )

			printer_cDs_new  = []
			for cA_new, cD_new in coefficients:
				cD = np.abs(cD_new)
				printer_cDs_new.append( [min(cD), np.median(cD), np.average(cD), max(cD), len(cD)] )

		plot_param = np.array(coefficients)[:,0], np.array(coefficients)[:,1]
		#data_recon = np.sum( np.sum(coefficients, axis=0), axis=0)
		#print(maximum)
		data_recon = pywt.iswt(coefficients[::-1], wavelet)


	###  DUAL  TREE  COMPLEX  WAVELET  TRANSFORM  (DTCWT)
	elif method.lower() == 'dtcwt':

		if not level:
			level = pywt.dwt_max_level( len(data), wavelet)
		biort     = 'near_sym_a'
		qshift    = 'qshift_a'
		transform = dtcwt.Transform1d(biort=biort, qshift=qshift)
		pyramid   = transform.forward(data, nlevels=level)

		if print_cDs:
			cDs  = []
			for cD in pyramid.highpasses[::-1]:
				cD = np.abs(cD)
				cDs.append( [min(cD), np.median(cD), np.average(cD), max(cD), len(cD)] )

		if set_to_zero.lower() == 'level':
			cD_level = pyramid.highpasses[level-1]
			cD_level[np.abs(cD_level) < threshold] = 0
			if smooth_denoise and threshold != 0:
				cD_level[cD_level > 0] -= threshold
				cD_level[cD_level < 0] += threshold
		elif set_to_zero.lower() == 'all':
			for cD in pyramid.highpasses:
				cD[np.abs(cD) < threshold] = 0
				if smooth_denoise and threshold != 0:
					cD[cD > 0] -= threshold
					cD[cD < 0] += threshold
		else:
			print('Option "%s" not implemented to set certain cDs (detailed coefficients) to 0. Skipped!' % set_to_zero)

		if print_cDs:
			cDs_new  = []
			for cD_new in pyramid.highpasses[::-1]:
				cD_new = np.abs(cD_new)
				cDs_new.append( [min(cD_new), np.median(cD_new), np.average(cD_new), max(cD_new), len(cD_new)] )

		data_recon   = transform.inverse(pyramid)


	###  OTHER  METHODS  NOT  IMPLEMENTED
	else:
		print('Method "%s" not implemented. Choose between "dwt", "swt", or "dtcwt". Returned empty array!' % method)
		return np.array([])


	### UNDO DATA-PADDING
	#plot_param = plot_param[0][:,pad_width:len(plot_param[0])-pad_width], plot_param[1][:,pad_width:len(plot_param[1])-pad_width]
	
	#data       = scipy.signal.detrend( data[pad_width:len(data)-pad_width] )
	#data_recon = scipy.signal.detrend( data_recon[pad_width:len(data_recon)-pad_width] )
	data       = data[pad_width:len(data)-pad_width]
	data_recon = data_recon[pad_width:len(data_recon)-pad_width]


	### Print wavelet parameters, if wished
	if print_parameters:
		parameters                    = {}
		parameters['method']          = method
		if method == 'dtcwt':
			parameters['wavelet']     = '%s, %s' % (biort, qshift)			
		else:
			parameters['wavelet']     = wavelet
		if method == 'dtcwt':
			parameters['mode']        = 'N.A.'
		else:
			parameters['mode']        = mode
		parameters['level']           = level
		if set_to_zero == 'level':
			parameters['set_to_zero'] = 'deepest ' + set_to_zero
		else:
			parameters['set_to_zero'] = set_to_zero			
		parameters['threshold']       = threshold
		parameters['smooth_denoise']  = smooth_denoise
		print()
		print('Wavelet transform parameters:')
		for key in parameters.keys():
			print('  %20s : %s' % (key, parameters[key]))


	### Print some approximation coefficient (cA) overview, if wished 
	if print_Cs:
		print()
		print('Approximation coefficients (before zeroing coeffs: %s, %s):' % (set_to_zero, threshold))
		for i, cA in enumerate(printer_cAs):
			label = 'cA%s' % (i+1)
			print('  %20s : min=%-12.2f  med=%-12.2f  ave=%-12.2f  max=%-12.2f  len=%-d' % (label, cA[0], cA[1], cA[2], cA[3], cA[4]))
		print('Approximation coefficients (after zeroing coeffs: %s, %s):' % (set_to_zero, threshold))
		for i, cA in enumerate(printer_cAs_new):
			label = 'cA%s' % (i+1)
			print('  %20s : min=%-12.2f  med=%-12.2f  ave=%-12.2f  max=%-12.2f  len=%-d' % (label, cA[0], cA[1], cA[2], cA[3], cA[4]))

		print()
		print('Detailed coefficients (before zeroing coeffs: %s, %s):' % (set_to_zero, threshold))
		for i, cD in enumerate(printer_cDs):
			label = 'cD%s' % (i+1)
			print('  %20s : min=%-12.2f  med=%-12.2f  ave=%-12.2f  max=%-12.2f  len=%-d' % (label, cD[0], cD[1], cD[2], cD[3], cD[4]))
		print('Detailed coefficients (after zeroing coeffs: %s, %s):' % (set_to_zero, threshold))
		for i, cD in enumerate(printer_cDs_new):
			label = 'cD%s' % (i+1)
			print('  %20s : min=%-12.2f  med=%-12.2f  ave=%-12.2f  max=%-12.2f  len=%-d' % (label, cD[0], cD[1], cD[2], cD[3], cD[4]))


	### Print quality of match & run time, if wished
	if print_xcorr:
		quality = (xcorr(data, data_recon, mode= 'valid')[1]+1)/2 * 100
		print()
		print('Quality of match between data and %s-reconstructed data:' % method.upper())
		print('  %20s : %-.2f %%'                                        % ('cross-correlation', quality))
		print('  %20s : %s s'                                            % ('run time', sec2hms(time.time()-now, digits=1)))
		print('  %20s : %-d'                                             % ('length data', len(data)))


	### Plot compiarson of data into one axis, if wished
	if plot_compare:
		if title:
			outfile = title[:-6] + '1.png' #os.path.join('/home/scholz/Desktop/glitch_test/plots', title+'_%s%s_s%s_m%s_%0.f%%.png' % (set_to_zero.upper(), threshold, sin_width, malm_width, quality))
		else:
			outfile = None

		quick_plot(data, data_recon, data_labels=('original', 'reconstructed'), title='%s (%s, %s)' % (title, set_to_zero, threshold), lw=1, show=True, outfile=outfile, verts=[(len(data)+2*pad_width)*0.681-pad_width])


	### Plot decomposition coefficients, if wished
	if plot_Cs:
		#quick_plot( *plot_param[0][:-1], title='%s-approximation coefficients after zeroing' % method.upper(), lw=1, show=False)
		quick_plot( *plot_param[1], title='%s-detailed coefficients after zeroing' % method.upper(), lw=1, show=True)


	return data_recon
def match_lines_2D(angles1, angles2, amplitudes1, amplitudes2, weight1=0.5, weight2=0.5, periodicity=90):

	"""
	Calculate a quality number (normalised to 0..1) that expresses
	the match of lines on a 2-D sphere. It is accounted for
	both the angles and the and amplitude of the lines.
	Usuful for example if you want to compare SKS splitting measurements
	with their respective SKS preditions.

	`weight1` is the weight of the angles (or 1st variable), meaning, how 
	much more important you think they are compared to the line amplitudes.

	Note that `angles1`, `angles2`, `amplitudes1`, `amplitudes2` are arrays (or tuples)
	and must have all the same array length.

	`periodicity` means that angles `periodicity` degree apart 
	are considered a perfect angle match)
	"""

	angles1     = np.array( angles1 ).astype('float')*M.pi/180
	angles2     = np.array( angles2 ).astype('float')*M.pi/180
	amplitudes1 = np.array( amplitudes1 ).astype('float')
	amplitudes2 = np.array( amplitudes2 ).astype('float')
	
	angles_diff = np.abs(angles1-angles2)%(periodicity/180*M.pi)
	summands    = weight1 * np.abs( np.cos( 180/periodicity*np.abs( angles_diff ) )) + weight2 * np.minimum(amplitudes1, amplitudes2)/np.maximum(amplitudes1, amplitudes2)

	return np.sum(summands)/len(summands)
def xyz_smoothing(xyz_input, mode='cubic', length=1, weights_length=[1], xyz_output='/Users/john/Dropbox/JOHNS_WORK/projects/rhum_rum/maps_images/surface_wave_tomo/LABeach5_1100/LABeach5_1100_smoothedLAB_cubic_l1_w1.txt'):

	"""
	IMPORTANT:
	----------

	This script only works if the input file is regular in x, and y, and has clear cutedges.
	For example:

		o  o  o  o  o
		o  o  o  o  o
		o  o  o  o  o

	works correctly, whilst

		o  o  o  o  o                  o  o  o  o  o 
		o  o  o  o  o  o      or       o  o  -  o  o 
		o  o  o  o  o                  o  o  o  o  o 

	would give incorrect results. This is because the input file is first sorted after x,y to then
	calculate the line indices of lines that will be used for smoothing (weighted means).
	For correct index calculations, the input file must hence match the stated requirements. 
	The idea behind is that index calculation is much faster than searching for values using for 
	example np.where(), which can take tremendously long for large files.

	---

	Take a file with three or more columns and iterate over x, and y and calculate the
	weighted mean of all other columns accordant to 'mode', 'length', and 'weights_length'.

	Example:
		
		- mode='cubic' & length=1:

				z 	z 	z
				z 	Z 	z
				z 	z 	z

		- mode='cubic' & length=2:

			z   z   z   z   z
			z   z   z   z   z
			z   z   Z   z   z
			z   z   z   z   z
			z   z   z   z   z

		- mode='plus' & length=1
	
				- 	z 	-
				z 	Z 	z
				- 	z 	-

		- mode='plus' & length=2:

			- 	- 	z 	- 	-
			-	- 	z 	- 	- 
			z	z 	Z 	z 	z
			-	- 	z 	- 	-
			- 	- 	z 	- 	-

		- ...


	Intuitively, length=0 will reproduce the input file for both modes.

	'weights_length' is a list of weights to be used. If it has only one value, this weight is
	assigned to all considered z_values, with the exception of the centre value, which always
	has weight=1! If more than one value is given, it must have as many entries as the value of
	'length'. The position of the weight in the list 'weights_length' controls to which z-values
	the weight is assgined. For example, weights_length=[1, 0.5] means that the first ring of 
	values (mode=cubic) or the first nearest neighboors (mode=plus), respectively, around the 
	center have weight=1, the second ring/nearest neighboors weight=0.5, and so on...

	Encountered nan and infinity values are not considered and may result in funny weighted means.

	Output is written as: x y weighted_mean_col3 weighted_mean_col4 ... 
	This script can be used for example to average z-values of a file that will then be used for plotting.
	"""


	now = time.time()


	### PARSE passed variables
	# prepare weights for values 
	if len(weights_length) == 1:
		weights_length = weights_length*length
	elif len(weights_length) != length:
		print("ERROR:   No valid list 'weights_length' given.")
		print("         Either give one entry as weight for all")
		print("         or as many entries as the value of 'length'.")
		sys.exit()
	weights_length.insert(0, 1)				# weight on main (centre) value is always 1!
	#print('weights_length:  %s' % weights_length)
	
	# ensure that code works fine even for negative 'length' (just in case)
	if length==0:						
		index_range = np.arange(-length, length+1, 1)
	else:
		index_range = np.arange(-length, length+1, np.sign(length)*1)		


	### SORT INPUT FILE: first X, then Y, then Z from small to large (important because script demands that)
	try:
		xyz    = np.loadtxt( xyz_input, dtype='float' )
	except ValueError:
		print("ERROR:   Couldn't assign data as floats.")
		print("         Please make sure that only floats are given")
		print("         in each column. 'Nan' and 'inf' values are ok.")		
		sys.exit()
	xyz_sorted = np.array( sorted( xyz, key=lambda x: (float(x[0]),float(x[1]),float(x[2])) ))


	### ASSIGN variables
	columns           = len(      xyz_sorted[0,:] )
	x_size            = len( set( xyz_sorted[:,0] ))
	y_size            = len( set( xyz_sorted[:,1] ))
	y_min             = np.amin(  xyz_sorted[:,1] )
	y_max             = np.amax(  xyz_sorted[:,1] )


	### loop linewise over sorted x,y. Sort is important, as indices will be calculated
	### to adress the lines that are considered for smoothing (rather than looking for the lines
	### via np.where for example, which can take tremendously long for large files).
	### append smoothed columns to results, this list will be printed to file at the end
	results = []

	for i in range(len( xyz_sorted )):
	
		values_to_smooth = np.array([]).reshape(0,columns-2)
		weights          = []

		try:
			if mode=='cubic':

				for index in index_range:
					for index2 in index_range:
						
						# IndexError if we are at the edge of file (adress conditions for x, index2 conditions for y)
						adress = i+index*y_size+index2			# index of line considered for smoothing
						if adress<0 or adress>x_size*y_size-1 or (index2<0 and xyz_sorted[adress][1]==y_max) or (index2>0 and xyz_sorted[adress][1]==y_min):
							raise IndexError

						# get vals (all columns but the x and y) of line and raise exception (skip smoothing for this point) if nan or inf
						vals = xyz_sorted[adress][2:]
						#if np.isnan(z_val) or np.isinf(z_val):
						#	raise IndexError
						values_to_smooth = np.append( values_to_smooth, [vals], axis=0 )		# matrix with lines for each line of file to smooth, columns are columns>2 of file

						# assign weight to adressed line
						ind    = max(abs(index), abs(index2))
						weight = weights_length[ind]
						weights.append( weight )


			if mode=='plus':
	
				for index in index_range:

					# IndexError if we are at the edge of file (adress conditions for x)
					adress = i+index*y_size
					if adress<0 or adress>x_size*y_size-1:
						raise IndexError
					
					# get vals (all columns but the x and y) of line  and raise exception (skip smoothing for this point) if nan or inf
					vals = xyz_sorted[adress][2:]
					#if np.isnan(z_val) or np.isinf(z_val):
					#	raise IndexError
					values_to_smooth = np.append( values_to_smooth, [vals], axis=0 )
	
					# assign weight to adressed line
					weight = weights_length[abs(index)]
					weights.append( weight )
	

				for index in index_range:
					if index!=0:

						# IndexError if we are at the edge of file (index2 conditions for y)
						adress = i+index
						if (index<0 and xyz_sorted[adress][1]==y_max) or (index>0 and xyz_sorted[adress][1]==y_min):
							raise IndexError

						# get vals (all columns but the x and y) of line and raise exception (skip smoothing for this point) if nan or inf
						vals = xyz_sorted[adress][2:]
						#if np.isnan(z_val) or np.isinf(z_val):
						#	raise IndexError
						values_to_smooth = np.append( values_to_smooth, [vals], axis=0 )
						
						# assign weight to adressed line
						weight = weights_length[abs(index)]
						weights.append( weight )


			### RESULTS						
			x_val          = xyz_sorted[i][0]
			y_val          = xyz_sorted[i][1]
			values_meaned  = np.average( values_to_smooth, weights=weights, axis=0 )
			
			results_line_i = np.insert( values_meaned, 0,  y_val )
			results_line_i = np.insert( results_line_i, 0, x_val )

			results.append( results_line_i )
			print( '\t'.join( [str(e) for e in results_line_i] ))


		# if IndexError, accordant line for smoothing is not there (e.g., at the edge of the file).
		# This point is skipped and hence not smoothed, neither written to file.
		except IndexError:
			#print( xyz_sorted[i][0], xyz_sorted[i][1], xyz_sorted[i][2], xyz_sorted[i][3] )
			continue


	### OUTPUT
	# to file
	header = 'This file was produced using python script "xyz_smoothing(%s, mode=%s, length=%s, weights_length=%s, xyz_output=%s)"\n' % (xyz_input, mode, length, weights_length, xyz_output)
	np.savetxt( xyz_output, results, fmt='%10.8s', header=header )

	# to shell print statement of needed time 
	print( u'Done in: %s (h:m:s). Check file & header:' % sec2hms( float( time.time()-now )) )
	print( xyz_output )
def xywz_smoothing(xywz_input, mode='cubic', length=1, weights_length=[1], xywz_output='/Users/john/Dropbox/JOHNS_WORK/projects/rhum_rum/maps_images/surface_wave_tomo/vs_model/vs_smoothed_plus_l1_w1_.xyz'):

	"""
	IMPORTANT:
	----------

	This script only works if the input file is regular in x, y, and w, and has clear cut edges.
	For example:

		o  o  o  o  o
		o  o  o  o  o
		o  o  o  o  o

	works correctly, whilst

		o  o  o  o  o                  o  o  o  o  o 
		o  o  o  o  o  o      or       o  o  -  o  o 
		o  o  o  o  o                  o  o  o  o  o 

	would give incorrect results. This is because the input file is first sorted after x,y,w to then
	calculate the line indices of lines that will be used for smoothing (weighted means).
	For correct index calculations, the input file must hence match the stated requirements. 
	The idea behind is that index calculation is much faster than searching for values using for 
	example np.where(), which can take tremendously long for large files.

	---

	Take a file with four or more columns and iterate over x, y and w, and calculate the
	weighted mean of all other columns accordant to 'mode', 'length', and 'weights_length'.

	Example:
		
		- mode='cubic' & length=1:

				z 	z 	z
				z 	Z 	z
				z 	z 	z

		- mode='cubic' & length=2:

			z   z   z   z   z
			z   z   z   z   z
			z   z   Z   z   z
			z   z   z   z   z
			z   z   z   z   z

		- mode='plus' & length=1
	
				- 	z 	-
				z 	Z 	z
				- 	z 	-

		- mode='plus' & length=2:

			- 	- 	z 	- 	-
			-	- 	z 	- 	- 
			z	z 	Z 	z 	z
			-	- 	z 	- 	-
			- 	- 	z 	- 	-

		- ...


	Intuitively, length=0 will reproduce the input file for both modes.

	'weights_length' is a list of weights to be used. If it has only one value, this weight is
	assigned to all considered z_values, with the exception of the centre value, which always
	has weight=1! If more than one value is given, it must have as many entries as the value of
	'length'. The position of the weight in the list 'weights_length' controls to which z-values
	the weight is assgined. For example, weights_length=[1, 0.5] means that the first ring of 
	values (mode=cubic) or the first nearest neighboors (mode=plus), respectively, around the 
	center have weight=1, the second ring/nearest neighboors weight=0.5, and so on...

	Encountered nan and infinity values are not considered and may result in funny weighted means.

	Output is written as: x y w weighted_mean_col4 weighted_mean_col5 ... 
	This script can be used for example to average z-values of a file that will then be used for plotting.
	"""


	now = time.time()


	### PARSE passed variables
	# prepare weights for values 
	if len(weights_length) == 1:
		weights_length = weights_length*length
	elif len(weights_length) != length:
		print("ERROR:   No valid list 'weights_length' given.")
		print("         Either give one entry as weight for all")
		print("         or as many entries as the value of 'length'.")
		sys.exit()
	weights_length.insert(0, 1)				# weight on main (centre) value is always 1!
	#print('weights_length:  %s' % weights_length)
	
	# ensure that code works fine even for negative 'length' (just in case)
	if length==0:						
		index_range = np.arange(-length, length+1, 1)
	else:
		index_range = np.arange(-length, length+1, np.sign(length)*1)		


	### SORT INPUT FILE: first X, then Y, then Z from small to large (important because script demands that)
	try:
		xywz        = np.loadtxt( xywz_input, dtype='float' )
	except ValueError:
		print("ERROR:   Couldn't assign data as floats.")
		print("         Please make sure that only floats are given")
		print("         in each column. 'Nan' and 'inf' values are ok.")		
		sys.exit()
	xywz_sorted = np.array( sorted( xywz, key=lambda x: (float(x[0]),float(x[1]),float(x[2])) ))


	### ASSIGN variables
	columns = len(      xywz_sorted[0,:] )
	x_size  = len( set( xywz_sorted[:,0] ))
	y_size  = len( set( xywz_sorted[:,1] ))
	w_size  = len( set( xywz_sorted[:,2] ))
	y_min   = np.amin(  xywz_sorted[:,1] )
	y_max   = np.amax(  xywz_sorted[:,1] )


	### loop linewise over sorted x,y,w. Sort is important, as indices will be calculated
	### to adress the lines that are considered for smoothing (rather than looking for the lines
	### via np.where for example, which can take tremendously long for large files).
	### append smoothed columns to results, this list will be printed to file at the end
	results = []

	for i in range(len( xywz_sorted )):
	
		values_to_smooth = np.array([]).reshape(0,columns-3)
		weights          = []

		try:
			if mode=='cubic':

				for index in index_range:
					for index2 in index_range:
						
						# IndexError if we are at the edge of file (adress conditions for x, index2 conditions for y)
						adress = i+index*y_size*w_size+index2*w_size			# index of line considered for smoothing
						if adress<0 or adress>x_size*y_size*w_size-1 or (index2<0 and xywz_sorted[adress][1]==y_max) or (index2>0 and xywz_sorted[adress][1]==y_min):
							raise IndexError

						# get vals (all columns but the x and y) of line and raise exception (skip smoothing for this point) if nan or inf
						vals = xywz_sorted[adress][3:]
						#if np.isnan(z_val) or np.isinf(z_val):
						#	raise IndexError
						values_to_smooth = np.append( values_to_smooth, [vals], axis=0 )

						# assign weight to adressed line
						ind    = max(abs(index), abs(index2))
						weight = weights_length[ind]
						weights.append( weight )


			if mode=='plus':
	
				for index in index_range:

					# IndexError if we are at the edge of file (adress conditions for x)
					adress = i+index*y_size*w_size
					if adress<0 or adress>x_size*y_size*w_size-1:
						raise IndexError
					
					# get vals (all columns but the x and y) of line and raise exception (skip smoothing for this point) if nan or inf
					vals = xywz_sorted[adress][3:]
					#if np.isnan(z_val) or np.isinf(z_val):
					#	raise IndexError
					values_to_smooth = np.append( values_to_smooth, [vals], axis=0 )
	
					# assign weight to adressed line
					weight = weights_length[abs(index)]
					weights.append( weight )
	

				for index in index_range:
					if index!=0:

						# IndexError if we are at the edge of file (index2 conditions for y)
						adress = i+index*w_size
						if (index<0 and xywz_sorted[adress][1]==y_max) or (index>0 and xywz_sorted[adress][1]==y_min):
							raise IndexError

						# get vals (all columns but the x and y) of line and raise exception (skip smoothing for this point) if nan or inf
						vals = xywz_sorted[adress][3:]
						#if np.isnan(z_val) or np.isinf(z_val):
						#	raise IndexError
						values_to_smooth = np.append( values_to_smooth, [vals], axis=0 )
						
						# assign weight to adressed line
						weight = weights_length[abs(index)]
						weights.append( weight )


			### RESULTS						
			x_val          = xywz_sorted[i][0]
			y_val          = xywz_sorted[i][1]
			w_val          = xywz_sorted[i][2]
			values_meaned  = np.average( values_to_smooth, weights=weights, axis=0 )

			results_line_i = np.insert( values_meaned, 0,  w_val )
			results_line_i = np.insert( results_line_i, 0, y_val )
			results_line_i = np.insert( results_line_i, 0, x_val )

			results.append( results_line_i )
			print( '\t'.join( [str(e) for e in results_line_i] ))


		# if IndexError, accordant line for smoothing is not there (e.g., at the edge of the file).
		# This point is skipped and hence not smoothed, neither written to file.
		except IndexError:
			continue


	### OUTPUT
	# to file
	header = 'This file was produced using python script "xywz_smoothing(%s, mode=%s, length=%s, weights_length=%s, xywz_output=%s)"\n' % (xywz_input, mode, length, weights_length, xywz_output)
	np.savetxt( xywz_output, np.array( results ), fmt='%10.8s', header=header )

	# to shell print statement of needed time 
	print( u'Done in: %s (h:m:s). Check file & header:' % sec2hms( float( time.time()-now )) )
	print( xywz_output )
def smooth(x, window_len=11, window='blackman'):

    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError( "smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError( "Input vector needs to be bigger than window size." )


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError( "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'" )


    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='same')
    return y
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
	0------> comp_2


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
		angle = 360-angle				# convert angle to as it would be clockwise rotation


	comp_1_new =  comp_1*M.cos(angle*M.pi/180) + comp_2*M.sin(angle*M.pi/180)
	comp_2_new = -comp_1*M.sin(angle*M.pi/180) + comp_2*M.cos(angle*M.pi/180)

	return comp_1_new, comp_2_new
def normalise(data, scale_to_between=[]):

	"""
	Normalise passed data (array-like).
	`scale_to_between` is list with 2 elements.
	"""

	data = np.array(data)

	
	if len(data) == 0:
		return data

	if scale_to_between:
		scale_to_between.sort()
		if len(data) == 1:
			data /= data * scale_to_between[1]
			return data
		scale  = abs(scale_to_between[-1]-scale_to_between[0]) / 2.
		drange = max(data)-min(data)
		data   = data * 2 / drange
		data   = data - max(data)+1				# data have y values between [-1,1]
		data  *= scale							# data have y values between [-1/scale,1/scale] 
		data  += scale_to_between[-1]-scale 	# eacht trace has y values filling range `scale_to_between`
	
	else:
		if len(data) == 1:
			data /= data
			return data
		data /= max(abs(data))

	return data
def coherence(data1, data2, fs=1.0, window='hann', freq_bounds=[], plot=False):

	"""
	La cohérence avec la pression avec la composante Z (plutôt la journée) entre 0.05Hz et 0.2 Hz.
	Et les horizontaux à plus basses fréquences (0.02-0.1).


	Check:
	  https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.welch.html#r34b375daf612-1
	  https://docs.scipy.org/doc/scipy-0.18.1/reference/generated/scipy.signal.csd.html


	Return: - one-sided frequencies (positive), 
			- coherence array, 
			- coherence value as average of coherence arraybetween `freq_bounds`
	"""

	data1     = np.array( data1 )
	data2     = np.array( data2 )

	f11, P11  = scipy.signal.welch(data1,      fs=fs, window=window, nperseg=None, noverlap=None, nfft=None, detrend='constant', return_onesided=True, scaling='density', axis=-1, average='mean')
	f22, P22  = scipy.signal.welch(data2,      fs=fs, window=window, nperseg=None, noverlap=None, nfft=None, detrend='constant', return_onesided=True, scaling='density', axis=-1, average='mean')
	f12, P12  = scipy.signal.csd(data1, data2, fs=fs, window=window, nperseg=None, noverlap=None, nfft=None, detrend='constant', return_onesided=True, scaling='density', axis=-1) 		#P12 is complex

	coh_array = P12 * np.conj(P12) / (P11 * P22)

	if freq_bounds:
		idx1  = ( np.abs( np.array(f11)-freq_bounds[0]) ).argmin()
		idx2  = ( np.abs( np.array(f11)-freq_bounds[1]) ).argmin()
		verts = (f11[idx1], f11[idx2])
	else:
		idx1  = None
		idx2  = None
		verts = ()

	coh_value = np.real( np.average(coh_array[idx1:idx2]) )

	if plot:
		quick_plot(coh_array, x=f11, data_labels=('Coherence(freq) = P12 * conj(P12) / (P11 * P22)',), xlabel='Frequency (Hz)', ylabel='Coherence', lw=1.5, ylim=[0,1], verts=verts)
	
	return f11, coh_array, coh_value
def euclidean_filter(data, window_length_in_samples=10, overlap_in_samples=1, step_in_samples=1):

	"""
	See: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4667093/#SM1
	"""


	data = np.array( data )
	
	new_data = []
	i = 0
	while True:

		start_seg1 =  i    * step_in_samples
		start_seg2 =  i    * step_in_samples+(window_length_in_samples-overlap_in_samples)
		end_seg1   = (i+1) * step_in_samples+window_length_in_samples
		end_seg2   = (i+1) * step_in_samples+window_length_in_samples+(window_length_in_samples-overlap_in_samples)
		segment1 = data[start_seg1:end_seg1]
		segment2 = data[start_seg2:end_seg2]

		i += 1
		if end_seg2>=len(data):
			break

		value = np.sqrt( np.sum((segment1-segment2)**2) )
		new_data.append(value)

	return np.array(new_data)
def Zscore_filter(data, threshold=3, window_length_in_samples=10, influence=1):

	"""
	Moving windows, step size 1 sample, window of length 
	`window_length_in_samples`, each for which Z-score is calulated. 

	Note that this Z-score of each new data point is calculated
	against a window that ends on the previous sample, allowing
	to introduce the `influence` parameter.

	`influence`:
	If put at 0, detections have no influence on the threshold, such that future detections 
	are detected based on a threshold that is calculated with a mean and standard 
	deviation that is not influenced by past detections. Another way to think about 
	this is that if you put the influence at 0, you implicitly assume stationarity 
	(i.e. no matter how many signals there are, the time series always returns to 
	the same average over the long term). If this is not the case, you should put 
	the influence parameter somewhere between 0 and 1, depending on the extent to 
	which signals can systematically influence the time-varying trend of the data. 
	E.g., if signals lead to a structural break of the long-term average of the time 
	series, the influence parameter should be put high (close to 1) so the threshold 
	can adjust to these changes quickly.


	Return:
	  Array of Z-scores. 


	 Note:
	  - influence=1 will return normal Z-scores of a moving window
	    (i.e., `threshold` has no impact at all)


	After:
	  http://stackoverflow.com/a/22640362/6029703
	"""

	data       = np.array(data)
	data_alter = np.array(data)
	z_scores   = np.zeros( len(data) )

	for i in range(window_length_in_samples, len(data)-1):
	
		mean    = np.mean(data_alter[i-window_length_in_samples:i])
		std     = np.std( data_alter[i-window_length_in_samples:i])
		z_score = (data[i] - mean) / std

		if abs(z_score) >= threshold:
			data_alter[i] = influence * data[i] + (1-influence) * data_alter[i-1]
		else:
			data_alter[i] = data[i]

		z_scores[i] = z_score	

	return np.asarray(z_scores)

# seismological
def x2ded(position, digits=6):

	"""
	Take any kind of gps data ('xx.', 'xx xx', 'xx xx xx', '-xx xx.xx', 'xx,xxxx')
	as string and return DEcimal Degree converted position as string.

	Keep in mind, cardinal directions such as 'W', 'E', 'S' and 'N' are accepted,
	but the actual sign must still be given, e.g. - for S/W. 
	
	:param position, contains GPS position as string in arbitrary format
	:return: final_position (decimal deg) as string. If == '-12345.0', input not valid.
	"""


	position = str(position)
	point, k = r"[.,]{1}", []
	regex = r"\A\s*[-+]?\s*(\d+)\s*[.,]?\s*(°|d|deg)?\s*[NS]?\s*\Z|" \
			r"\A\s*[-+]?\s*(\d+)\s*[.,]{1}\s*(\d+)\s*(°|d|deg)?\s*[NSEW]?\s*\Z|" \
			r"\A\s*[-+]?\s*(\d+)\s*(°|d|deg)?\s*(\d+)\s*('|`|m|min)?\s*[NSEW]?\s*\Z|" \
			r"\A\s*[-+]?\s*(\d+)\s*(°|d|deg)?\s*(\d+)\s*[.,]{1}\s*(\d*)\s*('|`|m|min)?\s*[NSEW]?\s*\Z|" \
			r"\A\s*[-+]?\s*(\d+)\s*(°|d|deg)?\s*(\d+)\s*('|`|m|min)?\s*(\d+)\s*(''|``|\"|s|sec)?\s*[NSEW]?\s*\Z"
	# decimal without comma digits
	# decimal with comma digits
	# decimal and minutes without seconds
	# decimal and decimal minutes
	# decimal and minute and seconds
	
	
	try:
	
		### decode the 'raw_input' string to 'UTF-8' to ensure correct coding
		position = position.decode("utf-8")

		### check passed string and append not none blocks to list 'k'
		# is 'regex' pattern int 'position' ?
		prog = re.compile(regex, flags=re.UNICODE)
		for i in [e for e in prog.search(position).groups() if e != None and re.search(r'\A\d*\Z', e)]:
			if i == '' : 
				i = '0'
			k.append(i)

		### actual conversions
		if len(k) == 1:											# decimal without comma
			final_position = '%s.0' % k[0]
		elif len(k) == 2 and re.search(point, position):		# decimal with comma
			final_position = '%s.%s' % ( k[0], k[1] )
		elif len(k) == 2:										# minute without seconds
			final_position = str( float( k[0] ) + float( k[1] )/60. )
		elif len(k) == 3 and re.search(point, position):		# decimal minute 
			minute = float( '%s.%s' % ( k[1], k[2] )) / 60
			final_position = str( float( k[0] ) + minute )
		elif len(k) == 3:										# minute with seconds
			second = float( k[2] )/60
			decimal = (float( k[1] )+second)/60
			final_position = str( float( k[0] ) + decimal )
										
		### consider sign / NS
		if (re.search(r"\A\s*-", position) or re.search(r"\s*S\s*\Z", position)) and not re.search(r"\s*N\s*\Z", position): 
			final_position = '-' + final_position

	except Exception:										# in case of invalid gps coordinates
		final_position = '-12345.0'

	return ("%%.%df" % digits) % float(final_position)
def x2dem(position, digits=6):

	"""
	Take any kind of gps data ('xx.', 'xx xx', 'xx xx xx', '-xx xx.xx', 'xx,xxxx')
	as string and return DEcimal Minutes converted position as string.
	
	Keep in mind, cardinal directions such as 'W', 'E', 'S' and 'N' are accepted,
	but the actual sign must still be given, e.g. - for S/W. 

	:param position, contains GPS position as string in arbitrary format
	:return: final_position (decimal deg) as string. If == '-12345.0', input not valid.
	"""


	position = str(position)
	point, k = r"[.,]{1}", []
	regex = r"\A\s*[-+]?\s*(\d+)\s*[.,]?\s*(°|d|deg)?\s*[NS]?\s*\Z|" \
			r"\A\s*[-+]?\s*(\d+)\s*[.,]{1}\s*(\d+)\s*(°|d|deg)?\s*[NSEW]?\s*\Z|" \
			r"\A\s*[-+]?\s*(\d+)\s*(°|d|deg)?\s*(\d+)\s*('|`|m|min)?\s*[NSEW]?\s*\Z|" \
			r"\A\s*[-+]?\s*(\d+)\s*(°|d|deg)?\s*(\d+)\s*[.,]{1}\s*(\d*)\s*('|`|m|min)?\s*[NSEW]?\s*\Z|" \
			r"\A\s*[-+]?\s*(\d+)\s*(°|d|deg)?\s*(\d+)\s*('|`|m|min)?\s*(\d+)\s*(''|``|\"|s|sec)?\s*[NSEW]?\s*\Z"
	# decimal without comma digits
	# decimal with comma digits
	# decimal and minutes without seconds
	# decimal and decimal minutes
	# decimal and minute and seconds
	
	
	try:
		### decode the 'raw_input' string to 'UTF-8' to ensure correct coding
		position = position.decode("utf-8")

		### check passed string and append not none blocks to list 'k'
		# is 'regex' pattern in 'position' ?
		prog = re.compile(regex, flags=re.UNICODE)
		for i in [e for e in prog.search(position).groups() if e != None and re.search(r'\A\d*\Z', e)]:
			if i == '' : 
				i = '0'
			k.append(i)

		### actual conversion
		if len(k) == 1:												# decimal without comma
			final_position = '%s 00.00' % k[0].zfill(2)
		elif len(k) == 2 and re.search(point, position):			# decimal with comma
			minute = str( float('.' + k[1])*60 ).zfill(2)
			if len( minute.split('.')[0] ) < 2:
				minute = minute.zfill( len(minute)+1 )
			if len( minute.split('.')[1] ) < 2:
				minute = minute + '0'
			final_position = k[0].zfill(2) + ' ' + minute.zfill(2)
		elif len(k) == 2:											# minute without seconds
			minute = str(float('.' + k[1])*60).zfill(2)
			final_position = '%s %s.00' % (k[0].zfill(2), k[1].zfill(2))
		elif len(k) == 3 and re.search(point, position):			# decimal minute
			degree = str( float(k[0]) + float(k[1] + '.' + k[2] )/60)
			minute = str( float( '.' + degree.split('.')[1] )*60).zfill(2)
			if len( minute.split('.')[0] ) < 2:
				minute = minute.zfill( len(minute)+1 )
			final_position = '%s %s' % (degree.split('.')[0].zfill(2), minute)
		elif len(k) == 3:											# minute with seconds
			second = float(k[2])/60
			decimal = (float(k[1]) + second)/60
			degree = str( float( k[0] ) + decimal )
			minute = str( float( '.' + degree.split('.')[1] )*60)
			if len( minute.split('.')[0] ) < 2:
				minute = minute.zfill( len(minute)+1 )
			final_position = '%s %s' % (degree.split('.')[0].zfill(2), minute.zfill(2))

		### consider sign / NS
		if (re.search(r"\A\s*-", position) or re.search(r"\s*S\s*\Z", position)) and not re.search(r"\s*N\s*\Z", position): 
			final_position = '-' + final_position

	except Exception:												# in case of invalid gps coordinates
		final_position = '-12 34.50'

	final_position = final_position.split()[0] + ' ' + ("%%.%df" % digits) % float(final_position.split()[1])
	return(final_position)
def x2ams(position, digits=6):

	"""
	Take any kind of gps data ('xx.', 'xx xx', 'xx xx xx', '-xx xx.xx', 'xx,xxxx')
	as string and return Arc Min Sec converted position as string.
	
	Keep in mind, cardinal directions such as 'W', 'E', 'S' and 'N' are accepted,
	but the actual sign must still be given, e.g. - for S/W. 

	:param: position, contains GPS position as string in arbitrary format
	:return: final_position (decimal deg) as string. If == '-12345.0', input not valid.
	"""


	position = str(position)
	point, k = r"[.,]{1}", []
	regex = r"\A\s*[-+]?\s*(\d+)\s*[.,]?\s*(°|d|deg)?\s*[NS]?\s*\Z|" \
			r"\A\s*[-+]?\s*(\d+)\s*[.,]{1}\s*(\d+)\s*(°|d|deg)?\s*[NSEW]?\s*\Z|" \
			r"\A\s*[-+]?\s*(\d+)\s*(°|d|deg)?\s*(\d+)\s*('|`|m|min)?\s*[NSEW]?\s*\Z|" \
			r"\A\s*[-+]?\s*(\d+)\s*(°|d|deg)?\s*(\d+)\s*[.,]{1}\s*(\d*)\s*('|`|m|min)?\s*[NSEW]?\s*\Z|" \
			r"\A\s*[-+]?\s*(\d+)\s*(°|d|deg)?\s*(\d+)\s*('|`|m|min)?\s*(\d+)\s*(''|``|\"|s|sec)?\s*[NSEW]?\s*\Z"
	#decimal without comma digits
	#decimal with comma digits
	#decimal and minutes without seconds
	#decimal and decimal minutes
	#decimal and minute and seconds
	

	try:
	
		### decode the 'raw_input' string to 'UTF-8' to ensure correct coding
		position = position.decode("utf-8")

		### check passed string and append not none blocks to list 'k'
		# is 'regex' pattern in 'position' ?
		prog = re.compile(regex, flags=re.UNICODE)
		for i in [e for e in prog.search(position).groups() if e != None and re.search(r'\A\d*\Z', e)]:
			if i == '' : 
				i = '0'
			k.append(i)

		### actual conversion
		if len(k) == 1:												# decimal without comma
			final_position = '%s 00 00' % k[0].zfill(2) 
		elif len(k) == 2 and re.search(point, position):			# decimal with comma
			minute = str( float('.' + k[1])*60 ).split('.')
			second = str( float('.' + minute[1])*60 )
			final_position = k[0].zfill(2) + ' ' + minute[0].zfill(2) + ' ' + second.zfill(2)
		elif len(k) == 2:											# minute without seconds
			d, m = divmod( int( k[1] ), 60)
			final_position = '%s %s 00' % ( int( k[0] ) + d, m ) 
		elif len(k) == 3 and re.search(point, position):			# decimal minute
			decimal = float(k[1] + '.' + k[2])/60
			degree = str(float(k[0]) + decimal)
			minute = str(float('.' + degree.split('.')[1])*60. ).split('.')
			second = str( float('.' + minute[1])*60. )
			final_position = degree.split('.')[0].zfill(2) + ' ' + minute[0].zfill(2) + ' ' + second.zfill(2)
		elif len(k) == 3:											# minute with seconds
			degree = str((float(k[0])+(float(k[1])+float(k[2])/60)/60)).split('.')
			minute = str(float('.' + degree[1])*60).split('.')
			second = str( float('.' + minute[1])*60. )
			final_position = degree[0] + ' ' + minute[0].zfill(2) + ' ' + second.zfill(2)

		### consider sign / NS
		if (re.search(r"\A\s*-", position) or re.search(r"\s*S\s*\Z", position)) and not re.search(r"\s*N\s*\Z", position): 
			final_position = '-' + final_position

	except Exception:												# in case of invalid gps coordinates
		final_position = '-12345.0'

	final_position = final_position.split()[0] + ' ' + final_position.split()[1] + ' ' + ("%%02.%df" % digits) % float(final_position.split()[2])
	return(final_position)
def read(file=None):
	# wrapper to make to return Stream2 objects instead of (ObsPy's) Stream object.
	st        = obspy.read(file)
	st        = Stream2(st, file=file)
	return st
class Stream2(Stream):

	"""	Extended class for ObsPy's stream object. """

	def __init__(self, traces=None, file=None):
		super().__init__(traces=traces)
		self.origin       = str(file)
		self.original     = None
		self.removed      = None
		self.gain_removed = False
		self.unit         = 'RAW'
		self.times        = self._get_times()
		self.inventory    = None
		self.filters      = {'0' : {'freqmin':None,  'freqmax':None, 'corners':3, 'zerophase':True, 'type_taper':'hann', 'max_percentage_taper':0.03, 'num':0, 'string':'0: None'},
		                     '1' : {'freqmin':1.0,   'freqmax':None, 'corners':3, 'zerophase':True, 'type_taper':'hann', 'max_percentage_taper':0.03, 'num':1, 'string':'1: >1 Hz'},
					         #'2' : {'freqmin':0.1,   'freqmax':None, 'corners':3, 'zerophase':True, 'type_taper':'hann', 'max_percentage_taper':0.03, 'num':2, 'string':'2: >0.1 Hz'},
		                     #'3' : {'freqmin':None,  'freqmax':None, 'corners':3, 'zerophase':True, 'type_taper':'hann', 'max_percentage_taper':0.03, 'num':3, 'string':'3: >0.01 Hz'},
		                     '2' : {'freqmin':None,  'freqmax':1,    'corners':3, 'zerophase':False,'type_taper':'hann', 'max_percentage_taper':0.03, 'num':2, 'string':'2: <1 Hz (0-phase=False)'},
		                     '3' : {'freqmin':None,  'freqmax':1,    'corners':3, 'zerophase':True, 'type_taper':'hann', 'max_percentage_taper':0.03, 'num':3, 'string':'3: <1 Hz (0-phase=True)'},
		                     '4' : {'freqmin':0.001, 'freqmax':None, 'corners':3, 'zerophase':True, 'type_taper':'hann', 'max_percentage_taper':0.03, 'num':4, 'string':'4: >0.001 Hz'},
		                     #'5' : {'freqmin':1./100,   'freqmax':None,  'corners':3, 'zerophase':True, 'type_taper':'hann', 'max_percentage_taper':0.03, 'num':5, 'string': '5: Lowpass 100s'},
		                     #'6' : {'freqmin':1./200,   'freqmax':None,  'corners':3, 'zerophase':True, 'type_taper':'hann', 'max_percentage_taper':0.03, 'num':5, 'string': '6: Lowpass 200s'},
		                     #'7' : {'freqmin':1./300,   'freqmax':None,  'corners':3, 'zerophase':True, 'type_taper':'hann', 'max_percentage_taper':0.03, 'num':5, 'string':'7: Lowpass 300s'},
		                     #'8' : {'freqmin':1./400,   'freqmax':None,  'corners':3, 'zerophase':True, 'type_taper':'hann', 'max_percentage_taper':0.03, 'num':5, 'string':'8: Lowpass 400s'},
		                     #'9' : {'freqmin':1./500,   'freqmax':None,  'corners':3, 'zerophase':True, 'type_taper':'hann', 'max_percentage_taper':0.03, 'num':5, 'string':'9: Lowpass 500s'},
		                     '5' : {'freqmin':0.1,   'freqmax':2.0,  'corners':3, 'zerophase':True, 'type_taper':'hann', 'max_percentage_taper':0.03, 'num':5, 'string':'5: RFs (0.1-2 Hz)'},
		                     '6' : {'freqmin':1/9.0, 'freqmax':1.0,  'corners':3, 'zerophase':True, 'type_taper':'hann', 'max_percentage_taper':0.03, 'num':6, 'string':'6: BP 1-9 s'},
		                     '7' : {'freqmin':1/6.0, 'freqmax':1.0,  'corners':3, 'zerophase':True, 'type_taper':'hann', 'max_percentage_taper':0.03, 'num':7, 'string':'7: BP 1-6 s'},
		                     '8' : {'freqmin':2.0,   'freqmax':3.0,  'corners':3, 'zerophase':True, 'type_taper':'hann', 'max_percentage_taper':0.03, 'num':8, 'string':'8: 2-3 Hz'},
		                     '9' : {'freqmin':2.0,   'freqmax':10,   'corners':3, 'zerophase':True, 'type_taper':'hann', 'max_percentage_taper':0.03, 'num':9, 'string':'9: BP 2-10 Hz'},
		                     }
	def is_unique(self):

		"""
		Takes a ObsPy stream object and returns its station code of the format: network.station.location.channel.
		If this is non-unique, i.e., either of them is non-unique, i.e., the stream object contains more than one
		specific component, this script returns `False`. If it is unique, returns `True`.

		This script can be used to make sure there is only one component within the passed stream.
		This is useful when checking for gaps/overlaps or for any other application
		where it's mandatory for correct processing to have only one component contained.
		"""

		networks  = []
		stations  = []
		locations = []
		channels  = []

		for tr in self:
			networks.append(tr.stats.network)
			stations.append(tr.stats.station)
			locations.append(tr.stats.location)
			channels.append(tr.stats.channel)
	
		if len(set(networks))==1 and len(set(stations))==1 and len(set(locations))==1 and len(set(channels))==1:
			return True
		else:
			return False
	def snr(self, axis=0, ddof=1):
		""" 
		Signal-to-noise ratio (SNR), as by Scipy.
		Retuns a dictionary with keys=station ID,
		values=SNR. 
		"""

		SNRs = {}
		for trace in self:
			data           = np.array(trace.data)
			SNRs[trace.id] = snr(data, axis=axis, ddof=ddof)

		for traceID in SNRs.keys():
			print('ID: %s  SNR: %8.2f' % (traceID, SNRs[traceID]))

		return SNRs
	def trim(self, *args, **kwargs):

		"""
		Use Obspy's trim fucntion and improve sliglhty.
		Improvement is that you don't have to type
		`UTCDateTime` every time you specify a time.
		"""

		# Try to convert positional arguments passed to `UTCDateTime` 
		args = list(args)
		for l, arg in enumerate(args):
			try:
				args[l] = UTCDateTime(args[l])
			except Exception as err:
				#print(err)
				pass

		# Try to convert keyword arguments passed to `UTCDateTime` 
		try:
			kwargs['starttime'] = UTCDateTime(kwargs['starttime'])
		except KeyError:
			pass
		try:
			kwargs['endtime'] = UTCDateTime(kwargs['endtime'])
		except KeyError:
			pass

		return super().trim(*args, **kwargs)
	def trim_common(self):

		"""
		Trim traces in stream to common start- and endtime, i.e.,
		maximum starttime and minimum endtime
		"""
	
		max_starttime = max([tr.stats.starttime for tr in self ])
		min_endtime   = min([tr.stats.endtime   for tr in self ])				
		self.trim(starttime=max_starttime, endtime=min_endtime)	
		self.times = self._get_times()
	def truncate(self, start_per=0, end_per=1):
	
		"""
		Truncate stream bei percentage.
		Default values mean nothing is trucated.
		"""
	
		start_per = float(start_per)
		end_per   = float(end_per)
	
		if start_per<0: start_per = 0
		if end_per>1: end_per = 1
		if start_per==0 and end_per==1:
			return
	
		mint       = min([tr.stats.starttime for tr in self])
		maxt       = max([tr.stats.endtime for tr in self])
		starttime  = mint + (maxt-mint)*start_per
		endtime    = maxt - (maxt-mint)*(1-end_per)
		self.trim(starttime=starttime, endtime=endtime)
		self.times = self._get_times()
	def normalise(self, scale_to_between=[]):

		"""
		Normalise each trace in stream.
		scale_to_between is list with 2 elements.
		"""

		for tr in self:
			tr.data = normalise(tr.data, scale_to_between=scale_to_between)
	def print_stats(self):

		"""
		Print stats for each each within the stream.
		"""

		print('\nStats:')
		for trace in self:
			print('  %s' % trace.id)
			for info in trace.stats.keys():
				if info == trace.stats._format.lower() or info == 'processing':
					continue
				print('%17s : %s' % (info, trace.stats[info]))
			print('%17s : %s' % ('.. and', 'processing_steps & %s-dictionary' % trace.stats._format))
	def print_process(self, truncate=100):
	
		"""
		Print Processing steps for all trace objects
		contained in the passed ObsPy stream object.
		"""
	
		print('\nProcessing steps:')
		for trace in self:
			print('  %s:' % trace.id)
			try:
				processed_steps = trace.stats.processing
			except AttributeError:
				print('     Nothing processed.')
			else:
				for processed_step in processed_steps:
					matches = re.findall(r'(\w*\(.*\))\Z',processed_step)						
					for match in matches:
						for i in range( M.ceil(len(match)/truncate) ):
							if i == 0:
								print('    %s' % match[i*truncate:(i+1)*truncate])
							else:
								print('        %s' % match[i*truncate:(i+1)*truncate])
	def print_filter(self):
		"""
		Print filter of stream, if it was fitlered.
		"""
		print( self.current_filter_str )
	def gain_correction(self): 
		if self.inventory:
			print()
			print(u'GAIN CORRECTION APPLIED')
			try:

				if not self.gain_removed:
					self.gain_removed = True
					for trace in self:
						response   = self.inventory.get_response(trace.id, trace.stats.starttime)
						gain       = response._get_overall_sensitivity_and_gain()[1]				
						trace.data = trace.data / gain
						print(u'  %15s : overall sensitivity and gain (division) %s' % (trace.id, gain))
				else:
					self.gain_removed = False
					for trace in self:
						response   = self.inventory.get_response(trace.id, trace.stats.starttime)
						gain       = response._get_overall_sensitivity_and_gain()[1]				
						trace.data = trace.data * gain
						print(u'  %15s : overall sensitivity and gain (multiplication) %s' % (trace.id, gain))

			except Exception as err:
				print(u'WARNING:  %s' % err)

		else:
			print()
			print(u'No matching response found for gain correction. Nothing done.')
	def filtering(self, filter):
	
		""" 
		Filter ObsPy's stream object has to the given filter specifications. 
		If only freqmax or freqmin is given, highpass or lowpass filter is used, respectively.
		By default the Butterworth filter implementation is used. For a different filter,
		change the code!
	
		See details:
		  https://docs.obspy.org/packages/autogen/obspy.core.stream.Stream.filter.html
		"""
	
		if isinstance(filter, (str, int)):
			filter = self.filters[str(filter)]


		#filter = {**taper, **filter}

		if filter['freqmin'] or filter['freqmax']:
			self.detrend('demean')
			self.detrend('linear')
			self.taper(type=filter['type_taper'], max_percentage=filter['max_percentage_taper'])
			if not filter['freqmin']:
				self.filter('lowpass', freq=filter['freqmax'], corners=filter['corners'], zerophase=filter['zerophase'])
			elif not filter['freqmax']:
				self.filter('highpass', freq=filter['freqmin'], corners=filter['corners'], zerophase=filter['zerophase'])
			else:
				self.filter('bandpass', freqmin=filter['freqmin'], freqmax=filter['freqmax'], corners=filter['corners'], zerophase=filter['zerophase'])

			self.current_filter_num = filter['num']
			self.current_filter_str = filter['string']
	def _get_filter_str(self):
		# return filter of stream object (use only first trace to check)
		try:
			filt = self.current_filter_str

		except AttributeError:
			self.current_filter_num = 0
			filt = 'no filter'
			try:
				for processed_step in self[0].stats.processing:
					if 'filter' in processed_step:
						filt = re.findall(r'filter\(options={(.*)}', processed_step)[0]
						break
			except AttributeError:
				pass

		return filt
	def _get_components(self):
		components = '*'
		try:
			components = sorted( list( set( [tr.stats.channel[-1] for tr in self] )), reverse=True)
			components = ''.join(components)
		except AttributeError:
			pass

		return components
	def _get_channels(self):
		channels = '*'
		try:
			channels = sorted( list( set( [tr.stats.channel for tr in self] )), reverse=False)
		except AttributeError:
			pass

		return channels
	def _get_ids(self):
		ids = '*'
		try:
			ids = sorted( list( set( [tr.id for tr in self] )), reverse=False)
		except AttributeError:
			pass

		return ids
	def _get_times(self):
		try:
			times = min([tr.stats.starttime for tr in self]), max([tr.stats.endtime for tr in self])
		except ValueError:  # empty stream object
			times = None, None

		return times
	def _get_inventory(self, file='/home/scholz/Downloads/dataless*.seed'):

		"""
		Gets latest `file` (with respect to download) and reads it as an Obspy inventory object.
		From this inventory, only the part is extracted that matches the start and end times of the stream `self`,
		as well as matches all networks, stations, locations and channels of the stream.

		The inventory file is then assigned to `self.inventory` and also returned in case
		further processing is needed.
		"""


		inv = Inventory()

		if re.search(r'BH|MH|SH', 'XXX'.join([tr.id for tr in self])):
			inv_file  = max( glob.glob( file ), key=os.path.getctime) 	# latest download
			inventory = read_inventory(inv_file)
			for trace in self:
				network, station, location, channel = trace.id.split('.')
				inv += inventory.select(network   = network, 
										station   = station, 
										location  = location, 
										channel   = channel, 
										starttime = self.times[0], 
										endtime   = self.times[1])

			#print('Read inventory file between %s - %s:' % (self.times[0], self.times[1]))
			#print(inv_file)

		return inv
	def _set_inventory(self):

		inventory      = self._get_inventory()
		self.inventory = inventory
	def rotate_2D(self, angle, components='NE', new_components='12', clockwise=False):

		"""
		"""

		if len(components) != 2:
			print('To rotate data in 2-D, you need to give two')
			print('component letters. You gave %s' % len(components))

		self.trim_common()
		comp_1 = self.select(component=components[0])
		comp_2 = self.select(component=components[1])
		if comp_1 and comp_2:
			comp_1[0].data, comp_2[0].data = rotate_2D(comp_1[0].data, comp_2[0].data, angle, clockwise=clockwise)
			comp_1[0].stats.channel        = comp_1[0].stats.channel[:-1] + new_components[0]
			comp_2[0].stats.channel        = comp_2[0].stats.channel[:-1] + new_components[1]
			print('Rotated `components` (%s) to `new_components` (%s)' % (components, new_components))
			print('`angle` (%s), `clockwise` (%s)' % (angle, clockwise))
		else:
			print('No rotation performed, as `components` (%s)' % components)
			print('are not contained in stream.')
	def rotate2zne(self):

		if not self.inventory:
			self._set_inventory()

		self.rotate('->ZNE', inventory=self.inventory, components=('UVW'))
	def rotation_correlation(self, max_splittime=4):
	
		"""
		"""

		# prep incoming stream
		stream_Q      = self.select(component='Q')
		stream_T      = self.select(component='T')
		if not (stream_Q and stream_T):
				print('')
				print('Return. No idea how to perform splitting on components: %s' % ', '.join([i.stats.channel for i in self]) )
				return
		st_cross_Q    = stream_Q[0]
		st_cross_T    = stream_T[0]		
		window_length = st_cross_Q.stats.endtime - st_cross_Q.stats.starttime
		
		stream2       = self.copy()
		try:
			stream2.trim(stream2[0].stats.starttime+max_splittime, stream2[0].stats.endtime-max_splittime)
		except:
			print('')
			print('Return. Need window length >%ss (>2*max_splittime). Given: %.1f s.' % (2*max_splittime, window_length))
			return
		st_cross_Q2    = stream2.select(component='Q')[0]
		st_cross_T2    = stream2.select(component='T')[0]

		# check entry-condition(s)
		#SNR_Q        = self.snr(st_cross_Q.data)
		#SNR_T        = self.snr(st_cross_T.data)
		#print('SNR_Q    : %s' % SNR_Q)v
		#print('SNR_T    : %s' % SNR_Q)
	
		# determine match wave
		if not self.current_filter:
			freqs = np.linspace(0.01, 0.5*st_cross_Q.stats.sampling_rate, 35)
		else:
			freqs = np.linspace(self.current_filter['freqmin'], self.current_filter['freqmax'], 35)

		coeffs       = []
		for freq in freqs:
			if 1./freq > window_length:
				full_periods = window_length * freq
			else:
				full_periods = 1
			x, sin   = sin_gen(freq, full_periods=full_periods, delta=st_cross_Q.stats.delta)
			sample_shift, cc_factor, cc_factor_abs, corr_array = xcorr(st_cross_Q.data, sin, mode='valid')
			coeffs.append( [sample_shift, cc_factor, cc_factor_abs, corr_array, freq, x, sin] )

		coeffs         = np.array(coeffs)
		arg            = np.argmax( coeffs[:,2] )	#sort after highest absolute cc_factor

		spl_shift_good = coeffs[arg,0]
		cc_good        = coeffs[arg,1]
		cc_abs_good    = coeffs[arg,2]
		corray_good    = coeffs[arg,3]
		freq_good      = coeffs[arg,4]
		x_good         = coeffs[arg,5]
		sin_good       = coeffs[arg,6]
	
		match_part     = st_cross_Q.data[spl_shift_good:spl_shift_good+len(sin_good)]
		#match_part     = sin_good
		match_wave     = np.concatenate(( st_cross_Q.data[:spl_shift_good], match_part, st_cross_Q.data[spl_shift_good+len(match_part):] ))
		#quick_plot(range(len(corray_good)), corray_good)

		###### rotation correlation
		coeffs2        = []
		for rot_angle in np.linspace(0,359,180):
			data_rot_T, data_rot_Q   = rotate_2D(st_cross_T.data, st_cross_Q.data, rot_angle, clockwise=True)
			data_rot_T2, data_rot_Q2 = rotate_2D(st_cross_T2.data, st_cross_Q2.data, rot_angle, clockwise=True)
			sample_shift, cc_factor, cc_factor_abs, corr_array = xcorr(data_rot_T, data_rot_Q2, mode='valid')
			coeffs2.append( [sample_shift, cc_factor, cc_factor_abs, corr_array, rot_angle] )
		
		coeffs2          = np.array(coeffs2)
		arg2             = np.argmax( coeffs2[:,2] )	#sort after highest absolute cc_factor
		
		spl_shift_good2  = coeffs2[arg2,0]
		cc_good2         = coeffs2[arg2,1]
		cc_abs_good2     = coeffs2[arg2,2]
		corray_good2     = coeffs2[arg2,3]
		rot_angle_good2  = coeffs2[arg2,4]
		######

		# plotting
		#data2_rot_T_good, data2_rot_Q_good                          = rotate_2D(st_cross_T2.data, st_cross_Q2.data, rot_angle_good2, clockwise=True)	
		#st_cross_T2.data[:sample_shift]                             = 0
		#st_cross_T2.data[sample_shift+len(match_part):]             = 0
		st_cross_T.data[spl_shift_good:spl_shift_good+len(match_part)] = sin_good*np.max(np.abs(st_cross_Q.data[spl_shift_good:spl_shift_good+len(match_part)]))
		st_cross_T.stats.channel                                   += '_M'
		self.append(st_cross_T2)
		self.plot()

		split_time = spl_shift_good2*st_cross_T2.stats.delta - max_splittime
		split_dire = 4
		
		# output
		print()
		print('SPLITTING ANALYSIS (ROT-CORR):' )
		print('       Chosen window length : %.1f s'   % window_length)
		print('       f of BHQ matched sin : %.3f Hz (T=%.3f s)' % (freq_good, (1./freq_good)))
		print('         Q-T rotation angle : %.0f'     % rot_angle_good2)
		print('        Match of sin in BHQ : %.1f %%'  % (cc_good*100))
		print('  Match (abs) of sin in BHQ : %.1f %%'  % (cc_abs_good*100))
		print('        Match of BHQ in BHT : %.1f %%'  % (cc_good2*100))
		print('  Match (abs) of BHQ in BHT : %.1f %%'  % (cc_abs_good2*100))
		print('                 Split time : %.2f s'   % abs(split_time))
	def _spectrogram(self, data, samp_rate, per_lap=0.9, wlen=None, log=False, outfile=None, fmt=None, axes=None, dbscale=False, mult=8.0, cmap=obspy_sequential, zorder=None, title=None, show=True, sphinx=False, normalize=True, clip=[0.0, 1.0]):

		"""
		Slightly modified from ObsPy.
		Computes and plots spectrogram of the input data.

			:param data: Input data
			:type samp_rate: float
			:param samp_rate: Samplerate in Hz
			:type per_lap: float
			:param per_lap: Percentage of overlap of sliding window, ranging from 0
			    to 1. High overlaps take a long time to compute.
			:type wlen: int or float
			:param wlen: Window length for fft in seconds. If this parameter is too
			    small, the calculation will take forever. If None, it defaults to
			    (samp_rate/100.0).
			:type log: bool
			:param log: Logarithmic frequency axis if True, linear frequency axis
			    otherwise.
			:type outfile: str
			:param outfile: String for the filename of output file, if None
			    interactive plotting is activated.
			:type fmt: str
			:param fmt: Format of image to save
			:type axes: :class:`matplotlib.axes.Axes`
			:param axes: Plot into given axes, this deactivates the fmt and
			    outfile option.
			:type dbscale: bool
			:param dbscale: If True 10 * log10 of color values is taken, if False the
			    sqrt is taken.
			:type mult: float
			:param mult: Pad zeros to length mult * wlen. This will make the
			    spectrogram smoother.
			:type cmap: :class:`matplotlib.colors.Colormap`
			:param cmap: Specify a custom colormap instance. If not specified, then the
			    default ObsPy sequential colormap is used.
			:type zorder: float
			:param zorder: Specify the zorder of the plot. Only of importance if other
			    plots in the same axes are executed.
			:type title: str
			:param title: Set the plot title
			:type show: bool
			:param show: Do not call `plt.show()` at end of routine. That way, further
			    modifications can be done to the figure before showing it.
			:type sphinx: bool
			:param sphinx: Internal flag used for API doc generation, default False
			:type clip: [float, float]
			:param clip: adjust colormap to clip at lower and/or upper end. The given
			    percentages of the amplitude range (linear or logarithmic depending
			    on option `dbscale`) are clipped.
		"""

		# enforce float for samp_rate
		samp_rate = float(samp_rate)

		# set wlen from samp_rate if not specified otherwise
		if not wlen:
			wlen = samp_rate / 100. 

		npts = len(data)
		# nfft needs to be an integer, otherwise a deprecation will be raised
		# XXX add condition for too many windows => calculation takes for ever
		nfft = int(next_power_of_2(int(wlen * samp_rate)))
		if nfft > npts:
		    nfft = int(next_power_of_2(npts / 8.0))

		if mult is not None:
		    mult = int(next_power_of_2(int(mult)))
		    mult = mult * nfft
		nlap = int(nfft * float(per_lap))

		data = data - data.mean()
		end = npts / samp_rate

		# Here we call not plt.specgram as this already produces a plot
		# matplotlib.mpl.specgram should be faster as it computes only the
		# arrays
		# XXX mlab.specgram uses fft, would be better and faster use rfft
		specgram, freq, time = mpl.mlab.specgram(data, Fs=samp_rate, NFFT=nfft, noverlap=nlap)
		if normalize:
			specgram = specgram/np.max(specgram)

		# db scale and remove zero/offset for amplitude
		if dbscale:
		    specgram = 10 * np.log10(specgram[1:, :])
		else:
		    specgram = np.sqrt(specgram[1:, :])
		freq = freq[1:]

		vmin, vmax = clip
		if vmin < 0 or vmax > 1 or vmin >= vmax:
		    msg = "Invalid parameters for clip option."
		    raise ValueError(msg)
		_range = float(specgram.max() - specgram.min())
		vmin = specgram.min() + vmin * _range
		vmax = specgram.min() + vmax * _range
		norm = Normalize(vmin, vmax, clip=True)

		if not axes:
		    fig = plt.figure()
		    ax = fig.add_subplot(111)
		else:
		    ax = axes

		# calculate half bin width
		halfbin_time = (time[1] - time[0]) / 2.0
		halfbin_freq = (freq[1] - freq[0]) / 2.0

		# argument None is not allowed for kwargs on matplotlib python 3.3
		kwargs = {k: v for k, v in (('cmap', cmap), ('zorder', zorder))
		          if v is not None}

		if log:
		    # pcolor expects one bin more at the right end
		    freq = np.concatenate((freq, [freq[-1] + 2 * halfbin_freq]))
		    time = np.concatenate((time, [time[-1] + 2 * halfbin_time]))
		    # center bin
		    time -= halfbin_time
		    freq -= halfbin_freq
		    # Log scaling for frequency values (y-axis)
		    ax.set_yscale('log')
		    # Plot times
		    ax.pcolormesh(time, freq, specgram, norm=norm, **kwargs)
		else:
		    # this method is much much faster!
		    specgram = np.flipud(specgram)
		    # center bin
		    extent = (time[0] - halfbin_time, time[-1] + halfbin_time,
		              freq[0] - halfbin_freq, freq[-1] + halfbin_freq)
		    ax.imshow(specgram, interpolation="nearest", extent=extent, **kwargs)

		# set correct way of axis, whitespace before and after with window
		# length
		ax.axis('tight')
		ax.set_xlim(0, end)
		ax.grid(False)

		if axes:
		    return ax

		ax.set_xlabel('Time [s]')
		ax.set_ylabel('Frequency [Hz]')
		if title:
		    ax.set_title(title)

		if not sphinx:
		    # ignoring all NumPy warnings during plot
		    with np.errstate(all='ignore'):
		        plt.draw()
		if outfile:
		    if fmt:
		        fig.savefig(outfile, format=fmt)
		    else:
		        fig.savefig(outfile)
		elif show:
		    plt.show()
		else:
		    return fig
	def obs_plot(self, **kwargs):
		# pass
		return super().plot(**kwargs)
	def obs_select(self, *args, **kwargs):
		#pass
		return super().select(*args, **kwargs)
	def select(self, *args, **kwargs):
		st_obs_select              = self.obs_select(*args, **kwargs)
		st_obs_select.origin       = self.origin   
		st_obs_select.original     = self.original
		st_obs_select.removed      = self.removed
		st_obs_select.gain_removed = self.gain_removed
		st_obs_select.unit         = self.unit
		st_obs_select.times        = self.times
		st_obs_select.inventory    = self.inventory
		return st_obs_select
	def plot_polarization(self, title='', return_results=False, show=True):
	
		"""
		Plot polarization of data.	
		"""

		### get data ready
		stream_N = self.select(component='N') or self.select(component='Q')  or self.select(component='U') or self.select(component='1') or self.select(component='T')
		stream_E = self.select(component='E') or self.select(component='T')  or self.select(component='V') or self.select(component='2') or self.select(component='R')
		stream_Z = self.select(component='Z') or self.select(component='L')  or self.select(component='W') or self.select(component='3')

		if not (stream_N and stream_E):
			print('')
			print('Return. No idea how to perform polarization on components: %s' % ', '.join([i.stats.channel for i in self]) )
			return

		# demeaning should be done
		trace_N  = stream_N[0].detrend('demean')
		trace_E  = stream_E[0].detrend('demean')
		if stream_Z:
			trace_Z = stream_Z[0].detrend('demean')
		else:
			trace_Z      = Trace()
			trace_Z.data = np.array([])

		# making sure for small numbers everything works out (precisely: arrow head of angle indicators)
		if max( [max(abs(trace_Z.data)), max(abs(trace_N.data)), max(abs(trace_E.data))] ) <= 1e-5:
			factor        = 1e9
			factor_str    = '%.0g' % (factor**-1)
			try:
				trace_Z.data *= factor
			except:
				pass
			trace_N.data *= factor
			trace_E.data *= factor

		else:
			factor        = 1
			factor_str    = ''

	
		### take care of labels
		label_R = '<--  %s  -->' % (trace_E.stats.channel[:-1] + 'R')
		label_N = '<--  %s  -->' % trace_N.stats.channel
		label_E = '<--  %s  -->' % trace_E.stats.channel


		### P-Pol 
		phi_2D, phi_3D, err_phi_2D, err_phi_3D, INCapp_2D, err_INCapp_2D, INCapp_3D, err_INCapp_3D, SNR_hori, SNR_2D_radZ, Rect_3D, Rect_H, Rect_RZ, comp_R, comp_T, eigvecs, eigvals = ppol_calc(trace_N.data, trace_E.data, trace_Z.data)


		### Printing to shell
		if trace_Z:
			print()
			print('P-POL CALCULATIONS (%s):' % title)
			print('   %11s : %-5.1f ± %3.1f' % ('baz 2D', phi_2D, err_phi_2D))
			print('   %11s : %-5.1f ± %3.1f' % ('baz 3D', phi_3D, err_phi_3D))
			print('   %11s : %-5.1f ± %3.1f' % ('incl. 2D', INCapp_2D, err_INCapp_2D))
			print('   %11s : %-5.1f ± %3.1f' % ('incl. 3D', INCapp_3D, err_INCapp_3D))
			print('   %11s : %-10.1f' % ('SNR hori', SNR_hori))
			print('   %11s : %-10.1f' % ('SNR radZ', SNR_2D_radZ))
			print('   %11s : %-10.2f' % ('Rect_hori', Rect_H))
			print('   %11s : %-10.2f' % ('Rect_radZ', Rect_RZ))
			print('   %11s : %-10.2f' % ('Rect_3D', Rect_3D))
	
			label_Z = '<--  %s  -->' % trace_Z.stats.channel
			ncols = 3
		else:
			print()
			print('P-POL CALCULATIONS (%s):' % title)
			print('   %11s : %-6.1f ± %3.1f' % ('baz 2D', phi_2D, err_phi_2D))
			print('   %11s : %-10.1f' % ('SNR hori', SNR_hori))
			print('   %11s : %-10.2f' % ('Rect_hori', Rect_H))
	
			ncols = 1


		### Rreturn 
		if return_results:
			return phi_2D, phi_3D, err_phi_2D, err_phi_3D, INCapp_2D, err_INCapp_2D, INCapp_3D, err_INCapp_3D, SNR_hori, SNR_2D_radZ, Rect_3D, Rect_H, Rect_RZ, comp_R, comp_T, eigvecs, eigvals


		### plot results
		if title in plt.get_figlabels():
			fig = plt.figure(title)
			fig.clf()
		fig     = plt.figure(figsize=plt.figaspect(1./ncols), num=title)
		colours = [[0, 0, 1-0.6*i/len(trace_E.data)] for i in range(len(trace_E.data))]
	
		maxi = max( list(np.abs(trace_N.data))+list(np.abs(trace_E.data)) ) * 1.05
		arr  = maxi/17
	
		ax   = fig.add_subplot(131)
		ax.grid(ls='-.', lw=0.5)
		ax.set(xlabel=label_E, ylabel=label_N)
		ax.set_ylim([-maxi, maxi])
		ax.set_xlim([-maxi, maxi])
		
		ax.spines['right'].set_color('none')
		ax.spines['top'].set_color('none')
		ax.spines['left'].set_color('none')
		ax.spines['bottom'].set_color('none')	
		
		xx      = np.linspace(-maxi, maxi, 100)
		yy      = 1/np.tan(phi_2D*M.pi/180) * xx
	
		# black lines
		ax.plot([-maxi, maxi], [0, 0], 'k', lw=1)
		ax.plot( [0, 0], [-maxi, maxi], 'k', lw=1)
		
		# data
		ax.scatter(trace_E.data, trace_N.data, s=7, c=colours)
		
		# arc for angle
		ax.plot( xx, yy, 'indianred', lw=1.5)
		ax.add_patch( Arc([0,0], maxi,  maxi, 90, -phi_2D, 0, color='indianred', lw=1.5, zorder=10 ))
		
		# arrows at end of arcs
		x  = maxi/2*M.sin( (phi_2D-7)*M.pi/180 )
		y  = maxi/2*M.cos( (phi_2D-7)*M.pi/180 )
		x2 = maxi/2*M.sin( (phi_2D+2)*M.pi/180 )
		y2 = maxi/2*M.cos( (phi_2D+2)*M.pi/180 )
		a  = FancyArrowPatch([x,y], [x2,y2], mutation_scale=20, lw=1.5, arrowstyle="-|>", color="indianred", zorder=15)
		ax.add_artist(a)		

		# others
		ax.set_aspect('equal', adjustable='box')
		ax.text(1, 0, title+' (%s points)' % len(trace_E.data), ha='right', color='red', transform=ax.transAxes, fontsize=6, bbox=dict(boxstyle='round,pad=0', fc='white', ec="white", lw=0))
		ax.text(0, 1, factor_str, ha='center', color='black', transform=ax.transAxes, fontsize=7, bbox=dict(boxstyle='round,pad=0', fc='white', ec="white", lw=0))
		
		# legend
		scatter_proxy = mpl.lines.Line2D([0],[0], c="indianred", marker='<')
		ax.legend([scatter_proxy], ['BAZ_2D=(%.1f\u00B1%.1f)\u00b0' % (phi_2D, err_phi_2D)], numpoints=1, loc='upper right', prop={'size': 8})
	
		
		if trace_Z:
			maxi2 = max( list(np.abs(comp_R))+list(np.abs(trace_Z.data)) ) * 1.05
			arr2  = maxi2/17
	
			ax = fig.add_subplot(132)
			ax.grid(ls='-.', lw=0.5)
			ax.set(xlabel=label_R, ylabel=label_Z)
			ax.set_ylim([-maxi2, maxi2])
			ax.set_xlim([-maxi2, maxi2])
		
			ax.spines['right'].set_color('none')
			ax.spines['top'].set_color('none')
			ax.spines['left'].set_color('none')
			ax.spines['bottom'].set_color('none')	
		
			xx2 = np.linspace(-maxi2/1.05*0.9, maxi2/1.05*0.9, 100)
			yy2 = 1/np.tan( INCapp_2D*M.pi/180) * xx2
		
			ax.plot([-maxi2, maxi2], [0, 0], 'k', lw=1)
			ax.plot( [0, 0], [-maxi2, maxi2], 'k', lw=1)
			ax.scatter(comp_R, trace_Z.data, s=5, c=colours)
			ax.plot( xx2, yy2, 'indianred', lw=1.5)
			ax.add_patch( Arc([0,0], maxi2,  maxi2, 90, -INCapp_2D, 0, color='indianred', lw=1.5, zorder=10 ))

			# arrows
			x  = maxi2/2*M.sin( (INCapp_2D-7)*M.pi/180 )
			y  = maxi2/2*M.cos( (INCapp_2D-7)*M.pi/180 )
			x2 = maxi2/2*M.sin( (INCapp_2D+2)*M.pi/180 )
			y2 = maxi2/2*M.cos( (INCapp_2D+2)*M.pi/180 )
			a  = FancyArrowPatch([x,y], [x2,y2], mutation_scale=20, lw=1.5, arrowstyle="-|>", color="indianred", zorder=15)
			ax.add_artist(a)

			# others
			ax.set_aspect('equal', adjustable='box')
			ax.text(1, 0, title+' (%s points)' % len(trace_E.data), horizontalalignment='right', color='red', transform=ax.transAxes, fontsize=6, bbox=dict(boxstyle='round,pad=0', fc='white', ec="white", lw=0))
			ax.text(0, 1, factor_str, ha='center', color='black', transform=ax.transAxes, fontsize=7, bbox=dict(boxstyle='round,pad=0', fc='white', ec="white", lw=0))
			
			#legend
			scatter_proxy = mpl.lines.Line2D([0],[0], c="indianred", marker='<')
			ax.legend([scatter_proxy], ['INCL_2D=(%.1f\u00B1%.1f)\u00b0' % (INCapp_2D, err_INCapp_2D)], loc='upper right', prop={'size': 8})
	

			## 3-D plot
			mean_E   = np.average(trace_E.data)
			mean_N   = np.average(trace_N.data)
			mean_Z   = np.average(trace_Z.data)
			maxi3    = max( list(np.abs(trace_E.data))+list(np.abs(trace_N.data))+list(np.abs(trace_Z.data)) ) * 1.05
			max_dist = M.sqrt(3*(maxi3*2)**2)

			ax = fig.add_subplot(133, projection='3d')
			ax.set_xlabel(trace_E.stats.channel[-1], labelpad=5, fontsize=11)
			ax.set_ylabel(trace_N.stats.channel[-1], labelpad=5, fontsize=11)
			ax.set_zlabel(trace_Z.stats.channel[-1], labelpad=5, fontsize=11)
			ax.scatter(trace_E.data, trace_N.data, trace_Z.data, c=colours, s=5, depthshade=False, zorder=-1)
			ax.scatter([mean_E], [mean_N], [mean_Z], s=20, c='white', edgecolors='k', depthshade=False, zorder=50)
			ax.set_xlim( [-maxi3,maxi3] )
			ax.set_ylim( [-maxi3,maxi3] )
			ax.set_zlim( [-maxi3,maxi3] )

			# eigenvectors
			eigvecs  = [ unit_vector(eigvecs[:,0]), unit_vector(eigvecs[:,1]), unit_vector(eigvecs[:,2]) ]
			for l, eigvec in enumerate( eigvecs ):
				length = (0.7-l*0.15) * max_dist/2
				a = Arrow3D(np.array( [mean_E,eigvec[1]] )*length,
							np.array( [mean_N,eigvec[0]] )*length,
					        np.array( [mean_Z,eigvec[2]] )*length,
							mutation_scale=20, lw=1.5, arrowstyle="-|>", color="indianred", zorder=20)
				ax.add_artist(a)

			# legend
			scatter_proxy = mpl.lines.Line2D([0],[0], c="indianred", marker='<')
			ax.legend([scatter_proxy], ['eigen vectors'], loc='upper right', prop={'size': 8})

		if show:
			plt.tight_layout()
			plt.show()
	def plot_spectrum(self, method='fft'):

		"""

		"""

		nrows     = len(self)
		fig, axes = plt.subplots(nrows=nrows, ncols=1, figsize=(5, nrows*3), num='Spectrum of Stream')

		if nrows < 2:
			axes = [axes]

		for i, tr in enumerate(self):
			freqs, ps = spectrum( tr.data, method=method, fs=tr.stats.sampling_rate, nfft=next_power_of_2(len(tr.data)) )
			axes[i].plot(np.abs(freqs), ps)
			axes[i].text(0.05, 0.95, tr.id, ha='left', color='k', transform=axes[i].transAxes, fontsize=7, bbox=dict(boxstyle='round,pad=0', fc='white', ec="white", lw=0))

		plt.tight_layout()
		plt.show()		
	def plot_spectro(self, times='utcdatetime', equal_scale=False, **spectro_kwargs):


		"""
		Take stream and plot first all seismograms, then respective spectrograms (spectrogram
		variables other than offered as variables are wrapped into "**spectro_kwargs" and 
		passed to ObsPy's spectrogram function).
	
		If times=relative, plot allows to zoom for all axes simultaneously, 
		if times=matplotlib, seismograms show dates and spectrogram seconds, so zoom won't
		work for both at the same time.
		
		Note that if you zoom, the sup-title will be updated with the new start- and endtimes 
		of the shown section.
		
	
		WARNING YLIM EXCEED minimum SAMPLE RATE ...
		dbscale
		when new parameter, keep xlims()
	
	
		POSSIBLE IMRPROVEMENTS:
		  - make zoom work also if you choose dates for seismogram, so spectrogram times
		   (seconds) would need to be converted somehow)
		  - when zooming, make each seismogram recentered on y-axis (canvas redraw?)
		"""

		try: spectro_kwargs['wlen']
		except: spectro_kwargs['wlen'] = 50

		try: spectro_kwargs['dbscale']
		except: spectro_kwargs['dbscale'] = False
	                                                                        
		try: spectro_kwargs['log']
		except: spectro_kwargs['log'] = False
	
		### callback and helpter functionsfunctions 
		def print_spectrogram_parameters():
			print()
			print('Spectrogram parameters:')
			for key in spectro_kwargs.keys():
				print('  %10s : %s' % (key, spectro_kwargs[key]))		
		def on_xlims_change(axes):		# update sub-title
			new_x = axes.get_xlim()
			fig.suptitle('%s - %s' % (mint+new_x[0],maxt+new_x[1]), fontsize=11)
		def on_click(evt):				# plot mouse position
			try:
				text = "mouse: x=%.2f,  y=%.2f" % (evt.xdata, evt.ydata)
				text = text + (35-len(text))*' '
			except:
				text = "mouse: None                                    "
			axes[0].text(0.99, 0.8, text, ha='right', color='red', transform=axes[0].transAxes, fontsize=7, bbox=dict(boxstyle='round,pad=0.4', fc='white', ec="white", lw=0.3))
			fig.canvas.draw()
			fig.canvas.flush_events()
		def new_zoom(evt):
			for ax in axes[:2*len(self)]:
				navmode = ax.get_navigate_mode()
				if navmode is not None:
					break
			if not navmode == 'ZOOM':
				return
			for i, ax in enumerate(axes[:2*len(self)]):
				ax.set_ylim( ylims[i] )
		def submit1(wlen):
			spectro_kwargs['wlen'] = float(wlen)
			text_box1.text_disp.set_color('red')
			_redraw()
		def submit2(ylim_min):
			#ylim[0] = float(ylim_min)
			text_box2.text_disp.set_color('red')
			#_redraw(ylim=ylim)
		def submit3(ylim_max):
			#ylim[1] = float(ylim_max)
			text_box3.text_disp.set_color('red')
			#_redraw(ylim=ylim)
		def submit4(dbscale):
			spectro_kwargs['dbscale'] = bool(dbscale)
			text_box4.text_disp.set_color('red')
			_redraw()
		def submit5(log):
			spectro_kwargs['log'] = bool(log)
			text_box5.text_disp.set_color('red')
			_redraw()
		def _redraw():
			pass
		#	print_spectrogram_parameters()
		#	for k in range(len(self)):
		#		axes[len(self)+k].clear()
		#		_spectro_axis( axes[len(self)+k], self[k], k, ylim=ylim, **spectro_kwargs )
		#	plt.draw()		
		def _spectro_axis(axis, trace, j, **spectro_kwargs):
			st=Stream2()
			st.append(trace)
			st._spectrogram(trace.data, trace.stats.sampling_rate, axes=axis, **spectro_kwargs)
			axis.yaxis.set_minor_formatter(plticker.FormatStrFormatter('%.2f'))
			axis.text(0.01, 0.8, trace.id, ha='left', color='orange', transform=axis.transAxes, fontsize=7, bbox=dict(boxstyle='round,pad=0.2', fc='none', ec="orange", lw=0.3))
			if j != len(self)-1:
				axis.get_xaxis().set_visible(False)
			else:
				axis.set_xlabel('Time (s)')
			if j == int(len(self)/2.):
				axis.yaxis.set_label_coords(-0.03, 0.5)
				axis.set_ylabel('Frequency (Hz)')
	
		### print the in ObsPy's stream object contained processing information
		#self.print_process(truncate=100)
		print_spectrogram_parameters()
		mint, maxt =  min([st.stats.starttime for st in self]), max([st.stats.endtime for st in self])
	
		### ONE figure for all
		fig        = plt.figure(figsize=(15,15))
		fig.suptitle('min %s - max %s' % (mint, maxt), fontsize=10)
		fig.canvas.mpl_connect('button_press_event', on_click)
		gs0        = gridspec.GridSpec(3, 1, height_ratios=[1, 1, 0.1])
		axes       = []
	
		### plot all seismograms
		gs00 = gridspec.GridSpecFromSubplotSpec(len(self), 1, subplot_spec=gs0[0], hspace=0)
		for i in range(len( self )):
	
			trace = self[i]
			if len(axes) > 0:
				if equal_scale:
					axes.append(fig.add_subplot(gs00[i], sharex=axes[-1], sharey=axes[-1]))
				else:
					axes.append(fig.add_subplot(gs00[i], sharex=axes[-1]))
			else:
				axes.append(fig.add_subplot(gs00[i]))
			axes[i].plot(trace.times(times), trace.data, lw=1)
			axes[i].margins(0)
			axes[i].text(0.01, 0.8, trace.id, ha='left', color='k', transform=axes[i].transAxes, fontsize=7, bbox=dict(boxstyle='round,pad=0.2', fc='none', ec="k", lw=0.3))
			#axes[i].callbacks.connect('xlim_changed', on_xlims_change)
			if i == int(len(self)/2.):
				axes[i].yaxis.set_label_coords(-0.03, 0.5)
				axes[i].set_ylabel('Amplitude')
			if i < len(self)-1:
				axes[i].get_xaxis().set_visible(False)
	
		### plot all spectrograms
		gs01 = gridspec.GridSpecFromSubplotSpec(len(self), 1, subplot_spec=gs0[1], hspace=0)
		for j in range(len( self )):
			trace = self[j]
			axes.append(fig.add_subplot(gs01[j]))
			#_spectro_axis( axes[j+i+1], self[j], j, **spectro_kwargs )
	
		### Textboxes to enter parameter
		gs02 = gridspec.GridSpecFromSubplotSpec(1, 7, subplot_spec=gs0[2], wspace=0.5, hspace=0)
		axes.append(fig.add_subplot(gs02[0]))
		axes.append(fig.add_subplot(gs02[1]))
		axes.append(fig.add_subplot(gs02[2]))
		axes.append(fig.add_subplot(gs02[3]))
		axes.append(fig.add_subplot(gs02[4]))
		axes.append(fig.add_subplot(gs02[5]))
		axes.append(fig.add_subplot(gs02[6]))
		for m in [-1, -2, -3, -4, -5, -6, -7]:
			axes[m].axis('off')
	
		#text_box1 = TextBox(axes[-6], 'wlen =', initial=str(spectro_kwargs['wlen']), color='1.00')
		#text_box1.on_submit(submit1)
		#text_box1.text_disp.set_color('red')
		#
		#text_box2 = TextBox(axes[-5], 'ylim min =', initial=str(min(ylim)), color='1.00')
		#text_box2.on_submit(submit2)
		#text_box2.text_disp.set_color('red')
		#
		#text_box3 = TextBox(axes[-4], 'ylim max =', initial=str(max(ylim)), color='1.00')
		#text_box3.on_submit(submit3)
		#text_box3.text_disp.set_color('red')
		#
		#text_box4 = TextBox(axes[-3], 'dbscale =', initial=str(spectro_kwargs['dbscale']), color='1.00')
		#text_box4.on_submit(submit4)
		#text_box4.text_disp.set_color('red')
		#
		#text_box5 = TextBox(axes[-2], 'log =', initial=str(spectro_kwargs['log']), color='1.00')
		#text_box5.on_submit(submit5)
		#text_box5.text_disp.set_color('red')
	
		### final statements
		gs0.tight_layout(fig, rect=[0, 0.03, 1, 0.95])
		ylims     = [ax.get_ylim() for ax in axes]
		plt.connect('button_release_event', new_zoom)
		plt.show()
	def plot(self, store_dir=pwd(), store_name='*', verticals=(), type='normal', method='full', save_and_no_show=False, xlim=[], ylim=[], **kwargs):

		"""
		Enhanced plot of stream as 'type', using ObsPy's plot function.
		  check: https://docs.obspy.org/packages/autogen/obspy.core.stream.Stream.plot.html
	  
		  'store_name' is the name of the plot when it shall be saved (default: .png) in 'store_dir'.
	  
		  'verticals' (a tuple) is needed to plot a vertical line in a specific axis of the seismogram.
		  Can be used for example to plot seismic phase arrival times or similar. The times, however,
		  must be given in seconds relative to starttime of axis. Example:
		  verticals = (
		              (time, label, colour, index),
		              (..), 
		              )
		  Note: - time must be passed as relative seconds with respect to axis beginning
		        - colour can be chosen from the matplotlib colours
		 		- index should be integer, but can be 'all' (to be plotted in all axes)
	  
		  If 'save_and_no_show'=True, no plot is shown but saved
		  if 'save_and_no_show'=False, plot is shown and further options will be displayed (e.g. saving)
	  
		  Following mouse and button events are available:
			- double left click, 's'         : save figure in pwd as `store_name`.png
			- double right click, 'q'  '     : close figure
			- left click                     : give one of two limits for polarization analysis
			- right click                    : delete all current picks done by left click
			- 'p', or enter                  : plot polarization if two picks have been chosen (left click)
			- 'C'                            : abord and exit python
			- '1' (or int number)            : filter seismogram (if wasn't filtered before), as specified in `self.filters` (see __init__)
		"""
	
		# helper functions
		def on_close(evt):
			if flags['saved']:
				print(u'Plot saved: %s' % outfile)
			else:
				pass
				#print(u'Plot not saved.')
		def on_click(evt):
			axis        = evt.inaxes
			try:
				stat_id = '.'.join( [c for c in axis.get_children() if isinstance(c, mpl.text.Text) ][0].get_text().split('.')[:3] )
				indexes = [i for i, ax in enumerate(axes) if stat_id in [c for c in ax.get_children() if isinstance(c, mpl.text.Text)][0].get_text() ]
			except AttributeError:
				#if clicking not on axes but on canvas ..
				return
	
			if not evt.dblclick:
				if axis.get_navigate_mode() == 'ZOOM':	
					return																		### single click

				try:
					mouse_clicks[stat_id]
				except:
					mouse_clicks[stat_id] = 0

				if evt.button==3:																# right-click
					for i in indexes:
						lines = axes[i].lines
						remove = []
						for l in range(len(lines)):
							xdata = lines[l].get_xdata()
							if xdata[0] == xdata[-1]:
								remove.append(l)
					
						remove.sort(reverse=True)
						print('remove', remove)
						for k in remove:
							lines[k].remove()
					mouse_clicks[stat_id] = 0


				if evt.button == 1: 								 							# left-click

					if mouse_clicks[stat_id]%2==0 and mouse_clicks[stat_id]>0:
						for i in indexes:
							axes[i].lines[-2].remove()
						mouse_clicks[stat_id] = 1

					mouse_clicks[stat_id] += 1	
					x = [evt.xdata, evt.xdata]
					for i in indexes:
						y = []
						for line in axes[i].lines:

							xdata = line.get_xdata()
							if xdata[0] == xdata[-1]:
								# avoid vertical lines that may be present already
								continue

							ymin, ymax = axes[i].get_ylim()
							y.append(ymin)
							y.append(ymax)

						margin = (max(y)-min(y))*1.05
						ylim   = [ min(y)-margin, max(y)+margin ]
						axes[i].plot(x, ylim, 'deepskyblue')
				fig.canvas.draw()
	

			else:																			### double click				
				if evt.button == 1:															# nothing right now
					pass
					#flags['saved'] = True
					#plt.savefig(outfile, dpi=200, bbox_inches='tight')
					#plt.close(fig)
				elif evt.button == 2:														# middle-click (set ylim for each axis to minimum of all axes and maximum of all axes)
					ylims = []
					for axis in axes:
						ylims.append( axis.get_ylim() )
					max_y = max( np.array(ylims)[:,0] )
					min_y = min( np.array(ylims)[:,1] )
					for axis in axes:
						axis.set_ylim( [max_y,min_y] )
					fig.canvas.draw()
				elif evt.button == 3:														# right-click (for each axis individually, set ylim to respective data minimum and maximum + margin)
					xlim   = axis.get_xlim()
					for axis in axes:
						y_mins = []
						y_maxs = []
						for line in axis.lines:
							x        = line.get_xdata()
							if x[0]==x[-1]:													# exclude vertical lines
								continue
							y        = line.get_ydata()
							i        = np.where( (x >= xlim[0]) & (x <= xlim[1]) )[0]		# all indexes of y_data according to xlim
							if not i.any():
								continue													# e.g. data gaps, that is, 2 lines within axes where one does not lie within chose xlims

							line_min = y[i].min()											# get minimum y within all data according to xlim
							line_max = y[i].max()											# get maximum y within all data according to xlim
							y_mins.append(line_min)
							y_maxs.append(line_max)
						
						y_min = min(y_mins)
						y_max = max(y_maxs)
						ylim  = [y_min-0.025*np.abs(y_max-y_min), y_max+0.025*np.abs(y_max-y_min)] # give 2.5% margin that works both for pos. & neg. values 
						axis.set_ylim( ylim )		
					fig.canvas.draw()
		def on_key(evt):
			if evt.inaxes:
				axis = evt.inaxes
			else:
				axis = axes[0]

			def instrument_removal(output, pre_filt=None, water_level=60): 
				xlim   = axis.get_xlim()
				output = output.upper()

				if self.current_filter_num == 0 and not self.removed:
					corr          = self.copy()
					corr.original = self.copy()
					corr.removed  = self.copy()
					corr.unit     = output

				elif self.current_filter_num == 0 and self.removed:
					corr          = self.removed.copy()
					corr.original = self.removed.copy()
					corr.removed  = None
					corr.unit     = 'RAW'
					corr.plot(xlim=xlim)
					return

				elif self.current_filter_num != 0 and not self.removed:
					corr          = self.original.copy()
					corr.original = self.original.copy()
					corr.removed  = self.original.copy()
					corr.unit     = output

				elif self.current_filter_num != 0 and self.removed:
					corr          = self.removed.copy()
					corr.original = self.removed.copy()
					corr.removed  = None
					corr.unit     = 'RAW'
					corr.filtering(self.current_filter_num)
					corr.plot(xlim=xlim)
					return


				if not mouse_clicks.keys():
					print()
					print('Please select at least one station location!')

				else:
					if not self.inventory:
						print('Reading inventory file for instrument response removal ..')
						self._set_inventory()
					

					components = self._get_components()
					for stat_id in mouse_clicks.keys():
						corr_part  = corr.select(id=stat_id+'*')
						corr_part.merge()

						corr_part2 = corr.original.select(id=stat_id+'*')
						corr_part2.merge()
				
						if re.search(r'[U|V|W]', components):
							corr_part.remove_response(inventory=self.inventory, output=output, pre_filt=pre_filt, water_level=water_level)
							corr_part2.remove_response(inventory=self.inventory, output=output, pre_filt=pre_filt, water_level=water_level)

						elif re.search(r'[Z|N|E]', components):
							corr_part  = rotate2VBBUVW(corr_part, inventory=self.inventory)
							corr_part.remove_response(inventory=self.inventory, output=output, pre_filt=pre_filt, water_level=water_level)
							corr_part.rotate('->ZNE', inventory=self.inventory, components='UVW')

							corr_part2 = rotate2VBBUVW(corr_part2, inventory=self.inventory)
							corr_part2.remove_response(inventory=self.inventory, output=output, pre_filt=pre_filt, water_level=water_level)
							corr_part2.rotate('->ZNE', inventory=self.inventory, components='UVW')

				
					corr.filtering(self.current_filter_num)
					corr.plot(xlim=xlim)


			if evt.key.lower()=='a':							# Acceleration as target unit after instrument response removal 
				pre_filt = (0.001, 0.002, 50, 60)
				instrument_removal('ACC')


			elif evt.key.lower()=='b':							# Before xlim were changed, meaning go back to full available data view 
				
				# set new xlim of all available data
				for axis in axes:
					x = []
					for line in axis.lines:
						x += list( line.get_xdata() )
					xlim = [min(x), max(x)]

				axis.set_xlim(xlim)

				# set best ylim for new xlim
				for axis in axes:
					y_mins = []
					y_maxs = []
					for line in axis.lines:
						x        = line.get_xdata()
						if x[0]==x[-1]:													# exclude vertical lines
							continue
						y        = line.get_ydata()
						i        = np.where( (x >= xlim[0]) & (x <= xlim[1]) )[0]		# all indexes of y_data according to xlim
						if not i.any():
							continue													# e.g. data gaps, that is, 2 lines within axes where one does not lie within chose xlims

						line_min = y[i].min()											# get minimum y within all data according to xlim
						line_max = y[i].max()											# get maximum y within all data according to xlim
						y_mins.append(line_min)
						y_maxs.append(line_max)
					
					y_min = min(y_mins)
					y_max = max(y_maxs)
					ylim  = [y_min-0.025*np.abs(y_max-y_min), y_max+0.025*np.abs(y_max-y_min)] # give 2.5% margin that works both for pos. & neg. values 
					axis.set_ylim( ylim )	
				
				plt.draw()


			elif evt.key.lower()=='c':							# Cease all 
				print('STOPPED ALL')
				plt.close('all')
				sys.exit()


			elif evt.key.lower()=='d':							# Displacement as target uni after instrument response removal 
				pre_filt = (0.001, 0.002, 50, 60)
				instrument_removal('DISP')


			elif evt.key.lower()=='g':							# Gain removal / add 
				self.gain_correction()
				self.plot(xlim=axis.get_xlim())


			elif evt.key.lower()=='h':							# Help summary 
				# TO BE IMRPOVED
				pass


			elif evt.key.lower()=='m':							# Narrow view 10% (zooming in) 
				xlim = axis.get_xlim()
				xlim_range = (xlim[1]-xlim[0])
				xlim = [xlim[0]+0.1*xlim_range, xlim[1]-0.1*xlim_range]
				for axis in axes:
					axis.set_xlim( xlim )
				plt.draw()	


			elif evt.key.lower()=='n':							# Move out 10% (zooming out) 
				xlim = axis.get_xlim()
				xlim_range = (xlim[1]-xlim[0])
				xlim = [xlim[0]-0.1*xlim_range, xlim[1]+0.1*xlim_range]
				for axis in axes:
					axis.set_xlim( xlim )
				plt.draw()	


			elif evt.key=='p':									# Polarisation plot 
				for l, stat_id in enumerate( mouse_clicks.keys() ):
					if mouse_clicks[stat_id] == 1:
						continue
					elif mouse_clicks[stat_id] != 2:
						print('Choose at least 2 time limits (left-click) per station to perform action.')
						continue

					if l==len(mouse_clicks.keys())-1:
						show=True
					else:
						show=False						
					indexes       = [i for i, ax in enumerate(axes) if stat_id in [c for c in ax.get_children() if isinstance(c, mpl.text.Text)][0].get_text() ]
					ax_1st        = axes[indexes[0]]
					xlim          = sorted( [ax_1st.lines[-1].get_xdata()[0], ax_1st.lines[-2].get_xdata()[0]] )
					if type == 'normal':
						xlim   = [ UTCDateTime(mdates.num2date(x)) for x in xlim ]
					elif type == 'relative':
						xlim = 	[ min_starttime+x for x in xlim ]
					
					polar = self.copy()
					polar = polar.select(id=stat_id+'*')
					polar.trim(starttime=xlim[0], endtime=xlim[1])
					polar.plot_polarization(title='%s, %s, %s, %s-%s'% (stat_id, self.unit, self._get_filter_str(), xlim[0].strftime('%H:%M:%S'), xlim[1].strftime('%H:%M:%S')), show=show)


			elif evt.key=='P':									# Polarisation plot as batch job 
				for l, stat_id in enumerate( mouse_clicks.keys() ):
					if mouse_clicks[stat_id] == 1:
						continue
					elif mouse_clicks[stat_id] != 2:
						print('Choose at least 2 time limits (left-click) per station to perform action.')
						continue
					print()
					print('Running P-pol batch job ..')

					if l==len(mouse_clicks.keys())-1:
						show=True
					else:
						show=False
					indexes       = [i for i, ax in enumerate(axes) if stat_id in [c for c in ax.get_children() if isinstance(c, mpl.text.Text)][0].get_text() ]
					ax_1st        = axes[indexes[0]]
					xlim          = sorted( [ax_1st.lines[-1].get_xdata()[0], ax_1st.lines[-2].get_xdata()[0]] )
					if type == 'normal':
						xlim   = [ UTCDateTime(mdates.num2date(x)) for x in xlim ]
					elif type == 'relative':
						xlim = 	[ min_starttime+x for x in xlim ]
					
					# create windows for batch job
					xlims    = []
					span     = xlim[1]-xlim[0]
					pro_side = 2
					fraction = 0.15
					for l in np.arange(-pro_side, pro_side+1):
						for i in np.arange(-pro_side, pro_side+1):
							window = [ xlim[0]+l*span*fraction, xlim[1]+i*span*fraction]
							xlims.append( window )

					# iterate over windows and filters to calculate polarisations
					results = []
					unit    = self.unit
					for xlim in xlims:
						for fil_num in self.filters.keys():
							polar   = self.original or self
							polar   = polar.copy()
							polar   = polar.select(id=stat_id+'*')
							polar.filtering(fil_num)
							polar.trim(starttime=xlim[0], endtime=xlim[1])
							
							fil_str = polar._get_filter_str()
							window  = '%s-%s' % (xlim[0].strftime('%H:%M:%S'), xlim[1].strftime('%H:%M:%S'))
							title   = '%s, %s, %s, %s' % (stat_id, unit, fil_str, window)						
							result  = polar.plot_polarization(title=title, return_results=True)
							results.append( list(result)+xlim+[fil_num] )

					# retrieve needed variables for best polarisation
					results      = np.array(results)
					ind_high     = np.argmax(results[:,11])		# highest Rect_H
					rect_Hs      = list( results[:,11] )
					fil_num_high = results[ind_high,-1] 
					xlim_high    = results[ind_high,-3:-1] 

					# plot best polarisationg
					polar_best   = self.copy()
					polar_best   = polar_best.select(id=stat_id+'*')
					polar_best.filtering(fil_num_high)
					polar_best.trim(starttime=xlim_high[0], endtime=xlim_high[1])

					fil_str      = polar_best._get_filter_str()
					window       = '%s-%s' % (xlim_high[0].strftime('%H:%M:%S'), xlim_high[1].strftime('%H:%M:%S'))
					title        = '%s, %s, %s, %s' % (stat_id, unit, fil_str, window)

					print('')
					print('B E S T')
					print('-' * 40)
					polar_best.plot_polarization(title='Best: %s' % title, show=False)
					print('-' * 40)

					# plot histogram
					hist_title = 'Histogram P-pol batch job (%s)' % self.unit
					if title in plt.get_figlabels():
						fig = plt.figure(hist_title)
						fig.clf()
					fig = plt.figure(num=hist_title)
					fig.suptitle('Histogram of %s P-wave polarisation calculations' % (len(rect_Hs)), fontsize=9)
					ax = fig.add_subplot(111)
					ax.grid(ls='-.', lw=0.5, zorder=0)
					ax.set(xlabel='Horizontal rectilinearity', ylabel='Number of occurences')
					ax.hist(rect_Hs, bins=100, range=(0,1), rwidth=1, ec='k', lw=0.5, zorder=3)
					if show:
						plt.show()


			elif evt.key.lower()=='q':							# Quit current figures
				# close
				print('Closed current figure window')
				plt.gcf()
				plt.close()


			elif evt.key.lower()=='r':							# Rotate data 
				
				self.merge(method=1)
				components = self._get_components()


				## rotate data to ZNE
				if re.search(r'[U|V|W]', components):

					# correct gain
					if not self.gain_removed:
						self.gain_correction()

						if self.original:
							self.original.gain_correction()
							self.original.gain_removed = True

						if self.removed:
							self.removed.gain_correction()
							self.removed.gain_removed = True

					# rotate
					self.rotate('->ZNE', inventory=self.inventory, components='UVW')
					
					if self.original:
						self.original.rotate('->ZNE', inventory=self.inventory, components='UVW')

					if self.removed:
						self.removed.rotate('->ZNE', inventory=self.inventory, components='UVW')

				
				## rotate data to UVW
				elif re.search(r'[Z|N|E]', components):
					rotate2VBBUVW(self, inventory=self.inventory)

					if self.original:
						rotate2VBBUVW(self.original, inventory=self.inventory)

					if self.removed:
						rotate2VBBUVW(self.removed, inventory=self.inventory)

				self.plot(xlim=axis.get_xlim())


			elif evt.key=='t':									# Trim file and save original data in viewed x-lim bounds 
				for stat_id in mouse_clicks.keys():

					indexes       = [i for i, ax in enumerate(axes) if stat_id in [c for c in ax.get_children() if isinstance(c, mpl.text.Text)][0].get_text() ]
					ax_1st        = axes[indexes[0]]
					
					if mouse_clicks[stat_id]==0:
						xlim = ax_1st.get_xlim()

					elif mouse_clicks[stat_id]==2:
						xlim = sorted( [ax_1st.lines[-1].get_xdata()[0], ax_1st.lines[-2].get_xdata()[0]] )

					else:
						print('Choose None or 2 time limits (left-click) per station to trim data & save data.')
						continue

					
					if type == 'normal':
						xlim   = [ UTCDateTime(mdates.num2date(x)) for x in xlim ]
					elif type == 'relative':
						xlim = 	[ min_starttime+x for x in xlim ]

					def submit(filename):
						if self.original:
							trim = self.original.copy()
						else:
							trim = self.copy()
						trim = trim.select(id=stat_id+'*')
						trim.trim(starttime=xlim[0], endtime=xlim[1])
						trim.write(filename)
						plt.close()
					
					fig      = plt.figure(figsize=(17,0.7), num='Save trimmed file')
					ax       = fig.add_subplot(111)
					
					propose  = os.path.join( os.path.dirname(self.origin), 'TRIMORIG_'+os.path.basename(self.origin) )
					text_box = TextBox(ax, 'Filename:', initial=propose)
					text_box.on_submit(submit)
					plt.show()


			elif evt.key=='T':									# Trim file and save present data in viewed x-lim bbounds 
				for stat_id in mouse_clicks.keys():

					indexes       = [i for i, ax in enumerate(axes) if stat_id in [c for c in ax.get_children() if isinstance(c, mpl.text.Text)][0].get_text() ]
					ax_1st        = axes[indexes[0]]
					
					if mouse_clicks[stat_id]==0:
						xlim = ax_1st.get_xlim()

					elif mouse_clicks[stat_id]==2:
						xlim = sorted( [ax_1st.lines[-1].get_xdata()[0], ax_1st.lines[-2].get_xdata()[0]] )

					else:
						print('Choose None or 2 time limits (left-click) per station to trim data & save data.')
						continue

					
					if type == 'normal':
						xlim   = [ UTCDateTime(mdates.num2date(x)) for x in xlim ]
					elif type == 'relative':
						xlim = 	[ min_starttime+x for x in xlim ]

					def submit(filename):
						trim = self.copy()
						trim = trim.select(id=stat_id+'*')
						trim.trim(starttime=xlim[0], endtime=xlim[1])
						trim.write(filename)
						plt.close()
					
					fig      = plt.figure(figsize=(17,0.7), num='Save trimmed file')
					ax       = fig.add_subplot(111)
					
					propose  = os.path.join( os.path.dirname(self.origin), 'TRIMPRES_'+os.path.basename(self.origin) )
					text_box = TextBox(ax, 'Filename:', initial=propose)
					text_box.on_submit(submit)
					plt.show()


			elif evt.key.lower()=='u':							# Update plot as to chosen xlim boundaries 
				for stat_id in mouse_clicks.keys():
					if mouse_clicks[stat_id] == 1:
						continue
					elif mouse_clicks[stat_id] == 2 :
						indexes       = [i for i, ax in enumerate(axes) if stat_id in [c for c in ax.get_children() if isinstance(c, mpl.text.Text)][0].get_text() ]
						ax_1st        = axes[indexes[0]]
						xlim          = sorted( [ax_1st.lines[-1].get_xdata()[0], ax_1st.lines[-2].get_xdata()[0]] )			
					else:
						xlim          = axes[0].get_xlim()

					for axis in axes:
						y_mins = []
						y_maxs = []
						for line in axis.lines:
							x        = line.get_xdata()
							if x[0]==x[-1]:													# exclude vertical lines
								continue
							y        = line.get_ydata()
							i        = np.where( (x >= xlim[0]) & (x <= xlim[1]) )[0]		# all indexes of y_data according to xlim
							if not i.any():
								continue													# e.g. data gaps, that is, 2 lines within axes where one does not lie within chose xlims
							
							line_min = y[i].min()											# get minimum y within all data according to xlim
							line_max = y[i].max()											# get maximum y within all data according to xlim
							y_mins.append(line_min)
							y_maxs.append(line_max)
			
						y_min = min(y_mins)
						y_max = max(y_maxs)
						ylim  = [y_min-0.025*np.abs(y_max-y_min), y_max+0.025*np.abs(y_max-y_min)] # give 2.5% margin that works both for pos. & neg. values 
						axis.set_ylim( ylim )		
						axis.set_xlim( xlim )

					plt.draw()		


			elif evt.key.lower()=='v':							# Velocity as target uni after instrument response removal 
				pre_filt = (0.001, 0.002, 50, 60)
				instrument_removal('VEL')


			elif evt.key.lower()=='w':							# Window closing (all) 
				# close
				print('Close all figure windows')
				plt.close('all')
			

			elif evt.key.lower()=='z':							# Demean current window viewed
				xlim = axis.get_xlim()

				for axis in axes:
					y_mins = []
					y_maxs = []
					for line in axis.lines:
						
						x_data = line.get_xdata()
						if not x_data[0] == x_data[-1]: #exclude vertical lines)

							# set ylim to best ylim of current viewing window
							i = np.where( (x_data >= xlim[0]) & (x_data <= xlim[1]) )[0]
							if not i.any():
								continue

							# demean window
							y_data = line.get_ydata()
							y_data = y_data - y_data[i].mean()
							line.set_ydata( y_data ) 

							line_min = y_data[i].min()
							line_max = y_data[i].max()
							y_mins.append(line_min)
							y_maxs.append(line_max)							
					y_min = min(y_mins)
					y_max = max(y_maxs)
					ylim  = [y_min-0.025*np.abs(y_max-y_min), y_max+0.025*np.abs(y_max-y_min)] # give 2.5% margin that works both for pos. & neg. values 
					axis.set_ylim( ylim )

				print(u'Demeaned current window.')	
				plt.draw()


			elif evt.key.lower()==',':							# Move 10% left 
				xlim = axis.get_xlim()
				xlim = xlim - 0.1*(xlim[1]-xlim[0])
				for axis in axes:
					axis.set_xlim( xlim )
				plt.draw()	


			elif evt.key.lower()=='.':							# Move 10% right 
				xlim = axis.get_xlim()
				xlim = xlim + 0.1*(xlim[1]-xlim[0])
				for axis in axes:
					axis.set_xlim( xlim )
				plt.draw()	


			elif evt.key.lower()=='up':							# Reduce ylim 10% 
				for axis in axes:
					ylim = axis.get_ylim()
					ylim = [ylim[0]+0.1*(ylim[1]-ylim[0]), ylim[1]-0.1*(ylim[1]-ylim[0])]
					axis.set_ylim( ylim )
				plt.draw()	


			elif evt.key.lower()=='down':						# Increase ylim 10% 
				for axis in axes:
					ylim = axis.get_ylim()
					ylim = [ylim[0]-0.1*(ylim[1]-ylim[0]), ylim[1]+0.1*(ylim[1]-ylim[0])]
					axis.set_ylim( ylim )
				plt.draw()	


			elif evt.key.isdigit():								# filter data, depending on choice 
				xlim = axis.get_xlim()

				if self._get_filter_str() == 'no filter':
					filt              = self.copy()
					filt.original     = self.copy()

				else:
					filt              = self.original.copy()
					filt.original     = self.original.copy()
				
				filt.removed = self.removed
				filt.filtering(evt.key)
				filt.plot(xlim=xlim)	
	

		# variables
		mouse_clicks  = {}
		flags         = {'saved' : False}
		colours       = 'bgcm'
		outfile       = os.path.join(store_dir, store_name+'.png')	# change format if desired (used if plot shall be saved)
		title         = '%s, %s' % (self.unit, self._get_filter_str())
		min_starttime = min( np.array( [tr.stats.starttime for tr in self] ))
	
	
		# plotting seismogram
		figs = [[manager.num, manager.canvas.figure.canvas.get_window_title()] for manager in mpl._pylab_helpers.Gcf.get_all_fig_managers()]
		for fignum, figtit in figs:
			if title == figtit:
				fig = mpl.pyplot.figure(fignum)
				plt.close()

		fig = self.obs_plot(type=type, show=False, handle=True, equal_scale=False, method='full', **kwargs)
		fig.canvas.draw()
		fig.canvas.set_window_title(title)
		fig.canvas.mpl_connect('button_press_event', on_click)
		fig.canvas.mpl_connect('key_press_event', on_key)
		fig.canvas.mpl_connect('close_event', on_close)
		fig.suptitle('Min starttime: %s (julday %d)' % (min_starttime.strftime('%Y-%m-%dT%H:%M:%S.%f').rstrip("0").rstrip("."), min_starttime.julday), fontsize=9)
		
		axes = fig.axes
	

		# loop over axes to add certain stuff
		for i in range(len(axes)):			
			if verticals:
				x_startpoint = axes[i].lines[0].get_xdata()[0]
				for time, label, colour, index in verticals:
					time = float(time)
					if type == 'normal':
						v_line_y     = UTCDateTime( mdates.num2date(x_startpoint) )
						v_line_y     = (v_line_y+time).datetime
					elif type == 'relative':
						v_line_y     = x_startpoint+time
					if index == 'all' or float(index)==i:
						axes[i].plot([v_line_y,v_line_y], axes[i].get_ylim(), color=colour, lw=1.0)
	
			if i == 0:
				axes[i].text(0.99, 0.92, self._get_filter_str(), horizontalalignment='right', color='red', transform=axes[i].transAxes, fontsize=8, bbox=dict(boxstyle='round,pad=0.2', fc='white', ec="red", lw=0.5))
	
			if i+1 == len(self):
				axes[i].text(0.99, 0.95, "Close all windows: 'w'.", horizontalalignment='right', color='red', transform=axes[i].transAxes, fontsize=7)
				if verticals:
					verts = np.unique( np.array(verticals)[:,:-1], axis=0 )
					verts = verts[verts[:,0].astype('float').argsort()]
					for j, verts_uni in enumerate(verts):
						axes[i].text(0.99, 0.95-j*0.06, "%s" % verts_uni[1], horizontalalignment='right', color=verts_uni[2], transform=axes[i].transAxes, fontsize=6)




		# set lims
		if list( xlim ):
			for axis in axes:
				axis.set_xlim(xlim)
		else:
			xlim = axes[0].get_xlim()
		if list( ylim ):
			for k, axis in enumerate(axes):
				if len(ylim)==len(axes) and not [y for y in ylim if not isinstance(y, (list, tuple, np.ndarray))]:
					axis.set_ylim( ylim[k] )
				else:
					axis.set_ylim(ylim)
		else:
			for axis in axes:
				y_mins = []
				y_maxs = []
				for line in axis.lines:
					x        = line.get_xdata()
					if x[0]==x[-1]:													# exclude vertical lines
						continue
					y        = line.get_ydata()
					i        = np.where( (x >= xlim[0]) & (x <= xlim[1]) )[0]		# all indexes of y_data according to xlim
					if not i.any():
						continue													# e.g. data gaps, that is, 2 lines within axes where one does not lie within chose xlims

					line_min = y[i].min()											# get minimum y within all data according to xlim
					line_max = y[i].max()											# get maximum y within all data according to xlim
					y_mins.append(line_min)
					y_maxs.append(line_max)
				
				y_min = min(y_mins)
				y_max = max(y_maxs)
				ylim  = [y_min-0.025*np.abs(y_max-y_min), y_max+0.025*np.abs(y_max-y_min)] # give 2.5% margin that works both for pos. & neg. values 
				axis.set_ylim( ylim )	


		# show / save figure
		if save_and_no_show:
			flags['saved'] = True
			plt.savefig(outfile, dpi=200, bbox_inches='tight')
			plt.close(fig)
		else:
			# show plot
			plt.show() 
def folder2sds(main_folder='/Users/jscholz/data/KNIPA_SUBSET', pattern='*.mseed', invers=False, seis_format=None):

	"""
	This script sorts all files within 'main_folder' (recursive file search) that contain 'pattern' with their actual names (not path to it)
	into seiscomp data structure (SDS), using the header information found within the files (so not based on actual file names).
	To read the accordant header information, ObsPy's 'read' function is used.

	Note: 
		- only seismlogical files of 'seis_format' are considered, e.g. SAC, MSEED, ... (None means no specification)
		- seismoglogical files should only contain one component (e.g. "BHZ", not multiple components)
		- non seis-files, files where start & end year mismatch, and files where channels within file have contradicting header information (e.g., different networks ..) are skipped
		- 'invers'==True means that all files from SDS (or in fact any random folder structure) are copied into 'main_folder' with no sub-directories.
		- files are only shifted into accordant folder structure, file names themselves don't change
		- folders that don't exist yet will be silently created (if permissions are sufficient, else error is thrown which is not caught by this script)
		- if a file already exists in the target (sub-)folder, then an errow is thrown that is not cauht by this script
		- empty folders remaining after re-organisation are removed
	"""

	now   = time.time()
	files = get_files(main_folder, pattern=pattern)
	no_mv = 0


	# loop over all files found under 'main_folder' (recursively and in combination with 'pattern')
	for file in files:


		# all files within main folder (also sub-dirs) organised into seiscomp strcture using seis file header information
		if not invers:
			try:
				st = read(file, headonly=True, format=seis_format)
			except Exception as err:
				print(u'Skipped, as likely no seismological file. File:  %s' % file)
				no_mv += 1
				continue

			# loop over each channel within stream object (can be same component but for instance multiple channels due to data gaps ..)
			networks, stations, channels = [], [], []
			for ch in st:

				start_year     = ch.stats.starttime.year
				end_year       = ch.stats.endtime.year
				if start_year != end_year:
					print(u'Skipped, as start & end year mismatch. File:     %s' % file)
					no_mv += 1
					break

				network        = ch.stats.network
				station        = ch.stats.station
				channel        = ch.stats.channel

				networks.append( network )
				stations.append( station )
				channels.append( channel )

			else:																						# executed if break statement within for-loop is not encountered (not exectued if break is encountered)

				if len(set(networks)) != 1 or len(set(stations)) != 1 or len(set(channels)) != 1:		# if all channels in one stream don't have same header information, they need to be treated manually
					print(u'Skipped, as different channels in file have contradicting header information. File:  %s' % file)
					no_mv += 1
					continue

				target_file    = os.path.join(main_folder, str(start_year), str(network), str(station), str(channel)+'.D', os.path.basename(file) )


		# SDS (or any other folder structure) into main folder without sub-directory structure
		else:
			target_file    = os.path.join(main_folder, os.path.basename(file) )


		# actual re-organisation done. Should also remove folders once they are empty
		os.renames(file, target_file)


	### OUTPUT
	print(u'')
	print(u'Finished "folder2sds()".')
	print(u'Re-organised %s/%s files!' % (len(files)-no_mv,len(files)) )
	if not invers:
		print(u'main folder (or random) --> seiscomp')
	if invers:
		print(u'seiscomp (or random) --> main folder')
	print(u"Main folder:                   '%s'"       % main_folder )
	print(u"File pattern:                  '%s'"       % pattern )
	print(u"Format seismological file:     '%s'"       % seis_format )
	print(u"Done in:                       %s (h:m:s)" % sec2hms( float( time.time()-now )) )
	print(u"Timestamp:                     %s"         % datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d, %H:%M:%S'))
def make_symlinks_to_data(target_dir=u'/Users/jscholz/data/KNIPA_SUBSET_MARCH', pattern=u'*', logfile=None, overwrite=True):

	"""
	Visualising data, e.g. with Pyrocko's 'snuffler' tool, can take tremendously long for large data sets.
	50.000 to 100.00 files are reasonably fast displayed, millions of files take just too long,
	making it impractical to work with.

	One work around is to create a new folder with symlinks to actual data of interest.
	This script does this. Just use (bash) wildcards for 'seiscomp_base_dir' and 'pattern',
	for these files symlinks are created in 'target_dir'. After that, you may do:

		snuffler target_dir

	And there you have your subset of data visualised without copying actual data. 
	If overwrite=True, target files that already exist will be deleted.
	If overwrite=False, target files that already exist cause a warning &
	                    no symlink creation.

	ATTENTION:
		- for some reason, symlinks that are stored somewhere on e.g. seisrv3 on not visible if 
		  an "ls" is done via the terminal if you just have mounted the accordant filesystem (also not via finder then)!
		  You actually need to connect to seissrv3, and do a "ls" an the specific folder the
		  symlinks are stored. No idea why ...
		- symlinks are named as original files, just stored in 'target_dir'
	"""


	now         = time.time()
	#pattern=u'1B.KNR??..???.?.20??-??-??_*.mseed'
	seiscomp_base_dir=[u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2016/1B/KNR01',
					   u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2016/1B/KNR02',
					   u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2016/1B/KNR03',
					   u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2016/1B/KNR04',
					   u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2016/1B/KNR05',
					   u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2016/1B/KNR06',
					   u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2016/1B/KNR07',
					   u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2016/1B/KNR08',
					   u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2016/1B/KNR09',
					   u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2016/1B/KNR10',
					   u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2016/1B/KNR11',
					   u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2016/1B/KNR12',
					   u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2016/1B/KNR13',
					   u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2016/1B/KNR14',
					   u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2016/1B/KNR15',
					   u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2016/1B/KNR16',
					   u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2016/1B/KNR17',
					   u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2016/1B/KNR18',
					   u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2016/1B/KNR19',
					   u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2016/1B/KNR20',
					   u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2016/1B/KNR21',
					   u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2016/1B/KNR22',
					   u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2016/1B/KNR23',
					   u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2016/1B/KNR28',
					   u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2016/1B/KNR29',
					   u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2016/1B/KNR30',
					   u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2017/1B/KNR01',
					   u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2017/1B/KNR02',
					   u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2017/1B/KNR03',
					   u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2017/1B/KNR04',
					   u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2017/1B/KNR05',
					   u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2017/1B/KNR06',
					   u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2017/1B/KNR07',
					   u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2017/1B/KNR08',
					   u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2017/1B/KNR09',
					   u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2017/1B/KNR10',
					   u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2017/1B/KNR11',
					   u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2017/1B/KNR12',
					   u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2017/1B/KNR13',
					   u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2017/1B/KNR14',
					   u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2017/1B/KNR15',
					   u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2017/1B/KNR16',
					   u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2017/1B/KNR17',
					   u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2017/1B/KNR18',
					   u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2017/1B/KNR19',
					   u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2017/1B/KNR20',
					   u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2017/1B/KNR21',
					   u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2017/1B/KNR22',
					   u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2017/1B/KNR23',
					   u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2017/1B/KNR28',
					   u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2017/1B/KNR29',
					   u'/Users/jscholz/data/KNIPA_NEW_SPLIT3600_DOWNSAMPLED50/2017/1B/KNR30',
					   ]


	### if logfile wished, do so
	if logfile:
		sys.stdout = Logger(logfile)


	### if target_dir not existing, create
	if not os.path.isdir(target_dir):
		os.makedirs( target_dir )


	### loop over all seis_comp_dir
	for dir in seiscomp_base_dir:

		source_files  = get_files(dir, pattern=pattern)
		created_links = 0
		if not source_files:
			print(u'Not ok! No files found matching given pattern "%s" in data directory: %s' % (pattern,dir))
			continue

		for source_file in source_files:
			target_file = os.path.join( target_dir,os.path.basename(source_file) )

			try:
				os.symlink(source_file, target_file)
				created_links += 1
			except OSError as err:
				if overwrite:
					os.remove( target_file )
					os.symlink(source_file, target_file)
					created_links += 1
				else:
					print(u'Skip, because file exists already: %s' % target_file)
					print(u'Check option: overwrite=True')
		print(u'OK! Found %s files with pattern "%s", made %s file-symlinks to data: %s' % (created_links,pattern,len(source_files),dir))


	### Final output
	print(u'Finished "make_symlinks_to_data()".')
	print(u'Target dir: %s '        % target_dir)
	print(u'Time:       %s (h:m:s)' % sec2hms( float( time.time()-now )) )
	print(u'Timestamp:  %s'         % datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d, %H:%M:%S'))
	if logfile:
		print(u'Logfile:    %s'     % logfile)
		sys.stdout = sys.__stdout__											# reset stdout to normal 
def delete_files(base_dir=u'/Users/jscholz/data/KNIPA_SUBSET', patterns=[]):

	"""
	Delete files found under 'base_dir', matching 'patterns'.
	E.g., files with certain station components that shall not be considered 
	after symlinks have been created..
	"""

	now         = time.time()
	if not patterns:
		patterns = [u'1B.KNR08..?H1.*.mseed',
		            u'1B.KNR09..?HZ.*.mseed',
		            u'1B.KNR13..?H1.*.mseed',
		            u'1B.KNR14..?HZ.*.mseed',
		            u'1B.KNR15..?H?.*.mseed',
		            u'1B.KNR20..?H2.*.mseed',
		            u'1B.KNR24..???.*.mseed',
		            u'1B.KNR25..???.*.mseed',
		            u'1B.KNR26..???.*.mseed',
		            u'1B.KNR27..???.*.mseed',
		            u'1B.KNR28..?DH.*.mseed',
		            u'1B.KNR29..?DH.*.mseed',
		            u'1B.KNR30..?DH.*.mseed',
		            ]

	total       = 0
	for pattern in patterns:
		files   = get_files(base_dir, pattern=pattern)

		for file in files:
			os.remove(file)
			total += 1

	print(u'')
	print(u'Finished "delete_files()". Deleted %s files in directory %s.' % (total,base_dir) )
	print(u'Time:       %s (h:m:s)'      % sec2hms( float( time.time()-now )) )
	print(u'Timestamp:  %s'              % datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d, %H:%M:%S'))
def diff_stream(stream1, stream2):

	"""
	Calculate difference of two stream objects.
	Only stream1 traces are subtracted by accordant stream2 traces 
	if they appear in both streams. Criteria is component letter.

	For stream1 traces not appearing in stream2, nothing is done
	and the returned stream contains these original traces.
	"""

	st_diff = stream1.copy()
	comps1  = ''.join( [tr.stats.channel[-1] for tr in stream1] )

	for comp in comps1:
		st_diff = stream1.select(component=comp)
		st_tmp1 = stream1.select(component=comp)
		st_tmp2 = stream2.select(component=comp)

		if st_tmp2:
			st_diff[0].data = st_tmp1[0].data - st_tmp2[0].data

	return st_diff
def get_station_metas(service='RESIF', network='YV', station='*', level='channel'):

	""" 
	Retrieve meda data of stations. Can read inventory files of catalog or
	STATIONXML files of single station. 
	"""

	###### shell style
	print(75*'-' + ' GO')
	print('Exec:\t\t' + sys._getframe().f_code.co_name + '()')
	print('..\n\n')

	###### retrieve actual data (inventory object), save inventory
	inventory = Client( service ).get_stations(network=network, station=station,level=level)
	#inventory.write('/home/john/Dropbox/JOHNS_WORK/programming/python/obspyDMT/obs_orientation/RR*_INSU_inventory.xml', 'STATIONXML')


	##### output
	counter = 0
	for network in inventory:
		print('\033[91m%s\033[0m === %s' % (network.code.upper(), network.description))

		for station in network:
			counter += 1
			print()
			print(str( counter) + '. ' + str( station.code.upper() + ':' ).ljust(6) + \
                  str( station.start_date )[:-8].ljust(21) + \
                  str( station.end_date )[:-8].ljust(21)+   \
                  "{:.4f}".format( station.longitude ).ljust(11) + \
                  "{:.4f}".format( station.latitude ).zfill(8).ljust(11) + \
                  "{:<7}".format( station.elevation ))
			print('cha   start                end                   sr       azi  dip')
			print(75 * '-')

			# sort the channels after channel name and then after time,
			# starting with newest
			k = [[channel, channel.code, channel.start_date] for channel in station]
			k.sort(key=lambda x: x[1]+str(x[2]), reverse=True)
			channels = [channel[0] for channel in k]
			station_sort = Station(station.code, station.latitude, \
                                   station.longitude, station.elevation, \
                                   channels=channels)
			for channel in station_sort:
				print(channel.code.ljust(6)+ str( channel.start_date )[:-8].ljust(21)+str( channel.end_date )[:-8].ljust(21)+"{:>6.1f}".format( channel.sample_rate ) +"{:>7.1f}".format( channel.azimuth ) + "{:>7.1f}".format( channel.dip ))
            
	print(75*'-' + ' DONE')
	return inventory
def sort_catalog(cat, chronologic=True):

	"""
	Sort ObsPy event catalog by origin times/dates 
	(True=recent events first), as e.g. two merged catalogs will not be date-time
	sorted anymore.
	
	:return: obspy catalogue object or False
	"""

	try:
			
		k      = [[event, event.origins[0].time] for event in cat]
		k.sort(key=lambda x: x[1], reverse=chronologic)
		events = [event[0] for event in k]
		cat    = Catalog(events=events)

	except Exception:
		cat = False

	return cat
def merge_cat(*cats, output='/Users/john/Desktop/catalog.ml', chronologic=True):

	"""
	Merge earthquake catalogs and sort them. Identic events will be trown out.
	"""

	events = []

	try:
		for cat in cats:
			cat = readEvents( cat )
			for event in cat:
				if event not in events:
					events.append( event )

		catalog = Catalog(events=events)
		catalog.write( output, format="QUAKEML")

		sort_catalog(catalog, chronologic=chronologic)

		print('Done !\n%s' % output )
	except Exception as err:
		message = 'ERROR:      %s' % err 
		sys.exit( message )
def produce_inventory_file():

	"""
	Produce stationXML inventory file using obspy (https://docs.obspy.org/tutorial/code_snippets/stationxml_file_from_scratch.html).
	The stationXML file is produced based on a 'stationfile' that containts all the important information, such as:
	NAME, LAT, LON, ELEVATION, DEPTH, STARTTIME, ENDTIME, REMARKS (no whitespaces), ...

	Fill in the variables at the top of the script and it should be fine.
	Note for channels, the naming, azimuths, dips is hardcoded (not position, starttime, endtime). 
	I assumed there are 4 channels: hydrophone, 2 horizontals and 1 vertical.

	Network time is assumed to be the starttime of first station and endtime of last station running.

	No response is attached to channels.

	Skew values are put into a comment for each station, along with a second comment stating remarks.

	By now, works only for one network. If you like this script to work for multiple networks, adapt it!
	"""


	# VARIABLES TO FILL
	author              = 'Alfred-Wegener-Institut: Alfons Eckstaller, Tanja Fromm & John-Robert Scholz'	# institution and operator whoever created the file.
	network_code        = 'AW'																				# network code according to the SEED standard
	network_description = 'Neumayer II seismic array from 1999 to 2009.'
	stationfile         = '/gsys/seismology1/MOVE_SWARM/JOHNS_WORK/NEUMAYER/needed_files/stations.txt'
	outfile             = '/gsys/seismology1/MOVE_SWARM/JOHNS_WORK/NEUMAYER/needed_files/stations.xml'


	# Read in data from file that was assembled by hand.
	stat_file   = np.loadtxt( stationfile, dtype=str )
	stats       = stat_file[:,0]
	lons        = stat_file[:,1].astype(np.float)
	lats        = stat_file[:,2].astype(np.float)
	eles        = stat_file[:,3].astype(np.float)
	starttimes  = stat_file[:,4]													# given in format: 2016-09-03T04:44:39 
	endtimes    = stat_file[:,5]													# given in format: 2016-09-03T04:44:39 
	samplerates = stat_file[:,6].astype(np.float)
	channels    = stat_file[:,8]													# three letters: band code, instrument code, orientation code (SEED standard)
	skews       = stat_file[:,9]													# in seconds
	remarks     = stat_file[:,10]													# no whitespaces here please



	# We'll first create all the various objects. These strongly follow the
	# hierarchy of StationXML files.
	inv = Inventory(
	    # We'll add networks later.
	    networks = [],
	    source   = author)
	
	net = Network(
	    stations                 = [],
	    total_number_of_stations = len(stats),
	    code                     = network_code,
	    description              = network_description,
	    end_date                 = UTCDateTime( sorted(starttimes)[0] ),
	    start_date               = UTCDateTime( sorted(endtimes)[-1] ))
	

	for i in range( len( stat_file )):

		bandcode          = channels[i][0]									# need because we had different band codes in our network
		comments          = [Comment("Skew is = %s s" % skews[i]),
							 Comment(remarks[i])]

		# not all stations have all components, this is one way catching this
		flag_Z = True
		flag_N = True
		flag_E = True
		flag_H = True


		# I used the remarks column in the station.txt file to state if the accordant station had one (Z) component or three (Z, E, N) components ...
		if remarks[i] == '1':
			flag_N = False
			flag_E = False
			flag_H = False
			total_number_of_channels=1
		elif remarks[i] == '3':
			flag_H = False
			total_number_of_channels=3

		sta = Station(
		    # This is the station code according to the SEED standard.
		    channels                 = [],
		    code                     = stats[i],
		    longitude                = lons[i],
		    latitude                 = lats[i],
		    elevation                = eles[i],
		    creation_date            = UTCDateTime(starttimes[i]),
		    start_date               = UTCDateTime(starttimes[i]),
		    end_date                 = UTCDateTime(endtimes[i]),
		    comments                 = comments,
		    total_number_of_channels = total_number_of_channels,
		    site                     = Site(name=stats[i]))					# description country, city, site location ...

		# change channel azimuths if you know them ...
		if flag_Z:
			chaZ = Channel(
			    # This is the channel code according to the SEED standard.
			    code          = "%sHZ" % bandcode,
			    location_code = "",
			    latitude      = lons[i],
			    longitude     = lats[i],
			    elevation     = eles[i],
			    depth         = 0.0,
			    azimuth       = 0.0,
			    dip           = -90.0,
			    start_date    = UTCDateTime(starttimes[i]),
			    end_date      = UTCDateTime(endtimes[i]),
			    sample_rate   = samplerates[i])
			sta.channels.append(chaZ)

		if flag_N:
			cha1 = Channel(
			    code          = "%sHN" % bandcode,
			    location_code = "",
			    latitude      = lons[i],
			    longitude     = lats[i],
			    elevation     = eles[i],
			    depth         = 0.0,
			    azimuth       = 0.0,
			    dip           = 0.0,
			    start_date    = UTCDateTime(starttimes[i]),
			    end_date      = UTCDateTime(endtimes[i]),
			    sample_rate   = samplerates[i])
			sta.channels.append(cha1)

		if flag_E:
			cha2 = Channel(
			    code          = "%sHE" % bandcode,
			    location_code = "",
			    latitude      = lons[i],
			    longitude     = lats[i],
			    elevation     = eles[i],
			    depth         = 0.0,
			    azimuth       = 0.0,
			    dip           = 0.0,
			    start_date    = UTCDateTime(starttimes[i]),
				end_date      = UTCDateTime(endtimes[i]),
			    sample_rate   = samplerates[i])
			sta.channels.append(cha2)

		if flag_H:
			chaH = Channel(
			    code          = "%sDH" % bandcode,
			    location_code = "",
			    latitude      = lons[i],
			    longitude     = lats[i],
			    elevation     = eles[i],
			    depth         = 0.0,
			    azimuth       = 0.0,
			    dip           = 0.0,
			    start_date    = UTCDateTime(starttimes[i]),
				end_date      = UTCDateTime(endtimes[i]),
			    sample_rate   = samplerates[i])
			sta.channels.append(chaH)

		net.stations.append(sta)

		## By default this accesses the NRL online. Offline copies of the NRL can
		## also be used instead
		#nrl = NRL()
		## The contents of the NRL can be explored interactively in a Python prompt,
		## see API documentation of NRL submodule:
		## http://docs.obspy.org/packages/obspy.clients.nrl.html
		## Here we assume that the end point of data logger and sensor are already
		## known:
		#response = nrl.get_response( # doctest: +SKIP
		#    sensor_keys=['Streckeisen', 'STS-1', '360 seconds'],
		#    datalogger_keys=['REF TEK', 'RT 130 & 130-SMA', '1', '200'])
		# Now tie it all together.
		#cha.response = response


	inv.networks.append(net)
	inv.write(outfile, validate=True)
	print('Done! Check:')
	print(outfile)
def get_station_metas(service='RESIF', network='YV', station='*', level='channel'):

	""" 
	Using ObsPy, retrieve station meta data for specified parameters.
	Returns a an ObsPy inventory file.
	"""

	### PRINTING
	print('Exec:\t\t' + sys._getframe().f_code.co_name + '()')
	print('..')
	print()
	print()


	### RETRIEVE actual data (inventory object), save inventory (outcomment second line)
	inventory = Client( service ).get_stations(network=network, station=station, level=level)
	#inventory.write('/home/john/Dropbox/JOHNS_WORK/programming/python/obspyDMT/obs_orientation/RR*_INSU_inventory.xml', 'STATIONXML')


	### OUTPUT
	counter = 0
	for network in inventory:
		print('\033[91m%s\033[0m === %s' % (network.code.upper(), network.description))

		for station in network:
			counter += 1
			print()
			print(str( counter) + '. ' + str( station.code.upper() + ':' ).ljust(6) + \
                  str( station.start_date )[:-8].ljust(21) + \
                  str( station.end_date )[:-8].ljust(21)+   \
                  "{:.4f}".format( station.longitude ).ljust(11) + \
                  "{:.4f}".format( station.latitude ).zfill(8).ljust(11) + \
                  "{:<7}".format( station.elevation ))
			print('cha   start                end                   sr       azi  dip')
			print(75 * '-')

			# sort the channels after channel name and then after time,
			# starting with newest
			k = [[channel, channel.code, channel.start_date] for channel in station]
			k.sort(key=lambda x: x[1]+str(x[2]), reverse=True)
			channels = [channel[0] for channel in k]
			station_sort = Station(station.code, station.latitude, \
                                   station.longitude, station.elevation, \
                                   channels=channels)
			for channel in station_sort:
				print(channel.code.ljust(6)+ str( channel.start_date )[:-8].ljust(21)+str( channel.end_date )[:-8].ljust(21)+"{:>6.1f}".format( channel.sample_rate ) +"{:>7.1f}".format( channel.azimuth ) + "{:>7.1f}".format( channel.dip ))
            
	print(75*'-' + ' DONE')
	return inventory
def seis2wav(abs_file, starttime=None, endtime=None, freqmin=0.01, freqmax=35, normalise_trace=False, speed=20):

	"""
	Create a wav-file from a seismogram.

	From the data, the mean and linear trend is removed, in addtion,
	a 'Hann'-taper is applied to both data ends (1%).

	Data are then Butterworth bandpass filtered (corners=2, zerophase=True)
	with frequencies given. Set both corner frequencies to 'None' to avoid
	filtering.

	The wav-file will be shorter by the factor
	of 'speed', meaning, the frequencies are raised by this factor.

	The output wav-file is stored in the same directory as the input file,
	with the same file-name but the file-ending '.wav'.

	If the trace is not normalised before, the output wav-file will be louder,
	sometimes leading to better results in terms of volume, but perhaps to 
	worse sound quality (?).


	:type abs_file:  string
	:param abs_file: absolute file path of seismological file to be converted

	:type starttime, endtime:  UTC based DateTime object
	:param starttime: define starttime of trace to be processed
	:param endtime: define endtime of trace to be processed

	:type freqmin:  float or int
	:param freqmin: lower corner frequency of Butterworth bandpass filter

	:type freqmax:  float or int
	:param freqmin: higher corner frequency of Butterworth bandpass filter

	:type normalise_trace:  Boolean
	:param normalise_trace: if True, seismological data array are normalised 
	                         by maximum absolute array value

	:type speed:  float or int
	:param speed: factor wav-file is speed up, thereby raising
	              contained frequencies by this factor
	"""

	t1 = time.time()



	### PRE-PROCESS STREAM
	tr = read(abs_file)[0]
	tr.detrend('demean')
	tr.detrend('linear')
	tr.taper(type='hann', max_percentage=0.02, side='both')
	sample_rate_original = tr.stats.sampling_rate

	if freqmin and freqmax:
		tr.filter('bandpass', freqmin=freqmin, freqmax=freqmax, corners=2, zerophase=True)

	if normalise_trace:
		tr.data   = tr.data/np.max( np.abs(tr.data) )


	### VARIABLES
	tr_time_in_s  = tr.stats.endtime - tr.stats.starttime
	wav_sps       = speed * sample_rate_original					# samples per s
	wav_file_name = '.'.join( abs_file.split('.')[:-1] ) + '_x%s.wav' % speed


	### WRITE WAV-FILE
	sci_wav.write(wav_file_name, int(wav_sps), tr.data)


	### OUTPUT SHELL
	print('')
	print('-' * 50)
	print('|  Output wav-file:')
	print('|    %s' % wav_file_name)
	print('|  Length wav-file:   %s s' % round(tr_time_in_s/speed, 1) )
	print('|')
	print('|  Length seis-file:  %s s' % tr_time_in_s)
	print('|  Filter applied:    %s Hz - %s Hz' % (freqmin,freqmax) )
	print('|  Trace normaised:   %s' % normalise_trace)
	print('|')
	print('|  Orignal frequencies raised:')
	print('|  x%s' % speed )
	print('-' * 50)
	print('')
	print('Finished in %s (h:m:s).' % sec2hms( time.time()-t1 ))
def seis2mp4(abs_file, starttime=None, endtime=None, freqmin=10, freqmax=35, normalise_traces=True, speed=20, frames=500, wlen=50, make_mp4=True):

	"""
	Create an matplotlib animation using a seismological-file.

	Three subplots are created.

	1. Unfiltered trace
	2. Filtered trace with corner frequencies given 
	   (Butterworth bandpass, corners=2, zerophase=True)
	3. Spectrogram of filtered trace 
	   Spectrogram parameters can be adjusted within the scipt.
	   Default: wlen=None, per_lap=0.90, mult=8, dbscale=True, mode='default'

	For both the unfiltered and filtered trace, 
	the mean and linear trend is removed, in addtion,
	a 'Hann'-taper is applied to both data ends (1%).

	The 'filtered' data trace is furthermore Butterworth bandpass
	filtered (corners=2, zerophase=True) with frequencies given. 
	Set both corner frequencies to 'None' to avoid fitlering of 
	this trace.

	The animation is a red line marking the actual trace time and
	running chronologically through	all 3 subplots in real time,
	meaning, one animation second refers to one second watching
	it run.

	If factor 'speed' != 1, the red line will run faster by this
	factor with respect to real time.
	Eg: speed=100, the animated red line will cover 100s whilst 
	    just one second passed watching it.

	(Attention: in dependence of frames and speed number given, 
				the animation will perhaps not play accordingly
				fast, however, if an mp4-file shall be produced, 
				it will be correct.)

	The variable 'frames' states, how many single frames shall be
	produced for the animation. Choose as much as needed but as
	least as wanted, as a high number is computational expensive.

	if make_mp4=True:
	  The animation is stored as mp4-file in the same directory as the
	  input file, with the same file-name but the file-ending '.mp4'.

	  You can adjust the mp4-parameter within the script.
	  Default: dpi=500, bitrate=5000, codec=ffmeg
	  (Those are high quality creteria, leading to big file-sizes.)

	  You can also pass some mp4-metadata information there.


	:type abs_file:  string
	:param abs_file: absolute file path of seismological file to be converted

	:type starttime, endtime:  UTC based DateTime object
	:param starttime: define starttime of trace to be processed
	:param endtime: define endtime of trace to be processed

	:type freqmin:  float or int
	:param freqmin: lower corner frequency of Butterworth bandpass filter

	:type freqmax:  float or int
	:param freqmin: higher corner frequency of Butterworth bandpass filter

	:type normalise_traces:  Boolean
	:param normalise_traces: if True, seismological data arrays are normalised 
	                         by maximum absolute array value

	:type speed:  float or int
	:param speed: factor wav-file is speed up, thereby raising
	              contained frequencies by this factor

	:type frames:  int
	:param frames: number of frames to be produced for the animation

	:type make_mp4:  Boolean
	:param make_mp4: if True, the animation will be stored as mp4-file
	"""

	t1 = time.time()

	### SEISMOLOGICAL TRACE
	try:
		st        = read(abs_file, format=None, headeronly=False)
		tr        = st[0]
		tr.trim(starttime=starttime, endtime=endtime)
		tr_or     = st[0]
	except IOError:
		print('')
		print('  File couldn\'t be read. Is it a')
		print('  seismological file? Exit ..')
		print('')
		sys.exit()

	# pre-processing
	tr_or.detrend('demean')
	tr_or.detrend('linear')
	tr_or.taper(type='hann', max_percentage=0.01, side='both')

	tr_filter = tr_or.copy()
	if freqmin and freqmax:
		tr_filter.filter('bandpass', freqmin=freqmin, freqmax=freqmax, corners=2, zerophase=True)

	if normalise_traces:
		tr_or.data        = tr_or.data/np.max( np.abs(tr_or.data) )
		tr_filter.data    = tr_filter.data/np.max( np.abs(tr_filter.data) )


	### VARIABLES
	# animation parameters
	station       = tr_or.stats.station
	tr_time_in_s  = tr_or.stats.endtime - tr_or.stats.starttime				# length trace in s
	spf           = tr_time_in_s / frames									# seconds per frame 
	interval      = (spf / speed) * 1000									# in ms
	fps           = (speed / spf)											# used for actual *mp4 file
	#print( tr_time_in_s, spf, interval, fps, frames)		
	
	# variables to adjust if wanted
	face_bg_color = (153/255., 173/255., 217/255.)							# back- and facecolour of plot
	title         = 'Records on OBS station %s' % station
	mp4_file_name = '.'.join( abs_file.split('.')[:-1] ) + '_x%s.mp4' % speed
	line_color    = 'red'

	# animation quality & metadata parameters to adjust if wanted	
	bitrate       = 5000													# if make_mp4=True
	dpi_mp4		  = 300														# if make_mp4=True
	codec         = 'ffmpeg'												# if make_mp4=True
	extra_args    = None 													# if make_mp4=True, eg: ['-vcodec', 'libx264']
	metadata = dict(title=title, artist='JRS', comment="Records on OBS")	# if make_mp4=True


	### FIGURE OBJECT & make background plots (not animated part)
	fig, axes = plt.subplots(nrows=2, ncols=1, sharex=True, dpi=None, facecolor=face_bg_color)
	fig.subplots_adjust(hspace=0.3)
	fig.canvas.set_window_title( title )

	x = np.arange( len( tr_or.data )) / tr_or.stats.sampling_rate			# in s
	lablex= -0.09															# adjustment label y-axis
	bbox=dict(boxstyle='round,pad=0.1', fc=face_bg_color, ec="black", lw=0)	# bounding box text labels

	minorLocator = MultipleLocator(500)										# minor ticks subplots
	cbarLocator  = MultipleLocator(50)										# ticks colorbar


	#SUBPLOT 1
	axes[0].plot(x, tr_or.data, lw=1)
	#axes[0].set_axis_bgcolor(face_bg_color)
	#axes[0].set_xlabel('Trace Time / [s]')
	axes[0].set_ylabel('Counts norm.')
	axes[0].ticklabel_format(style='sci', axis='y', scilimits=(0,3))
	axes[0].yaxis.set_label_coords(lablex, 0.5)
	axes[0].text(0.02, 0.85, 'Unfiltered Data', horizontalalignment='left', color='red', transform=axes[0].transAxes, bbox=bbox)
	#axes[0].set_title('Unfiltered Displacement Data')
	#axes[0].set_xlim([])
	#axes[0].set_ylim([])


	#SUBPLOT 2
	axes[1].plot(x, tr_filter.data, lw=1)
	#axes[1].set_axis_bgcolor(face_bg_color)
	#axes[1].set_xlabel('Trace Time / [s]')
	axes[1].set_ylabel('Counts norm.')
	axes[1].ticklabel_format(style='sci', axis='y', scilimits=(0,3))
	axes[1].yaxis.set_label_coords(lablex, 0.5)
	axes[1].text(0.02, 0.85, 'Filtered: %s Hz - %s Hz' % (freqmin,freqmax), horizontalalignment='left', color='red', transform = axes[1].transAxes, bbox=bbox)
	#axes[1].set_title('Filtered: %s Hz - %s Hz' % (freqmin,freqmax) )
	#axes[1].set_xlim()
	#axes[1].set_ylim([])


	#SUBPLOT 3
	#axes[2], img = spectro(tr_filter.data, tr_filter.stats.sampling_rate, wlen=wlen, per_lap=0.95, mult=8, show=False, axes=axes[2], ylim=(freqmin,freqmax), dbscale=True, mode='default' )
	##axes[2].set_axis_bgcolor(face_bg_color)
	#axes[2].set_xlabel('Trace Time (s)')
	#axes[2].set_ylabel('Frequency (Hz)')
	##axes[2].ticklabel_format(style='sci', axis='y', scilimits=(0,3))
	#axes[2].yaxis.set_label_coords(lablex, 0.5)
	#axes[2].xaxis.set_minor_locator(minorLocator)
	#axes[2].text(0.02, 0.85, 'Spectrogram filt. (db)', horizontalalignment='left', color='black', transform=axes[2].transAxes, bbox=None)
	#axes[2].set_title('Filtered: %s Hz - %s Hz' % (freqmin,freqmax) )
	#axes[2].set_xlim()
	#axes[2].set_ylim([])

	#COLOURBAR of subplot 3
	#fig.subplots_adjust(right=0.89)
	#cbar_ax = fig.add_axes([0.95, 0.1, 0.02, 0.22])
	#cbar_ax.tick_params(labelsize=8) 
	#cbar = fig.colorbar(img, cax=cbar_ax, ticklocation='left', spacing='uniform', ticks=cbarLocator)


	### ANIMATION
	anim_objects = []
	for j in range( len( axes )):
		anim_object, = axes[j].plot([], [], lw=1, color=line_color )
		ylim = axes[j].get_ylim()
		anim_objects.append( [anim_object,ylim] )


	# initialization function
	def init(anim_objects=anim_objects):

		container = ()
		for anim_object, ylim in anim_objects:
			anim_object.set_data( [],[] )
			anim_objects = container + ( anim_object, )
		
		return container


	# animation function. This is called sequentially
	def animate(i, anim_objects=anim_objects, spf=spf):

		container = ()
		for anim_object, ylim in anim_objects:
			anim_object.set_data( i*spf, ylim )
			anim_objects = container + ( anim_object, )
		
		#print(i, i*spf )
		return container


	### ANIMATOR blit=True means only re-draw the parts that have changed.
	anim = animation.FuncAnimation(fig, animate, init_func=init, frames=frames, interval=interval, blit=False)
	# save the animation as an mp4. This requires ffmpeg or mencoder to be
	# installed. The extra_args ensure that the x264 codec is used, so that
	# the video can be embedded in html5. You may need to adjust this for.
	# --> extra_args=['-vcodec', 'libx264']


	### SHOW / SAVE
	if make_mp4:
		Writer = animation.writers[codec]
		writer=Writer(fps=fps, bitrate=bitrate, metadata=metadata)
		savefig_kwargs = dict(facecolor=face_bg_color, edgecolor='none')			# dict of savefig commands

		anim.save(mp4_file_name, writer=writer, dpi=dpi_mp4, savefig_kwargs=savefig_kwargs, extra_args=extra_args)
	else:
		plt.show()


	### OUTPUT SHELL
	print('')
	print('-' * 50)
	print('|  If "make_mp4=True":')
	print('|    Output mp4-file:')
	print('|      %s' % mp4_file_name)
	print('|    Codec used:       %s' % codec)
	print('|    Bitrate used:     %s' % bitrate)
	print('|    Dots per Inch:    %s' % dpi_mp4)
	print('|    Length mp4-file:  %s s' % round(tr_time_in_s/speed, 1) )
	print('|')
	print('|  Length seis-file:   %s s' % tr_time_in_s)
	print('|  Filter applied:     %s Hz - %s Hz' % (freqmin,freqmax) )
	print('|  Traces normaised:   %s' % normalise_traces)
	print('|')
	print('|  Orignal frequencies raised:')
	print('|  x%s' % speed )
	print('|')
	print('|  Frames produced:    %s' % frames)
	print('-' * 50)
	print('')
	print('Finished in %s (h:m:s).' % sec2hms( time.time()-t1 ))
def interpolate_model(model_in, outfile, delimiter=',', fmt='%.8f', max_depth=7000, depth_column=1, depth_frac_to_km=1000):

	"""
	Allows to interpolate row-wise values of a velocity model.
	The idea is that often it's useful to have the values of such models stated for each
	x kilometers in depth (instead of whatever depth the models state), as this 
	simplifies usage for other scripts that may just look up entries for a certain depth. 
	
	dilimiter : used for model_in & output as delimiter
	max_depth : maximum depth you care for to be processed
	fmt : typical numpy format to output stuff
	depth_column : column of depth in model_in, starting to count from 1
	depth_frac_to_km : if model_in depths are in km, set to 1
					   if model_in depths are in m, set to 1000

	NOTE:
		- if depth_column is not 1, outfile will have different column order
		  than model_in, as the depth_column is always the first in outfile (for now)
		- By now, works only if model_in contains only numbers (i.e., no letters).
		  If strings shall work too, you must take care on conversions
	"""

	model  = np.loadtxt(model_in, delimiter=delimiter)
	model[:,depth_column-1] = model[:,depth_column-1]/depth_frac_to_km		# depth to km
	output = []

	### loop over lines of model file
	for i in range(len( model )):

		# keep current line in any case
		line = model[i,:]
		output.append( list(line) )

		if model[i,depth_column-1] >= max_depth:
			break
		 
		try: 
			if model[i+1,depth_column-1] - model[i,depth_column-1] > 1:
				
				x_interp = [x for x in range( int(model[i,depth_column-1])+1, int(model[i+1,depth_column-1]) ) ]

				if model[i,depth_column] == model[i+1,depth_column]:
					line = model[i,:]
					for x in x_interp:
						line[depth_column-1] = x
						output.append( list(line) )

				else:
					lines = [x_interp]
					for j in range(1, len(model[depth_column-1,:])):
						xp = [ model[i,depth_column-1],model[i+1,depth_column-1] ]
						fp = [ model[i,j],model[i+1,j] ]
						values = np.interp(x_interp, xp, fp)
						lines.append(values)
					lines = np.transpose( np.array(lines) )
					for line in lines:
						output.append( list(line) )

		except:		# catches case @ EOF
			pass 


	### SAVE TO FILE or PRINT TO SHELL
	if outfile:
		np.savetxt(outfile, output, fmt=fmt, delimiter=delimiter)
	else:
		for line in output:
			print(','.join( [str(e) for e in line] ))

	print('Finished. :)')
def read_obspyDMT_events(this, pattern='*', state='raw'):

	"""
	Read given obspyDMT path `this`, matching `pattern`.
	`state` is the kind of data you like, e.g. 'processed'
	or 'raw'.

	Return tuple of shape:

			( (event1, stream_of_all_event1_files), 
			  (event2, stream_of_all_event1_files), 
			  .. )

	where `event` is a dictionary with the event parameter,
	and stream object is ObsPy stream object.
	"""

	### read obspyDMT pickle file
	try:
		event_list_pickle = get_files( os.path.dirname(this), pattern='event_list_pickle')[0]
		events_dict       = pickle.load( open( event_list_pickle, "rb" ) )
		if not events_dict:
			raise Exception
	except:
		print("ERROR:  Could not find 'event_list_pickle' file. Or it is empty.")
		return

	### all files that match `pattern` in given directory `this` 
	all_seismic_files  = [f for f in get_files( this, pattern=pattern ) if re.search(state, os.path.dirname(f))]
	for event_dict in events_dict:

		event_seismic_files = [f for f in all_seismic_files if re.search(event_dict['event_id'], os.path.dirname(f))]
		if not event_seismic_files:
			continue

		stream = Stream()
		for event_seismic_file in event_seismic_files:
			stream += read(event_seismic_file)

		yield event_dict, stream
def view_events(this, pattern='*BH?', archive='obspyDMT', filter={'freqmin':0.02, 'freqmax':0.2, 'corners':2, 'zerophase':True}, rotate_NE=False, rotate_LQT=False, final_process=None, plot_station_wise=True, v_lines=True, **kwargs):

	"""
	  - kwargs passed to ObsPy's plot function
	"""

	### further parameters needed
	try:
		# stat : [cw angle from N to BH1, lon, lat, elevation in m]
		stats_parameters = {'RR29' : [266, 51.7488, -24.9657, -4829],
	                        'MAYO' : [0,   045.1868, -12.8456, 0]}
		for stat in stats_parameters.keys():
			if len(stats_parameters[stat]) != 4:
				raise Exception
	except:
		print('Please provide 4 station parameters:')
		print('orientation cw N-->BH1 or BHN, lon, lat, elevation in m!')
		return                	   

	try:	      
		# model : [phase1, phase2, ...]
		model_phases     = { 'prem' : ['Pdiff', 'P', 'S', 'SKS'],
	    	                'ak135' : ['Pdiff', 'P', 'S', 'SKS'],
	        	           'iasp91' : ['Pdiff', 'P', 'S', 'SKS'],}
		if not model_phases:
			raise Exception
	except:
		print('Please provide (global) v-models and respective phases to be considered,')
		print('even if you do not want to rotate into the L-Q-T frame.')
		print('This is important for the inclination calculations.')
		print('(if rotate_LQT=True, used inclination is average of all inclinations of all v-models & phases.)')
		return        	          


	### retrieve event_stream tuple of structure: ((event1_dict, stream_all_files_event1), (..), ..)
	if archive == 'obspyDMT':
		function = read_obspyDMT_events
	else:
		print("Method '%s' for archive reading not yet implemented. Please check." % archive)
		return


	### loop over events and accordant streams
	for data_event, stream_event in function(this, pattern=pattern):
		stream_event       = Stream2( stream_event )

		# retrieve event information
		source_latitude    = data_event['latitude']
		source_longitude   = data_event['longitude']
		source_depth_in_km = data_event['depth']
		source_time        = data_event['datetime']
		source_magnitude   = data_event['magnitude']
		source_mag_type    = data_event['magnitude_type']
		source_id          = data_event['event_id']
		if (source_id!='20130405_130002.a'):# and source_id!='20130904_001824.a'):
			continue

		# print event information
		print()
		print('---------------------------------------')
		print('EVENT: %s' % source_time)
		print('   %11s : %-10.4f'  % ('lon', source_longitude))
		print('   %11s : %-10.4f'  % ('lat', source_latitude))
		print('   %11s : %-10.1f'  % ('dep', source_depth_in_km))
		print('   %11s : %-.1f %s' % ('mag', source_magnitude, source_mag_type, ))
		print('   %11s : %-10s'    % ('id',  source_id))
		
		# loop over stations within each event
		stream_final = Stream2()
		stats        = sorted( list( set( [tr.stats.station for tr in stream_event] )))
		verticals    = []
		for stat in stats:

			#retrieve station & phase information
			station_orientation  = stats_parameters[stat][0]
			station_longitude    = stats_parameters[stat][1]
			station_latitude     = stats_parameters[stat][2]
			receiver_depth_in_km = stats_parameters[stat][3]/-1000.
			dist_in_m, baz, az   = gps2dist_azimuth(station_latitude, station_longitude, source_latitude, source_longitude, a=6378137.0, f=0.0033528106647474805)
			distance_in_degree   = locations2degrees(station_latitude, station_longitude, source_latitude, source_longitude)
			inclinations         = []
			all_arrivals         = []
			for model in model_phases.keys():
				velo_model = tau.TauPyModel(model=model.lower(), verbose=False, planet_flattening=0.0, cache=None)
				phase_list = model_phases[model]
				arrivals   = velo_model.get_travel_times(source_depth_in_km, distance_in_degree=distance_in_degree, phase_list=phase_list, receiver_depth_in_km=receiver_depth_in_km)
				all_arrivals.append( arrivals )
				for arrival in arrivals:
					inclinations.append(arrival.incident_angle)
			inclination = sum( inclinations )/len(inclinations) 		# average over all inclinations of all models and all phases

			# print station & seismic information
			print()		
			print('STATION: %s' % stat)
			print('   %11s : %-10.4f' % ('lon', station_longitude))
			print('   %11s : %-10.4f' % ('lat', station_latitude))
			print('   %11s : %-10.1f' % ('dep', receiver_depth_in_km))
			print('   %11s : %-10.1f' % ('orient', station_orientation))
			print('SEISMIC PHASE:')
			print('   %11s : %-10.1f' % ('baz', baz))
			print('   %11s : %-10.1f' % ('az', az))
			print('   %11s : %-10.1f' % ('dist deg', distance_in_degree))
			print('   %11s : %-10.1f' % ('dist km', dist_in_m/1000.))
			print('   %11s : %-10.1f' % ('incl.', inclination))
			for j, model in enumerate(model_phases.keys()):
				for arrival in all_arrivals[j]:
					print('   %11s : tt=%-.1f, incl=%-10.1f' % ((model+':'+arrival.name), arrival.time, arrival.incident_angle))

			#select stream as to station & filter stream
			stream_select = stream_event.select(station=stat)
			stream_select.filtering(filter=filter)
			
			if rotate_NE or rotate_LQT:
				# get components needed
				stream_select.trim_common()
				comp_1 = stream_select.select(component='1') or stream_select.select(component='N')
				comp_2 = stream_select.select(component='2') or stream_select.select(component='E')
				comp_Z = stream_select.select(component='Z')

				# check condition
				condition = len(comp_Z)==1 and len(comp_1)==1 and len(comp_2)==1
				if not condition:
					print('ERROR:  Either one or more component not found or data gaps for this event & station.')
					print('        Event: %s, stat: %s left unrotated.' % (event_id, stat))
					continue
	
      		    # rotate using self written function (note that `stream` data are changed as select creates only reference)
				comp_1[0].data, comp_2[0].data = rotate_2D(comp_1[0].data, comp_2[0].data, station_orientation, clockwise=False)
				comp_1[0].stats.channel        = comp_1[0].stats.channel[:-1] + 'N'
				comp_2[0].stats.channel        = comp_2[0].stats.channel[:-1] + 'E'

				# no need to tell obspy order of components, it chooses automatically based on channel naming
				if rotate_LQT:
					stream_select.rotate('ZNE->LQT', baz, inclination)

				# do further processing of stream, if wished (one can pass a function containing the steps, but it must return a stream object)
				if final_process:
					stream_select = final_process(stream_select)

			# verticals (pass times as relative seconds with respect to beginning of respective axis)
			if v_lines:
				stream_select.trim_common()
				for i, tr in enumerate(stream_select):
					verticals.append( [source_time-tr.stats.starttime, 'origin', 'teal', i] )
					for arrival in arrivals:
						vertical = [source_time+arrival.time-tr.stats.starttime, '%s arrival' % arrival.name, 'yellowgreen', i]
						verticals.append( vertical )
						
			if plot_station_wise:
				stream_select.plot(store_name='Event ID: %s' % source_id, verticals=verticals, save_and_no_show=False, **kwargs)
			else:
				stream_final += stream_select

		if not plot_station_wise:
			stream_final.plot( store_name='Event ID: %s' % source_id, verticals=verticals, save_and_no_show=False, **kwargs)
		print('---------------------------------------')
def pierce_stream(stream): 

	"""
	Does NOT work if there are any overlaps in the original stream
	for one particular channel. That is why ObsPy's merge function
	is first applied to stream.

	See options:
	https://docs.obspy.org/packages/autogen/obspy.core.trace.Trace.__add__.html#obspy.core.trace.Trace.__add__

	After that:
	Pierce stream into X streams, in which each one contains all unique channels of the input stream,
	all with same start and end times. For example useful, if you need 3-D component, equally trimmed data
	that are contained within one gappy stream ..

	NOTE: all unique traces in stream are considered. so if there are 5 different traces, the
	resulting streams will have the start and end times where all 5 traces are present, nothing more.
	So make sure stream contains only those traces that you'd like.
	"""


	### handle overlaps first (if there any gaps, data will be masked array), that is why then split masked array contiguous unmasked traces to apply piercing
	stream_copy = stream.copy()
	stream_copy.merge(method=1, fill_value=None, interpolation_samples=0)
	stream_copy = stream_copy.split()


	### piercing traces 
	ids           = list(set( [tr.id for tr in stream_copy] ))
	array         = np.array([[len(stream_copy.select(id=id)),id] for id in ids])
	id_max_traces = array[np.argmax(array[:,0].astype('float')),1]

	new_streams   = []
	for trace in stream_copy.select(id=id_max_traces):


		#traces     = [tr for tr in stream_copy if tr.id != trace.id and tr.stats.endtime>=trace.stats.starttime and tr.stats.starttime<=trace.stats.endtime]
		traces = []
		for tr in stream_copy:
			if tr.id != trace.id:
				if tr.stats.endtime>=trace.stats.starttime and tr.stats.starttime<=trace.stats.endtime:
					traces.append(tr)

		traces     = [trace] + traces
		tmp_stream = Stream2(traces=traces).copy()
		tmp_ids    = list(set( [tr.id for tr in tmp_stream] ))


		if len(tmp_ids)<len(ids):
			continue
		
		
		if len(tmp_stream) != len(ids):
			counter = 0
			while True:
				try:
					tmp_stream2 = tmp_stream.copy()
					for _ in range( len(tmp_stream)-len(ids) ):
							for _ in range(len(tmp_stream2)):
								index = random.randint(0,len(tmp_stream2))
								try:
									tmp_stream2.pop( index )
									break
								except IndexError:
									continue
		
					tmp_stream2.trim_common()
					new_streams.append(tmp_stream2)
				except ValueError:
					continue
		
				counter += 1
				if counter==len(tmp_stream)-len(ids)+1:
					break
		
		else:
			tmp_stream.trim_common()
			new_streams.append(tmp_stream)


	new_streams = sorted(new_streams, key=lambda x: x[0].stats.starttime)
	return new_streams

# illustrational
def quick_plot(*y, x=[], data_labels=[], lw=1.5, win_title='', title='', xlabel='x-axis', ylabel='y-axis', x_invert=False, y_invert=False, xlim=(), ylim=(), verts=(), outfile=None, show=True):

	"""
	This function allows for some convenient, quick 2-D plotting.
	
	It requires loaded classes `Trace`, `Stream`, `Stream2`,
	`UTCDateTime`, and `datetime.datetime`.
	
	Otherwise, it should be very user friendly. Times will
	be displayed as '%d/%m %H:%M:%S', which is harcoded.


	PARAMETERS:
	-----------

	`y` : - y-data to be plotted, as many as you like
	      - if your data are `Trace` object(s), it is plotted
	        with date strings (if no `x` is given).
	      - if your data are a `Stream` or `Stream2` object, it is plotted
	        with their internal plot functions and further parameters
	        are not passed.
	      - if none of above, then x-data are determined using np.arange
	`x` : - this overrides x-data (dates) when passing `Trace` objects as `y`,
	        so you can assign anything as x-data and still conveniently use
	        `Trace` objects.
	      - `UTCDateTime`, `datetime.datetime`, and `str` objects will
	         be converted to `UTCDateTime` and then matplotlib times,
	         then plotted as date strings '%d/%m %H:%M:%S'
	`data_labels` : - will appear in legend
	                - should have as many elements as you pass in `y`,
	                  otherwise it will not used but data simply numbered
	`lw`       : - linewidth of plots (lines are harcoded, no option to plot points only)
	`win_title`: - title of figures canvas (does not appear in plot)
	`title`    : - sup-title of figure instance (appears as title in plot)
	`xlabel    : - will appear beneath x-axis
	`ylabel    : - will appear next to y-axis
	`x_invert` : - will invert direction of x-axis, if True
	`y_invert` : - will invert direction of y-axis, if True
	`xlim`     : - list (or similar) with two elements. Choose same type as `y`. If 'xlim' is set, `ylim` is set automatically as to the y-range within these xlim
	`ylim`     : - list (or similar) with two elements.
	`verts`    : - either list or list of lists
	             - elements will be interpreted as x-data where vertical lines shall be plotted
	             - each list gets an own colour for all elements
	             - when elements are: - `UTCDateTime`, `datetime.datetime`, and `str` 
	                                     --> conversion to matplotlib times (so x-data should be times, too)
	                                  - `float` or `int`
	                                     --> where you would expect it to happen (each sample on x-axis numbered).
	                                     	 If x-data are times (e.g. matplotlib), you can pass here matplotlib times as well
	                                  - `complex` (e.g. 3+0j)
	                                     --> real part will be interpreted as relative position on x-axis,
	                                         so (0.5+0j) will plot a vertical line in the middle of the x-axis.
	                                         Imaginery part is neglected.
	`outfile`  : - (absolut) path where plot shall be saved. Use endings like '.pdf' or '.png'
	`show`     : - if True (default), plot will be shown, otherwise not until the next show=True command ;)
	             - if False, figure instance of matplotlib will be returned so it can be further used
	             - no matter the status of `show`, if an `outfile` is specified a plot will be saved

	NOTES:
	------
	  - Matplotlib internally uses matplotlib times, always.
	    (example: 733643.01392361, type=numpy.float64)

	    This holds true also when extracting times, for example via
	    ax.get_ylim() or similar, no matter the original input format. Therefore,
	    this script converts `str`, `UTCDateTime` and `datetime.datetime` objects all to matplotb times
	    before the plot command (although matplotlib can plot datetime.datetime objects natively), so
	    these converted times later can be used to adjust xlim and ylim (if specified).
	  - natively, matplotlib plots `datetime.datetime` objects with labels
	    hour, minute, and second. This is true for both ax.plot and ax.plot_date ..
	  - as I wished a personal format including month, day, hour, minute and second, and
	    as times are converted to matplotlib times regardless (which matplotlib then does not
	    plot with a nice string but as simples numbers), the following two commands are used:
		
		myFmt = mdates.DateFormatter('%d/%m %H:%M:%S')
		ax.xaxis.set_major_formatter(myFmt)
	"""


	### Figure instance & settings
	fig = plt.figure(figsize=(16,10), num=title)
	fig.canvas.set_window_title(win_title)
	ax = fig.add_subplot(111)
	ax.set_title(title, fontsize=13)
	ax.set_xlabel(xlabel, fontsize=12)
	ax.set_ylabel(ylabel, fontsize=12)
	ax.grid(ls='-.', lw=0.5)
	if y_invert:
		ax.invert_yaxis()
	if x_invert:
		ax.invert_xaxis()


	### Labels
	if len(data_labels)!=len(y):	# not same amount of entries (also true if none given)
		
		if all(isinstance(ele, Trace) for ele in y):
			labels = [trace.id for trace in y]

		else:
			labels = ['Data %s' % (i+1) for i in range(len(y))]

	else:
		labels = data_labels


	### Data plotting
	for j, data in enumerate(y):

		if isinstance(data, Trace):
			if list(x):
				if all(isinstance(ele, (UTCDateTime, datetime.datetime, str)) for ele in x):
					x     = [UTCDateTime(e).datetime for e in x] 	# convert all to datetime.datetime objects
					x     = mdates.date2num(x)						# convert all to matplotlib times (wouldn't need to, but for later it is need as ax.get_xlim() retrieves matplotlib times!)
					ax.plot(x, data.data, '-', lw=lw, label=labels[j])
					myFmt = mdates.DateFormatter()
					ax.xaxis.set_major_formatter(myFmt)				# because we want dates in customised way, not matplotlib times
				else:	# normal array, relative times, matlab times, or POSIX timestamps
					ax.plot(x, data, lw=lw, label=labels[j])		
			else:
				myFmt = mdates.DateFormatter('%d/%m %H:%M:%S')
				xdata = data.times(type="matplotlib")
				ax.plot(xdata, data.data, '-', lw=lw, label=labels[j])
				ax.xaxis.set_major_formatter(myFmt)

		elif isinstance(data, (Stream, Stream2)):
			print(u'Using plot function of %s object and then return.' % type(data))
			print(u'No further variables passed.')
			data.plot()		# use Stream or Stream2, respectively, object's plot function
			return

		else:
			if list(x):
				if all(isinstance(ele, (UTCDateTime, datetime.datetime, str)) for ele in x):
					xdata = [UTCDateTime(e).datetime for e in x] 	# convert all to datetime.datetime objects
					xdata = mdates.date2num(xdata)						# convert all to matplotlib times (wouldn't need to, but for later it is need as ax.get_xlim() retrieves matplotlib times!)
					myFmt = mdates.DateFormatter('%d/%m %H:%M:%S')
					ax.plot(xdata, data, '-', lw=lw, label=labels[j])
					ax.xaxis.set_major_formatter(myFmt)				# because we want dates in customised way, not matplotlib times
				else:	# normal array, relative times, matlab times, or POSIX timestamps
					ax.plot(x, data, lw=lw, label=labels[j])
			else:
				xdata = np.arange(len(data))
				ax.plot(xdata, data, lw=lw, label=labels[j])


	### Limits of x- and y-axis
	if list(xlim):
		if all(isinstance(ele, (UTCDateTime, datetime.datetime, str)) for ele in x):
			xlim = [UTCDateTime(e).datetime for e in xlim]
			xlim = [mdates.date2num(e) for e in xlim]
		ax.set_xlim( xlim )
	else:
		xlim = ax.get_xlim()
	
	if list(ylim):
		ax.set_ylim( ylim )
	else:			# make sure ylim is according to newly set xlim
		y_mins = []
		y_maxs = []
		for line in ax.lines:
			x = line.get_xdata()
			y = line.get_ydata()
			i = np.where( (x >= xlim[0]) & (x <= xlim[1]) )[0]		 # all indexes of y_data according to xlim
			line_min = y[i].min()									 # get minimum y within all data according to xlim
			line_max = y[i].max()									 # get maximum y within all data according to xlim
			y_mins.append(line_min)
			y_maxs.append(line_max)

		y_min = min(y_mins)
		y_max = max(y_maxs)
		ylim  = [y_min-0.025*np.abs(y_max-y_min), y_max+0.025*np.abs(y_max-y_min)] # give 2.5% margin that works both for pos. & neg. values 
		ax.set_ylim( ylim )		

	
	### Vertical lines for data indications
	xlim    = ax.get_xlim()
	ylim    = ax.get_ylim()
	colours = ['k', 'grey', 'lightgrey', 'red', 'green']
	if list(verts):

		try:
			if not all(isinstance(ele, (np.ndarray, tuple, list)) for ele in verts): 	#empty list `()` or `[]` is True ..
				raise TypeError
		except TypeError:	# UTCDateTime object is not iterable (case when verts= list of UTCDateTimes only)
			verts = [verts]																# make list of lists so fo-loops work

		for k, vert in enumerate(verts):
			colour_index = k % len(colours)
			for vline in vert:
					
				if isinstance(vline, (UTCDateTime, datetime.datetime, str)):
					vline = UTCDateTime( vline )
					vline = mdates.date2num(vline.datetime)
					
				elif isinstance(vline, (complex)):
					vline = vline.real
					vline = xlim[0] + (xlim[1]-xlim[0])*vline 	# relative to beginning in %

				else:
					pass

				if vline>=xlim[0] and vline<=xlim[1]:
					ax.plot( [vline, vline], ylim, lw=1, color=colours[colour_index])


	### Saving & showing ..
	ax.legend()
	plt.tight_layout()

	if outfile:
		plt.savefig(outfile, bbox_inches='tight', dpi=300)

	if show:		
		plt.show()
	else:
		return fig


################  _ _ N A M E _ _ = = " _ _ M A I N _ _ "  ################
if __name__ == "__main__":
	#argues = sys.argv
	#eval(argues[1])
	pass