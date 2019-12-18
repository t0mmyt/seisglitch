#!/usr/bin/env ipython3
# -*- coding: utf-8 -*-

#----------------------------------------------------------------------
#   Filename:  glitch_detector.py
""" Functions for seismological applications, using Python & ObsPy. """
#   Author:    John - Robert Scholz, Anna Horleston
#   Date:      July 2019 - present
#   Email:     john.robert.scholz@gmail.com
#   License:   GPLv3
#---------------------------------------------------------------------


###############  python 2.7.x & 3.x modules import  ######################
from __future__ import (absolute_import, unicode_literals, division, print_function) 

import os
import re
import sys
import copy
import glob
import time
import scipy
import fnmatch
import datetime
import math as M
import numpy as np
import shutil


###############  mtplotlib modules import  ######################
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as mticker



###############  obspy modules import  ######################
import obspy
from obspy import read_inventory, UTCDateTime
from obspy.core.inventory import Inventory, Network, Station, Channel, Site, Comment
from obspy.imaging.cm import obspy_sequential
from obspy.core.trace import Trace
from obspy.core.stream import Stream
from obspy.signal import rotate


######################  code  ######################
# mathematical
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


	comp_1_new =  comp_1*np.cos(angle*np.pi/180) + comp_2*np.sin(angle*np.pi/180)
	comp_2_new = -comp_1*np.sin(angle*np.pi/180) + comp_2*np.cos(angle*np.pi/180)

	return comp_1_new, comp_2_new
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
def normalise(data, scale_to_between=[]):

	"""
	Normalise passed data (array-like).
	scale_to_between is list with 2 elements.
	"""

	data = np.array(data)

	
	if len(data) == 0:
		return data

	if scale_to_between:
		scale_to_between.sort()
		if len(data)== 1:
			data /= data * scale_to_between[1]
			return data
		scale  = abs(scale_to_between[-1]-scale_to_between[0]) / 2.
		drange = max(data)-min(data)
		data  *= 2 / drange
		data   = data - max(data)+1				# data have y values between [-1,1]
		data  *= scale							# data have y values between [-1/scale,1/scale] 
		data  += scale_to_between[-1]-scale 	# eacht trace has y values filling range `scale_to_between`
	
	else:
		if len(data)== 1:
			data /= data
			return data
		data /= max(abs(data))

	return data
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
	phi_2D                           = (np.arctan2( eig_vec_2D_hori_1[1], eig_vec_2D_hori_1[0] ) * 180/np.pi )%360
	err_phi_2D                       = np.arctan( np.sqrt( eig_val_2D_hori[1]/eig_val_2D_hori[0] )) * 180/np.pi 	    # ...
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
		phi_3D                           = (np.arctan2( eig_vec_3D_1[1], eig_vec_3D_1[0] ) * 180/np.pi) % 360
		err_phi_3D                       = np.arctan( np.sqrt( eig_val_3D[2]/(eig_val_3D[1]+eig_val_3D[0]) )) * 180/np.pi 		# ...
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
		INCapp_3D                        = np.arctan( eig_vec_2D_radZ_1[1] / eig_vec_2D_radZ_1[0] ) * 180/np.pi
		err_INCapp_3D                    = np.arctan( np.sqrt( eig_val_2D_radZ[1]/eig_val_2D_radZ[0] )) * 180/np.pi


		### 2-D phi, radial & Z-comp plane
		comp_R, comp_T                   = rotate.rotate_ne_rt(comp_1, comp_2, phi_2D)
		covariance_matrix_2D_radZ        = covariance_matrix(comp_Z, comp_R)
		eig_val_2D_radZ, eig_vec_2D_radZ = np.linalg.eig(covariance_matrix_2D_radZ)
		index_array_descending_eig_val   = np.argsort( eig_val_2D_radZ )[::-1]
		eig_val_2D_radZ                  = eig_val_2D_radZ[index_array_descending_eig_val]							# E-values descending
		eig_vec_2D_radZ                  = eig_vec_2D_radZ[:,index_array_descending_eig_val]						# E-vectors sorted acc. to E-values
		eig_vec_2D_radZ_1                = eig_vec_2D_radZ[:,0]

		# derived
		INCapp_2D                        = np.arctan( eig_vec_2D_radZ_1[1] / eig_vec_2D_radZ_1[0] ) * 180/np.pi
		err_INCapp_2D                    = np.arctan( np.sqrt( eig_val_2D_radZ[1]/eig_val_2D_radZ[0] )) * 180/np.pi
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
	print(stream)


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

# InSight stuff
def pierce_stream(stream, method=1, fill_value=None, interpolation_samples=0): 

	"""
	Does NOT work if there are any overlaps in the original stream
	for one particular channel. That is why ObsPy's merge function
	is first applied to stream.

	See options:
	https://docs.obspy.org/packages/autogen/obspy.core.trace.Trace.__add__.html#obspy.core.trace.Trace.__add__

	After that:
	Pierce stream into X streams, in which each trace represents a unique channel of the input stream,
	all with same start and end times. For example useful, if you need 3-D component, equally trimmed data
	that are contained within one gappy stream ..
	"""


	### handle overlaps first (if there any gaps, data will be masked array), that is why then split masked array contiguous unmasked traces to apply piercing
	stream_copy = stream.copy()
	stream_copy.merge(method=method, fill_value=fill_value, interpolation_samples=interpolation_samples)
	stream_copy = stream_copy.split()


	### piercing traces 
	ids           = list(set( [tr.id for tr in stream_copy] ))
	array         = np.array([[len(stream_copy.select(id=id)),id] for id in ids])
	id_max_traces = array[np.argmax(array[:,0].astype('float')),1]

	new_streams   = []
	for trace in stream_copy.select(id=id_max_traces):


		traces     = [tr for tr in stream_copy if tr.id != trace.id and tr.stats.endtime>=trace.stats.starttime and tr.stats.starttime<=trace.stats.endtime]
		traces     = [trace] + traces
		tmp_stream = Stream2(traces=traces)
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
def rotate2VBBUVW(stream, inventory=read_inventory( max( glob.glob( '/home/scholz/Downloads/dataless*.seed'), key=os.path.getctime)), is_UVW=False, plot=False):

	"""
	Returns data (no copy of data)!
	No gain correction applied, must be done beforehand!

	https://docs.obspy.org/packages/autogen/obspy.signal.rotate.rotate2zne.html
	"""


	# VALUES EXTRACTED from INVENTORY FILE
	VBBU_azimuth = 135.10
	VBBV_azimuth = 15.00
	VBBW_azimuth = 255.00
	
	VBBU_dip     = -29.40	# defined as down from horizontal!
	VBBV_dip     = -29.20
	VBBW_dip     = -29.70


	# case data (SP or VBB) are in UVW, first go to ZNE using inventory file
	if is_UVW:
		stream_new = stream.select(component='[UVW]')
		stream_new.sort()
		stream_new.rotate('->ZNE', inventory=inventory, components=('UVW'))		# in case SP data are passed (or similar)

	else:
		stream_new = stream.select(component='[ZNE]')
		stream_new.sort(reverse=True)


	# rotate from ZNE data to VBB UVW
	stream_new[0].data, stream_new[1].data, stream_new[2].data = rotate.rotate2zne(	data_1     = stream_new[0].data,
																					azimuth_1  = VBBU_azimuth, 
																					dip_1      = VBBU_dip, 

																					data_2     = stream_new[1].data, 
																					azimuth_2  = VBBV_azimuth, 
																					dip_2      = VBBV_dip, 

																					data_3     = stream_new[2].data, 
																					azimuth_3  = VBBW_azimuth, 
																					dip_3      = VBBW_dip, 

																					inverse    = True)

	# because of sorting, can just index new component
	stream_new[0].stats.channel = stream_new[0].stats.channel[:-1] + 'U'
	stream_new[1].stats.channel = stream_new[1].stats.channel[:-1] + 'V'
	stream_new[2].stats.channel = stream_new[2].stats.channel[:-1] + 'W'


	# plot, if desired
	if plot:
		if is_UVW:
			data_labels=('U old', 'U new', 'V old', 'V new', 'W old', 'W new')
		else:
			data_labels=('Z old', 'U new', 'N old', 'V new', 'E old', 'W new')
		quick_plot(stream[0], stream_new[0], stream[1], stream_new[1], stream[2], stream_new[2], data_labels=data_labels)
	
	return stream_new
def glitch_detector(seismic_file, correct_gain=True, window=5, step_in_samples=10, input_polar=None, threshold=0.97, glitch_length=25, plot=False):

	"""

	--------------------------------------------------------
	|                                                      |
	|   Hello, are you sure about getting into glitches?   |
	|                                                      |
	--------------------------------------------------------


	This script can be run in two different ways.
	For both ways, the script expects `seismic_file` (absolute path)
	to contain at least three seismic components (UVW or ZNE).
	No matter which one is passed, two streams are created; one with ZNE &
	one with UVW. Decide if you want `correct_gain` (True or False).
	Note: the file is read as an ObsPy `Stream2` object - a self-written class.
	Note: for rotating the components, the `_get_inventory()` method
	(see Stream2 class) attempts to read an inventory file that follows
	are certain naming in a certain folder. Please adjust line 401!

	Way 1): `input_polar` = None

		Using a 3-D PCA (function `ppol_moving`), the rectinilinearity (polarisation)
		and other measures are calculated for each moving window of length `window` 
		(in seconds), stepping `step_in_samples`. These polarisations are written into 
		`output_polar`.

		BE AWARE: For many data (few days), this step may take a while!


	Way 2): `input_polar` = output of way 1)

		The calculated polarisations for each window (see way 1) are 
		looped over and glitches are extracted.
		The glitch condition is: `POLRZ` > `threshold` and
		glitch start time is later than the end time of previous glitch.

		==> Here one can perhaps find a better condition!
		(`POLRZ` is 13th column of output file of way 1)

		Once a glitch start time is detected, a window around it is extracted:
		  --> 60s before glitch start + `glitch_length` + 60 s after. This window is demeaned
		and then another 3-D PCA (same code) is run over ONLY the glitch to detect all glitch parameters
		(amplitudes of ZNEUVW, BAZ+err, INC+err, rectinilineairty, SNR ...).
		All glitch parameters are written to `output_glitch`.
		You may want to plot the result if you wish, `plot`=True.


	SUGGESTIONS:
	  - For VBB, I worked with: window=5, step_in_samples=10, threshold=0.97
	  - Since I saw (too) many polarised signals on SP, maybe window (way 1) and threshold (way 2)
	    could be increased
	  - data could be filtered prior to running this code, perhaps this can improve
	    the performance, too
	  - A more fancy and elaborated glitch condition can probably be found - why not?
	  - the code doesn't handle data gaps, so run the code only on files having no gaps
	    (see e.g. function `pierce_stream`)
	"""


	### ENTRY OF SCRIPT
	now           = time.time()
	streamVBB     = read(seismic_file)
	output_polar  = os.path.join( os.path.dirname(streamVBB.origin), 'polarization_' + '.'.join(os.path.basename(streamVBB.origin).split('.')[:-1]) + '.txt' )
	output_glitch = os.path.join( os.path.dirname(streamVBB.origin), 'glitches_'     + '.'.join(os.path.basename(streamVBB.origin).split('.')[:-1]) + '.txt' )


	### GET STREAMS of ZNE and UVW components, depending on input orientation and if gain correction is wished
	streamVBB_components = list( set( [tr.stats.channel[-1].upper() for tr in streamVBB] ))
	streamVBB.trim_common()


	if 'Z' in streamVBB_components and 'N' in streamVBB_components and 'E' in streamVBB_components:
		streamZNE = streamVBB.select(component='[ZNE]')
		streamUVW = rotate2VBBUVW(streamZNE.copy(), inventory=streamZNE.inventory, is_UVW=False)
		if correct_gain:
			streamUVW.gain_correction()
			streamZNE = streamUVW.copy()
			streamZNE.rotate('->ZNE', inventory=streamZNE.inventory, components='UVW')			

	elif 'U' in streamVBB_components and 'V' in streamVBB_components and 'W' in streamVBB_components:
		streamUVW = streamVBB.select(component='[UVW]')
		if correct_gain:
			streamUVW.gain_correction()
		streamZNE = streamUVW.copy()
		streamZNE.rotate('->ZNE', inventory=streamUVW.inventory, components='UVW')

	else:
		print(u'Did you pass a stream with either `ZNE` or `UVW` components?')
		print(u'Leave.')
		return

	streamZNE.sort(reverse=True)	# sorted: Z, N, E
	streamUVW.sort(reverse=False)	# sorted: U, V, W
	print()
	print('STREAMS:')
	print(streamZNE)
	print(streamUVW)


	### LOOP (MOVING WINDOWS)
	if not input_polar:
		print()
		print(u'DATA ANALYSIS:  %s - %s' % (streamVBB[0].stats.starttime, streamVBB[0].stats.endtime))
		print(u'Running moving windows of P-wave polarization ..')

		# MOVING P-WAVE POLARIZATION
		ppol_moving(streamZNE, window=window, step_in_samples=step_in_samples, output_polar=output_polar, print_info=True)

		# ADDING TO SHELL OUTPUT
		print(u'Done in:          %s (h:m:s).' % sec2hms( float( time.time()-now )) )
		print(u'Output polar.:    %s'          % output_polar)


	else:
		# assign needed variables
		data_start                 = str(streamVBB[0].stats.starttime)
		data_end                   = str(streamVBB[0].stats.endtime)
		data_delta                 = streamVBB[0].stats.delta
		
		polarizations              = np.loadtxt(input_polar, dtype='str')		
		polarizations_start        = UTCDateTime(polarizations[0][-2])
		polarizations_end          = UTCDateTime(polarizations[-1][-1])
		polarizations_num          = len(polarizations)
		
		pol_start_sorted           = sorted(polarizations[:,-2])
		pol_end_sorted             = sorted(polarizations[:,-1])
		pol_min_start, pol_max_end = min(pol_start_sorted), max(pol_end_sorted)

		if data_end<=pol_min_start or data_start>=pol_max_end:
			# if data do not at all overlap with polarization windows, return
			print(u'Data start/end times do not overlap with polarization start/end times. Return.')
			return
		else:
			print()
			print(u'DATA ANALYSIS:  %s - %s' % (data_start, data_end))
			print(u'Running glitch analyses on polarizations ..')

		window                     = UTCDateTime(polarizations[0][-1])-UTCDateTime(polarizations[0][-2])
		window_overlap             = UTCDateTime(polarizations[0][-1])-UTCDateTime(polarizations[1][-2]) 
		window_step                = round(window-window_overlap, 6)

		# prepare data outputs
		header                     = 'NUM     POL                START                  END      Z-AMP      N-AMP      E-AMP      U-AMP      V-AMP      W-AMP      SNR    DIR  DERR    INC  IERR'
		line_formatter             = '%05d:  %5.3f  %15s  %15s  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %7.1f  %5.1f  %4.1f  %5.1f  %4.1f'
		print('  '+header)
		
		# prepare loop variables (loop over polarization windows)
		glitch_starts, glitch_ends = [], []
		glitches                   = []
		glitch_end, glitch_counter = 0, 0

		for i in range(1, len(polarizations)):

			if not (polarizations[i-1][-2]>=data_start and polarizations[i-1][-1]<=data_end):
				# if data do not fully overlap with polarization window, continue
				continue
			
			pol_value_b = float(polarizations[i-1][12])
			pol_start_b = UTCDateTime(polarizations[i-1][-2])
			pol_end_b   = UTCDateTime(polarizations[i-1][-1])
			
			pol_value   = float(polarizations[i][12])
			pol_start   = UTCDateTime(polarizations[i][-2])
			pol_end     = UTCDateTime(polarizations[i][-1])
			pol_mid     = pol_start+(pol_end-pol_start)/2.
			
			try: 
				# GLITCH CONDITION
				#print(pol_value, threshold, pol_value_b)
				if not (pol_value_b<threshold and pol_value>=threshold and pol_start>=glitch_end):
					raise Exception
			except Exception:
				continue 


			# PREPARE DATA FOR GLITCH PARAMETERS
			glitch_start     = pol_mid
			
			streamZNE_window = streamZNE.copy()
			streamUVW_window = streamUVW.copy()
			streamZNE_window.trim(starttime=glitch_start-60, endtime=glitch_start+glitch_length+60)
			streamUVW_window.trim(starttime=glitch_start-60, endtime=glitch_start+glitch_length+60)
			 
			streamZNE_window.detrend('demean')
			streamUVW_window.detrend('demean')
			
			streamZNE_window.trim(starttime=glitch_start, endtime=glitch_start+glitch_length)
			streamUVW_window.trim(starttime=glitch_start, endtime=glitch_start+glitch_length)
			
			data_Z          = streamZNE_window[0].data
			data_N          = streamZNE_window[1].data
			data_E          = streamZNE_window[2].data
			data_U          = streamUVW_window[0].data
			data_V          = streamUVW_window[1].data
			data_W          = streamUVW_window[2].data

			ppol_results    = ppol_calc(data_N, data_E, data_Z, fix_angles='AMP', Xoffset_in_samples_for_amplitude=int(0.5*window/data_delta))[:-4] 	# no data, eigen vectors or eigen values return
			phi_2D          = ppol_results[0]
			phi_3D          = ppol_results[1]
			err_phi_2D      = ppol_results[2]
			err_phi_3D      = ppol_results[3]
			INCapp_2D       = ppol_results[4]
			err_INCapp_2D   = ppol_results[5]
			INCapp_3D       = ppol_results[6]
			err_INCapp_3D   = ppol_results[7]
			SNR_hori        = ppol_results[8]
			SNR_2D_radZ     = ppol_results[9]
			Rect_3D         = ppol_results[10]
			Rect_H          = ppol_results[11]
			Rect_RZ         = ppol_results[12]


			# GLITCH PARAMETERS
			glitch_counter += 1

			glitch_polar    = Rect_3D
			glitch_start    = glitch_start
			glitch_end      = glitch_start+glitch_length
			glitch_amp_Z    = data_Z[int(0.5*window/data_delta):][np.argmax(np.abs(data_Z[int(0.5*window/data_delta):]))]		# skip first half window for amplitude determination due to
			glitch_amp_N    = data_N[int(0.5*window/data_delta):][np.argmax(np.abs(data_N[int(0.5*window/data_delta):]))]		# occasional FIR precursor ringings at beginning of glitch!
			glitch_amp_E    = data_E[int(0.5*window/data_delta):][np.argmax(np.abs(data_E[int(0.5*window/data_delta):]))]
			glitch_amp_U    = data_U[int(0.5*window/data_delta):][np.argmax(np.abs(data_U[int(0.5*window/data_delta):]))]
			glitch_amp_V    = data_V[int(0.5*window/data_delta):][np.argmax(np.abs(data_V[int(0.5*window/data_delta):]))]
			glitch_amp_W    = data_W[int(0.5*window/data_delta):][np.argmax(np.abs(data_W[int(0.5*window/data_delta):]))]
			glitch_SNR      = SNR_2D_radZ
			glitch_pol      = phi_3D
			glitch_pol_err  = err_phi_3D
			glitch_inc      = INCapp_3D
			glitch_inc_err  = err_INCapp_3D

			glitch = glitch_counter, \
					 glitch_polar, \
					 glitch_start.strftime('%Y-%m-%dT%H:%M:%S'), \
					 glitch_end.strftime('%Y-%m-%dT%H:%M:%S'), \
					 glitch_amp_Z, \
					 glitch_amp_N, \
					 glitch_amp_E, \
					 glitch_amp_U, \
					 glitch_amp_V, \
					 glitch_amp_W, \
					 glitch_SNR, \
					 glitch_pol, \
					 glitch_pol_err, \
					 glitch_inc, \
					 glitch_inc_err  

			glitches.append(line_formatter % glitch)
			print(line_formatter % glitch)


			# for plotting, if wished
			if plot:
				glitch_starts  += [glitch_start]
				glitch_ends    += [glitch_end]


		# OUTPUT FILE
		try:
			np.savetxt(output_glitch, glitches, fmt='%s', header=header)
		except ValueError: # no glitches found
			pass


		# OUTPUT SHELL
		print()
		print(u'- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')
		print(u'Polarizations start: %s'               % polarizations_start)
		print(u'Polarizations end:   %s'               % polarizations_end)
		print(u'Polarizations #:     %s'               % polarizations_num)
		print()		
		print(u'Data start:          %s'               % data_start)
		print(u'Data end:            %s'               % data_end)
		print(u'Data length:         %s (h:m:s)'       % sec2hms( UTCDateTime(data_end)-UTCDateTime(data_start) ))
		print(u'Data SPS:            %sHz (%ss)'       % (1./data_delta, data_delta))
		print()
		print(u'Window length:       %ss'              % window)
		print(u'Window overlap:      %ss'              % window_overlap)
		print(u'Window step:         %ss (%s samples)' % (window_step,int(window_step/data_delta)))
		print(u'Number runs:         %s'               % len(polarizations))
		print()
		print(u'Threshold:           %s'               % threshold)
		print(u'Glitch length:       %ss'              % glitch_length)
		print(u'Glitches #:          %s'               % glitch_counter)
		print(u'- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')
		print(u'Done in:             %s (h:m:s).'      % sec2hms( float( time.time()-now )) )
		print(u'Input polarizations: %s'               % input_polar)
		print(u'Outfile glitches:    %s'               % output_glitch)


		### PLOT
		if plot:

			x    = np.linspace(0, 10, num=len(polarizations))
			y    = polarizations[:,12].astype('float')
			
			#f    = scipy.interpolate.interp1d(x, y, kind='linear')
			#xnew = np.linspace(0, 10, num=len(polarizations)*window_step/data_delta-window/data_delta)

			# hack, assign results to dummy stream object (so I have the datetimes), plotting data then works
			stream         = streamZNE[0:1].copy()
			#stream[0].data = f(xnew)
			stream[0].data = np.repeat(polarizations[:,12].astype('float'), window_step/data_delta)
			stream[0].stats.starttime = polarizations_start+window/2.

			# plotting
			streamZNE.normalise([-1,0])
			stream.normalise([0,1])
			quick_plot(*streamZNE, *stream, data_labels=(streamZNE[0].id, streamZNE[1].id, streamZNE[2].id, 'Rectinilinearity'), verts=glitch_starts, xlabel='Time', ylabel='Normalized Amplitude')


################  _ _ N A M E _ _ = = " _ _ M A I N _ _ "  ################
if __name__ == "__main__":

	# adjust line 401 for inventory file (needed for rotations!)
	# `seismic_file` is absolute path to waveforms containing (UVW or ZNE)


	## print documentation
	print( glitch_detector.__doc__ )


	## variables
	#seismic_file = '/home/scholz/Downloads/XB.68.SH.155.merge.mseed'
	#input_polar  = '/home/scholz/Downloads/polarization_XB.68.SH.155.merge.txt'	 # output of way 1 (hardcoded name in top of script)	

	seismic_file = '/home/scholz/Downloads/filt-test-merge.mseed'
	input_polar  = '/home/scholz/Downloads/polarization_filt-test-merge.txt'		 # output of way 1 (hardcoded name in top of script)
	

	## actual function
	glitch_detector(seismic_file, correct_gain=True, window=20, step_in_samples=10, input_polar=None)							# produces output `polarization_...`
	glitch_detector(seismic_file, correct_gain=True, input_polar=input_polar, threshold=0.97, glitch_length=25, plot=True)		# produces output `glitches_...`