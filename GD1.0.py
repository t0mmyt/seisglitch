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
def read(file=None):
	# wrapper to return Stream2 objects instead of (ObsPy's) Stream object.
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
		self.inventory    = self._get_inventory()
		self.filters      = {'0' : {'freqmin':None,  'freqmax':None, 'corners':3, 'zerophase':True, 'type_taper':'hann', 'max_percentage_taper':0.03, 'num':0, 'string':'0: None'},
		                     '1' : {'freqmin':1.0,   'freqmax':None, 'corners':3, 'zerophase':True, 'type_taper':'hann', 'max_percentage_taper':0.03, 'num':1, 'string':'1: >1 Hz'},
					         '2' : {'freqmin':0.1,   'freqmax':None, 'corners':3, 'zerophase':True, 'type_taper':'hann', 'max_percentage_taper':0.03, 'num':2, 'string':'2: >0.1 Hz'},
		                     '3' : {'freqmin':0.01,  'freqmax':None, 'corners':3, 'zerophase':True, 'type_taper':'hann', 'max_percentage_taper':0.03, 'num':3, 'string':'3: >0.01 Hz'},
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
	def trim_common(self):

		"""
		Trim traces in stream to common start- and endtime, i.e.,
		maximum starttime and minimum endtime
		"""
	
		max_starttime = max([tr.stats.starttime for tr in self ])
		min_endtime   = min([tr.stats.endtime   for tr in self ])				
		self.trim(starttime=max_starttime, endtime=min_endtime)	
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
	
		mint      = min([tr.stats.starttime for tr in self])
		maxt      = max([tr.stats.endtime for tr in self])
		starttime = mint + (maxt-mint)*start_per
		endtime   = maxt - (maxt-mint)*(1-end_per)
		self.trim(starttime=starttime, endtime=endtime)
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
	
		if isinstance(filter, str) or isinstance(filter, int):
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
	def _get_inventory(self):
		inv = Inventory()

		if re.search(r'BH|MH|SH', 'XXX'.join([tr.id for tr in self])):
			inv_file  = max( glob.glob( '/home/scholz/Downloads/dataless*.seed'), key=os.path.getctime) 	# latest download
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
	def rotate_2D(self, angle, components='NE', new_components='12', clockwise=False):

		"""
		"""

		if len(components) != 2:
			print('To rotate data in 2-D, you need to give two')
			print('component letters. You gave %s' % components)

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
		yy      = 1/np.tan(phi_2D*seismology.py/180) * xx
	
		# black lines
		ax.plot([-maxi, maxi], [0, 0], 'k', lw=1)
		ax.plot( [0, 0], [-maxi, maxi], 'k', lw=1)
		
		# data
		ax.scatter(trace_E.data, trace_N.data, s=7, c=colours)
		
		# arc for angle
		ax.plot( xx, yy, 'indianred', lw=1.5)
		ax.add_patch( Arc([0,0], maxi,  maxi, 90, -phi_2D, 0, color='indianred', lw=1.5, zorder=10 ))
		
		# arrows at end of arcs
		x  = maxi/2*np.sin( (phi_2D-7)*np.pi/180 )
		y  = maxi/2*np.cos( (phi_2D-7)*np.pi/180 )
		x2 = maxi/2*np.sin( (phi_2D+2)*np.pi/180 )
		y2 = maxi/2*np.cos( (phi_2D+2)*np.pi/180 )
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
			yy2 = 1/np.tan( INCapp_2D*np.pi/180) * xx2
		
			ax.plot([-maxi2, maxi2], [0, 0], 'k', lw=1)
			ax.plot( [0, 0], [-maxi2, maxi2], 'k', lw=1)
			ax.scatter(comp_R, trace_Z.data, s=5, c=colours)
			ax.plot( xx2, yy2, 'indianred', lw=1.5)
			ax.add_patch( Arc([0,0], maxi2,  maxi2, 90, -INCapp_2D, 0, color='indianred', lw=1.5, zorder=10 ))

			# arrows
			x  = maxi2/2*np.sin( (INCapp_2D-7)*np.pi/180 )
			y  = maxi2/2*np.cos( (INCapp_2D-7)*np.pi/180 )
			x2 = maxi2/2*np.sin( (INCapp_2D+2)*np.pi/180 )
			y2 = maxi2/2*np.cos( (INCapp_2D+2)*np.pi/180 )
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
			max_dist = np.sqrt(3*(maxi3*2)**2)

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
	
				if (evt.button==3 or mouse_clicks[stat_id]%2==0) and mouse_clicks[stat_id]!=0:	# right-click
					for i in indexes:
						for l in range(mouse_clicks[stat_id]):
							axes[i].lines[-1].remove()
							#axes[i].annotations[-1].remove()
					mouse_clicks[stat_id] = 0

				if evt.button == 1: 								 							# left-click
					mouse_clicks[stat_id] += 1	
					x = [evt.xdata, evt.xdata]
					for i in indexes:
						y = []
						for line in axes[i].lines:

							xdata = line.get_xdata()
							if xdata[0] == xdata[-1]:
								# avoid vertical lines that may be present already
								continue

							ydata      = line.get_ydata()
							ymin, ymax = min(ydata), max(ydata)
							y.append(ymin)
							y.append(ymax)

						margin = (max(y)-min(y))*1.05
						ylim   = [ min(y)-margin, max(y)+margin ]
						axes[i].plot(x, ylim, 'deepskyblue')
				fig.canvas.draw()
	

			else:																			### double click				
				if evt.button == 1:															# left-click (save and close figure)
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

			def instrument_removal(output, pre_filt=(0.001, 0.002, 50, 60)): 
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
						self.inventory = read_inventory('/home/scholz/Downloads/dataless.XB.ELYSE (4).seed')
					

					components = self._get_components()
					for stat_id in mouse_clicks.keys():
						corr_part  = corr.select(id=stat_id+'*')
						corr_part.merge()

						corr_part2 = corr.original.select(id=stat_id+'*')
						corr_part2.merge()
				
						if re.search(r'[U|V|W]', components):
							corr_part.remove_response(inventory=self.inventory, output=output, pre_filt=pre_filt)
							corr_part2.remove_response(inventory=self.inventory, output=output, pre_filt=pre_filt)

						elif re.search(r'[Z|N|E]', components):
							corr_part  = rotate2VBBUVW(corr_part, inventory=self.inventory)
							corr_part.remove_response(inventory=self.inventory, output=output, pre_filt=pre_filt)
							corr_part.rotate('->ZNE', inventory=self.inventory, components='UVW')

							corr_part2 = rotate2VBBUVW(corr_part2, inventory=self.inventory)
							corr_part2.remove_response(inventory=self.inventory, output=output, pre_filt=pre_filt)
							corr_part2.rotate('->ZNE', inventory=self.inventory, components='UVW')

				
					corr.filtering(self.current_filter_num)
					corr.plot(xlim=xlim)


			if evt.key.lower()=='a':							# Acceleration as target uni after instrument response removal 
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


			elif evt.key.lower()=='t':							# Trim file and save 
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
						if self.removed:
							trim = self.removed.copy()
						elif self.original:
							trim = self.original.copy()
						else:
							trim = self.copy()
						trim = trim.select(id=stat_id+'*')
						trim.trim(starttime=xlim[0], endtime=xlim[1])
						trim.write(filename)
						plt.close()
					
					fig      = plt.figure(figsize=(17,0.7), num='Save trimmed file')
					ax       = fig.add_subplot(111)
					
					propose  = os.path.join( os.path.dirname(self.origin), 'TRIM_'+os.path.basename(self.origin) )
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

# Plotting
def quick_plot(*y, x=[], data_labels=[], lw=1.5, win_title='', title='', xlabel='x-axis', ylabel='y-axis', x_invert=False, y_invert=False, xlim=(), ylim=(), verts=(), outfile=None, show=True):

	"""
	This function allows for some convenient, quick plotting.
	
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
	         then plotted as dates
	`data_labels` : - will appear in legend
	                - should have as many elements as you pass as `y`,
	                  otherwise it will not used but data simply numbered
	`lw`       : - linewidth of plots (lines are harcoded, no option to plot points only)
	`win_title`: - title of figures canvas (does not appear in plot)
	`title`    :     - supo-title of figure (appears in plot)
	`xlabel    : - will appear beneath x-axis
	`ylabel    : - will appear next to y-axis
	`x_invert` : - will invert direction of x-axis, if True
	`y_invert` : - will invert direction of y-axis, if True
	`verts`    : - either list or list of lists
	             - elements will be interpreted as x-data where vertical lines shall be plotted
	             - each list gets an own colour for all elements
	             - when elements are: - `UTCDateTime`, `datetime.datetime`, and `str` 
	                                     --> conversion to matplotlib times (so x-data should be times, too)
	                                  - `float` or `int`
	                                     --> where you would expect it to happen. If x-data are matplotlib,
	                                         you can pass here matplotlib times as well
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
	if not len(data_labels)==len(y):			# not same amount of entries (also true if none given)
		labels = []
		for i in range(len(y)):
			labels.append( 'Data %s' % (i+1) )
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