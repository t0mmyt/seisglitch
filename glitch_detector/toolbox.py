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


#####  python modules import  #####
import os
import re
import sys
import copy
import time
import glob
import scipy
import datetime
import numpy as np

#####  matplotlib modules import  #####
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.switch_backend('TKAgg')
import matplotlib.dates as mdates
import matplotlib.ticker as mticker
from matplotlib.widgets import TextBox
from matplotlib.patches import Arc, FancyArrowPatch
from mpl_toolkits.mplot3d import Axes3D, proj3d

#####  obspy modules import  #####
from obspy import read, read_inventory, UTCDateTime
from obspy.core.stream import Stream
from obspy.core.trace import Trace
from obspy.core.inventory import Inventory
from obspy.signal import rotate



# Extended Obspy Stream class (adding some more functionality)
def read2(file=None):
	# wrapper to make to return Stream2 objects instead of (ObsPy's) Stream object.
	st = read(file)
	st = Stream2(st, file=file)
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
		specific component, this script returns `False`. If it is unique, it returns `True`.

		This script can be used to make sure there is only one component within the passed stream.
		This is useful when checking for gaps/overlaps or for any other application
		where it's mandatory for correct processing to have only one component within the stream.
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
		Use Obspy's trim fucntion and improve slightly.
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

		super().trim(*args, **kwargs)
		self.times = self._get_times()
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
						for i in range( np.ceil(len(match)/truncate) ):
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
			print(u'GAIN CORRECTION APPLIED:')
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
			ids = sorted( set( [tr.id for tr in self] ), reverse=False)
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
	def _set_inventory(self, inventory=None):

		if not inventory:
			inventory = self._get_inventory()

		self.inventory = inventory
	def rotate_2D(self, angle, components='NE', new_components='12', clockwise=False):

		"""
		"""

		if len(components) != 2:
			print('To rotate data in 2-D, you need to give two')
			print('component letters. You gave %s' % len(components))

		self.trim_common()
		comp_1 = self.select2(component=components[0])
		comp_2 = self.select2(component=components[1])
		if comp_1 and comp_2:
			comp_1[0].data, comp_2[0].data = rotate_2D(comp_1[0].data, comp_2[0].data, angle, clockwise=clockwise)
			comp_1[0].stats.channel        = comp_1[0].stats.channel[:-1] + new_components[0]
			comp_2[0].stats.channel        = comp_2[0].stats.channel[:-1] + new_components[1]
			print('Rotated `components` (%s) to `new_components` (%s)' % (components, new_components))
			print('`angle` (%s), `clockwise` (%s)' % (angle, clockwise))
		else:
			print('No rotation performed, as `components` (%s)' % components)
			print('are not contained in stream.')
	def select2(self, *args, **kwargs):
		st_obs_select              = self.select(*args, **kwargs)
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
		stream_N = self.select(component='N') or self.select(component='Q') or self.select(component='U') or self.select(component='1') or self.select(component='T')
		stream_E = self.select(component='E') or self.select(component='T') or self.select(component='V') or self.select(component='2') or self.select(component='R')
		stream_Z = self.select(component='Z') or self.select(component='L') or self.select(component='W') or self.select(component='3')

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


		### Return, if wished
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
		yy      = 1/np.tan(phi_2D*np.pi/180) * xx
	
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
	def plot2(self, store_dir=os.getcwd(), store_name='*', verticals=(), type='normal', method='full', save_and_no_show=False, xlim=[], ylim=[], **kwargs):

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

# Convenient helpers
def sec2hms(seconds, digits=0):

	"""
	Convert seconds given as float into hours, minutes and seconds.
	Optional 'digits' determines position after decimal point.

	Returns string.
	"""


	string = str( datetime.timedelta(seconds=seconds) )
	parts  = string.split('.')

	try:
		frac_sec = parts[1]
		parts[1] = frac_sec[:digits]

	except IndexError:
		if digits>0:
			parts += [digits*'0']

	string = '.'.join(parts)
	string = string.rstrip('.')

	return string
def moving_window(*data, window_length_in_samples=100, step_in_samples=50, equal_end=True):

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


		yield_data = (array[index_window_start:index_window_end] for array in data)

		yield (index_window_start, index_window_end, *yield_data)

# Mathematical
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

		scale  = abs(scale_to_between[-1]-scale_to_between[0]) / 2.
		drange = max(data)-min(data)
		data   = data * 2 / drange
		data   = data - max(data)+1				# data have y values between [-1,1]
		data  *= scale							# data have y values between [-1/scale,1/scale] 
		data  += scale_to_between[-1]-scale 	# eacht trace has y values filling range `scale_to_between`
	
	else:
		data /= max(abs(data))

	return data
def snr(data, axis=0, ddof=1):

	"""
	Signal-to-noise ratio (SNR), as by Scipy.
	"""

	data = np.array(data)
	mean = data.mean(axis)
	sd_d = data.std(axis=axis, ddof=ddof)
	return np.where(sd_d==0, 0, mean/sd_d)

# InSight time conversions
def solify(UTC_time, sol0=UTCDateTime(2018, 11, 26, 5, 10, 50.33508)):

	"""
	:copyright:
	    Simon Stähler (mail@simonstaehler.com), 2018
	    Martin van Driel (Martin@vanDriel.de), 2018
	:license:
	    None
	"""

	SEC_PER_DAY_EARTH = 86400
	SEC_PER_DAY_MARS  = 88775.2440 #before: 88775.244147	

	MIT = (UTC_time - sol0) / SEC_PER_DAY_MARS
	t   = UTCDateTime((MIT - 1) * SEC_PER_DAY_EARTH)

	return t
def UTCify(LMST_time, sol0=UTCDateTime(2018, 11, 26, 5, 10, 50.33508)):
	"""
	:copyright:
	    Simon Stähler (mail@simonstaehler.com), 2018
	    Martin van Driel (Martin@vanDriel.de), 2018
	:license:
	    None
	"""
	SEC_PER_DAY_EARTH = 86400
	SEC_PER_DAY_MARS  = 88775.2440 #before: 88775.244147

	MIT      = float(LMST_time) / SEC_PER_DAY_EARTH + 1
	UTC_time = UTCDateTime(MIT * SEC_PER_DAY_MARS + float(sol0))

	return UTC_time
def sol2UTC(sol):
	# Convert a float, interpreted as InSight sol, to UTC.
	return UTCify(UTCDateTime('1970-01-01T00:00:00.000000Z')+datetime.timedelta(days=sol-1))
def ptime(time=None):

	"""
	PRINT TIME

	Small script to display current UTC time in most common ways.
	Pass a terrestrial time if you pass one.
	"""

	if time:
		print(u'TIME GIVEN')
		print(u'----------')
		if isinstance(time, float):										# MATLAB
			time = mdates.num2date(time)
		elif isinstance(time, int):										# sol
			time = sol2UTC(time).datetime
		elif isinstance(time, (UTCDateTime,datetime.datetime,str)):		# typical obspy str
			time = UTCDateTime(time).datetime
		else:															# others
			print(u'Passed time variable has no valid format (float, str, datetime.datetime, or UTCDateTime.)')
	else:
		print(u'TIME UTC NOW')
		print(u'------------')
		time = datetime.datetime.utcnow()

	print(u'MatLab:         %s'                % mdates.date2num(time))
	print(u'DateTime:       %s'                % time.__repr__())
	print(u"UTCDateTime:    UTCDateTime('%s')" % UTCDateTime(time))
	print(u'LocalMeanSolar: %sM%s'             % (solify(UTCDateTime(time)).julday, solify(UTCDateTime(time)).strftime('%H:%M:%S')))
def ltime(hms, sols_range=[], is_UTC=False):

	"""
	LIST TIME
	
	Small script to display a range of time. Useful e.g.
	when having to check data at a specific UTC or LMST
	each day / sol.

	`sol_range` can be a list or single numer, e.g.:
	   - sols_range=[17,193] (displays these sols)
	   - sols_range=17       (displays last 17 sols)
	If no `sols_range` specified, the last 10 sols are displayed by default.

	`is_UTC`=True means `hms` is given as UTC and you would like to
	see the corresponding LMST with respect to `sols_range`.
	`is_UTC`=False means `hms` is given as LMST and you would like to
	see the corresponding UTC with respect to `sols_range`.
	"""

	## VARIABLES
	hms = '%06d' % ( int(hms)*10**(6-len(str(hms))) )

	if isinstance(sols_range, (float, int)):
		sols_range = [sols_range]
	if not sols_range or len(sols_range)==1:
		UTC_now    = UTCDateTime( time.time() )
		LMST_now   = solify(UTC_now)

		if len(sols_range)==1:
			sols_range = [LMST_now.julday-sols_range[0], LMST_now.julday-1]
		else:
			sols_range = [LMST_now.julday-10,            LMST_now.julday-1]


	## PRINTS
	print('UTC                    LMST')
	print('---                    ----')	
	for sol in range(sols_range[0], sols_range[1]+1):

		if is_UTC:
			time_str_ymd = sol2UTC(sol).strftime('%Y-%m-%d')
			time_str_HMS = '%s:%s:%s' % (hms[0:2], hms[2:4], hms[4:6])

			UTC_time     = UTCDateTime( time_str_ymd + 'T' + time_str_HMS)
			LMST_time    = solify(UTC_time)


		else:
			LMST_time    = UTCDateTime('1970-01-01T%s:%s:%s.000000Z' % (hms[:2], hms[2:4], hms[4:])) + datetime.timedelta(days=sol)
			UTC_time     = UTCify(LMST_time)


		print('%s    %sS%s' % (UTC_time.strftime('%Y-%m-%dT%H:%M:%S'), LMST_time.julday, LMST_time.strftime('%H:%M:%S')))

# Illustrational
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
class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)


################  _ _ N A M E _ _ = = " _ _ M A I N _ _ "  ################
if __name__ == "__main__":
	argues = sys.argv
	eval(argues[1])
	#print('Define Testing')