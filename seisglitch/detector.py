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


#####  python modules import  #####
import os
import re
import sys
import time
import numpy as np


#####  obspy modules import  #####
from obspy import UTCDateTime


#####  seisglitch modules import  #####
from seisglitch.util import read2, Stream2, moving_window, normalise
from ppol.core import ppol_calc
from ppol.util import quick_plot, sec2hms



### GLITCH DETECTION: based on a step function convoluted with low-pass filtered acceleration data
def glitch_detector(*RAW_UVW, stepwin_length_total=80, ACCfilter=3, threshold=1e-7, glitch_mindist=5, glitch_length=25, plot=False):
	

	"""
	SETTINGS:
		VBB: - stepwin_length_total = 100
		     - threshold = 1e-7
		 SP: - stepwin_length_total = 150
		     - threshold = 5e-6

	https://stackoverflow.com/questions/48000663/step-detection-in-one-dimensional-data

	Hard-coded FIR pre-cursor offset of 1 s.
	"""


	now = time.time()



	### SMALL OUTPUT
	print()
	print(u'  -----------------------')
	print(u'  RUNNING GLITCH DETECTOR')
	print(u'  -----------------------')



	### DOING THE SAME ANALYSIS / OUTPUT FOR ALL PASSED FILES
	for i in range(len(RAW_UVW)):
	

		### READ DATA and assign output file
		stRAW         = read2(RAW_UVW[i])
		stRAW.truncate(0.4,0.7)
		output_glitch = os.path.join( os.path.dirname(RAW_UVW[i]), 'glitches_' + '.'.join(os.path.basename(RAW_UVW[i]).split('.')[:-1]) + '.txt' )
			

		### SMALL OUTPUT
		print()
		print()
		print(u'Analysing file: %s/%s' % (i+1,len(RAW_UVW)))
		print()
		print(u'INPUT STREAM:')
		print( str(stRAW).replace('\n','\n  '))


		### PREPARE DATA, YA' KNOW, RAW & DIS & VEL & ACC as UVW & ZNE ..
		stRAW.sort(reverse=False)	# if you change this, god may stand by your side
		stRAW.trim_common() 		# all start with very same time (not necessary, but still desirable)
		stRAW._set_inventory()		# get & set station inventory for stream (response information)
		data_start, data_end = stRAW.times
		data_delta           = stRAW[0].stats.delta

		stGAI = stRAW.copy()
		stGAI.gain_correction()
		stZNE = stGAI.copy()
		stZNE.rotate('->ZNE', inventory=stZNE.inventory, components=('UVW'))
		stZNE.sort(reverse=True)
		stGAI += stZNE  			# having in GAI (overall sensitivy and gain corrected): UVWZNE (in that order)

		stDIS = stRAW.copy()
		stDIS.remove_response(inventory=stRAW.inventory, output='DISP', pre_filt=None, water_level=60)
		stZNE = stDIS.copy()
		stZNE.rotate('->ZNE', inventory=stZNE.inventory, components=('UVW'))
		stZNE.sort(reverse=True)
		stDIS += stZNE  			# having in DIS: UVWZNE (in that order)

		stVEL = stRAW.copy()
		stVEL.remove_response(inventory=stRAW.inventory, output='VEL',  pre_filt=None, water_level=60)
		stZNE = stVEL.copy()
		stZNE.rotate('->ZNE', inventory=stZNE.inventory, components=('UVW'))
		stZNE.sort(reverse=True)
		stVEL += stZNE  			# having in VEL: UVWZNE (in that order)

		stACC = stRAW.copy()
		stACC.remove_response(inventory=stRAW.inventory, output='ACC',  pre_filt=None, water_level=60)
		stACC_fil = stACC.copy()
		stACC_fil.filtering(ACCfilter)	# LP <1 Hz to better see step in acceleration
		stACC_fil.filtering(4)		    # DOES IT MAKE A DIFFERENCE FOR GLITCH DETECTIONS?
		stZNE = stACC.copy()
		stZNE.rotate('->ZNE', inventory=stZNE.inventory, components=('UVW'))
		stZNE.sort(reverse=True)
		stACC += stZNE  			# having in ACC: UVWZNE (in that order)

		stCON = stRAW.copy()  		# dummy container having to put convolution calculations into 


		st = Stream2(stRAW+stGAI+stDIS+stVEL+stACC+stACC_fil+stCON)
		st.taper(0.05) 				# SHOULD PERHAPS MAKE AT BEGINNING FOR stRAW (to ease response effects at edges of data?!)


		### CONVOLUTE WITH STEP FUNCTION FOR EACH OF UVW-COMPONENTS, EXTRACT GLITCH-STARTS FOR EACH COMPONENT INDIVIDUALLY (!)
		glitch_start_times = [[] for i in range(len(stRAW))]
		convolutions       = [[] for i in range(len(stRAW))]
		for l, trace in enumerate(stACC_fil):

			step        = np.hstack(( np.ones(stepwin_length_total//2), -1*np.ones(stepwin_length_total//2) )) # 
			convolve    = np.convolve(trace.data, step, mode='valid')
			
	        #indices     = np.hstack((scipy.signal.argrelextrema(convolve, np.greater)[0], scipy.signal.argrelextrema(convolve, np.less)[0]))
			indices     = np.where((convolve >= threshold) | (convolve <= -threshold))
			indices_pos = [index for index in indices[0] if convolve[index-1]<=convolve[index] and convolve[index]>=convolve[index+1]]
			indices_neg = [index for index in indices[0] if convolve[index-1]>=convolve[index] and convolve[index]<=convolve[index+1]]
			
			new_indices = []
			for index in indices_pos:
				for i in range(index, len(convolve)):
					try:
						condition = (convolve[i-1]>=convolve[i] and convolve[i]<=convolve[i+1])
					except:
						pass
					if condition:
						new_indices.append( i )
						break
			for index in indices_neg:
				for i in range(index, len(convolve)):
					try:
						condition = (convolve[i-1]<=convolve[i] and convolve[i]>=convolve[i+1])
					except:
						pass
					if condition:
						new_indices.append( i )
						break
			new_indices.sort()

			try:
				start_times           = trace.times('utcdatetime')[new_indices]
				start_times           = [start_times[0]-(glitch_mindist+0.1)] + list(start_times)
				start_times           = [start_times[i] for i in range(1, len(start_times)) if start_times[i-1]+glitch_mindist<=start_times[i]]
				
				glitch_start_times[l] = start_times
				stCON[l].data         = convolve

			except IndexError:	# =empty lists ==> no glitches found
				print(u'No glitches found for given parameters.')
				sys.exit()


		### UNIFY GLITCH START TIMES ACROSS UVW-COMPONENTS (using parameter: `glitch_mindist` given in seconds)
		glitch_start_times = np.array(glitch_start_times)
		flat_starts        = np.sort( np.hstack(glitch_start_times) )	# np.flatten wouldn't work because sublists contain different number of elements ..
		
		glitch_starts      = {}
		shifter_avoid      = 0

		for l in range( len( flat_starts )):

			k = l + shifter_avoid

			try:
				indices_same_glitch = np.where( (flat_starts < flat_starts[k]+glitch_mindist) & (flat_starts > flat_starts[k]-glitch_mindist) )
			except IndexError:		# composed index `k` >= len(flat_starts) --> end of data basically
				break

			same_glitch_starts   = flat_starts[indices_same_glitch]
			unified_glitch_start = min( same_glitch_starts )
			unified_glitch_comps = []
			for gs in same_glitch_starts:
				for m, comp_starts in enumerate(glitch_start_times):
					if gs in comp_starts and not str(m) in unified_glitch_comps:
						unified_glitch_comps.append( str(m) )
						break

			unified_glitch_comps = ''.join( sorted( unified_glitch_comps ) )
			unified_glitch_comps = unified_glitch_comps.replace('0','U').replace('1','V').replace('2','W')

			try:
				glitch_starts[str(unified_glitch_start)] += unified_glitch_comps
			except KeyError: 	# not yet in dictionary
				glitch_starts[str(unified_glitch_start)]  = unified_glitch_comps

			#print(k, shifter_avoid)
			#print(same_glitch_starts)
			#print(unified_glitch_start)
			#print(unified_glitch_comps)
			#print()
			shifter_avoid += len(same_glitch_starts)-1

		glitch_starts       = np.array( [[time, ''.join(sorted(set(glitch_starts[time]))) ] for time in glitch_starts.keys()] )
		glitch_starts       = glitch_starts[np.argsort(glitch_starts[:,0])]
		glitch_starts_all   = glitch_starts[:,0]
		glitch_starts_U     = [time_comp[0] for time_comp in glitch_starts if 'U' in time_comp[1]]
		glitch_starts_V     = [time_comp[0] for time_comp in glitch_starts if 'V' in time_comp[1]]
		glitch_starts_W     = [time_comp[0] for time_comp in glitch_starts if 'W' in time_comp[1]]
		glitch_starts_Uonly = [time_comp[0] for time_comp in glitch_starts if 'U' in time_comp[1] and len(time_comp[1])==1]
		glitch_starts_Vonly = [time_comp[0] for time_comp in glitch_starts if 'V' in time_comp[1] and len(time_comp[1])==1]
		glitch_starts_Wonly = [time_comp[0] for time_comp in glitch_starts if 'W' in time_comp[1] and len(time_comp[1])==1]


		### CALCULATE GLITCH PARAMETERS (e.g. TIMES, AMPLITUDES, BACKAZIMUTHS, INCIDENCE ANGLES, SNRs, POLARZATIONS)
		header         = u'   0                    1                    2         3         4         5          6          7          8          9         10         11         12         13         14         15         16         17         18         19         20         21         22         23         24         25         26         27         28         29         30         31         32          33          34          35          36          37          38          39          40          41          42          43          44          45          46          47          48\n' \
		                 u' NUM         GLITCH-START           GLITCH-END  U-GLITCH  V-GLITCH  W-GLITCH      U-RAW      V-RAW      W-RAW      U-GAI      V-GAI      W-GAI      Z-GAI      N-GAI      E-GAI      U-DIS      V-DIS      W-DIS      Z-DIS      N-DIS      E-DIS      U-VEL      V-VEL      W-VEL      Z-VEL      N-VEL      E-VEL      U-ACC      V-ACC      W-ACC      Z-ACC      N-ACC      E-ACC  PHI_3D_GAI  INC_3D_GAI  SNR_3D_GAI  POL_3D_GAI  PHI_3D_DIS  INC_3D_DIS  SNR_3D_DIS  POL_3D_DIS  PHI_3D_VEL  INC_3D_VEL  SNR_3D_VEL  POL_3D_VEL  PHI_3D_ACC  INC_3D_ACC  SNR_3D_ACC  POL_3D_ACC\n' \
		                 u'----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'
		line_formatter = u'%05d:  %15s  %15s  %8d  %8d  %8d  %9d  %9d  %9d  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %10.1f  %10.1f  %10.3f  %10.3g  %10.1f  %10.1f  %10.3f  %10.3g  %10.1f  %10.1f  %10.3g  %10.3f  %10.1f  %10.1f  %10.3g  %10.3f'

		glitches       = []
		glitch_counter = 0
		for glitch_start, comps in glitch_starts:

			# prepare data outputs
			glitch_start    = UTCDateTime(glitch_start) 
			glitch_end      = glitch_start+glitch_length
			glitch_counter += 1
		
			# prepare different streams
			stRAW_window = stRAW.copy()
			stGAI_window = stGAI.copy()
			stDIS_window = stDIS.copy()
			stVEL_window = stVEL.copy()
			stACC_window = stACC.copy()
			
			stRAW_window.trim(starttime=glitch_start-60, endtime=glitch_start+glitch_length+60)
			stGAI_window.trim(starttime=glitch_start-60, endtime=glitch_start+glitch_length+60)
			stDIS_window.trim(starttime=glitch_start-60, endtime=glitch_start+glitch_length+60)
			stVEL_window.trim(starttime=glitch_start-60, endtime=glitch_start+glitch_length+60)
			stACC_window.trim(starttime=glitch_start-60, endtime=glitch_start+glitch_length+60)
			 
			stRAW_window.detrend('demean')
			stGAI_window.detrend('demean')
			stDIS_window.detrend('demean')
			stVEL_window.detrend('demean')
			stACC_window.detrend('demean')

			# trim streams further to glitch window for amplitude extractions
			stRAW_window.trim(starttime=glitch_start, endtime=glitch_start+glitch_length)
			stGAI_window.trim(starttime=glitch_start, endtime=glitch_start+glitch_length)
			stDIS_window.trim(starttime=glitch_start, endtime=glitch_start+glitch_length)
			stVEL_window.trim(starttime=glitch_start, endtime=glitch_start+glitch_length)
			stACC_window.trim(starttime=glitch_start, endtime=glitch_start+glitch_length)

			data_U_RAW = stRAW_window[0].data
			data_V_RAW = stRAW_window[1].data
			data_W_RAW = stRAW_window[2].data
			data_U_GAI = stGAI_window[0].data
			data_V_GAI = stGAI_window[1].data
			data_W_GAI = stGAI_window[2].data
			data_Z_GAI = stGAI_window[3].data
			data_N_GAI = stGAI_window[4].data
			data_E_GAI = stGAI_window[5].data
			data_U_DIS = stDIS_window[0].data
			data_V_DIS = stDIS_window[1].data
			data_W_DIS = stDIS_window[2].data
			data_Z_DIS = stDIS_window[3].data
			data_N_DIS = stDIS_window[4].data
			data_E_DIS = stDIS_window[5].data
			data_U_VEL = stVEL_window[0].data
			data_V_VEL = stVEL_window[1].data
			data_W_VEL = stVEL_window[2].data
			data_Z_VEL = stVEL_window[3].data
			data_N_VEL = stVEL_window[4].data
			data_E_VEL = stVEL_window[5].data
			data_U_ACC = stACC_window[0].data
			data_V_ACC = stACC_window[1].data
			data_W_ACC = stACC_window[2].data
			data_Z_ACC = stACC_window[3].data
			data_N_ACC = stACC_window[4].data
			data_E_ACC = stACC_window[5].data


			# glitch flag per UVW-component 
			if 'U' in comps:
				U_GLITCH = 1
			else:
				U_GLITCH = 0
			
			if 'V' in comps:
				V_GLITCH = 1
			else:
				V_GLITCH = 0

			if 'W' in comps:
				W_GLITCH = 1
			else:
				W_GLITCH = 0


			# get UVWZNE amplitudes of glitches
			offset = int(1./data_delta)	# 1 second (in samples) offset to avoid FIR-precursors 
			
			U_RAW  = data_U_RAW[offset:][np.argmax(np.abs(data_U_RAW[offset:]))]
			U_GAI  = data_U_GAI[offset:][np.argmax(np.abs(data_U_GAI[offset:]))]
			U_DIS  = data_U_DIS[offset:][np.argmax(np.abs(data_U_DIS[offset:]))]
			U_VEL  = data_U_VEL[offset:][np.argmax(np.abs(data_U_VEL[offset:]))]
			U_ACC  = data_U_ACC[offset:][np.argmax(np.abs(data_U_ACC[offset:]))]
			
			V_RAW  = data_V_RAW[offset:][np.argmax(np.abs(data_V_RAW[offset:]))]
			V_GAI  = data_V_GAI[offset:][np.argmax(np.abs(data_V_GAI[offset:]))]
			V_DIS  = data_V_DIS[offset:][np.argmax(np.abs(data_V_DIS[offset:]))]
			V_VEL  = data_V_VEL[offset:][np.argmax(np.abs(data_V_VEL[offset:]))]
			V_ACC  = data_V_ACC[offset:][np.argmax(np.abs(data_V_ACC[offset:]))]
			
			W_RAW  = data_W_RAW[offset:][np.argmax(np.abs(data_W_RAW[offset:]))]
			W_GAI  = data_W_GAI[offset:][np.argmax(np.abs(data_W_GAI[offset:]))]
			W_DIS  = data_W_DIS[offset:][np.argmax(np.abs(data_W_DIS[offset:]))]
			W_VEL  = data_W_VEL[offset:][np.argmax(np.abs(data_W_VEL[offset:]))]
			W_ACC  = data_W_ACC[offset:][np.argmax(np.abs(data_W_ACC[offset:]))]
			
			Z_GAI  = data_Z_GAI[offset:][np.argmax(np.abs(data_Z_GAI[offset:]))]
			Z_DIS  = data_Z_DIS[offset:][np.argmax(np.abs(data_Z_DIS[offset:]))]
			Z_VEL  = data_Z_VEL[offset:][np.argmax(np.abs(data_Z_VEL[offset:]))]
			Z_ACC  = data_Z_ACC[offset:][np.argmax(np.abs(data_Z_ACC[offset:]))]
			
			N_GAI  = data_N_GAI[offset:][np.argmax(np.abs(data_N_GAI[offset:]))]
			N_DIS  = data_N_DIS[offset:][np.argmax(np.abs(data_N_DIS[offset:]))]
			N_VEL  = data_N_VEL[offset:][np.argmax(np.abs(data_N_VEL[offset:]))]
			N_ACC  = data_N_ACC[offset:][np.argmax(np.abs(data_N_ACC[offset:]))]
			
			E_GAI  = data_E_GAI[offset:][np.argmax(np.abs(data_E_GAI[offset:]))]
			E_DIS  = data_E_DIS[offset:][np.argmax(np.abs(data_E_DIS[offset:]))]
			E_VEL  = data_E_VEL[offset:][np.argmax(np.abs(data_E_VEL[offset:]))]
			E_ACC  = data_E_ACC[offset:][np.argmax(np.abs(data_E_ACC[offset:]))]

			# get glitch BAZ, INC, and POL (=rectinilinearity)
			PHI_3D_GAI, INC_3D_GAI, SNR_3D_GAI, POL_3D_GAI = np.asarray( ppol_calc(data_N_GAI, data_E_GAI, data_Z_GAI, fix_angles='AMP', Xoffset_in_samples_for_amplitude=offset) )[[2,6,9,13]]
			PHI_3D_DIS, INC_3D_DIS, SNR_3D_DIS, POL_3D_DIS = np.asarray( ppol_calc(data_N_DIS, data_E_DIS, data_Z_DIS, fix_angles='AMP', Xoffset_in_samples_for_amplitude=offset) )[[2,6,9,13]]
			PHI_3D_VEL, INC_3D_VEL, SNR_3D_VEL, POL_3D_VEL = np.asarray( ppol_calc(data_N_VEL, data_E_VEL, data_Z_VEL, fix_angles='AMP', Xoffset_in_samples_for_amplitude=offset) )[[2,6,9,13]]
			PHI_3D_ACC, INC_3D_ACC, SNR_3D_ACC, POL_3D_ACC = np.asarray( ppol_calc(data_N_ACC, data_E_ACC, data_Z_ACC, fix_angles='AMP', Xoffset_in_samples_for_amplitude=offset) )[[2,6,9,13]]

			# get all glitch measures together
			glitch = glitch_counter, \
					 glitch_start.strftime('%Y-%m-%dT%H:%M:%S'), \
					 glitch_end.strftime('%Y-%m-%dT%H:%M:%S'),\
					 U_GLITCH, \
					 V_GLITCH, \
					 W_GLITCH, \
					 U_RAW, \
					 V_RAW, \
					 W_RAW, \
					 U_GAI, \
					 V_GAI, \
					 W_GAI, \
					 Z_GAI, \
					 N_GAI, \
					 E_GAI, \
					 U_DIS, \
					 V_DIS, \
					 W_DIS, \
					 Z_DIS, \
					 N_DIS, \
					 E_DIS, \
					 U_VEL, \
					 V_VEL, \
					 W_VEL, \
					 Z_VEL, \
					 N_VEL, \
					 E_VEL, \
					 U_ACC, \
					 V_ACC, \
					 W_ACC, \
					 Z_ACC, \
					 N_ACC, \
					 E_ACC, \
					 PHI_3D_GAI, \
					 INC_3D_GAI, \
					 SNR_3D_GAI, \
					 POL_3D_GAI, \
					 PHI_3D_DIS, \
					 INC_3D_DIS, \
					 SNR_3D_DIS, \
					 POL_3D_DIS, \
					 PHI_3D_VEL, \
					 INC_3D_VEL, \
					 SNR_3D_VEL, \
					 POL_3D_VEL, \
					 PHI_3D_ACC, \
					 INC_3D_ACC, \
					 SNR_3D_ACC, \
					 POL_3D_ACC

			# prapare file output and print into terminal
			glitches.append(line_formatter % glitch)
			#print( (line_formatter % glitch))


		### FILING GLITCHES
		np.savetxt(output_glitch, glitches, fmt='%s', header=header)


		### SHELLING SUMMARY
		print()
		print(u'- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')
		print(u'Data start:          %s'               % data_start)
		print(u'Data end:            %s'               % data_end)
		print(u'Data length:         %s (h:m:s)'       % sec2hms( data_end-data_start ))
		print(u'Data SPS:            %sHz (%ss)'       % (1./data_delta, data_delta))
		print(u'Data path:           %s'               % RAW_UVW[i])
		print()
		print(u'Step window length:  %s samples'       % len(step))
		print(u'Acceleration filter: %s'               % stACC_fil._get_filter_str())
		print(u'Threshold:           %s'               % threshold)
		print(u"Glitches' min dist:  %ss"              % glitch_mindist)
		print(u"Glitch length (fix): %ss"              % glitch_length)
		print(u'Glitches ALL:        %s'               % glitch_counter)
		print(u'        on U:        %s'               % len(glitch_starts_U))
		print(u'        on V:        %s'               % len(glitch_starts_V))
		print(u'        on W:        %s'               % len(glitch_starts_W))
		print(u'   only on U:        %s'               % len(glitch_starts_Uonly))
		print(u'   only on V:        %s'               % len(glitch_starts_Vonly))
		print(u'   only on W:        %s'               % len(glitch_starts_Wonly))
		print(u'- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')
		print(u'Done in:             %s (h:m:s).'      % sec2hms( time.time()-now ))
		print(u'Outfile glitches:    %s'               % output_glitch)


		### PLOTTING, if wished
		if plot:
			#stCON.normalise([0,1])
			#stRAW.normalise([mi, ma])
			#stRAW.filtering(1)
			#stDIS.filtering(1)
			#stVEL.filtering(1)
			#stACC.filtering(1)

			#stCON.filtering(4)
			stRAW.filtering(4)
			stDIS.filtering(4)
			stVEL.filtering(4)
			stACC.filtering(4)
			#stACC_fil.filtering(4)

			mi, ma = min(stCON[2].data), max(stCON[2].data)
			stACC_fil.normalise([mi+mi, ma+mi])
			stCON.normalise([mi, ma])
			stRAW.normalise([mi-1*mi, ma-1*mi])
			stDIS.normalise([mi-2*mi, ma-2*mi])
			stVEL.normalise([mi-3*mi, ma-3*mi])
			stACC.normalise([mi-4*mi, ma-4*mi])

			quick_plot(stCON[2], stACC_fil[2], stRAW[2], stDIS[2], stVEL[2], stACC[2], data_labels=('W CON','W ACC-FIL','W RAW','W DIS','W VEL','W ACC'), title='%s glitches' % glitch_counter, verts=glitch_starts_W, xlabel='Time', ylabel='Amplitude', outfile=output_plot)
			#quick_plot(stCON[2], stACC_fil[2], stACC_fil2[2],  stRAW[2], data_labels=('W CON', 'W ACC FIL', 'W ACC FIL2', 'W RAW'), title='%s glitches' % glitch_counter, verts=glitch_starts_W, xlabel='Time', ylabel='Amplitude')
			#quick_plot(stACC_fil[0], stACC_fil2[0],  stRAW[0], data_labels=('U ACC FIL', 'U ACC FIL2', 'U RAW'), title='%s glitches' % glitch_counter, verts=glitch_starts_U, xlabel='Time', ylabel='Amplitude')
			#quick_plot(*stCON, *stACC_fil, *stRAW, data_labels=('U CON', 'V CON', 'W CON', 'U ACC FIL', 'V ACC FIL', 'W ACC FIL', 'U RAW', 'V RAW', 'W RAW'), title='%s glitches' % glitch_counter, verts=glitch_starts_all, xlabel='Time', ylabel='Amplitude' )
			print(u'Outfile plot:        %s' % output_plot)



### _ _ N A M E _ _ = = " _ _ M A I N _ _ "  
if __name__ == "__main__":
	pass