#!/usr/bin/env python
# -*- coding: utf-8 -*-

#----------------------------------------------------------------------
#   Filename:  main.py
""" Useful functions for glitch detection, including plot scripts. """
#   Author:    John - Robert Scholz
#   Email:     john.robert.scholz@gmail.com
#   Date:      Jan 2019 - present
#   License:   TBD
#---------------------------------------------------------------------


from toolbox import Stream2, moving_window, snr, ppol_calc, solify, UTCify, quick_plot


### two different glitch detectors
def glitch_detector1(RAW_UVW, window=5, step_in_samples=10, input_polar=None, threshold=0.97, glitch_length=25, plot=False):

	"""

	---------------------------------------------------------
	|                                                       |
	|   Hello ! Are you sure about getting into glitches?   |
	|                                                       |
	---------------------------------------------------------


	This script can be run in two different ways.
	For both ways, the script expects `seismic_file` (absolute path)
	to contain at least three seismic components (UVW or ZNE).
	No matter which one is passed, two streams are created; one with ZNE &
	one with UVW. Decide if you want `correct_gain` (True or False).
	Note: the file is read as an ObsPy `Stream2` object - a self-written class.
	Note: for rotating the components, the `_get_inventory()` method
	(see Stream2 class) attempts to read an inventory file that follows
	are certain naming in a certain folder. Please adjust line 421!

	Way 1): `input_polar` = None

		Using a 3-D PCA (function `ppol_moving`), the rectinilinearity (polarisation)
		and other measures are calculated for each moving window of length `window`, 
		in seconds, stepping `step_in_samples`. These polarisations are written into 
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
	stUVW         = read2(RAW_UVW)
	stUVW._set_inventory()
	stUVW.filtering(4)
	output_polar  = os.path.join( os.path.dirname(stUVW.origin), 'polarization_' + '.'.join(os.path.basename(stUVW.origin).split('.')[:-1]) + '.txt' )
	output_glitch = os.path.join( os.path.dirname(stUVW.origin), 'glitches_'     + '.'.join(os.path.basename(stUVW.origin).split('.')[:-1]) + '.txt' )
	output_plot   = os.path.join( os.path.dirname(stUVW.origin), 'plot_'         + '.'.join(os.path.basename(stUVW.origin).split('.')[:-1]) + '.png' )



	stZNE = stUVW.copy()
	streamZNE.rotate('->ZNE', inventory=stZNE.inventory, components='UVW')
	stZNE.sort(reverse=True)	# sorted: Z, N, E
	stUVW.sort(reverse=False)	# sorted: U, V, W
	print()
	print('INPUT STREAMS:')
	print(streamUVW)
	print(streamZNE)


	### CREATE POLARIZATION DATA TO LATER DETECT GLITCHES ON 
	if not input_polar:


		# SMALL SHELL OUTPUT
		print()
		print(u'DATA ANALYSIS:  %s - %s' % (streamVBB[0].stats.starttime, streamVBB[0].stats.endtime))
		print(u'Running moving windows of P-wave polarization ..')


		# MOVING P-WAVE POLARIZATION
		stZNE.trim_common()
		stream_string  = stZNE[0].id[:-1] + '?'


		# EXTRACT VARIABLES NEEDED FROM STREAM OBJECT
		starttime      = stZNE[0].stats.starttime
		endtime        = stZNE[0].stats.endtime
		duration       = endtime - starttime
		samples        = int( (endtime-starttime)/delta+1)
		delta          = stZNE[0].stats.delta
		SPS            = stZNE[0].stats.sampling_rate
		
		window_step    = step_in_samples*delta
		window_overlap = window-window_step
		window_counter = 0
		
		#data_Z = stZNE.select(component='[Z]'  )[0].data
		#data_N = stZNE.select(component='[N1X]')[0].data
		#data_E = stZNE.select(component='[E2Y]')[0].data
		data_Z = stZNE[0].data
		data_N = stZNE[1].data
		data_E = stZNE[2].data	
		

		# LOOP (MOVING WINDOWS)
		pol_values     = []
		line_formatter = '%6.2f %6.2f %5.2f %5.2f %6.2f %5.2f %6.2f %5.2f %7.1f %7.1f %5.3f %5.3f %5.3f %23s %23s'

		for mwindow in moving_window(data_Z, data_N, data_E, window_length_in_samples=window*SPS, step_in_samples=step_in_samples, equal_end=True):

			window_index_start = mwindow[0]
			window_index_end   = mwindow[1]
			window_Z           = mwindow[2]
			window_N           = mwindow[3]
			window_E           = mwindow[4]
			
			start              = starttime + window_index_start*delta
			end                = start + window

			# P-POL RESULT FOR WINDOW
			window_pol = ppol_calc(window_N, window_E, window_Z)[:13] + (start, end)
			pol_values.append( line_formatter % window_pol )

			# small output into shell
			window_counter += 1
			if window_counter%100==0:
				print(u'Window %7d:  %s - %s' % (window_counter, start, end))		


		# OUTPUT FILE
		header  = 'Polarizations  %s  (DATA_RANGE:  %s  %s)\n' % (stream_string, starttime, endtime)
		header += 'PH2D   PH3D 2DERR 3DERR  INC2D 2DERR  INC3D 3DERR SNR_HOR  SNR_RZ POL3D  POLH POLRZ                       START                         END'
		np.savetxt(output_polar, pol_values, header=header, fmt='%s', delimiter=None)


		# OUTPUT SHELL, IF WISHED
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
		print(u'Windows total:    %s'          % window_counter)
		print(u'- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')
		print(u'Done in:          %s (h:m:s).' % sec2hms( float( time.time()-now )) )
		print(u'Output polar.:    %s'          % output_polar)



	### DETECT GLITCHES ON ALREADY CALCULATED POLARIZATIONS
	else:


		# assign needed variables
		data_start                 = str(stZNE[0].stats.starttime)
		data_end                   = str(stZNE[0].stats.endtime)
		data_delta                 = stZNE[0].stats.delta
		
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
		header                     = ' NUM    POL                START                  END      Z-AMP      N-AMP      E-AMP      U-AMP      V-AMP      W-AMP      SNR    DIR  DERR    INC  IERR'
		line_formatter             = '%05d:  %5.3f  %15s  %15s  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %7.1f  %5.1f  %4.1f  %5.1f  %4.1f'
		print('# '+header)
		

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
				if not (pol_value_b<threshold and pol_value>=threshold and pol_start>=glitch_end and pol_mid+glitch_length<=polarizations_end):
					raise Exception
			except Exception:
				continue 


			# PREPARE DATA FOR GLITCH PARAMETERS
			glitch_start     = pol_mid
			
			streamZNE_window = stZNE.copy()
			streamUVW_window = stUVW.copy()
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

			ppol_results    = ppol_calc(data_N, data_E, data_Z, fix_angles='AMP', Xoffset_in_samples_for_amplitude=int(.25*glitch_length/data_delta))[:-4] 	# no data, eigen vectors or eigen values return
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
			glitch_amp_Z    = data_Z[int(.25*glitch_length/data_delta):][np.argmax(np.abs(data_Z[int(.25*glitch_length/data_delta):]))]		# skip first quarter of glitch length for amplitude determination due to
			glitch_amp_N    = data_N[int(.25*glitch_length/data_delta):][np.argmax(np.abs(data_N[int(.25*glitch_length/data_delta):]))]		# occasional FIR precursor ringings at beginning of glitch!
			glitch_amp_E    = data_E[int(.25*glitch_length/data_delta):][np.argmax(np.abs(data_E[int(.25*glitch_length/data_delta):]))]
			glitch_amp_U    = data_U[int(.25*glitch_length/data_delta):][np.argmax(np.abs(data_U[int(.25*glitch_length/data_delta):]))]
			glitch_amp_V    = data_V[int(.25*glitch_length/data_delta):][np.argmax(np.abs(data_V[int(.25*glitch_length/data_delta):]))]
			glitch_amp_W    = data_W[int(.25*glitch_length/data_delta):][np.argmax(np.abs(data_W[int(.25*glitch_length/data_delta):]))]
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
			quick_plot(*streamZNE, *stream, data_labels=(streamZNE[0].id, streamZNE[1].id, streamZNE[2].id, 'Rectinilinearity'), verts=glitch_starts, title='%s glitches' % len(glitch_starts), xlabel='Time', ylabel='Normalized Amplitude', show=False, outfile=output_plot)
def glitch_detector2(RAW_UVW, stepwin_length_total=80, ACCfilter=3, threshold=1e-7, glitch_mindist=5, glitch_length=25, plot=True):
	

	"""
	https://stackoverflow.com/questions/48000663/step-detection-in-one-dimensional-data
	"""

	now = time.time()



	### PREPARE DATA, YA' KNOW, RAW & DIS & VEL & ACC as UVW & ZNE ..
	stRAW = read(RAW_UVW)
	stRAW.sort(reverse=False)	# if you change this, god may stand by your side
	stRAW.trim_common() 		# all start with very same time (not necessary, but still desirable)
	stRAW.truncate(0.4,0.7)
	stRAW._set_inventory()		# get & set station inventory for stream (response information)

	stDIS = stRAW.copy()
	stDIS.remove_response(inventory=stRAW.inventory, output='DISP', pre_filt=None, water_level=60)
	stZNE = stDIS.copy()
	stZNE.rotate('->ZNE', inventory=stZNE.inventory, components=('UVW'))
	stZNE.sort(reverse=True)
	stDIS = stDIS + stZNE  		# having in DIS: UVWZNE (in that order)

	stVEL = stRAW.copy()
	stVEL.remove_response(inventory=stRAW.inventory, output='VEL',  pre_filt=None, water_level=60)
	stZNE = stVEL.copy()
	stZNE.rotate('->ZNE', inventory=stZNE.inventory, components=('UVW'))
	stZNE.sort(reverse=True)
	stVEL = stVEL + stZNE  		# having in VEL: UVWZNE (in that order)

	stACC = stRAW.copy()
	stACC.remove_response(inventory=stRAW.inventory, output='ACC',  pre_filt=None, water_level=60)
	stACC_fil = stACC.copy()
	stACC_fil.filtering(ACCfilter)	# LP <1 Hz to better see step in acceleration
	stACC_fil.filtering(4)	# LP <1 Hz to better see step in acceleration
	stZNE = stACC.copy()
	stZNE.rotate('->ZNE', inventory=stZNE.inventory, components=('UVW'))
	stZNE.sort(reverse=True)
	stACC = stACC + stZNE  		# having in ACC: UVWZNE (in that order)

	stCON = stRAW.copy()  		# dummy container having to put convolution into 

	stZNE = stRAW.copy()
	stZNE.rotate('->ZNE', inventory=stZNE.inventory, components=('UVW'))
	stZNE.sort(reverse=True)
	stRAW = stRAW + stZNE
	stRAW.origin = stZNE.origin  # having in RAW: UVWZNE (in that order)



	# DOES IT MAKES SENSE --> ZNE in RAW? Should rather remove gain?
	st = Stream2(stRAW+stDIS+stVEL+stACC+stACC_fil+stCON)
	st.taper(5) 				# in per cent

	output_glitch        = os.path.join( os.path.dirname(stRAW.origin), 'glitches2_' + '.'.join(os.path.basename(stRAW.origin).split('.')[:-1]) + '.txt' )
	data_start, data_end = stRAW.times
	data_delta           = stRAW[0].stats.delta	



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



	### UNIFY GLITCH START TIMES ACROSS UVW-COMPONENTS
	glitch_start_times = np.array(glitch_start_times)
	flat_starts        = np.sort( np.hstack(glitch_start_times) )		# np.flatten wouldn't work because sublists contain different of elements
	
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
		except KeyError: 	# no yet in dictionary
			glitch_starts[str(unified_glitch_start)]  = unified_glitch_comps

		#print(k, shifter_avoid)
		#print(same_glitch_starts)
		#print(unified_glitch_start)
		#print(unified_glitch_comps)
		#print()
		shifter_avoid += len(same_glitch_starts)-1

	glitch_starts     = np.array( [[time, ''.join(sorted(set(glitch_starts[time]))) ] for time in glitch_starts.keys()] )
	glitch_starts     = glitch_starts[np.argsort(glitch_starts[:,0])]
	glitch_starts_all = glitch_starts[:,0]
	glitch_starts_U   = [time_comp[0] for time_comp in glitch_starts if 'U' in time_comp[1]]
	glitch_starts_V   = [time_comp[0] for time_comp in glitch_starts if 'V' in time_comp[1]]
	glitch_starts_W   = [time_comp[0] for time_comp in glitch_starts if 'W' in time_comp[1]]



	### SHELLING GLITCHES & CALCULATE THEIR PARAMETERS / AMPLITUDES
	header         = u'   0                    1                    2         3         4         5          6          7          8          9         10         11         12         13         14         15         16         17         18         19         20         21         22         23         24         25         26         27         28         29         30         31         32         33         34         35         36         37         38         39         40         41         42        43          44         45         46         47         48         49         50         51         52         53       54       55       56       57       58       59       60       61       62       63       64       65\n' \
	                 u' NUM         GLITCH-START           GLITCH-END  U-GLITCH  V-GLITCH  W-GLITCH      U-RAW      V-RAW      W-RAW      Z-RAW      N-RAW      E-RAW      U-DIS      V-DIS      W-DIS      Z-DIS      N-DIS      E-DIS      U-VEL      V-VEL      W-VEL      Z-VEL      N-VEL      E-VEL      U-ACC      V-ACC      W-ACC      Z-ACC      N-ACC      E-ACC  SNR_U-RAW  SNR_V-RAW  SNR_W-RAW  SNR_Z-RAW  SNR_N-RAW  SNR_E-RAW  SNR_U-DIS  SNR_V-DIS  SNR_W-DIS  SNR_Z-DIS  SNR_N-DIS  SNR_E-DIS  SNR_U-VEL  SNR_V-VEL  SNR_W-VEL  SNR_Z-VEL  SNR_N-VEL  SNR_E-VEL  SNR_U-ACC  SNR_V-ACC  SNR_W-ACC  SNR_Z-ACC  SNR_N-ACC  SNR_E-ACC  BAZ-RAW  INC-RAW  POL-RAW  BAZ-DIS  INC-DIS  POL-DIS  BAZ-VEL  INC-VEL  POL-VEL  BAZ-ACC  INC-ACC  POL-ACC\n' \
	                 u'----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'
	line_formatter = u'%05d:  %15s  %15s  %8d  %8d  %8d  %9d  %9d  %9d  %9d  %9d  %9d  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %7.1f  %7.1f  %7.3f  %7.1f  %7.1f  %7.3f  %7.1f  %7.1f  %7.3f  %7.1f  %7.1f  %7.3f'

	print(u'INPUT STREAM:')
	print(stRAW)
	#print()
	#print(u'# '+header)

	glitches       = []
	glitch_counter = 0
	for glitch_start, comps in glitch_starts:

		# prepare data outputs
		glitch_start    = UTCDateTime(glitch_start) 
		glitch_end      = glitch_start+glitch_length
		glitch_counter += 1
	
		# prepare different streams
		stRAW_window = stRAW.copy()
		stDIS_window = stDIS.copy()
		stVEL_window = stVEL.copy()
		stACC_window = stACC.copy()
		
		stRAW_window.trim(starttime=glitch_start-60, endtime=glitch_start+glitch_length+60)
		stDIS_window.trim(starttime=glitch_start-60, endtime=glitch_start+glitch_length+60)
		stVEL_window.trim(starttime=glitch_start-60, endtime=glitch_start+glitch_length+60)
		stACC_window.trim(starttime=glitch_start-60, endtime=glitch_start+glitch_length+60)
		 
		stRAW_window.detrend('demean')
		stDIS_window.detrend('demean')
		stVEL_window.detrend('demean')
		stACC_window.detrend('demean')

		# calculate SNRs on demeaned window (SNR=0 means the mean of the data window is zero!)
		SNR_U_RAW = snr(stDIS_window[0].data, ddof=1)
		SNR_V_RAW = snr(stDIS_window[1].data, ddof=1)
		SNR_W_RAW = snr(stDIS_window[2].data, ddof=1)
		SNR_Z_RAW = snr(stDIS_window[3].data, ddof=1)
		SNR_N_RAW = snr(stDIS_window[4].data, ddof=1)
		SNR_E_RAW = snr(stDIS_window[5].data, ddof=1)
		SNR_U_DIS = snr(stDIS_window[0].data, ddof=1)
		SNR_V_DIS = snr(stDIS_window[1].data, ddof=1)
		SNR_W_DIS = snr(stDIS_window[2].data, ddof=1)
		SNR_Z_DIS = snr(stDIS_window[3].data, ddof=1)
		SNR_N_DIS = snr(stDIS_window[4].data, ddof=1)
		SNR_E_DIS = snr(stDIS_window[5].data, ddof=1)
		SNR_U_VEL = snr(stVEL_window[0].data, ddof=1)
		SNR_V_VEL = snr(stVEL_window[1].data, ddof=1)
		SNR_W_VEL = snr(stVEL_window[2].data, ddof=1)
		SNR_Z_VEL = snr(stVEL_window[3].data, ddof=1)
		SNR_N_VEL = snr(stVEL_window[4].data, ddof=1)
		SNR_E_VEL = snr(stVEL_window[5].data, ddof=1)
		SNR_U_ACC = snr(stACC_window[0].data, ddof=1)
		SNR_V_ACC = snr(stACC_window[1].data, ddof=1)
		SNR_W_ACC = snr(stACC_window[2].data, ddof=1)
		SNR_Z_ACC = snr(stACC_window[3].data, ddof=1)
		SNR_N_ACC = snr(stACC_window[4].data, ddof=1)
		SNR_E_ACC = snr(stACC_window[5].data, ddof=1)

		# trim streams further to glitch window for amplitude extractions
		stRAW_window.trim(starttime=glitch_start, endtime=glitch_start+glitch_length)
		stDIS_window.trim(starttime=glitch_start, endtime=glitch_start+glitch_length)
		stVEL_window.trim(starttime=glitch_start, endtime=glitch_start+glitch_length)
		stACC_window.trim(starttime=glitch_start, endtime=glitch_start+glitch_length)

		data_U_RAW = stRAW_window[0].data
		data_V_RAW = stRAW_window[1].data
		data_W_RAW = stRAW_window[2].data
		data_Z_RAW = stRAW_window[3].data
		data_N_RAW = stRAW_window[4].data
		data_E_RAW = stRAW_window[5].data
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
		versatz      = int(.1*glitch_length/data_delta)	# in samples, to avoid FIR-precursors 
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
		U_RAW = data_U_RAW[versatz:][np.argmax(np.abs(data_U_RAW[versatz:]))]
		U_DIS = data_U_DIS[versatz:][np.argmax(np.abs(data_U_DIS[versatz:]))]
		U_VEL = data_U_VEL[versatz:][np.argmax(np.abs(data_U_VEL[versatz:]))]
		U_ACC = data_U_ACC[versatz:][np.argmax(np.abs(data_U_ACC[versatz:]))]
		
		V_RAW = data_V_RAW[versatz:][np.argmax(np.abs(data_V_RAW[versatz:]))]
		V_DIS = data_V_DIS[versatz:][np.argmax(np.abs(data_V_DIS[versatz:]))]
		V_VEL = data_V_VEL[versatz:][np.argmax(np.abs(data_V_VEL[versatz:]))]
		V_ACC = data_V_ACC[versatz:][np.argmax(np.abs(data_V_ACC[versatz:]))]
		
		W_RAW = data_W_RAW[versatz:][np.argmax(np.abs(data_W_RAW[versatz:]))]
		W_DIS = data_W_DIS[versatz:][np.argmax(np.abs(data_W_DIS[versatz:]))]
		W_VEL = data_W_VEL[versatz:][np.argmax(np.abs(data_W_VEL[versatz:]))]
		W_ACC = data_W_ACC[versatz:][np.argmax(np.abs(data_W_ACC[versatz:]))]
		
		Z_RAW = data_Z_RAW[versatz:][np.argmax(np.abs(data_Z_RAW[versatz:]))]
		Z_DIS = data_Z_DIS[versatz:][np.argmax(np.abs(data_Z_DIS[versatz:]))]
		Z_VEL = data_Z_VEL[versatz:][np.argmax(np.abs(data_Z_VEL[versatz:]))]
		Z_ACC = data_Z_ACC[versatz:][np.argmax(np.abs(data_Z_ACC[versatz:]))]
		
		N_RAW = data_N_RAW[versatz:][np.argmax(np.abs(data_N_RAW[versatz:]))]
		N_DIS = data_N_DIS[versatz:][np.argmax(np.abs(data_N_DIS[versatz:]))]
		N_VEL = data_N_VEL[versatz:][np.argmax(np.abs(data_N_VEL[versatz:]))]
		N_ACC = data_N_ACC[versatz:][np.argmax(np.abs(data_N_ACC[versatz:]))]
		
		E_RAW = data_E_RAW[versatz:][np.argmax(np.abs(data_E_RAW[versatz:]))]
		E_DIS = data_E_DIS[versatz:][np.argmax(np.abs(data_E_DIS[versatz:]))]
		E_VEL = data_E_VEL[versatz:][np.argmax(np.abs(data_E_VEL[versatz:]))]
		E_ACC = data_E_ACC[versatz:][np.argmax(np.abs(data_E_ACC[versatz:]))]


		# get glitch BAZ, INC, and POL
		phi_3D_RAW, INC_3D_RAW, Rec_3D_RAW = np.array( ppol_calc(data_N_RAW, data_E_RAW, data_Z_RAW, fix_angles='AMP', Xoffset_in_samples_for_amplitude=versatz) )[[1,6,10]]
		phi_3D_DIS, INC_3D_DIS, Rec_3D_DIS = np.array( ppol_calc(data_N_DIS, data_E_DIS, data_Z_DIS, fix_angles='AMP', Xoffset_in_samples_for_amplitude=versatz) )[[1,6,10]]
		phi_3D_VEL, INC_3D_VEL, Rec_3D_VEL = np.array( ppol_calc(data_N_VEL, data_E_VEL, data_Z_VEL, fix_angles='AMP', Xoffset_in_samples_for_amplitude=versatz) )[[1,6,10]]
		phi_3D_ACC, INC_3D_ACC, Rec_3D_ACC = np.array( ppol_calc(data_N_ACC, data_E_ACC, data_Z_ACC, fix_angles='AMP', Xoffset_in_samples_for_amplitude=versatz) )[[1,6,10]]


		# assing all to glitch
		glitch = glitch_counter, \
				 glitch_start.strftime('%Y-%m-%dT%H:%M:%S'), \
				 glitch_end.strftime('%Y-%m-%dT%H:%M:%S'),\
				 U_GLITCH, \
				 V_GLITCH, \
				 W_GLITCH, \
				 U_RAW, \
				 V_RAW, \
				 W_RAW, \
				 Z_RAW, \
				 N_RAW, \
				 E_RAW, \
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
				 SNR_U_RAW, \
				 SNR_V_RAW, \
				 SNR_W_RAW, \
				 SNR_Z_RAW, \
				 SNR_N_RAW, \
				 SNR_E_RAW, \
				 SNR_U_DIS, \
				 SNR_V_DIS, \
				 SNR_W_DIS, \
				 SNR_Z_DIS, \
				 SNR_N_DIS, \
				 SNR_E_DIS, \
				 SNR_U_VEL, \
				 SNR_V_VEL, \
				 SNR_W_VEL, \
				 SNR_Z_VEL, \
				 SNR_N_VEL, \
				 SNR_E_VEL, \
				 SNR_U_ACC, \
				 SNR_V_ACC, \
				 SNR_W_ACC, \
				 SNR_Z_ACC, \
				 SNR_N_ACC, \
				 SNR_E_ACC, \
				 phi_3D_RAW, \
				 INC_3D_RAW, \
				 Rec_3D_RAW, \
				 phi_3D_DIS, \
				 INC_3D_DIS, \
				 Rec_3D_DIS, \
				 phi_3D_VEL, \
				 INC_3D_VEL, \
				 Rec_3D_VEL, \
				 phi_3D_ACC, \
				 INC_3D_ACC, \
				 Rec_3D_ACC


		# prapare file output and print into terminal
		glitches.append(line_formatter % glitch)
		#print(line_formatter % glitch)



	### FILING GLITCHES
	np.savetxt(output_glitch, glitches, fmt='%s', header=header)



	### SHELLING SUMMARY
	print()
	print(u'- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')
	print(u'Data start:          %s'               % data_start)
	print(u'Data end:            %s'               % data_end)
	print(u'Data length:         %s (h:m:s)'       % sec2hms( data_end-data_start ))
	print(u'Data SPS:            %sHz (%ss)'       % (1./data_delta, data_delta))
	print()
	print(u'Step window length:  %s samples'       % stepwin_length_total)
	print(u'Threshold:           %s'               % threshold)
	print(u"Glitches' min dist:  %ss"              % glitch_mindist)
	print(u"Glitch length (fix): %ss"              % glitch_length)
	print()
	print(u'Glitches ALL:        %s'               % glitch_counter)
	print(u'        on U:        %s'               % len(glitch_starts_U))
	print(u'        on V:        %s'               % len(glitch_starts_V))
	print(u'        on W:        %s'               % len(glitch_starts_W))
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

		quick_plot(stCON[2], stACC_fil[2], stRAW[2], stDIS[2], stVEL[2], stACC[2], data_labels=('W CON','W ACC-FIL','W RAW','W DIS','W VEL','W ACC'), title='%s glitches' % glitch_counter, verts=glitch_starts_W, xlabel='Time', ylabel='Amplitude')
		#quick_plot(stCON[2], stACC_fil[2], stACC_fil2[2],  stRAW[2], data_labels=('W CON', 'W ACC FIL', 'W ACC FIL2', 'W RAW'), title='%s glitches' % glitch_counter, verts=glitch_starts_W, xlabel='Time', ylabel='Amplitude')
		#quick_plot(stACC_fil[0], stACC_fil2[0],  stRAW[0], data_labels=('U ACC FIL', 'U ACC FIL2', 'U RAW'), title='%s glitches' % glitch_counter, verts=glitch_starts_U, xlabel='Time', ylabel='Amplitude')
		#quick_plot(*stCON, *stACC_fil, *stRAW, data_labels=('U CON', 'V CON', 'W CON', 'U ACC FIL', 'V ACC FIL', 'W ACC FIL', 'U RAW', 'V RAW', 'W RAW'), title='%s glitches' % glitch_counter, verts=glitch_starts_all, xlabel='Time', ylabel='Amplitude' )


### plot scripts to visualise detected glitches



### _ _ N A M E _ _ = = " _ _ M A I N _ _ "  
if __name__ == "__main__":	
	print('Define Testing')