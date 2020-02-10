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
import scipy
import datetime
import numpy as np


#####  obspy modules import  #####
from obspy import UTCDateTime, read


#####  seisglitch modules import  #####
from seisglitch.math import moving_window, normalise
from seisglitch.util import read2, marstime, sec2hms, quick_plot, ppol

#TEST IF TAPERING TIMES CAN BE USED ANYWAY
#CHECK WHY SOME FILES DON'T WORK


### GLITCH DETECTION
def glitch_detector(*RAW_UVW, 
    inventory_file            = 'IRIS',
    taper_length_per_side     = 0.05,
    pre_filt                  = None,
    water_level               = 60,
    ACCfilter                 = {'type' : 'bandpass',  'options' : {'freqmin':0.001, 'freqmax':0.1}, 'string':'bp: 0.001 < f < 0.1'}, 
    window_length_minutes     = 5, 
    average_peak_height_times = 3, 
    plot_triggering           = False,
    glitch_min_dist           = 5, 
    glitch_length             = 20,
    glitch_min_polarization   = 0.9, 
    plot_individual           = False,
    plot_all                  = False):
    

    """
    Hard-coded FIR pre-cursor offset of 1 s.
    """



    ### OUTPUT
    print()
    print(u'  -----------------------')
    print(u'  RUNNING GLITCH DETECTOR')
    print(u'  -----------------------')



    ### SOME VARIABLES
    header         = u'    0                    1                    2              3              4         5         6         7          8          9         10           11         12         13         14         15         16           17         18         19         20         21         22           23         24         25         26         27         28           29         30         31         32         33         34            35          36          37          38            39          40          41          42            43          44          45          46            47          48          49          50\n' \
                     u'  NUM         GLITCH-UTC-S         GLITCH-UTC-E  GLITCH-LMST-S  GLITCH-LMST-E  U-GLITCH  V-GLITCH  W-GLITCH      U-RAW      V-RAW      W-RAW        U-GAI      V-GAI      W-GAI      Z-GAI      N-GAI      E-GAI        U-DIS      V-DIS      W-DIS      Z-DIS      N-DIS      E-DIS        U-VEL      V-VEL      W-VEL      Z-VEL      N-VEL      E-VEL        U-ACC      V-ACC      W-ACC      Z-ACC      N-ACC      E-ACC    PHI_3D_GAI  INC_3D_GAI  SNR_3D_GAI  POL_3D_GAI    PHI_3D_DIS  INC_3D_DIS  SNR_3D_DIS  POL_3D_DIS    PHI_3D_VEL  INC_3D_VEL  SNR_3D_VEL  POL_3D_VEL    PHI_3D_ACC  INC_3D_ACC  SNR_3D_ACC  POL_3D_ACC\n' \
                     u'------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'
    line_formatter = u'%06d:  %19s  %19s  %13s  %13s  %8d  %8d  %8d  %9d  %9d  %9d    %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g    %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g    %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g    %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g    %10.1f  %10.1f  %10.1g  %10.3f    %10.1f  %10.1f  %10.1g  %10.3f    %10.1f  %10.1f  %10.1g  %10.3f    %10.1f  %10.1f  %10.1g  %10.3f'
        


    ### DOING THE SAME ANALYSIS / OUTPUT FOR ALL PASSED FILES
    for o in range(len(RAW_UVW)):
    

        ### READ DATA
        now    = time.time()        
        stream = read2(RAW_UVW[o])
        stream.sort(reverse=False)                 # Sorted: U then V then W: if you change this, god may stand by your side
        #stream.trim2(0.5,0.55)        
        stream._set_inventory(source=inventory_file)  # get & set station inventory for stream (response information)
        data_delta           = stream[0].stats.delta
        data_start, data_end = [marstime(time) for time in stream.times]


        ### OUTPUT
        print()
        print(u'Analysing file: %s/%s' % (o+1,len(RAW_UVW)))
        print()
        print(u'INPUT STREAM:')
        print( str(stream).replace('\n','\n  '))


        ### QUICK SANITY CHECK
        num_ids = len(stream._get_ids())
        if num_ids != 3:
            print()
            print(u'ERROR: Found %s instead of 3 trace ids in file: %s' % (num_ids, ', '.join(stream._get_ids())))
            print(u'       Cannot perform analysis.')
            continue

        if data_end.UTC_time-data_start.UTC_time<60:
            print()
            print(u'WARNING: File length <60 s. Too little for proper analysis. Skipped.')
            continue


        ### VARIABLES
        output_txt             = os.path.join( os.path.dirname(RAW_UVW[o]), 'glitches_' + '.'.join(os.path.basename(RAW_UVW[o]).split('.')[:-1]) + '.txt' )
        glitch_dict            = {}
        glitch_time_comps_keep = []  
        glitches_text          = []
        glitch_counter         = 0
        print_counter          = 0


        ### LOOP OVER SPLIT STREAM OBJECT (split each file into 12h parts - avoiding timing & RAM issues)
        length = 3600 * 12
        step   = length*(1-2*taper_length_per_side) - 3600*0.5
        for stream_slide in stream.slide(length, step, include_partial_windows=True):


            ## PREPARE DATA: RAW, DIS, VEL, and ACC for both UVW & ZNE ..
            print()
            print(u'Info: Prepare data (tapering, remove instrument response, filtering, ..)')

            stRAW = stream_slide.copy()
            stRAW.detrend('demean')        
            stRAW.taper(taper_length_per_side)       
            stRAW.detrend('simple')        
            
            stGAI = stRAW.copy()
            stGAI.gain_correction()
            stZNE = stGAI.copy()
            stZNE.rotate('->ZNE', inventory=stream.inventory, components=('UVW'))
            stZNE.sort(reverse=True)
            stGAI += stZNE                                            # having in GAI (overall sensitivity and gain corrected): UVWZNE (in that order)
            
            stDIS = stRAW.copy()
            stDIS.remove_response(inventory=stream.inventory, output='DISP', pre_filt=eval(str(pre_filt)), water_level=eval(str(water_level)))
            stZNE = stDIS.copy()
            stZNE.rotate('->ZNE', inventory=stream.inventory, components=('UVW'))
            stZNE.sort(reverse=True)
            stDIS += stZNE                                            # having in DIS: UVWZNE (in that order)
            
            stVEL = stRAW.copy()
            stVEL.remove_response(inventory=stream.inventory, output='VEL',  pre_filt=eval(str(pre_filt)), water_level=eval(str(water_level)))
            stZNE = stVEL.copy()
            stZNE.rotate('->ZNE', inventory=stream.inventory, components=('UVW'))
            stZNE.sort(reverse=True)
            stVEL += stZNE                                            # having in VEL: UVWZNE (in that order)

            stACC = stRAW.copy()
            stACC.remove_response(inventory=stream.inventory, output='ACC',  pre_filt=eval(str(pre_filt)), water_level=eval(str(water_level)))
            stACC_fil = stACC.copy()
            stACC_fil.detrend('demean')        
            stACC_fil.taper(taper_length_per_side)                    # ease response effects at edges of data
            stACC_fil.detrend('simple')        
            stACC_fil.filter(ACCfilter['type'], **ACCfilter['options'])  
            stACC_fil.trim2(taper_length_per_side, 1-taper_length_per_side)
 
            stACC_fil_der = stACC_fil.copy()
            for trace in stACC_fil_der:
                trace.data = np.gradient(trace.data)

            stZNE = stACC.copy()
            stZNE.rotate('->ZNE', inventory=stream.inventory, components=('UVW'))
            stZNE.sort(reverse=True)
            stACC += stZNE                                             # having in ACC: UVWZNE (in that order)


            ## CALCULATE THRESHOLD FOR GLITCH TRIGGERING FOR EACH COMPONENT
            print(u'Info: Calculate moving average-peak-height of filtered acceleration data.')

            stAPH  = stACC_fil_der.copy()         # dummy container to put AVERAGE PEAK HEIGHT calculations into 
            for h in range(len( stAPH )):

                average_peak_heights      = []
                window_length_in_samples  = int(window_length_minutes*60*stAPH[h].stats.sampling_rate)    # length of moving window to retrieve `average_peak_height` from 
                step_in_samples           = window_length_in_samples//4                               # quarter of `window_length_minutes`

                for window_start, window_end, window_data in moving_window(stAPH[h].data, window_length_in_samples=window_length_in_samples, step_in_samples=step_in_samples, equal_end=True):
                    peaks_window, _       = scipy.signal.find_peaks(np.abs( window_data ))
                    average_peak_height   = np.average( np.abs( window_data[peaks_window] ))
                    average_peak_heights.append(average_peak_height)

                average_peak_heights      = np.repeat(average_peak_heights, step_in_samples)           # now number of energies equal to "len(data)-step_in_samples"
                
                shift                     = int(window_length_in_samples//2)
                try:
                    average_peak_heights  = np.hstack(( np.ones(shift)*average_peak_heights[0], average_peak_heights, np.ones(shift)*average_peak_heights[-1])) # pre- and append due to shift of window so its corresponding value is middle of window
                except IndexError:
                    print(u'No peaks detected on component %s. Moving on.' % stAPH[h].stats.channel)
                    print()
                    continue
                stAPH[h].data = average_peak_heights[:len(stAPH[h].data)] * average_peak_height_times  # these are the thresholds used to detect glitches
            

            ## FIND PEAKS FOR EACH OF UVW-COMPONENTS = GLITCH-STARTS FOR EACH COMPONENT INDIVIDUALLY !
            print(u'Info: Find peaks in time derivative of filtered acceleration, using `average_peak_height_times`.')

            glitch_start_times = [[] for i in range(len(stRAW))]
            for l in range(len( stACC_fil_der )):
        
                peaks, _              = scipy.signal.find_peaks(np.abs(stACC_fil_der[l].data), height=stAPH[l].data, distance=glitch_min_dist*stACC_fil_der[l].stats.sampling_rate)
                trigger_times         = stACC_fil_der[l].times('utcdatetime')[peaks]
                glitch_start_times[l] = trigger_times
            

            ### SHOW TRIGGERING, if wished
            if plot_triggering:

                print(u'Info: Triggering on derivative of filtered acceleration data is only done on parts not affected by taper.')
                print(u'Info: Amplitudes are scaled so plots are readable, they do not reflect true amplitude values.')
                print(u'Info: Glitches will still be unified across components according to `glitch_min_dist`.')
                print(u'Info: Glitches will still be excluded based on their polarization, according to `glitch_min_polarization`.')

                for p, comp in enumerate('UVW'):

                    ACC_fil_der = stACC_fil_der.copy()
                    RAW         = stRAW.copy()
                    ACC_fil     = stACC_fil.copy()
                    APH1        = stAPH.copy()
                    APH2        = stAPH.copy()
                    
                    mini        = min(stACC_fil_der[p].data)
                    maxi        = max(stACC_fil_der[p].data) 

                    ACC_fil_der.normalise( [mini              , maxi              ] )
                    RAW.filter('highpass', **{'freq':0.001, 'corners':3, 'zerophase':True})
                    RAW.normalise(         [mini-0.5*(maxi-mini), maxi-0.5*(maxi-mini)] )
                    ACC_fil.normalise(     [mini-1.0*(maxi-mini), maxi-1.0*(maxi-mini)] )
                    APH2[p].data *= -1     
                
                    title       = '%s-component, %s glitches' % (comp, len(glitch_start_times[p]))
                    data_labels = ('Derivative of filtered ACC', 'RAW original (HP > 0.001Hz)', 'filtered ACC', ' %s * average-peak-height' % average_peak_height_times, '-%s * average-peak-height' % average_peak_height_times)
                    verticals   = glitch_start_times[p]
                    if p!=2:
                        quick_plot(stACC_fil_der[p], RAW[p], ACC_fil[p], APH1[p], APH2[p], win_title=title, title=title, data_labels=data_labels, verts=[verticals], show=False, keep_open=True)
                    else:
                        quick_plot(stACC_fil_der[p], RAW[p], ACC_fil[p], APH1[p], APH2[p], win_title=title, title=title, data_labels=data_labels, verts=[verticals], show=True)


            ## UNIFY GLITCH START TIMES ACROSS UVW-COMPONENTS (using parameter: `glitch_min_dist` given in seconds)
            print(u'Info: Unifying glitches across components, using `glitch_min_dist`.')
            print()

            if all(not list(liste) for liste  in glitch_start_times):       # no possible glitches detected
                continue

            glitch_start_times = np.array(glitch_start_times)
            flat_starts        = np.sort( np.hstack(glitch_start_times) )   # np.flatten wouldn't work because sublists contain different number of elements ..
            
            same_glitch_starts = []
            shifter_avoid      = 0

            for l in range( len( flat_starts )):

                # unify glitch start times: loop glitches, find all +/- `glitch_min_dist`, take minimum, go to next not unified glitch (`shifter_avoid`)
                k = l + shifter_avoid
                try:
                    indices_same_glitch = np.where( (flat_starts < flat_starts[k]+glitch_min_dist) & (flat_starts > flat_starts[k]-glitch_min_dist))
                except IndexError:      # composed index `k` >= len(flat_starts) --> end of data basically
                    break

                same_glitch_starts   = [flat_starts[i] for i in indices_same_glitch[0] if flat_starts[i] not in same_glitch_starts]
                shifter_avoid       += len(same_glitch_starts)-1
                unified_glitch_start = min( same_glitch_starts )

                # from `same_glitch_starts`, collect respective components
                unified_glitch_comps = []
                for gs in same_glitch_starts:
                    for m, comp_starts in enumerate(glitch_start_times):
                        if gs in comp_starts and not str(m) in unified_glitch_comps:
                            unified_glitch_comps.append( str(m) )
                            break

                unified_glitch_comps = ''.join( sorted( unified_glitch_comps ) )
                unified_glitch_comps = unified_glitch_comps.replace('0','U').replace('1','V').replace('2','W')

                # assign `unified_glitch_comps` to keys `unified_glitch_start`
                try:
                    glitch_dict[str(unified_glitch_start)] += unified_glitch_comps
                except KeyError:    # not yet in dictionary
                    glitch_dict[str(unified_glitch_start)]  = unified_glitch_comps

            # from dictionary containing glitches, create sorted list
            glitch_time_comps = np.array( [[UTCDateTime(time), ''.join(sorted(set(glitch_dict[time]))) ] for time in glitch_dict.keys()] )
            glitch_time_comps = glitch_time_comps[np.argsort(glitch_time_comps[:,0])]         # time sorted, chronologic


            ## GET / CALCULATE GLITCH PARAMETERS (e.g. AMPLITUDES, BACKAZIMUTHS, INCIDENCE ANGLES, SNRs, POLARIZATIONS)
            for glitch_UTC_start, comps in glitch_time_comps:

                # avoiding re-iteratring over glitches in `glitch_time_comps` already processed im previous `stream_slide` window
                if glitch_UTC_start<stRAW[0].stats.starttime:
                    continue

                # prepare data outputs
                glitch_time_start = marstime(glitch_UTC_start)
                glitch_time_end   = marstime(glitch_UTC_start+glitch_length)

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

                # prepare different streams
                stRAW_window = stRAW.copy()
                stGAI_window = stGAI.copy()
                stDIS_window = stDIS.copy()
                stVEL_window = stVEL.copy()
                stACC_window = stACC.copy()

                # trim streams further to glitch window for amplitude extractions
                stRAW_window.trim(starttime=glitch_UTC_start, endtime=glitch_UTC_start+glitch_length)
                stGAI_window.trim(starttime=glitch_UTC_start, endtime=glitch_UTC_start+glitch_length)
                stDIS_window.trim(starttime=glitch_UTC_start, endtime=glitch_UTC_start+glitch_length)
                stVEL_window.trim(starttime=glitch_UTC_start, endtime=glitch_UTC_start+glitch_length)
                stACC_window.trim(starttime=glitch_UTC_start, endtime=glitch_UTC_start+glitch_length)
                
                # demean
                stRAW_window.detrend('demean')
                stGAI_window.detrend('demean')
                stDIS_window.detrend('demean')
                stVEL_window.detrend('demean')
                stACC_window.detrend('demean')

                # assign data for amplitude determination
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

                # get UVWZNE amplitudes of glitches
                offset = int(1./data_delta) # 1 second (in samples) offset to avoid FIR-precursors 
                
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

                # get glitch BAZ, INC, SNR, and POL (=rectinilinearity)
                GAI_measurement = ppol(stream=stGAI_window, demean=False, fix_angles='AMP', Xoffset_samples_for_amplitude=offset)
                DIS_measurement = ppol(stream=stDIS_window, demean=False, fix_angles='AMP', Xoffset_samples_for_amplitude=offset)
                VEL_measurement = ppol(stream=stVEL_window, demean=False, fix_angles='AMP', Xoffset_samples_for_amplitude=offset)
                ACC_measurement = ppol(stream=stACC_window, demean=False, fix_angles='AMP', Xoffset_samples_for_amplitude=offset)

                PHI_3D_GAI, INC_3D_GAI, SNR_3D_GAI, POL_3D_GAI = np.asarray( GAI_measurement.results )[[2,6,9,13]]
                PHI_3D_DIS, INC_3D_DIS, SNR_3D_DIS, POL_3D_DIS = np.asarray( DIS_measurement.results )[[2,6,9,13]]
                PHI_3D_VEL, INC_3D_VEL, SNR_3D_VEL, POL_3D_VEL = np.asarray( VEL_measurement.results )[[2,6,9,13]]
                PHI_3D_ACC, INC_3D_ACC, SNR_3D_ACC, POL_3D_ACC = np.asarray( ACC_measurement.results )[[2,6,9,13]]
                

                # exclude glitches with polarization < `glitch_min_polarization` 
                # (glitches have high polarization, thus this condition intends to throw out non-glitches)
                if POL_3D_GAI < float( glitch_min_polarization ):
                    continue
               
                # exclude glitches that already were detected 
                # (this is the case when a glitch was detected in 2 subsequent 12 h data windows that overlap a 1 hour + 2 * taper length)
                try:
                    condition = glitch_UTC_start >= glitch_time_comps_keep[-1][0] + glitch_min_dist
                    if not condition:
                        continue
                except IndexError:              # happens only if `glitch_time_comps_keep` is still empty
                    pass
                glitch_time_comps_keep.append( [glitch_UTC_start, comps] )

                # collect all glitch measures
                glitch_counter += 1                
                glitch = glitch_counter, \
                         glitch_time_start.UTC,  \
                         glitch_time_end.UTC,    \
                         glitch_time_start.LMST, \
                         glitch_time_end.LMST,   \
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

                glitches_text.append(line_formatter % glitch)

                # save individual ppol measurement, if wished
                if plot_individual:
                    GAI_measurement.plot(show=False, outfile='%06d_GAI.png' % glitch_counter)
                    DIS_measurement.plot(show=False, outfile='%06d_DIS.png' % glitch_counter)
                    VEL_measurement.plot(show=False, outfile='%06d_VEL.png' % glitch_counter)
                    ACC_measurement.plot(show=False, outfile='%06d_ACC.png' % glitch_counter)

                # small output to user
                if glitch_counter%5 == 0:

                    if print_counter == 0:
                        print_counter += 1
                        for header_line in header.split('\n')[1:]:
                            print('# ' + header_line[:107])

                    print( (line_formatter % glitch)[:109] )


        ### STATISTICS
        glitch_starts_U        = [time_comp[0] for time_comp in glitch_time_comps_keep if 'U' in time_comp[1]]
        glitch_starts_V        = [time_comp[0] for time_comp in glitch_time_comps_keep if 'V' in time_comp[1]]
        glitch_starts_W        = [time_comp[0] for time_comp in glitch_time_comps_keep if 'W' in time_comp[1]]
        glitch_starts_Uonly    = [time_comp[0] for time_comp in glitch_time_comps_keep if 'U' in time_comp[1] and len(time_comp[1])==1]
        glitch_starts_Vonly    = [time_comp[0] for time_comp in glitch_time_comps_keep if 'V' in time_comp[1] and len(time_comp[1])==1]
        glitch_starts_Wonly    = [time_comp[0] for time_comp in glitch_time_comps_keep if 'W' in time_comp[1] and len(time_comp[1])==1]
        try:
            glitch_indivdual_ratio = (len(glitch_starts_Uonly)+len(glitch_starts_Vonly)+len(glitch_starts_Wonly))/len(glitch_time_comps_keep)*100
        except ZeroDivisionError:
            glitch_indivdual_ratio = np.nan


        ### SUMMARY TEXT
        output    = [u'DATA:',
                     u'  path:                      %s'           % RAW_UVW[o],
                     u'  IDs:                       %s'           % ', '.join(stream._get_ids()),
                     u'  SPS:                       %s Hz (%s s)' % (1./data_delta, data_delta),
                     u'  UTC start:                 %s'           % data_start.UTC_time,
                     u'  UTC end:                   %s'           % data_end.UTC_time,
                     u'  UTC length:                %s (h:m:s)'   % sec2hms( data_end.UTC_time-data_start.UTC_time ),
                     u'  LMST start:                %s'           % data_start.LMST,
                     u'  LMST end:                  %s'           % data_end.LMST,
                     u'  LMST length:               %s (h:m:s)'   % sec2hms( data_end.LMST_time-data_start.LMST_time ),
                     u'',      
                     u'PARAMETERS:',      
                     u'  Taper length per side %%:   %s (%.1f s)' % (taper_length_per_side, (data_end.UTC_time-data_start.UTC_time)*taper_length_per_side),
                     u'  Pre filtering:             %s'           % str(pre_filt),
                     u'  Water level:               %s'           % str(water_level),
                     u'  Acceleration filter:       %s'           % ACCfilter['string'],
                     u'  Window length:             %s min'       % window_length_minutes,
                     u'  Average_peak_height_times: %s'           % average_peak_height_times,
                     u"  Show triggering:           %s"           % str(plot_triggering),
                     u"  Glitches' min dist:        %s s"         % glitch_min_dist,
                     u"  Glitch length (fix):       %s s"         % glitch_length,
                     u'  Glitch min polarization:   %s'           % glitch_min_polarization,
                     u"  Plot individual:           %s "          % str(plot_individual),
                     u'',      
                     u'RESULTS:',      
                     u'    Glitches all:            %s'           % len(glitch_time_comps_keep),
                     u'       all / sol:            %.0f'         % (len(glitch_time_comps_keep)*86400/(data_end.LMST_time-data_start.LMST_time)),
                     u'            on U:            %s'           % len(glitch_starts_U),
                     u'            on V:            %s'           % len(glitch_starts_V),
                     u'            on W:            %s'           % len(glitch_starts_W),
                     u'       only on U:            %s'           % len(glitch_starts_Uonly),
                     u'       only on V:            %s'           % len(glitch_starts_Vonly),
                     u'       only on W:            %s'           % len(glitch_starts_Wonly),
                     u'    indi./ all %%:            %.1f'        % glitch_indivdual_ratio,
                     u'',
                     u'Done in:   %s (h:m:s).'                    % sec2hms( time.time()-now ),
                     u'Timestamp: %s'                             % datetime.datetime.now().strftime('%Y-%m-%dT%H:%M:%S')]

        # create good-looking output for later
        len_max_line  = max( [len(line) for line in output] )
        formatter_str = "\n#   | %-" + str(len_max_line) + 's |'

        output = [formatter_str % line for line in output]
        output.insert(0,           '\n\n#   +-' +  '-' * len_max_line + '-+'   )
        output.insert(len(output),   '\n#   +-' +  '-' * len_max_line + '-+\n' )


        ### FILING GLITCHES + OUTPUT
        np.savetxt(output_txt, glitches_text, fmt='%s', header=header)
        
        with open(output_txt, 'r') as fp_in:
            lines = fp_in.readlines()

        with open(output_txt, 'w') as fp_out:

            for line in lines:
                fp_out.write(line)

            for line in output:
                fp_out.write(line)


        ### FINAL OUTPUT TO SHELL
        print()
        print(u'OUTPUT GLITCH DETECTOR:')
        for line in output:
            print(line.strip('\n'))
        print()
        print(u'Output file:')
        print(output_txt)


        ### FINAL PLOT, if wished
        if plot_all:
            title = '%s glitches' % len(glitch_time_comps_keep)
            verts = np.array( glitch_time_comps_keep )[:,0]
            quick_plot(*stream, title=title, verts=[verts], xlabel='Time', ylabel='Amplitude')



### _ _ N A M E _ _ = = " _ _ M A I N _ _ "  
if __name__ == "__main__": 
    pass