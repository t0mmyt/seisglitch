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
from seisglitch.math import normalise
from seisglitch.util import read2, Stream2, glitch_statistic, pierce_stream, decimate_SEIS, marstime, sec2hms, quick_plot, ppol



### GLITCH DETECTION
def detect(*RAW_UVW, 
    inventory_file          = 'IPGP',
    taper_length_per_side   = 0.05,
    pre_filt                = None,
    water_level             = 60,
    ACCfilter               = {'type' : 'bandpass',  'options' : {'freqmin':0.001, 'freqmax':0.1}, 'string':'0.001 < f < 0.1'},
    threshold               = 1e-10,
    plot_triggering         = False,
    glitch_min_length       = 5, 
    glitch_length           = 20,
    glitch_min_polarization = 0.9, 
    **kwargs):
    

    """
    WRITE DOCSTRING.
    """


    ### OUTPUT
    print()
    print(u'  -----------------------')
    print(u'  RUNNING GLITCH DETECTOR')
    print(u'  -----------------------')



    ### FIXED PARAMETERS
    DATA_MIN_LENGTH = 300      # in s
    OFFSET_AMP      = 1        # in s, avoiding determining glitch amplitude on glitch precursors for majority of cases



    ### SOME VARIABLES
    header         = u'   0                    1                    2              3              4           5         6         7             8          9         10           11         12         13         14         15         16            17              18              19              20              21          22\n' \
                     u' NUM         GLITCH-UTC-S         GLITCH-UTC-E  GLITCH-LMST-S  GLITCH-LMST-E    U-GLITCH  V-GLITCH  W-GLITCH         U-RAW      V-RAW      W-RAW        U-GAI      V-GAI      W-GAI      Z-GAI      N-GAI      E-GAI    AZI_3D_GAI  AZI_ERR_3D_GAI      INC_3D_GAI  INC_ERR_3D_GAI      SNR_3D_GAI  POL_3D_GAI\n' \
                     u'--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'
    line_formatter = u'%06d  %19s  %19s  %13s  %13s    %8d  %8d  %8d     %9d  %9d  %9d    %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g    %10.1f  %14.1f  %14.1f  %14.1f  %14.3g  %10.3f'
        


    ### DOING THE SAME ANALYSIS / OUTPUT FOR ALL PASSED FILES
    for o in range(len(RAW_UVW)):
    

        ### READ DATA
        now           = time.time()
        try: 
            stream = read2(RAW_UVW[o])
        except TypeError:       # could not read file, skip
            continue
        stream_select = stream.select(channel='?[LMH]?')
        stream_ids    = sorted(set( stream_select._get_ids() ))
        stream_comps  = sorted(set( stream_select._get_components() ))
        stream_times  = [marstime(time) for time in stream_select.times]


        ### OUTPUT
        print()
        print()
        print(u'INFO: Analysing seismic traces of file: %s/%s' % (o+1,len(RAW_UVW)))
        print(RAW_UVW[o])
        for trace in stream:
            print('  %s' % trace)


        ### SANITY CHECK
        if len(stream_ids) != 3:
            print()
            print(u'ERROR: Found %s instead of 3 seismic trace ids in file: %s' % (len(stream_ids), ', '.join(stream_ids)))
            print(u'       Cannot perform analysis.')
            continue
        if 'U' not in stream_comps or 'V' not in stream_comps or 'W' not in stream_comps:
            print()
            print(u'ERROR: Found components %s instead of `U`, `V`, and `W`.' % ', '.join(stream_ids))
            print(u'       Cannot perform analysis.')
            continue


        ### VARIABLES
        output_txt           = os.path.join( os.path.dirname(RAW_UVW[o]), 'glitches_' + '.'.join(os.path.basename(RAW_UVW[o]).split('.')[:-1]) + '.txt' )
        glitches_text        = []
        glitches             = []
        total_data_time_UTC  = 0
        total_data_time_LMST = 0
        glitch_counter       = 0
        print_counter        = 0

        for stream2 in stream_select.slide(12*3600, (1-2*taper_length_per_side)*12*3600-60, offset=0, include_partial_windows=True, nearest_sample=True):

            start_read                = min([str(tr.stats.starttime) for tr in stream2]) 
            end_read                  = max([str(tr.stats.endtime)   for tr in stream2])
            stWORK                    = read2(RAW_UVW[o], starttime=UTCDateTime(start_read), endtime=UTCDateTime(end_read))
            streamU, streamV, streamW = pierce_stream(stWORK, minimum_sample_length=10, ids=stream_ids)

            for traceU, traceV, traceW in zip(streamU, streamV, streamW):

                glitch_dict = {}


                ## PRINT INFORMATION
                print()
                print(u'INFO: Preparing data.')
                print('  %s' % traceU)
                print('  %s' % traceV)
                print('  %s' % traceW)


                ## DATA assignment
                stRAW = Stream2(traces=[traceU, traceV, traceW])
                stRAW.sort(reverse=False)                    # Sorted: U then V then W: if you change this, god may stand by your side
                stGAI = stRAW.copy()
                data_start, data_end  = [marstime(time) for time in stRAW.times]


                ## SANITY CHECK
                if data_end.UTC_time-data_start.UTC_time<DATA_MIN_LENGTH:
                    print(u'WARNING: Time length <%s s. Too little for proper analysis. Probably related to data overlaps. Skipped.' % DATA_MIN_LENGTH)
                    continue

                total_data_time_UTC  += data_end.UTC_time  - data_start.UTC_time
                total_data_time_LMST += data_end.LMST_time - data_start.LMST_time


                ## PREPARE RAW DATA
                # downsampling
                for trace in stRAW:
                    if trace.stats.sampling_rate==100:
                        decimate_SEIS(trace, 5, verbose=False)      # SPS then 20
                        decimate_SEIS(trace, 5, verbose=False)      # SPS then 4
                        decimate_SEIS(trace, 2, verbose=False)      # SPS then 2
                    elif trace.stats.sampling_rate==25:
                        print(u'WARNING: Input sampling rate %s. Cannot go to 2 SPS but only 2.5 SPS. Slight differences in detection may be expected.' % trace.stats.sampling_rate)
                        decimate_SEIS(trace, 5, verbose=False)      # SPS then 5
                        decimate_SEIS(trace, 2, verbose=False)      # SPS then 2.5
                    elif trace.stats.sampling_rate==20:
                        decimate_SEIS(trace, 5, verbose=False)      # SPS then 4
                        decimate_SEIS(trace, 2, verbose=False)      # SPS then 2                        
                    elif trace.stats.sampling_rate==10:
                        decimate_SEIS(trace, 5, verbose=False)      # SPS then 2
                    elif trace.stats.sampling_rate==5:
                        print(u'WARNING: Input sampling rate %s. Cannot go to 2 SPS but only 2.5 SPS. Slight differences in detection may be expected.' % trace.stats.sampling_rate)
                        decimate_SEIS(trace, 2, verbose=False)      # SPS then 2.5
                    elif trace.stats.sampling_rate==4:
                        decimate_SEIS(trace, 2, verbose=False)      # SPS then 2
                    elif trace.stats.sampling_rate<2:
                        print(u'WARNING: Input sampling rate %s < 2 SPS. Slight differences in detection may be expected.' % trace.stats.sampling_rate)                
                    else:
                        pass


                ## SANITY CHECK 2
                stRAW_sampling_rates  = sorted(set( stRAW._get_sampling_rates() ))
                if len(stRAW_sampling_rates) != 1:
                    print()
                    print(u'ERROR: Found %s different sampling rates after decimation: %s Hz. They all need to be equal.' % (len(stRAW_sampling_rates), ', Hz'.join(len(stRAW_sampling_rates))))
                    print(u'       Cannot perform analysis.')
                    continue                


                ## PREPARE RAW DATA 2
                stRAW.set_inventory(source=inventory_file)      # get & set station inventory for stream (response information)
                if not stRAW.inventory:                          # Inventory file empty, i.e., for given trace no information available
                    print(u'WARNING: Could not find meta-data for given trace(s). Skipped.')
                    continue

                stGAI.set_inventory(stRAW.inventory)
                if stGAI.gain_correction() == -1:                # Could not gain correct although inventory file given for traces. ==> Inv file maybe faulty
                    log_file = os.path.join(os.getcwd(), 'inventory_error.log')
                    with open(log_file, 'a') as fp:
                        fp.write('Trace %s:  %s - %s\n' % (stGAI[0].id[:-1]+'?', stGAI.times[0], stGAI.times[1]))
                        print(u"WARNING: Skipped. Trace times added to '%s'." % log_file)
                    continue


                ## PREPARE ACCELERATION DATA
                print(u'INFO: Removing instrument response, using `pre_filt` and/or `water_level`.')

                stACC = stRAW.copy()
                stACC.remove_response(inventory=stACC.inventory, output='ACC',  pre_filt=eval(str(pre_filt)), water_level=eval(str(water_level)))
                stACC.detrend('demean')        
                stACC.detrend('simple')        
                stACC.taper(taper_length_per_side)   # ease response / filter effects at edges of data
                if ACCfilter['type']:
                    stACC.filter(ACCfilter['type'], **ACCfilter['options'])  
                stACC.trim2(taper_length_per_side, 1-taper_length_per_side)

                stACC_grad = stACC.copy()
                for trace in stACC_grad:
                    trace.data = np.gradient(trace.data)


                ## FIND PEAKS FOR EACH OF UVW-COMPONENTS = GLITCH-STARTS FOR EACH OF THE COMPONENTS INDIVIDUALLY !
                print(u'INFO: Finding peaks in time derivative of filtered acceleration, using `ACCfilter` and `threshold`.')

                glitch_start_times = [[] for i in range(len(stACC_grad))]
                for l in range(len( stACC_grad )):

                    peaks_all,  _ = scipy.signal.find_peaks(np.abs(stACC_grad[l].data))            
                    peaks_trig, _ = scipy.signal.find_peaks(np.abs(stACC_grad[l].data), height=threshold, distance=glitch_min_length*stRAW_sampling_rates[0])
                    peaks_keep    = []
                    
                    shifter = 0
                    for j in range(len( peaks_all )):

                        try:
                            i    = j + shifter
                            peak = peaks_all[i]

                            if peak in peaks_trig:

                                    peaks_tmp = [peak]
                                    runner    = 0
                                    while True:
                                        runner += 1
                                        next_peak = peaks_all[i+runner]

                                        if next_peak in peaks_trig:
                                            peaks_tmp.append(next_peak)

                                        else:
                                            break

                                    amps_tmp  = [np.abs(stACC_grad[l].data[peak_tmp]) for peak_tmp in peaks_tmp]
                                    best_peak = peaks_tmp[np.argmax(amps_tmp)]
                                    peaks_keep.append(best_peak)

                                    shifter  += len(peaks_tmp) - 1

                        except IndexError:
                            continue

                    peaks_keep            = sorted(set( peaks_keep ))
                    trigger_times         = stACC_grad[l].times('utcdatetime')[peaks_keep]
                    glitch_start_times[l] = trigger_times


                ## SHOW TRIGGERING, if wished
                if plot_triggering:

                    print(u'INFO: Triggering plot: Triggering on derivative of filtered acceleration data is only done on parts not affected by taper.')
                    print(u'INFO: Triggering plot: Amplitudes are scaled so plots are readable, they do not reflect true amplitude values.')
                    print(u'INFO: Triggering plot: Glitches will still be unified across components according to `glitch_min_length`.')
                    print(u'INFO: Triggering plot: Glitches will still be excluded based on their polarization, according to `glitch_min_polarization`.')

                    for p, comp in enumerate(stream_comps):
                        ACC_der = stACC_grad.copy()
                        RAW     = stRAW.copy()
                        RAW.filter('highpass', **{'freq':0.001, 'corners':3, 'zerophase':True})
                        ACC     = stACC.copy()

                        mini    = np.min(ACC_der[p].data)
                        maxi    = np.max(ACC_der[p].data) 

                        ACC_der.normalise( [mini                , maxi                ] )
                        RAW.normalise(     [mini-0.5*(maxi-mini), maxi-0.5*(maxi-mini)] )
                        ACC.normalise(     [mini-1.0*(maxi-mini), maxi-1.0*(maxi-mini)] )  
                    
                        title       = '%s-component, %s preliminary glitches' % (comp, len(glitch_start_times[p]))
                        xlabel      = 'Time (UTC)'
                        data_labels = ('Derivative of filtered ACC', 'RAW original (HP 0.001Hz)', 'filtered ACC (%s)' % ACCfilter['string'])
                        verticals   = glitch_start_times[p]
                        if p!=2:
                            quick_plot(ACC_der[p], RAW[p], ACC[p], win_title=title, title=title, xlabel=xlabel, data_labels=data_labels, horis=[[threshold, -threshold]], verts=[verticals], show=False, keep_open=True)
                        else:
                            quick_plot(ACC_der[p], RAW[p], ACC[p], win_title=title, title=title, xlabel=xlabel, data_labels=data_labels, horis=[[threshold, -threshold]], verts=[verticals], show=True)


                ## UNIFY GLITCH START TIMES ACROSS UVW-COMPONENTS (using parameter: `glitch_min_length` given in seconds)
                print(u'INFO: Unifying glitches across components, using `glitch_min_length`.')

                if all(not list(liste) for liste in glitch_start_times):       # no possible glitches detected
                    continue                    

                glitch_start_times = np.array(glitch_start_times)
                flat_starts        = np.sort( np.hstack(glitch_start_times) )   # np.flatten wouldn't work because sublists contain different number of elements ..
                
                same_glitch_starts = []
                shifter_avoid      = 0
                for l in range( len( flat_starts )):

                    # unify glitch start times: loop glitches, find all +/- `glitch_min_length`, take minimum, go to next not unified glitch (`shifter_avoid`)
                    k = l + shifter_avoid
                    try:
                        indices_same_glitch = np.where( (flat_starts < flat_starts[k]+glitch_min_length) & (flat_starts > flat_starts[k]-glitch_min_length))
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


                ## CALCULATE GLITCH PARAMETERS (e.g. AMPLITUDES, AZIMUTHS, INCIDENCE ANGLES, SNRs, POLARIZATIONS)
                glitch_time_comps = np.array( [[UTCDateTime(time), ''.join(sorted(set(glitch_dict[time]))) ] for time in glitch_dict.keys()] )
                if glitch_time_comps.size==0:
                    print(u'INFO: No glitches found.')
                    continue
                else:
                    print(u'INFO: Verifying glitches, using `glitch_min_polarization`.')
                glitch_time_comps = glitch_time_comps[np.argsort(glitch_time_comps[:,0])]         # time sorted, chronologic

                disregard_counter = 0
                for glitch_UTC_start, comps in glitch_time_comps:

                    ## Prepare data outputs
                    glitch_time_start = marstime(glitch_UTC_start)
                    glitch_time_end   = marstime(glitch_UTC_start+glitch_length)
                    offset            = int(OFFSET_AMP * stRAW_sampling_rates[0])                # 1 second (in samples) offset to avoid glitch precursors 

                    ## Glitch flag per UVW-component 
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

                    ## Polarization check
                    stGAI_window = stGAI.copy()
                    stGAI_window.trim(starttime=glitch_UTC_start, endtime=glitch_UTC_start+glitch_length)
                    stGAI_window.rotate('->ZNE', inventory=stRAW.inventory, components=('UVW'))
                    GAI_measurement = ppol(stream=stGAI_window, demean=True, fix_angles='AMP', Xoffset_samples_for_amplitude=offset)
                    AZI_3D_GAI, AZI_err_3D_GAI, INC_3D_GAI, INC_err_3D_GAI, SNR_3D_GAI, POL_3D_GAI = np.asarray( GAI_measurement.results )[[2,3,6,7,9,13]]

                    # exclude glitches with polarization < `glitch_min_polarization` 
                    # (glitches have high polarization, thus this condition intends to throw out non-glitches)
                    if POL_3D_GAI < float( glitch_min_polarization ):
                        disregard_counter += 1
                        if disregard_counter==len(glitch_time_comps):
                            print(u'INFO: No glitches found.')
                        continue

                    ## Glitch amplitudes
                    stRAW_window = stRAW.copy()
                    stRAW_window.trim(starttime=glitch_UTC_start, endtime=glitch_UTC_start+glitch_length)

                    stZNE = stGAI_window.copy()
                    stZNE.rotate('->ZNE', inventory=stRAW.inventory, components=('UVW'))
                    stZNE.sort(reverse=True)
                    stGAI_window += stZNE    # having in GAI (overall sensitivity and gain corrected): UVWZNE (in that order)    
                    
                    # demean
                    stRAW_window.detrend('demean')
                    stGAI_window.detrend('demean')

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

                    # get UVWZNE peak2peak amplitudes of glitches                
                    U_RAW  = np.abs (data_U_RAW[offset:][np.argmax(data_U_RAW[offset:])] - data_U_RAW[offset:][np.argmin(data_U_RAW[offset:])] )
                    V_RAW  = np.abs (data_V_RAW[offset:][np.argmax(data_V_RAW[offset:])] - data_V_RAW[offset:][np.argmin(data_V_RAW[offset:])] )
                    W_RAW  = np.abs (data_W_RAW[offset:][np.argmax(data_W_RAW[offset:])] - data_W_RAW[offset:][np.argmin(data_W_RAW[offset:])] )

                    U_GAI  = np.abs (data_U_GAI[offset:][np.argmax(data_U_GAI[offset:])] - data_U_GAI[offset:][np.argmin(data_U_GAI[offset:])] )
                    V_GAI  = np.abs (data_V_GAI[offset:][np.argmax(data_V_GAI[offset:])] - data_V_GAI[offset:][np.argmin(data_V_GAI[offset:])] )
                    W_GAI  = np.abs (data_W_GAI[offset:][np.argmax(data_W_GAI[offset:])] - data_W_GAI[offset:][np.argmin(data_W_GAI[offset:])] )
                    Z_GAI  = np.abs (data_Z_GAI[offset:][np.argmax(data_Z_GAI[offset:])] - data_Z_GAI[offset:][np.argmin(data_Z_GAI[offset:])] )
                    N_GAI  = np.abs (data_N_GAI[offset:][np.argmax(data_N_GAI[offset:])] - data_N_GAI[offset:][np.argmin(data_N_GAI[offset:])] )
                    E_GAI  = np.abs (data_E_GAI[offset:][np.argmax(data_E_GAI[offset:])] - data_E_GAI[offset:][np.argmin(data_E_GAI[offset:])] )

                    # Exclude glitches already detected due to small overlap of sliding windows (should not happen often)                    
                    if glitches:
                        if not all(glitch_time_start.UTC_time>=UTCDateTime(glitch_already[1])+glitch_min_length or 
                                   glitch_time_start.UTC_time<=UTCDateTime(glitch_already[1])-glitch_min_length for glitch_already in glitches):
                            continue

                    # Collect all glitch measures                                        
                    glitch_counter += 1                
                    glitch = glitch_counter, \
                             glitch_time_start.UTC_string,  \
                             glitch_time_end.UTC_string,    \
                             glitch_time_start.LMST_string, \
                             glitch_time_end.LMST_string,   \
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
                             AZI_3D_GAI,     \
                             AZI_err_3D_GAI, \
                             INC_3D_GAI,     \
                             INC_err_3D_GAI, \
                             SNR_3D_GAI,     \
                             POL_3D_GAI

                    glitches_text.append(line_formatter % glitch)
                    glitches.append(glitch)

                    # Output to user
                    if glitch_counter%5 == 0:

                        if print_counter == 0:
                            print_counter += 1
                            print()
                            for header_line in header.split('\n')[1:]:
                                print('# ' + header_line[:108])

                        print( (line_formatter % glitch)[:110] )


        ### STATISTICS
        statistic = glitch_statistic(glitches, total_data_time_UTC)
        try:
            cum_time_sols = 1/(total_data_time_LMST/(3600*24.))
        except ZeroDivisionError:
            cum_time_sols = np.nan


        ### SUMMARY TEXT
        output    = [u'DATA:',
                     u'  path:                       %s'           % RAW_UVW[o],
                     u'  IDs:                        %s'           % ', '.join(stream_ids),
                     u'  UTC start:                  %s'           % stream_times[0].UTC_string,
                     u'  UTC end:                    %s'           % stream_times[1].UTC_string,
                     u'  UTC length (without gaps):  %s (h:m:s)'   % sec2hms( total_data_time_UTC ),
                     u'  LMST start:                 %s'           % stream_times[0].LMST_string,
                     u'  LMST end:                   %s'           % stream_times[1].LMST_string,
                     u'  LMST length (without gaps): %s (h:m:s)'   % sec2hms( total_data_time_LMST ).replace('days', 'sols').replace('day','sol'),
                     u'',      
                     u'PARAMETERS:',      
                     u'  Taper length per side %%:    %s'          % taper_length_per_side,
                     u'  Pre filtering:              %s'           % str(pre_filt),
                     u'  Water level:                %s'           % str(water_level),
                     u'  ACCfilter:                  %s'           % ACCfilter['type'],
                     u'                              %s'           % ACCfilter['string'],
                     u'  Threshold:                  %s m/s**3'    % float(threshold),
                     u"  Show triggering:            %s"           % str(plot_triggering),
                     u"  Glitches' min dist:         %s s"         % glitch_min_length,
                     u"  Glitch length (fix):        %s s"         % glitch_length,
                     u'  Glitch min polarization:    %s'           % glitch_min_polarization,
                     u'',      
                     u'RESULTS:',      
                     u'    Glitches all:             %s'           % statistic[0],
                     u'       all / sol:             %.0f'         % statistic[1],
                     u'            on U:             %s'           % statistic[2],
                     u'            on V:             %s'           % statistic[3],
                     u'            on W:             %s'           % statistic[4],
                     u'       only on U:             %s'           % statistic[5],
                     u'       only on V:             %s'           % statistic[6],
                     u'       only on W:             %s'           % statistic[7],
                     u'    indi./ all %%:             %.1f'        % statistic[8],
                     u'',
                     u'Done in:   %s (h:m:s), %s / sol.'           % (sec2hms( time.time()-now ), sec2hms( (time.time()-now)/cum_time_sols )),
                     u'Timestamp: %s'                              % datetime.datetime.now().strftime('%Y-%m-%dT%H:%M:%S')]

        # create good-looking output for later
        len_max_line  = max( [len(line) for line in output] )
        formatter_str = "\n#   | %-" + str(len_max_line) + 's |'

        output = [formatter_str % line for line in output]
        output.insert(0,           '\n\n#   +-' +  '-' * len_max_line + '-+'   )
        output.insert(len(output),   '\n#   +-' +  '-' * len_max_line + '-+\n' )


        ### FILING GLITCHES + OUTPUT
        try:
            sort_indices  = np.array( glitches )[:,1].argsort()         # glitches will always be sorted according to starttime !
            glitches_text = np.asarray(glitches_text)[sort_indices]
        except IndexError:                                              # case e.g., if there are no glitches detected.
            pass
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
        print(u'OUTFILE GLITCH DETECTOR:')
        print(output_txt)



### _ _ N A M E _ _ = = " _ _ M A I N _ _ "  
if __name__ == "__main__": 
    pass