#!/usr/bin/env ipython

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
import time
import scipy
import datetime
import numpy as np


#####  obspy modules import  #####
from obspy.core.trace import Trace


#####  toolbox modules import  #####
from seisglitch.util import read2, Stream2, marstime, quick_plot, sec2hms



### GLTICH REMOVAL
def synthetic_step(component_response, sampling_period, total_length, step_unit='ACC', step_amp=1e-9):


    # Variables
    num_samples  = 10000
    total_length = min(total_length, num_samples*sampling_period)       # final signal as long as possible but not longer than distance between two steps (box-car usage)
    signal       = np.zeros(4*num_samples)
    step         = np.hstack(( np.zeros(num_samples), 
                               np.ones(num_samples),
                              -np.ones(num_samples),
                               np.zeros(num_samples) ))


    # FFT
    step_fft    = np.fft.fft(step)
    step_freqs  = np.fft.fftfreq(step.size, d=sampling_period)

    if step_unit.upper()=='ACC':
        # glitch
        resp_fft = component_response.get_evalresp_response_for_frequencies(step_freqs, output='ACC')
    elif step_unit.upper()=='VEL':
        # no terminology defined
        resp_fft = component_response.get_evalresp_response_for_frequencies(step_freqs, output='VEL')
    elif step_unit.upper()=='DIS' or step_unit.upper()=='DISP':
        # precursor
        resp_fft = component_response.get_evalresp_response_for_frequencies(step_freqs, output='DISP')


    # Synthezising signal
    signal_fft  = step_fft * resp_fft
    signal     += np.fft.ifft(signal_fft).real
    signal     *= float(step_amp)


    # Cut only away zeroes in front
    peak2peak_signal = np.max(signal)-np.min(signal)
    first_index      = np.where(signal>peak2peak_signal/1000.)[0][0]  # factor 1000 has influence on how many non-"zero" values are before actual glitch starts. For 2 SPS data these are more than for higher SPS. Thus influences precursor fits somewhat. 1000 is an okay trade-off.
    if first_index>0:
        first_index -= 1
    signal           = signal[first_index:]

    # Cut again to fixed "total_length" as we need non-varying synthetic signal length
    signal = signal[:int(total_length/sampling_period)]

    return signal
def remove(*glitch_detector_files, 
    waveform_files         = [], 
    inventory_file         = 'IRIS', 
    glitch_length          = 40,
    glitch_shift_time_s    = 10, 
    precursor_fit          = False,
    var_reduction_min      = 85, 
    show_glitch_fit        = False,
    plot_removal_statistic = False, 
    **kwargs):


    now = time.time()



    ### OUTPUT
    print()
    print(u'  ----------------------')
    print(u'  RUNNING GLITCH REMOVER')
    print(u'  ----------------------')
    print()



    ### FIXED PARAMETERS
    PREPEND_STRAIGHT                 = min(1,glitch_shift_time_s)  # in s
    ACC_STEP                         = 1e-9                        # in m/s**2, for glitches
    DIS_STEP                         = 1e-12                       # in m, for glitch precursors
    PRECURSOR_LENGTH_SAMPLES         = 35                          # should be >25 and smaller than glitch_length*sampling_period
    PRECURSOR_SHIFT_SAMPLES_PER_SIDE = 7                           # maximum of precursor shifted left and right w.r.t. determined fitted glitch onset (should be larger than samples modeled glitch signal before real glitch starts)
    VAR_REDUCTION_MIN_PRECURSOR      = 2                           # in %



    ### READ GLITCH-FILE
    all_glitches = []

    for glitch_detector_file in glitch_detector_files:
        glitches      = np.loadtxt(glitch_detector_file, dtype='str')
        all_glitches += list(glitches)
    all_glitches = np.array(all_glitches)
    


    ### OTHER VARIABLES
    assign      = {'U':5, 'V':6, 'W':6}
    var_red_U   = []
    var_red_V   = []
    var_red_W   = []
    acc_steps_U = []
    acc_steps_V = []
    acc_steps_W = []



    ### READ WAVEFORM FILES
    for o, waveform_file in enumerate(waveform_files):
        removed  = []

        # read streams
        stream = read2(waveform_file)
        stream.sort(reverse=False)

        # small output
        print(u'Info: Analysing file: %s/%s' % (o+1,len(waveform_files)))
        print(waveform_file)
        for trace in stream:
            print('  %s' % trace)
        print()
        print(u'Info: Handling: %s' % waveform_file)

        stream.set_inventory(inventory_file)

        # loop traces
        for trace in stream.select(channel='?[LMH]?'):

            print()

            # data prep
            component = trace.stats.channel[-1]
            glitches  = all_glitches[ (all_glitches[:,1]>=str(trace.stats.starttime)) & (all_glitches[:,2]<=str(trace.stats.endtime)) ]
            inv = stream.inventory.select(network   = trace.stats.network, 
                                          station   = trace.stats.station, 
                                          location  = trace.stats.location, 
                                          channel   = trace.stats.channel, 
                                          starttime = trace.stats.starttime, 
                                          endtime   = trace.stats.endtime)
            response         = inv[0][0][0].response
            sampling_period  = trace.stats.delta
            sampling_rate    = trace.stats.sampling_rate

            # synthetic glitch generation, for each new trace once
            prepend    = int(PREPEND_STRAIGHT/sampling_period)
            syn_glitch = synthetic_step(response, sampling_period, glitch_length, step_unit='ACC', step_amp=ACC_STEP)
            syn_glitch = np.hstack(( np.zeros(prepend), syn_glitch ))
            def glitch_model(x, m, n, o):
                """
                Three fit-variables.
                """
                return x * m + n + syn_glitch * o

            # synthetic glitch precursor generation, for each new trace once
            if precursor_fit:
                syn_precur          = synthetic_step(response, sampling_period, PRECURSOR_LENGTH_SAMPLES*sampling_period, step_unit='DIS', step_amp=DIS_STEP)
                precursor_index_max = np.argmax(np.abs(syn_precur))
                def precursor_model(x, m, n, o):
                    """
                    Three fit-variables.
                    """
                    return x * m + n + syn_precur * o

            # looping over glitches to be corrected
            for g in range(len( glitches )):

                # glitch variables
                glitch_number  = int(glitches[g][0].replace(':',''))
                glitch_start   = marstime(glitches[g][1])
                glitch_end     = marstime(glitches[g][2])
                glitch_len     = prepend+int(glitch_length/sampling_period)

                # variables needed
                x_range_glitch = np.arange(glitch_len)
                raw_slice      = trace.slice(starttime=glitch_start.UTC_time-glitch_shift_time_s, endtime=glitch_end.UTC_time+glitch_shift_time_s)
                original_slice = raw_slice.copy()
                data_len       = len(raw_slice.data)
                residuals      = []
                residuals2     = []
                fits           = []
                fits2          = []
                popts          = []
                popts2         = []

                # looping over each data point for glitch correction
                for i in range(data_len-glitch_len):

                    # DATA
                    data_shifted = raw_slice.data[i:i+glitch_len]

                    # FIT with variables: x, m, n, o
                    p0       = [0,data_shifted[0],1]
                    #bounds   = (-np.inf, np.inf)               # is Scipy default
                    bounds   = ([-1, np.min(data_shifted), -np.inf],[1, np.max(data_shifted), np.inf])
                    popt, _  = scipy.optimize.curve_fit(glitch_model, x_range_glitch, data_shifted, p0=p0, bounds=bounds)
                    fit      = glitch_model(x_range_glitch, *popt)
                    residual = np.linalg.norm(data_shifted - fit)

                    # storing needed results
                    residuals.append(residual)
                    popts.append(popt)
                    fits.append(fit)

                # best fit glitch
                best_index    = np.array( residuals ).argmin()
                best_popt     = popts[best_index]
                best_fit      = fits[best_index]
                shift         = best_index
                scaled_glitch = best_fit - x_range_glitch * best_popt[0] - best_popt[1]

                # actual glitch correction!
                raw_slice.data[shift:shift+len(scaled_glitch)] = raw_slice.data[shift:shift+len(scaled_glitch)]-scaled_glitch

                # looping over each data point for glitch precursor correction
                tag_precursor = False
                if precursor_fit:
                    precur_len     = PRECURSOR_LENGTH_SAMPLES
                    x_range_precur = np.arange(precur_len)

                    for i in range(2*PRECURSOR_SHIFT_SAMPLES_PER_SIDE+1):

                        # DATA
                        start_index  = max(shift+prepend-precursor_index_max-PRECURSOR_SHIFT_SAMPLES_PER_SIDE+i,0)      # avoid negative numbers as start index, which might happen
                        end_index    = start_index+precur_len
                        data_shifted = raw_slice.data[start_index:end_index]

                        # FIT with variables: x, m, n, o
                        p0 = [0,data_shifted[0],1]

                        #bounds   = (-np.inf, np.inf)               # is Scipy default
                        bounds   = ([-1, np.min(data_shifted), -np.inf],[1, np.max(data_shifted), np.inf])
                        popt, _  = scipy.optimize.curve_fit(precursor_model, x_range_precur, data_shifted, p0=p0, bounds=bounds)
                        fit      = precursor_model(x_range_precur, *popt)
                        residual = np.linalg.norm(data_shifted - fit)

                        # storing needed results
                        residuals2.append(residual)
                        popts2.append(popt)
                        fits2.append(fit)

                    # best fit precursor
                    best_index2   = np.array( residuals2 ).argmin()
                    best_popt2    = popts2[best_index2]
                    best_fit2     = fits2[best_index2]
                    shift2        = shift+prepend-precursor_index_max-PRECURSOR_SHIFT_SAMPLES_PER_SIDE+best_index2
                    scaled_precur = best_fit2 - x_range_precur * best_popt2[0] - best_popt2[1]

                    # actual glitch precursor correction!
                    deglitched_slice = raw_slice.copy()
                    raw_slice.data[shift2:shift2+len(scaled_precur)] = raw_slice.data[shift2:shift2+len(scaled_precur)]-scaled_precur

                    # variance reduction between deglitched and deglitched+deprecursored data
                    var_data           = np.var(deglitched_slice)
                    var_data_corrected = np.var(raw_slice)
                    var_red            = (1 - var_data_corrected/var_data) * 100
                    if var_red>100:
                        var_red = 100
                    if var_red<0:
                        var_red = 0

                    # precursor correction undone or not
                    if var_red>=VAR_REDUCTION_MIN_PRECURSOR:
                        tag_precursor = True
                    else:
                        tag_precursor = False
                        raw_slice.data = deglitched_slice.data

                # variance reduction between original and deglitched or deglitched+deprecursored data
                var_data           = np.var(original_slice)
                var_data_corrected = np.var(raw_slice)
                var_red            = (1 - var_data_corrected/var_data) * 100
                if var_red>100:
                    var_red = 100
                if var_red<0:
                    var_red = 0

                # glitch correction undone or not
                if var_red>=var_reduction_min:
                    removed.append(component)
                    if precursor_fit:
                        if tag_precursor:
                            print(u'Glitch %06d,  %s,  %s,  var_red = %4.1f %%,  Correction is done.  Acc. step = %6.1f nm/s**2.  Dis. step = %6.1f pm'  % (glitch_number, glitch_start.UTC_string, component, var_red, best_popt[2], best_popt2[2]))
                            label = 'glitch+precursor corrected'
                        else:
                            print(u'Glitch %06d,  %s,  %s,  var_red = %4.1f %%,  Correction is done.  Acc. step = %6.1f nm/s**2'  % (glitch_number, glitch_start.UTC_string, component, var_red, best_popt[2]))
                            label = 'glitch corrected'
                    else:
                        print(u'Glitch %06d,  %s,  %s,  var_red = %4.1f %%,  Correction is done.  Acc. step = %6.1f nm/s**2'  % (glitch_number, glitch_start.UTC_string, component, var_red, best_popt[2]))
                        label = 'glitch corrected'
                else:
                    print(u'Glitch %06d,  %s,  %s,  var_red = %4.1f %%.  No correction done.'  % (glitch_number, glitch_start.UTC_string, component, var_red))
                    raw_slice.data = original_slice.data                    
                    label = 'not corrected'

                # only to convey information to user
                if show_glitch_fit:
                    glitch_slice = original_slice.copy()
                    glitch       = np.hstack(( np.nan*np.ones(shift), best_fit, np.nan*np.ones(data_len-glitch_len-shift) ))          # glitch
                    if precursor_fit and tag_precursor:
                        precursor = np.hstack(( np.nan*np.ones(shift2), best_fit2-best_fit2[1], np.nan*np.ones(data_len-precur_len-shift2) ))   # precursor
                    else:
                        precursor = np.nan*np.ones(data_len)

                    data_plot = []
                    for i in range(len(glitch)):
                        if np.isnan(glitch[i]) and np.isnan(precursor[i]):
                            value = np.nan
                        elif not np.isnan(glitch[i]) and np.isnan(precursor[i]):
                            value = glitch[i]
                        elif np.isnan(glitch[i]) and not np.isnan(precursor[i]):
                            value = precursor[i]+best_fit[1]
                        else:
                            value = glitch[i] + precursor[i]
                        data_plot.append( value )

                    glitch_slice.data  = np.array( data_plot )
                    title              = 'Glitch %s on component %s (var_red=%.1f%%, %s)' % (glitch_number, component, var_red, label)
                    data_labels        = ('original', 'corrected', 'fit (on original)')
                    quick_plot(original_slice, raw_slice, glitch_slice, title=title, xlabel='Time UTC (s)', ylabel='Raw Units', lw=[1,1,1.5], lc=[None, None, 'k'], data_labels=data_labels)

                # results
                if component.upper()=='U':
                    if var_red>=var_reduction_min or glitches[g][5]=='1':
                        var_red_U.append( var_red )
                        acc_steps_U.append( best_popt[2]*ACC_STEP )
                elif component.upper()=='V':
                    if var_red>=var_reduction_min or glitches[g][6]=='1':
                        var_red_V.append( var_red )
                        acc_steps_V.append( best_popt[2]*ACC_STEP )
                elif component.upper()=='W':
                    if var_red>=var_reduction_min or glitches[g][7]=='1':
                        var_red_W.append( var_red )
                        acc_steps_W.append( best_popt[2]*ACC_STEP )
                else:
                    pass


        ## WRITE DEGLITCHED FILE
        outfile = '.'.join(waveform_file.split('.')[:-1]) + '_deglitched.' + waveform_file.split('.')[-1]
        stream.write(outfile)


        ## OUTPUT
        string       = ''
        removed      = np.array(removed)
        for comp in sorted(set(removed)):
            string += '%s:%s, ' % (comp, len(removed[removed==comp]))
        string = '(' + string[:-2] + ')'

        print()
        print(u'Removed a total of %s individual glitches on all traces.' % len(removed))
        print(string)
        print()
        print(u'DEGLITCHED FILE:')
        print(outfile)
        print()
        print(u'Done in:   %s (h:m:s).' % sec2hms( time.time()-now ))
        print(u'Timestamp: %s'          % datetime.datetime.now().strftime('%Y-%m-%dT%H:%M:%S'))


    ### FINAL STATISTIC
    if plot_removal_statistic:
        # variables
        var_red_U   = np.array(var_red_U)
        var_red_V   = np.array(var_red_V)
        var_red_W   = np.array(var_red_W)
        acc_steps_U = np.array(acc_steps_U)
        acc_steps_V = np.array(acc_steps_V)
        acc_steps_W = np.array(acc_steps_W)
        
        # Plot 1
        sum_acc_U                 = np.sum(np.abs( acc_steps_U ))
        sum_acc_V                 = np.sum(np.abs( acc_steps_V ))
        sum_acc_W                 = np.sum(np.abs( acc_steps_W ))
        
        norm_cumul_accerl_U       = [ np.sum( np.abs(acc_steps_U[np.where(var_red_U<=var_red)] ))/sum_acc_U for var_red in range(0,101)]
        norm_cumul_accerl_V       = [ np.sum( np.abs(acc_steps_V[np.where(var_red_V<=var_red)] ))/sum_acc_V for var_red in range(0,101)]
        norm_cumul_accerl_W       = [ np.sum( np.abs(acc_steps_W[np.where(var_red_W<=var_red)] ))/sum_acc_W for var_red in range(0,101)]
        
        percent_below_threshold_U = np.sum( np.abs(acc_steps_U[np.where(var_red_U>=var_red)] ))/sum_acc_U
        percent_below_threshold_V = np.sum( np.abs(acc_steps_V[np.where(var_red_V>=var_red)] ))/sum_acc_V
        percent_below_threshold_W = np.sum( np.abs(acc_steps_W[np.where(var_red_W>=var_red)] ))/sum_acc_W
        
        quick_plot(norm_cumul_accerl_U, norm_cumul_accerl_V, norm_cumul_accerl_W, win_title='Norm. Cum. Acc.', data_labels=['U (%3.1f of total acceleration above threshold)' % percent_below_threshold_U, 'V (%3.1f of total acceleration above threshold)' % percent_below_threshold_V,'W (%3.1f of total acceleration above threshold)' % percent_below_threshold_W], x=np.arange(0,101), xlim=[0,100], ylim=[0,1], xlabel='Variance Reduction (%)', ylabel='Normalized Cumulative Accleration', verts=[[85]])
        
        # Plot 2
        var_red_U_sorted          = var_red_U[np.argsort(var_red_U)]
        var_red_V_sorted          = var_red_V[np.argsort(var_red_V)]
        var_red_W_sorted          = var_red_W[np.argsort(var_red_W)]
        
        percent_below_threshold_U = len(var_red_U[var_red_U>=85])/len(var_red_U) * 100
        percent_below_threshold_V = len(var_red_V[var_red_V>=85])/len(var_red_V) * 100
        percent_below_threshold_W = len(var_red_W[var_red_W>=85])/len(var_red_W) * 100
        
        quick_plot(var_red_U_sorted, var_red_V_sorted, var_red_W_sorted, win_title='Variance Reduction', data_labels=['U (%4.1f%% above threshold)' % percent_below_threshold_U, 'V (%4.1f%% above threshold)' % percent_below_threshold_V,'W (%4.1f%% above threshold)' % percent_below_threshold_W], xlim=[1,np.max([len(var_red_U_sorted), len(var_red_V_sorted), len(var_red_W_sorted)])], ylim=[0,100], xlabel='Glitch index', ylabel='Variance Reduction (%)', horis=[[85]])


### _ _ N A M E _ _ = = " _ _ M A I N _ _ "  
if __name__ == "__main__":


    files = ['/home/scholz/Desktop/data/XB.ELYSE.02.MH?_2019-03-15T06:28:06.779000Z-2019-03-15T08:32:53.279000Z_raw.MSEED',
             '/home/scholz/Desktop/data/XB.ELYSE.02.BH?_raw_S0325a_QB.MSEED',
             '/home/scholz/Desktop/data/XB.ELYSE.00.HH?_sol423night.mseed',
             '/home/scholz/Desktop/data/XB.ELYSE.67.MH?_2019-03-08T03:56:39.154000Z-2019-03-08T04:30:59.154000Z_raw.MSEED',
             '/home/scholz/Desktop/data/XB.ELYSE.67.SH?_2019-03-01T00:00:09.731000Z-2019-03-01T12:00:21.214000Z_raw.MSEED',
             '/home/scholz/Desktop/data/XB.ELYSE.65.EH?_2019-03-08T03:59:59.314000Z-2019-03-08T04:30:01.944000Z_raw.MSEED']

    traces = []
    for file in files:
        stream = read2(file, headonly=True)
        #stream.trim2(0,0.01)
        stream.set_inventory('IRIS')
        response        = stream.inventory[0][0][0].response
        sampling_period = stream[0].stats.delta
        header          = {'delta'    : sampling_period, 
                           'network'  : stream[0].stats.network, 
                           'station'  : stream[0].stats.station, 
                           'location' : stream[0].stats.location, 
                           'channel'  : stream[0].stats.channel}

        syn_glitch = Trace(data=synthetic_step(response, sampling_period, 25, step_unit='ACC', step_amp=1e-9), header=header)
        syn_glitch = Trace(data=syn_glitch.data+synthetic_step(response, sampling_period, 25, step_unit='DIS', step_amp=1e-12), header=header)
        traces.append(syn_glitch)
    
    glitch_stream = Stream2(traces=traces)
    print(glitch_stream)
    #quick_plot(*[trace.data for trace in glitch_stream], y_label='Digital Units', data_labels=('VBB 2 SPS', 'VBB 20 SPS', 'VBB 100 SPS', 'SP 2SPS', 'SP 20 SPS', 'SP 100 SPS', ))
    quick_plot(*glitch_stream, x_label='Relative Time (s)', y_label='Digital Units')    