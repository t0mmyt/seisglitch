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
def synthetic_step(component_response, sampling_period, total_length, step_unit='ACC', step_amp=1e-9, prepend_straight_in_s=1):


    # Variables
    num_samples = 10000
    signal      = np.zeros(4*num_samples)
    step        = np.hstack(( np.zeros(num_samples), 
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
    first_index = np.where(signal>1e-2)[0][0]
    signal      = signal[first_index:]


    # Prepend to signal (can be better for fitting actual data, especially for glitches)
    if prepend_straight_in_s:
        prepend_ones = np.ones( int(prepend_straight_in_s/sampling_period) )
        signal       = np.hstack(( prepend_ones*0, signal )) 


    # Cut again to fixed "total_length" as we need non-varying synthetic signal length
    signal = signal[:int(total_length/sampling_period)]

    return signal
def remove(*glitch_detector_files, 
    waveform_files        = [], 
    inventory_file        = 'IRIS', 
    glitch_length         = 40, 
    var_reduction_min     = 85, 
    shift_time_s          = 10, 
    lowpass_Hz            = 1, 
    show_removed_glitch   = False,
    **kwargs):


    now = time.time()



    ### OUTPUT
    print()
    print(u'  ----------------------')
    print(u'  RUNNING GLITCH REMOVER')
    print(u'  ----------------------')
    print()



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
        for trace in stream:

            def glitch_model(x, m, n, o, dis_step, prepend):

                """
                Five fit-variables.
                #glitch_model(x, m, n, o, dis_step, prepend)
                """
                syn_glitch  = synthetic_glitch(response, sampling_period, glitch_length, acceleration_step=o, displacement_step=dis_step, prepend_straight_in_s=prepend)
                return x * m + n + syn_glitch            
            print()

            # synthetic glitch
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

            # more data prep
            if isinstance(lowpass_Hz, (float,int)) and not isinstance(lowpass_Hz, bool) and 0.5*sampling_rate > float(lowpass_Hz) :
                trace_filter = trace.copy()
                trace_filter.detrend('demean')        
                trace_filter.taper(0.03)       
                trace_filter.detrend('simple') 
                trace_filter.filter('lowpass', freq=float(lowpass_Hz), corners=3, zerophase=True)
                filtered = True
            else:
                trace_filter = trace.copy()
                filtered     = False

            # looping over glitches to be corrected
            for g in range(len( glitches )):

                glitch_number      = int(glitches[g][0].replace(':',''))
                glitch_start       = marstime(glitches[g][1])
                glitch_end         = marstime(glitches[g][2])
                #if component!='V':
                #    continue     

                
                len_glitch         = int(glitch_length/sampling_period)
                x_range_glitch     = np.arange(len_glitch)
                trace_filter_slice = trace_filter.slice(starttime=glitch_start.UTC_time-shift_time_s, endtime=glitch_end.UTC_time+shift_time_s)
                trace_slice        = trace.slice(       starttime=glitch_start.UTC_time-shift_time_s, endtime=glitch_end.UTC_time+shift_time_s)

                original           = trace_slice.copy()
                len_data           = len(trace_filter_slice.data)
                
                residuals          = []
                fits               = []
                popts              = []

                # looping over each data point for glitch correction
                for i in range(len_data-len_glitch):

                    # DATE
                    data_shifted = trace_filter_slice.data[i:i+len_glitch]

                    # FIT-VARIABLES
                    # x, m, n, o, dis_step, prepend
                    #p0       = [0,data_shifted[0],1]
                    #bounds   = ([-np.inf, np.min(data_shifted),0],[np.inf, np.max(data_shifted),np.inf])
                    p0       = [0,data_shifted[0],1e-9,1e-11,1]
                    #bounds   = ([-np.inf, np.min(data_shifted),0,0,0],[np.inf, np.max(data_shifted),np.inf,np.inf,np.min([shift_time_s,5])])
                    popt, _  = scipy.optimize.curve_fit(glitch_model, x_range_glitch, data_shifted, p0=p0)

                    # FIT
                    #residual = np.linalg.norm(data_shifted - (fit - x_range_glitch * popt[0] - popt[1]))
                    fit      = glitch_model(x_range_glitch, *popt)
                    residual = np.linalg.norm(data_shifted - (fit))
                    residuals.append(residual)
                    popts.append(popt)
                    fits.append(fit)

                # best fit
                best_index    = np.array( residuals ).argmin()
                best_popt     = popts[best_index]
                best_fit      = fits[best_index]
                shift         = best_index
                scaled_glitch = best_fit - x_range_glitch * best_popt[0] - best_popt[1]

                # variance reduction
                var_data           = np.var(trace_slice)
                var_data_corrected = np.var( np.hstack((trace_slice.data[:shift], trace_slice.data[shift:shift+len(scaled_glitch)] - scaled_glitch, trace_slice.data[shift+len(scaled_glitch):])) )
                var_red            = (1 - var_data_corrected/var_data) * 100
                if var_red>100:
                    var_red = 100
                if var_red<0:
                    var_red = 0

                # results
                if component.upper()=='U':
                    if var_red>=var_reduction_min or glitches[g][5]=='1':
                        var_red_U.append( var_red )
                        acc_steps_U.append( best_popt[2] )
                elif component.upper()=='V':
                    if var_red>=var_reduction_min or glitches[g][6]=='1':
                        var_red_V.append( var_red )
                        acc_steps_V.append( best_popt[2] )
                elif component.upper()=='W':
                    if var_red>=var_reduction_min or glitches[g][7]=='1':
                        var_red_W.append( var_red )
                        acc_steps_W.append( best_popt[2] )
                else:
                    pass

                # actual data correction
                if var_red>=var_reduction_min:

                    print(u'Glitch %06d,  %s,  %s,  var_red = %7.1f %%,  Correction is done.  Acceleration step = %6.1f nm/s**2'  % (glitch_number, glitch_start.UTC, component, var_red, best_popt[2]))

                    trace_slice.data[shift:shift+len(scaled_glitch)] = trace_slice.data[shift:shift+len(scaled_glitch)]-scaled_glitch
                    removed.append(component)
                    glitch_trace = original.copy()
                    glitch_trace.data = np.hstack(( np.nan*np.ones(shift), best_fit ))

                    if show_removed_glitch:
                        if filtered:
                            title             = 'Glitch %s on component %s (var_red=%.1f%%)' % (glitch_number, component, var_red)
                            xlabel            = 'Time (UTC)'
                            data_labels       = ('original', 'original filtered (lowpass %s Hz)' % lowpass_Hz, 'fitted glitch (on filtered)', 'deglitched')
                            quick_plot(original, trace_filter_slice, glitch_trace, trace_slice, title=title, xlabel=xlabel, data_labels=data_labels,)
                        else:
                            title             = 'Glitch %s on component %s (var_red=%.1f%%)' % (glitch_number, component, var_red)
                            xlabel            = 'Time (UTC)'
                            data_labels       = ('original', 'fitted glitch (on original)', 'deglitched')
                            quick_plot(original, glitch_trace, trace_slice, title=title, xlabel=xlabel, data_labels=data_labels)

                else:
                    print(u'Glitch %06d,  %s,  %s,  var_red = %7.1f %%.  No correction done.  Acceleration step = %6.1f nm/s**2'  % (glitch_number, glitch_start.UTC, component, var_red, best_popt[2]))


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

    stream          = read2('/home/scholz/Desktop/data/XB.ELYSE.00.HH?_sol423night.mseed', headonly=True)
    #stream.trim2(0,0.01)
    
    stream.set_inventory('IRIS')
    response        = stream.inventory[0][0][0].response
    sampling_period = stream[0].stats.delta

    syn_glitch1     = Trace(data=synthetic_step(response, sampling_period, 10, step_unit='ACC', step_amp=1e-9)+10, header={'delta':sampling_period})
    syn_glitch2     = Trace(data=synthetic_step(response, sampling_period, 20, step_unit='ACC', step_amp=1e-9)+20, header={'delta':sampling_period})
    syn_glitch3     = Trace(data=synthetic_step(response, sampling_period, 30, step_unit='ACC', step_amp=1e-9)+30, header={'delta':sampling_period})
    syn_glitch4     = Trace(data=synthetic_step(response, sampling_period, 40, step_unit='ACC', step_amp=1e-9)+40, header={'delta':sampling_period})
    syn_glitch5     = Trace(data=synthetic_step(response, sampling_period, 50, step_unit='ACC', step_amp=1e-9)+50, header={'delta':sampling_period})
    
    glitch_stream   = Stream2(traces=[syn_glitch1, 
                                      syn_glitch2, 
                                      syn_glitch3, 
                                      syn_glitch4, 
                                      syn_glitch5, 
                                      ])
        
    quick_plot(*glitch_stream, x_label='Time', y_label='Digital Units', data_labels=['10 s', '20 s', '30 s', '40 s', '50 s'])