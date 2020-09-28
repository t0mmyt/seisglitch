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
import matplotlib.pyplot as plt


#####  obspy modules import  #####
import obspy
from obspy.core.trace import Trace


#####  toolbox modules import  #####
from seisglitch.util import read2, Stream2, marstime, quick_plot, sec2hms



### GLTICH REMOVAL
def sourceFFT(num_samples, sampling_period, int_samples=None, gauss=False, abkling=None, amp=None, freq=None):

    # Instantaneous step or finite ramp
    if int_samples:
        if int_samples%2 != 0:
            int_samples += 1
        if not gauss:
            interp = np.array( [(i+1)*1/(int_samples+1) for i in range(int_samples)] )
        else:
            interp = scipy.signal.gaussian(2*int_samples, std=int_samples/4)[:int_samples]
        step = np.hstack(( np.zeros(num_samples-int(int_samples/2)),
                           interp,
                           np.ones(num_samples-int(int_samples/2)) ))
        #step = np.hstack(( np.zeros(num_samples-int(int_samples)),
        #                   interp[::-1],
        #                   interp,
        #                   np.zeros(num_samples-int(int_samples)) ))  
        #quick_plot(step)

    else:
        step = np.hstack(( np.zeros(num_samples), 
                           np.ones(num_samples) ))

    # Additional osscillation
    if abkling:
        x        = np.arange(num_samples)
        omega    = np.pi*2*freq
        func     = amp * np.exp(-x*sampling_period*abkling) * np.cos(2*np.pi*freq*sampling_period * x)
        step[:] += np.hstack(( func, func ))

    # FFT
    step_freqs  = np.fft.fftfreq(step.size, d=sampling_period)  # independent of actual input, only of its length
    source_FFT  = np.fft.fft(step)

    return source_FFT, step_freqs
def responseFFT(component_response, freqs, step_unit='ACC'):

    # FFT depending on 'step_unit'
    if step_unit.upper()=='ACC':
        # glitch
        resp_fft = component_response.get_evalresp_response_for_frequencies(freqs, output='ACC')
    elif step_unit.upper()=='VEL':
        # no terminology defined
        resp_fft = component_response.get_evalresp_response_for_frequencies(freqs, output='VEL')
    elif step_unit.upper()=='DIS' or step_unit.upper()=='DISP':
        # spike
        resp_fft = component_response.get_evalresp_response_for_frequencies(freqs, output='DISP')

    return resp_fft, freqs
def fft2signal(fft, freqs, tau=0, amp_factor=1, real=True):

    # 
    omega   = 2*np.pi*freqs
    fft    *= np.exp(-1j*omega*tau)         # positive 'tau' produces shift to the right
    signal  = np.fft.ifft(fft)
    signal *= float(amp_factor)
    if real:
        signal = signal.real

    return signal

def remove(*glitch_detector_files, 
    waveform_files               = [], 
    inventory_file               = 'IRIS',
    glitch_window_leftright      = 2,
    glitch_prepend_zeros         = 0,
    glitch_interpolation_samples = 0,
    glitch_subsample_factor      = 1,
    spike_fit                    = False,
    spike_subsample_factor       = 2,
    spike_fit_samples_leftright  = 7,
    var_reduction_spike          = 2,
    var_reduction_total          = 80, 
    show_fit                     = False,
    store_glitches               = False,
    plot_removal_statistic       = False, 
    **kwargs):


    now = time.time()



    ### FIXED PARAMETERS
    GLITCH_WINDOW_LEFTRIGHT      = float(glitch_window_leftright)  # in s
    PREPEND_ZEROS                = float(glitch_prepend_zeros)     # in s, this is added additionally to 'GLITCH_WINDOW_LEFTRIGHT' as detected glitch onset is pretty accurate
    NUM_SAMPLES_GLITCH           = 10000   # must be > max_glitch_length_s * max_sampling period, so larger 10000 (= max_glitch_length_s=100 * max_sampling_period=100)
    CUT_SOURCE_BEFORE            = 20      # in samples for both glitch and spike, reason: acausal FIR-Filters)

    ACC_STEP                     = 1e-9    # in m/s**2, for glitches
    TRY_DECAY                    = False   # Source time function maniuplation. Here, osscilating decay after steo simulating some swinging og of interal structure
    TRY_INTER                    = int(glitch_interpolation_samples) # Source time function maniuplation. Here, interpolation samples to deviate from instantaneous step to finite ramp. Will be samples every 10 samples. Must be samller than 2*NUM_SAMPLES
    GAUSS_SOURCE_INTER           = False   # If True and 'TRY_INTER' != 0, the interpolation samples of the non-zero acc. rise will not be following a straight line but the first half of a Gaussian
    SUBSAMPLE_FACTOR_GLITCH      = int(glitch_subsample_factor)

    DIS_STEP                     = 1e-12   # in m, for glitch spikes
    SPIKE_SHIFT_SAMPLES_PER_SIDE = int(spike_fit_samples_leftright) # maximum of spike shifted left and right w.r.t. determined fitted glitch onset (should be larger than samples modeled glitch signal before real glitch starts)
    SUBSAMPLE_FACTOR_SPIKE       = int(spike_subsample_factor)
    VAR_REDUCTION_MIN_SPIKE      = float(var_reduction_spike)       # in %



    ### OUTPUT
    print()
    print(u'  ----------------------')
    print(u'  RUNNING GLITCH REMOVER')
    print(u'  ----------------------')



    ### FUNCTIONS TO MODEL GLITCHES AND SPIKES, DEPENDING ON SOME EXPERIMENTAL STUFF
    def residual_model(x, a,b,c,d,e):

        return a*x**4 + b*x**3 + c*x**2 + d*x + e
    def syn_glitch(interpolation_samples=0, amp=None, abkling=None, freq=None):

        source_FFT, step_freqs = sourceFFT(NUM_SAMPLES_GLITCH, sampling_period, int_samples=interpolation_samples, gauss=GAUSS_SOURCE_INTER, abkling=abkling, amp=amp, freq=freq)
        signal                 = fft2signal(source_FFT*ACC_fft, step_freqs, tau=tau, amp_factor=ACC_STEP)
        signal                 = signal[cut_index:cut_index+len(data_shifted)-prepend]
        signal                 = np.hstack(( np.zeros(prepend), signal )) 

        return signal
    def glitch_model_nodecay(x, m, n, o):
        """
        Three fit-variables.
        """
        return x * m + n + syn_glitch(interpolation_samples=interpolation_samples) * o
    def glitch_model_decay(x, m, n, o, amp, abkling, freq):
        """
        Six fit-variables.
        """
        return x * m + n + syn_glitch(abkling=abkling, amp=amp, freq=freq) * o
    def spike_model(x, m, n, o):
        """
        Three fit-variables.
        """
        return x * m + n + syn_spike * o



    ### READ GLITCH-FILE
    all_glitches = []

    for glitch_detector_file in glitch_detector_files:
        glitches      = np.loadtxt(glitch_detector_file, dtype='str')
        all_glitches += list(glitches)
    all_glitches = np.array(all_glitches)
    


    ### OTHER VARIABLES
    # for evaluation
    assign      = {'U':5, 'V':6, 'W':6}
    var_red_U   = []
    var_red_V   = []
    var_red_W   = []
    acc_steps_U = []
    acc_steps_V = []
    acc_steps_W = []



    ### READ WAVEFORM FILES
    for o, waveform_file in enumerate(waveform_files):
        removed              = []
        total_data_time_UTC  = 0
        total_data_time_LMST = 0

        # read streams
        try: 
            stream    = read2(waveform_file)
            format_in = stream[0].stats._format
        except TypeError:       # could not read file, skip
            continue   

        # small output
        print()        
        print(u'INFO: Analysing file: %s/%s' % (o+1,len(waveform_files)))
        print(waveform_file)
        for trace in stream:
            print('  %s' % trace)

        stream.set_inventory(inventory_file)

        # loop traces
        for trace in stream.select(channel='?[LMH]?'):

            print()
            trace.data = trace.data.astype(dtype=np.float32, copy=False)        # important as glitch removal introduces float RAW values
            total_data_time_UTC  += marstime(trace.stats.endtime).UTC_time  - marstime(trace.stats.starttime).UTC_time
            total_data_time_LMST += marstime(trace.stats.endtime).LMST_time - marstime(trace.stats.starttime).LMST_time

            # data prep
            component        = trace.stats.channel[-1]
            sampling_period  = trace.stats.delta
            prepend          = int(max(0,PREPEND_ZEROS)/sampling_period)
            glitches         = all_glitches[ (all_glitches[:,1]>=str(trace.stats.starttime)) & (all_glitches[:,2]<=str(trace.stats.endtime)) ]
            glitches         = glitches[np.argsort(glitches[:,0])]
            inv              = stream.inventory.select(network   = trace.stats.network, 
                                                       station   = trace.stats.station, 
                                                       location  = trace.stats.location, 
                                                       channel   = trace.stats.channel, 
                                                       starttime = trace.stats.starttime, 
                                                       endtime   = trace.stats.endtime)
            response         = inv[0][0][0].response
            optimum_interval = int(np.sqrt( (1+prepend+2*GLITCH_WINDOW_LEFTRIGHT/sampling_period)/(2*SUBSAMPLE_FACTOR_GLITCH) ))

            # output
            print(u'INFO: Handling %s glitches.' % len(glitches))

            # synthetic glitch+spike generation, for each new trace once
            step_source_FFT, source_freqs = sourceFFT(NUM_SAMPLES_GLITCH, sampling_period, int_samples=0)
            ACC_fft, _                    = responseFFT(response, source_freqs, step_unit='ACC')
            DIS_fft, _                    = responseFFT(response, source_freqs, step_unit='DIS')

            # looping over glitches to be corrected
            for g in range(len( glitches )):

                # glitch variables
                glitch_number  = glitches[g][0].replace(':','')
                glitch_start   = marstime(glitches[g][1])
                glitch_end     = marstime(glitches[g][2])
                try:
                    if marstime(glitches[g+1][1]).UTC_time<marstime(glitches[g][2]).UTC_time:
                        glitch_end = marstime(glitches[g+1][1])              # case of poly-glitches, fit window only until the start of next glitch
                except IndexError:
                    pass

                glitch_len = prepend+int((glitch_end.UTC_time-glitch_start.UTC_time)/sampling_period)
                #if not (component=='V' and glitch_number=='000162'):
                #    continue

                # variables needed
                fit_slice          = trace.slice(starttime=glitch_start.UTC_time-GLITCH_WINDOW_LEFTRIGHT-PREPEND_ZEROS, endtime=glitch_end.UTC_time+GLITCH_WINDOW_LEFTRIGHT)
                ori_slice          = fit_slice.copy()
                data_len_fit       = len(fit_slice.data)
                cor_slice          = trace.slice(starttime=glitch_start.UTC_time-GLITCH_WINDOW_LEFTRIGHT-PREPEND_ZEROS)
                uncor_slice        = cor_slice.copy()
                int_samples_array  = np.arange(0,TRY_INTER+1,10)



                ### G L I T C H
                # some lists to be filled
                residuals          = [[[np.nan  for _ in range(SUBSAMPLE_FACTOR_GLITCH)] for _ in range(data_len_fit-glitch_len)] for _ in range(len(int_samples_array))]
                fits               = [[[np.nan  for _ in range(SUBSAMPLE_FACTOR_GLITCH)] for _ in range(data_len_fit-glitch_len)] for _ in range(len(int_samples_array))]
                popts              = [[[np.nan  for _ in range(SUBSAMPLE_FACTOR_GLITCH)] for _ in range(data_len_fit-glitch_len)] for _ in range(len(int_samples_array))]
                                
                for k, interpolation_samples in enumerate(int_samples_array):
                    
                    # this bit fits the modeled glitch against the data in a coarse search and then determines index of residual minimum (using a polynomial 4th order)
                    tau              = 0
                    cut_index        = NUM_SAMPLES_GLITCH - max(CUT_SOURCE_BEFORE, int(interpolation_samples/2))
                    residuals_fit    = []
                    indices_fit      = []

                    for m in range(0,data_len_fit-glitch_len,optimum_interval):

                        # DATA
                        data_shifted = cor_slice.data[m:m+glitch_len]
                        
                        # FIT
                        p0           = [0, data_shifted[0], 1]
                        bounds       = ([-np.inf, np.min(data_shifted), -np.inf],[np.inf, np.max(data_shifted), np.inf])
                        popt, _      = scipy.optimize.curve_fit(glitch_model_nodecay, np.arange(len(data_shifted)), data_shifted, p0=p0, bounds=bounds)
                        fit          = glitch_model_nodecay(np.arange(len(data_shifted)), *popt)
                        residual     = np.linalg.norm(data_shifted - fit)
                        residuals_fit.append(residual)
                        indices_fit.append(m)

                    # fit of determined residuals
                    popt, _   = scipy.optimize.curve_fit(residual_model, np.arange(len(residuals_fit)), residuals_fit)
                    fit       = residual_model(np.arange(len(residuals_fit)), *popt)
                    index_min = indices_fit[np.argmin(fit)]

                    # around minimum of fitted residual, fit modeled glitch in finer grids against data to retrieve (including sub-sample shifts now)
                    for i in range(index_min-optimum_interval,index_min+optimum_interval+1,1):

                        if i < 0:
                            continue

                        if i >= data_len_fit-glitch_len:
                            break

                        # DATA
                        data_shifted = cor_slice.data[i:i+glitch_len]


                        for j in range(SUBSAMPLE_FACTOR_GLITCH):

                            # More signal prep
                            tau = sampling_period * j/SUBSAMPLE_FACTOR_GLITCH
                            #print(k,i,j)

                            # FITS: start values and their boundaries can improve fits, but not significantly whilst increasing calculation time significantly
                            if not TRY_DECAY:   # FIT variables: m, n, o ; but with fixed 'interpolation_samples' passed as global variable
                                p0      = [0, data_shifted[0], 1]
                                bounds  = ([-np.inf, np.min(data_shifted), -np.inf],[np.inf, np.max(data_shifted), np.inf])
                                popt, _ = scipy.optimize.curve_fit(glitch_model_nodecay, np.arange(len(data_shifted)), data_shifted, p0=p0, bounds=bounds)
                                fit     = glitch_model_nodecay(np.arange(len(data_shifted)), *popt)
                                
                            else:               # FIT variables: m, n, o, + amp, abkling, freq (for decay)
                                p0      = [0, data_shifted[0], 1, 0, 0, 4]
                                bounds  = ([-np.inf, np.min(data_shifted), -np.inf, 0, 0, 0],[np.inf, np.max(data_shifted), np.inf, 20, 0.1, 1/(2*sampling_period)])
                                popt, _ = scipy.optimize.curve_fit(glitch_model_decay, np.arange(len(data_shifted)), data_shifted, p0=p0, bounds=bounds)
                                fit     = glitch_model_decay(np.arange(len(data_shifted)), *popt)
                            
                            residual           = np.linalg.norm(data_shifted - fit)
                            residuals[k][i][j] = residual
                            popts[k][i][j]     = popt
                            fits[k][i][j]      = fit


                # best fit glitch
                residuals    = np.array(residuals)
                fits         = np.array(fits)
                popts        = np.array(popts)
                #quick_plot(*(residuals[i,:,0] for i in range(len(residuals)) ), data_labels=['# ramp inter. samples=%s' % i for i in int_samples_array])
                
                best_indices = np.where(residuals==np.nanmin(residuals))
                best_popt    = popts[best_indices][0]
                best_fit     = fits[best_indices][0]
                best_inter   = int_samples_array[best_indices[0][0]]           
                best_shift   = best_indices[1][0]
                best_tau     = best_indices[2][0] * sampling_period/SUBSAMPLE_FACTOR_GLITCH

                # from data, we do not subtact the fit which is as long as glitch stated in glitch detector file, but subtract a longer glitch of sample length ´NUM_SAMPLES_GLITCH´
                length_glitch = np.min([NUM_SAMPLES_GLITCH,len(cor_slice.data[best_shift:])])-prepend
                cut_index     = NUM_SAMPLES_GLITCH - max(CUT_SOURCE_BEFORE, int(best_inter/2))
                if TRY_INTER:
                    source_FFT, step_freqs = sourceFFT(NUM_SAMPLES_GLITCH, sampling_period, int_samples=best_inter, gauss=GAUSS_SOURCE_INTER)
                    print('int_samples=%s' % best_inter)
                elif TRY_DECAY:
                    source_FFT, step_freqs = sourceFFT(NUM_SAMPLES_GLITCH, sampling_period, abkling=best_popt[4], amp=best_popt[3], freq=best_popt[5])
                    print('amp=%s, decay=%s, freq=%s' % (best_popt[3], best_popt[4], best_popt[5]))
                else:
                    source_FFT, step_freqs = sourceFFT(NUM_SAMPLES_GLITCH, sampling_period, int_samples=0)

                scaled_glitch = fft2signal(source_FFT*ACC_fft, source_freqs, tau=best_tau, amp_factor=ACC_STEP*best_popt[2])
                scaled_glitch = scaled_glitch[cut_index:cut_index+length_glitch]   
                scaled_glitch = np.hstack(( np.zeros(prepend), scaled_glitch ))

                # actual glitch correction!
                cor_slice.data[best_shift:best_shift+len(scaled_glitch)] = cor_slice.data[best_shift:best_shift+len(scaled_glitch)]-scaled_glitch



                ### S P I K E
                # some lists to be filled
                residuals2         = []
                fits2              = []
                popts2             = []
                correction_indices = []   

                # looping over each data point for glitch spike correction
                tag_spike = False
                source_FFT, step_freqs = sourceFFT(NUM_SAMPLES_GLITCH, sampling_period)
                if spike_fit:

                    for i in range(2*SPIKE_SHIFT_SAMPLES_PER_SIDE+1+best_inter):

                        # DATA
                        start_index = i+best_shift+prepend-SPIKE_SHIFT_SAMPLES_PER_SIDE-int(best_inter/2)

                        if start_index<0:                           # case: best_shift+prepend+i<SPIKE_SHIFT_SAMPLES_PER_SIDE which creates negative index
                            start_index = 0
                            correction_index = abs(start_index)     # if that happens, cut accordant samples from model to fit data
                        else:
                            correction_index = 0
                        correction_indices.append(correction_index)
                        end_index    = start_index+glitch_len
                        data_shifted = fit_slice.data[start_index:end_index]

                        for j in range(SUBSAMPLE_FACTOR_SPIKE):

                            # More signal prep
                            tau        = sampling_period * j/SUBSAMPLE_FACTOR_SPIKE
                            syn_spike = fft2signal(step_source_FFT*DIS_fft, step_freqs, tau=tau, amp_factor=DIS_STEP)
                            syn_spike = syn_spike[cut_index+correction_index:cut_index+correction_index+len(data_shifted)]

                            # FIT with variables: m, n, o
                            # start values and their boundaries can improve fits, but not significantly whilst increasing calculation time significantly
                            #p0       = [0,data_shifted[0],1]
                            #bounds   = ([-1, np.min(data_shifted), -np.inf],[1, np.max(data_shifted), np.inf])
                            popt, _  = scipy.optimize.curve_fit(spike_model, np.arange(len(data_shifted)), data_shifted)
                            fit      = spike_model(np.arange(len(data_shifted)), *popt)
                            residual = np.linalg.norm(data_shifted - fit)

                            # storing needed results
                            residuals2.append(residual)
                            popts2.append(popt)
                            fits2.append(fit)

                    # best fit spike
                    best_index2   = np.array( residuals2 ).argmin()
                    best_popt2    = popts2[best_index2]
                    best_fit2     = fits2[best_index2]
                    best_shift2   = max(0, best_index2//SUBSAMPLE_FACTOR_SPIKE+best_shift+prepend-SPIKE_SHIFT_SAMPLES_PER_SIDE-int(best_inter/2))
                    best_cindex2  = correction_indices[best_index2//SUBSAMPLE_FACTOR_SPIKE]            # cutting  of slightly in case not enough data in beginning
                    best_tau2     = (best_index2%SUBSAMPLE_FACTOR_SPIKE) * sampling_period/SUBSAMPLE_FACTOR_SPIKE
                    #print(best_shift,best_index2,best_shift2,best_tau2,data_len_fit,len(best_fit2),best_cindex2)

                    # from data, we do not subtact the fit which is as long as glitch stated in glitch detector file, but subtract a longer glitch of sample length ´NUM_SAMPLES_GLITCH´
                    length_spike = np.min([NUM_SAMPLES_GLITCH,len(cor_slice.data[best_shift2:])])
                    scaled_spike = fft2signal(source_FFT*DIS_fft, step_freqs, tau=best_tau2, amp_factor=DIS_STEP*best_popt2[2])
                    scaled_spike = scaled_spike[cut_index+best_cindex2:cut_index+length_spike]

                    # variance reduction between deglitched and deglitched+despikeed data
                    var_data           = np.var(fit_slice)
                    deg_slice          = cor_slice.copy()
                    cor_slice.data[best_shift2:best_shift2+len(scaled_spike)] = cor_slice.data[best_shift2:best_shift2+len(scaled_spike)]-scaled_spike      # actual spike correction!           
                    var_data_corrected = np.var(fit_slice)
                    var_red_spi        = (1 - var_data_corrected/var_data) * 100

                    if var_red_spi>100:
                        var_red_spi = 100
                    if var_red_spi<0:
                        var_red_spi = 0

                    # spike correction
                    if var_red_spi>=VAR_REDUCTION_MIN_SPIKE:
                        tag_spike = True
                        time_diff_GP  = (best_shift - best_shift2 + prepend) * sampling_period + best_tau - best_tau2
                    else:
                        tag_spike = False
                        cor_slice.data[:] = deg_slice.data[:]           # undo spike correction as not good enough

                # variance reduction between original and deglitched or deglitched+despikeed data
                var_data           = np.var(ori_slice)
                var_data_corrected = np.var(fit_slice)
                var_red            = (1 - var_data_corrected/var_data) * 100
                if var_red>100:
                    var_red = 100
                if var_red<0:
                    var_red = 0

                # glitch correction undone or not
                if var_red>=var_reduction_total:
                    removed.append(component)
                    if spike_fit:
                        if tag_spike:
                            print(u'Glitch %6s,  %s,  %s,  var_red=%4.1f %%.  Correction is done.  Acc. step = %6.1f nm/s**2.  Dis. step = %6.1f pm (var_red_spi=%4.1f %%, Tgli-Tspi=%6.3f s)'  % (glitch_number, glitch_start.UTC_string, component, var_red, best_popt[2], best_popt2[2], var_red_spi, time_diff_GP))
                            label = 'glitch+spike corrected'
                        else:
                            print(u'Glitch %6s,  %s,  %s,  var_red=%4.1f %%.  Correction is done.  Acc. step = %6.1f nm/s**2'  % (glitch_number, glitch_start.UTC_string, component, var_red, best_popt[2]))
                            label = 'glitch corrected'
                    else:
                        print(u'Glitch %6s,  %s,  %s,  var_red=%4.1f %%.  Correction is done.  Acc. step = %6.1f nm/s**2'  % (glitch_number, glitch_start.UTC_string, component, var_red, best_popt[2]))
                        label = 'glitch corrected'
                else:

                    print(u'Glitch %6s,  %s,  %s,  var_red=%4.1f %%.  No correction done.'  % (glitch_number, glitch_start.UTC_string, component, var_red))
                    cor_slice.data[:] = uncor_slice.data[:]           # undo glitch correction as not good enough
                    label = 'not corrected'

                # only to convey information to user
                if show_fit:
                    glitch_slice = fit_slice.copy()
                    glitch       = np.hstack(( np.nan*np.ones(best_shift), best_fit, np.nan*np.ones(data_len_fit-len(best_fit)-best_shift) ))                   # glitch
                    if spike_fit and tag_spike:
                        spike = np.hstack(( np.nan*np.ones(best_shift2), scaled_spike[:len(best_fit2)], np.nan*np.ones(data_len_fit-len(best_fit2)-best_shift2) ))   # spike
                    else:
                        spike = np.nan*np.ones(data_len_fit)

                    data_plot = []
                    for i in range(len(glitch)):
                        if np.isnan(glitch[i]) and np.isnan(spike[i]):
                            value = np.nan
                        elif not np.isnan(glitch[i]) and np.isnan(spike[i]):
                            value = glitch[i]
                        else:
                            value = glitch[i] + spike[i]
                        data_plot.append( value )

                    glitch_slice.data  = np.array( data_plot )
                    title              = 'Glitch %s on component %s (var_red=%.1f%%, %s)' % (glitch_number, component, var_red, label)
                    data_labels        = ('original', 'corrected', 'fit (on original)')
                    quick_plot(ori_slice, fit_slice, glitch_slice, title=title, xlabel='Time UTC (s)', ylabel='Raw Units', lw=[1,1,1.5], lc=[None, None, 'k'], data_labels=data_labels)

                # results
                if component.upper()=='U':
                    if var_red>=var_reduction_total or glitches[g][5]=='1':
                        var_red_U.append( var_red )
                        acc_steps_U.append( best_popt[2]*ACC_STEP )
                elif component.upper()=='V':
                    if var_red>=var_reduction_total or glitches[g][6]=='1':
                        var_red_V.append( var_red )
                        acc_steps_V.append( best_popt[2]*ACC_STEP )
                elif component.upper()=='W':
                    if var_red>=var_reduction_total or glitches[g][7]=='1':
                        var_red_W.append( var_red )
                        acc_steps_W.append( best_popt[2]*ACC_STEP )


        ## WRITE DEGLITCHED FILE
        outfile_degl = '.'.join(waveform_file.split('.')[:-1]) + '_deglitched.' + format_in.lower()
        stream.write(outfile_degl, format=format_in)
        if store_glitches:                              # write out only fitted glitches as time series
            stream_orig = read2(waveform_file)
            for tr_orig, tr_degl in zip(stream_orig, stream):
                tr_orig.data = tr_orig.data - tr_degl.data
            outfile_glit = '.'.join(waveform_file.split('.')[:-1]) + '_glitches.' + format_in.lower()
            stream_orig.write(outfile_glit, format=format_in)


        ## OUTPUT
        string  = ''
        removed = np.array(removed)
        for comp in sorted(set(removed)):
            string += '%s:%s, ' % (comp, len(removed[removed==comp]))
        string = '(' + string[:-2] + ')'

        output    = [u'GLITCH DECTECTOR FILE(S):',
                     u'',      
                     u'DATA:',
                     u'  path:                         %s'           % waveform_file,
                     u'  IDs:                          %s'           % ', '.join( sorted(set(stream._get_ids())) ),
                     u'  UTC start:                    %s'           % marstime(stream.times[0]).UTC_string,
                     u'  UTC end:                      %s'           % marstime(stream.times[1]).UTC_string,
                     u'  UTC length (without gaps):    %s (h:m:s)'   % sec2hms( total_data_time_UTC ),
                     u'  LMST start:                   %s'           % marstime(stream.times[0]).LMST_string,
                     u'  LMST end:                     %s'           % marstime(stream.times[1]).LMST_string,
                     u'  LMST length (without gaps):   %s (h:m:s)'   % sec2hms( total_data_time_LMST ).replace('days', 'sols').replace('day','sol'),
                     u'',      
                     u'PARAMETERS:',      
                     u'  Glitch window left/right:     %s s'         % glitch_window_leftright,
                     u'  Glitch prepend zeros:         %s s'         % glitch_prepend_zeros,
                     u'  Glitch interpolation samples: %s'           % glitch_interpolation_samples,
                     u'  Glitch subsample factor:      %s'           % glitch_subsample_factor,
                     u'  Spike fit                     %s'           % spike_fit,
                     u'  Spike subsample factor:       %s m/s**3'    % spike_subsample_factor,
                     u"  Spike_samples left/right:     %s"           % spike_fit_samples_leftright,
                     u"  Variance reduction spike:     %s %%"        % var_reduction_spike,
                     u"  Variance reduction total:     %s %%"        % var_reduction_total,
                     u'  Show fit:                     %s'           % show_fit,
                     u"  Store glitches:               %s"           % store_glitches,
                     u'  Plot removal statistic:       %s'           % plot_removal_statistic,
                     u'',      
                     u'RESULTS:',      
                     u'    Removed a total of %s individual glitches on all traces.' % removed.size,
                     u'    %s' % string,
                     u'']

        # create good-looking output for later
        for f, file in enumerate(glitch_detector_files):
            output.insert(1+f, u'  %s' % file)
        if removed.size:                                # at least one glitch removed 
            output.insert(len(output), u'Done in:   %s (h:m:s), %.1f s per glitch' % (sec2hms( time.time()-now ), (time.time()-now)/removed.size ))
        else:
            output.insert(len(output), u'Done in:   %s (h:m:s)'                    %  sec2hms( time.time()-now ))
        output.insert(len(output), u'Timestamp: %s'   % datetime.datetime.now().strftime('%Y-%m-%dT%H:%M:%S'))

        len_max_line  = max( [len(line) for line in output] )
        formatter_str = "\n#   | %-" + str(len_max_line) + 's |'

        output = [formatter_str % line for line in output]
        output.insert(0,           '\n\n#   +-' +  '-' * len_max_line + '-+'   )
        output.insert(len(output),   '\n#   +-' +  '-' * len_max_line + '-+\n' )


        ### FINAL OUTPUT TO SHELL
        print()
        print(u'OUTPUT GLITCH REMOVER:')
        for line in output:
            print(line.strip('\n'))
        print()

        print(u'INPUT FILE:')
        print(waveform_file)
        print()
        print(u'DEGLITCHED FILE:')
        print(outfile_degl)
        print()

        if store_glitches:
            print(u'GLITCH FILE:')
            print(outfile_glit)


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

    # Variables
    PLOT_TRACES = False
    PLOT_LENGTH = 100    # in s

    files = ['Example/XB.ELYSE.03.BH?_2019-05-19T21:20:24.619000Z-2019-05-27T01:58:22.782000Z_raw.mseed']
    traces_glitch = []
    traces_spike  = []

    # Data prep
    for file in files:
        stream = read2(file, headonly=True)
        stream.set_inventory('IRIS')
        response        = stream.inventory[0][0][0].response
        sampling_period = stream[0].stats.delta
        header          = {'delta'    : sampling_period, 
                           'network'  : stream[0].stats.network, 
                           'station'  : stream[0].stats.station, 
                           'location' : stream[0].stats.location, 
                           'channel'  : stream[0].stats.channel}
        cut_index = 10000 - 20    # fixed. '20' samples before modeled glitch starts model is taken and fit against data (for both glitch and spike)


        source_FFT, step_freqs = sourceFFT(10000, sampling_period)
        ACC_fft, _ = responseFFT(response, step_freqs, step_unit='ACC')
        syn_glitch = fft2signal(source_FFT*ACC_fft, step_freqs, tau=0, amp_factor=1e-9)
        syn_glitch = syn_glitch[cut_index:cut_index+int(PLOT_LENGTH/sampling_period)]
        
        DIS_fft, _ = responseFFT(response, step_freqs, step_unit='DIS')
        syn_spike  = fft2signal(source_FFT*DIS_fft, step_freqs, tau=0, amp_factor=1e-12)
        syn_spike  = syn_spike[cut_index:cut_index+int(PLOT_LENGTH/sampling_period)]
        
        glitch     = Trace(data=syn_glitch, header=header)
        spike      = Trace(data=syn_spike, header=header)
        traces_glitch.append(glitch)
        traces_spike.append(spike)

    stream_glitch = Stream2(traces=traces_glitch)
    stream_spike  = Stream2(traces=traces_spike)

    # Plotting
    fig, axes = plt.subplots(2, figsize=(10,15), sharex=True)
    title     = 'VBB & SP modeled glitches and their spikes'
    fig.suptitle(title, fontsize=13)
    fig.align_ylabels()
    fig.subplots_adjust(hspace=0.2)
    axes[0].xaxis.set_ticks_position('none')


    if not PLOT_TRACES:
        axes[0] = quick_plot(*[trace.data for trace in stream_glitch],
                        axis=axes[0],
                        title='Glitches (step: 1e-9 m/s**2)',
                        ylabel='Digital Units',
                        xlabel=None,
                        data_labels=['VBB 10 SPS'],
                        legend_loc='upper right',
                        show=False)
        axes[1] = quick_plot(*[trace.data for trace in stream_spike],
                        axis=axes[1],
                        title='Precursors (step: 1e-12 m)',
                        ylabel='Digital Units',
                        xlabel='Data points',
                        data_labels=['VBB 10 SPS'],
                        legend_loc='upper right',
                        show=False)
    else:
        axes[0] = quick_plot(*stream_glitch,
                        axis=axes[0],
                        title='Glitches (step: 1e-9 m/s**2)',
                        ylabel='Digital Units',
                        xlabel=None,
                        legend_loc='upper right',
                        show=False)
        axes[1] = quick_plot(*stream_spike,
                        axis=axes[1],
                        title='Spikes (step: 1e-12 m)',
                        ylabel='Digital Units',
                        xlabel='Relative Time (s)',
                        legend_loc='upper right',
                        show=False)
    plt.show()