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
def stepFFT(sampling_period, num_samples=10000, decay=False, a=0.1, freq=0.01, abkling=0.01):

    # Step
    step        = np.hstack(( np.zeros(num_samples), 
                              np.ones(num_samples),
                             -np.ones(num_samples),
                              np.zeros(num_samples) ))
    if decay:
        t     = np.arange(num_samples)
        omega = np.pi*2*freq
        step[:num_samples]              += a * np.exp(-t*sampling_period*abkling) * np.cos(2*np.pi*freq*sampling_period * t)
        step[num_samples:2*num_samples] += a * np.exp(-t*sampling_period*abkling) * np.cos(2*np.pi*freq*sampling_period * t)

    # FFT
    step_fft    = np.fft.fft(step)
    step_freqs  = np.fft.fftfreq(step.size, d=sampling_period)

    # FFT multiplication
    return step_fft, step_freqs
def responseFFT(component_response, freqs, step_unit='ACC'):

    # FFT depending on 'step_unit'
    if step_unit.upper()=='ACC':
        # glitch
        resp_fft = component_response.get_evalresp_response_for_frequencies(freqs, output='ACC')
    elif step_unit.upper()=='VEL':
        # no terminology defined
        resp_fft = component_response.get_evalresp_response_for_frequencies(freqs, output='VEL')
    elif step_unit.upper()=='DIS' or step_unit.upper()=='DISP':
        # precursor
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
def _cut_index(component_response=None, offset=None):

    # Cut only away zeroes in front
    if not component_response:

        peak2peak_signal = np.max(signal)-np.min(signal)
        first_index      = np.where(signal>peak2peak_signal/1000.)[0][0]  # factor 1000 has influence on how many non-"zero" values are before actual glitch starts. For 2 SPS data these are more than for higher SPS. Thus influences precursor fits somewhat. 1000 is an okay trade-off.
        if first_index>0:
            first_index -= 1

    else:

        if offset is not None:
            pre_cut = 0
            for FIRstage in component_response.response_stages:
                
                if isinstance(FIRstage, obspy.core.inventory.response.FIRResponseStage):                        # inventory from IPGP
                    len_FIRcoeffs = len(FIRstage.coefficients)

                elif isinstance(FIRstage, obspy.core.inventory.response.CoefficientsTypeResponseStage):         # inventory from IRIS
                    len_FIRcoeffs = len(FIRstage.numerator)

                else:
                    continue

                decimation_factor = FIRstage.decimation_factor
                pre_cut           = (pre_cut+int(len_FIRcoeffs/2))/decimation_factor

            first_index = int(offset-pre_cut)
        else:
            print(u'WARNING: Cannot determine cut index without ´offset´. Set to 0.')
            first_index = 0

    return first_index

def remove(*glitch_detector_files, 
    waveform_files         = [], 
    inventory_file         = 'IRIS', 
    precursor_fit          = False,
    var_reduction_min      = 80, 
    show_glitch_fit        = False,
    store_glitches         = False,
    plot_removal_statistic = False, 
    **kwargs):


    now = time.time()



    ### OUTPUT
    print()
    print(u'  ----------------------')
    print(u'  RUNNING GLITCH REMOVER')
    print(u'  ----------------------')
    print()



    ### EQUATIONS TO MODEL GLTICHES AND PRECURSORS
    def glitch_model(x, m, n, o):
        """
        Three fit-variables.
        """
        return x * m + n + syn_glitch * o
    def precursor_model(x, m, n, o):
        """
        Three fit-variables.
        """
        return x * m + n + syn_precur * o



    ### FIXED PARAMETERS
    GLITCH_WINDOW_LEFTRIGHT          = 3       # in s
    PREPEND_ZEROS                    = 0       # in s
    NUM_SAMPLES_GLITCH               = 15000   # must be > max_glitch_length_s * max_sampling period, so larger max_glitch_length_s=100 * max_sampling_period=100 = 10000
    CUT_MODEL_BEFORE                 = 20      # in samples
    ACC_STEP                         = 1e-9    # in m/s**2, for glitches
    DIS_STEP                         = 1e-12   # in m, for glitch precursors
    SUBSAMPLE_FACTOR_GLITCH          = 5
    SUBSAMPLE_FACTOR_PRECUR          = 10
    PRECURSOR_SHIFT_SAMPLES_PER_SIDE = 5       # maximum of precursor shifted left and right w.r.t. determined fitted glitch onset (should be larger than samples modeled glitch signal before real glitch starts)
    VAR_REDUCTION_MIN_PRECURSOR      = 2       # in %



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
        print(u'INFO: Analysing file: %s/%s' % (o+1,len(waveform_files)))
        print(waveform_file)
        for trace in stream:
            print('  %s' % trace)

        stream.set_inventory(inventory_file)

        # loop traces
        for trace in stream.select(channel='?[LMH]?'):

            print()
            if store_glitches:
                trace_copy = trace.copy()

            # data prep
            component = trace.stats.channel[-1]
            sampling_period  = trace.stats.delta
            prepend          = int(PREPEND_ZEROS/sampling_period)
            glitches         = all_glitches[ (all_glitches[:,1]>=str(trace.stats.starttime)) & (all_glitches[:,2]<=str(trace.stats.endtime)) ]
            glitches         = glitches[np.argsort(glitches[:,0])]
            inv              = stream.inventory.select(network   = trace.stats.network, 
                                                       station   = trace.stats.station, 
                                                       location  = trace.stats.location, 
                                                       channel   = trace.stats.channel, 
                                                       starttime = trace.stats.starttime, 
                                                       endtime   = trace.stats.endtime)
            response         = inv[0][0][0].response
            #cut_index        = _cut_index(component_response=response, offset=NUM_SAMPLES_GLITCH)
            cut_index        = NUM_SAMPLES_GLITCH - CUT_MODEL_BEFORE      # fixed. '20' samples before modeled step starts, model is taken and fit against data (for both glitch and precursor, reason: acausal FIR-Filters)

            # ouput
            print(u'INFO: Handling %s glitches.' % len(glitches))

            # synthetic glitch generation, for each new trace once
            step_fft, step_freqs = stepFFT(sampling_period, num_samples=NUM_SAMPLES_GLITCH)
            glit_fft, _          = responseFFT(response, step_freqs, step_unit='ACC')
            prec_fft, _          = responseFFT(response, step_freqs, step_unit='DIS')

            # synthetic glitch precursor generation, for each new trace once
            if precursor_fit:
                prec_fft, _ = responseFFT(response, step_freqs, step_unit='DIS')

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
                #if component!='W':
                #    continue

                # variables needed
                fit_slice          = trace.slice(starttime=glitch_start.UTC_time-GLITCH_WINDOW_LEFTRIGHT-PREPEND_ZEROS, endtime=glitch_end.UTC_time+GLITCH_WINDOW_LEFTRIGHT)
                ori_slice          = fit_slice.copy()
                data_len_fit       = len(fit_slice.data)
                cor_slice          = trace.slice(starttime=glitch_start.UTC_time-GLITCH_WINDOW_LEFTRIGHT-PREPEND_ZEROS)
                uncor_slice        = cor_slice.copy()
                residuals          = []
                residuals2         = []
                fits               = []
                fits2              = []
                popts              = []
                popts2             = []
                correction_indices = []

                # looping over each data point for glitch correction
                for i in range(data_len_fit-glitch_len):

                    # DATA
                    data_shifted = cor_slice.data[i:i+glitch_len]

                    for j in range(SUBSAMPLE_FACTOR_GLITCH):

                        # More signal prep
                        tau        = sampling_period * j/SUBSAMPLE_FACTOR_GLITCH
                        syn_glitch = fft2signal(step_fft*glit_fft, step_freqs, tau=tau, amp_factor=ACC_STEP)
                        syn_glitch = syn_glitch[cut_index:cut_index+len(data_shifted)-prepend]
                        syn_glitch = np.hstack(( np.zeros(prepend), syn_glitch ))     

                        # FIT with variables: m, n, o
                        # start values and their boundaries can improve fits, but not significantly whilst increasing calculation time significantly                        
                        p0         = [0,data_shifted[0],1]
                        bounds     = ([-np.inf, np.min(data_shifted), -np.inf],[np.inf, np.max(data_shifted), np.inf])
                        popt, _    = scipy.optimize.curve_fit(glitch_model, np.arange(len(data_shifted)), data_shifted, p0=p0, bounds=bounds)
                        fit        = glitch_model(np.arange(len(data_shifted)), *popt)
                        residual   = np.linalg.norm(data_shifted - fit)

                        # storing needed results
                        residuals.append(residual)
                        popts.append(popt)
                        fits.append(fit)

                # best fit glitch
                best_index    = np.array( residuals ).argmin()
                best_popt     = popts[best_index]
                best_fit      = fits[best_index]
                best_shift    = best_index//SUBSAMPLE_FACTOR_GLITCH
                best_tau      = (best_index%SUBSAMPLE_FACTOR_GLITCH) * sampling_period/SUBSAMPLE_FACTOR_GLITCH
                #print(best_index,best_shift,best_tau)

                # from data, we do not subtact the fit which is as long as glitch stated in glitch detector file, but subtract a longer glitch of sample length ´NUM_SAMPLES_GLITCH´
                length_glitch = np.min([NUM_SAMPLES_GLITCH,len(cor_slice.data[best_shift:])])-prepend
                scaled_glitch = fft2signal(step_fft*glit_fft, step_freqs, tau=best_tau, amp_factor=ACC_STEP*best_popt[2])
                scaled_glitch = scaled_glitch[cut_index:cut_index+length_glitch]   
                scaled_glitch = np.hstack(( np.zeros(prepend), scaled_glitch ))

                # actual glitch correction!
                cor_slice.data[best_shift:best_shift+len(scaled_glitch)] = cor_slice.data[best_shift:best_shift+len(scaled_glitch)]-scaled_glitch

                # looping over each data point for glitch precursor correction
                tag_precursor = False
                if precursor_fit:

                    for i in range(2*PRECURSOR_SHIFT_SAMPLES_PER_SIDE+1):

                        # DATA
                        start_index  = best_shift+prepend-PRECURSOR_SHIFT_SAMPLES_PER_SIDE+i
                        if start_index<0:                                                       # case: best_shift+prepend+i<PRECURSOR_SHIFT_SAMPLES_PER_SIDE which creates negative index
                            start_index = 0
                            correction_index = abs(start_index)                                 # if that happens, cut accordant samples from model to fit data
                        else:
                            correction_index = 0
                        correction_indices.append(correction_index)
                        end_index    = start_index+glitch_len
                        data_shifted = fit_slice.data[start_index:end_index]

                        for j in range(SUBSAMPLE_FACTOR_PRECUR):

                            # More signal prep
                            tau        = sampling_period * j/SUBSAMPLE_FACTOR_PRECUR
                            syn_precur = fft2signal(step_fft*prec_fft, step_freqs, tau=tau, amp_factor=DIS_STEP)
                            syn_precur = syn_precur[cut_index+correction_index:cut_index+correction_index+len(data_shifted)]

                            # FIT with variables: m, n, o
                            # start values and their boundaries can improve fits, but not significantly whilst increasing calculation time significantly
                            #p0       = [0,data_shifted[0],1]
                            #bounds   = ([-1, np.min(data_shifted), -np.inf],[1, np.max(data_shifted), np.inf])
                            popt, _  = scipy.optimize.curve_fit(precursor_model, np.arange(len(data_shifted)), data_shifted)
                            fit      = precursor_model(np.arange(len(data_shifted)), *popt)
                            residual = np.linalg.norm(data_shifted - fit)

                            # storing needed results
                            residuals2.append(residual)
                            popts2.append(popt)
                            fits2.append(fit)

                    # best fit precursor                   
                    best_index2   = np.array( residuals2 ).argmin()
                    best_popt2    = popts2[best_index2]
                    best_fit2     = fits2[best_index2]
                    best_shift2   = best_index2//SUBSAMPLE_FACTOR_PRECUR+best_shift+prepend-PRECURSOR_SHIFT_SAMPLES_PER_SIDE
                    if best_shift2<0:
                        best_shift2 = 0
                    best_cindex2  = correction_indices[best_index2//SUBSAMPLE_FACTOR_PRECUR]
                    best_tau2     = (best_index2%SUBSAMPLE_FACTOR_PRECUR) * sampling_period/SUBSAMPLE_FACTOR_PRECUR
                    #print(best_index2,best_shift2,best_tau2)

                    # from data, we do not subtact the fit which is as long as glitch stated in glitch detector file, but subtract a longer glitch of sample length ´NUM_SAMPLES_GLITCH´
                    length_precur = np.min([NUM_SAMPLES_GLITCH,len(cor_slice.data[best_shift2:])])
                    scaled_precur = fft2signal(step_fft*prec_fft, step_freqs, tau=best_tau2, amp_factor=DIS_STEP*best_popt2[2])
                    scaled_precur = scaled_precur[cut_index+best_cindex2:cut_index+length_precur]

                    # variance reduction between deglitched and deglitched+deprecursored data
                    var_data           = np.var(fit_slice)
                    deg_slice          = cor_slice.copy()
                    cor_slice.data[best_shift2:best_shift2+len(scaled_precur)] = cor_slice.data[best_shift2:best_shift2+len(scaled_precur)]-scaled_precur      # actual precursor correction!           
                    var_data_corrected = np.var(fit_slice)
                    var_red_pre        = (1 - var_data_corrected/var_data) * 100

                    if var_red_pre>100:
                        var_red_pre = 100
                    if var_red_pre<0:
                        var_red_pre = 0

                    # precursor correction
                    if var_red_pre>=VAR_REDUCTION_MIN_PRECURSOR:
                        tag_precursor = True
                        time_diff_GP  = (best_shift - best_shift2) * sampling_period + best_tau - best_tau2 - prepend
                    else:
                        tag_precursor = False
                        cor_slice.data[:] = deg_slice.data[:]           # undo precursor correction as not good enough

                # variance reduction between original and deglitched or deglitched+deprecursored data
                var_data           = np.var(ori_slice)
                var_data_corrected = np.var(fit_slice)
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
                            print(u'Glitch %6s,  %s,  %s,  var_red=%4.1f %%,  tau=%.3f s.  Correction is done.  Acc. step = %6.1f nm/s**2.  Dis. step = %6.1f pm (var_red_pre=%4.1f %%, tau=%.3f s, Tgli-Tpre=%6.3f s)'  % (glitch_number, glitch_start.UTC_string, component, var_red, best_tau, best_popt[2], best_popt2[2], var_red_pre, best_tau2, time_diff_GP))
                            label = 'glitch+precursor corrected'
                        else:
                            print(u'Glitch %6s,  %s,  %s,  var_red=%4.1f %%,  tau=%.3f s.  Correction is done.  Acc. step = %6.1f nm/s**2'  % (glitch_number, glitch_start.UTC_string, component, var_red, best_tau, best_popt[2]))
                            label = 'glitch corrected'
                    else:
                        print(u'Glitch %6s,  %s,  %s,  var_red=%4.1f %%,  tau=%.3f s.  Correction is done.  Acc. step = %6.1f nm/s**2'  % (glitch_number, glitch_start.UTC_string, component, var_red, best_tau, best_popt[2]))
                        label = 'glitch corrected'
                else:
                    print(u'Glitch %6s,  %s,  %s,  var_red=%4.1f %%,  tau=%.3f s.  No correction done.'  % (glitch_number, glitch_start.UTC_string, component, var_red, best_tau))
                    cor_slice.data[:] = uncor_slice.data[:]           # undo glitch correction as not good enough
                    label = 'not corrected'

                # only to convey information to user
                if show_glitch_fit:
                    glitch_slice = fit_slice.copy()
                    glitch       = np.hstack(( np.nan*np.ones(best_shift), best_fit, np.nan*np.ones(data_len_fit-len(best_fit)-best_shift) ))          # glitch
                    if precursor_fit and tag_precursor:
                        precursor = np.hstack(( np.nan*np.ones(best_shift2), best_fit2-best_fit2[1], np.nan*np.ones(data_len_fit-len(best_fit2)-best_shift2) ))   # precursor
                    else:
                        precursor = np.nan*np.ones(data_len_fit)

                    data_plot = []
                    for i in range(len(glitch)):
                        if np.isnan(glitch[i]) and np.isnan(precursor[i]):
                            value = np.nan
                        elif not np.isnan(glitch[i]) and np.isnan(precursor[i]):
                            value = glitch[i]
                        else:
                            value = glitch[i] + precursor[i]
                        data_plot.append( value )

                    glitch_slice.data  = np.array( data_plot )
                    title              = 'Glitch %s on component %s (var_red=%.1f%%, %s)' % (glitch_number, component, var_red, label)
                    data_labels        = ('original', 'corrected', 'fit (on original)')
                    quick_plot(ori_slice, fit_slice, glitch_slice, title=title, xlabel='Time UTC (s)', ylabel='Raw Units', lw=[1,1,1.5], lc=[None, None, 'k'], data_labels=data_labels)

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

            # if True, assign trace data to only what was removed
            if store_glitches:
                trace.data[:] = trace_copy.data[:] - trace.data[:]


        ## WRITE DEGLITCHED FILE
        if store_glitches:
            outfile = '.'.join(waveform_file.split('.')[:-1]) + '_glitches.' + waveform_file.split('.')[-1]
        else:
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
        if store_glitches:
            print(u'GLITCH FILE:')
        else:
            print(u'DEGLITCHED FILE:')
        print(outfile)
        print()
        print(u'Done in:   %s (h:m:s), %.1f s per glitch per component.' % (sec2hms( time.time()-now ),(time.time()-now)/(len(stream)*len(glitches))))
        print(u'Timestamp: %s'                                           % datetime.datetime.now().strftime('%Y-%m-%dT%H:%M:%S'))


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


    return stream


### _ _ N A M E _ _ = = " _ _ M A I N _ _ "  
if __name__ == "__main__":

    # Variables
    PLOT_TRACES = True
    PLOT_LENGTH = 100    # in s

    files = ['/home/scholz/Desktop/data/XB.ELYSE.02.MH?_2019-03-15T06:28:06.779000Z-2019-03-15T08:32:53.279000Z_raw.MSEED',
             '/home/scholz/Desktop/data/XB.ELYSE.02.BH?_raw_S0325a_QB.mseed',
             '/home/scholz/Desktop/data/XB.ELYSE.00.HH?_sol423night.mseed',
             '/home/scholz/Desktop/data/XB.ELYSE.67.MH?_2019-03-08T03:56:39.154000Z-2019-03-08T04:30:59.154000Z_raw.mseed',
             '/home/scholz/Desktop/data/XB.ELYSE.67.SH?_2019-03-01T00:00:09.731000Z-2019-03-01T12:00:21.214000Z_raw.mseed',
             '/home/scholz/Desktop/data/XB.ELYSE.65.EH?_2019-03-08T03:59:59.314000Z-2019-03-08T04:30:01.944000Z_raw.mseed']
    traces_glitch = []
    traces_precur = []

    # Data prep
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
        #cut_index = _cut_index(component_response=response, offset=15000)
        cut_index = 15000 - 20          # fixed. '20' samples before modeled glitch starts model is taken and fit against data (for both glitch and precursor)


        step_fft, step_freqs = stepFFT(sampling_period, num_samples=15000)
        glit_fft, _          = responseFFT(response, step_freqs, step_unit='ACC')
        syn_glitch           = fft2signal(step_fft*glit_fft, step_freqs, tau=0, amp_factor=1e-9)
        syn_glitch           = syn_glitch[cut_index:cut_index+int(PLOT_LENGTH/sampling_period)]

        prec_fft, _          = responseFFT(response, step_freqs, step_unit='DIS')
        syn_precur           = fft2signal(step_fft*prec_fft, step_freqs, tau=0, amp_factor=1e-12)
        syn_precur           = syn_precur[cut_index:cut_index+int(PLOT_LENGTH/sampling_period)]

        glitch = Trace(data=syn_glitch, header=header)
        precur = Trace(data=syn_precur, header=header)
        traces_glitch.append(glitch)
        traces_precur.append(precur)

    stream_glitch = Stream2(traces=traces_glitch)
    stream_precur = Stream2(traces=traces_precur)
    print(stream_glitch)

    # Plotting
    fig, axes = plt.subplots(2, figsize=(10,15), sharex=True)
    title     = 'VBB & SP modeled glitches and their precursors'
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
                        data_labels=['VBB 2 SPS', 'VBB 20 SPS', 'VBB 100 SPS', 'SP 2 SPS', 'SP 20 SPS', 'SP 100 SPS'],
                        legend_loc='upper right',
                        show=False)
        axes[1] = quick_plot(*[trace.data for trace in stream_precur],
                        axis=axes[1],
                        title='Precursors (step: 1e-12 m)',
                        ylabel='Digital Units',
                        xlabel='Data points',
                        data_labels=['VBB 2 SPS', 'VBB 20 SPS', 'VBB 100 SPS', 'SP 2 SPS', 'SP 20 SPS', 'SP 100 SPS'],
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
        axes[1] = quick_plot(*stream_precur,
                        axis=axes[1],
                        title='Precursors (step: 1e-12 m)',
                        ylabel='Digital Units',
                        xlabel='Relative Time (s)',
                        legend_loc='upper right',
                        show=False)
    plt.show()