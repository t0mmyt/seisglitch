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
import sys
import time
import datetime
import numpy as np
from collections import Counter


#####  matplotlib modules import  #####
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as mticker
import matplotlib.cm as cmx


#####  obspy modules import  #####
from obspy import read, UTCDateTime


#####  seisglitch util import  #####
from seisglitch.util import read2, marstime, time_funnel, sec2hms, ppol, quick_plot, on_click
from seisglitch.math import normalise



### SELECT FUNCTION
def select_glitches(glitch_files,
    components          = '',
    sols                = [], 
    LMST_range          = ['0','24'], 
    Amplitudes          = {},
    AZs                 = [], 
    AZ_Errs             = [], 
    INCs                = [],
    INC_Errs            = [], 
    glitch_polarization = [],
    glitch_SNR          = [],
    inverse_selection   = False,
    verbose             = True,
    **kwargs):


    variables = locals().copy()
    del variables['glitch_files']
    del variables['kwargs']


    ### READ GLITCH-FILE
    all_glitches = []
    for glitch_file in glitch_files:
        glitches      = np.loadtxt(glitch_file, dtype=str)
        all_glitches += list(glitches)
    all_glitches = np.array(all_glitches)



    ### SELECT GLITCHES ON GIVEN `components`
    components = str(components).upper().strip()
    if components=='1':
        keep = [g for g in range(len(all_glitches)) if (all_glitches[g,5]=='1' and all_glitches[g,6]=='0' and all_glitches[g,7]=='0') or (all_glitches[g,5]=='0' and all_glitches[g,6]=='1' and all_glitches[g,7]=='0') or (all_glitches[g,5]=='0' and all_glitches[g,6]=='0' and all_glitches[g,7]=='1')]
    elif components=='2':
        keep = [g for g in range(len(all_glitches)) if (all_glitches[g,5]=='1' and all_glitches[g,6]=='1' and all_glitches[g,7]=='0') or (all_glitches[g,6]=='1' and all_glitches[g,6]=='0' and all_glitches[g,7]=='1') or (all_glitches[g,5]=='0' and all_glitches[g,6]=='1' and all_glitches[g,7]=='1')]
    elif components=='3':
        keep = [g for g in range(len(all_glitches)) if (all_glitches[g,5]=='1' and all_glitches[g,6]=='1' and all_glitches[g,7]=='1')]
    elif components=='U':
        keep = [g for g in range(len(all_glitches)) if all_glitches[g,5]=='1' and all_glitches[g,6]=='0' and all_glitches[g,7]=='0']
    elif components=='V':
        keep = [g for g in range(len(all_glitches)) if all_glitches[g,5]=='0' and all_glitches[g,6]=='1' and all_glitches[g,7]=='0']
    elif components=='W':
        keep = [g for g in range(len(all_glitches)) if all_glitches[g,5]=='0' and all_glitches[g,6]=='0' and all_glitches[g,7]=='1']
    elif sorted(components)=='UV':
        keep = [g for g in range(len(all_glitches)) if all_glitches[g,5]=='1' and all_glitches[g,6]=='1' and all_glitches[g,7]=='0']
    elif sorted(components)=='UW':
        keep = [g for g in range(len(all_glitches)) if all_glitches[g,5]=='1' and all_glitches[g,6]=='0' and all_glitches[g,7]=='1']
    elif sorted(components)=='VW':
        keep = [g for g in range(len(all_glitches)) if all_glitches[g,5]=='0' and all_glitches[g,6]=='1' and all_glitches[g,7]=='1']
    elif sorted(components)=='UVW':
        keep = [g for g in range(len(all_glitches)) if all_glitches[g,5]=='1' and all_glitches[g,6]=='1' and all_glitches[g,7]=='1']
    else:
        keep = [g for g in range(len(all_glitches))]



    ### SELECT GLITCHES ON GIVEN `SOLS`
    if sols:
        sols = sorted(sols)
        keep = np.array( [g for g in keep if marstime(all_glitches[g][1]).sol>=sols[0] and marstime(all_glitches[g][1]).sol<=sols[1]] )



    ### SELECT GLITCHES IN SPECIFIED `LMST RANGE`
    if LMST_range and len(LMST_range)==2:
        LMST_range        = [(str(t).replace(':','') + '000000')[:6] for t in LMST_range]
        LMST_range        = [[t[i:i+2] for i in range(0,len(t),2)] for t in LMST_range]
        LMST_range        = [':'.join(t) for t in LMST_range]
        glitch_start_LMST = np.array( [all_glitches[g,3].split('M')[1] for g in keep] )
        try:
            if LMST_range[0]<=LMST_range[1]:
                indices = np.where( (glitch_start_LMST >= LMST_range[0]) & (glitch_start_LMST <= LMST_range[1]) )
            else:
                indices = np.where( (glitch_start_LMST >= LMST_range[0]) | (glitch_start_LMST <= LMST_range[1]) )
            keep = list(np.array(keep)[indices[0]])
        except:
            pass


    ### SELECT GLITCHES ON GIVEN `Amplitudes`
    if Amplitudes:
        for comp in Amplitudes.keys():
            if comp.upper()=='U':
                column = 11
            elif comp.upper()=='V':
                column = 12
            elif comp.upper()=='W':
                column = 13
            elif comp.upper()=='Z':
                column = 14
            elif comp.upper()=='N':
                column = 15
            elif comp.upper()=='E':
                column = 16
            else:
                continue    

            if not isinstance(Amplitudes[comp], (list, tuple, np.ndarray)):
                Amplitudes[comp] = [Amplitudes[comp]]

            if len(Amplitudes[comp])==1:
                keep = [g for g in keep if float(all_glitches[g,column])<=float(Amplitudes[comp][0])]
            elif len(Amplitudes[comp])==2:
                keep = [g for g in keep if float(all_glitches[g,column])>=float(Amplitudes[comp][0]) and float(all_glitches[g,column])<=float(Amplitudes[comp][1])]



    ### SELECT GLITCHES ON GIVEN `AZs`
    if AZs:
        AZs = [float(AZ) for AZ in AZs]

        if len(AZs)==2:
            keep = [g for g in keep if float(all_glitches[g,17])>=AZs[0] and float(all_glitches[g,17])<=AZs[1]]

        elif len(AZs)==4:
            keep = [g for g in keep if (float(all_glitches[g,17])>=AZs[0] and float(all_glitches[g,17])<=AZs[1]) or 
                                       (float(all_glitches[g,17])>=AZs[2] and float(all_glitches[g,17])<=AZs[3])]

        elif len(AZs)==6:
            keep = [g for g in keep if (float(all_glitches[g,17])>=AZs[0] and float(all_glitches[g,17])<=AZs[1]) or 
                                       (float(all_glitches[g,17])>=AZs[2] and float(all_glitches[g,17])<=AZs[3]) or
                                       (float(all_glitches[g,17])>=AZs[4] and float(all_glitches[g,17])<=AZs[5])]

        elif len(AZs)==8:
            keep = [g for g in keep if (float(all_glitches[g,17])>=AZs[0] and float(all_glitches[g,17])<=AZs[1]) or 
                                       (float(all_glitches[g,17])>=AZs[2] and float(all_glitches[g,17])<=AZs[3]) or
                                       (float(all_glitches[g,17])>=AZs[4] and float(all_glitches[g,17])<=AZs[5]) or
                                       (float(all_glitches[g,17])>=AZs[6] and float(all_glitches[g,17])<=AZs[7])]

        elif len(AZs)==10:
            keep = [g for g in keep if (float(all_glitches[g,17])>=AZs[0] and float(all_glitches[g,17])<=AZs[1]) or 
                                       (float(all_glitches[g,17])>=AZs[2] and float(all_glitches[g,17])<=AZs[3]) or
                                       (float(all_glitches[g,17])>=AZs[4] and float(all_glitches[g,17])<=AZs[5]) or
                                       (float(all_glitches[g,17])>=AZs[6] and float(all_glitches[g,17])<=AZs[7]) or
                                       (float(all_glitches[g,17])>=AZs[8] and float(all_glitches[g,17])<=AZs[9])]

        elif len(AZs)==12:
            keep = [g for g in keep if (float(all_glitches[g,17])>=AZs[0] and float(all_glitches[g,17])<=AZs[1]) or 
                                       (float(all_glitches[g,17])>=AZs[2] and float(all_glitches[g,17])<=AZs[3]) or
                                       (float(all_glitches[g,17])>=AZs[4] and float(all_glitches[g,17])<=AZs[5]) or
                                       (float(all_glitches[g,17])>=AZs[6] and float(all_glitches[g,17])<=AZs[7]) or
                                       (float(all_glitches[g,17])>=AZs[8] and float(all_glitches[g,17])<=AZs[9]) or
                                       (float(all_glitches[g,17])>=AZs[10] and float(all_glitches[g,17])<=AZs[11])]   
        
        else:
            print(u'WARNING: AZs list must contain 2, 4, 6, 8, 10, or 12 elements. Given: %s' % AZs)



    ### SELECT GLITCHES ON GIVEN `AZ_Errs`
    if AZ_Errs:
        if len(AZ_Errs)==1:
            keep = [g for g in range(len(keep)) if float(all_glitches[g,18])<=float(AZ_Errs[0])]
        elif len(AZ_Errs)==2:
            keep = [g for g in range(len(keep)) if float(all_glitches[g,18])>=float(AZ_Errs[0]) and float(all_glitches[g,18])<=float(AZ_Errs[1])]



    ### SELECT GLITCHES ON GIVEN `INCs`
    if INCs:
        INCs = [float(INC) for INC in INCs]

        if len(INCs)==2:
            keep = [g for g in keep if float(all_glitches[g,19])>=INCs[0] and float(all_glitches[g,19])<=INCs[1]]

        elif len(INCs)==4:
            keep = [g for g in keep if (float(all_glitches[g,19])>=INCs[0] and float(all_glitches[g,19])<=INCs[1]) or 
                                       (float(all_glitches[g,19])>=INCs[2] and float(all_glitches[g,19])<=INCs[3])]

        elif len(INCs)==6:
            keep = [g for g in keep if (float(all_glitches[g,19])>=INCs[0] and float(all_glitches[g,19])<=INCs[1]) or 
                                       (float(all_glitches[g,19])>=INCs[2] and float(all_glitches[g,19])<=INCs[3]) or
                                       (float(all_glitches[g,19])>=INCs[4] and float(all_glitches[g,19])<=INCs[5])]

        else:
            print(u'WARNING: INCs list must contain 2, 4, or 6 elements. Given: %s' % INCs)



    ### SELECT GLITCHES ON GIVEN `INC_Errs`
    if INC_Errs:
        if len(INC_Errs)==1:
            keep = [g for g in range(len(keep)) if float(all_glitches[g,20])<=float(INC_Errs[0])]
        elif len(INC_Errs)==2:
            keep = [g for g in range(len(keep)) if float(all_glitches[g,20])>=float(INC_Errs[0]) and float(all_glitches[g,20])<=float(INC_Errs[1])]



    ### SELECT GLITCHES ON GIVEN `glitch_SNR`
    if glitch_SNR:
        if len(glitch_SNR)==1:
            keep = [g for g in range(len(keep)) if float(all_glitches[g,21])>=float(glitch_SNR[0])]
        elif len(glitch_SNR)==2:
            keep = [g for g in range(len(keep)) if float(all_glitches[g,21])>=float(glitch_SNR[0]) and float(all_glitches[g,21])<=float(glitch_SNR[1])]



    ### SELECT GLITCHES ON GIVEN `glitch_polarization`
    if glitch_polarization:
        if len(glitch_polarization)==1:
            keep = [g for g in range(len(keep)) if float(all_glitches[g,22])>=float(glitch_polarization[0])]
        elif len(glitch_polarization)==2:
            keep = [g for g in range(len(keep)) if float(all_glitches[g,22])>=float(glitch_polarization[0]) and float(all_glitches[g,22])<=float(glitch_polarization[1])]



    ### SELECT INVERSE SELCTION
    if not inverse_selection:
        new_glitches = all_glitches[keep]
    else:
        inverse      = [i for i in range(len(all_glitches)) if i not in keep]
        new_glitches = all_glitches[inverse]



    ### PRINTING VARIABLES:
    if verbose:
        print(u'SELECTED %s/%s GLITCHES ON THE FOLLOWING PARAMETERS:' % (len(new_glitches),len(all_glitches)))
        for key, value in sorted(variables.items(), key=lambda item: len(item[0])):
            print(u'  %25s = %s' % (key, value))
        print()

    if new_glitches.size == 0:
        print(u'WARNING: With given parameters no glitches were selected.')
        sys.exit()
    else:
        return new_glitches



### PLOT SCRIPTS
def plot_glitch_remover(*glitch_files, run=True, original_data=None, deglitch_data=None, window=None, starttime=None, endtime=None, show=True, outfile='', **kwargs):
    
    """
    Plot comparison of deglitched data to original data.
    Plot also detected glitches per component.

    All input traces will be plotted e.g., also if they are not seismic.
    """


    ### RUN OR NOT:
    if not run:
        return

    now = time.time()



    ### OUTPUT
    print()
    print(u'  ------------------------------')
    print(u'  RUNNING GLITCH REMOVER PLOTTER')
    print(u'  ------------------------------')
    print()



    ### CHECK IF WAVEFORM FILES PASSED
    try:
        original_wf = read2(original_data)
    except:
        print(u'ERROR: Cannot read `original_data` via ObsPy.')
        sys.exit()

    try:
        deglitch_wf = read2(deglitch_data)
    except:
        print(u'ERROR: Cannot read `deglitch_data` via ObsPy.')
        sys.exit()

    combined_wf = original_wf + deglitch_wf


    ### WINDOW LENGTH
    if not window:
        window_length = 25*3600
    else:
        window_length = float(window)*3600



    ### SELECT GLITCHES
    all_glitches = select_glitches(glitch_files, **kwargs)
    


    ### CUT STREAMS
    if starttime:
        starttime = time_funnel(starttime).UTC_time
        original_wf.trim(starttime=starttime)
        deglitch_wf.trim(starttime=starttime)
    if endtime:
        endtime = time_funnel(endtime).UTC_time
        original_wf.trim(endtime=endtime)
        deglitch_wf.trim(endtime=endtime)



    ### SANITY CHECK
    if not original_wf or not deglitch_wf:
        print(u'WARNING: Times %s, and %s not present in %s. No plotting done.' % (starttime, endtime, waveform_file))
        sys.exit()



    ### PLOT DATA + GLITCHES
    j = 0
    combined_ids = sorted( set( combined_wf._get_ids() ), key=combined_wf._get_ids().index)

    for stream_window in combined_wf.slide(window_length, window_length*0.95, offset=0, include_partial_windows=True, nearest_sample=True):

        j += 1
        starttime_window = min([tr.stats.starttime for tr in stream_window])
        endtime_window   = max([tr.stats.endtime   for tr in stream_window])
        original_window  = original_wf.slice(starttime=starttime_window, endtime=endtime_window)
        deglitch_window  = deglitch_wf.slice(starttime=starttime_window, endtime=endtime_window)

        # glitches in shown window
        indices = sorted( [i for i in range(len(all_glitches)) if all_glitches[i][1]>=str(starttime_window) and all_glitches[i][1]<=str(endtime_window)] )
        if indices:
            for i in indices:
                print('  '+'   '.join(all_glitches[i,[0,1,5,6,7]]))
            print()

        # make plot       
        def on_xlims_change(event_ax):
            ax1Xs = axes[-1].get_xticks()
            ax2Xs = ['S'+marstime(mdates.num2date(time)).LMST_string.replace('M','\n') for time in ax1Xs]
            axis_twin.set_xticks(ax1Xs)
            axis_twin.set_xticklabels(ax2Xs)
        fig, _    = plt.subplots(len(combined_ids), figsize=(14,8), sharex=True)
        axes      = fig.axes
        win_title = 'Glitch remover plot (window %s)' % j
        title     = '%s glitches' % len(indices)
        fig.canvas.set_window_title(win_title)
        fig.canvas.mpl_connect('button_press_event', lambda event: on_click(event, fig))

        for i, id in enumerate(combined_ids):

            orientation_code = stream_window.select(id=id)[0].stats.channel[-1].upper()
            if orientation_code=='U':
                glitches_comp = all_glitches[indices][all_glitches[indices][:,5]=='1']
            elif orientation_code=='V':
                glitches_comp = all_glitches[indices][all_glitches[indices][:,6]=='1']
            elif orientation_code=='W':
                glitches_comp = all_glitches[indices][all_glitches[indices][:,7]=='1']
            else:
                glitches_comp = all_glitches[indices]

            title   = '%s glitches' % len(glitches_comp)
            if i==len(combined_ids)-1:
                xlabel = 'Time (UTC)'
            else:
                xlabel = None

            # actual data plot
            axes[i] = quick_plot(*original_window.select(id=id), 
                                 *deglitch_window.select(id=id), 
                                 title       = title, 
                                 verts       = [glitches_comp[:,1]], 
                                 lc          = ['C0'] * len(original_window.select(id=id)) + ['C1'] * len(deglitch_window.select(id=id)), 
                                 xlabel      = xlabel, 
                                 data_labels = ['original: %s' % id] * len(original_window.select(id=id)) + ['deglitch: %s' % id] * len(deglitch_window.select(id=id)), 
                                 legend_loc  = 'upper right', 
                                 axis        = axes[i])

        # twiny axis
        axis_twin = axes[-1].twiny()
        axes[-1].callbacks.connect('xlim_changed', on_xlims_change)
        ax1Xs = axes[-1].get_xticks()
        ax2Xs = ['S'+marstime(mdates.num2date(time)).LMST_string.replace('M','\n') for time in ax1Xs]
        axis_twin.set_xticks(ax1Xs)
        axis_twin.set_xbound(axes[-1].get_xbound())
        axis_twin.set_xticklabels(ax2Xs)
        axis_twin.set_xlabel(u"Time (LMST)", fontsize=11)
        plt.tight_layout()

        if show:
            plt.show()

        if outfile:
            outfile = '.'.join( outfile.split('.')[:-1]) + '%s.' % j + outfile.split('.')[-1]
            plt.savefig(outfile)
            print(outfile)
def plot_glitch_detector(*glitch_files, run=True, waveform_files=[], window=None, starttime=None, endtime=None, show=True, outfile='', **kwargs):
    
    """
    Plot original data.
    Plot also detected glitches per component.

    Only seismic traces will be plotted.
    """


    ### RUN OR NOT:
    if not run:
        return

    now = time.time()



    ### OUTPUT
    print()
    print(u'  -------------------------------')
    print(u'  RUNNING GLITCH DETECTOR PLOTTER')
    print(u'  -------------------------------')
    print()



    ### CHECK IF WAVEFORM FILES PASSED
    if not waveform_files:
        print(u'ERROR: You need to specify `waveform_files`.')
        sys.exit()



    ### WINDOW LENGTH
    if not window:
        window_length = 25*3600
    else:
        window_length = float(window)*3600



    ### SELECT GLITCHES
    all_glitches = select_glitches(glitch_files, **kwargs)



    ### LOOP OVER ALL FILES PASSED
    for waveform_file in waveform_files:

        # read / select stream
        try: 
            st = read2(waveform_file)
        except TypeError:       # could not read file, skip
            continue        
        st_select = st.select(channel='?[LMH]?')

        # sanity check
        if not st_select:
            print(u'WARNING: No seismic components found in %s. No plotting done.' % waveform_file)
            continue

        # if times specified
        if starttime:
            starttime = time_funnel(starttime).UTC_time
            st_select.trim(starttime=starttime)
        if endtime:
            endtime = time_funnel(endtime).UTC_time
            st_select.trim(endtime=endtime)
        
        # another sanity check
        if not st_select:
            print(u'WARNING: Times %s, and %s not present in %s. No plotting done.' % (starttime, endtime, waveform_file))
            continue
        else:
            print(u'Plotting: %s' % waveform_file)

        # sliding
        j = 0
        stream_ids = sorted( set( st_select._get_ids() ), key=st_select._get_ids().index)

        for stream_plot in st_select.slide(window_length, window_length*0.95, offset=0, include_partial_windows=True, nearest_sample=True):

            j += 1
            starttime_window = min([tr.stats.starttime for tr in stream_plot])
            endtime_window   = max([tr.stats.endtime   for tr in stream_plot])

            # glitches
            indices = sorted( [i for i in range(len(all_glitches)) if all_glitches[i][1]>=str(starttime_window) and all_glitches[i][1]<=str(endtime_window)] )
            if indices:
                for i in indices:
                    print('  '+'   '.join(all_glitches[i,[0,1,5,6,7]]))
                print()

            # make plot       
            def on_xlims_change(event_ax):
                ax1Xs = axes[-1].get_xticks()
                ax2Xs = ['S'+marstime(mdates.num2date(time)).LMST_string.replace('M','\n') for time in ax1Xs]
                axis_twin.set_xticks(ax1Xs)
                axis_twin.set_xticklabels(ax2Xs)
            fig, _    = plt.subplots(len(stream_ids), figsize=(14,8), sharex=True)
            axes      = fig.axes
            win_title = 'Glitch remover plot (window %s)' % j
            title     = '%s glitches' % len(indices)
            fig.canvas.set_window_title(win_title)
            fig.canvas.mpl_connect('button_press_event', lambda event: on_click(event, fig))

            for i, id in enumerate(stream_ids):

                orientation_code = stream_plot.select(id=id)[0].stats.channel[-1].upper()
                if orientation_code=='U':
                    glitches_comp = all_glitches[indices][all_glitches[indices][:,5]=='1']
                elif orientation_code=='V':
                    glitches_comp = all_glitches[indices][all_glitches[indices][:,6]=='1']
                elif orientation_code=='W':
                    glitches_comp = all_glitches[indices][all_glitches[indices][:,7]=='1']
                else:
                    glitches_comp = np.array( [] )    

                title   = '%s glitches' % len(glitches_comp)
                if i==len(stream_ids)-1:
                    xlabel = 'Time (UTC)'
                else:
                    xlabel = None

                # actual data plot
                axes[i] = quick_plot(*stream_plot.select(id=id), 
                                     title      = title, 
                                     verts      = [glitches_comp[:,1]], 
                                     lc         = ['C%s' % i], 
                                     xlabel     = xlabel, 
                                     legend_loc = 'upper right', 
                                     axis       = axes[i])

            # twiny axis
            axis_twin = axes[-1].twiny()
            axes[-1].callbacks.connect('xlim_changed', on_xlims_change)
            ax1Xs = axes[-1].get_xticks()
            ax2Xs = ['S'+marstime(mdates.num2date(time)).LMST_string.replace('M','\n') for time in ax1Xs]
            axis_twin.set_xticks(ax1Xs)
            axis_twin.set_xbound(axes[-1].get_xbound())
            axis_twin.set_xticklabels(ax2Xs)
            axis_twin.set_xlabel(u"Time (LMST)", fontsize=11)
            plt.tight_layout()

            if show:
                plt.show()

            if outfile:
                outfile = '.'.join( outfile.split('.')[:-1]) + '%s.' % j + outfile.split('.')[-1]
                plt.savefig(outfile)
                print(outfile)               
def plot_glitch_overview(*glitch_files, run=True, waveform_files=[], outfile='', **kwargs):
    
    """
    Plot glitches, based on glitch file produced by function `glitch_detector()`.

    amp parameters apply to all components: UVWZNE
    min_amp: at least one component amplitude must be larger
    max_amp: all component amplitudes must be smaller

    Matplotlib colorscale:
    https://matplotlib.org/3.1.1/tutorials/colors/colormaps.html
    """



    ### RUN OR NOT:
    if not run:
        return



    ### HELPER FUNCTION
    picked_already = False
    def onpick(event):

        indices = event.ind

        # print all picked glitches
        nonlocal picked_already
        if not picked_already:
            picked_already = True
            header         = u'#Number                  UTC          LMST    U    V    W      U-AMP      V-AMP      W-AMP      Z-AMP      N-AMP      E-AMP        AZ     INC     SNR_3D   POL_3D'
            print(header)

        for index in indices:
            print(u'%7s  %s  %s  %3s  %3s  %3s  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g    %6.1f   %-6.1f    %6.1g   %6.3f' % (all_glitches[index,0],
                                                                                                              marstime(all_glitches[index,1]).UTC_string, 
                                                                                                              marstime(all_glitches[index,1]).LMST_string, 
                                                                                                              int(all_glitches[index,5]),
                                                                                                              int(all_glitches[index,6]),
                                                                                                              int(all_glitches[index,7]),
                                                                                                              float(all_glitches[index,11]),
                                                                                                              float(all_glitches[index,12]),
                                                                                                              float(all_glitches[index,13]),
                                                                                                              float(all_glitches[index,14]),
                                                                                                              float(all_glitches[index,15]),
                                                                                                              float(all_glitches[index,16]),
                                                                                                              float(all_glitches[index,17]),
                                                                                                              float(all_glitches[index,19]),
                                                                                                              float(all_glitches[index,21]),
                                                                                                              float(all_glitches[index,22]),
                                                                                                              ))
        print(u'- - -')

        # print first picked glitches
        glitch_start = all_glitches[indices[0],1]
        plot_glitch_ppol(glitch_start, waveform_files=waveform_files, glitch_length=marstime(all_glitches[index,2]).UTC_time-marstime(all_glitches[index,1]).UTC_time, **kwargs)



    ### OUTPUT
    print()
    print(u'  -------------------------------')
    print(u'  RUNNING GLITCH OVERVIEW PLOTTER')
    print(u'  -------------------------------')
    print()



    ### SELECT GLITCHES
    all_glitches = select_glitches(glitch_files, **kwargs)
    LMST_range   = kwargs['LMST_range']



    ### ASSIGN NEEDED VARIABLES
    sols_range   = sorted( [marstime(e).sol for e in all_glitches[:,1]] )
    print(u'Detected sols: %s .. %s ' % (sols_range[0], sols_range[-1]))

    ymax_hist    = Counter( sols_range ).most_common()[0][1]

    try:
        glitch_starts     = np.array( [marstime(e) for e in all_glitches[:,1]] )

        Z_amps            = all_glitches[:,14].astype('float')
        N_amps            = all_glitches[:,15].astype('float')
        E_amps            = all_glitches[:,16].astype('float')
        phis              = all_glitches[:,17].astype('float') * np.pi/180
        incs              = all_glitches[:,19].astype('float') * np.pi/180
        sols              = np.array( [e.sol for e in glitch_starts] )
        

    except IndexError:  #no glitches for chosen conditions
        glitch_starts     = np.array( [] )

        Z_amps            = np.array( [] )
        N_amps            = np.array( [] )
        E_amps            = np.array( [] )
        phis              = np.array( [] )
        incs              = np.array( [] )
        sols              = np.array( [] )



    ### FIGURE 1
    fig1 = plt.figure(figsize=(15,7))
    fig1.canvas.set_window_title('Glitch overview plot 1')
    fig1.suptitle('Overview plotter: Sols=%s..%s, %s glitches' % (sols_range[0], sols_range[-1], len(all_glitches)), fontsize=12, y=0.99)
    fig1.subplots_adjust(wspace=0.4, hspace=0.4)
    fig1.canvas.mpl_connect('pick_event', onpick)
    if outfile:
        outfile1 = '.'.join( outfile.split('.')[:-1]) + '1.'+outfile.split('.')[-1]
        outfile2 = '.'.join( outfile.split('.')[:-1]) + '2.'+outfile.split('.')[-1]


    ## COLORMAP
    cmap      = plt.get_cmap('hsv')
    #norm      = Normalize(vmin=0, vmax=24)
    norm      = mpl.colors.BoundaryNorm(np.linspace(0,24,25), cmap.N)
    scalarMap = cmx.ScalarMappable(norm=norm, cmap=cmap)
    colours   = [scalarMap.to_rgba(glitch_starts[i].LMST_time.hour) for i in range(len(glitch_starts))]

    cax1 = fig1.add_axes([0.01, 0.55, 0.01, 0.4]) # left, bottom, width, height
    cax1.yaxis.set_ticks_position('left')
    cb1 = mpl.colorbar.ColorbarBase(cax1, drawedges=True, cmap=cmap, norm=norm, orientation='vertical', ticks=np.linspace(0,24,9), boundaries=np.linspace(0,24,25))
    cb1.set_label('LMST hour')

    # create hist data and colours
    glitch_start_hist_SOL_mat = [[] for _ in range(24)]
    for start in glitch_starts:
        glitch_start_hist_SOL_mat[start.LMST_time.hour] += [start.sol]
    #print(np.array(glitch_start_hist_SOL_mat)[0])
    
    if np.array(glitch_start_hist_SOL_mat).any():
        colors_hist = [scalarMap.to_rgba(i) for i in range(len(glitch_start_hist_SOL_mat))]
    else:   # if no glitches, make it still work for plotting
        colors_hist = [(0,0,0,0)]



    ## HORIZONTAL POLARIZATIONS (Fig 1, left)
    sizes = np.sqrt( N_amps**2 + E_amps**2 )
    sizes = normalise(sizes, scale_to_between=[5,20])

    ax   = fig1.add_subplot(121, polar=True)
    ax.set_title('Back-azimuths', pad=20, size=11)
    ax.spines['polar'].set_visible(False)
    ax.tick_params(pad=5)

    ax.set_theta_direction(-1)
    ax.set_theta_zero_location('N')
    ax.set_thetagrids([0,15,37,70,90,105,114,135,157,180,251,255,270,277,345.5], labels=['','LSA/VBBV\n(15°)','LVL1\n(37°)', r'HP$^{3}$'+'\n(70°)','','SP2\n(105°)','WTS E\n(114°)',r'VBBU'+'\n(135°)','LVL2\n(157°)','','\nWTS W\n(251°)','VBBW\n(255°)','','LVL3\n(277°)','SP3/WTS N\n(345°)'], size=10)

    #ax.set_rgrids([0.1,0.4,0.7,1], labels=['Hmin','','','Hmax'], size=6)
    ax.set_ylim( [min(sols_range)-0.1*(max(sols_range)-min(sols_range)), max(sols_range)+0.1*(max(sols_range)-min(sols_range))] )
    ax.tick_params(axis='y', labelsize=6)   # because with thetamin/max command, fontzise in set_rgrids doesn't work ...
    ax.set_yticklabels(['Sol %d' % x for x in ax.get_yticks()])
    ax.set_rlabel_position(315)

    ax.scatter(phis, sols, s=sizes, linewidths=.3, c=colours, rasterized=False, picker=1)
    ax.grid(ls='-.', lw=0.4, c='dimgray')



    ## VERTICAL POLARIZATIONS (Fig 1, right)
    sizes = np.sqrt(Z_amps**2)
    sizes = normalise(sizes , scale_to_between=[5,20])

    ax = fig1.add_subplot(122, polar=True)
    ax.set_title('Incidence angles', pad=20, size=11)
    ax.spines['polar'].set_visible(False)
    ax.tick_params(pad=5)

    ax.set_theta_direction(-1)
    ax.set_theta_zero_location('N')
    ax.set_thetamin(0)
    ax.set_thetamax(180)
    ax.set_thetagrids([0, 48.5, 60.5, 90, 131.5, 119.5, 180], labels=['0°','48.5°','VBB\n(60.5°)','90°','131.5°','VBB\n(119.5°)', '180°'], size=10)

    #ax.set_rgrids([0.1,0.4,0.7,1], labels=['Zmin','','','Zmax'], fontsize=6)
    ax.set_ylim( [min(sols_range)-0.1*(max(sols_range)-min(sols_range)), max(sols_range)+0.1*(max(sols_range)-min(sols_range))] )
    ax.tick_params(axis='y',labelsize=6)    # because with thetamin/max command, fontzise in set_rgrids doesn't work ...
    ax.set_yticklabels(['Sol %d' % x for x in ax.get_yticks()])
    ax.set_rlabel_position(315)

    ax.scatter(incs, sols, s=sizes, color=colours, linewidths=.3, rasterized=False, picker=1)
    ax.grid(ls='-.', lw=0.4, c='dimgray')

    if outfile:
        plt.savefig(outfile1)
        print(outfile1)


    ### FIGURE 2
    fig2 = plt.figure(figsize=(15,7))
    fig2.canvas.set_window_title('Glitch overview plot 2')
    fig2.suptitle('Overview plotter: Sols=%s..%s, %s glitches' % (sols_range[0], sols_range[-1], len(all_glitches)), fontsize=12, y=0.99)
    fig2.subplots_adjust(wspace=0.4, hspace=0.4)
    fig2.canvas.mpl_connect('pick_event', onpick)


    ## COLORMAP
    cax2 = fig2.add_axes([0.01, 0.55, 0.01, 0.4]) # left, bottom, width, height
    cax2.yaxis.set_ticks_position('left')
    cb2 = mpl.colorbar.ColorbarBase(cax2, drawedges=True, cmap=cmap, norm=norm, orientation='vertical', ticks=np.linspace(0,24,9), boundaries=np.linspace(0,24,25))
    cb2.set_label('LMST hour')


    ## 3-D PLOT (Fig 2, left)
    sizes = np.sqrt(Z_amps**2+N_amps**2+E_amps**2) 
    sizes = normalise(sizes, scale_to_between=[5,20])

    ax = fig2.add_subplot(121, projection='3d')
    ax.set_title('3-D Polarizations', y=1.10, size=11)
    ax.view_init(elev=15., azim=-50)
    ax.set_xlabel('E', labelpad=2, fontsize=8)
    ax.set_ylabel('N', labelpad=2, fontsize=8)
    ax.set_zlabel('Z', labelpad=2, fontsize=8)
    ax.set_xticks([-1,-0.5,0,0.5,1])
    ax.set_yticks([-1,-0.5,0,0.5,1])
    ax.set_zticks([-1,-0.5,0,0.5,1])
    ax.scatter(np.sin(phis)*abs(np.sin(incs)), np.cos(phis)*abs(np.sin(incs)), np.cos(incs), s=sizes, color=colours, linewidths=.3, depthshade=False, rasterized=True)
    ax.scatter(0, 0, 0, s=30, c='white', edgecolor='k', depthshade=False)
    ax.set_xlim( [-1,1] )
    ax.set_ylim( [-1,1] )
    ax.set_zlim( [-1,1] )



    ## HISTOGRAM (Fig 2, right)
    ax = fig2.add_subplot(122)
    #ax.set_title('Sol Histogram Glitches', pad=10)
    ax.grid(ls='-.', lw=0.5, zorder=0)
    ax.set(xlabel='Sol', ylabel=' Number of glitches / Sol')
    ax.autoscale(enable=True, axis='both', tight=False)
    ax.hist(glitch_start_hist_SOL_mat, bins=np.arange(min(sols_range), max(sols_range)+2,1), color=colors_hist, align='left', lw=0.5, stacked=True, zorder=3, rasterized=True) # "+2" instead of "+1" because otherwise there is a bug in the last bin
    #ax.set_xlim([min(sols_range), max(sols_range)])
    ax.set_ylim([0,ymax_hist])
    ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True, nbins='auto'))
    ax2   = ax.twiny()
    ax1Xs = ax.get_xticks()
    ax2Xs = [marstime(LMST=UTCDateTime('1969-12-31T00:00:00.000000Z')+datetime.timedelta(days=e)).UTC_time.strftime('%m-%d') for e in ax1Xs]
    #print(ax1Xs)
    #print(ax2Xs)
    
    ax2.set_xticks(ax1Xs)
    ax2.set_xbound(ax.get_xbound())
    ax2.set_xticklabels(ax2Xs)
    ax2.set_xlabel(u"Day")



    ### FINAL
    if outfile:
        plt.savefig(outfile2)
        print(outfile2)
def glitch_SOLoverLMST_plot(*glitch_files, run=True, mode='AZ', outfile='', **kwargs):

    """

    """


    ### RUN OR NOT:
    if not run:
        return

    now = time.time()



    ### HELPER FUNCTION
    picked_already = False
    def onpick(event):
        nonlocal picked_already
        if not picked_already:
            picked_already = True
            header         = u'#Number                  UTC          LMST      U-AMP      V-AMP      W-AMP      Z-AMP      N-AMP      E-AMP       AZ      INC'
            print(header)

        indices = event.ind
        for index in indices:
            print(u'%7s  %s  %s  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %7.1f  %7.1f' % (all_glitches[index,0],
                                                                                            glitch_times[index].UTC_time.strftime('%Y-%m-%dT%H:%M:%S'), 
                                                                                            glitch_times[index].LMST_string, 
                                                                                            float(all_glitches[index,11]),
                                                                                            float(all_glitches[index,12]),
                                                                                            float(all_glitches[index,13]),
                                                                                            float(all_glitches[index,14]),
                                                                                            float(all_glitches[index,15]),
                                                                                            float(all_glitches[index,16]),
                                                                                            float(all_glitches[index,17]),
                                                                                            float(all_glitches[index,19]),
                                                                                            ))
        print(u'- - -')



    ### OUTPUT
    print()
    print(u'  -------------------------------')
    print(u'  RUNNING GLITCH LMST sol PLOTTER')
    print(u'  -------------------------------')
    print()



    ### SELECT GLITCHES
    all_glitches = select_glitches(glitch_files, **kwargs)
    LMST_range   = kwargs['LMST_range']



    ### VARIABLES
    glitch_times = [marstime(glitch[1]) for glitch in all_glitches]
    LMSTs        = np.array( [time.LMST_time.datetime-datetime.timedelta(days=time.sol) for time in glitch_times] )
    sols         = np.array( [time.sol for time in glitch_times] )
    amps         = np.array( np.abs( [glitch[11:14].astype('float') for glitch in all_glitches] ))
    AZs          = np.array( [glitch[17].astype('float') for glitch in all_glitches] )
    mini_amp     = np.min( np.abs( amps ))
    maxi_amp     = np.max( np.abs( amps ))
    sols_range   = sorted( set( sols ))
    titles       = ['U-component','V-component','W-component']



    ### OUTPUT
    print()
    print(u'Detected sols: %s .. %s '         % (min(sols), max(sols)))

    

    ### PLOT
    if mode.upper() == 'AMP':
        fig, axes = plt.subplots(1, 3, figsize=(18, 7), sharex=True, sharey=True)
        titles    = ['U-component','V-component','W-component']
        cmap      = plt.get_cmap('viridis')
        norm      = mpl.colors.LogNorm(vmin=mini_amp, vmax=maxi_amp)
        label     = 'Glitch amplitude (gain corrected)'
        scalarMap = cmx.ScalarMappable(norm=norm, cmap=cmap)
        colours   = np.array( [[scalarMap.to_rgba(abs(amp[0])),scalarMap.to_rgba(abs(amp[1])),scalarMap.to_rgba(abs(amp[2]))] for amp in amps] )

    else:
        fig, axes = plt.subplots(1, 1, figsize=(10, 8), sharex=True, sharey=True)
        axes      = [axes]
        titles    = ['']
        cmap      = plt.get_cmap('hsv')
        norm      = mpl.colors.Normalize(vmin=0, vmax=360)
        label     = r'Glitch azimuth ($^{\circ}$ from N)'
        scalarMap = cmx.ScalarMappable(norm=norm, cmap=cmap)
        colours   = np.array( [scalarMap.to_rgba(AZ) for AZ in AZs] )

    cax1 = fig.add_axes([0.05, 0.5, 0.005, 0.4])             # left, bottom, width, height
    cb1  = mpl.colorbar.ColorbarBase(cax1, drawedges=False, cmap=cmap, norm=norm, orientation='vertical')
    cb1.set_label(label)
    cax1.yaxis.set_ticks_position('left')
    if not mode.upper() == 'AMP':
        cax1.yaxis.set_ticks(np.arange(0, 361, 45))

    fig.canvas.set_window_title('Glitch LMST sol plot')
    fig.suptitle('LMST sol plotter: Sols=%s..%s, %s glitches' % (sols_range[0], sols_range[-1], len(all_glitches)), fontsize=11, y=0.99)
    fig.subplots_adjust(wspace=0.4, hspace=0.4)
    fig.canvas.mpl_connect('pick_event', onpick)
    fig.autofmt_xdate()

    for l, ax in enumerate(axes):

        ax.set_title(titles[l], size=10)
        ax.autoscale(enable=True, axis='both', tight=False)
        ax.grid(ls='-.', lw=0.5, zorder=1)
        ax.set_xlabel('LMST')
        ax.xaxis.labelpad = 3
        ax.xaxis.set_major_formatter( mdates.DateFormatter('%H:%M:%S') )
        ax.xaxis.set_major_locator(plt.MaxNLocator(5))
        ax.set_ylabel('Sol')
        ax.set_ylim([min(sols_range),max(sols_range)+1])

        if mode.upper() == 'AMP':
            c = colours[:,l]
        else:
            c = colours
        ax.scatter(LMSTs, sols, s=5, linewidths=.3, c=c, edgecolor='k', rasterized=True, picker=2)


    if outfile:
        plt.savefig(outfile)
        print(outfile)
def plot_glitch_gutenberg(*glitch_files, run=True, outfile='', **kwargs):

    """

    """



    ### RUN OR NOT:
    if not run:
        return



    ### OUTPUT
    print()
    print(u'  --------------------------------')
    print(u'  RUNNING GLITCH GUTENBERG PLOTTER')
    print(u'  --------------------------------')
    print()



    ### SELECT GLITCHES
    all_glitches = select_glitches(glitch_files, **kwargs)



    ### PLOT HISTOGRAMS
    label_comps = 'UVWZNE'
    range_hist  = (1e-10,1e-5)

    fig, axes = plt.subplots(2,3, figsize=(12,7), sharex=True, sharey=True)
    fig.canvas.set_window_title('Glitch Gutenberg-Richter plot')
    fig.suptitle('Gutenberg plotter: %s glitches' % len(all_glitches), fontsize=11, y=0.99)
    plt.xscale('log')

    for l, ax in enumerate(axes.flatten()):

        # set-up plot
        ax.grid(ls='-.', lw=0.5, zorder=0)
        ax.set(aspect='equal', xlabel='Glitch amplitudes (m/s)', ylabel='Number of occurences')
        ax.set_xlim(range_hist)
        ax.set_ylim(1,1e5)
        ax.xaxis.set_tick_params(labelbottom=True)
        ax.yaxis.set_tick_params(labelbottom=True)

        # glitches for that happen on U, V, and/or W; and on ZNE 
        if l<=2:
            glitches = all_glitches[all_glitches[:,5+l]=='1']
        else:
            glitches = all_glitches

        # create bins but only needed once
        if l==0:
            hist, bins = np.histogram( np.abs(glitches[:,l+11].astype('float')), bins=30)
            logbins    = np.logspace(np.log10(bins[0]), np.log10(bins[-1]),len(bins))

        # plot histogram
        y, x, _         = ax.hist( np.abs(glitches[:,l+11].astype('float')), log=True, range=range_hist, color='k', bins=logbins, ec='w', lw=0.5, label=label_comps[l]+'-comp')
        
        # calculating fit
        x               = x[1:][y!=0]
        y               = y[y!=0]
        y_max_index     = np.argmax(y)
        b, a            = np.polyfit( np.log10(x[y_max_index:]), np.log10(y[y_max_index:]), 1)
        
        # plot fit
        x_int           = x[y_max_index:]
        y_b             = 10**(a+np.log10(x_int)*b)
        ax.plot(x_int, y_b, lw =2, c='darkorange', label='b=%.2f' % b)
        
        # swap order of legends
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[::-1], labels[::-1], loc='upper right')
        #print(u'%s: b=%.2f, a=%.2f' % (label_comps[l], b, a))



    ### PRINT INFORMATION
    #print(u'Please treat `a` values with caution.')
    #print(u'Fits have been done with negative exponents, so the `a` are not quite right!')



    ### FINAL
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    if outfile:
        plt.savefig(outfile)
        print(outfile)
def plot_glitch_ppol(glitch_start=None, glitch_length=30, run=True, waveform_files=[], inventory_file='IRIS', show=True, outfile='', **kwargs):

    """

    """



    ### RUN OR NOT:
    if not run:
        return



    ### OUTPUT
    print()
    print(u'  ---------------------------')
    print(u'  RUNNING GLITCH PPOL PLOTTER')
    print(u'  ---------------------------')
    print()



    ### READ GLITCH-FILE
    if glitch_start is None:
        print(u'No starttime given of glitch, cannot know which one you would like to plot.')
        return
    else:
        glitch_start = UTCDateTime(glitch_start)

    glitch_end = glitch_start + glitch_length
    title      = u'Glitch at %s' % marstime(glitch_start).UTC_string


    ### LOOP OVER ALL WAVEFORM FILES AND FIND CORRECT ONE
    for waveform_file in waveform_files:
        
        # reading stream
        try: 
            stream = read(waveform_file, starttime=glitch_start-30, endtime=glitch_end+30, headonly=True)
        except TypeError:       # could not read file, skip
            continue             

        if not stream:
            continue
        else:
            stream = read2(waveform_file)
            stream_select = stream.select(channel='?[LMH]?')
            stream_select.trim(starttime=glitch_start-10, endtime=glitch_end+10)

        # data pre processing
        stream_select.set_inventory(inventory_file)
        stream_select.gain_correction(verbose=False)
        stream2 = stream_select.copy()
        stream2.rotate('->ZNE', inventory=stream_select.inventory, components=('UVW'))
        #stream2.detrend('linear')
        #stream2.detrend('demean')
        #stream2.taper(0.05)
        #stream2.filter('highpass', freq=0.1, corners=2, zerophase=True)
        #stream2.filter('lowpass', freq=2, corners=2, zerophase=True)

        # ppol measurement + plot
        ppol_measurement = ppol(stream=stream2, demean=True, fix_angles='EQ', starttime=glitch_start, endtime=glitch_end, Xoffset_samples_for_amplitude=1*stream2[0].stats.sampling_rate)
        ppol_measurement.display(tag=title)
        print()
        print(u'Glitch times:')
        print(u'    UTC: %20s - %s' % (marstime(glitch_start).UTC_string,  marstime(glitch_end).UTC_string))
        print(u'   LMST: %20s - %s' % (marstime(glitch_start).LMST_string, marstime(glitch_end).LMST_string))
        print(u'Info: Results may slightly deviate from those of the glitch detector.')
        ppol_measurement.plot(title=title, show=show, outfile=outfile, vertsc=['k'], original_data=stream_select.detrend('demean'))

        break

    else:
        print(u'Could not find any waveform data corresponding to glitch time.')



### _ _ N A M E _ _ = = " _ _ M A I N _ _ "  
if __name__ == "__main__":
    pass