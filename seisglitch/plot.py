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


# CHANGE DETECTOR PLOT SO IT SPLITS LONG RECORDS

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
from seisglitch.util import read2, marstime, sec2hms, ppol, quick_plot
from seisglitch.math import normalise



### PLOT SCRIPTS TO VISUALIZE

#    header         = u'   0                    1                    2              3              4           5         6         7             8          9         10           11         12         13         14         15         16            17              18              19              20              21          22\n' \
#                     u' NUM         GLITCH-UTC-S         GLITCH-UTC-E  GLITCH-LMST-S  GLITCH-LMST-E    U-GLITCH  V-GLITCH  W-GLITCH         U-RAW      V-RAW      W-RAW        U-GAI      V-GAI      W-GAI      Z-GAI      N-GAI      E-GAI    PHI_3D_GAI  PHI_ERR_3D_GAI      INC_3D_GAI  INC_ERR_3D_GAI      SNR_3D_GAI  POL_3D_GAI\n' \
#                     u
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
    **kwargs):


    UNIT      = UNIT.upper()
    variables = locals().copy()
    del variables['glitch_files']
    del variables['kwargs']


    ### READ GLITCH-FILE
    all_glitches = []
    for gf, glitch_file in enumerate(glitch_files):
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
                keep = [g for g in range(len(keep)) if float(all_glitches[g,column])<=float(Amplitudes[comp][0])]
            elif len(Amplitudes[comp])==2:
                keep = [g for g in range(len(keep)) if float(all_glitches[g,column])>=float(Amplitudes[comp][0]) and float(all_glitches[g,column])<=float(Amplitudes[comp][1])]
            


    ### SELECT GLITCHES ON GIVEN `AZs`
    if AZs:
        AZs = [float(BAZ) for BAZ in AZs]

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
            keep = [g for g in range(len(keep)) if float(all_glitches[g,21])<=float(glitch_SNR[0])]
        elif len(glitch_SNR)==2:
            keep = [g for g in range(len(keep)) if float(all_glitches[g,21])>=float(glitch_SNR[0]) and float(all_glitches[g,21])<=float(glitch_SNR[1])]



    ### SELECT GLITCHES ON GIVEN `glitch_polarization`
    if glitch_polarization:
        if len(glitch_polarization)==1:
            keep = [g for g in range(len(keep)) if float(all_glitches[g,22])<=float(glitch_polarization[0])]
        elif len(glitch_polarization)==2:
            keep = [g for g in range(len(keep)) if float(all_glitches[g,22])>=float(glitch_polarization[0]) and float(all_glitches[g,22])<=float(glitch_polarization[1])]



    ### SELECT INVERSE SELCTION
    if not inverse_selection:
        all_glitches = all_glitches[keep]
    else:
        inverse      = [i for i in range(len(all_glitches)) if i not in keep]
        all_glitches = all_glitches[inverse]



    ### PRINTING VARIABLES:
    print(u'SELECTED %s GLITCHES ON THE FOLLOWING PARAMETERS:' % len(all_glitches))
    for key, value in sorted(variables.items(), key=lambda item: len(item[0])):
        print(u'  %25s = %s' % (key, value))
    print()

    if all_glitches.size == 0:
        print(u'WARNING: With given parameters no glitches were selected.')
        sys.exit()
    else:
        return all_glitches

def plot_glitch_detector(*glitch_files, run=True, waveform_files=[], starttime=None, endtime=None, components='*', show=True, outfile='', **kwargs):
    
    """
    Plot glitches, based on glitch file produced by function `glitch_detector()`.t.
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
        print(u'ERROR: You need to specify `waveform_files` in the config.yml so waveform data can be read.')
        sys.exit()



    ### SELECT GLITCHES
    all_glitches = select_glitches(glitch_files, **kwargs)



    ### LOOP OVER ALL FILES PASSED
    for k, waveform_file in enumerate(waveform_files):

        # read / select stream
        st = read2(waveform_file)
        if components != '*':
            comps = '['+components+']'
        else:
            comps = components
        st_select = st.select(component=comps)

        # sanity check
        if not st_select:
            print(u'WARNING: No `components` %s found in %s. No plotting done.' % (comps, waveform_file))
            continue


        # if times specified
        if starttime or endtime:
            st_select.trim2(starttime=starttime, endtime=endtime)
        
        # another sanity check
        if not st_select:
            print(u'WARNING: Times %s, and %s not present in %s. No plotting done.' % (starttime, endtime, waveform_file))
            continue
        else:
            print(u'Plotting: %s' % waveform_file)

        # get times
        glitch_times = []
        stream_times = st_select.times

        if 'U' in comps.upper() or 'V' in comps.upper() or 'W' in comps.upper():

            if 'U' in comps.upper():
                glitch_times += [glitch[1] for glitch in all_glitches if glitch[5]=='1']

            if 'V' in comps.upper():
                glitch_times += [glitch[1] for glitch in all_glitches if glitch[6]=='1']

            if 'W' in comps.upper():
                glitch_times += [glitch[1] for glitch in all_glitches if glitch[7]=='1']

        else:
            glitch_times = all_glitches[:,1]

        glitch_times = [marstime(time).UTC_time for time in glitch_times if time>=str(stream_times[0]) and time<=str(stream_times[1])]

        # make plot
        title     = '%s glitches (components=%s)' % (len(glitch_times),comps)
        win_title = 'Glitch dectector plot'
        if outfile:
            fileout = '.'.join( outfile.split('.')[:-1]) + '%d.png' % (k+1)
        else:
            fileout = None
        quick_plot(*st_select, title=title, win_title=win_title, verts=[glitch_times], xlabel='Time', show=show, outfile=fileout)
def plot_glitch_overview(*glitch_files, run=True, waveform_files=[], glitch_length=None, outfile='', **kwargs):
    
    """
    Plot glitches, based on glitch file produced by function `glitch_detector()`.

    amp parameters apply to all components: UVWZNE
    min_amp: at least one component amplitude must be larger
    max_amp: all component amplitudes must be smaller

    UNIT allows to choose which polarizations you wish to plot.

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
            header         = u'#Number                  UTC          LMST    U    V    W      U-AMP      V-AMP      W-AMP      Z-AMP      N-AMP      E-AMP      BAZ-{0}    INC-{0}  SNR_3D-{0}  POL_3D-{0}'.format(UNIT)
            print(header)

        for index in indices:
            print(u'%7s  %s  %s  %3s  %3s  %3s  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g    %9.1f  %9.1f   %9.3g   %9.3f' % (all_glitches[index,0],
                                                                                                              marstime(all_glitches[index,1]).UTC_time.strftime('%Y-%m-%dT%H:%M:%S'), 
                                                                                                              marstime(all_glitches[index,1]).LMST, 
                                                                                                              int(all_glitches[index,5]),
                                                                                                              int(all_glitches[index,6]),
                                                                                                              int(all_glitches[index,7]),
                                                                                                              float(all_glitches[index,columns[UNIT][0]]),
                                                                                                              float(all_glitches[index,columns[UNIT][1]]),
                                                                                                              float(all_glitches[index,columns[UNIT][2]]),
                                                                                                              float(all_glitches[index,columns[UNIT][3]]),
                                                                                                              float(all_glitches[index,columns[UNIT][4]]),
                                                                                                              float(all_glitches[index,columns[UNIT][5]]),
                                                                                                              float(all_glitches[index,columns[UNIT][6]]),
                                                                                                              float(all_glitches[index,columns[UNIT][7]]),
                                                                                                              float(all_glitches[index,columns[UNIT][8]]),
                                                                                                              float(all_glitches[index,columns[UNIT][9]]),
                                                                                                              ))
        print(u'- - -')

        # print first picked glitches
        glitch_start = all_glitches[indices[0],1]
        plot_glitch_ppol(glitch_start, waveform_files=waveform_files, glitch_length=glitch_length, **kwargs)



    ### OUTPUT
    print()
    print(u'  -------------------------------')
    print(u'  RUNNING GLITCH OVERVIEW PLOTTER')
    print(u'  -------------------------------')
    print()



    ### SELECT GLITCHES
    all_glitches = select_glitches(glitch_files, **kwargs)
    UNIT         = kwargs['UNIT']
    LMST_range   = kwargs['LMST_range']



    ### ASSIGN NEEDED VARIABLES
    sols_range   = sorted( [marstime(e).sol for e in all_glitches[:,1]] )
    print(u'Detected sols: %s .. %s ' % (sols_range[0], sols_range[-1]))

    ymax_hist    = Counter( sols_range ).most_common()[0][1]

    try:
        glitch_starts     = np.array( [marstime(e) for e in all_glitches[:,1]] )

        Z_amps            = all_glitches[:,columns[UNIT][3]].astype('float')
        N_amps            = all_glitches[:,columns[UNIT][4]].astype('float')
        E_amps            = all_glitches[:,columns[UNIT][5]].astype('float')
        phis              = all_glitches[:,columns[UNIT][6]].astype('float') * np.pi/180
        incs              = all_glitches[:,columns[UNIT][7]].astype('float') * np.pi/180
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
    fig1 = plt.figure(figsize=(17,8))
    fig1.canvas.set_window_title('Glitch overview plot 1')
    fig1.suptitle('%s overview plotter: Sols=%s..%s, LMST %s-%s, %s glitches' % (UNIT, sols_range[0], sols_range[-1], LMST_range[0], LMST_range[1], len(all_glitches)), fontsize=11, y=0.99)
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
    ax.set_title('Back-azimuths', pad=20, size=10)
    ax.spines['polar'].set_visible(False)
    ax.tick_params(pad=5)

    ax.set_theta_direction(-1)
    ax.set_theta_zero_location('N')
    ax.set_thetagrids([0,15,37,70,90,105,114,135,157,180,251,255,270,277,345.5], labels=['','LSA/VBBV\n(15°)','LVL1\n(37°)', r'HP$^{3}$'+'\n(70°)','','SP2\n(105°)','WTS E\n(114°)',r'VBBU'+'\n(135°)','LVL2\n(157°)','','WTS W (251°)','VBBW\n(255°)','','LVL3\n(277°)','SP3/WTS N\n(345°)'], size=8)

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
    ax.set_title('Incidence angles', pad=20, size=10)
    ax.spines['polar'].set_visible(False)
    ax.tick_params(pad=5)

    ax.set_theta_direction(-1)
    ax.set_theta_zero_location('N')
    ax.set_thetamin(0)
    ax.set_thetamax(180)
    ax.set_thetagrids([0,48.5,60.5,90,131.5,119.5,180], labels=['0°','48.5°','VBB\n(60.5°)','90°','131.5°','VBB\n(119.5°)'], size=8)

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
    fig2 = plt.figure(figsize=(17,8))
    fig2.canvas.set_window_title('Glitch overview plot 2')
    fig2.suptitle('%s overview plotter: Sols=%s..%s, LMST %s-%s, %s glitches' % (UNIT, sols_range[0], sols_range[-1], LMST_range[0], LMST_range[1], len(all_glitches)), fontsize=11, y=0.99)
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
    ax.set_title('3-D Polarizations', y=1.10, size=10)
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
    ax2Xs = [marstime(LMST_time=UTCDateTime('1970-01-01T00:00:00.000000Z')+datetime.timedelta(days=e)).UTC_time.strftime('%m-%d') for e in ax1Xs]
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
def plot_glitch_LMST_sol(*glitch_files, run=True, mode='BAZ', outfile='', **kwargs):



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
            header         = u'#Number                  UTC          LMST      U-{0}      V-{0}      W-{0}      Z-{0}      N-{0}      E-{0}  BAZ-{0}  INC-{0}'.format(UNIT)
            print(header)

        indices = event.ind
        for index in indices:
            print(u'%7s  %s  %s  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %7.1f  %7.1f' % (all_glitches[index,0],
                                                                                            glitch_times[index].UTC_time.strftime('%Y-%m-%dT%H:%M:%S'), 
                                                                                            glitch_times[index].LMST, 
                                                                                            float(all_glitches[index,columns[UNIT][0]]),
                                                                                            float(all_glitches[index,columns[UNIT][1]]),
                                                                                            float(all_glitches[index,columns[UNIT][2]]),
                                                                                            float(all_glitches[index,columns[UNIT][3]]),
                                                                                            float(all_glitches[index,columns[UNIT][4]]),
                                                                                            float(all_glitches[index,columns[UNIT][5]]),
                                                                                            float(all_glitches[index,columns[UNIT][6]]),
                                                                                            float(all_glitches[index,columns[UNIT][7]]),
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
    UNIT         = kwargs['UNIT']
    LMST_range   = kwargs['LMST_range']
    min_amp      = kwargs['min_amp']
    max_amp      = kwargs['max_amp']



    ### VARIABLES
    glitch_times = [marstime(glitch[1]) for glitch in all_glitches]
    LMSTs        = np.array( [time.LMST_time.datetime-datetime.timedelta(days=time.sol) for time in glitch_times] )
    sols         = np.array( [time.sol for time in glitch_times] )
    amps         = np.array( np.abs( [glitch[columns[UNIT][0:3]].astype('float') for glitch in all_glitches] ))
    AZs         = np.array( [glitch[columns[UNIT][6]].astype('float') for glitch in all_glitches] )
    mini_amp     = np.min( np.abs( amps ))
    maxi_amp     = np.max( np.abs( amps ))
    sols_range   = sorted( set( sols ))
    titles       = ['U-component','V-component','W-component']



    ### OUTPUT
    print()
    print(u'Detected sols: %s .. %s '         % (min(sols), max(sols)))
    print()

    

    ### PLOT
    if mode.upper() == 'AMP':
        fig, axes = plt.subplots(1, 3, figsize=(18, 7), sharex=True, sharey=True)
        titles    = ['U-component','V-component','W-component']
        cmap      = plt.get_cmap('viridis')
        norm      = mpl.colors.LogNorm(vmin=mini_amp, vmax=maxi_amp)
        label     = '%s amplitude' % UNIT
        scalarMap = cmx.ScalarMappable(norm=norm, cmap=cmap)
        colours   = np.array( [[scalarMap.to_rgba(abs(amp[0])),scalarMap.to_rgba(abs(amp[1])),scalarMap.to_rgba(abs(amp[2]))] for amp in amps] )

    else:
        fig, axes = plt.subplots(1, 1, figsize=(10, 8), sharex=True, sharey=True)
        axes      = [axes]
        titles    = ['']
        cmap      = plt.get_cmap('hsv')
        norm      = mpl.colors.Normalize(vmin=0, vmax=360)
        label     = 'Glitch azimuth'
        scalarMap = cmx.ScalarMappable(norm=norm, cmap=cmap)
        colours   = np.array( [scalarMap.to_rgba(BAZ) for BAZ in AZs] )

    cax1 = fig.add_axes([0.005, 0.5, 0.005, 0.4]) # left, bottom, width, height
    cax1.yaxis.set_ticks_position('left')
    cb1 = mpl.colorbar.ColorbarBase(cax1, drawedges=False, cmap=cmap, norm=norm, orientation='vertical')
    cb1.set_label(label)

    fig.canvas.set_window_title('Glitch LMST sol plot')
    fig.suptitle('%s LMST sol plotter: Sols=%s..%s, LMST %s-%s, %s glitches (UVW: %g ≤ 1comp ≤ %g)' % (UNIT, sols_range[0], sols_range[-1], LMST_range[0], LMST_range[1], len(all_glitches), min_amp, max_amp), fontsize=11, y=0.99)
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
        #ax.xaxis.set_major_locator( mdates.HourLocator(interval=4) )
        ax.set_ylabel('Sol')
        ax.set_ylim([min(sols_range),max(sols_range)+1])

        if mode.upper() == 'AMP':
            c = colours[:,l]
        else:
            c = colours
        ax.scatter(LMSTs, sols, s=5, linewidths=.3, c=c, edgecolor='k', rasterized=True, picker=2)


    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
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
    #all_glitches = all_glitches[all_glitches[:,11].astype('float')<1e-6]
    UNIT         = kwargs['UNIT']
    LMST_range   = kwargs['LMST_range']



    ### PLOT HISTOGRAMS
    label_comps = 'UVWZNE'
    range_hist  = (1e-10,1e-5)

    fig, axes = plt.subplots(2,3, figsize=(12,7), sharex=True, sharey=True)
    fig.canvas.set_window_title('Glitch Gutenberg-Richter plot')
    fig.suptitle('%s Gutenberg plotter: LMST %s-%s, %s glitches' % (UNIT, LMST_range[0], LMST_range[1], len(all_glitches)), fontsize=11, y=0.99)
    #fig.subplots_adjust(wspace=0.5, hspace=0.0)    
    plt.xscale('log')

    for l, ax in enumerate(axes.flatten()):

        # set-up plot
        ax.grid(ls='-.', lw=0.5, zorder=0)
        ax.set(aspect='equal', xlabel='Glitch amplitudes (m/s)', ylabel='Number of occurences')
        ax.set_xlim(range_hist)
        ax.set_ylim(1,1e5)
        ax.xaxis.set_tick_params(labelbottom=True)
        ax.yaxis.set_tick_params(labelbottom=True)

        # create bins, but only needed for first one
        if l<=2:
            glitches = all_glitches[all_glitches[:,5+l]=='1']
        else:
            glitches = all_glitches

        if l==0:
            hist, bins = np.histogram( np.abs(glitches[:,columns[UNIT][l]].astype('float')), bins=30)
            logbins    = np.logspace(np.log10(bins[0]), np.log10(bins[-1]),len(bins))

        # plot histogram
        y, x, _         = ax.hist( np.abs(glitches[:,columns[UNIT][l]].astype('float')), log=True, range=range_hist, color='k', bins=logbins, ec='w', lw=0.5, label=label_comps[l]+'-comp')
        
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
    print(u'Please treat `a` values with caution.')
    print(u'Fits have been done with negative exponents, so they the `a` are not quite right!')



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
    title      = 'Glitch at %s' % marstime(glitch_start).UTC


    ### LOOP OVER ALL WAVEFORM FILES AND FIND CORRECT ONE
    for waveform_file in waveform_files:
        
        # reading stream
        stream = read(waveform_file, starttime=glitch_start-30, endtime=glitch_end+30, headonly=True)
        if not stream:
            continue
        else:
            stream = read2(waveform_file)
            stream.trim(starttime=glitch_start-10, endtime=glitch_end+10)

        # data pre processing
        stream._set_inventory(inventory_file)
        #stream.write('/home/scholz/Desktop/data/VBB_test_glitch.mseed', starttime=glitch_start, endtime=glitch_end)
        stream.gain_correction(verbose=False)
        #stream.filter('lowpass', freq=1, zerophase=True)
        stream2 = stream.copy()
        stream2.rotate('->ZNE', inventory=stream.inventory, components=('UVW'))

        # ppol measurement + plot
        ppol_measurement = ppol(stream=stream2, demean=False, fix_angles='AMP', starttime=glitch_start, endtime=glitch_end, Xoffset_samples_for_amplitude=1*stream2[0].stats.sampling_rate)
        ppol_measurement.display(tag=title)
        print()
        print(u'Glitch times:')
        print(u'    UTC: %20s - %s' % (marstime(glitch_start).UTC,  marstime(glitch_end).UTC))
        print(u'   LMST: %20s - %s' % (marstime(glitch_start).LMST, marstime(glitch_end).LMST))
        print(u'Info: Results may slightly differ from those of the glitch detector.')
        ppol_measurement.plot(title=title, show=show, outfile=outfile, original_data=stream.detrend('demean'))

        break

    else:
        print(u'Could not find any waveform data corresponding to glitch time.')

def plot_glitch_waveform(*glitch_files, run=True, waveform_files=[], sols=[], LMST_range=['0', '24'], min_amp=None, max_amp=None, AZs=[], scale=1, outfile='', **kwargs):



    ### RUN OR NOT:
    if not run:
        return

    now = time.time()



    ### HELPER FUNCTION
    def onclick(event):
        if event.dblclick:
            sol       = int(event.ydata)
            LMST_time = UTCDateTime(mdates.num2date(event.xdata) + datetime.timedelta(days=sol))
            time_obj  = marstime(LMST_time=LMST_time)
            print(time_obj)



    ### OUTPUT
    print()
    print(u'  -------------------------------')
    print(u'  RUNNING GLITCH WAVEFORM PLOTTER')
    print(u'  -------------------------------')
    print()
    print(u"Have to convert each data point's UTC time to LMST, this takes a minute ..")
    print()



    ### SANITY CHECK:
    if not waveform_files:
        print()
        print(u'ERROR: Cannot plot without `waveform_files`.')
        sys.exit()



    ### READ GLITCH-FILE
    all_glitches = []
    for glitch_file in glitch_files:
        glitches      = np.loadtxt(glitch_file, dtype='str')
        #glitches      = glitch_exclude(glitches, verbose=False)
        all_glitches += list(glitches)
    all_glitches = np.array(all_glitches)



    ### SELECT GLITCHES ON GIVEN `SOLS`
    if sols:
        sols = sorted(sols)
        all_glitches = np.array( [glitch for glitch in all_glitches if marstime(glitch[1]).sol>=sols[0] and marstime(glitch[1]).sol<=sols[1]] )



    ### SELECT GLITCHES IN SPECIFIED LMST RANGE
    LMST_range        = [(str(t).replace(':','') + '000000')[:6] for t in LMST_range]
    LMST_range        = [[t[i:i+2] for i in range(0,len(t),2)] for t in LMST_range]
    LMST_range        = [':'.join(t) for t in LMST_range]
    
    glitch_start_LMST = np.array( [glitch[3].split('M')[1] for glitch in all_glitches] )
    if LMST_range[0]<=LMST_range[1]:
        indices       = np.where( (glitch_start_LMST >= LMST_range[0]) & (glitch_start_LMST <= LMST_range[1]) )
    else:
        indices       = np.where( (glitch_start_LMST >= LMST_range[0]) | (glitch_start_LMST <= LMST_range[1]) )
    all_glitches      = all_glitches[ indices ]



    ### SELECT GLITCHES TO `MIN_AMP` AND `MAX_AMP`
    if eval(str(min_amp).capitalize()) is None or eval(str(min_amp).capitalize()) == False:
        min_amp = 0
    else:
        min_amp = float( min_amp )
    if eval(str(max_amp).capitalize()) is None or eval(str(max_amp).capitalize()) == False:
        max_amp = 1e9
    else:
        max_amp = float( max_amp )   

    all_glitches = all_glitches[((abs(all_glitches[:,columns['GAI'][0]].astype('float'))>=min_amp)  | \
                                 (abs(all_glitches[:,columns['GAI'][1]].astype('float'))>=min_amp)  | \
                                 (abs(all_glitches[:,columns['GAI'][2]].astype('float'))>=min_amp)  | \
                                 (abs(all_glitches[:,columns['GAI'][3]].astype('float'))>=min_amp)  | \
                                 (abs(all_glitches[:,columns['GAI'][4]].astype('float'))>=min_amp)  | \
                                 (abs(all_glitches[:,columns['GAI'][5]].astype('float'))>=min_amp)) & \

                                 (abs(all_glitches[:,columns['GAI'][0]].astype('float'))<=max_amp)  & \
                                 (abs(all_glitches[:,columns['GAI'][1]].astype('float'))<=max_amp)  & \
                                 (abs(all_glitches[:,columns['GAI'][2]].astype('float'))<=max_amp)  & \
                                 (abs(all_glitches[:,columns['GAI'][3]].astype('float'))<=max_amp)  & \
                                 (abs(all_glitches[:,columns['GAI'][4]].astype('float'))<=max_amp)  & \
                                 (abs(all_glitches[:,columns['GAI'][5]].astype('float'))<=max_amp)]



    ### SELECT GLITCHES ON GIVEN `AZs`
    if AZs:
        AZs = [float(BAZ) for BAZ in AZs]
        if AZs[0]<=AZs[1]:    # e.g.: 40 - 300 deg
            all_glitches = np.array( [glitch for glitch in all_glitches if float(glitch[columns['GAI'][6]])>=AZs[0] and float(glitch[columns['GAI'][6]])<=AZs[1]] )
        else:                   # e.g.: 300 - 40 deg
            all_glitches = np.array( [glitch for glitch in all_glitches if float(glitch[columns['GAI'][6]])>=AZs[0] or  float(glitch[columns['GAI'][6]])<=AZs[1]] )



    ### VARIABLES
    streams         = [read(file) for file in waveform_files]
    glitch_max_amp  = np.max( np.abs( all_glitches[:,8:11].astype('float') ))*1
    comps           = []
    times           = []
    datas           = []
    sols_range      = []
    counter_counter = 0



    ### LOOP OVER GLITCHES, PREPARE DATE FOR PLOTTING
    for i, glitch in enumerate(all_glitches):

        if (i+1)%10 == 0:
            print(u'Handling glitch %06d / %06d' % (i+1, len(all_glitches)))

        glitches_start_UTC         = UTCDateTime(glitch[1])
        glitches_end_UTC           = UTCDateTime(glitch[2])
        glitches_start_LMST_string = glitch[3]
        glitches_end_LMST_string   = glitch[4]
        glitch_start_sol           = int(glitch[3].split('M')[0])
        glitch_end_sol             = int(glitch[4].split('M')[0])
        glitch_U                   = glitch[5]
        glitch_V                   = glitch[6]
        glitch_W                   = glitch[7]


        for stream in streams:

            stream_glitch = stream.slice(starttime=glitches_start_UTC, endtime=glitches_end_UTC)


            if stream_glitch:

                glitch_times_UTC = stream_glitch[0].times('utcdatetime')
                glitch_times_obj = np.array( [marstime(glitch_time) for glitch_time in glitch_times_UTC] )

                if glitch_start_sol == glitch_end_sol:

                    glitch_comps = ''                  
                    if glitch_U == '1':
                        glitch_comps += 'U'
                    if glitch_V == '1':
                        glitch_comps += 'V'
                    if glitch_W == '1':
                        glitch_comps += 'W'

                    glitch_times  = [glitch_time_obj.LMST_time.datetime-datetime.timedelta(days=glitch_start_sol) for glitch_time_obj in glitch_times_obj]
                    
                    glitch_data_U = stream_glitch.select(component='U')[0]
                    glitch_data_U.detrend('demean')
                    glitch_data_U = glitch_data_U.data
                    glitch_data_U = normalise(glitch_data_U, scale_to_between=[np.min(glitch_data_U)/glitch_max_amp*0.5*scale, np.max(glitch_data_U)/glitch_max_amp*0.5*scale])
                    glitch_data_U = glitch_data_U + glitch_start_sol + 0.5
                    glitch_data_V = stream_glitch.select(component='V')[0]
                    glitch_data_V.detrend('demean')
                    glitch_data_V = glitch_data_V.data
                    glitch_data_V = normalise(glitch_data_V, scale_to_between=[np.min(glitch_data_V)/glitch_max_amp*0.5*scale, np.max(glitch_data_V)/glitch_max_amp*0.5*scale])
                    glitch_data_V = glitch_data_V + glitch_start_sol + 0.5
                    glitch_data_W = stream_glitch.select(component='W')[0]
                    glitch_data_W.detrend('demean')
                    glitch_data_W = glitch_data_W.data
                    glitch_data_W = normalise(glitch_data_W, scale_to_between=[np.min(glitch_data_W)/glitch_max_amp*0.5*scale, np.max(glitch_data_W)/glitch_max_amp*0.5*scale])
                    glitch_data_W = glitch_data_W + glitch_start_sol + 0.5

                    times.append( glitch_times )
                    datas.append( [glitch_data_U, glitch_data_V, glitch_data_W] )
                    comps.append( glitch_comps )

                    sols_range.append(glitch_start_sol)


                # case glitch starts before midnight but ends after midnight
                else:
                    counter_counter += 1

                    glitch_comps = ''                  
                    if glitch_U == '1':
                        glitch_comps += 'U'
                    if glitch_V == '1':
                        glitch_comps += 'V'
                    if glitch_W == '1':
                        glitch_comps += 'W'

                    glitches_times_sol = np.array( [glitch_time_obj.sol for glitch_time_obj in glitch_times_obj] )
                    indices_sol_start  = np.where(glitches_times_sol==glitch_start_sol)

                    glitch_times       = [glitch_time_obj.LMST_time.datetime-datetime.timedelta(days=glitch_start_sol) for glitch_time_obj in glitch_times_obj[indices_sol_start]]
                    
                    glitch_data_U      = stream_glitch.select(component='U')[0]
                    glitch_data_U.detrend('demean')
                    glitch_data_U      = glitch_data_U.data[indices_sol_start]
                    glitch_data_U      = normalise(glitch_data_U, scale_to_between=[-np.max(abs(glitch_data_U))/glitch_max_amp*0.5*scale, np.max(abs(glitch_data_U))/glitch_max_amp*0.5*scale])
                    glitch_data_U      = glitch_data_U + glitch_start_sol + 0.5
                    glitch_data_V      = stream_glitch.select(component='U')[0]
                    glitch_data_V.detrend('demean')
                    glitch_data_V      = glitch_data_V.data[indices_sol_start]
                    glitch_data_V      = normalise(glitch_data_V, scale_to_between=[-np.max(abs(glitch_data_V))/glitch_max_amp*0.5*scale, np.max(abs(glitch_data_V))/glitch_max_amp*0.5*scale])
                    glitch_data_V      = glitch_data_V + glitch_start_sol + 0.5
                    glitch_data_W      = stream_glitch.select(component='U')[0]
                    glitch_data_W.detrend('demean')
                    glitch_data_W      = glitch_data_W.data[indices_sol_start]
                    glitch_data_W      = normalise(glitch_data_W, scale_to_between=[-np.max(abs(glitch_data_W))/glitch_max_amp*0.5*scale, np.max(abs(glitch_data_W))/glitch_max_amp*0.5*scale])
                    glitch_data_W      = glitch_data_W + glitch_start_sol + 0.5
                    
                    times.append( glitch_times )
                    datas.append( [glitch_data_U, glitch_data_V, glitch_data_W] )
                    comps.append( glitch_comps )
                    
                    
                    indices_sol_end    = np.where(glitches_times_sol==glitch_end_sol)
                    glitch_times       = [glitch_time_obj.LMST_time.datetime-datetime.timedelta(days=glitch_start_sol) for glitch_time_obj in glitch_times_obj[indices_sol_end]]
                    
                    glitch_data_U      = stream_glitch.select(component='U')[0]
                    glitch_data_U.detrend('demean')
                    glitch_data_U      = glitch_data_U.data[indices_sol_end]
                    glitch_data_U      = normalise(glitch_data_U, scale_to_between=[-np.max(abs(glitch_data_U))/glitch_max_amp*0.5*scale, np.max(abs(glitch_data_U))/glitch_max_amp*0.5*scale])
                    glitch_data_U      = glitch_data_U + glitch_start_sol + 0.5
                    glitch_data_V      = stream_glitch.select(component='U')[0]
                    glitch_data_V.detrend('demean')
                    glitch_data_V      = glitch_data_V.data[indices_sol_end]
                    glitch_data_V      = normalise(glitch_data_V, scale_to_between=[-np.max(abs(glitch_data_V))/glitch_max_amp*0.5*scale, np.max(abs(glitch_data_V))/glitch_max_amp*0.5*scale])
                    glitch_data_V      = glitch_data_V + glitch_start_sol + 0.5
                    glitch_data_W      = stream_glitch.select(component='U')[0]
                    glitch_data_W.detrend('demean')
                    glitch_data_W      = glitch_data_W.data[indices_sol_end]
                    glitch_data_W      = normalise(glitch_data_W, scale_to_between=[-np.max(abs(glitch_data_W))/glitch_max_amp*0.5*scale, np.max(abs(glitch_data_W))/glitch_max_amp*0.5*scale])
                    glitch_data_W      = glitch_data_W + glitch_start_sol + 0.5

                    times.append( glitch_times )
                    datas.append( [glitch_data_U, glitch_data_V, glitch_data_W] )
                    comps.append( glitch_comps )

                    sols_range.append(glitch_start_sol)
                    sols_range.append(glitch_end_sol)



    ### OUTPUT
    print()
    print(u'Detected sols: %s .. %s '         % (sols_range[0], sols_range[-1]))
    print(u'Found waveforms for %s glitches.' % (len(comps)-counter_counter))
    print()
    print(u'Done in:   %s (h:m:s).'           % sec2hms( time.time()-now ))
    print(u'Timestamp: %s'                    % datetime.datetime.now().strftime('%Y-%m-%dT%H:%M:%S'))



    ### PLOT
    sols_range     = sorted( set( sols_range ))
    titles         = ['U-component','V-component','W-component']

    fig, axes = plt.subplots(1, 3, figsize=(15, 7), sharex=True, sharey=True)
    fig.canvas.set_window_title('Glitch waveform plot')
    fig.suptitle('RAW glitch waveform plotter: LMST %s-%s, %s glitches (black=glitch, gray=not)' % (LMST_range[0], LMST_range[1], (len(comps)-counter_counter)), fontsize=11, y=0.99)
    fig.subplots_adjust(wspace=0.4, hspace=0.4)
    fig.canvas.mpl_connect('button_press_event', onclick)
    fig.autofmt_xdate()

    for l, ax in enumerate(axes):

        ax.set_title(titles[l], size=10)
        ax.autoscale(enable=True, axis='both', tight=False)
        ax.grid(ls='-.', lw=0.5, zorder=1)
        ax.set_xlabel('LMST')
        ax.xaxis.labelpad = 3
        ax.xaxis.set_major_formatter( mdates.DateFormatter('%H:%M:%S') )
        ax.xaxis.set_major_locator(plt.MaxNLocator(5))
        #ax.xaxis.set_major_locator( mdates.HourLocator(interval=4) )
        ax.set_ylabel('Sol')
        ax.set_ylim([min(sols_range),max(sols_range)+1])

        # data
        for m in range(len(comps)):

            if 'UVW'[l] in comps[m]:
                colour = 'k'
            else:
                colour = 'darkgray'

            xdata = times[m]
            ydata = datas[m][l]
            ax.plot(xdata, ydata, lw=0.75, c=colour, zorder=3)


    # gray lines middle of sol (base line)
    xlim = ax.get_xlim()
    for l, ax in enumerate(axes):
        ax.set_xlim( xlim )     
        for sol in sols_range[::-1]:
            ax.plot( xlim, [sol+0.5, sol+0.5], lw=0.5, c='gray', zorder=2)

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    if outfile:
        plt.savefig(outfile)
        print(outfile)
def plot_glitch_align(*glitch_files, run=True, waveform_files=[], sols=[], LMST_range=['0', '24'], min_amp=None, max_amp=None, mode='largest', scale=1, align='maximum', outfile='', **kwargs):

    """

    """



    ### RUN OR NOT:
    if not run:
        return

    now = time.time()



    ### HELPER FUNCTION
    def onpick(event):
        index  = int(round( np.average( event.artist.get_ydata() ))) - 1
        glitch = all_glitches[glitches_found[index]]
        print( marstime(glitch[1] ))



    ### OUTPUT
    print()
    print(u'  ---------------------')
    print(u'  RUNNING ALIGN PLOTTER')
    print(u'  ---------------------')
    print()



    ### SANITY CHECK:
    if not waveform_files:
        print()
        print(u'ERROR: Cannot plot without `waveform_files`.')
        sys.exit()



    ### READ GLITCH-FILE
    all_glitches = []

    for glitch_file in glitch_files:
        glitches      = np.loadtxt(glitch_file, dtype='str')
        #glitches      = glitch_exclude(glitches, verbose=False)
        all_glitches += list(glitches)
    all_glitches = np.array(all_glitches)



    ### SELECT GLITCHES ON GIVEN `SOLS`
    if sols:
        sols = sorted(sols)
        all_glitches = np.array( [glitch for glitch in all_glitches if marstime(glitch[1]).sol>=sols[0] and marstime(glitch[1]).sol<=sols[1]] )



    ### SELECT GLITCHES IN SPECIFIED LMST RANGE
    LMST_range        = [(str(t).replace(':','') + '000000')[:6] for t in LMST_range]
    LMST_range        = [[t[i:i+2] for i in range(0,len(t),2)] for t in LMST_range]
    LMST_range        = [':'.join(t) for t in LMST_range]
    
    glitch_start_LMST = np.array( [glitch[3].split('M')[1] for glitch in all_glitches] )
    if LMST_range[0]<=LMST_range[1]:
        indices       = np.where( (glitch_start_LMST >= LMST_range[0]) & (glitch_start_LMST <= LMST_range[1]) )
    else:
        indices       = np.where( (glitch_start_LMST >= LMST_range[0]) | (glitch_start_LMST <= LMST_range[1]) )
    all_glitches      = all_glitches[ indices ]



    ### SELECT GLITCHES TO `MIN_AMP` AND `MAX_AMP`
    if eval(str(min_amp).capitalize()) is None or eval(str(min_amp).capitalize()) == False:
        min_amp = 0
    else:
        min_amp = float( min_amp )
    if eval(str(max_amp).capitalize()) is None or eval(str(max_amp).capitalize()) == False:
        max_amp = 1e9
    else:
        max_amp = float( max_amp )   

    all_glitches = all_glitches[((abs(all_glitches[:,8].astype('float'))>=min_amp)   | \
                                 (abs(all_glitches[:,9].astype('float'))>=min_amp)   | \
                                 (abs(all_glitches[:,10].astype('float'))>=min_amp)) & \

                                 (abs(all_glitches[:,8].astype('float'))<=max_amp)   & \
                                 (abs(all_glitches[:,9].astype('float'))<=max_amp)   & \
                                 (abs(all_glitches[:,10].astype('float'))<=max_amp)]



    ### SELECT GLITCHES ACCORDING TO SPECIFIED MODE
    if mode.lower()=='all':
        pass

    elif mode.lower()=='largest':
        all_glitches_new = []

        times            = np.array( [marstime(time) for time in all_glitches[:,1]] )
        sols             = sorted( set( [time.sol for time in times] ))

        for sol in sols:

            indices  = np.where( np.array([time.sol for time in times]) == sol)
            max_amps = []
            for glitch in all_glitches[indices]:

                U = np.abs(float(glitch[8]))
                V = np.abs(float(glitch[9]))
                W = np.abs(float(glitch[10]))
                max_amps.append( np.max( [U, V, W] ))

            max_index = np.argmax( max_amps )
            all_glitches_new.append( all_glitches[indices][max_index] )

        all_glitches = np.array(all_glitches_new)

    else:
        pass



    ### VARIABLES
    streams        = [read(file) for file in waveform_files]
    glitch_max_amp = np.max( np.abs( all_glitches[:,8:11].astype('float') ))*1
    glitches_found = []
    comps          = []
    times          = []
    datas          = []
    container      = []
    counter        = 0



    ### LOOP OVER GLITCHES, PREPARE DATE FOR PLOTTING
    for i, glitch in enumerate(all_glitches):

        if (i+1)%10 == 0:
            print(u'Handling glitch %06d / %06d' % (i+1, len(all_glitches)))

        glitches_start_UTC         = UTCDateTime(glitch[1])
        glitches_end_UTC           = UTCDateTime(glitch[2])
        glitch_U                   = glitch[5]
        glitch_V                   = glitch[6]
        glitch_W                   = glitch[7]


        for m, stream in enumerate(streams):

            stream_glitch = stream.slice(starttime=glitches_start_UTC, endtime=glitches_end_UTC)

            if stream_glitch:

                counter += 1
                glitches_found.append(i)

                glitch_times_UTC = stream_glitch[0].times('relative')

                glitch_comps = ''                  
                if glitch_U == '1':
                    glitch_comps += 'U'
                if glitch_V == '1':
                    glitch_comps += 'V'
                if glitch_W == '1':
                    glitch_comps += 'W'

                if 'maximum' in align.lower() or 'max' in align.lower():
                    U  = stream_glitch.select(component='U')[0]
                    V  = stream_glitch.select(component='V')[0]
                    W  = stream_glitch.select(component='W')[0]

                    if 'u' in align.lower():
                        index_comp  = 0
                    if 'v' in align.lower():
                        index_comp  = 1                    
                    if 'w' in align.lower():
                        index_comp  = 2
                    else:
                        index_comp = np.argmax( [max(U), max(U), max(W)] )

                    offset       = int(2*U.stats.sampling_rate)
                    index_shift  = np.argmax( np.abs( [U.data, V.data, W.data][index_comp][offset:] ))
                    glitch_times = stream_glitch[0].times('relative') - (index_shift+offset) * [U, V, W][index_comp].stats.delta
                else:
                    glitch_times = stream_glitch[0].times('relative')
                
                glitch_data_U = stream_glitch.select(component='U')[0]
                glitch_data_U.detrend('demean')
                glitch_data_U = glitch_data_U.data
                glitch_data_U = normalise(glitch_data_U, scale_to_between=[np.min(glitch_data_U)/glitch_max_amp*0.5*scale, np.max(glitch_data_U)/glitch_max_amp*0.5*scale])
                glitch_data_U = glitch_data_U + counter
                glitch_data_V = stream_glitch.select(component='V')[0]
                glitch_data_V.detrend('demean')
                glitch_data_V = glitch_data_V.data
                glitch_data_V = normalise(glitch_data_V, scale_to_between=[np.min(glitch_data_V)/glitch_max_amp*0.5*scale, np.max(glitch_data_V)/glitch_max_amp*0.5*scale])
                glitch_data_V = glitch_data_V + counter
                glitch_data_W = stream_glitch.select(component='W')[0]
                glitch_data_W.detrend('demean')
                glitch_data_W = glitch_data_W.data
                glitch_data_W = normalise(glitch_data_W, scale_to_between=[np.min(glitch_data_W)/glitch_max_amp*0.5*scale, np.max(glitch_data_W)/glitch_max_amp*0.5*scale])
                glitch_data_W = glitch_data_W + counter

                times.append( glitch_times )
                datas.append( [glitch_data_U, glitch_data_V, glitch_data_W] )
                comps.append( glitch_comps )


    ### OUTPUT
    print()
    print(u'Found waveforms for %s glitches.' % counter)
    print()
    print(u'Done in:   %s (h:m:s).'           % sec2hms( time.time()-now ))
    print(u'Timestamp: %s'                    % datetime.datetime.now().strftime('%Y-%m-%dT%H:%M:%S'))



    ### PLOT
    titles    = ['U-component','V-component','W-component']

    fig, axes = plt.subplots(1, 3, figsize=(15, 7), sharex=True, sharey=True)
    fig.canvas.set_window_title('Glitch align plot')
    fig.suptitle('RAW glitch align plotter: LMST %s-%s, %s glitches (black=glitch, gray=not), mode=%s, align=%s' % (LMST_range[0], LMST_range[1], len(comps), mode, align), fontsize=11, y=0.99)
    fig.subplots_adjust(wspace=0.4, hspace=0.4)
    fig.canvas.mpl_connect('pick_event', onpick)

    for l, ax in enumerate(axes):

        ax.set_title(titles[l], size=10)
        ax.set_xlabel('Relative time of glitch / (s)')
        ax.xaxis.labelpad = 3
        ax.yaxis.set_major_locator(mticker.MaxNLocator(integer=True))
        ax.set_ylabel('Number of glitch')

        # data
        for m in range(len(comps)):

            if 'UVW'[l] in comps[m]:
                colour = 'k'
            else:
                colour = 'darkgray'

            xdata     = times[m]
            ydata     = datas[m][l]
            ax.plot(xdata, ydata, lw=0.75, c=colour, zorder=3, picker=3)

            ydata_red = ydata - (m+3)
            ax.plot(xdata, ydata_red, lw=0.75, c='lightcoral', alpha=0.3, zorder=3)

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    if outfile:
        plt.savefig(outfile)
        print(outfile)


### _ _ N A M E _ _ = = " _ _ M A I N _ _ "  
if __name__ == "__main__":
    pass