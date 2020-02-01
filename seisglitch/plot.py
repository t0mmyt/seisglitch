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
from seisglitch.util import marstime, sec2hms
from seisglitch.math import normalise



### PLOT SCRIPTS TO VISUALIZE
# first 6 indices are for UVWZNE amplitudes, 7th for backazimuth angle, 8th for incidence angle 
# (w.r.t. glitch list as created by glitch_detector)
columns = {'GAI' : [11,12,13,14,15,16,35,36],
           'DIS' : [17,18,19,20,21,22,39,40],
           'VEL' : [23,24,25,26,27,28,43,44],
           'ACC' : [29,30,31,32,33,34,47,48]}
def glitch_overview_plot(*glitch_files, run=True, UNIT='GAI', sols=None, LMST_range=['0', '24'], min_amp=None, max_amp=None, show=True, outfile=None):
    
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
        nonlocal picked_already
        if not picked_already:
            picked_already = True
            header         = u'#               UTC          LMST      U-{0}      V-{0}      W-{0}      Z-{0}      N-{0}      E-{0}  BAZ-{0}  INC-{0}'.format(UNIT)
            print(header)

        indices = event.ind
        for index in indices:
            print(u'%s  %s  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %7.1f  %7.1f' % (glitch_starts[index].UTC_time.strftime('%Y-%m-%dT%H:%M:%S'), 
                                                                                       glitch_starts[index].LMST, 
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
    print(u'  RUNNING GLITCH OVERVIEW PLOTTER')
    print(u'  -------------------------------')
    print()



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
    
    glitch_start_LMST = np.array( [glitch[3].split('S')[1] for glitch in all_glitches] )
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

    all_glitches = all_glitches[((abs(all_glitches[:,columns[UNIT][0]].astype('float'))>=min_amp)  | \
                                 (abs(all_glitches[:,columns[UNIT][1]].astype('float'))>=min_amp)  | \
                                 (abs(all_glitches[:,columns[UNIT][2]].astype('float'))>=min_amp)  | \
                                 (abs(all_glitches[:,columns[UNIT][3]].astype('float'))>=min_amp)  | \
                                 (abs(all_glitches[:,columns[UNIT][4]].astype('float'))>=min_amp)  | \
                                 (abs(all_glitches[:,columns[UNIT][5]].astype('float'))>=min_amp)) & \
   
                                 (abs(all_glitches[:,columns[UNIT][0]].astype('float'))<=max_amp)  & \
                                 (abs(all_glitches[:,columns[UNIT][1]].astype('float'))<=max_amp)  & \
                                 (abs(all_glitches[:,columns[UNIT][2]].astype('float'))<=max_amp)  & \
                                 (abs(all_glitches[:,columns[UNIT][3]].astype('float'))<=max_amp)  & \
                                 (abs(all_glitches[:,columns[UNIT][4]].astype('float'))<=max_amp)  & \
                                 (abs(all_glitches[:,columns[UNIT][5]].astype('float'))<=max_amp)]


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



    ### FIGURE
    fig = plt.figure(figsize=(9,11))
    fig.canvas.set_window_title('Glitch overview plot')
    fig.suptitle('%s overview plotter: LMST %s-%s, %s glitches (UVW: %g ≤ 1comp ≤ %g)' % (UNIT, LMST_range[0], LMST_range[1], len(all_glitches), min_amp, max_amp), fontsize=11, y=0.99)
    fig.subplots_adjust(wspace=0.4, hspace=0.4)
    fig.canvas.mpl_connect('pick_event', onpick)



    ## COLORMAP
    cmap      = plt.get_cmap('hsv')
    #norm      = Normalize(vmin=0, vmax=24)
    norm      = mpl.colors.BoundaryNorm(np.linspace(0,24,25), cmap.N)
    scalarMap = cmx.ScalarMappable(norm=norm, cmap=cmap)
    colours   = [scalarMap.to_rgba(glitch_starts[i].LMST_time.hour) for i in range(len(glitch_starts))]

    cax = fig.add_axes([0.01, 0.55, 0.01, 0.4]) # left, bottom, width, height
    cax.yaxis.set_ticks_position('left')
    cb1 = mpl.colorbar.ColorbarBase(cax, drawedges=True, cmap=cmap, norm=norm, orientation='vertical', ticks=np.linspace(0,24,9), boundaries=np.linspace(0,24,25))
    cb1.set_label('LMST hour')

    # create hist data and colours
    glitch_start_hist_SOL_mat = [[] for _ in range(24)]
    for start in glitch_starts:
        glitch_start_hist_SOL_mat[start.LMST_time.hour] += [start.LMST_time.julday]
    if np.array(glitch_start_hist_SOL_mat).any():
        colors_hist = [scalarMap.to_rgba(i) for i in range(len(glitch_start_hist_SOL_mat))]
    else:   # if no glitches, make it still work for plotting
        colors_hist = [(0,0,0,0)]



    ## HORIZONTAL POLARIZATIONS (upper left)
    sizes = np.sqrt( N_amps**2 + E_amps**2 )
    sizes = normalise(sizes, scale_to_between=[5,20])

    ax   = fig.add_subplot(221, polar=True)
    ax.set_title('Back-azimuths', pad=20, size=10)
    ax.spines['polar'].set_visible(False)
    ax.tick_params(pad=5)

    ax.set_theta_direction(-1)
    ax.set_theta_zero_location('N')
    ax.set_thetagrids([0,15,37,70,90,105,114,135,157,180,225,251,255,270,277,345.5], labels=['','LSA/VBBV\n(15°)','LVL1\n(37°)', r'HP$^{3}$'+'\n(70°)','','SP2\n(105°)','WTS E\n(114°)',r'VBBU/Y$_{asm}$'+'\n(135°)','LVL2\n(157°)','',r'X$_{asm}$','WTS W (251°)','VBBW (255°)','','LVL3\n(277°)','SP3/WTS N\n(345°)'], size=8)

    #ax.set_rgrids([0.1,0.4,0.7,1], labels=['Hmin','','','Hmax'], size=6)
    ax.set_ylim( [min(sols_range)-0.1*(max(sols_range)-min(sols_range)), max(sols_range)+0.1*(max(sols_range)-min(sols_range))] )
    ax.tick_params(axis='y', labelsize=6)   # because with thetamin/max command, fontzise in set_rgrids doesn't work ...
    ax.set_yticklabels(['Sol %d' % x for x in ax.get_yticks()])
    ax.set_rlabel_position(315)

    ax.scatter(phis, sols, s=sizes, linewidths=.3, c=colours, edgecolor='k', rasterized=True, picker=2)
    ax.grid(ls='-.', lw=0.4, c='dimgray')



    ## VERTICAL POLARIZATIONS (upper right)
    sizes = np.sqrt(Z_amps**2)
    sizes = normalise(sizes , scale_to_between=[5,20])

    ax = fig.add_subplot(222, polar=True)
    ax.set_title('Incidence angles', pad=20, size=10)
    ax.spines['polar'].set_visible(False)
    ax.tick_params(pad=5)

    ax.set_theta_direction(-1)
    ax.set_theta_zero_location('N')
    ax.set_thetamin(0)
    ax.set_thetamax(180)
    ax.set_thetagrids([0,45,60.5,90,135,180], labels=['0°','45°','VBB\n(60.5°)','90°','135°','180°'], size=8)

    #ax.set_rgrids([0.1,0.4,0.7,1], labels=['Zmin','','','Zmax'], fontsize=6)
    ax.set_ylim( [min(sols_range)-0.1*(max(sols_range)-min(sols_range)), max(sols_range)+0.1*(max(sols_range)-min(sols_range))] )
    ax.tick_params(axis='y',labelsize=6)    # because with thetamin/max command, fontzise in set_rgrids doesn't work ...
    ax.set_yticklabels(['Sol %d' % x for x in ax.get_yticks()])
    ax.set_rlabel_position(315)

    ax.scatter(incs, sols, s=sizes, color=colours, edgecolor='k', linewidths=.3, rasterized=True, picker=2)
    ax.grid(ls='-.', lw=0.4, c='dimgray')



    ## 3-D PLOT (lower left)
    sizes = np.sqrt(Z_amps**2+N_amps**2+E_amps**2) 
    sizes = normalise(sizes, scale_to_between=[5,20])

    ax = fig.add_subplot(223, projection='3d')
    ax.set_title('3-D Polarizations', y=1.10, size=10)
    ax.view_init(elev=15., azim=-50)
    ax.set_xlabel('E', labelpad=2, fontsize=8)
    ax.set_ylabel('N', labelpad=2, fontsize=8)
    ax.set_zlabel('Z', labelpad=2, fontsize=8)
    ax.set_xticks([-1,-0.5,0,0.5,1])
    ax.set_yticks([-1,-0.5,0,0.5,1])
    ax.set_zticks([-1,-0.5,0,0.5,1])
    ax.scatter(np.sin(phis)*abs(np.sin(incs)), np.cos(phis)*abs(np.sin(incs)), -np.cos(incs), s=sizes, color=colours, edgecolor='k', linewidths=.3, depthshade=False, rasterized=True)
    ax.scatter(0, 0, 0, s=30, c='white', depthshade=False)
    ax.set_xlim( [-1,1] )
    ax.set_ylim( [-1,1] )
    ax.set_zlim( [-1,1] )



    ## HISTOGRAM (lower right)
    ax = fig.add_subplot(224)
    #ax.set_title('Sol Histogram Glitches', pad=10)
    ax.grid(ls='-.', lw=0.5, zorder=0)
    ax.set(xlabel='Sol', ylabel='Glitches / Sol')
    ax.autoscale(enable=True, axis='both', tight=False)
    ax.hist(glitch_start_hist_SOL_mat, bins=np.arange(min(sols_range), max(sols_range)+2,1), color=colors_hist, align='left', lw=0.5, stacked=True, zorder=3, rasterized=True) # "+2" instead of "+1"  because otherwise there is a bug in the last bin
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
        plt.savefig(outfile)
        print(outfile)
    if show:
        plt.show()
    plt.close()
def glitch_gutenberg_plot(*glitch_files, run=True, UNIT='GAI', sols=None, LMST_range=['0', '24'], show=True, outfile=None):

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
    
    glitch_start_LMST = np.array( [glitch[3].split('S')[1] for glitch in all_glitches] )
    if LMST_range[0]<=LMST_range[1]:
        indices       = np.where( (glitch_start_LMST >= LMST_range[0]) & (glitch_start_LMST <= LMST_range[1]) )
    else:
        indices       = np.where( (glitch_start_LMST >= LMST_range[0]) | (glitch_start_LMST <= LMST_range[1]) )
    all_glitches      = all_glitches[ indices ]



    ### SOME VARIABLE DECLARATIONS
    label_comps = 'UVWZNE'
    range_hist  = (1e-10,1e-6)



    ### PLOT HISTOGRAMS
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
        ax.set_ylim(1,1e4)
        ax.xaxis.set_tick_params(labelbottom=True)
        ax.yaxis.set_tick_params(labelbottom=True)

        # create bins, but only needed for first one
        if l==0:
            hist, bins = np.histogram( np.abs(all_glitches[:,columns[UNIT][l]].astype('float')), bins=25)
            logbins    = np.logspace(np.log10(bins[0]), np.log10(bins[-1]),len(bins))

        # plot histogram
        y, x, _         = ax.hist( np.abs(all_glitches[:,columns[UNIT][l]].astype('float')), log=True, range=range_hist, color='k', bins=logbins, ec='w', lw=0.5, label=label_comps[l]+'-comp')
        
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
    if show:
        plt.show()
    plt.close()
def glitch_waveform_plot(*glitch_files, run=True, waveform_files=[], sols=None, LMST_range=['0', '24'], min_amp=None, max_amp=None, scale=1, show=True, outfile=None):



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



    ### CHECK IF WAVEFORM FILES PASSED
    if not waveform_files:
        print('ERROR: You need to specify `waveform_files` in the config.yaml so waveform data can be read.')
        sys.exit()



    ### OUTPUT
    print()
    print(u'  -------------------------------')
    print(u'  RUNNING GLITCH WAVEFORM PLOTTER')
    print(u'  -------------------------------')
    print()
    print(u"Have to convert each data point's UTC time to LMST, this takes a bit ..")
    print()



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
    
    glitch_start_LMST = np.array( [glitch[3].split('S')[1] for glitch in all_glitches] )
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

                                 (abs(all_glitches[:,8].astype('float'))<=max_amp)  & \
                                 (abs(all_glitches[:,9].astype('float'))<=max_amp)  & \
                                 (abs(all_glitches[:,10].astype('float'))<=max_amp)]



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
        glitch_start_sol           = int(glitch[3].split('S')[0])
        glitch_end_sol             = int(glitch[4].split('S')[0])
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
        ax.set_ylabel('Sols')
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
    if show:
        plt.show()
    plt.close()
def glitch_XoverBAZ_plot(*glitch_files, run=True, comp='Z', UNIT='GAI', sols=None, LMST_range=['0', '24'], min_amp=None, max_amp=None, show=True, outfile=None):

    """

    """



    ### RUN OR NOT:
    if not run:
        return



    ### OUTPUT
    print()
    print(u'  -------------------------------')
    print(u'  RUNNING GLITCH %soverBAZ PLOTTER' % comp)
    print(u'  -------------------------------')
    print()



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
    
    glitch_start_LMST = np.array( [glitch[3].split('S')[1] for glitch in all_glitches] )
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

    all_glitches = all_glitches[((abs(all_glitches[:,columns[UNIT][0]].astype('float'))>=min_amp)  | \
                                 (abs(all_glitches[:,columns[UNIT][1]].astype('float'))>=min_amp)  | \
                                 (abs(all_glitches[:,columns[UNIT][2]].astype('float'))>=min_amp)  | \
                                 (abs(all_glitches[:,columns[UNIT][3]].astype('float'))>=min_amp)  | \
                                 (abs(all_glitches[:,columns[UNIT][4]].astype('float'))>=min_amp)  | \
                                 (abs(all_glitches[:,columns[UNIT][5]].astype('float'))>=min_amp)) & \

                                 (abs(all_glitches[:,columns[UNIT][0]].astype('float'))<=max_amp)  & \
                                 (abs(all_glitches[:,columns[UNIT][1]].astype('float'))<=max_amp)  & \
                                 (abs(all_glitches[:,columns[UNIT][2]].astype('float'))<=max_amp)  & \
                                 (abs(all_glitches[:,columns[UNIT][3]].astype('float'))<=max_amp)  & \
                                 (abs(all_glitches[:,columns[UNIT][4]].astype('float'))<=max_amp)  & \
                                 (abs(all_glitches[:,columns[UNIT][5]].astype('float'))<=max_amp)]



    ### SOME VARIABLE DECLARATIONS
    glitch_starts = np.array( [marstime(e).LMST_time.hour for e in all_glitches[:,1]] )
    U             = abs(all_glitches[:,columns[UNIT][0]].astype('float'))
    V             = abs(all_glitches[:,columns[UNIT][1]].astype('float'))
    W             = abs(all_glitches[:,columns[UNIT][2]].astype('float'))
    Z             = abs(all_glitches[:,columns[UNIT][3]].astype('float'))
    N             = abs(all_glitches[:,columns[UNIT][4]].astype('float'))
    E             = abs(all_glitches[:,columns[UNIT][5]].astype('float'))
    BAZ           = abs(all_glitches[:,columns[UNIT][6]].astype('float'))*np.pi/180



    ### FIGURE
    fig = plt.figure(figsize=(8,8))
    fig.canvas.set_window_title('Glitch %soverBaz plot' % comp)
    fig.suptitle('%s %soverBaz plotter: LMST %s-%s, %s glitches (UVW: %g ≤ 1comp ≤ %g)' % (UNIT, comp, LMST_range[0], LMST_range[1], len(all_glitches), min_amp, max_amp), fontsize=11, y=0.99)
    fig.subplots_adjust(wspace=0.4, hspace=0.4)
    

    ## COLORMAP
    cmap      = plt.get_cmap('hsv')
    #norm      = Normalize(vmin=0, vmax=24)
    norm      = mpl.colors.BoundaryNorm(np.linspace(0,24,25), cmap.N)
    scalarMap = cmx.ScalarMappable(norm=norm, cmap=cmap)
    colours   = [scalarMap.to_rgba(time) for time in glitch_starts]

    cax = fig.add_axes([0.01, 0.55, 0.01, 0.4]) # left, bottom, width, height
    cax.yaxis.set_ticks_position('left')
    cb1 = mpl.colorbar.ColorbarBase(cax, drawedges=True, cmap=cmap, norm=norm, orientation='vertical', ticks=np.linspace(0,24,9), boundaries=np.linspace(0,24,25))
    cb1.set_label('LMST hour')


    ## POLAR PLOT
    ax = fig.add_axes((0.1, 0.15, 0.8, 0.7), projection='polar')
    ax.spines['polar'].set_visible(False)
    ax.tick_params(pad=5)
    ax.grid(ls='-.', lw=0.4, c='dimgray')
    ax.set_theta_direction(-1)
    ax.set_theta_zero_location('N')
    ax.set_thetagrids([0,15,37,70,90,105,114,135,157,180,225,251,255,270,277,345.5], labels=['','LSA/VBBV\n(15°)','LVL1\n(37°)', r'HP$^{3}$'+'\n(70°)','','SP2\n(105°)','WTS E\n(114°)',r'VBBU/Y$_{asm}$'+'\n(135°)','LVL2\n(157°)','',r'X$_{asm}$','WTS W (251°)','VBBW (255°)','','LVL3\n(277°)','SP3/WTS N\n(345°)'], size=8)
    ax.tick_params(axis='y', labelsize=6)   # because with thetamin/max command, fontzise in set_rgrids doesn't work ...
    #ax.set_yticklabels(['%d' % x for x in ax.get_yticks()])
    ax.set_ylim([-1e-8,1.25e-7])
    ax.set_rlabel_position(315)
    ax.scatter(BAZ, eval(comp.upper()), s=10, c=colours, edgecolor='k', linewidths=.3)



    ### FINAL
    if outfile:
        plt.savefig(outfile)
        print(outfile)
    if show:
        plt.show()
    plt.close()
def glitch_align_plot(*glitch_files, run=True, waveform_files=[], sols=None, LMST_range=['0', '24'], min_amp=None, max_amp=None, mode='largest', scale=1, align='maximum', show=True, outfile=None):

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
    
    glitch_start_LMST = np.array( [glitch[3].split('S')[1] for glitch in all_glitches] )
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


        for stream in streams:

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
    if show:
        plt.show()
    plt.close()



### _ _ N A M E _ _ = = " _ _ M A I N _ _ "  
if __name__ == "__main__":
    pass