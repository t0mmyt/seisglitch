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
from obspy import UTCDateTime


#####  seisglitch util import  #####
from seisglitch.util import normalise, solify, UTCify, ptime



### PLOT SCRIPTS TO VISUALIZE
# first 6 indices are for UVWZNE amplitudes, 7th for backazimuth angle, 8th for incidence angle 
# (w.r.t. glitch list as created by glitch_detector)
columns = {'GAI' : [9, 10,11,12,13,14,33,34],
           'DIS' : [15,16,17,18,19,20,37,38],
           'VEL' : [21,22,23,24,25,26,41,42],
           'ACC' : [27,28,29,30,31,31,45,46]}
def glitch_overview_plot(*glitch_files, LMST_hour=None, min_amp=1e-10, max_amp=1e-5, UNIT='DIS', show=True, outfile=None):
    
    """
    Plot glitches, based on glitch file produced by function `glitch_detector()`.

    amp parameters apply to all components: UVWZNE
    min_amp: at least one component amplitude must be larger
    max_amp: all component amplitudes must be smaller

    UNIT allows to choose which polarizations you wish to plot.

    Matplotlib colorscale:
    https://matplotlib.org/3.1.1/tutorials/colors/colormaps.html
    """


    ### PICK FUNCTION
    picked_already = False
    def onpick(event):
        nonlocal picked_already
        if not picked_already:
            picked_already = True
            header         = u'#               UTC          LMST      U-{0}      V-{0}      W-{0}      Z-{0}      N-{0}      E-{0}  BAZ-{0}  INC-{0}'.format(UNIT)
            print(header)

        indices = event.ind
        for index in indices:
            print(u'%s  %sS%s  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %7.1f  %7.1f' % (glitch_starts_UTC[index].strftime('%Y-%m-%dT%H:%M:%S'), 
                                                                                          glitch_starts_SOL[index].julday, 
                                                                                          glitch_starts_SOL[index].strftime('%H:%M:%S'),
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



    ### READ GLITCH-FILE
    all_glitches = []

    for glitch_file in glitch_files:
        glitches      = np.loadtxt(glitch_file, dtype='str')
        #glitches      = glitch_exclude(glitches, verbose=False)
        all_glitches += list(glitches)

    all_glitches = np.array(all_glitches)
    if LMST_hour or LMST_hour==0:
        all_glitches = np.array( [e for e in all_glitches if solify(UTCDateTime(e[1])).hour == LMST_hour] )



    ### PREPARE DATA
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
    sols_range   = np.array( [solify(UTCDateTime(e[1])).julday for e in all_glitches]  )
    ymax_hist    = Counter( [solify(UTCDateTime(e[1])).julday for e in all_glitches] ).most_common()[0][1]



    ### LMST HOUR EXCLUSION, if wished
    if LMST_hour != False and LMST_hour!=None:
        LMST_hours   = np.arange(float(str(LMST_hour).split('-')[0]), float(str(LMST_hour).split('-')[-1])+1)
        all_glitches = np.array( [e for e in all_glitches if solify(UTCDateTime(e[1])).hour in LMST_hours] )



    ### ASSIGN NEEDED VARIABLES
    try:
        glitch_starts_UTC = np.array( [UTCDateTime(e) for e in all_glitches[:,1]] )
        glitch_starts_SOL = np.array( [solify(e) for e in glitch_starts_UTC] )

        Z_amps            = np.array( [                float(e)  for e in all_glitches[:,columns[UNIT][3]]] )
        N_amps            = np.array( [                float(e)  for e in all_glitches[:,columns[UNIT][4]]] )
        E_amps            = np.array( [                float(e)  for e in all_glitches[:,columns[UNIT][5]]] )
        phis              = np.array( [    np.pi/180 * float(e)  for e in all_glitches[:,columns[UNIT][6]]] )
        incs              = np.array( [abs(np.pi/180 * float(e)) for e in all_glitches[:,columns[UNIT][7]]] )
        sols              = np.array( [solify(e).julday for e in glitch_starts_UTC] )
        

    except IndexError:  #no glitches for chosen conditions
        glitch_starts_UTC = np.array( [] )
        glitch_starts_SOL = np.array( [] )

        Z_amps            = np.array( [] )
        N_amps            = np.array( [] )
        E_amps            = np.array( [] )
        phis              = np.array( [] )
        incs              = np.array( [] )
        sols              = np.array( [] )



    ### FIGURE
    if LMST_hour  or LMST_hour==0:
        title = 'LMST %02s, %s, %s glitches (UVW: %g ≤ 1comp ≤ %g)' % (LMST_hour, UNIT, len(glitches), min_amp, max_amp)
    else:
        title = 'LMST 00:00-24:00, %s, %s glitches (UVW: %g ≤ 1comp ≤ %g)' % (UNIT, len(glitches), min_amp, max_amp)        
    fig = plt.figure(figsize=(10,10))
    fig.canvas.set_window_title('Glitch overview plot')
    fig.suptitle(title, fontsize=11, y=0.99)
    fig.subplots_adjust(wspace=0.4, hspace=0.4)
    fig.canvas.mpl_connect('pick_event', onpick)


    ## COLORMAP
    cmap      = plt.get_cmap('hsv')
    #norm      = Normalize(vmin=0, vmax=24)
    norm      = mpl.colors.BoundaryNorm(np.linspace(0,24,25), cmap.N)
    scalarMap = cmx.ScalarMappable(norm=norm, cmap=cmap)
    colours   = [scalarMap.to_rgba(glitch_starts_SOL[i].hour) for i in range(len(glitch_starts_SOL))]

    cax = fig.add_axes([0.01, 0.55, 0.01, 0.4]) # left, bottom, width, height
    cax.yaxis.set_ticks_position('left')
    cb1 = mpl.colorbar.ColorbarBase(cax, drawedges=True, cmap=cmap, norm=norm, orientation='vertical', ticks=np.linspace(0,24,9), boundaries=np.linspace(0,24,25))
    cb1.set_label('LMST hour')

    # create hist data and colours
    glitch_start_hist_SOL_mat = [[] for _ in range(24)]
    for start in glitch_starts_SOL:
        glitch_start_hist_SOL_mat[start.hour] += [start.julday]
    if np.array(glitch_start_hist_SOL_mat).any():
        colors_hist = [scalarMap.to_rgba(i) for i in range(len(glitch_start_hist_SOL_mat))]
    else:   # if no glitches, make it still work for plotting
        colors_hist = [(0,0,0,0)]

    if LMST_hour or LMST_hour==0:   # arrow next to colorbar
        ax = fig.add_axes([0.017, 0.55, 0.01, 0.4])
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        #ax.text(1, LMST_hour/24., r'$\rightarrow$' , ha='right', color='k', transform=ax.transAxes, fontsize=20)


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

    ax.scatter(phis, sols, s=sizes, edgecolor='k', linewidths=.3, c=colours, rasterized=True, picker=2)
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
    ax.scatter(np.sin(phis)*abs(np.sin(incs)), np.cos(phis)*abs(np.sin(incs)), np.cos(incs), s=sizes, color=colours, edgecolor='k', linewidths=.3, depthshade=False, rasterized=True)
    ax.scatter(0, 0, 0, s=30, c='white', edgecolors='k', depthshade=False)
    ax.set_xlim( [-1,1] )
    ax.set_ylim( [-1,1] )
    ax.set_zlim( [-1,1] )


    ## HISTOGRAM (lower right)
    ax = fig.add_subplot(224)
    #ax.set_title('Sol Histogram Glitches', pad=10)
    ax.grid(ls='-.', lw=0.5, zorder=0)
    ax.set(xlabel='Sol', ylabel='Glitches / Sol')
    ax.autoscale(enable=True, axis='both', tight=False)
    ax.hist(glitch_start_hist_SOL_mat, bins=np.arange(min(sols_range), max(sols_range)+2,1), color=colors_hist, align='left', ec='k', lw=0.5, stacked=True, zorder=3, rasterized=True) # "+2" instead of "+1"  because otherwise there is a bug in the last bin
    ax.set_xlim([min(sols_range), max(sols_range)])
    ax.set_ylim([0,ymax_hist])
    ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True, nbins='auto'))
    ax2   = ax.twiny()
    ax1Xs = ax.get_xticks()
    ax2Xs = [UTCify(UTCDateTime('1970-01-01T00:00:00.000000Z')+datetime.timedelta(days=e)).strftime('%m-%d') for e in ax1Xs]
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
def glitch_gutenberg_plot(*glitch_files, LMST_hour=None, UNIT='DIS', show=True, outfile=None):

    """

    """


    ### READ GLITCH-FILE
    all_glitches = []

    for glitch_file in glitch_files:
        glitches      = np.loadtxt(glitch_file, dtype='str')
        #glitches      = glitch_exclude(glitches, verbose=False)
        all_glitches += list(glitches)

    all_glitches = np.array(all_glitches)
    if LMST_hour or LMST_hour==0:
        all_glitches = np.array( [e for e in all_glitches if solify(UTCDateTime(e[1])).hour == LMST_hour] )



    ### SOME VARIABLES DECLARATION
    label_comps = 'UVWZNE'
    range_hist  = (1e-9,1e-5)

    if LMST_hour != False and LMST_hour!=None:
        title = 'LMST %02d:00-%02d:00, %s, %s glitches' % (LMST_hour, LMST_hour+1, UNIT, len(all_glitches))
    else:
        title = 'LMST 00:00-24:00, %s, %s glitches' % (UNIT, len(all_glitches))



    ### PLOT HISTOGRAMS
    fig, axes = plt.subplots(2,3, figsize=(12,7), sharex=True, sharey=True)
    fig.canvas.set_window_title('Glitch Gutenberg-Richter plot')
    fig.suptitle(title, fontsize=12)
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
            hist, bins = np.histogram( abs(all_glitches[:,columns[UNIT][l]].astype('float')), bins=25)
            logbins    = np.logspace(np.log10(bins[0]), np.log10(bins[-1]),len(bins))

        # plot histogram
        y, x, _         = ax.hist( abs(all_glitches[:,columns[UNIT][l]].astype('float')), log=True, range=range_hist, color='k', bins=logbins, ec='w', lw=0.5, label=label_comps[l]+'-comp')
        
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
    #print(u'Fits have been done with negative exponents, so they the `a` are not quite right!')



    ### FINAL
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    if outfile:
        plt.savefig(outfile)
        print(outfile)
    if show:
        plt.show()
    plt.close()
def glitch_envelope_plot(*glitch_files, min_amp=1e-10, max_amp=1e-5, UNIT='DIS', show=True, outfile=None):



    ### HELPER FUNCTION
    def onclick(event):
        sol       = int(event.ydata)
        sol_delta = sol - int(min(sols_range))
        UTC_time  = UTCify(UTCDateTime(mdates.num2date(event.xdata) + datetime.timedelta(days=sol_delta)))
        ptime(UTC_time)



    ### READ GLITCH-FILE
    all_glitches = []

    for glitch_file in glitch_files:
        glitches      = np.loadtxt(glitch_file, dtype='str')
        #glitches      = glitch_exclude(glitches, verbose=False)
        all_glitches += list(glitches)

    all_glitches = np.array(all_glitches)



    ### PREPARE DATA
    all_glitches = all_glitches[((abs(all_glitches[:,columns[UNIT][0]].astype('float'))>=min_amp)  | \
                                 (abs(all_glitches[:,columns[UNIT][1]].astype('float'))>=min_amp)  | \
                                 (abs(all_glitches[:,columns[UNIT][2]].astype('float'))>=min_amp)) & \

                                 (abs(all_glitches[:,columns[UNIT][0]].astype('float'))<=max_amp)  & \
                                 (abs(all_glitches[:,columns[UNIT][1]].astype('float'))<=max_amp)  & \
                                 (abs(all_glitches[:,columns[UNIT][2]].astype('float'))<=max_amp)]



    ### PLOT PREPARATION
    sols         = [solify(UTCDateTime(e[1])).julday for e in all_glitches]
    sols_range   = np.arange(min(sols), max(sols)+1, 1)
    data         = []
    glitches_UVW = np.array( list(all_glitches[:,columns[UNIT][0]].astype('float')) + \
                             list(all_glitches[:,columns[UNIT][1]].astype('float')) + \
                             list(all_glitches[:,columns[UNIT][2]].astype('float')) )
    glitches_UVW = normalise(np.abs(glitches_UVW), scale_to_between=[0,1])

    for k, sol in enumerate(sols_range):
        times             = glitches[:,2]
        indices_sol       = np.array( [i for i in range(len(times)) if solify(UTCDateTime(times[i])).julday==sol] )
        if not indices_sol.any():
            continue
        
        sol_times         = times[indices_sol]
        Matlab_lmst_times = [mdates.date2num(solify(UTCDateTime(time)).datetime-datetime.timedelta(days=k)) for time in sol_times]
        
        y_data_U          = glitches_UVW[indices_sol+0*len(times)] + sol
        y_data_V          = glitches_UVW[indices_sol+1*len(times)] + sol
        y_data_W          = glitches_UVW[indices_sol+2*len(times)] + sol
        
        data.append( [Matlab_lmst_times, y_data_U, y_data_V, y_data_W] )



    ### PLOT
    titles = ['U-component','V-component','W-component']

    fig, axes = plt.subplots(1, 3, figsize=(15,7), sharex=True, sharey=True)
    fig.canvas.set_window_title('Glitch envelope plot')
    fig.suptitle('Envelope plotter for %s glitches (%s)' % (len(glitches), UNIT), fontsize=11, y=0.99)
    fig.subplots_adjust(wspace=0.4, hspace=0.4)
    fig.canvas.mpl_connect('button_press_event', onclick)
    fig.autofmt_xdate()

    for l, ax in enumerate(axes):

        ax.set_title(titles[l], size=10)
        ax.autoscale(enable=True, axis='both', tight=False)
        ax.grid(ls='-.', lw=0.5, zorder=0)
        ax.set_xlabel('LMST')
        ax.xaxis.labelpad = 3
        ax.xaxis.set_major_formatter( mdates.DateFormatter('%H:%M:%S') )
        ax.xaxis.set_major_locator(plt.MaxNLocator(5))
        #ax.xaxis.set_major_locator( mdates.HourLocator(interval=4) )
        ax.set_ylabel('Sols')
        ax.set_ylim([min(sols_range),max(sols_range)+1])

        for line in data:
            x = line[0]
            y = line[l+1]
            ax.plot(x, y, lw=0.75, c='gray')

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    if outfile:
        plt.savefig(outfile)
        print(outfile)
    if show:
        plt.show()
    plt.close()
def glitch_XoverBAZ_plot(*glitch_files, LMST_hour=None, min_amp=1e-10, max_amp=1e-5, comp='Z', UNIT='DIS', show=True, outfile=None):

    """

    """



    ### READ GLITCH-FILE
    all_glitches = []

    for glitch_file in glitch_files:
        glitches      = np.loadtxt(glitch_file, dtype='str')
        #glitches      = glitch_exclude(glitches, verbose=False)
        all_glitches += list(glitches)

    all_glitches = np.array(all_glitches)
    if LMST_hour or LMST_hour==0:
        all_glitches = np.array( [e for e in all_glitches if solify(UTCDateTime(e[1])).hour == LMST_hour] )


    ### PREPARE DATA
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



    ### PREPARE PLOT
    ### SOME VARIABLES DECLARATION
    if LMST_hour != False and LMST_hour!=None:
        title = '%s over BAZ plotter (LMST %02d:00-%02d:00, %s, %s glitches)' % (LMST_hour, LMST_hour+1, comp.upper(), UNIT, len(all_glitches))
    else:
        title = '%s over BAZ plotter (LMST 00:00-24:00, %s, %s glitches)' % (comp.upper(), UNIT, len(all_glitches))

    glitch_starts = np.array( [solify(UTCDateTime(e)).hour for e in all_glitches[:,1]] )
    U             = abs(all_glitches[:,columns[UNIT][0]].astype('float'))
    V             = abs(all_glitches[:,columns[UNIT][1]].astype('float'))
    W             = abs(all_glitches[:,columns[UNIT][2]].astype('float'))
    Z             = abs(all_glitches[:,columns[UNIT][3]].astype('float'))
    N             = abs(all_glitches[:,columns[UNIT][4]].astype('float'))
    E             = abs(all_glitches[:,columns[UNIT][5]].astype('float'))
    BAZ           = abs(all_glitches[:,columns[UNIT][6]].astype('float'))*np.pi/180



    ### FIGURE
    fig = plt.figure(figsize=(8,8))
    fig.canvas.set_window_title('Glitch %s over BAZ plotter'    %  comp.upper())
    fig.suptitle(title, fontsize=11, y=0.99)
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



### _ _ N A M E _ _ = = " _ _ M A I N _ _ "  
if __name__ == "__main__":
    pass