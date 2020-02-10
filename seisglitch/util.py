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
import copy
import time
import glob
import yaml
import fnmatch
import datetime
import itertools
import warnings
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
import numpy as np


#####  matplotlib modules import  #####
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.switch_backend('TKAgg')
import matplotlib.dates as mdates
from matplotlib.patches import Arc, FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
from matplotlib.widgets import TextBox


#####  obspy modules import  #####
from obspy import read, read_inventory, UTCDateTime
from obspy.core.stream import Stream, Trace
from obspy.core.inventory import Inventory
from obspy.signal import rotate
from obspy.clients.fdsn import Client


##### seisglitch modules import #####
from seisglitch.math import normalise, covariance_matrix, rotate_2D, unit_vector


# Extended Obspy Stream class (adding some more functionality)
def read2(file=None, headonly=False):
    # wrapper to make to return Stream2 objects instead of (ObsPy's) Stream object.
    st = read(file, headonly=headonly)
    st = Stream2(st, file=file)
    return st
class Stream2(Stream):

    """ Extended class for ObsPy's stream object. """

    def __init__(self, traces=None, file=None):
        super().__init__(traces=traces)
        self.origin       = str(file)
        self.original     = None
        self.removed      = None
        self.gain_removed = False
        self.unit         = 'RAW'
        self.times        = self._get_times()
        self.inventory    = None
        self.filters      = {'0' : {'type' : False,       'options' : {                                                               },   'string':'0: no filter'},
                             '1' : {'type' : 'highpass',  'options' : {'freq':1.0,                       'corners':3, 'zerophase':True},   'string':'1: >1 Hz'},
                            #'2' : {'type' : 'highpass',  'options' : {'freq':0.1,                       'corners':3, 'zerophase':True},   'string':'2: >0.1 Hz'},
                            #'3' : {'type' : False,       'options' : {                                                               },   'string':'3: >0.01 Hz'},
                             '2' : {'type' : 'bandpass',  'options' : {'freqmin':0.001,  'freqmax':1.0,  'corners':3, 'zerophase':True},   'string':'2: 0.001s < f < 1'},
                             '3' : {'type' : 'bandpass',  'options' : {'freqmin':0.001,  'freqmax':0.1,  'corners':3, 'zerophase':True},   'string':'3: 0.001s < f < 0.1'},
                             '4' : {'type' : 'bandpass',  'options' : {'freqmin':0.01,   'freqmax':0.1,  'corners':3, 'zerophase':True},   'string':'4: 0.01s < f < 0.1'},
                             '5' : {'type' : 'highpass',  'options' : {'freq':0.001,                     'corners':3, 'zerophase':True},   'string':'5: 0.001s < f'},
                            #'6' : {'type' : 'highpass',  'options' : {'freq':1./200,                    'corners':3, 'zerophase':True},   'string': '6: Lowpass 200s'},
                            #'7' : {'type' : 'highpass',  'options' : {'freq':1./300,                    'corners':3, 'zerophase':True},   'string':'7: Lowpass 300s'},
                            #'8' : {'type' : 'highpass',  'options' : {'freq':1./400,                    'corners':3, 'zerophase':True},   'string':'8: Lowpass 400s'},
                            #'9' : {'type' : 'highpass',  'options' : {'freq':1./500,                    'corners':3, 'zerophase':True},   'string':'9: Lowpass 500s'},
                            #'5' : {'type' : 'bandpass',  'options' : {'freqmin':0.1,    'freqmax':2.0,  'corners':3, 'zerophase':True},   'string':'5: RFs (0.1-2 Hz)'},
                             '6' : {'type' : 'bandpass',  'options' : {'freqmin':1/9.0,  'freqmax':1.0,  'corners':3, 'zerophase':True},   'string':'6: BP 1-9 s'},
                             '7' : {'type' : 'bandpass',  'options' : {'freqmin':1/6.0,  'freqmax':1.0,  'corners':3, 'zerophase':True},   'string':'7: BP 1-6 s'},
                             '8' : {'type' : 'bandpass',  'options' : {'freqmin':2.0,    'freqmax':3.0,  'corners':3, 'zerophase':True},   'string':'8: 2-3 Hz'},
                             '9' : {'type' : 'bandpass',  'options' : {'freqmin':2.0,    'freqmax':10,   'corners':3, 'zerophase':True},   'string':'9: BP 2-10 Hz'},
                             }
    def is_unique(self):

        """
        Takes a ObsPy stream object and returns its station code of the format: network.station.location.channel.
        If this is non-unique, i.e., either of them is non-unique, i.e., the stream object contains more than one
        specific component, this script returns `False`. If it is unique, it returns `True`.

        This script can be used to make sure there is only one component within the passed stream.
        This is useful when checking for gaps/overlaps or for any other application
        where it's mandatory for correct processing to have only one component within the stream.
        """

        networks  = []
        stations  = []
        locations = []
        channels  = []

        for tr in self:
            networks.append(tr.stats.network)
            stations.append(tr.stats.station)
            locations.append(tr.stats.location)
            channels.append(tr.stats.channel)
    
        if len(set(networks))==1 and len(set(stations))==1 and len(set(locations))==1 and len(set(channels))==1:
            return True
        else:
            return False
    def trim2(self, *args, samples=None, common=False, **kwargs):

        """
        Use Obspy's trim function and improve slightly.
        Improvement is that you don't have to type
        `UTCDateTime` every time you specify a time.
        """


        ## minimum and maximum times of stream
        mint = min([tr.stats.starttime for tr in self])
        maxt = max([tr.stats.endtime   for tr in self])


        ## Try to convert positional arguments passed to `UTCDateTime` 
        args  = list(args)
        args2 = np.array(args)

        if args2 is not None:
            for l, arg in enumerate(args2[:2]):

                if isinstance(arg, (str,int,float)):
                    try:
                        args2[l] = float(args2[l])
                    except ValueError:
                        pass

                if l==0:
                    kwargs['starttime'] = args2[l]
                    del args[0]
                elif l==1:
                    kwargs['endtime']   = args2[l]
                    del args[0]


        ## Try to convert keyword arguments passed to `UTCDateTime` 
        try:
            kwargs['starttime']
        except KeyError:
            kwargs['starttime'] = mint
        
        try:
            kwargs['endtime']
        except KeyError:
            kwargs['endtime']   = maxt


        ## if start or endtimes or floats, convert to percentage of respective start and ent times
        if isinstance(kwargs['starttime'], float):
            kwargs['starttime'] = max(mint, mint + (maxt-mint)*kwargs['starttime'])
        else:
            kwargs['starttime'] = UTCDateTime(kwargs['starttime'])

        if isinstance(kwargs['endtime'], float):
            kwargs['endtime']   = min(maxt,  maxt - (maxt-mint)*(1-kwargs['endtime']))
        else:
            kwargs['endtime']   = UTCDateTime(kwargs['endtime'])


        ## cut to samples length, if wished:
        if samples:
            for trace in self:
                kwargs['endtime'] = kwargs['starttime'] + (samples-1)*trace.stats.delta
                trace.trim(*args, **kwargs)
        else:
            self.trim(*args, **kwargs)


        ## trim common, if wished
        if common:
            ids = self._get_ids()

            minis, maxis = [], []
            for id in ids:
                mini = min([trace.stats.starttime for trace in self.select(id=id)])
                maxi = max([trace.stats.endtime   for trace in self.select(id=id)])
                minis.append(mini)
                maxis.append(maxi)

            kwargs['starttime'] = max(minis)
            kwargs['endtime']   = min(maxis)
            #print(kwargs['starttime'])
            #print(kwargs['endtime'])
            
            self.trim(*args, **kwargs)


        ## update Stream2 times
        self.times = self._get_times()
    def normalise(self, scale_to_between=[]):

        """
        Normalise each trace in stream.
        scale_to_between is list with 2 elements.
        """

        for tr in self:
            tr.data = normalise(tr.data, scale_to_between=scale_to_between)
    def print_stats(self):

        """
        Print stats for each each within the stream.
        """

        print('\nStats:')
        for trace in self:
            print('  %s' % trace.id)
            for info in trace.stats.keys():
                if info == trace.stats._format.lower() or info == 'processing':
                    continue
                print('%17s : %s' % (info, trace.stats[info]))
            print('%17s : %s' % ('.. and', 'processing_steps & %s-dictionary' % trace.stats._format))
    def print_process(self, truncate=100):
    
        """
        Print Processing steps for all trace objects
        contained in the passed ObsPy stream object.
        """
    
        print('\nProcessing steps:')
        for trace in self:
            print('  %s:' % trace.id)
            try:
                processed_steps = trace.stats.processing
            except AttributeError:
                print('     Nothing processed.')
            else:
                for processed_step in processed_steps:
                    matches = re.findall(r'(\w*\(.*\))\Z',processed_step)                       
                    for match in matches:
                        for i in range( np.ceil(len(match)/truncate) ):
                            if i == 0:
                                print('    %s' % match[i*truncate:(i+1)*truncate])
                            else:
                                print('        %s' % match[i*truncate:(i+1)*truncate])
    def print_filter(self):
        """
        Print filter of stream, if it was fitlered.
        """
        print( self.current_filter_str )
    def gain_correction(self, verbose=False): 
        if self.inventory:

            if verbose:
                print()
                print(u'GAIN CORRECTION APPLIED:')


            try:

                if not self.gain_removed:
                    self.gain_removed = True
                    for trace in self:
                        response   = self.inventory.get_response(trace.id, trace.stats.starttime)
                        gain       = response._get_overall_sensitivity_and_gain()[1]                
                        trace.data = trace.data / gain

                        if verbose:
                            print(u'  %15s : overall sensitivity and gain (division) %s' % (trace.id, gain))
                else:
                    self.gain_removed = False
                    for trace in self:
                        response   = self.inventory.get_response(trace.id, trace.stats.starttime)
                        gain       = response._get_overall_sensitivity_and_gain()[1]                
                        trace.data = trace.data * gain

                        if verbose:
                            print(u'  %15s : overall sensitivity and gain (multiplication) %s' % (trace.id, gain))

            except Exception as err:
                print(u'WARNING:  %s' % err)

        else:
            print()
            print(u'No matching response found for gain correction. Nothing done.')
    def filtering(self, filter_num):
    
        """ 
        Filter ObsPy's stream object has to the given filter specifications. 
        If only freqmax or freqmin is given, highpass or lowpass filter is used, respectively.
        By default the Butterworth filter implementation is used. For a different filter,
        change the code!
    
        See details:
          https://docs.obspy.org/packages/autogen/obspy.core.stream.Stream.filter.html
        """
    
        if isinstance(filter_num, (str, int)):
            filt = self.filters[str(filter_num)]

        if filt['type']:
            self.taper(max_percentage=0.03)
            self.filter(filt['type'], **filt['options'])

        self.current_filter_num = filter_num
        self.current_filter_str = filt['string']
    def _get_filter_str(self):
        # return filter of stream object (use only first trace to check)
        try:
            filt = self.current_filter_str

        except AttributeError:
            self.current_filter_num = 0
            filt = 'no filter'
            try:
                for processed_step in self[0].stats.processing:
                    if 'filter' in processed_step:
                        filt = re.findall(r'filter\(options={(.*)}', processed_step)[0]
                        break
            except AttributeError:
                pass

        return filt
    def _get_components(self):
        components = '*'
        try:
            components = sorted( set( [tr.stats.channel[-1] for tr in self] ))
        except AttributeError:
            pass

        return components
    def _get_channels(self):
        channels = '*'
        try:
            channels = sorted( list( set( [tr.stats.channel for tr in self] )), reverse=False)
        except AttributeError:
            pass

        return channels
    def _get_ids(self):
        ids = '*'
        try:
            ids = sorted( set( [tr.id for tr in self] ), reverse=False)
        except AttributeError:
            pass

        return ids
    def _get_times(self):
        try:
            times = min([tr.stats.starttime for tr in self]), max([tr.stats.endtime for tr in self])
        except ValueError:  # empty stream object
            times = None, None

        return times
    def _get_inventory(self, source):

        """
        Gets latest `source` and reads it as an Obspy inventory object.
        From this inventory, only the part is extracted that matches the start and end times of the stream `self`,
        as well as matches all networks, stations, locations and channels of the stream.

        The inventory is then returned.

        https://docs.obspy.org/packages/obspy.clients.fdsn.html
        """


        inv = Inventory()

        try:

            inv_file  = max( glob.glob( source ), key=os.path.getctime)   # most recent one, in case there multiple
            inventory = read_inventory(inv_file)

            for trace in self:
                network, station, location, channel = trace.id.split('.')
                inv += inventory.select(network   = network, 
                                        station   = station, 
                                        location  = location, 
                                        channel   = channel, 
                                        starttime = self.times[0], 
                                        endtime   = self.times[1])
            
            print()
            print(u'Info: Inventory read from source %s' % source)


        except Exception:      # `source` couldn't be read, because it was e.g. something like 'IRIS' or 'IPGP'
            
            print()
            print(u'Info: Inventory retrieval online from %s' % source)
           
            client = Client(source)

            for trace in self:
                network, station, location, channel = trace.id.split('.')
                inv += client.get_stations(network = network, 
                                         station   = station, 
                                         location  = location, 
                                         channel   = channel,
                                         starttime = self.times[0], 
                                         endtime   = self.times[1],
                                         level     = 'response')                

        return inv
    def _set_inventory(self, source='IRIS', inventory=None):

        """
        
        """


        # small output
        if inventory is not None:
            print()
            print(u'Info: Inventory used that was passed.')
        

        # if not `inventory` is passed, get it from `source` 
        else:
            inventory = self._get_inventory(source=source)


        # assign final inventory to object
        self.inventory = inventory
    def rotate_2D(self, angle, components='NE', new_components='12', clockwise=False):

        """
        """

        if len(components) != 2:
            print('To rotate data in 2-D, you need to give two')
            print('component letters. You gave %s' % len(components))

        self.trim_common()
        comp_1 = self.select2(component=components[0])
        comp_2 = self.select2(component=components[1])
        if comp_1 and comp_2:
            comp_1[0].data, comp_2[0].data = rotate_2D(comp_1[0].data, comp_2[0].data, angle, clockwise=clockwise)
            comp_1[0].stats.channel        = comp_1[0].stats.channel[:-1] + new_components[0]
            comp_2[0].stats.channel        = comp_2[0].stats.channel[:-1] + new_components[1]
            print('Rotated `components` (%s) to `new_components` (%s)' % (components, new_components))
            print('`angle` (%s), `clockwise` (%s)' % (angle, clockwise))
        else:
            print('No rotation performed, as `components` (%s)' % components)
            print('are not contained in stream.')
    def select2(self, *args, **kwargs):
        st_obs_select              = self.select(*args, **kwargs)
        st_obs_select.origin       = self.origin   
        st_obs_select.original     = self.original
        st_obs_select.removed      = self.removed
        st_obs_select.gain_removed = self.gain_removed
        st_obs_select.unit         = self.unit
        st_obs_select.times        = self.times
        st_obs_select.inventory    = self.inventory
        return st_obs_select
    def plot2(self, store_dir=os.getcwd(), store_name='*', verticals=(), time='normal', method='full', xlim=[], ylim=[], store_plot=False, show=True, **kwargs):


        """
        Enhanced plot of stream as 'time', using ObsPy's plot function.
          check: https://docs.obspy.org/packages/autogen/obspy.core.stream.Stream.plot.html
      
          'store_name' is the name of the plot when it shall be saved (default: .png) in 'store_dir'.
      
          'verticals' (a tuple) is needed to plot a vertical line in a specific axis of the seismogram.
          Can be used for example to plot seismic phase arrival times or similar. The times, however,
          must be given in seconds relative to starttime of axis. Example:
          verticals = (
                      (time, label, colour, index),
                      (..), 
                      )
          Note: - time must be passed as relative seconds with respect to axis beginning
                - colour can be chosen from the matplotlib colours
                - index should be integer, but can be 'all' (to be plotted in all axes)
      
          If 'save_and_no_show'=True, no plot is shown but saved
          if 'save_and_no_show'=False, plot is shown and further options will be displayed (e.g. saving)
      
          Following mouse and button events are available:
            - double left click, 's'         : save figure in pwd as `store_name`.png
            - double right click, 'q'  '     : close figure
            - left click                     : give one of two limits for polarization analysis
            - right click                    : delete all current picks done by left click
            - 'p', or enter                  : plot polarization if two picks have been chosen (left click)
            - 'C'                            : abord and exit python
            - '1' (or int number)            : filter seismogram (if wasn't filtered before), as specified in `self.filters` (see __init__)
        """
    
        # helper functions
        def on_close(evt):
            if flags['saved']:
                print(u'Plot saved: %s' % outfile)
            else:
                pass
                #print(u'Plot not saved.')
        def on_click(evt):
            axis        = evt.inaxes
            try:
                stat_id = '.'.join( [c for c in axis.get_children() if isinstance(c, mpl.text.Text) ][0].get_text().split('.')[:3] )
                indexes = [i for i, ax in enumerate(axes) if stat_id in [c for c in ax.get_children() if isinstance(c, mpl.text.Text)][0].get_text() ]
            except AttributeError:
                #if clicking not on axes but on canvas ..
                return
    
            if not evt.dblclick:
                if axis.get_navigate_mode() == 'ZOOM':  
                    return                                                                      ### single click

                try:
                    mouse_clicks[stat_id]
                except:
                    mouse_clicks[stat_id] = 0

                if evt.button==3:                                                               # right-click
                    for i in indexes:
                        lines = axes[i].lines
                        remove = []
                        for l in range(len(lines)):
                            xdata = lines[l].get_xdata()
                            if xdata[0] == xdata[-1]:
                                remove.append(l)
                    
                        remove.sort(reverse=True)
                        for k in remove:
                            lines[k].remove()
                    mouse_clicks[stat_id] = 0


                if evt.button == 1:                                                             # left-click

                    if mouse_clicks[stat_id]%2==0 and mouse_clicks[stat_id]>0:
                        for i in indexes:
                            axes[i].lines[-2].remove()
                        mouse_clicks[stat_id] = 1

                    mouse_clicks[stat_id] += 1  
                    x = [evt.xdata, evt.xdata]
                    for i in indexes:
                        y = []
                        for line in axes[i].lines:

                            xdata = line.get_xdata()
                            if xdata[0] == xdata[-1]:
                                # avoid vertical lines that may be present already
                                continue

                            ymin, ymax = axes[i].get_ylim()
                            y.append(ymin)
                            y.append(ymax)

                        margin = (max(y)-min(y))*1.05
                        ylim   = [ min(y)-margin, max(y)+margin ]
                        axes[i].plot(x, ylim, 'deepskyblue')
                fig.canvas.draw()
    

            else:                                                                           ### double click                
                if evt.button == 1:                                                         # nothing right now
                    pass
                    #flags['saved'] = True
                    #plt.savefig(outfile, dpi=200, bbox_inches='tight')
                    #plt.close(fig)
                elif evt.button == 2:                                                       # middle-click (set ylim for each axis to minimum of all axes and maximum of all axes)
                    ylims = []
                    for axis in axes:
                        ylims.append( axis.get_ylim() )
                    max_y = max( np.array(ylims)[:,0] )
                    min_y = min( np.array(ylims)[:,1] )
                    for axis in axes:
                        axis.set_ylim( [max_y,min_y] )
                    fig.canvas.draw()
                elif evt.button == 3:                                                       # right-click (for each axis individually, set ylim to respective data minimum and maximum + margin)
                    xlim   = axis.get_xlim()
                    for axis in axes:
                        y_mins = []
                        y_maxs = []
                        for line in axis.lines:
                            x        = line.get_xdata()
                            if x[0]==x[-1]:                                                 # exclude vertical lines
                                continue
                            y        = line.get_ydata()
                            i        = np.where( (x >= xlim[0]) & (x <= xlim[1]) )[0]       # all indexes of y_data according to xlim
                            if not i.any():
                                continue                                                    # e.g. data gaps, that is, 2 lines within axes where one does not lie within chose xlims

                            line_min = y[i].min()                                           # get minimum y within all data according to xlim
                            line_max = y[i].max()                                           # get maximum y within all data according to xlim
                            y_mins.append(line_min)
                            y_maxs.append(line_max)
                        
                        y_min = min(y_mins)
                        y_max = max(y_maxs)
                        ylim  = [y_min-0.025*np.abs(y_max-y_min), y_max+0.025*np.abs(y_max-y_min)] # give 2.5% margin that works both for pos. & neg. values 
                        axis.set_ylim( ylim )       
                    fig.canvas.draw()
        def on_key(evt):
            if evt.inaxes:
                axis = evt.inaxes
            else:
                axis = axes[0]

            def instrument_removal(output, pre_filt=None, water_level=0):
                xlim   = axis.get_xlim()
                output = output.upper()

                if self.current_filter_num == 0 and not self.removed:
                    corr          = self.copy()
                    corr.original = self.copy()
                    corr.removed  = self.copy()
                    corr.unit     = output

                elif self.current_filter_num == 0 and self.removed:
                    corr          = self.removed.copy()
                    corr.original = self.removed.copy()
                    corr.removed  = None
                    corr.unit     = 'RAW'
                    corr.unit    +='wt%s' % water_level
                    corr.plot2(xlim=xlim)
                    return

                elif self.current_filter_num != 0 and not self.removed:
                    corr          = self.original.copy()
                    corr.original = self.original.copy()
                    corr.removed  = self.original.copy()
                    corr.unit     = output

                elif self.current_filter_num != 0 and self.removed:
                    corr          = self.removed.copy()
                    corr.original = self.removed.copy()
                    corr.removed  = None
                    corr.unit     = 'RAW'
                    corr.unit    +='wt%s' % water_level
                    corr.filtering(self.current_filter_num)
                    corr.plot2(xlim=xlim)
                    return


                if not mouse_clicks.keys():
                    print()
                    print('Please select at least one station location!')

                else:
                    if not self.inventory:
                        self._set_inventory()
                    

                    components = ''.join( self._get_components() )

                    for stat_id in mouse_clicks.keys():
                        corr_part  = corr.select(id=stat_id+'*')
                        corr_part.merge()

                        corr_part2 = corr.original.select(id=stat_id+'*')
                        corr_part2.merge()
                
                        if re.search(r'[U|V|W]', components):
                            corr_part.remove_response(inventory=self.inventory, output=output, pre_filt=pre_filt, water_level=water_level)
                            corr_part2.remove_response(inventory=self.inventory, output=output, pre_filt=pre_filt, water_level=water_level)

                        elif re.search(r'[Z|N|E]', components):
                            corr_part  = rotate2VBBUVW(corr_part, inventory=self.inventory)
                            corr_part.remove_response(inventory=self.inventory, output=output, pre_filt=pre_filt, water_level=water_level)
                            corr_part.rotate('->ZNE', inventory=self.inventory, components='UVW')

                            corr_part2 = rotate2VBBUVW(corr_part2, inventory=self.inventory)
                            corr_part2.remove_response(inventory=self.inventory, output=output, pre_filt=pre_filt, water_level=water_level)
                            corr_part2.rotate('->ZNE', inventory=self.inventory, components='UVW')

                    corr.filtering(self.current_filter_num)
                    #corr.unit    +='wt%s' % water_level
                    corr.plot2(xlim=xlim)


            if evt.key.lower()=='a':                            # Acceleration as target unit after instrument response removal 
                pre_filt    = (0.001, 0.002, 50, 60)
                instrument_removal('ACC')


            elif evt.key.lower()=='b':                          # Before xlim were changed, meaning go back to full available data view 
                
                # set new xlim of all available data
                for axis in axes:
                    x = []
                    for line in axis.lines:
                        x += list( line.get_xdata() )
                    xlim = [min(x), max(x)]

                axis.set_xlim(xlim)

                # set best ylim for new xlim
                for axis in axes:
                    y_mins = []
                    y_maxs = []
                    for line in axis.lines:
                        x        = line.get_xdata()
                        if x[0]==x[-1]:                                                 # exclude vertical lines
                            continue
                        y        = line.get_ydata()
                        i        = np.where( (x >= xlim[0]) & (x <= xlim[1]) )[0]       # all indexes of y_data according to xlim
                        if not i.any():
                            continue                                                    # e.g. data gaps, that is, 2 lines within axes where one does not lie within chose xlims

                        line_min = y[i].min()                                           # get minimum y within all data according to xlim
                        line_max = y[i].max()                                           # get maximum y within all data according to xlim
                        y_mins.append(line_min)
                        y_maxs.append(line_max)
                    
                    y_min = min(y_mins)
                    y_max = max(y_maxs)
                    ylim  = [y_min-0.025*np.abs(y_max-y_min), y_max+0.025*np.abs(y_max-y_min)] # give 2.5% margin that works both for pos. & neg. values 
                    axis.set_ylim( ylim )   
                
                plt.draw()


            elif evt.key.lower()=='c':                          # Cease all 
                print('STOPPED ALL')
                plt.close('all')
                sys.exit()


            elif evt.key.lower()=='d':                          # Displacement as target uniy after instrument response removal 
                pre_filt    = (0.001, 0.002, 50, 60)
                instrument_removal('DISP')


            elif evt.key.lower()=='e':                          # Velocity as target uniy after instrument response removal 
                pre_filt    = (0.001, 0.002, 50, 60)
                water_level = 0
                instrument_removal('VEL', water_level=water_level)
           

            elif evt.key.lower()=='j':                          # Acceleration as target uniy after instrument response removal 
                pre_filt    = (0.001, 0.002, 50, 60)
                water_level = 0
                instrument_removal('ACC', water_level=water_level)
            

            elif evt.key.lower()=='g':                          # Gain removal / add 
                self.gain_correction(verbose=True)
                self.plot2(xlim=axis.get_xlim())


            elif evt.key.lower()=='h':                          # Help summary 
                # TO BE IMRPOVED
                pass


            elif evt.key.lower()=='m':                          # Narrow view 10% (zooming in) 
                xlim = axis.get_xlim()
                xlim_range = (xlim[1]-xlim[0])
                xlim = [xlim[0]+0.1*xlim_range, xlim[1]-0.1*xlim_range]
                for axis in axes:
                    axis.set_xlim( xlim )
                plt.draw()  


            elif evt.key.lower()=='n':                          # Move out 10% (zooming out) 
                xlim = axis.get_xlim()
                xlim_range = (xlim[1]-xlim[0])
                xlim = [xlim[0]-0.1*xlim_range, xlim[1]+0.1*xlim_range]
                for axis in axes:
                    axis.set_xlim( xlim )
                plt.draw()  


            elif evt.key=='p':                                  # Polarisation plot 
                for l, stat_id in enumerate( mouse_clicks.keys() ):
                    if mouse_clicks[stat_id] == 1:
                        continue
                    elif mouse_clicks[stat_id] != 2:
                        print('Choose at least 2 time limits (left-click) per station to perform action.')
                        continue

                    if l==len(mouse_clicks.keys())-1:
                        show=True
                    else:
                        show=False                      
                    indexes       = [i for i, ax in enumerate(axes) if stat_id in [c for c in ax.get_children() if isinstance(c, mpl.text.Text)][0].get_text() ]
                    ax_1st        = axes[indexes[0]]
                    xlim          = sorted( [ax_1st.lines[-1].get_xdata()[0], ax_1st.lines[-2].get_xdata()[0]] )
                    if time == 'normal':
                        xlim   = [ UTCDateTime(mdates.num2date(x)) for x in xlim ]
                    elif time == 'relative':
                        xlim =  [ min_starttime+x for x in xlim ]
                    
                    polar = self.copy()
                    polar = polar.select(id=stat_id+'*')
                    polar.trim(starttime=xlim[0], endtime=xlim[1])

                    title = '%s, %s, %s, %s-%s'% (stat_id, self.unit, self._get_filter_str(), xlim[0].strftime('%H:%M:%S'), xlim[1].strftime('%H:%M:%S'))
                    measurement = ppol(stream=polar)
                    measurement.display()
                    measurement.plot(title=title)


            elif evt.key=='P':                                  # Polarisation plot as batch job 
                for l, stat_id in enumerate( mouse_clicks.keys() ):
                    if mouse_clicks[stat_id] == 1:
                        continue
                    elif mouse_clicks[stat_id] != 2:
                        print('Choose at least 2 time limits (left-click) per station to perform action.')
                        continue
                    print()
                    print('Running Ppol batch job ..')


                    indexes       = [i for i, ax in enumerate(axes) if stat_id in [c for c in ax.get_children() if isinstance(c, mpl.text.Text)][0].get_text() ]
                    ax_1st        = axes[indexes[0]]
                    xlim          = sorted( [ax_1st.lines[-1].get_xdata()[0], ax_1st.lines[-2].get_xdata()[0]] )
                    if time == 'normal':
                        xlim   = [ UTCDateTime(mdates.num2date(x)) for x in xlim ]
                    elif time == 'relative':
                        xlim =  [ min_starttime+x for x in xlim ]
                    
                    # create windows for batch job
                    xlims    = []
                    span     = xlim[1]-xlim[0]
                    per_side = 2
                    fraction = 0.15
                    for l in np.arange(-per_side, per_side+1):
                        for i in np.arange(-per_side, per_side+1):
                            window = [ xlim[0]+l*span*fraction, xlim[1]+i*span*fraction]
                            xlims.append( window )

                    # iterate over windows and filters to calculate polarisations
                    results = []
                    unit    = self.unit
                    for xlim in xlims:
                        for fil_num in self.filters.keys():
                            polar = self.original or self
                            polar = polar.copy()
                            polar = polar.select(id=stat_id+'*')
                            polar.filtering(fil_num)
                            polar.trim(starttime=xlim[0], endtime=xlim[1])
                            measurement = ppol(stream=polar)
                            results.append( measurement.results )

                    # retrieve needed variables for best polarisation
                    results      = np.array(results)
                    POLs_2D      = results[:,12]
                    ind_best     = np.argmax(POLs_2D)     # highest POL_2D

                    fil_num_best =       int( ind_best %  len(self.filters)) 
                    xlim_best    = xlims[int( ind_best // len(xlims))       ]

                    # plot best polarization
                    polar_best   = self.original or self
                    polar_best   = polar_best.select(id=stat_id+'*')
                    polar_best.filtering(fil_num_best)
                    polar_best.trim(starttime=xlim_best[0], endtime=xlim_best[1])

                    measurement_best = ppol(stream=polar_best)
                    fil_str      = polar_best._get_filter_str()
                    window       = '%s-%s' % (xlim_best[0].strftime('%H:%M:%S'), xlim_best[1].strftime('%H:%M:%S'))
                    title        = '%s, %s, %s, %s' % (stat_id, unit, fil_str, window)

                    # plot histogram
                    hist_title = 'Histogram P-pol batch job (%s)' % self.unit
                    if title in plt.get_figlabels():
                        fig = plt.figure(hist_title)
                        fig.clf()
                    fig = plt.figure(num=hist_title)
                    fig.suptitle('Histogram of %s P-wave polarisation calculations' % (len(POLs_2D)), fontsize=9)
                    ax = fig.add_subplot(111)
                    ax.grid(ls='-.', lw=0.5, zorder=0)
                    ax.set(xlabel='Horizontal polarization (POL_2D)', ylabel='Number of occurences')
                    ax.hist(POLs_2D, bins=100, range=(0,1), rwidth=1, ec='k', lw=0.5, zorder=3)

                    # print and plot ppol_results
                    print('')
                    print('B E S T')
                    print('-' * 8)
                    measurement_best.display(title)
                    print('-' * 8)
                    measurement.plot(title=title)


            elif evt.key.lower()=='q':                          # Quit current figures
                # close
                print('Closed current figure window')
                plt.gcf()
                plt.close()


            elif evt.key.lower()=='r':                          # Rotate data 
                

                ## inventory & components
                if not self.inventory:
                    self._set_inventory()

                components = ''.join( self._get_components() )

                print(components)

                ## rotate data to ZNE
                if re.search(r'[U|V|W]', components):

                    # correct gain
                    if not self.gain_removed and self.unit=='RAW':
                        self.gain_correction(verbose=True)

                        if self.original:
                            self.original.gain_correction(verbose=False)
                            self.original.gain_removed = True

                        if self.removed:
                            self.removed.gain_correction(verbose=False)
                            self.removed.gain_removed = True

                    # rotate
                    self.rotate('->ZNE', inventory=self.inventory, components='UVW')
                    
                    if self.original:
                        self.original.rotate('->ZNE', inventory=self.inventory, components='UVW')

                    if self.removed:
                        self.removed.rotate('->ZNE', inventory=self.inventory, components='UVW')

                
                ## rotate data to UVW
                elif re.search(r'[Z|N|E]', components):

                    rotate2VBBUVW(self, inventory=self.inventory)

                    if self.original:
                        print(self.original[0].data)
                        rotate2VBBUVW(self.original, inventory=self.inventory)
                    
                    if self.removed:
                        rotate2VBBUVW(self.removed, inventory=self.inventory)

                self.plot2(xlim=axis.get_xlim())


            elif evt.key.lower()=='t':                          # Trim file and save present data ("as-is") in current x-bounds (vertical blue lines)
                for stat_id in mouse_clicks.keys():

                    indexes       = [i for i, ax in enumerate(axes) if stat_id in [c for c in ax.get_children() if isinstance(c, mpl.text.Text)][0].get_text() ]
                    ax_1st        = axes[indexes[0]]
                    
                    if mouse_clicks[stat_id]==0:
                        xlim = ax_1st.get_xlim()

                    elif mouse_clicks[stat_id]==2:
                        xlim = sorted( [ax_1st.lines[-1].get_xdata()[0], ax_1st.lines[-2].get_xdata()[0]] )

                    else:
                        print('Choose None or 2 time limits (left-click) per station to trim data & save data.')
                        continue

                    xlim   = [ UTCDateTime(mdates.num2date(x)) for x in xlim ]

                    def submit(filename):
                        trim = self.copy()
                        trim = trim.select(id=stat_id+'*')
                        trim.trim2(starttime=xlim[0], endtime=xlim[1])
                        trim.write(filename)
                        plt.close()
                    
                    fig      = plt.figure(figsize=(17,0.7), num='Save trimmed file')
                    ax       = fig.add_subplot(111)
                    
                    propose  = os.path.join( os.path.dirname(self.origin), 'AsIs_'+os.path.basename(self.origin) )
                    text_box = TextBox(ax, 'Filename:', initial=propose)
                    text_box.on_submit(submit)
                    plt.show()


            elif evt.key.lower()=='u':                          # Update plot as to chosen xlim boundaries 
                for stat_id in mouse_clicks.keys():
                    if mouse_clicks[stat_id] == 1:
                        continue
                    elif mouse_clicks[stat_id] == 2 :
                        indexes       = [i for i, ax in enumerate(axes) if stat_id in [c for c in ax.get_children() if isinstance(c, mpl.text.Text)][0].get_text() ]
                        ax_1st        = axes[indexes[0]]
                        xlim          = sorted( [ax_1st.lines[-1].get_xdata()[0], ax_1st.lines[-2].get_xdata()[0]] )            
                    else:
                        xlim          = axes[0].get_xlim()

                    for axis in axes:
                        y_mins = []
                        y_maxs = []
                        for line in axis.lines:
                            x        = line.get_xdata()
                            if x[0]==x[-1]:                                                 # exclude vertical lines
                                continue
                            y        = line.get_ydata()
                            i        = np.where( (x >= xlim[0]) & (x <= xlim[1]) )[0]       # all indexes of y_data according to xlim
                            if not i.any():
                                continue                                                    # e.g. data gaps, that is, 2 lines within axes where one does not lie within chose xlims
                            
                            line_min = y[i].min()                                           # get minimum y within all data according to xlim
                            line_max = y[i].max()                                           # get maximum y within all data according to xlim
                            y_mins.append(line_min)
                            y_maxs.append(line_max)
            
                        y_min = min(y_mins)
                        y_max = max(y_maxs)
                        ylim  = [y_min-0.025*np.abs(y_max-y_min), y_max+0.025*np.abs(y_max-y_min)] # give 2.5% margin that works both for pos. & neg. values 
                        axis.set_ylim( ylim )       
                        axis.set_xlim( xlim )

                    plt.draw()      


            elif evt.key.lower()=='v':                          # Velocity as target uniy after instrument response removal 
                pre_filt = (0.001, 0.002, 50, 60)
                instrument_removal('VEL')


            elif evt.key.lower()=='w':                          # Window closing (all) 
                # close
                print('Close all figure windows')
                plt.close('all')
            

            elif evt.key.lower()=='z':                          # Demean current window viewed
                xlim = axis.get_xlim()

                for axis in axes:
                    y_mins = []
                    y_maxs = []
                    for line in axis.lines:
                        
                        x_data = line.get_xdata()
                        if not x_data[0] == x_data[-1]: #exclude vertical lines)

                            # set ylim to best ylim of current viewing window
                            i = np.where( (x_data >= xlim[0]) & (x_data <= xlim[1]) )[0]
                            if not i.any():
                                continue

                            # demean window
                            y_data = line.get_ydata()
                            y_data = y_data - y_data[i].mean()
                            line.set_ydata( y_data ) 

                            line_min = y_data[i].min()
                            line_max = y_data[i].max()
                            y_mins.append(line_min)
                            y_maxs.append(line_max)                         
                    y_min = min(y_mins)
                    y_max = max(y_maxs)
                    ylim  = [y_min-0.025*np.abs(y_max-y_min), y_max+0.025*np.abs(y_max-y_min)] # give 2.5% margin that works both for pos. & neg. values 
                    axis.set_ylim( ylim )

                print(u'Demeaned current window.')  
                plt.draw()


            elif evt.key.lower()==',':                          # Move 10% left 
                xlim = axis.get_xlim()
                xlim = xlim - 0.1*(xlim[1]-xlim[0])
                for axis in axes:
                    axis.set_xlim( xlim )
                plt.draw()  


            elif evt.key.lower()=='.':                          # Move 10% right 
                xlim = axis.get_xlim()
                xlim = xlim + 0.1*(xlim[1]-xlim[0])
                for axis in axes:
                    axis.set_xlim( xlim )
                plt.draw()  


            elif evt.key.lower()=='up':                         # Reduce ylim 10% 
                for axis in axes:
                    ylim = axis.get_ylim()
                    ylim = [ylim[0]+0.1*(ylim[1]-ylim[0]), ylim[1]-0.1*(ylim[1]-ylim[0])]
                    axis.set_ylim( ylim )
                plt.draw()  


            elif evt.key.lower()=='down':                       # Increase ylim 10% 
                for axis in axes:
                    ylim = axis.get_ylim()
                    ylim = [ylim[0]-0.1*(ylim[1]-ylim[0]), ylim[1]+0.1*(ylim[1]-ylim[0])]
                    axis.set_ylim( ylim )
                plt.draw()  


            elif evt.key.isdigit():                             # filter data, depending on choice 
                xlim = axis.get_xlim()

                if self._get_filter_str() == 'no filter':
                    filt              = self.copy()
                    filt.original     = self.copy()

                else:
                    filt          = self.original.copy()
                    filt.original = self.original.copy()
                
                filt.removed = self.removed
                filt.filtering(evt.key)
                filt.plot2(xlim=xlim)    
    
        # variables
        mouse_clicks  = {}
        flags         = {'saved' : False}
        colours       = 'bgcm'
        outfile       = os.path.join(store_dir, store_name+'.png')  # change format if desired (used if plot shall be saved)
        title         = '%s, %s' % (self.unit, self._get_filter_str())
        min_starttime = min( np.array( [tr.stats.starttime for tr in self] ))
    
    
        # plotting seismogram
        figs = [[manager.num, manager.canvas.figure.canvas.get_window_title()] for manager in mpl._pylab_helpers.Gcf.get_all_fig_managers()]
        for fignum, figtit in figs:
            if title == figtit:
                fig = mpl.pyplot.figure(fignum)
                plt.close()

        fig = self.plot(type='normal', show=False, handle=True, equal_scale=False, method='full', **kwargs)
        fig.canvas.draw()
        fig.canvas.set_window_title(title)
        fig.canvas.mpl_connect('button_press_event', on_click)
        fig.canvas.mpl_connect('key_press_event', on_key)
        fig.canvas.mpl_connect('close_event', on_close)
        fig.suptitle('Min starttime: %s (julday %d)' % (min_starttime.strftime('%Y-%m-%dT%H:%M:%S.%f').rstrip("0").rstrip("."), min_starttime.julday), fontsize=9)
        
        axes = fig.axes
    

        # loop over axes to add certain stuff
        for i in range(len(axes)):          
            if verticals:
                x_startpoint = axes[i].lines[0].get_xdata()[0]
                for time, label, colour, index in verticals:
                    time = float(time)
                    if type == 'normal':
                        v_line_y     = UTCDateTime( mdates.num2date(x_startpoint) )
                        v_line_y     = (v_line_y+time).datetime
                    elif type == 'relative':
                        v_line_y     = x_startpoint+time
                    if index == 'all' or float(index)==i:
                        axes[i].plot([v_line_y,v_line_y], axes[i].get_ylim(), color=colour, lw=1.0)
    
            if i == 0:
                axes[i].text(0.99, 0.92, self._get_filter_str(), horizontalalignment='right', color='red', transform=axes[i].transAxes, fontsize=8, bbox=dict(boxstyle='round,pad=0.2', fc='white', ec="red", lw=0.5))
    
            if i+1 == len(self):
                axes[i].text(0.99, 0.93, "Close all windows: 'w'.", horizontalalignment='right', color='red', transform=axes[i].transAxes, fontsize=7)
                if verticals:
                    verts = np.unique( np.array(verticals)[:,:-1], axis=0 )
                    verts = verts[verts[:,0].astype('float').argsort()]
                    for j, verts_uni in enumerate(verts):
                        axes[i].text(0.99, 0.95-j*0.06, "%s" % verts_uni[1], horizontalalignment='right', color=verts_uni[2], transform=axes[i].transAxes, fontsize=6)


        # set lims
        if list( xlim ):
            for axis in axes:
                axis.set_xlim(xlim)
        else:
            xlim = axes[0].get_xlim()
        if list( ylim ):
            for k, axis in enumerate(axes):
                if len(ylim)==len(axes) and not [y for y in ylim if not isinstance(y, (list, tuple, np.ndarray))]:
                    axis.set_ylim( ylim[k] )
                else:
                    axis.set_ylim(ylim)
        else:
            for axis in axes:
                y_mins = []
                y_maxs = []
                for line in axis.lines:
                    x        = line.get_xdata()
                    if x[0]==x[-1]:                                                 # exclude vertical lines
                        continue
                    y        = line.get_ydata()
                    i        = np.where( (x >= xlim[0]) & (x <= xlim[1]) )[0]       # all indexes of y_data according to xlim
                    if not i.any():
                        continue                                                    # e.g. data gaps, that is, 2 lines within axes where one does not lie within chose xlims

                    line_min = y[i].min()                                           # get minimum y within all data according to xlim
                    line_max = y[i].max()                                           # get maximum y within all data according to xlim
                    y_mins.append(line_min)
                    y_maxs.append(line_max)
                
                y_min = min(y_mins)
                y_max = max(y_maxs)
                ylim  = [y_min-0.025*np.abs(y_max-y_min), y_max+0.025*np.abs(y_max-y_min)] # give 2.5% margin that works both for pos. & neg. values 
                axis.set_ylim( ylim )   


        # show / save figure
        if store_plot:
            flags['saved'] = True
            plt.savefig(outfile, dpi=200, bbox_inches='tight')

        if show:
            plt.show()
            plt.close()


# InSight related (e.g. time conversions, rotation)
class marstime():

    def __init__(self, UTC_time=None, LMST_time=None, sec_per_day_mars=88775.2440, sec_per_day_earth=86400, sol0=UTCDateTime(2018, 11, 26, 5, 10, 50.33508)):

        self.sec_per_day_mars  = sec_per_day_mars
        self.sec_per_day_earth = sec_per_day_earth
        self.sol0              = sol0
        
        self.UTC_time          = UTC_time
        self.LMST_time         = LMST_time
        self.UTC               = None
        self.LMST              = None
        self.sol               = None

        if UTC_time is None and LMST_time is None:
            self.UTC_time = UTCDateTime( datetime.datetime.utcnow() )       # now

        if self.UTC_time is not None:
            self.LMSTify()

        if self.LMST_time is not None:
            self.UTCify()
    def __str__(self):

        """
        Display current UTC time in most common ways.
        """

        output = []

        output.append(u'TIME GIVEN\n')
        output.append(u'----------\n')
        output.append(u'MatLab UTC:     %s\n' % mdates.date2num(self.UTC_time.datetime))
        output.append(u'DateTime UTC:   %s\n' % self.UTC_time.datetime.__repr__())
        output.append(u'UTCDateTime:    %s\n' % self.UTC_time.__repr__())
        output.append(u'LocalMeanSolar: %s'   % self.LMST)

        string = ''.join( output )
        return string
    def UTCify(self, LMST_time=None):

        """
        :copyright:
            Simon Stähler (mail@simonstaehler.com), 2018
            Martin van Driel (Martin@vanDriel.de), 2018
        :license:
            None
        """

        if LMST_time is not None:
            LMST_time  = UTCDateTime(LMST_time)
        elif self.LMST_time is not None:
            LMST_time  = self.LMST_time
        else:
            sys.exit(u'WARNING: No time specified.')
        
        MIT            = float(LMST_time) / self.sec_per_day_earth + 1
        UTC_time       = UTCDateTime(MIT  * self.sec_per_day_mars  + float(self.sol0))
        
        self.LMST_time = UTCDateTime(LMST_time)
        self.LMST      = self._LMST_string()
        self.sol       = int(self.LMST.split('M')[0])

        self.UTC_time  = UTC_time
        self.UTC       = self._UTC_string()
        self.julday    = self.UTC_time.julday

        return UTC_time
    def LMSTify(self, UTC_time=None):

        """
        :copyright:
            Simon Stähler (mail@simonstaehler.com), 2018
            Martin van Driel (Martin@vanDriel.de), 2018
        :license:
            None
        """

        if UTC_time is not None:
            UTC_time  = UTCDateTime(UTC_time)
        elif self.UTC_time is not None:
            UTC_time  = self.UTC_time
        else:
            sys.exit(u'WARNING: No time specified.')


        UTC_time       = UTCDateTime(UTC_time)
        
        MIT            = (UTC_time - self.sol0) / self.sec_per_day_mars
        LMST_time      = UTCDateTime((MIT - 1)  * self.sec_per_day_earth)

        self.UTC_time  = UTC_time
        self.UTC       = self._UTC_string()
        self.julday    = self.UTC_time.julday

        self.LMST_time = LMST_time
        self.LMST      = self._LMST_string()
        self.sol       = int(self.LMST.split('M')[0])

        return LMST_time
    def _UTC_string(self, strftime=('%Y-%m-%dT%H:%M:%S')):

        if self.UTC_time is None:
            print(u'WARNING: Need to pass time first.')
            sys.exit()

        string = self.UTC_time.strftime(strftime)

        return string
    def _LMST_string(self):

        sol    = (self.LMST_time.datetime - UTCDateTime('1969-12-31T00:00:00.000000Z').datetime).days
        string = '%03dM%s' % (sol, self.LMST_time.strftime('%H:%M:%S'))

        return string
def mars_list(sols_range=[], hms='120000', is_UTC=False):

    """
    LIST TIME
    
    Small script to display a range of time. Useful e.g.
    when having to check data at a specific UTC or LMST
    each day / sol.

    `sol_range` can be a list or single numer, e.g.:
       - sols_range=[17,193] (displays these sols)
       - sols_range=17       (displays last 17 sols)
    If no `sols_range` specified, the last 10 sols are displayed by default.

    `is_UTC`=True means `hms` is given as UTC and you would like to
    see the corresponding LMST with respect to `sols_range`.
    `is_UTC`=False means `hms` is given as LMST and you would like to
    see the corresponding UTC with respect to `sols_range`.
    """


    ## VARIABLES
    hms = '%06d' % ( int(hms)*10**(6-len(str(hms))) )

    if isinstance(sols_range, (float, int)):
        sols_range = [sols_range]

    elif not sols_range:
        now        = time.time()
        time_mars  = marstime( UTC_time=now )
        sols_range = [time_mars.sol-10, time_mars.sol]

    else:
        pass


    flow = np.sign(sols_range[-1]-sols_range[0])
    if flow == 0:
        flow = 1
    sols_range = np.arange(sols_range[0], sols_range[-1]+1, flow)


    ## PRINTS
    print('UTC                    LMST')
    print('---                    ----')    
    for sol in sols_range:

        if is_UTC:
            
            time_mars_sol = marstime( LMST_time=UTCDateTime('1970-01-01T00:00:00.000000Z')+datetime.timedelta(days=int(sol)) )
            time_str_ymd = time_mars_sol.UTC_time.strftime('%Y-%m-%d')
            time_str_HMS = '%s:%s:%s' % (hms[0:2], hms[2:4], hms[4:6])
            time_mars    = marstime( time_str_ymd + 'T' + time_str_HMS )

        else:

            LMST_time    = UTCDateTime('1970-01-01T%s:%s:%s.000000Z' % (hms[:2], hms[2:4], hms[4:])) + datetime.timedelta(days=int(sol)-1)
            time_mars    = marstime( LMST_time=LMST_time )

        print('%s    %s' % (time_mars.UTC_time.strftime('%Y-%m-%dT%H:%M:%S'), time_mars.LMST))
def rotate2VBBUVW(stream, inventory, is_UVW=False, plot=False):

    """
    Returns data (no copy of data)!
    No gain correction applied, must be done beforehand!

    https://docs.obspy.org/packages/autogen/obspy.signal.rotate.rotate2zne.html
    """
 

    # VALUES EXTRACTED from INVENTORY FILE
    VBBU_azimuth = 135.11
    VBBV_azimuth = 15.04
    VBBW_azimuth = 254.96
    
    VBBU_dip     = -29.28   # defined as down from horizontal!
    VBBV_dip     = -29.33
    VBBW_dip     = -29.61


    # case data (SP or VBB) are in UVW, first go to ZNE using inventory file
    if is_UVW:
        stream_new = stream.select(component='[UVW]')
        stream_new.sort()
        stream_new.rotate('->ZNE', inventory=inventory, components=('UVW'))     # in case SP data are passed (or similar)

    else:
        stream_new = stream.select(component='[ZNE]')
        stream_new.sort(reverse=True)


    # rotate from ZNE data to VBB UVW
    stream_new[0].data, stream_new[1].data, stream_new[2].data = rotate.rotate2zne( data_1     = stream_new[0].data,
                                                                                    azimuth_1  = VBBU_azimuth, 
                                                                                    dip_1      = VBBU_dip, 

                                                                                    data_2     = stream_new[1].data, 
                                                                                    azimuth_2  = VBBV_azimuth, 
                                                                                    dip_2      = VBBV_dip, 

                                                                                    data_3     = stream_new[2].data, 
                                                                                    azimuth_3  = VBBW_azimuth, 
                                                                                    dip_3      = VBBW_dip, 

                                                                                    inverse    = True)

    # because of sorting, can just index new component
    stream_new[0].stats.channel = stream_new[0].stats.channel[:-1] + 'U'
    stream_new[1].stats.channel = stream_new[1].stats.channel[:-1] + 'V'
    stream_new[2].stats.channel = stream_new[2].stats.channel[:-1] + 'W'


    # plot, if desired
    if plot:
        if is_UVW:
            data_labels=('U old', 'U new', 'V old', 'V new', 'W old', 'W new')
        else:
            data_labels=('Z old', 'U new', 'N old', 'V new', 'E old', 'W new')
        quick_plot(stream[0], stream_new[0], stream[1], stream_new[1], stream[2], stream_new[2], data_labels=data_labels)
    
    return stream_new


# Other
def pierce_stream(file, format_out='mseed'): 

    """
    ObsPy's merge function is first applied to stream.

    See options:
    https://docs.obspy.org/packages/autogen/obspy.core.trace.Trace.__add__.html#obspy.core.trace.Trace.__add__

    After that:
    Pierce stream into X streams, in which each stream contains all unique channels of the input stream,
    all with same start and end times. For example useful, if you need 3-D component, equally trimmed data
    that are contained within one gappy stream ..

    NOTE: all unique traces in stream are considered. so if there are 5 different traces, the
    resulting streams will have the start and end times where all 5 traces are present, nothing more.
    So make sure stream contains only those traces that you'd like.
    """



    ### handle overlaps first (if there any gaps, data will be masked array), that is why then split masked array contiguous unmasked traces to apply piercing
    stream = read2(file)
    stream.merge(method=1, fill_value=None, interpolation_samples=0)
    stream = stream.split()



    ### piercing traces 
    ids           = list(set( [tr.id for tr in stream] ))
    array         = np.array([[len(stream.select(id=id)),id] for id in ids])
    id_max_traces = array[np.argmax(array[:,0].astype('float')),1]

    new_streams   = []
    for l, trace in enumerate( stream.select(id=id_max_traces) ): 

        traces = []
        for tr in stream:
            if tr.id != trace.id:
                if tr.stats.endtime>=trace.stats.starttime and tr.stats.starttime<=trace.stats.endtime:
                    traces.append(tr)

        traces     = [trace] + traces
        tmp_stream = Stream2(traces=traces).copy()
        tmp_ids    = list( set( [tr.id for tr in tmp_stream] ))


        # not enough to pierce, moving on to next one
        if len(tmp_ids)<len(ids):
            continue
        

        # containing multiple traces of one ID
        if len(tmp_stream) != len(ids):

            for possibility in itertools.combinations(tmp_stream, len(ids)):

                possibility_ids = len( set( [trace.id for trace in possibility] ))

                if possibility_ids != len(ids):
                    continue

                else:
                    new_stream = Stream2(traces=possibility)
                    new_stream = new_stream.copy()
                    try:
                        new_stream.trim2(common=True)
                    except ValueError:
                        continue
                    new_streams.append(new_stream)
        
        else:
            tmp_stream.trim2(common=True)
            new_streams.append(tmp_stream)

    new_streams = sorted(new_streams, key=lambda x: x[0].stats.starttime)



    ### WRITING FILES OUT
    for stream in new_streams:
        print( str(stream).replace('\n','\n  '))

    print()
    print(u'STREAMS PIERCED AND WRITTEN:')
    for l, stream in enumerate(new_streams):
        start_str = stream[0].stats.starttime.strftime('%Y%m%dT%H%M%S')
        end_str   = stream[0].stats.endtime.strftime('%Y%m%dT%H%M%S')
        ids       = stream._get_ids()
        comps     = ''.join( stream._get_components() )

        filepath  = os.path.dirname( file )
        filename  = '%s_%s-%s.%s' % (ids[0][:-1] + comps, start_str, end_str, format_out)
        fileout   = os.path.join(filepath, filename)

        stream.write(fileout, format=format_out)
        print(fileout)
def merge_glitch_detector_files(outfile, path='', pattern='*', starttime_sort=True):

    """
    Merging together glitch extration files produced by `glitch_detector`.
    This helps combining into one file representing a certain period, e.g.
    """


    now                   = time.time()
    glitch_detector_files = get_files(path, pattern=pattern)


    ### SANITY CHECK
    if not glitch_detector_files:
        print(u'ERROR: No glitch detector files for merging found.')
        sys.exit()



    ### OUTPUT
    print(u'Merging the following files:')
    for file in glitch_detector_files:
        print(file)



    ### RETRIEVE HEADER AND ALL DATA
    data    = []
    times   = []
    with open(outfile, 'w') as fp: 
        counter = 0

        for file in glitch_detector_files: 

            with open(file, 'r') as fp2:
                lines    = fp2.readlines() 
                header   = lines[0:3]
                data    += [line for line in lines[3:] if line.strip() and not line.startswith('#')]
                
                start    = [marstime(line.split()[4]) for line in lines[3:] if line.strip() and line.startswith('#') and 'UTC start:' in line][0]
                end      = [marstime(line.split()[4]) for line in lines[3:] if line.strip() and line.startswith('#') and 'UTC end:' in line][0]
                times.append( [start,end] )
                    

    ### SORTING
    if starttime_sort:
        data_to_sort = np.array( [line.split() for line in data] )
        sort_indices = data_to_sort[:,1].argsort()
        data         = np.array(data)
        data         = data[sort_indices]


    ### STATISTICS    
    data_split     = np.array( [line.split() for line in data] )
    total_LMST     = sum( [time[1].LMST_time-time[0].LMST_time for time in times] )

    glitch_all     = len(data)
    glitch_per_sol = len(data)*86400/total_LMST
    glitch_U       = len(data_split[ data_split[:,5]=='1'])
    glitch_V       = len(data_split[ data_split[:,6]=='1'])
    glitch_W       = len(data_split[ data_split[:,7]=='1'])
    glitch_Uonly   = len(data_split[(data_split[:,5]=='1') & (data_split[:,6]=='0') & (data_split[:,7]=='0')])
    glitch_Vonly   = len(data_split[(data_split[:,5]=='0') & (data_split[:,6]=='1') & (data_split[:,7]=='0')])
    glitch_Wonly   = len(data_split[(data_split[:,5]=='0') & (data_split[:,6]=='0') & (data_split[:,7]=='1')]) 
    try:
        glitch_indivdual_ratio = (glitch_Uonly+glitch_Vonly+glitch_Wonly)/glitch_all*100
    except ZeroDivisionError:
        glitch_indivdual_ratio = np.nan

    output1   = [u'MERGED:']
    output2   = ['  %s' % file for file in glitch_detector_files]
    output3   = [u'',
                 u'RESULTS:',      
                 u'  Glitches total:            %s'           % glitch_all,
                 u'       all / sol:            %.0f'         % glitch_per_sol,
                 u'            on U:            %s'           % glitch_U,
                 u'            on V:            %s'           % glitch_V,
                 u'            on W:            %s'           % glitch_W,
                 u'       only on U:            %s'           % glitch_Uonly,
                 u'       only on V:            %s'           % glitch_Vonly,
                 u'       only on W:            %s'           % glitch_Wonly,
                 u'    indi./ all %%:            %.1f'        % glitch_indivdual_ratio,
                 u'',
                 u'Done in:   %s (h:m:s).'                    % sec2hms( time.time()-now ),
                 u'Timestamp: %s'                             % datetime.datetime.now().strftime('%Y-%m-%dT%H:%M:%S')]

    # create good-looking output for later
    len_max_line  = max( [len(line) for line in output1+output2+output3] )
    formatter_str = "\n#   | %-" + str(len_max_line) + 's |'

    output1 = [formatter_str % line for line in output1]
    output2 = [formatter_str % line for line in output2]
    output3 = [formatter_str % line for line in output3]
    output1.insert(0,           '\n\n#   +-' +  '-' * len_max_line + '-+'   )
    output3.insert(len(output3),  '\n#   +-' +  '-' * len_max_line + '-+\n' )


    ### WRITING OUT
    counter = 0
    with open(outfile, 'w') as fp: 
                    
        for header_line in header:
            fp.write(header_line)

        for line in data:
            counter += 1 
            number   = '%06d' % counter 
            fp.write(number+line[6:])

        for line in output1+output2+output3:
            fp.write(line)


    ### OUTPUT
    print()
    print(u'Merged glitch file to:')
    print(outfile)
def read_config(config_file):

    """
    Read yaml config file via `full_load`. See details in:
    https://github.com/yaml/pyyaml/wiki/PyYAML-yaml.load(input)-Deprecation

    Included is also a hack to read numbers like `1e-7` as floats and not
    as strings. Shouldn't be like this, but is: See details in:
    https://stackoverflow.com/questions/30458977/yaml-loads-5e-6-as-string-and-not-a-number
    """

    print()
    print(u'Reading config file:')
    print(os.path.abspath(config_file))





    try:
        with open(config_file, 'r') as yamlfile:

            loader = yaml.FullLoader
            loader.add_implicit_resolver(
                u'tag:yaml.org,2002:float',
                re.compile(u'''^(?:
                 [-+]?(?:[0-9][0-9_]*)\\.[0-9_]*(?:[eE][-+]?[0-9]+)?
                |[-+]?(?:[0-9][0-9_]*)(?:[eE][-+]?[0-9]+)
                |\\.[0-9_]+(?:[eE][-+][0-9]+)?
                |[-+]?[0-9][0-9_]*(?::[0-5]?[0-9])+\\.[0-9_]*
                |[-+]?\\.(?:inf|Inf|INF)
                |\\.(?:nan|NaN|NAN))$''', re.X),
                list(u'-+0123456789.'))

            params = yaml.load(yamlfile, Loader=loader)

    except FileNotFoundError:
        print(u'ERROR: read_config: Config file not found.')
        sys.exit()

    return params
def sec2hms(seconds, digits=0):

    """
    Convert seconds given as float into hours, minutes and seconds.
    Optional 'digits' determines position after decimal point.

    Returns string.
    """


    string = str( datetime.timedelta(seconds=seconds) )
    parts  = string.split('.')

    try:
        frac_sec = parts[1]
        parts[1] = frac_sec[:digits]

    except IndexError:
        if digits>0:
            parts += [digits*'0']

    string = '.'.join(parts)
    string = string.rstrip('.')

    return string
def get_files(this, pattern='*'):

    """
    Return all files in dir 'this' (and sub-directories) if 'this' is dir,
    or return file (if 'this' is file). You may use 'pattern' to search only for
    certain file name patterns (eg: file endings, station names, ..).

    Script uses os.walk and thus lists all files in a directory, meaning also
    those in sub-directories. Bash wildcards are understood, python regex not. 

    :type this: str
    :param this: file or dir (bash wildcards work)

    :type pattern: str
    :param pattern: pattern to filter file(s) (bash wildcards work)
    
    :return: list of files (or empty list, if no matching files were found)
    """
    

    ### get files from 'this' to loop over
    uber  = glob.glob(this)
    files = []

    try:
        for u in uber:
            if os.path.isfile(u):
                files.append(u)

            else:
                # walk trough (sub)-directory and put ALL files into 'files'
                for root, dirs, files_walk in os.walk(u):
                    for file_walk in files_walk:

                        # take care on 'pattern'
                        file = os.path.join( root,file_walk )
                        if fnmatch.fnmatch( os.path.basename(file), pattern ) and not file in files:
                            files.append( file )

    except Exception:
        files = []

    return sorted(files)


# Ppol (P-wave Polarization)
class ppol():

    """
    TO DO:
        - data gaps
    """

    def __init__(self, comp_N=None, comp_E=None, comp_Z=None, stream=None, demean=True, **kwargs):
        traces       = self._assign_traces(comp_N=comp_N, comp_E=comp_E, comp_Z=comp_Z, stream=stream)
        self.trace_N = traces[0].copy()
        self.trace_E = traces[1].copy()
        self.trace_Z = traces[2].copy()

        if demean:
            self.trace_N.detrend('demean')
            self.trace_E.detrend('demean')
            self.trace_Z.detrend('demean')
            self.demeaned = True
        else:
            self.demeaned = False

        self._sanity_check()
        self.results = self.calc(**kwargs)
    def __str__(self, tag='', print_3D=True):

        """
        When calling print(...)
        """

        if print_3D:
            string =    u'\n' \
                      + u'PPOL CALCULATIONS: %s   \n' % tag \
                      + u'   %11s : %s            \n' % ('demeaned',   self.demeaned) \
                      + u'   %11s : %s            \n' % ('NPTS',       len(self.trace_N.data)) \
                      + u'   %11s : %-5.1f ± %3.1f\n' % ('BAZ_2D',     self.results[0], self.results[1]) \
                      + u'   %11s : %-5.1f ± %3.1f\n' % ('BAZ_3D',     self.results[2], self.results[3]) \
                      + u'   %11s : %-5.1f ± %3.1f\n' % ('INC_2D',     self.results[4], self.results[5]) \
                      + u'   %11s : %-5.1f ± %3.1f\n' % ('INC_3D',     self.results[6], self.results[7]) \
                      + u'   %11s : %-10.1g       \n' % ('SNR_HOR_2D', self.results[8]) \
                      + u'   %11s : %-10.1g       \n' % ('SNR_3D',     self.results[9]) \
                      + u'   %11s : %-10.1g       \n' % ('SNR_RZp_2D', self.results[10]) \
                      + u'   %11s : %-10.1g       \n' % ('SNR_RZp_3D', self.results[11]) \
                      + u'   %11s : %-10.3f       \n' % ('POL_HOR_2D', self.results[12]) \
                      + u'   %11s : %-10.3f       \n' % ('POL_3D',     self.results[13]) \
                      + u'   %11s : %-10.3f       \n' % ('POL_RZp_2D', self.results[14]) \
                      + u'   %11s : %-10.3f       \n' % ('POL_RZp_3D', self.results[15]) \
                      + u'   %11s : %s            \n' % ('EigVecs_2D', self.results[20][0]) \
                      + u'   %11s : %s            \n' % (''          , self.results[20][1]) \
                      + u'   %11s : %s            \n' % ('EigVals_2D', self.results[21]) \
                      + u'   %11s : %s            \n' % ('EigVecs_3D', self.results[22][0]) \
                      + u'   %11s : %s            \n' % (''          , self.results[22][1]) \
                      + u'   %11s : %s            \n' % (''          , self.results[22][2]) \
                      + u'   %11s : %s              ' % ('EigVals_3D', self.results[23])
        else:
            string =    u'\n' \
                      + u'PPOL CALCULATIONS: %s   \n' % tag \
                      + u'   %11s : %s            \n' % ('demeaned',   self.demeaned) \
                      + u'   %11s : %s            \n' % ('NPTS',       len(self.trace_N.data)) \
                      + u'   %11s : %-5.1f ± %3.1f\n' % ('BAZ_2D',     self.results[0], self.results[1]) \
                      + u'   %11s : %-10.1g       \n' % ('SNR_HOR_2D', self.results[8]) \
                      + u'   %11s : %-10.3f       \n' % ('POL_HOR_2D', self.results[12]) \
                      + u'   %11s : %s            \n' % ('EigVecs_2D', self.results[20][0]) \
                      + u'   %11s : %s            \n' % (''          , self.results[20][1]) \
                      + u'   %11s : %s              ' % ('EigVals_2D', self.results[21])

        return string
    def _assign_traces(self, comp_N=None, comp_E=None, comp_Z=None, stream=None):

        """
        """

        if stream is not None:
            stream_N = stream.select(component='N') or stream.select(component='Q') or stream.select(component='U') or stream.select(component='1') or stream.select(component='T')
            stream_E = stream.select(component='E') or stream.select(component='T') or stream.select(component='V') or stream.select(component='2') or stream.select(component='R')
            stream_Z = stream.select(component='Z') or stream.select(component='L') or stream.select(component='W') or stream.select(component='3')

            if not stream_N and not stream_E:
                print(u'No idea how to perform polarization analysis on components: %s' % ', '.join( [tr.id for tr in stream] ))
                sys.exit()

            else:
                trace_N = stream_N[0]
                trace_E = stream_E[0]
                if stream_Z is not None:
                    trace_Z = stream_Z[0]
                else:
                    trace_Z = Trace()


        elif comp_N is not None and comp_E is not None:

            if isinstance(comp_N, Trace):
                trace_N = comp_N
            elif isinstance(comp_N, (list, tuple, np.ndarray)):
                trace_N = Trace(data=np.asarray(comp_N), header={'channel':'N'})
            else:
                print(u'`comp_N` must be either a ObsPy `Trace` object or a `list`-like object containing your waveform data.')
                sys.exit()

            if isinstance(comp_E, Trace):
                trace_E = comp_E
            elif isinstance(comp_E, (list, tuple, np.ndarray)):
                trace_E = Trace(data=np.asarray(comp_E), header={'channel':'E'})
            else:
                print(u'`comp_E` must be either a ObsPy `Trace` object or a `list`-like object containing your waveform data.')
                sys.exit()

            if comp_Z is not None:
                if isinstance(comp_Z, Trace):
                    trace_Z = comp_Z
                elif isinstance(comp_Z, (list, tuple, np.ndarray)):
                    trace_Z = Trace(data=np.asarray(comp_Z), header={'channel':'Z'})
                else:
                    print(u'`comp_Z` must be either a ObsPy `Trace` object or a `list`-like object containing your waveform data.')
                    sys.exit()                    
            else:
                trace_Z = Trace()


        else:
            print(u'You must either specify an ObsPy `stream` object containing at least two horizontal components')
            print(u'or `comp_N` and `comp_E` (both either as ObsPy `Trace` objects or lists).')
            sys.exit()


        return trace_N , trace_E, trace_Z
    def _sanity_check(self):

        
        len_N   = self.trace_N.stats.npts
        len_E   = self.trace_E.stats.npts
        len_Z   = self.trace_Z.stats.npts

        start_N = self.trace_N.stats.starttime
        start_E = self.trace_E.stats.starttime
        start_Z = self.trace_Z.stats.starttime

        if len_Z != 0:
            if len_N!=len_E or len_N!=len_Z or len_E!=len_Z:
                print(u'ERROR: Data lengths of components are not equal. Cannot perform ppol analysis.')
                sys.exit()

        else:
            if len_N!=len_E:
                print(u'ERROR: Data lengths of components are not equal. Cannot perform ppol analysis.')
                sys.exit()


        if len_Z != 0:
            if start_N!=start_E or start_N!=start_Z or start_E!=start_Z:
                print(u'WARNING: Data do not start at the same time. Analysis is performed nevertheless.')
                print('  '+self.trace_N.__str__())
                print('  '+self.trace_E.__str__())
                print('  '+self.trace_Z.__str__())

        else:
            if start_N!=start_E:
                print(u'WARNING: Data do not start at the same time. Analysis is performed nevertheless.')
                print('  '+self.trace_N.__str__())
                print('  '+self.trace_E.__str__())
    def calc(self, bias=False, fix_angles='EQ', Xoffset_samples_for_amplitude=None, **kwargs):

        r"""
        DO INDIVIDUAL PPOL MEASUREMENTS

        Take data arrays and perform 2-D and 3-D principle component
        analysis (2-D and 3-D PCA). Components must be of equal length
        and correspond to same start and thus end times.


        Useful for seismology, for example.

        |       ^  
        |       \| data_1
        |       \|
        |       \|
        |       \|
        |       x------------> data_2
        |     data_Z 
        |  (pointing to you)


        180° ambiguity:
        
          P-waves  -->  2-D data  -->  BAZ unknown  -->  180° ambiguity    (but if one knows expected first motion)
                   -->  3-D data  -->  BAZ unknown  -->  no 180° ambiguity (because P-wave must arrive from below, i.e. INC must >0°, or because of known, expected first motion)
          R-waves  -->  3-D data  -->  BAZ unknown  -->  no 180° ambiguity (retrogradicity, however, this algoirthmus does not treat surface / Rayleigh waves).

        Note
        ----

            An unknown event BAZ is the same as an unknown station orientation
            with a known event BAZ.

        This function does not demean data before running.
        """



        ### assign needed data
        data_1 = self.trace_N.data
        data_2 = self.trace_E.data
        data_Z = self.trace_Z.data



        ### IF NO Z-COMPONENT GIVEN
        if data_Z is None:
     

            ### 2-D phi, horizontal plane
            covariance_matrix_2D_hori      = covariance_matrix(data_1, data_2, bias=bias)                                  # does be default no demeaning internally!
            eig_val_2D, eig_vec_2D         = np.linalg.eig(covariance_matrix_2D_hori)
            index_array_descending_eig_val = eig_val_2D.argsort()[::-1]
            eig_val_2D                     = eig_val_2D[index_array_descending_eig_val]                                    # E-values descending
            eig_vec_2D                     = eig_vec_2D[:,index_array_descending_eig_val]                                  # E-vectors sorted acc. to E-values
            eig_vec_2D_1                   = eig_vec_2D[:,0]
            
            # Derived
            PHI_2D                         = (np.arctan2( eig_vec_2D_1[1], eig_vec_2D_1[0] ) * 180/np.pi )%360
            PHI_err_2D                     = np.arctan( np.sqrt( eig_val_2D[1]/eig_val_2D[0] )) * 180/np.pi
            SNR_HOR_2D                     = (eig_val_2D[0] - eig_val_2D[1]) / eig_val_2D[1]                               # De Meersman et al. (2006)
            POL_HOR_2D                     = 1 - eig_val_2D[1]/eig_val_2D[0]                                               # rectilinearity in horizontal plane (1 for linearised, 0 for circular polarisation). Jurkevics (1988)
            data_R_2D, data_T_2D           = rotate.rotate_ne_rt(data_1, data_2, PHI_2D)
            
            # Others
            PHI_3D                         = np.nan
            PHI_err_3D                     = np.nan
            INC_2D                         = np.nan
            INC_err_2D                     = np.nan
            INC_3D                         = np.nan
            INC_err_3D                     = np.nan
            SNR_3D                         = np.nan
            SNR_RZp_2D                     = np.nan
            SNR_RZp_3D                     = np.nan
            POL_3D                         = np.nan
            POL_RZp_2D                     = np.nan
            POL_RZp_3D                     = np.nan
            data_R_3D                      = np.nan*data_R_2D
            data_T_3D                      = np.nan*data_T_2D
            eig_val_3D                     = np.nan*np.ones(3)                                                            # E-values descending
            eig_vec_3D                     = np.nan*np.ones(3)                                     


            ### AMBIGUITY 180°
            if fix_angles.upper()=='EQ':         # Nothing we can do solving the 180° ambiguity in 2D (this code doesn't make use of first motion information)
                pass


            elif fix_angles.upper()=='AMP':
                data_1_offest = data_1[Xoffset_samples_for_amplitude:]
                data_2_offest = data_2[Xoffset_samples_for_amplitude:]
                
                amp_1         = data_1[np.argmax(np.abs(data_1_offest))]
                amp_2         = data_2[np.argmax(np.abs(data_2_offest))]
                
                PHI_2D_OLD    = PHI_2D

                # 2-D PHI
                if abs(amp_1)>=abs(amp_2):              # `data_1` more significant than `data_2` 
                    if amp_1>=0:                        # `data_1` positive
                        if PHI_2D>=90 and PHI_2D<180 or PHI_2D>=270 and PHI_2D<360:
                            PHI_2D = PHI_2D%180
                        else:
                            PHI_2D = PHI_2D%180 + 180
                    else:                               # `data_1` negative
                        if PHI_2D>=90 and PHI_2D<180 or PHI_2D>=270 and PHI_2D<360:
                            PHI_2D = PHI_2D%180 + 180
                        else:
                            PHI_2D = PHI_2D%180
                else:                                   # `data_2` more significant than `data_1` 
                    if amp_2>=0:                        # `data_2` positive
                        PHI_2D = PHI_2D%180 + 180
                    else:                               # `data_2` negative
                        PHI_2D = PHI_2D%180             

                # correct radial and transverse data
                data_R_2D, data_T_2D = rotate.rotate_ne_rt(data_1, data_2, PHI_2D)

                # 2-D INC
                if PHI_2D_OLD != PHI_2D:
                    INC_2D *= -1
                INC_2D = INC_2D % 180


            else:
                pass



        ### IF Z-COMPONENT GIVEN
        else:


            ### 2-D phi, horizontal plane
            covariance_matrix_2D_hori      = covariance_matrix(data_1, data_2, bias=bias)                               # does be default no demeaning internally!
            eig_val_2D, eig_vec_2D         = np.linalg.eig(covariance_matrix_2D_hori)
            index_array_descending_eig_val = eig_val_2D.argsort()[::-1]
            eig_val_2D                     = eig_val_2D[index_array_descending_eig_val]                                 # E-values descending
            eig_vec_2D                     = eig_vec_2D[:,index_array_descending_eig_val]                               # E-vectors sorted acc. to E-values
            eig_vec_2D_1                   = eig_vec_2D[:,0]
            
            # Derived
            PHI_2D                         = (np.arctan2( eig_vec_2D_1[1], eig_vec_2D_1[0] ) * 180/np.pi )%360
            PHI_err_2D                     = np.arctan( np.sqrt( eig_val_2D[1]/eig_val_2D[0] )) * 180/np.pi
            SNR_HOR_2D                     = (eig_val_2D[0] - eig_val_2D[1]) / eig_val_2D[1]                            # De Meersman et al. (2006)
            POL_HOR_2D                     = 1 - eig_val_2D[1]/eig_val_2D[0]                                            # rectilinearity in horizontal plane (1 for linearised, 0 for circular polarisation). Jurkevics (1988)
            data_R_2D, data_T_2D           = rotate.rotate_ne_rt(data_1, data_2, PHI_2D)


            ### 2-D phi, radial-vertical plane
            data_R_2D, data_T_2D             = rotate.rotate_ne_rt(data_1, data_2, PHI_2D)
            covariance_matrix_2D_radZ        = covariance_matrix(data_Z, data_R_2D, bias=bias)
            eig_val_2D_radZ, eig_vec_2D_radZ = np.linalg.eig(covariance_matrix_2D_radZ)
            index_array_descending_eig_val   = np.argsort( eig_val_2D_radZ )[::-1]
            eig_val_2D_radZ                  = eig_val_2D_radZ[index_array_descending_eig_val]                          # E-values descending
            eig_vec_2D_radZ                  = eig_vec_2D_radZ[:,index_array_descending_eig_val]                        # E-vectors sorted acc. to E-values
            eig_vec_2D_radZ_1                = eig_vec_2D_radZ[:,0]

            # derived
            INC_2D                           = np.arctan( eig_vec_2D_radZ_1[1] / eig_vec_2D_radZ_1[0] ) * 180/np.pi
            INC_err_2D                       = np.arctan( np.sqrt( eig_val_2D_radZ[1]/eig_val_2D_radZ[0] )) * 180/np.pi
            SNR_RZp_2D                       = (eig_val_2D_radZ[0] - eig_val_2D_radZ[1]) / eig_val_2D_radZ[1]           # De Meersman et al. (2006)
            POL_RZp_2D                       = 1 - eig_val_2D_radZ[1]/eig_val_2D_radZ[0]                                # rectilinearity in radial-vertical plane (1 for linearised, 0 for circular polarisation). Jurkevics (1988)


            ### 3-D
            covariance_matrix_3D             = covariance_matrix(data_1, data_2, data_Z, bias=bias)                     # does be default no demeaning internally!
            eig_val_3D, eig_vec_3D           = np.linalg.eig(covariance_matrix_3D)
            index_array_descending_eig_val   = np.argsort( eig_val_3D )[::-1]
            eig_val_3D                       = eig_val_3D[index_array_descending_eig_val]                               # E-values descending
            eig_vec_3D                       = eig_vec_3D[:,index_array_descending_eig_val]                             # E-vectors sorted acc. to E-values
            eig_vec_3D_1                     = eig_vec_3D[:,0]

            # derived
            PHI_3D                           = (np.arctan2( eig_vec_3D_1[1], eig_vec_3D_1[0] ) * 180/np.pi) % 360
            PHI_err_3D                       = np.arctan( np.sqrt( eig_val_3D[2]/(eig_val_3D[1]+eig_val_3D[0]) )) * 180/np.pi
            SNR_3D                           = abs((eig_val_3D[0] - (eig_val_3D[1] + eig_val_3D[2]))) / (eig_val_3D[1] + eig_val_3D[2])  # De Meersman et al. (2006)        
            POL_3D                           = 1 - ( eig_val_3D[1]+eig_val_3D[2] )/( 2*eig_val_3D[0] )                  # rectilinearity in 3-D. (1 for linearised, 0 for circular polarisation). Jurkevics (1988)
            data_R_3D, data_T_3D             = rotate.rotate_ne_rt(data_1, data_2, PHI_3D)


            ### 3-D phi, radial & Z-data plane
            data_R_3D, data_T_3D             = rotate.rotate_ne_rt(data_1, data_2, PHI_3D)
            covariance_matrix_3D_radZ        = covariance_matrix(data_Z, data_R_3D, bias=bias)                          # does be default no demeaning internally!
            eig_val_3D_radZ, eig_vec_3D_radZ = np.linalg.eig(covariance_matrix_3D_radZ)
            index_array_descending_eig_val   = np.argsort( eig_val_3D_radZ )[::-1]
            eig_val_3D_radZ                  = eig_val_3D_radZ[index_array_descending_eig_val]                          # E-values descending
            eig_vec_3D_radZ                  = eig_vec_3D_radZ[:,index_array_descending_eig_val]                        # E-vectors sorted acc. to E-values
            eig_vec_3D_radZ_1                = eig_vec_3D_radZ[:,0]

            # derived
            INC_3D                           = np.arctan( eig_vec_3D_radZ_1[1] / eig_vec_3D_radZ_1[0] ) * 180/np.pi
            INC_err_3D                       = np.arctan( np.sqrt( eig_val_3D_radZ[1]/eig_val_3D_radZ[0] )) * 180/np.pi
            SNR_RZp_3D                       = (eig_val_3D_radZ[0] - eig_val_3D_radZ[1]) / eig_val_3D_radZ[1]           # De Meersman et al. (2006)
            POL_RZp_3D                       = 1 - eig_val_3D_radZ[1]/eig_val_3D_radZ[0]                                # rectilinearity in radial-vertical plane (1 for linearised, 0 for circular polarisation). Jurkevics (1988)


            ### AMBIGUITY 180°
            if fix_angles.upper()=='EQ':    # However, correct baz must deliver incidence angle>0, therefore can solve ambiguity
                if INC_2D < 0:
                    PHI_2D               = (PHI_2D+180)%360
                    INC_2D               = abs(INC_2D)
                    data_R_2D, data_T_2D = rotate_2D(data_R_2D, data_T_2D, 180)
                if INC_3D < 0:
                    PHI_3D               = (PHI_3D+180)%360
                    INC_3D               = abs(INC_3D)
                    data_R_3D, data_T_3D = rotate_2D(data_R_3D, data_T_3D, 180)


            elif fix_angles.upper()=='AMP':
                data_1_offest = data_1[Xoffset_samples_for_amplitude:]
                data_2_offest = data_2[Xoffset_samples_for_amplitude:]
                data_Z_offest = data_Z[Xoffset_samples_for_amplitude:]
                
                amp_1         = data_1[np.argmax(np.abs(data_1_offest))]
                amp_2         = data_2[np.argmax(np.abs(data_2_offest))]
                amp_Z         = data_Z[np.argmax(np.abs(data_Z_offest))]
                
                PHI_2D_OLD    = PHI_2D
                PHI_3D_OLD    = PHI_3D

                # 2-D PHI
                if abs(amp_1)>=abs(amp_2):              # `data_1` more significant than `data_2` 
                    if amp_1>=0:                        # `data_1` positive
                        if PHI_2D>=90 and PHI_2D<180 or PHI_2D>=270 and PHI_2D<360:
                            PHI_2D = PHI_2D%180
                        else:
                            PHI_2D = PHI_2D%180 + 180
                    else:                               # `data_1` negative
                        if PHI_2D>=90 and PHI_2D<180 or PHI_2D>=270 and PHI_2D<360:
                            PHI_2D = PHI_2D%180 + 180
                        else:
                            PHI_2D = PHI_2D%180
                else:                                   # `data_2` more significant than `data_1` 
                    if amp_2>=0:                        # `data_2` positive
                        PHI_2D = PHI_2D%180 + 180
                    else:                               # `data_2` negative
                        PHI_2D = PHI_2D%180             
    
                # 3-D PHI
                if abs(amp_1)>=abs(amp_2):              # `data_1` more significant than `data_2` 
                    if amp_1>=0:                        # `data_1` positive
                        if PHI_3D>=90 and PHI_3D<180 or PHI_3D>=270 and PHI_3D<360:
                            PHI_3D = PHI_3D%180
                        else:
                            PHI_3D = PHI_3D%180 + 180
                    else:                               # `data_1` negative
                        if PHI_3D>=90 and PHI_3D<180 or PHI_3D>=270 and PHI_3D<360:
                            PHI_3D = PHI_3D%180 + 180
                        else:
                            PHI_3D = PHI_3D%180
                else:                                   # `data_2` more significant than `data_1` 
                    if amp_2>=0:                        # `data_2` positive
                        PHI_3D = PHI_3D%180 + 180
                    else:                               # `data_2` negative
                        PHI_3D = PHI_3D%180   

                # correct radial and transverse data
                data_R_2D, data_T_2D = rotate.rotate_ne_rt(data_1, data_2, PHI_2D)
                data_R_3D, data_T_3D = rotate.rotate_ne_rt(data_1, data_2, PHI_3D)

                # 2-D INC
                if PHI_2D_OLD != PHI_2D:
                    INC_2D *= -1
                INC_2D = INC_2D % 180

                # 3-D INC
                if PHI_3D_OLD != PHI_3D:
                    INC_3D *= -1
                INC_3D = INC_3D % 180


            else:
                pass
        


        ### RESULTS
        results = PHI_2D,     PHI_err_2D, PHI_3D,     PHI_err_3D, \
                  INC_2D,     INC_err_2D, INC_3D,     INC_err_3D, \
                  SNR_HOR_2D, SNR_3D,     SNR_RZp_2D, SNR_RZp_3D, \
                  POL_HOR_2D, POL_3D,     POL_RZp_2D, POL_RZp_3D, \
                  data_R_2D,  data_T_2D,  data_R_3D,  data_T_3D,  \
                  eig_vec_2D, eig_val_2D, eig_vec_3D, eig_val_3D



        ### ASSIGN / RETURN RESULTS
        self.results = results
        return self.results
    def plot(self, title='', verticals=(), outfile=None, show=True, **kwargs):

        """
        PLOT INDIVIDUAL PPOL MEASUREMENT


        Either provide num,py lists or Obspy traces.
        Always plots 2-D angle
        Expects that all traces have same start time and same amount of samples!
        """


        ### assign needed labels and data
        label_1 = '<--  %s  -->' %  self.trace_N.stats.channel
        label_2 = '<--  %s  -->' %  self.trace_E.stats.channel
        label_Z = '<--  %s  -->' %  self.trace_Z.stats.channel
        label_R = '<--  %s  -->' % (self.trace_E.stats.channel[:-1] + 'R')
        
        data_1  = self.trace_N.data
        data_2  = self.trace_E.data
        data_Z  = self.trace_Z.data


        ### Ppol result (some needed for plotting)
        #   only measures based on BAZ_2D printed, because they are always available
        BAZ_2D     = self.results[0]
        BAZ_err_2D = self.results[1]
        INC_2D     = self.results[4]
        INC_err_2D = self.results[5]
        data_R     = self.results[16]   # comp_R_2D
        eigvecs    = self.results[22]   # eig_vec_3D


        ### PLOT WAVEFORMS
        # Figure instance
        fig = plt.figure(figsize=(12,8))
        fig.canvas.set_window_title('Ppol plot individual measurement: %s' % title)
        fig.suptitle(title, fontsize=11)
        gs  = fig.add_gridspec(2, 3)

        ax1 = fig.add_subplot(gs[0, :])
        ax1 = quick_plot(self.trace_N, self.trace_E, self.trace_Z, verts=verticals, ylabel='Amplitudes', xlabel='Time', legend_loc='upper right', axis=ax1)


        ### SMALL HACKS to make sure for small amplitudes everything works out (precisely: arrow head of angle indicators)
        factor      = 1
        factor_str  = ''

        if data_Z is not None:
            if max( [max(abs(data_Z)), max(abs(data_1)), max(abs(data_2)), max(abs(data_R))] ) <= 1e-3:
                factor      = 1e9
                factor_str  = '%.0g' % (factor**-1)
                data_Z     *= factor
                data_1     *= factor
                data_2     *= factor
                data_R     *= factor

        else:
            if max( [max(abs(data_1)), max(abs(data_2)), max(abs(data_R))] ) <= 1e-3:
                factor      = 1e9
                factor_str  = '%.0g' % (factor**-1)
                data_1     *= factor
                data_2     *= factor


        ### PLOT PPOL
        # Variables needed for (beautiful) plotting
        colours = [[0, 0, 1-0.6*i/len(data_1)] for i in range(len(data_2))]
        
        maxi    = max( list(np.abs(data_1))+list(np.abs(data_2)) ) * 1.05
        
        xx      = np.linspace(-maxi, maxi, 100)
        yy      = 1/np.tan(BAZ_2D*np.pi/180) * xx
        yy_be   = 1/np.tan((BAZ_2D-BAZ_err_2D)*np.pi/180) * xx
        yy_af   = 1/np.tan((BAZ_2D+BAZ_err_2D)*np.pi/180) * xx
        
        x       = maxi/2*np.sin( (BAZ_2D-7)*np.pi/180 )
        y       = maxi/2*np.cos( (BAZ_2D-7)*np.pi/180 )
        x2      = maxi/2*np.sin( (BAZ_2D+2)*np.pi/180 )
        y2      = maxi/2*np.cos( (BAZ_2D+2)*np.pi/180 )

        # Set-up
        ax2 = fig.add_subplot(gs[1, 0])
        ax2.grid(ls='-.', lw=0.5)
        ax2.set(xlabel=label_2, ylabel=label_1)
        ax2.set_ylim([-maxi, maxi])
        ax2.set_xlim([-maxi, maxi])
        ax2.set_aspect('equal', adjustable='box')
        ax2.text(1, 0, title+' (%s points)' % len(data_1), ha='right', color='red', transform=ax2.transAxes, fontsize=6, bbox=dict(boxstyle='round,pad=0', fc='white', ec="white", lw=0))
        ax2.text(0, 1, factor_str, ha='center', color='black', transform=ax2.transAxes, fontsize=9, bbox=dict(boxstyle='round,pad=0', fc='white', ec="white", lw=0))
        
        ax2.spines['right'].set_color('none')
        ax2.spines['top'].set_color('none')
        ax2.spines['left'].set_color('none')
        ax2.spines['bottom'].set_color('none')  
        
        # Plot commands; data & Phi angle + erorrs
        if (BAZ_2D+BAZ_err_2D)//180 != BAZ_2D//180.001 or BAZ_2D == 0:
            ax2.fill_betweenx(yy_af, xx, facecolor='red', alpha=0.075)
            ax2.fill_betweenx(yy,    xx, facecolor='red', alpha=0.075)
        else:
            ax2.fill_between(xx, yy, yy_af, facecolor='red', alpha=0.075)       # error area of PHI
        if (BAZ_2D-BAZ_err_2D)//180 != BAZ_2D//180.001 or BAZ_2D == 0:
            ax2.fill_betweenx(yy_be, xx, facecolor='red', alpha=0.075)
            ax2.fill_betweenx(yy,    xx, facecolor='red', alpha=0.075)
        else:
            ax2.fill_between(xx, yy_be, yy, facecolor='red', alpha=0.075)       # error area of PHI
        ax2.plot([-maxi, maxi], [0, 0], 'k', lw=1)                              # centered horizontal line
        ax2.plot([0, 0], [-maxi, maxi], 'k', lw=1)                              # centered vertical line
        ax2.scatter(data_2, data_1, s=7, c=colours, zorder=3)                             # data
        ax2.plot( xx, yy, 'indianred', lw=1.5, zorder=4)                                  # PHI results
                
        # Angle arc + arrow head
        ax2.add_patch( Arc([0,0], maxi,  maxi, 90, -BAZ_2D, 0, color='indianred', lw=1.5, zorder=5))
        a = FancyArrowPatch([x,y], [x2,y2], mutation_scale=20, lw=1.5, arrowstyle="-|>", color="indianred", zorder=6)
        ax2.add_artist(a)       

        # Legend
        scatter_proxy = mpl.lines.Line2D([0],[0], c="indianred", marker='>')
        ax2.legend([scatter_proxy], ['BAZ_2D=(%.1f\u00B1%.1f)\u00b0' % (BAZ_2D, BAZ_err_2D)], numpoints=1, loc='upper right', prop={'size': 8})


        ## plot axis (part) 2 and 3, if there's Z-data
        if data_Z is not None:

            # Variables needed for (beautiful) plotting
            maxi  = max( list(np.abs(data_R))+list(np.abs(data_Z)) ) * 1.05
            
            xx    = np.linspace(-maxi, maxi, 100)
            yy    = 1/np.tan( INC_2D*np.pi/180) * xx
            yy_be = 1/np.tan((INC_2D-INC_err_2D)*np.pi/180) * xx
            yy_af = 1/np.tan((INC_2D+INC_err_2D)*np.pi/180) * xx
            
            x     = maxi/2*np.sin( (INC_2D-7)*np.pi/180 )
            y     = maxi/2*np.cos( (INC_2D-7)*np.pi/180 )
            x2    = maxi/2*np.sin( (INC_2D+2)*np.pi/180 )
            y2    = maxi/2*np.cos( (INC_2D+2)*np.pi/180 )
            
            # Set-up
            ax3 = fig.add_subplot(gs[1, 1])
            ax3.grid(ls='-.', lw=0.5, zorder=-1)
            ax3.set(xlabel=label_R, ylabel=label_Z)
            ax3.set_ylim([-maxi, maxi])
            ax3.set_xlim([-maxi, maxi])
            ax3.set_aspect('equal', adjustable='box')
            ax3.text(1, 0, title+' (%s points)' % len(data_2), horizontalalignment='right', color='red', transform=ax3.transAxes, fontsize=6, bbox=dict(boxstyle='round,pad=0', fc='white', ec="white", lw=0))
            ax3.text(0, 1, factor_str, ha='center', color='black', transform=ax3.transAxes, fontsize=9, bbox=dict(boxstyle='round,pad=0', fc='white', ec="white", lw=0))
        
            ax3.spines['right'].set_color('none')
            ax3.spines['top'].set_color('none')
            ax3.spines['left'].set_color('none')
            ax3.spines['bottom'].set_color('none')  
        
            # Plot commands; data & Phi angle + erorrs
            if (INC_2D+INC_err_2D)//180 != INC_2D//180.001 or INC_2D == 0:
                ax3.fill_betweenx(yy_af, xx, facecolor='red', alpha=0.075)
                ax3.fill_betweenx(yy,    xx, facecolor='red', alpha=0.075)
            else:
                ax3.fill_between(xx, yy, yy_af, facecolor='red', alpha=0.075)           # error area of INC
            if (INC_2D-INC_err_2D)//180 != INC_2D//180.001 or INC_2D == 0:
                ax3.fill_betweenx(yy_be, xx, facecolor='red', alpha=0.075)
                ax3.fill_betweenx(yy,    xx, facecolor='red', alpha=0.075)
            else:
                ax3.fill_between(xx, yy_be, yy, facecolor='red', alpha=0.075)           # error area of INC         
            ax3.plot([-maxi, maxi], [0, 0],  'k', lw=1)
            ax3.plot( [0, 0], [-maxi, maxi], 'k', lw=1)
            ax3.scatter(data_R, data_Z, s=5, c=colours, zorder=3)
            ax3.plot( xx, yy, 'indianred', lw=1.5, zorder=4)

            # Angle arc + arrow head
            ax3.add_patch( Arc([0,0], maxi,  maxi, 90, -INC_2D, 0, color='indianred', lw=1.5, zorder=5))
            a = FancyArrowPatch([x,y], [x2,y2], mutation_scale=20, lw=1.5, arrowstyle="-|>", color="indianred", zorder=6)
            ax3.add_artist(a)
            
            # Legend
            scatter_proxy = mpl.lines.Line2D([0],[0], c="indianred", marker='>')
            ax3.legend([scatter_proxy], ['INC_2D=(%.1f\u00B1%.1f)\u00b0' % (INC_2D, INC_err_2D)], loc='upper right', prop={'size': 8})



            ## 3-D plot
            # Variables needed for (beautiful) plotting
            mean_N   = np.average(data_1)
            mean_E   = np.average(data_2)
            mean_Z   = np.average(data_Z)
            maxi     = max( list(np.abs(data_Z))+list(np.abs(data_1))+list(np.abs(data_2)) ) * 1.05
            max_dist = np.sqrt(3*(maxi*2)**2)

            # Set-up
            ax4 = fig.add_subplot(gs[1, 2], projection='3d')
            ax4.set_xlabel(label_2, labelpad=5, fontsize=10)
            ax4.set_ylabel(label_1, labelpad=5, fontsize=10)
            ax4.zaxis.set_rotate_label(False)    # disable automatic rotation
            ax4.set_zlabel(label_Z, labelpad=5, fontsize=10, rotation=90)
            ax4.set_xlim( [-maxi,maxi] )
            ax4.set_ylim( [-maxi,maxi] )
            ax4.set_zlim( [-maxi,maxi] )
            ax4.view_init(elev=15, azim=-135)
            ax4.text(0, 1, 1, factor_str, ha='center', color='black', transform=ax4.transAxes, fontsize=9, bbox=dict(boxstyle='round,pad=0', fc='white', ec="white", lw=0))


            # Plot commands
            ax4.scatter(data_2, data_1, data_Z, c=colours, s=5, depthshade=False, zorder=3)
            ax4.scatter([mean_E], [mean_N], [mean_Z], s=20, c='white', edgecolors='k', depthshade=False, zorder=4)

            # eigenvectors
            eigvecs  = [ unit_vector(eigvecs[:,0]), unit_vector(eigvecs[:,1]), unit_vector(eigvecs[:,2]) ]
            for l, eigvec in enumerate( eigvecs ):
                length = (0.7-l*0.15) * max_dist/2      # arbitrary length, not related to actual eigen values
                a = _Arrow3D(np.asarray( [mean_E, eigvec[1]*length] ),
                             np.asarray( [mean_N, eigvec[0]*length] ),
                             np.asarray( [mean_Z, eigvec[2]*length] ),
                             mutation_scale=20, lw=1.5, arrowstyle="-|>", color="indianred", zorder=5)
                ax4.add_artist(a)

            # legend
            scatter_proxy = mpl.lines.Line2D([0],[0], c="indianred", marker='>')
            ax4.legend([scatter_proxy], ['eigen vectors'], loc='upper right', prop={'size': 8})


        ## Figure save / close
        plt.tight_layout(rect=[0, 0.02, 1, 0.98])
        if outfile:
            plt.savefig(outfile)
        if show:
            plt.show()
        plt.close()
    def display(self, tag='', print_3D=True):

        """
        Print the results of `calc`.
        """

        print(self.__str__(tag=tag, print_3D=print_3D))


# Plot related
def quick_plot(*y, x=None, data_labels=(), lw=1.5, win_title='', title='', xlabel='Data points', ylabel='Amplitude', x_invert=False, y_invert=False, xlim=None, ylim=None, verts=None, horis=None, legend_loc='best', axis=False, outfile=None, show=True, keep_open=False):

    """
    This function allows for some convenient, quick 2-D plotting.
    
    It requires loaded classes `Trace`, `Stream`,
    `UTCDateTime`, and `datetime.datetime`.
    
    Otherwise, it should be very user friendly. Times will
    be displayed as '%d/%m %H:%M:%S', which is harcoded.


    PARAMETERS
    ----------

    `y` : - y-data to be plotted, as many as you like
          - if your data are `Trace` object(s), it is plotted
            with date strings (if no `x` is given).
          - if your data are a `Stream` object, it is plotted
            with their internal plot functions and further parameters
            are not passed.
          - if none of above, then x-data are determined using np.arange
    `x` : - this overrides x-data (dates) when passing `Trace` objects as `y`,
            so you can assign anything as x-data and still conveniently use
            `Trace` objects.
          - `UTCDateTime`, `datetime.datetime`, and `str` objects will
             be converted to `UTCDateTime` and then matplotlib times,
             then plotted as date strings '%d/%m %H:%M:%S'
    `data_labels` : - will appear in legend
                    - should have as many elements as you pass in `y`,
                      otherwise it will not used but data simply numbered
    `lw`       : - linewidth of plots (lines are harcoded, no option to plot points only)
    `win_title`: - title of figures canvas (does not appear in plot)
    `title`    : - sup-title of figure instance (appears as title in plot)
    `xlabel    : - will appear beneath x-axis
    `ylabel    : - will appear next to y-axis
    `x_invert` : - will invert direction of x-axis, if True
    `y_invert` : - will invert direction of y-axis, if True
    `xlim`     : - list (or similar) with two elements. Choose same type as `y`. If 'xlim' is set, `ylim` is set automatically as to the y-range within these xlim
    `ylim`     : - list (or similar) with two elements.
    `verts`    : - either list or list of lists
                 - elements will be interpreted as x-data where vertical lines shall be plotted
                 - each list gets an own colour for all elements
                 - when elements are: - `UTCDateTime`, `datetime.datetime`, and `str` 
                                         --> conversion to matplotlib times (so x-data should be times, too)
                                      - `float` or `int`
                                         --> where you would expect it to happen (each sample on x-axis numbered).
                                             If x-data are times (e.g. matplotlib), you can pass here matplotlib times as well
                                      - `complex` (e.g. 3+0j)
                                         --> real part will be interpreted as relative position on x-axis,
                                             so (0.5+0j) will plot a vertical line in the middle of the x-axis.
                                             Imaginery part is neglected.
    `legend_loc`: - position of legend, by default automaically determined. Choose e.g. 'upper left' or 'lower right' or ...
    `axis`     : - if passed, all will be plotted into this axis and returned and next two options are ignored.                                            
    `outfile`  : - (absolut) path where plot shall be saved. Use endings like '.pdf' or '.png'
    `show`     : - if True (default), plot will be shown, otherwise not until the next show=True command ;)
                 - if False, matplotlib figure instance of is returned so it can be further used
                 - no matter the status of `show`, if an `outfile` is specified a plot will be saved

    NOTES
    -----
      - Matplotlib internally uses matplotlib times, always.
        (example: 733643.01392361, type=numpy.float64)

        This holds true also when extracting times, for example via
        ax.get_ylim() or similar, no matter the original input format. Therefore,
        this script converts `str`, `UTCDateTime` and `datetime.datetime` objects all to matplotb times
        before the plot command (although matplotlib can plot datetime.datetime objects natively), so
        these converted times later can be used to adjust xlim and ylim (if specified).
      - natively, matplotlib plots `datetime.datetime` objects with labels
        hour, minute, and second. This is true for both ax.plot and ax.plot_date ..
      - as I wished a personal format including month, day, hour, minute and second, and
        as times are converted to matplotlib times regardless (which matplotlib then does not
        plot with a nice string but as simples numbers), the following two commands are used:
        
        myFmt = mdates.DateFormatter('%d/%m %H:%M:%S')
        ax.xaxis.set_major_formatter(myFmt)
    """



    ### Figure instance & settings
    if not axis:
        fig = plt.figure(figsize=(16,10), num=title)
        fig.canvas.set_window_title(win_title)
        ax = fig.add_subplot(111)
    else:
        ax = axis

    ax.ticklabel_format(axis='y', style='sci', scilimits=(-3, 3))
    ax.set_title(title, fontsize=13)
    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.grid(ls='-.', lw=0.5)
    if y_invert:
        ax.invert_yaxis()
    if x_invert:
        ax.invert_xaxis()



    ### Labels
    if data_labels is not None:
        if len(data_labels)!=len(y):    # not same amount of entries (also true if none given)
            
            if all(isinstance(ele, Trace) for ele in y):
                labels = [trace.id if trace.data.any() else None for trace in y]    # empty data don't get a label

            else:
                labels = ['Data %s' % (i+1) for i in range(len(y))]

        else:
            labels = data_labels

    else:
        labels = [None for i in range(len(y))]



    ### Data plotting (empty data are not plotted)
    for j, data in enumerate(y):

        if isinstance(data, Trace):
            if x is not None:
                if all(isinstance(ele, (UTCDateTime, datetime.datetime, str)) for ele in x):
                    xdata = [UTCDateTime(e).datetime for e in x]    # convert all to datetime.datetime objects
                    xdata = mdates.date2num(x)                      # convert all to matplotlib times (wouldn't need to, but for later it is need as ax.get_xlim() retrieves matplotlib times!)
                    ax.plot_date(x, data.data, '-', lw=lw, label=labels[j])
                    #myFmt = mdates.DateFormatter()
                    #ax.xaxis.set_major_formatter(myFmt)                # because we want dates in customised way, not matplotlib times
                else:   # normal array, relative times, matlab times, or POSIX timestamps
                    ax.plot(xdata, data, lw=lw, label=labels[j])        
            else:
                xdata = data.times(type="matplotlib")
                ax.plot_date(xdata, data.data, '-', lw=lw, label=labels[j])
                #myFmt = mdates.DateFormatter('%Y-%m-%d %H:%M:%S')
                #ax.xaxis.set_major_formatter(myFmt)

        elif isinstance(data, (Stream)):
            print(u'Using plot function of %s object and then return.' % type(data))
            print(u'No further variables passed.')
            data.plot()     # use Stream, respectively, object's plot function
            return

        else:
            if x is not None:
                if all(isinstance(ele, (UTCDateTime, datetime.datetime, str)) for ele in x):
                    xdata = [UTCDateTime(e).datetime for e in x]    # convert all to datetime.datetime objects
                    #xdata = mdates.date2num(xdata)                     # convert all to matplotlib times (wouldn't need to, but for later it is need as ax.get_xlim() retrieves matplotlib times!)
                    #myFmt = mdates.DateFormatter('%Y-%m-%d %H:%M:%S')
                    ax.plot(xdata, data, '-', lw=lw, label=labels[j])
                    #ax.xaxis.set_major_formatter(myFmt)                # because we want dates in customised way, not matplotlib times
                else:   # normal array, relative times, matlab times, or POSIX timestamps
                    ax.plot(x, data, lw=lw, label=labels[j])
            else:
                xdata = np.arange(len(data))
                ax.plot(xdata, data, lw=lw, label=labels[j])



    ### Limits of x- and y-axis
    if xlim is not None:
        if all(isinstance(ele, (UTCDateTime, datetime.datetime, str)) for ele in xlim):
            xlim = [UTCDateTime(e).datetime for e in xlim]
            xlim = [mdates.date2num(e) for e in xlim]
        ax.set_xlim( xlim )
    

    if ylim is not None:
        ax.set_ylim( ylim )
    else:           # make sure ylim is according to newly set xlim
        y_mins = []
        y_maxs = []
        for line in ax.lines:
            x = line.get_xdata()
            y = line.get_ydata()
            i = np.where( (x >= ax.get_xlim()[0]) & (x <= ax.get_xlim()[1]) )[0]          # all indexes of y_data according to xlim
            try:                                                        # e.g. if one of the `y` is empty
                line_min = y[i].min()                                   # get minimum y within all data according to xlim
                line_max = y[i].max()                                   # get maximum y within all data according to xlim
                y_mins.append(line_min)
                y_maxs.append(line_max)
            except ValueError:
                continue

        y_min = min(y_mins, default=None)
        y_max = max(y_maxs, default=None)
        if y_min is not None and y_max is not None:
            ylim = [y_min-0.025*np.abs(y_max-y_min), y_max+0.025*np.abs(y_max-y_min)] # give 2.5% margin that works both for pos. & neg. values 
            ax.set_ylim( ylim )     
        else:
            print(u'WARNING: quick_plot: specified `xlim` outside data boundaries.')
    xlim    = ax.get_xlim()
    ylim    = ax.get_ylim()

    

    ### Vertical lines for data indications
    colours = ['k', 'grey', 'lightgrey', 'red', 'green']
    if verts is not None:
                                                        # make list of lists so fo-loops work

        for k, vert in enumerate(verts):
            colour_index = k % len(colours)
                    
            if isinstance(vert, (np.ndarray, tuple, list)): 
                pass

            else:
                vert = [vert]


            for vline in vert:

                if isinstance(vline, (UTCDateTime, datetime.datetime, str)):
                    vline = UTCDateTime( vline )
                    vline = mdates.date2num(vline.datetime)
                    
                elif isinstance(vline, (complex)):
                    vline = vline.real
                    vline = xlim[0] + (xlim[1]-xlim[0])*vline   # relative to beginning in %

                else:
                    pass

                if vline>=xlim[0] and vline<=xlim[1]:
                    ax.plot( [vline, vline], ylim, lw=1, color=colours[colour_index])



    ### Horinzontal lines for data indications
    colours = ['deepskyblue','darkorange','olivedrab', 'indianred', 'orchid', 'red', 'sienna']
    if horis is not None:
                                                          # make list of lists so fo-loops work
        for k, hori in enumerate(horis):
            colour_index = k % len(colours)

            if isinstance(hori, (np.ndarray, tuple, list)):
                pass

            else:
                hori = [hori]


            for hline in hori:
                ax.plot( xlim, [hline, hline], lw=1, color=colours[colour_index])

    ### Saving & showing ..
    ax.legend(loc=legend_loc)
    if axis:
        return ax

    plt.tight_layout()

    if outfile:
        plt.savefig(outfile)

    if show:
        plt.show()

    if not keep_open:
        plt.close()
class _Arrow3D(FancyArrowPatch):
    """
    Needed for nice arrow heads in ppol plot.
    Hack found online.
    """
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)    


################  _ _ N A M E _ _ = = " _ _ M A I N _ _ "  ################
if __name__ == "__main__":
    pass