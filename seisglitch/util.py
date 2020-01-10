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
import datetime
import numpy as np


#####  matplotlib modules import  #####
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.widgets import TextBox


#####  obspy modules import  #####
from obspy import read, read_inventory, UTCDateTime
from obspy.core.stream import Stream
from obspy.core.inventory import Inventory
from obspy.signal import rotate
from obspy.clients.fdsn import Client



# Extended Obspy Stream class (adding some more functionality)
def read2(file=None):
    # wrapper to make to return Stream2 objects instead of (ObsPy's) Stream object.
    st = read(file)
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
        self.filters      = {'0' : {'freqmin':None,  'freqmax':None, 'corners':3, 'zerophase':True, 'type_taper':'hann', 'max_percentage_taper':0.03, 'num':0, 'string':'0: None'},
                             '1' : {'freqmin':1.0,   'freqmax':None, 'corners':3, 'zerophase':True, 'type_taper':'hann', 'max_percentage_taper':0.03, 'num':1, 'string':'1: >1 Hz'},
                             #'2' : {'freqmin':0.1,   'freqmax':None, 'corners':3, 'zerophase':True, 'type_taper':'hann', 'max_percentage_taper':0.03, 'num':2, 'string':'2: >0.1 Hz'},
                             #'3' : {'freqmin':None,  'freqmax':None, 'corners':3, 'zerophase':True, 'type_taper':'hann', 'max_percentage_taper':0.03, 'num':3, 'string':'3: >0.01 Hz'},
                             '2' : {'freqmin':None,  'freqmax':1,    'corners':3, 'zerophase':False,'type_taper':'hann', 'max_percentage_taper':0.03, 'num':2, 'string':'2: <1 Hz (0-phase=False)'},
                             '3' : {'freqmin':None,  'freqmax':1,    'corners':3, 'zerophase':True, 'type_taper':'hann', 'max_percentage_taper':0.03, 'num':3, 'string':'3: <1 Hz (0-phase=True)'},
                             '4' : {'freqmin':0.001, 'freqmax':None, 'corners':3, 'zerophase':True, 'type_taper':'hann', 'max_percentage_taper':0.03, 'num':4, 'string':'4: >0.001 Hz'},
                             #'5' : {'freqmin':1./100,   'freqmax':None,  'corners':3, 'zerophase':True, 'type_taper':'hann', 'max_percentage_taper':0.03, 'num':5, 'string': '5: Lowpass 100s'},
                             #'6' : {'freqmin':1./200,   'freqmax':None,  'corners':3, 'zerophase':True, 'type_taper':'hann', 'max_percentage_taper':0.03, 'num':5, 'string': '6: Lowpass 200s'},
                             #'7' : {'freqmin':1./300,   'freqmax':None,  'corners':3, 'zerophase':True, 'type_taper':'hann', 'max_percentage_taper':0.03, 'num':5, 'string':'7: Lowpass 300s'},
                             #'8' : {'freqmin':1./400,   'freqmax':None,  'corners':3, 'zerophase':True, 'type_taper':'hann', 'max_percentage_taper':0.03, 'num':5, 'string':'8: Lowpass 400s'},
                             #'9' : {'freqmin':1./500,   'freqmax':None,  'corners':3, 'zerophase':True, 'type_taper':'hann', 'max_percentage_taper':0.03, 'num':5, 'string':'9: Lowpass 500s'},
                             '5' : {'freqmin':0.1,   'freqmax':2.0,  'corners':3, 'zerophase':True, 'type_taper':'hann', 'max_percentage_taper':0.03, 'num':5, 'string':'5: RFs (0.1-2 Hz)'},
                             '6' : {'freqmin':1/9.0, 'freqmax':1.0,  'corners':3, 'zerophase':True, 'type_taper':'hann', 'max_percentage_taper':0.03, 'num':6, 'string':'6: BP 1-9 s'},
                             '7' : {'freqmin':1/6.0, 'freqmax':1.0,  'corners':3, 'zerophase':True, 'type_taper':'hann', 'max_percentage_taper':0.03, 'num':7, 'string':'7: BP 1-6 s'},
                             '8' : {'freqmin':2.0,   'freqmax':3.0,  'corners':3, 'zerophase':True, 'type_taper':'hann', 'max_percentage_taper':0.03, 'num':8, 'string':'8: 2-3 Hz'},
                             '9' : {'freqmin':2.0,   'freqmax':10,   'corners':3, 'zerophase':True, 'type_taper':'hann', 'max_percentage_taper':0.03, 'num':9, 'string':'9: BP 2-10 Hz'},
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
    def snr(self, axis=0, ddof=1):
        """ 
        Signal-to-noise ratio (SNR), as by Scipy.
        Retuns a dictionary with keys=station ID,
        values=SNR. 
        """

        SNRs = {}
        for trace in self:
            data           = np.array(trace.data)
            SNRs[trace.id] = snr(data, axis=axis, ddof=ddof)

        for traceID in SNRs.keys():
            print('ID: %s  SNR: %8.2f' % (traceID, SNRs[traceID]))

        return SNRs
    def trim(self, *args, **kwargs):

        """
        Use Obspy's trim fucntion and improve slightly.
        Improvement is that you don't have to type
        `UTCDateTime` every time you specify a time.
        """

        # Try to convert positional arguments passed to `UTCDateTime` 
        args = list(args)
        for l, arg in enumerate(args):
            try:
                args[l] = UTCDateTime(args[l])
            except Exception as err:
                #print(err)
                pass

        # Try to convert keyword arguments passed to `UTCDateTime` 
        try:
            kwargs['starttime'] = UTCDateTime(kwargs['starttime'])
        except KeyError:
            pass
        try:
            kwargs['endtime'] = UTCDateTime(kwargs['endtime'])
        except KeyError:
            pass

        super().trim(*args, **kwargs)
        self.times = self._get_times()
    def trim_common(self):

        """
        Trim traces in stream to common start- and endtime, i.e.,
        maximum starttime and minimum endtime
        """
    
        max_starttime = max([tr.stats.starttime for tr in self ])
        min_endtime   = min([tr.stats.endtime   for tr in self ])               
        self.trim(starttime=max_starttime, endtime=min_endtime) 
        self.times = self._get_times()
    def truncate(self, start_per=0, end_per=1):
    
        """
        Truncate stream by percentag with respect
        to earliest and latest trace times found in
        stream.

        Default values mean nothing is truncated.
        """
    
        start_per = float(start_per)
        end_per   = float(end_per)
    
        if start_per<0: 
            start_per = 0

        if end_per>1: 
            end_per = 1

        if start_per==0 and end_per==1:
            return
    
        mint       = min([tr.stats.starttime for tr in self])
        maxt       = max([tr.stats.endtime for tr in self])
        starttime  = mint + (maxt-mint)*start_per
        endtime    = maxt - (maxt-mint)*(1-end_per)
        self.trim(starttime=starttime, endtime=endtime)
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
    def gain_correction(self): 
        if self.inventory:
            print()
            print(u'GAIN CORRECTION APPLIED:')
            try:

                if not self.gain_removed:
                    self.gain_removed = True
                    for trace in self:
                        response   = self.inventory.get_response(trace.id, trace.stats.starttime)
                        gain       = response._get_overall_sensitivity_and_gain()[1]                
                        trace.data = trace.data / gain
                        print(u'  %15s : overall sensitivity and gain (division) %s' % (trace.id, gain))
                else:
                    self.gain_removed = False
                    for trace in self:
                        response   = self.inventory.get_response(trace.id, trace.stats.starttime)
                        gain       = response._get_overall_sensitivity_and_gain()[1]                
                        trace.data = trace.data * gain
                        print(u'  %15s : overall sensitivity and gain (multiplication) %s' % (trace.id, gain))

            except Exception as err:
                print(u'WARNING:  %s' % err)

        else:
            print()
            print(u'No matching response found for gain correction. Nothing done.')
    def filtering(self, filter):
    
        """ 
        Filter ObsPy's stream object has to the given filter specifications. 
        If only freqmax or freqmin is given, highpass or lowpass filter is used, respectively.
        By default the Butterworth filter implementation is used. For a different filter,
        change the code!
    
        See details:
          https://docs.obspy.org/packages/autogen/obspy.core.stream.Stream.filter.html
        """
    
        if isinstance(filter, (str, int)):
            filter = self.filters[str(filter)]


        #filter = {**taper, **filter}

        if filter['freqmin'] or filter['freqmax']:
            self.detrend('demean')
            self.detrend('linear')
            self.taper(type=filter['type_taper'], max_percentage=filter['max_percentage_taper'])
            if not filter['freqmin']:
                self.filter('lowpass', freq=filter['freqmax'], corners=filter['corners'], zerophase=filter['zerophase'])
            elif not filter['freqmax']:
                self.filter('highpass', freq=filter['freqmin'], corners=filter['corners'], zerophase=filter['zerophase'])
            else:
                self.filter('bandpass', freqmin=filter['freqmin'], freqmax=filter['freqmax'], corners=filter['corners'], zerophase=filter['zerophase'])

            self.current_filter_num = filter['num']
            self.current_filter_str = filter['string']
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
            components = sorted( list( set( [tr.stats.channel[-1] for tr in self] )), reverse=True)
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
    def _get_inventory(self, file=None, online='IPGP'):

        """
        Gets latest `file` (with respect to download) and reads it as an Obspy inventory object.
        From this inventory, only the part is extracted that matches the start and end times of the stream `self`,
        as well as matches all networks, stations, locations and channels of the stream.

        The inventory file is then assigned to `self.inventory` and also returned in case
        further processing is needed.


        https://docs.obspy.org/packages/obspy.clients.fdsn.html
        """


        inv = Inventory()

        # retrieve from online
        if not file:

            client = Client(online)

            for trace in self:
                network, station, location, channel = trace.id.split('.')
                inv += client.get_stations(network = network, 
                                        station    = station, 
                                        location   = location, 
                                        channel    = channel,
                                        starttime  = self.times[0], 
                                        endtime    = self.times[1],
                                        level      = 'response')


        # retrieve from file
        else:

            inv_file  = max( glob.glob( file ), key=os.path.getctime)   # most recent one, in case there multiple
            inventory = read_inventory(inv_file)

            for trace in self:
                network, station, location, channel = trace.id.split('.')
                inv += inventory.select(network   = network, 
                                        station   = station, 
                                        location  = location, 
                                        channel   = channel, 
                                        starttime = self.times[0], 
                                        endtime   = self.times[1])

        return inv
    def _set_inventory(self, inventory=None, file=None, online='IPGP'):

        if not inventory:
            inventory = self._get_inventory(file=file, online=online)

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
    def plot_polarization(self, title='', return_results=False, show=True):
    
        """
        Plot polarization of data.

        TO BE REDONE!
        """

        ### get data ready
        stream_N = self.select(component='N') or self.select(component='Q') or self.select(component='U') or self.select(component='1') or self.select(component='T')
        stream_E = self.select(component='E') or self.select(component='T') or self.select(component='V') or self.select(component='2') or self.select(component='R')
        stream_Z = self.select(component='Z') or self.select(component='L') or self.select(component='W') or self.select(component='3')

        if not (stream_N and stream_E):
            print('')
            print('Return. No idea how to perform polarization on components: %s' % ', '.join([i.stats.channel for i in self]) )
            return
    def plot2(self, store_dir=os.getcwd(), store_name='*', verticals=(), type='normal', method='full', save_and_no_show=False, xlim=[], ylim=[], **kwargs):

        """
        Enhanced plot of stream as 'type', using ObsPy's plot function.
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
                        print('remove', remove)
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

            def instrument_removal(output, pre_filt=None, water_level=60): 
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
                    corr.plot(xlim=xlim)
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
                    corr.filtering(self.current_filter_num)
                    corr.plot(xlim=xlim)
                    return


                if not mouse_clicks.keys():
                    print()
                    print('Please select at least one station location!')

                else:
                    if not self.inventory:
                        print('Reading inventory file for instrument response removal ..')
                        self._set_inventory()
                    

                    components = self._get_components()
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
                    corr.plot(xlim=xlim)


            if evt.key.lower()=='a':                            # Acceleration as target unit after instrument response removal 
                pre_filt = (0.001, 0.002, 50, 60)
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


            elif evt.key.lower()=='d':                          # Displacement as target uni after instrument response removal 
                pre_filt = (0.001, 0.002, 50, 60)
                instrument_removal('DISP')


            elif evt.key.lower()=='g':                          # Gain removal / add 
                self.gain_correction()
                self.plot(xlim=axis.get_xlim())


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
                    if type == 'normal':
                        xlim   = [ UTCDateTime(mdates.num2date(x)) for x in xlim ]
                    elif type == 'relative':
                        xlim =  [ min_starttime+x for x in xlim ]
                    
                    polar = self.copy()
                    polar = polar.select(id=stat_id+'*')
                    polar.trim(starttime=xlim[0], endtime=xlim[1])
                    polar.plot_polarization(title='%s, %s, %s, %s-%s'% (stat_id, self.unit, self._get_filter_str(), xlim[0].strftime('%H:%M:%S'), xlim[1].strftime('%H:%M:%S')), show=show)


            elif evt.key=='P':                                  # Polarisation plot as batch job 
                for l, stat_id in enumerate( mouse_clicks.keys() ):
                    if mouse_clicks[stat_id] == 1:
                        continue
                    elif mouse_clicks[stat_id] != 2:
                        print('Choose at least 2 time limits (left-click) per station to perform action.')
                        continue
                    print()
                    print('Running P-pol batch job ..')

                    if l==len(mouse_clicks.keys())-1:
                        show=True
                    else:
                        show=False
                    indexes       = [i for i, ax in enumerate(axes) if stat_id in [c for c in ax.get_children() if isinstance(c, mpl.text.Text)][0].get_text() ]
                    ax_1st        = axes[indexes[0]]
                    xlim          = sorted( [ax_1st.lines[-1].get_xdata()[0], ax_1st.lines[-2].get_xdata()[0]] )
                    if type == 'normal':
                        xlim   = [ UTCDateTime(mdates.num2date(x)) for x in xlim ]
                    elif type == 'relative':
                        xlim =  [ min_starttime+x for x in xlim ]
                    
                    # create windows for batch job
                    xlims    = []
                    span     = xlim[1]-xlim[0]
                    pro_side = 2
                    fraction = 0.15
                    for l in np.arange(-pro_side, pro_side+1):
                        for i in np.arange(-pro_side, pro_side+1):
                            window = [ xlim[0]+l*span*fraction, xlim[1]+i*span*fraction]
                            xlims.append( window )

                    # iterate over windows and filters to calculate polarisations
                    results = []
                    unit    = self.unit
                    for xlim in xlims:
                        for fil_num in self.filters.keys():
                            polar   = self.original or self
                            polar   = polar.copy()
                            polar   = polar.select(id=stat_id+'*')
                            polar.filtering(fil_num)
                            polar.trim(starttime=xlim[0], endtime=xlim[1])
                            
                            fil_str = polar._get_filter_str()
                            window  = '%s-%s' % (xlim[0].strftime('%H:%M:%S'), xlim[1].strftime('%H:%M:%S'))
                            title   = '%s, %s, %s, %s' % (stat_id, unit, fil_str, window)                       
                            result  = polar.plot_polarization(title=title, return_results=True)
                            results.append( list(result)+xlim+[fil_num] )

                    # retrieve needed variables for best polarisation
                    results      = np.array(results)
                    ind_high     = np.argmax(results[:,11])     # highest Rect_H
                    rect_Hs      = list( results[:,11] )
                    fil_num_high = results[ind_high,-1] 
                    xlim_high    = results[ind_high,-3:-1] 

                    # plot best polarisationg
                    polar_best   = self.copy()
                    polar_best   = polar_best.select(id=stat_id+'*')
                    polar_best.filtering(fil_num_high)
                    polar_best.trim(starttime=xlim_high[0], endtime=xlim_high[1])

                    fil_str      = polar_best._get_filter_str()
                    window       = '%s-%s' % (xlim_high[0].strftime('%H:%M:%S'), xlim_high[1].strftime('%H:%M:%S'))
                    title        = '%s, %s, %s, %s' % (stat_id, unit, fil_str, window)

                    print('')
                    print('B E S T')
                    print('-' * 40)
                    polar_best.plot_polarization(title='Best: %s' % title, show=False)
                    print('-' * 40)

                    # plot histogram
                    hist_title = 'Histogram P-pol batch job (%s)' % self.unit
                    if title in plt.get_figlabels():
                        fig = plt.figure(hist_title)
                        fig.clf()
                    fig = plt.figure(num=hist_title)
                    fig.suptitle('Histogram of %s P-wave polarisation calculations' % (len(rect_Hs)), fontsize=9)
                    ax = fig.add_subplot(111)
                    ax.grid(ls='-.', lw=0.5, zorder=0)
                    ax.set(xlabel='Horizontal rectilinearity', ylabel='Number of occurences')
                    ax.hist(rect_Hs, bins=100, range=(0,1), rwidth=1, ec='k', lw=0.5, zorder=3)
                    if show:
                        plt.show()


            elif evt.key.lower()=='q':                          # Quit current figures
                # close
                print('Closed current figure window')
                plt.gcf()
                plt.close()


            elif evt.key.lower()=='r':                          # Rotate data 
                
                self.merge(method=1)
                components = self._get_components()


                ## rotate data to ZNE
                if re.search(r'[U|V|W]', components):

                    # correct gain
                    if not self.gain_removed:
                        self.gain_correction()

                        if self.original:
                            self.original.gain_correction()
                            self.original.gain_removed = True

                        if self.removed:
                            self.removed.gain_correction()
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
                        rotate2VBBUVW(self.original, inventory=self.inventory)

                    if self.removed:
                        rotate2VBBUVW(self.removed, inventory=self.inventory)

                self.plot(xlim=axis.get_xlim())


            elif evt.key=='t':                                  # Trim file and save original data in viewed x-lim bounds 
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

                    
                    if type == 'normal':
                        xlim   = [ UTCDateTime(mdates.num2date(x)) for x in xlim ]
                    elif type == 'relative':
                        xlim =  [ min_starttime+x for x in xlim ]

                    def submit(filename):
                        if self.original:
                            trim = self.original.copy()
                        else:
                            trim = self.copy()
                        trim = trim.select(id=stat_id+'*')
                        trim.trim(starttime=xlim[0], endtime=xlim[1])
                        trim.write(filename)
                        plt.close()
                    
                    fig      = plt.figure(figsize=(17,0.7), num='Save trimmed file')
                    ax       = fig.add_subplot(111)
                    
                    propose  = os.path.join( os.path.dirname(self.origin), 'TRIMORIG_'+os.path.basename(self.origin) )
                    text_box = TextBox(ax, 'Filename:', initial=propose)
                    text_box.on_submit(submit)
                    plt.show()


            elif evt.key=='T':                                  # Trim file and save present data in viewed x-lim bbounds 
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

                    
                    if type == 'normal':
                        xlim   = [ UTCDateTime(mdates.num2date(x)) for x in xlim ]
                    elif type == 'relative':
                        xlim =  [ min_starttime+x for x in xlim ]

                    def submit(filename):
                        trim = self.copy()
                        trim = trim.select(id=stat_id+'*')
                        trim.trim(starttime=xlim[0], endtime=xlim[1])
                        trim.write(filename)
                        plt.close()
                    
                    fig      = plt.figure(figsize=(17,0.7), num='Save trimmed file')
                    ax       = fig.add_subplot(111)
                    
                    propose  = os.path.join( os.path.dirname(self.origin), 'TRIMPRES_'+os.path.basename(self.origin) )
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


            elif evt.key.lower()=='v':                          # Velocity as target uni after instrument response removal 
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
                    filt              = self.original.copy()
                    filt.original     = self.original.copy()
                
                filt.removed = self.removed
                filt.filtering(evt.key)
                filt.plot(xlim=xlim)    
    

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

        fig = self.obs_plot(type=type, show=False, handle=True, equal_scale=False, method='full', **kwargs)
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
                axes[i].text(0.99, 0.95, "Close all windows: 'w'.", horizontalalignment='right', color='red', transform=axes[i].transAxes, fontsize=7)
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
        if save_and_no_show:
            flags['saved'] = True
            plt.savefig(outfile, dpi=200, bbox_inches='tight')
            plt.close(fig)
        else:
            # show plot
            plt.show() 

# Glitch related
def merge_glitch_files(outfile, *glitch_files, starttime_sort=True):

    """
    Merging together glitch extration files produced by `glitch_detector`.
    This helps combining into one file representing a certain period, e.g.
    """


    ### OUTPUT
    print(u'Merging the following files:')
    for file in glitch_files:
        print(file)


    ### RETRIEVE HEADER AND ALL DATA
    data = []
    with open(outfile, 'w') as fp: 
        counter = 0

        for file in glitch_files: 

            with open(file, 'r') as fp2:
                lines  = fp2.readlines() 
                header = lines[0:3]
                data  += [line for line in lines[2:] if line and not line.startswith('#')]
                    

    ### SORTING
    if starttime_sort:
        data_to_sort = np.array( [line.split() for line in data] )
        sort_indices = data_to_sort[:,1].argsort()
        data         = np.array(data)
        data         = data[sort_indices]


    ### WRITING OUT
    counter = 0
    with open(outfile, 'w') as fp: 
                    
        for header_line in header:
            fp.write(header_line)

        for line in data:
            counter += 1 
            number   = '%06d' % counter 
            fp.write(number+line[6:])


    ### OUTPUT
    print()
    print(u'Merged glitch file to:')
    print(outfile)

# Convenient helpers
def moving_window(*data, window_length_in_samples=100, step_in_samples=50, equal_end=True):

    """
    Yield data of moving windows according to given parameters.
    """

    i = 0
    while True:

        # window indices
        index_window_start = i * step_in_samples
        index_window_end   = index_window_start + window_length_in_samples

        # latest when start time of window is >= number samples, then no more windows
        i += 1
        if index_window_start >= len(data):
            break

        # insist last window is of length `window_length_in_samples` 
        if equal_end:
            if len(data[index_window_start:index_window_end]) < window_length_in_samples:
                break


        yield_data = (array[index_window_start:index_window_end] for array in data)

        yield (index_window_start, index_window_end, *yield_data)

# Mathematical
def rotate_2D(comp_1, comp_2, angle, clockwise=False):

    """
    This function rotates 2 data traces (i.e. numpy arrays, e.g., seismological data) 
    into the desired coordinate system using 'angle' (clockwise direction by default).

    The given algorithm works for 'comp_2' being oriented 90° clockwise to 'comp_1'.

    ^  
    | comp_1
    |
    |
    |
    0------> comp_2


    :type comp_1:  np.array
    :param comp_1: List with floats representing data.
                   Note 'comp_1' is oriented 90° counter-clockwise to 'comp_2'.

    :type comp_2:  np.array
    :param comp_2: List with floats representing data.
                   Note 'comp_2' is oriented 90° clockwise to 'comp_1'.

    :type angle:  float
    :param angle: Angle about which both data traces are rotated. Can be negative.

    :type clockwise:  bool
    :param clockwise: if True  :         clockwise rotation of both data traces 
                      if False : counter-clockwise rotation of both data traces 

    :type return:  np.array, np.array
    :param return: rotated data lists
    """

    if not clockwise:
        angle = 360-angle               # convert angle to as it would be clockwise rotation


    comp_1_new =  comp_1*np.cos(angle*np.pi/180) + comp_2*np.sin(angle*np.pi/180)
    comp_2_new = -comp_1*np.sin(angle*np.pi/180) + comp_2*np.cos(angle*np.pi/180)

    return comp_1_new, comp_2_new
def normalise(data, scale_to_between=[]):

    """
    Normalise passed data (array-like).
    `scale_to_between` is list with 2 elements.

    Returns data as was if length of data is one or two.
    """

    data = np.asarray( data )
    
    if len(data)==0 or len(data)==1:
        return data

    if isinstance(scale_to_between, (int, float)):
        scale_to_between = [scale_to_between]

    if scale_to_between:
        if len(scale_to_between) == 1:
            scale_to_between = [0, scale_to_between[0]]
        scale_to_between.sort()

        scale  = abs(scale_to_between[-1]-scale_to_between[0]) / 2.
        drange = max(data)-min(data)
        data   = data * 2 / drange
        data   = data - max(data)+1             # data have y values between [-1,1]
        data  *= scale                          # data have y values between [-1/scale,1/scale] 
        data  += scale_to_between[-1]-scale     # eacht trace has y values filling range `scale_to_between`
    
    else:
        data /= max(abs(data))

    return data
def snr(data, axis=0, ddof=1):

    """
    Signal-to-noise ratio (SNR), as by Scipy.
    """

    data = np.array(data)
    mean = data.mean(axis)
    sd_d = data.std(axis=axis, ddof=ddof)
    return np.where(sd_d==0, 0, mean/sd_d)

# InSight time conversions / ways of displaying
def solify(UTC_time, sol0=UTCDateTime(2018, 11, 26, 5, 10, 50.33508)):

    """
    :copyright:
        Simon Stähler (mail@simonstaehler.com), 2018
        Martin van Driel (Martin@vanDriel.de), 2018
    :license:
        None
    """

    SEC_PER_DAY_EARTH = 86400
    SEC_PER_DAY_MARS  = 88775.2440 #before: 88775.244147    

    MIT = (UTC_time - sol0) / SEC_PER_DAY_MARS
    t   = UTCDateTime((MIT - 1) * SEC_PER_DAY_EARTH)

    return t
def UTCify(LMST_time, sol0=UTCDateTime(2018, 11, 26, 5, 10, 50.33508)):
    """
    :copyright:
        Simon Stähler (mail@simonstaehler.com), 2018
        Martin van Driel (Martin@vanDriel.de), 2018
    :license:
        None
    """
    SEC_PER_DAY_EARTH = 86400
    SEC_PER_DAY_MARS  = 88775.2440 #before: 88775.244147

    MIT      = float(LMST_time) / SEC_PER_DAY_EARTH + 1
    UTC_time = UTCDateTime(MIT * SEC_PER_DAY_MARS + float(sol0))

    return UTC_time
def sol2UTC(sol):
    # Convert a float, interpreted as InSight sol, to UTC.
    return UTCify(UTCDateTime('1970-01-01T00:00:00.000000Z')+datetime.timedelta(days=sol-1))
def ptime(time=None):

    """
    PRINT TIME

    Small script to display current UTC time in most common ways.
    Pass a terrestrial time if you pass one.
    """

    if time:
        print(u'TIME GIVEN')
        print(u'----------')
        if isinstance(time, float):                                     # MATLAB
            time = mdates.num2date(time)
        elif isinstance(time, int):                                     # sol
            time = sol2UTC(time).datetime
        elif isinstance(time, (UTCDateTime,datetime.datetime,str)):     # typical obspy str
            time = UTCDateTime(time).datetime
        else:                                                           # others
            print(u'Passed time variable has no valid format (float, str, datetime.datetime, or UTCDateTime.)')
    else:
        print(u'TIME UTC NOW')
        print(u'------------')
        time = datetime.datetime.utcnow()

    print(u'MatLab:         %s'                % mdates.date2num(time))
    print(u'DateTime:       %s'                % time.__repr__())
    print(u"UTCDateTime:    UTCDateTime('%s')" % UTCDateTime(time))
    print(u'LocalMeanSolar: %sM%s'             % (solify(UTCDateTime(time)).julday, solify(UTCDateTime(time)).strftime('%H:%M:%S')))
def ltime(hms, sols_range=[], is_UTC=False):

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
    if not sols_range or len(sols_range)==1:
        UTC_now    = UTCDateTime( time.time() )
        LMST_now   = solify(UTC_now)

        if len(sols_range)==1:
            sols_range = [LMST_now.julday-sols_range[0], LMST_now.julday-1]
        else:
            sols_range = [LMST_now.julday-10,            LMST_now.julday-1]


    ## PRINTS
    print('UTC                    LMST')
    print('---                    ----')    
    for sol in range(sols_range[0], sols_range[1]+1):

        if is_UTC:
            time_str_ymd = sol2UTC(sol).strftime('%Y-%m-%d')
            time_str_HMS = '%s:%s:%s' % (hms[0:2], hms[2:4], hms[4:6])

            UTC_time     = UTCDateTime( time_str_ymd + 'T' + time_str_HMS)
            LMST_time    = solify(UTC_time)


        else:
            LMST_time    = UTCDateTime('1970-01-01T%s:%s:%s.000000Z' % (hms[:2], hms[2:4], hms[4:])) + datetime.timedelta(days=sol)
            UTC_time     = UTCify(LMST_time)


        print('%s    %sS%s' % (UTC_time.strftime('%Y-%m-%dT%H:%M:%S'), LMST_time.julday, LMST_time.strftime('%H:%M:%S')))


################  _ _ N A M E _ _ = = " _ _ M A I N _ _ "  ################
if __name__ == "__main__":
    pass