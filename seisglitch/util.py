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
import io
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
import urllib
import numpy as np
import scipy



#####  matplotlib modules import  #####
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.switch_backend('TKAgg')
import matplotlib.dates as mdates
from matplotlib.patches import Arc, FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
from matplotlib.widgets import TextBox


#####  obspy modules import  #####
import obspy
from obspy import read, read_inventory, UTCDateTime
from obspy.core.stream import Stream, Trace
from obspy.core.inventory import Inventory
from obspy.signal import rotate
from obspy.clients.fdsn import Client
from obspy.clients.fdsn.client import CustomRedirectHandler, NoRedirectionHandler, build_url
from obspy.signal.filter import envelope
from obspy.io.mseed.core import ObsPyMSEEDFilesizeTooLargeError


##### seisglitch modules import #####
from seisglitch.math import normalise, covariance_matrix, rotate_2D, unit_vector


# Extended Obspy Stream class (adding some more functionality)
def read2(file=None, reclen=512, chunksize=1000000, **kwargs):
    # wrapper to return Stream2 objects instead of ObsPy's Stream object.
    try:
        st = read(file, **kwargs)
    except ObsPyMSEEDFilesizeTooLargeError:
        st        = Stream()
        reclen    = reclen
        chunksize = chunksize * reclen        # around 500 MB by default values
        with io.open(file, "rb") as fh:
            while True:
                with io.BytesIO() as buf:
                    c = fh.read(chunksize)
                    if not c:
                        break
                    buf.write(c)
                    buf.seek(0, 0)
                    stream = read(buf, **kwargs)
                    st += stream

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
        self.filters      = {'0' : {'type' : False,       'options' : {                                                               },   'string':'no filter'},
                             '1' : {'type' : 'highpass',  'options' : {'freq':1.0,                       'corners':3, 'zerophase':True},   'string':'1: >1 Hz'},
                            #'2' : {'type' : 'highpass',  'options' : {'freq':0.1,                       'corners':3, 'zerophase':True},   'string':'2: >0.1 Hz'},
                            #'3' : {'type' : False,       'options' : {                                                               },   'string':'3: >0.01 Hz'},
                             '2' : {'type' : 'bandpass',  'options' : {'freqmin':0.01,  'freqmax':0.1,   'corners':1, 'zerophase':False},   'string':'2: 0.01s < f < 0.1'},
                             '3' : {'type' : 'lowpass',   'options' : {'freq':0.5,                       'corners':3, 'zerophase':True},   'string':'3: <0.5Hz'},
                             '4' : {'type' : 'bandpass',  'options' : {'freqmin':0.001,  'freqmax':0.1,  'corners':1, 'zerophase':True},   'string':'4: 0.001s < f < 0.1'},
                            #'5' : {'type' : 'highpass',  'options' : {'freq':0.001,                     'corners':3, 'zerophase':True},   'string':'5: 0.001s < f'},
                            #'6' : {'type' : 'highpass',  'options' : {'freq':1./200,                    'corners':3, 'zerophase':True},   'string': '6: Lowpass 200s'},
                            #'7' : {'type' : 'highpass',  'options' : {'freq':1./300,                    'corners':3, 'zerophase':True},   'string':'7: Lowpass 300s'},
                            #'8' : {'type' : 'highpass',  'options' : {'freq':1./400,                    'corners':3, 'zerophase':True},   'string':'8: Lowpass 400s'},
                            #'9' : {'type' : 'highpass',  'options' : {'freq':1./500,                    'corners':3, 'zerophase':True},   'string':'9: Lowpass 500s'},
                            #'5' : {'type' : 'bandpass',  'options' : {'freqmin':0.1,    'freqmax':2.0,  'corners':3, 'zerophase':True},   'string':'5: RFs (0.1-2 Hz)'},
                             '5' : {'type' : 'bandpass',  'options' : {'freqmin':1/6.0,  'freqmax':1.0,  'corners':3, 'zerophase':True},   'string':'5: BP 1-6 s'},
                             '6' : {'type' : 'bandpass',  'options' : {'freqmin':1/8.0,  'freqmax':1.0,  'corners':3, 'zerophase':True},   'string':'6: BP 1-8 s'},
                             '7' : {'type' : 'highpass',  'options' : {'freq':5,                         'corners':3, 'zerophase':True},   'string':'7: >5 Hz'},
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
    def gain_correction(self, verbose=False, channel='?[LMH]?'): 
        if self.inventory:

            if verbose:
                print()
                print(u'GAIN CORRECTION APPLIED:')


            try:

                if not self.gain_removed:
                    for trace in self.select(channel=channel):
                        response   = self.inventory.get_response(trace.id, trace.stats.starttime)
                        gain       = response._get_overall_sensitivity_and_gain()[1]                
                        trace.data = trace.data / gain

                        if verbose:
                            print(u'  %15s : overall sensitivity and gain (division) %s' % (trace.id, gain))
                    self.gain_removed = True

                else:
                    for trace in self.select(channel=channel):
                        response   = self.inventory.get_response(trace.id, trace.stats.starttime)
                        gain       = response._get_overall_sensitivity_and_gain()[1]                
                        trace.data = trace.data * gain

                        if verbose:
                            print(u'  %15s : overall sensitivity and gain (multiplication) %s' % (trace.id, gain))
                    self.gain_removed = False
                
                return 0

            except Exception as err:
                print(u'WARNING: %s No gain correction done.' % err)
                return -1

        else:
            print(u'INFO: No inventory found. No gain correction done.')
            return 0
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
    def splitwrite(self, num=2, folder=None, file_format='mseed'):

        """
        When having a long file, this method can help.
        It splits the stream object in `num` streams that have no overlapping times
        and writes out each new stream out. The new start and end times are
        aquidistant, however, actual data contained may be distributed unevenly due
        to gaps / overlaps.
        """

        print(u'Written file(s):')

        for i in range(num):

            st=self.copy()
            st.trim2(i/num, (1+i)/num)

            if not folder: 
                if self.origin:
                    folder = os.path.dirname( self.origin )
                else:
                    folder = os.getcwd()
            name = '%s_%s_raw.%s' % (st.times[0].strftime('%Y%m%dT%H%M'),st.times[1].strftime('%Y%m%dT%H%M'),file_format)
            outfile = os.path.join(folder, name)

            st.write(outfile)
            print(outfile)
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
            except Exception:
                pass

        return filt
    def _get_components(self):
        components = []
        try:
            components = [tr.stats.channel[-1] for tr in self]
        except AttributeError:
            pass

        return components
    def _get_channels(self):
        channels = []
        try:
            channels = [tr.stats.channel for tr in self]
        except AttributeError:
            pass

        return channels
    def _get_sampling_rates(self):
        sampling_rates = []
        try:
            sampling_rates = [tr.stats.sampling_rate for tr in self]
        except AttributeError:
            pass

        return sampling_rates
    def _get_ids(self):
        ids = []
        try:
            ids = [tr.id for tr in self]
        except AttributeError:
            pass

        return ids
    def _get_times(self):
        try:
            times = min([tr.stats.starttime for tr in self]), max([tr.stats.endtime for tr in self])
        except ValueError:  # empty stream object
            times = None, None

        return times
    def _get_inventory(self, source, verbose=True):

        """
        Gets latest `inventory` and reads it as an Obspy inventory object.
        From this inventory, only the part is extracted that matches the start and end times of the stream `self`,
        as well as matches all networks, stations, locations and channels of the stream.

        The inventory is then returned.

        https://docs.obspy.org/packages/obspy.clients.fdsn.html
        """


        inv = Inventory()

        if isinstance(source, Inventory): 
            for trace in self:
                network, station, location, channel = trace.id.split('.')
                inv += source.select(network   = network, 
                                     station   = station, 
                                     location  = location, 
                                     channel   = channel, 
                                     starttime = trace.stats.starttime, 
                                     endtime   = trace.stats.endtime)
            if verbose:
                print(u'INFO: Inventory used that was passed')


        elif os.path.isfile(source):
            inventory = read_inventory(source)
            for trace in self:
                network, station, location, channel = trace.id.split('.')
                inv += inventory.select(network   = network, 
                                        station   = station, 
                                        location  = location, 
                                        channel   = channel, 
                                        starttime = trace.stats.starttime, 
                                        endtime   = trace.stats.endtime)
            if verbose:            
                print(u"INFO: Inventory read from file '%s'" % source)


        else:
            client = Client(source)
            for trace in self:
                network, station, location, channel = trace.id.split('.')
                inv += client.get_stations(network   = network, 
                                           station   = station, 
                                           location  = location, 
                                           channel   = channel,
                                           starttime = trace.stats.starttime, 
                                           endtime   = trace.stats.endtime,
                                           level     = 'response')                
            if verbose:
                print(u"INFO: Inventory retrieved online from '%s'" % source)

        return inv
    def set_inventory(self, source='IPGP', verbose=True):

        """
        
        """

        inventory      = self._get_inventory(source=source, verbose=verbose)
        self.inventory = inventory
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
    def envelope_smooth(self, envelope_window_in_sec=10, mode='same'):

        for tr in self:

            tr.data = envelope(tr.data)
            w       = np.ones(int(envelope_window_in_sec * tr.stats.sampling_rate))
            w      /= w.sum()
            tr.data = np.convolve(tr.data, w, mode=mode)
    def plot2(self, store_dir=os.getcwd(), store_name='*', verticals=[], time='normal', method='full', xlim=[], ylim=[], store_plot=False, show=True, **kwargs):


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
                        for line in lines[::-1][:mouse_clicks[stat_id]]:
                            line.remove()
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
                elif evt.button == 2:                                                       # middle-click (set ylim for each axis to maximum available ylim of all axes whilst centering data of each axis == equal scale)
                    xlim  = axis.get_xlim()
                    ylims = []

                    for axis in axes:

                        y_mins = []
                        y_maxs = []     

                        for line in axis.lines:
                            x = line.get_xdata()
                            if x[0]==x[-1]:                                                 # exclude vertical lines
                                continue
                            y = line.get_ydata()
                            i = np.where( (x >= xlim[0]) & (x <= xlim[1]) )[0]              # all indexes of y_data according to xlim
                            if not i.any():
                                continue                                                    # e.g. data gaps, that is, 2 lines within axes where one does not lie within chose xlims

                            line_min = y[i].min()                                           # get minimum y within all data according to xlim
                            line_max = y[i].max()                                           # get maximum y within all data according to xlim
                            y_mins.append(line_min)
                            y_maxs.append(line_max)
                        
                        y_min = min(y_mins)
                        y_max = max(y_maxs)
                        ylims.append([y_min, y_max])

                    ylims        = np.array(ylims)
                    ylims_vals   = [ylim[1]-ylim[0] for ylim in ylims]
                    ylim_max_ind = np.argmax(ylims_vals)

                    for r, ylim in enumerate(ylims):
                        if r == ylim_max_ind:
                            pass
                        else:
                            ylim = [ylim[0]-ylims_vals[ylim_max_ind]/2+(ylim[1]-ylim[0])/2, ylim[1]+ylims_vals[ylim_max_ind]/2-(ylim[1]-ylim[0])/2]
                        ylim_set = [ylim[0]-0.025*np.abs(ylim[1]-ylim[0]), ylim[1]+0.025*np.abs(ylim[1]-ylim[0])] # give 2.5% margin that works both for pos. & neg. values
                        axes[r].set_ylim( ylim_set ) 
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
                        self.set_inventory()
                    

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
                        xlim = [ UTCDateTime(mdates.num2date(x)) for x in xlim ]
                    elif time == 'relative':
                        xlim =  [ min_starttime+x for x in xlim ]
                    polar = self.copy()
                    polar = polar.select(id=stat_id+'*')

                    title = '%s, %s, %s, %s-%s'% (stat_id, self.unit, self._get_filter_str(), xlim[0].strftime('%H:%M:%S'), xlim[1].strftime('%H:%M:%S'))
                    measurement = ppol(stream=polar, starttime=xlim[0], endtime=xlim[1], **kwargs)
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
                    self.set_inventory()

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
                filt.plot2(xlim=xlim, verticals=verticals, time=time, method=method, show=True, **kwargs)    
    
        # variables
        mouse_clicks  = {}
        flags         = {'saved' : False}
        colours       = 'bgcm'
        outfile       = os.path.join(store_dir, store_name+'.png')  # change format if desired (used if plot shall be saved)
        title         = '%s, %s' % (self.unit, self._get_filter_str())
        min_starttime = min( np.array( [tr.stats.starttime for tr in self] ), default=np.nan)
    
    
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
                for vert in verticals:
                    axes[i].plot([vert.datetime,vert.datetime], axes[i].get_ylim(), color='red', lw=1.0)
    
            if i == 0:
                axes[i].text(0.99, 0.92, self._get_filter_str(), horizontalalignment='right', color='red', transform=axes[i].transAxes, fontsize=8, bbox=dict(boxstyle='round,pad=0.2', fc='white', ec="red", lw=0.5))
    
            if i+1 == len(self):
                axes[i].text(0.99, 0.93, "Close all windows: 'w'.", horizontalalignment='right', color='red', transform=axes[i].transAxes, fontsize=7)


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
# FANCY DETRENDING ?

# InSight related
class marstime():

    def __init__(self, UTC=None, LMST=None, sec_per_day_mars=88775.2440, sec_per_day_earth=86400, sol0=UTCDateTime(2018, 11, 26, 5, 10, 50.33508), strftime=('%Y-%m-%dT%H:%M:%S')):

        self.sec_per_day_mars  = sec_per_day_mars
        self.sec_per_day_earth = sec_per_day_earth

        self.UTC_time          = None
        self.UTC_string        = None
        self.julday            = None

        self.LMST_time         = None
        self.LMST_string       = None
        self.sol0              = sol0

        if UTC is None and LMST is None:
            self.UTC_time  = UTCDateTime( datetime.datetime.utcnow() )       # now
            self._LMSTify(strftime)

        if UTC is not None:
            self.UTC_time = UTCDateTime(UTC)
            self._LMSTify(strftime)

        if LMST is not None:
            try:
                self.LMST_time = UTCDateTime(LMST)
            except:
                sol            = int(LMST.split('M')[0])
                hms            = '000000.000000'
                if len(LMST.split('M'))>1:
                    passed_hms = LMST.split('M')[1].replace(':','')
                    hms        = passed_hms + hms[len(passed_hms):]
                self.LMST_time = UTCDateTime('1969-12-31T00:00:00.000000Z')     + \
                                    datetime.timedelta(days = sol)              + \
                                    datetime.timedelta(hours = int(hms[:2]))    + \
                                    datetime.timedelta(minutes = int(hms[2:4])) + \
                                    datetime.timedelta(seconds = float(hms[4:]))
            self._UTCify(strftime)    
    def __str__(self):

        """
        Display current UTC time in most common ways.
        """

        output = []

        output.append(u'TIMES\n')
        output.append(u'-----\n')
        output.append(u'UTC string:     %s\n' % self.UTC_string.__repr__())
        output.append(u'DateTime UTC:   %s\n' % self.UTC_time.datetime.__repr__())
        output.append(u'UTCDateTime:    %s\n' % self.UTC_time.__repr__())
        output.append(u'MatLab UTC:     %s\n' % mdates.date2num(self.UTC_time.datetime))
        output.append(u'LocalMeanSolar: %s'   % self.LMST_string.__repr__())

        string = ''.join( output )
        return string
    def _UTCify(self, strftime):

        """
        :copyright:
            Simon Stähler (mail@simonstaehler.com), 2018
            Martin van Driel (Martin@vanDriel.de), 2018
        :license:
            None
        """

        MIT              = float(self.LMST_time) / self.sec_per_day_earth + 1
        UTC_time         = UTCDateTime(MIT  * self.sec_per_day_mars  + float(self.sol0))
        
        self.LMST_time   = self.LMST_time        
        self.LMST_string = self._LMST_string()
        self.sol         = int(self.LMST_string.split('M')[0])
        
        self.UTC_time    = UTC_time
        self.UTC_string  = self._UTC_string(strftime)
        self.julday      = self.UTC_time.julday

        return UTC_time
    def _LMSTify(self, strftime):

        """
        :copyright:
            Simon Stähler (mail@simonstaehler.com), 2018
            Martin van Driel (Martin@vanDriel.de), 2018
        :license:
            None
        """

        MIT              = (self.UTC_time - self.sol0) / self.sec_per_day_mars
        LMST_time        = UTCDateTime((MIT - 1)  * self.sec_per_day_earth)
        
        self.UTC_time    = self.UTC_time
        self.UTC_string  = self._UTC_string(strftime)
        self.julday      = self.UTC_time.julday
        
        self.LMST_time   = LMST_time
        self.LMST_string = self._LMST_string()
        self.sol         = int(self.LMST_string.split('M')[0])

        return LMST_time
    def _UTC_string(self, strftime):

        string = self.UTC_time.strftime(strftime)

        return string
    def _LMST_string(self):

        sol    = (self.LMST_time.datetime - UTCDateTime('1969-12-31T00:00:00.000000Z').datetime).days
        string = '%03dM%s' % (sol, self.LMST_time.strftime('%H:%M:%S'))

        return string
def marstime_list(sols_range=[10], hms='120000', hms_in_UTC=False):

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

    if len(sols_range)==1:
        now        = time.time()
        time_mars  = marstime(UTC=now)
        sols_range = [time_mars.sol-sols_range[0], time_mars.sol]

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

        if hms_in_UTC:
            time_mars_sol = marstime( LMST=UTCDateTime('1969-12-31T00:00:00.000000Z')+datetime.timedelta(days=int(sol)) )
            time_str_ymd  = time_mars_sol.UTC_time.strftime('%Y-%m-%d')
            time_str_HMS  = '%s:%s:%s' % (hms[0:2], hms[2:4], hms[4:6])
            time_mars     = marstime( time_str_ymd + 'T' + time_str_HMS )

        else:
            time_mars     = marstime( LMST='%03dM%s' % (sol,hms))

        print('%s    %s' % (time_mars.UTC_string, time_mars.LMST_string))
def time_funnel(time):

    """
    Small wrapper to convert times and catch errors.
    If no error occurs due to bad input, returns marstime object 
    within correct UTC / LMST times accesible via its attributes.
    """

    time     = str(time)
    time_new = None

    try:
        time_new = marstime(UTC=time)

    except ValueError:          # e.g. something like 25 hours etc ..
        pass

    except:                     # UTC conversion didn't work, so let's try if input was fiven as LMST
        try:
            time_new = marstime(LMST=time)

        except:                 # also didn't work, no matter why
            pass

    return time_new
def UVW2ZNE(stream, minimum_sample_length=2, inventory_file=None):

    """
    Rotate any UVW stream with arbitrary amount of traces per component 
    into ZNE-system.

    Passed stream object must contain UVW traces and no other traces.
    A new stream2 object is returned.
    """
    
    stream = stream.copy()
    stream = Stream2(stream)

    try:
        if stream.inventory:
            inv = stream.inventory
        else:
            raise Exception
    except:
        if inventory_file:
            stream.set_inventory(verbose=False, source=inventory_file)
        else:
            stream.set_inventory(verbose=False, source='IRIS')
        inv = stream.inventory

    pierced_streams = pierce_stream(stream, minimum_sample_length=minimum_sample_length)

    stRET = Stream2()
    for trU, trV, trW in zip(*pierced_streams):

        stWORK = Stream2(traces=[trU, trV, trW])
        stWORK.rotate('->ZNE', inventory=inv, components=('UVW')) 
        stRET += stWORK

    stRET.set_inventory(source=inv, verbose=False)
    stRET.sort()
    return stRET
def decimate_SEIS(trace, decimation_factor, verbose=True):

    """
    Downsample (decimate) original trace! That is,
    no copy is retured, the input trace is changed!
    """

    decimation_factor = int(decimation_factor)

    if decimation_factor is 1:
        return trace

    coeffs, delay = SEIS_FIR_coefficients(decimation_factor)
    if not coeffs:
        return trace

    trace_mean        = np.mean(trace)
    trace.data        = trace.data - trace_mean
    
    filtered_signal   = np.convolve(trace.data, coeffs, mode = "full")
    old_delta         = trace.stats.delta
    new_delta         = decimation_factor * old_delta
    index_time_sample = np.arange(trace.stats.npts)
    itimes            = index_time_sample[delay::int(decimation_factor)]
    mask              = (itimes < trace.stats.npts)
    itimes            = itimes[mask]
    
    trace.data        = filtered_signal[itimes]
    trace.data        = trace.data + trace_mean
    trace.stats.delta = new_delta

    if verbose:
        print(u'INFO: Decimated trace by factor %s.' % decimation_factor)
def SEIS_FIR_coefficients(decimation_factor):

    """
    Valid for VBB and SP.
    """

    if int(decimation_factor) is 2 :
        coeffs = [0.47312629, 0.31697387, 0.026678406, -0.10212776, -0.026100356
                , 0.057145566, 0.025159817, -0.036581896, -0.023883324, 0.024321463
                , 0.022318326, -0.016004425, -0.020517055, 0.0099688675, 0.018537275
                , -0.0054459628, -0.016433518, 0.0020369552, 0.014277564, 0.00050962088
                , -0.012123583, -0.0023524212, 0.010030696, 0.0036127516, -0.0080487542
                , -0.0043881042, 0.0062203603, 0.0047640391, -0.0045795627, -0.0048192143
                , 0.0031499979, 0.0046267062, -0.0019449759, -0.004253204, 0.00096769247
                , 0.0037597087, -0.00021371132, -0.0032004546, -0.00032992219, 0.0026274207
                , 0.00067918003, -0.0020938036, -0.00085202418, 0.0016962707, 0.00087003352
                , -0.0019099547, -0.0033559431, -0.00033131917]
    elif int(decimation_factor) is 3 :
        coeffs = [0.31535816, 0.26616585, 0.14575759, 0.017844718, -0.057989098
                                    , -0.061381415, -0.017455403, 0.027205415, 0.038973324, 0.016825307
                                    , -0.014348116, -0.027822647, -0.015969098, 0.0071774079, 0.0207357
                                    , 0.014920589, -0.0026327996, -0.015616039, -0.013714714, -0.00041046151
                                    , 0.01165056, 0.012385793, 0.0024693119, -0.0084637869, -0.010977356
                                    , -0.0038227537, 0.0058662295, 0.0095316246, 0.0046434188, -0.0037460821
                                    , -0.0080886949, -0.0050521158, 0.0020347144, 0.0066878777, 0.0051394254
                                    , -0.00068287796, -0.0053611984, -0.0049797446, -0.00035139307, 0.0041367821
                                    , 0.0046373457, 0.0011029486, -0.0030388678, -0.0041676573, -0.0016067328
                                    , 0.0020830208, 0.0036205212, 0.0018995546, -0.0012775299, -0.0030384576
                                    , -0.0020179669, 0.00062483293, 0.0024575219, 0.0019959896, -0.00012281968
                                    , -0.0019090844, -0.0018706999, -0.00023879687, 0.0014155055, 0.0016748102
                                    , 0.00047103246, -0.0010020733, -0.0014473109, -0.00058572053, 0.00071034848
                                    , 0.0012577202, 0.00060853921, -0.00079453795, -0.001507015, -0.0027735722
                                    , -0.0005772619]
    elif int(decimation_factor) is 4 :
        coeffs = [0.23653504, 0.2153236, 0.15848082, 0.083793312, 0.013365105, -0.034419119
                                    , -0.051046878, -0.039812922, -0.013070323, 0.013739644, 0.028548554
                                    , 0.026923735, 0.012592813, -0.0055519193, -0.018260375, -0.020213991
                                    , -0.011949001, 0.0011730746, 0.012122922, 0.015773855, 0.011154728
                                    , 0.0014588218, -0.0079624243, -0.012446117, -0.010243686, -0.0030995123
                                    , 0.0049446635, 0.0097706541, 0.0092410743, 0.0040912377, -0.0026885252
                                    , -0.0075423168, -0.0081806388, -0.0046254322, 0.00098837493, 0.0056574643
                                    , 0.0070931828, 0.0048204092, 0.00027697056, -0.0040607406, -0.0060094716
                                    , -0.0047603725, -0.0011885283, 0.0027206815, 0.0049573286, 0.0045121713
                                    , 0.0018078703, -0.0016159373, -0.0039624022, -0.0041275406, -0.002183456
                                    , 0.00072517502, 0.0030510179, 0.0036549442, 0.002359204, -3.6519799e-05
                                    , -0.0022335923, -0.0031325626, -0.0023753699, -0.0004719717, 0.0015254335
                                    , 0.002596593, 0.00226994, 0.00081899471, -0.00093285192, -0.0020744186
                                    , -0.002075599, -0.001027084, 0.00045596727, 0.0015875453, 0.0018238572
                                    , 0.0011186297, -9.3354625e-05, -0.0011566642, -0.0015448434, -0.0011173668
                                    , -0.00015881853, 0.00079965324, 0.0012694919, 0.0010489, 0.00030336727
                                    , -0.00054487761, -0.0010457698, -0.00094576413, -0.00029340666, 0.00060871453
                                    , 0.0013765346, 0.0017459453, 0.0017342302, 0.0017703171, -0.00078610331]
    elif int(decimation_factor) is 5 :
        coeffs = [0.18922436, 0.17825589, 0.14762825, 0.10361496, 0.054943368, 0.010695269
                                    , -0.021801252, -0.038539425, -0.039482988, -0.028172743, -0.010461548
                                    , 0.007229303, 0.019629899, 0.023857031, 0.019867271, 0.010076012, -0.0016730814
                                    , -0.011430031, -0.016343743, -0.015406003, -0.009560192, -0.0012011472
                                    , 0.0067212107, 0.01168697, 0.012363423, 0.0089212917, 0.0028455369, -0.003654296
                                    , -0.00841631, -0.010022298, -0.0081948247, -0.0038015433, 0.0015175294
                                    , 0.0059415763, 0.0080774855, 0.0073840916, 0.0042904904, -1.4044337e-05
                                    , -0.0040227529, -0.0064289104, -0.0065423511, -0.0044770911, -0.0010602982
                                    , 0.0024869707, 0.004979196, 0.0056554032, 0.0044051381, 0.0017638088
                                    , -0.00130429, -0.003753541, -0.0048074462, -0.0042014252, -0.002237522
                                    , 0.00034752686, 0.0026567061, 0.0039369017, 0.0038236496, 0.0024387422
                                    , 0.00031702407, -0.0017899896, -0.0031829374, -0.0034412947, -0.0025477768
                                    , -0.000868481, 0.00098881521, 0.0023946902, 0.002907305, 0.0024096994
                                    , 0.001126837, -0.00046738447, -0.0018248872, -0.0025087805, -0.0023346422
                                    , -0.0014139037, -0.00010080911, 0.0011430616, 0.0019060248, 0.001963608
                                    , 0.0013379015, 0.00027951988, -0.00083641545, -0.0016343265, -0.0018772145
                                    , -0.0015157261, -0.00071785948, 0.00022151141, 0.00098129292, 0.001313176
                                    , 0.001132468, 0.00053038984, -0.00026810885, -0.00099059893, -0.0014061579
                                    , -0.0014018975, -0.0010119241, -0.00039880717, 0.00021137166, 0.00060708891
                                    , 0.00066196767, 0.00036935217, -0.00016124418, -0.00074757589, -0.0011960072
                                    , -0.0013651208, -0.0012058818, -0.00077180623, -0.00019336015, 0.00036705163
                                    , 0.00076842657, 0.00092948624, 0.0022739838]
    elif int(decimation_factor) is 6 :
        coeffs = [0.16570362, 0.15829441, 0.13725868, 0.10594265, 0.06921123, 0.032528229
                                    , 0.00095806678, -0.021721859, -0.033611827, -0.034888856, -0.027580921
                                    , -0.014997372, -0.00094359554, 0.011083508, 0.018581118, 0.020429686
                                    , 0.016964789, 0.0097183716, 0.00092194974, -0.0070848679, -0.012422763
                                    , -0.014064534, -0.011991728, -0.0071014166, -0.00089137035, 0.0049562268
                                    , 0.00900786, 0.01040828, 0.0090432446, 0.0054955026, 0.00085260952
                                    , -0.0036210883, -0.006805297, -0.0079971515, -0.0070541361, -0.0043872511
                                    , -0.00080762105, 0.0027008946, 0.0052497657, 0.0062641921, 0.0056016045
                                    , 0.003556136, 0.00075419003, -0.0020313612, -0.004091803, -0.0049545085
                                    , -0.0044863801, -0.002905136, -0.00069879252, 0.0015243914, 0.0031950874
                                    , 0.0039263805, 0.0035993042, 0.0023769522, 0.00063838286, -0.0011335427
                                    , -0.0024846792, -0.0030999891, -0.0028772806, -0.0019367724, -0.00057383557
                                    , 0.00083092984, 0.0019156958, 0.0024286588, 0.0022825976, 0.0015665065
                                    , 0.00050730573, -0.00059655984, -0.0014599159, -0.0018816034, -0.0017901543
                                    , -0.0012505304, -0.00043752801, 0.00041982584, 0.0010982298, 0.0014395756
                                    , 0.0013833372, 0.00097796321, 0.00035366282, -0.00030664937, -0.00083971326
                                    , -0.0010991644, -0.0010712864, -0.00071174919, -0.00025581248, 0.00043810415
                                    , 0.00071178772, 0.0014833882, 0.00072977401, 0.0022810167, 0.0011542854]
    elif int(decimation_factor) is 7 :
        coeffs = [0.14208198, 0.13739452, 0.12389062, 0.10316056, 0.077608004, 0.050108083
                                    , 0.023616169, 0.00077113276, -0.016440101, -0.02690912, -0.030497886
                                    , -0.027981848, -0.020865075, -0.011100596, -0.00076022197, 0.0082820896
                                    , 0.014618799, 0.017472565, 0.016758278, 0.013016257, 0.007257143, 0.00074292615
                                    , -0.0052608857, -0.0097016394, -0.011911346, -0.011686133, -0.009284161
                                    , -0.0053405212, -0.00071840314, 0.003663674, 0.0070070671, 0.0087675247
                                    , 0.0087357759, 0.0070546297, 0.0041608354, 0.00068841781, -0.0026673293
                                    , -0.0052824821, -0.0067139948, -0.0067736702, -0.0055438336, -0.00334271
                                    , -0.00065325515, 0.0019840142, 0.0040728515, 0.005251891, 0.0053568939
                                    , 0.0044368859, 0.0027307249, 0.0006137155, -0.0014870598, -0.0031737122
                                    , -0.0041507659, -0.0042769266, -0.003582113, -0.0022478071, -0.00056846417
                                    , 0.0011144583, 0.0024819272, 0.0032929941, 0.0034265392, 0.0029004803
                                    , 0.0018556728, 0.00052273611, -0.00082722993, -0.0019346501, -0.0026047612
                                    , -0.0027381419, -0.0023428579, -0.001526694, -0.00047217472, 0.00060416642
                                    , 0.0014966615, 0.0020483169, 0.0021752114, 0.0018809054, 0.0012483622
                                    , 0.00042107934, -0.00043112988, -0.0011449105, -0.0015943006, -0.0017113318
                                    , -0.0014960549, -0.001011861, -0.0003693403, 0.0002983216, 0.00086243439
                                    , 0.0012237863, 0.0013277864, 0.0011733547, 0.00080783467, 0.0003159871
                                    , -0.00019982469, -0.00063952431, -0.00092532497, -0.0010140305, -0.00090339431
                                    , -0.00062933657, -0.00025537529, 0.00013938319, 0.00047917303, 0.00070110161
                                    , 0.00077338412, 0.00068326876, 0.00047101121, 0.00014812217, -0.00015542324
                                    , -0.0005331845, -0.00062111754, -0.0010503677, -0.00043352076, -0.001660049
                                    , -0.00079061883]
    else:
        print(u'WARNING: Only decimation factors 2 .. 7 are implmented currently. Choose another decimation factor or provide FIRs within the code.')
        return False, False

    delay  = len(coeffs)-1
    coeffs = coeffs[-1:0:-1] + coeffs
    return coeffs, delay


# Other
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
        print()
        print(u'ERROR: read_config: config file not found.')
        sys.exit()

    return params
def pierce_stream(stream, merge_before=False, ids=None, minimum_sample_length=2, sort_starttime=True, outfile=None):

    """
    Pierce stream into X streams, each of which contains all traces of ONE unique channel found in the input stream.
    All traces corresponding to each other across these X stream have the same start and end times. 
    For example useful, if you need 3-D component, equally trimmed data
    that are contained within one gappy stream ..

    NOTE: all unique traces in the stream are considered. So if there are 5 unique traces, the
    resulting streams will have the start and end times where all 5 traces are present, nothing more.
    So make sure the stream contains from the beginning only those traces that you'd like.

    New stream is sorted by start time (if wished) and stored to file (if wished).
    Returned is a tuple of streams each containing the equalled cut traces with unique channel ID.
    """


    ### PARAMETERS
    MERGE_METHOD                = 1
    MERGE_FILL_VALUE            = None
    MERGE_INTERPOLATION_SAMPLES = 0
    NEAREST_SAMPLE              = True


    ### handle overlaps first (if there any gaps, data will be masked array, that is why then split is applied to unmask the traces before piercing)
    if merge_before:
        stream = read2(file)
        stream.merge(method=MERGE_METHOD, fill_value=MERGE_FILL_VALUE, interpolation_samples=MERGE_INTERPOLATION_SAMPLES)
        stream = stream.split()

    return_stream = Stream()


    ### PIERCING TRACES
    if ids is None:
        ids = sorted(set( [tr.id for tr in stream] ))
    array         = np.array([[len(stream.select(id=id)),id] for id in ids])
    id_max_traces = array[np.argmax(array[:,0].astype('float')),1]

    
    for l, trace in enumerate( stream.select(id=id_max_traces) ):
        traces = []
        for tr in stream:
            if tr.id != trace.id:
                if tr.stats.endtime>=trace.stats.starttime and tr.stats.starttime<=trace.stats.endtime:
                    traces.append(tr)

        traces     = [trace] + traces
        tmp_stream = Stream2(traces=traces)
        tmp_stream = tmp_stream.copy()
        tmp_ids    = list( set( [tr.id for tr in tmp_stream] ))


        # not enough to pierce, moving on to next one
        if len(tmp_ids)<len(ids):
            continue
        

        # containing multiple traces of one ID
        if len(tmp_stream)>len(ids):

            for possibility in itertools.combinations(tmp_stream, len(ids)):

                possibility_ids = len( set( [trace.id for trace in possibility] ))

                if possibility_ids != len(ids):
                    continue

                else:
                    pierced_stream = Stream2(traces=possibility)        
                    pierced_stream = pierced_stream.copy()          # necessary as otherwise you also trim `possibilities` (which one doesn't want)
                    try:
                        pierced_stream.trim2(common=True, nearest_sample=NEAREST_SAMPLE)
                        if pierced_stream[0].stats.npts >= minimum_sample_length:
                            if all(pierced_stream[0].stats.starttime!=trace.stats.starttime or pierced_stream[0].stats.endtime!=trace.stats.endtime for trace in return_stream):
                                return_stream += pierced_stream
                    except ValueError:                              # can happen that `pierced_stream` has traces not overlapping at all, then continue
                        continue

        else:
            tmp_stream.trim2(common=True, nearest_sample=NEAREST_SAMPLE)
            if tmp_stream[0].stats.npts >= minimum_sample_length:
                return_stream += tmp_stream


    ### SORT STARTTIME EQUALLY FOR EACH ID (SO ORDER WITH RESPECT TO EACH OTHER IS KEPT)
    if sort_starttime:    
        sort_indices = np.argsort( [trace.stats.starttime for trace in return_stream.select(id=id_max_traces)])
        final_stream = Stream()
        for id in ids:
            traces        = return_stream.select(id=id).traces
            traces        = [traces[i] for i in sort_indices]
            final_stream += Stream(traces=traces)


    ### WRITE / RETURN
    if outfile is not None:
        final_stream.write(outfile)


    return tuple(final_stream.select(id=id) for id in ids)
def glitch_statistic(glitches, total_time_s_UTC):

    glitches       = np.array( glitches )

    if glitches.size != 0:
        t1                = UTCDateTime( datetime.datetime.utcnow() )
        t2                = t1 + float(total_time_s_UTC)
        total_LMST        = marstime(t2).LMST_time - marstime(t1).LMST_time
        
        glitch_number     = len(glitches)
        glitch_per_sol    = len(glitches)*86400/total_LMST
        glitch_U          = len(glitches[ glitches[:,5]=='1'])
        glitch_V          = len(glitches[ glitches[:,6]=='1'])
        glitch_W          = len(glitches[ glitches[:,7]=='1'])
        glitch_Uonly      = len(glitches[(glitches[:,5]=='1') & (glitches[:,6]=='0') & (glitches[:,7]=='0')])
        glitch_Vonly      = len(glitches[(glitches[:,5]=='0') & (glitches[:,6]=='1') & (glitches[:,7]=='0')])
        glitch_Wonly      = len(glitches[(glitches[:,5]=='0') & (glitches[:,6]=='0') & (glitches[:,7]=='1')])
        try:
            glitch_indivdual_ratio = (glitch_Uonly+glitch_Vonly+glitch_Wonly)/glitch_number*100
        except ZeroDivisionError:
            glitch_indivdual_ratio = np.nan

        return glitch_number, glitch_per_sol, glitch_U, glitch_V, glitch_W, glitch_Uonly, glitch_Vonly, glitch_Wonly, glitch_indivdual_ratio

    else:
        return 0, 0, 0, 0, 0, 0, 0, 0, np.nan
def merge_glitch_detector_files(glitch_detector_files, outfile, starttime_sort=True, multiples_out=None):

    """
    Merging together glitch extration files produced by `glitch_detector`.
    This helps combining into one file representing a certain period, e.g.
    """


    now                   = time.time()
    print()


    ### SANITY CHECK
    if not glitch_detector_files:
        print(u'ERROR: No glitch detector files for merging found.')
        sys.exit()

    if not outfile:
        print(u'WARNING: Must specify file name in config file to write out combined glitch detector file.')
        return


    ### OUTPUT
    print(u'Merging the following files:')
    for file in glitch_detector_files:
        print(file)


    ### RETRIEVE HEADER AND ALL DATA
    data    = []
    lengths = []
    with open(outfile, 'w') as fp: 
        counter = 0

        for file in glitch_detector_files:

            if file == outfile:
                print(u'WARNING: File name of target file and one input file is identic. This should not be.')
                return

            with open(file, 'r') as fp2:
                lines    = fp2.readlines() 
                header   = lines[0:3]
                data    += [line for line in lines[3:] if line.strip() and not line.startswith('#')]

                length   = [re.search(r'.*UTC length.*:\s*( \d*)\s*\w{0,4}.*\s*( \d{1,2}):(\d{2}):(\d{2})', line).groups() for line in lines[3:] if line.strip() and line.startswith('#') and 'UTC length' in line][0]
                lengths.append( length )


    ### SORTING
    if starttime_sort:
        data_sort    = np.array( [line.split() for line in data] )
        sort_indices = data_sort[:,1].argsort()
        data         = np.array(data)
        data         = data[sort_indices]


    ### MULTIPLES OUT
    if multiples_out:

        multiples_out = float(multiples_out)
        data_keep     = np.array( [line.split() for line in data] )
        discard       = []

        if starttime_sort:
            for i in range(len(data_keep)):
                j = i+1
                while True:
                    try:
                        if UTCDateTime(data_keep[i][1])+multiples_out>UTCDateTime(data_keep[j][1]):
                            discard.append(j)
                            j += 1
                        else:
                            break
                    except IndexError:
                        break

        else:
            for i in range(len(data_keep)):
                for j in range(len(data_keep)):
                    if i==j:
                        continue
                    if UTCDateTime(data_keep[i][1])+multiples_out>UTCDateTime(data_keep[j][1]):
                        discard.append(j)

        print(u'Discarded %s glitches that occured within %s s to another glitch.' % (len(discard),multiples_out))
        data = np.array(data)
        keep = [i for i in range(len(data_keep)) if i not in discard]
        data = data[keep]


    ### STATISTICS    
    data_split = np.array( [line.split() for line in data] )
    total_UTC = 0

    for days, hours, minutes, seconds in lengths:
        #print(days, hours, minutes, seconds)
        if days.rstrip():
            total_UTC += int(days) * 24*3600
        total_UTC += int(hours) * 3600
        total_UTC += int(minutes) * 60
        total_UTC += int(seconds) 


    statistic  = glitch_statistic(data_split, total_UTC)

    output1    = [u'MERGED GLITCH DETECTOR FILES:']
    output2    = ['  %s' % file for file in glitch_detector_files]
    output3    = [u'',
                  u'RESULTS:',      
                  u'  Glitches total:            %s'           % statistic[0],
                  u'       all / sol:            %.0f'         % statistic[1],
                  u'            on U:            %s'           % statistic[2],
                  u'            on V:            %s'           % statistic[3],
                  u'            on W:            %s'           % statistic[4],
                  u'       only on U:            %s'           % statistic[5],
                  u'       only on V:            %s'           % statistic[6],
                  u'       only on W:            %s'           % statistic[7],
                  u'    indi./ all %%:            %.1f'        % statistic[8],
                  u'',
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
    print(u'Merged glitch file:')
    print(outfile)
def sec2hms(seconds, digits=0):

    """
    Convert seconds given as float into hours, minutes and seconds.
    Optional 'digits' determines position after decimal point.

    Returns string.
    """

    if np.isnan(seconds):
        return np.nan


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
def download_data(outdir=os.getcwd(), 
    starttime='2019-04-10T00:00:00', 
    endtime='2019-04-11T00:00:00', 
    network='XB', 
    station='ELYSE', 
    location='02', 
    channel='BH?', 
    source='IRIS',
    username='', 
    password='', 
    format_DATA='MSEED', 
    format_INV='STATIONXML'):

    """
    :copyright:
        Simon Staehler (simon.staehler@erdw.ethz.ch), 2018
        with a fair amount of changes by John-Robert Scholz, 2020.
    """

    class WrappedClient(Client):

        """
        Wrapped client to make data request work for restricted InSight data
        at the SEIS data portal. Data are restricted usually for three 3 months.

        :copyright:
        Martin van Driel (vandriel@erdw.ethz.ch), 2018
        """

        # hacking the Basic authentication into the fdsn client
        def _set_opener(self, user, password, digest=False):
            # Only add the authentication handler if required.
            handlers = []
            if user is not None and password is not None:
                # Create an OpenerDirector for HTTP Digest Authentication
                password_mgr = urllib.request.HTTPPasswordMgrWithDefaultRealm()
                password_mgr.add_password(None, self.base_url, user, password)

                handlers.append(urllib.request.HTTPBasicAuthHandler(password_mgr))

            if (user is None and password is None) or self._force_redirect is True:
                # Redirect if no credentials are given or the force_redirect
                # flag is True.
                handlers.append(CustomRedirectHandler())
            else:
                handlers.append(NoRedirectionHandler())

            # Don't install globally to not mess with other codes.
            self._url_opener = urllib.request.build_opener(*handlers)
            if self.debug:
                print('Installed new opener with handlers: {!s}'.format(handlers))

        # MSDS uses query also for authenticated queries
        def _build_url(self, service, resource_type, parameters={}):
            return build_url(self.base_url, service, self.major_versions[service],
                             resource_type, parameters,
                             service_mappings=self._service_mappings)

    # Variables
    if not outdir:
        outdir = os.getcwd()
    os.makedirs(outdir, exist_ok=True)

    starttime       = time_funnel(starttime).UTC_time
    endtime         = time_funnel(endtime).UTC_time
    network         = str(network)
    station         = str(station)
    location        = str(location)
    channel         = str(channel)
    source          = str(source)

    # Client
    if username and password:
        if source.upper()=='IPGP':
            client = WrappedClient('https://ws.seis-insight.eu', user=str(username), password=str(password))
        else:
            client = Client(source, user=str(username), password=str(password))
    else:
        client = Client(source)

    # Waveform data
    st = client.get_waveforms(network   = network,
                              station   = station,
                              location  = location,
                              channel   = channel,
                              starttime = starttime,
                              endtime   = endtime)
    st = Stream2(st)
    st.set_inventory(source=source, verbose=False)
    print()

    # Paths
    times            = st.times
    request          = '%s.%s.%s.%s' % (network, station, location, channel)
    outfile_inv      = os.path.join(outdir, 'inventory_%s-%s.xml' % (st.times[0].strftime('%Y%m%dT%H%M'),st.times[1].strftime('%Y%m%dT%H%M')))
    outfile_raw      = os.path.join(outdir, '%s_%s_%s_raw.%s'     % (request,st.times[0].strftime('%Y%m%dT%H%M'),st.times[1].strftime('%Y%m%dT%H%M'),format_DATA))

    # Processing
    st.inventory.write(outfile_inv, format=format_INV)
    print(u'INFO: written inventory file to:    %s' % outfile_inv)

    st.write(outfile_raw, format=format_DATA)
    print(u'INFO: written raw waveform data to: %s' % outfile_raw)
    print(u'Finished.')
def process_data(*waveform_files,
    inventory_file='IPGP',
    gain_correction=False, 
    remove_response={'unit':None, 'pre_filt':None, 'water_level':60},
    rotate2zne=False,
    decimate={'decimation_factors':[], 'final_location':'', 'final_bandcode':''}):


    ### Variables
    unit        = remove_response['unit']
    pre_filt    = remove_response['pre_filt']
    water_level = remove_response['water_level']


    ### Processing
    for file in waveform_files:

        stream    = read2(file)
        format_in = stream[0].stats._format
        stream.set_inventory(inventory_file, verbose=False)
        print()
        print(u'File: %s' % file)
        print(u'Found %s traces. Handling ..' % len(stream))


        # gain correction
        if gain_correction:
            stream.gain_correction(channel='?[LMH]?')
            print(u'INFO: gain corrected data')
        
            if rotate:
                stream = UVW2ZNE(stream, minimum_sample_length=2)
                print(u'INFO: rotated gain corrected data')


        # Removal instrument response
        elif remove_response['unit'] and remove_response['unit'].lower()!='none':
            stream.detrend(type='demean')
            stream.taper(0.1)
            stream.detrend(type='simple')
            for tr in stream:
                try:
                    tr.remove_response(stream.inventory, output=unit, pre_filt=pre_filt, water_level=water_level)
                except ValueError:
                    print(u'WARNING: Could not remove response to %s for channel: %s' % (remove_response, tr))
            print(u'INFO: instrument corrected data')

            if rotate:
                stream = UVW2ZNE(stream, minimum_sample_length=2)
                print(u'INFO: rotated instrument corrected data')


        # Decimation
        for trace in stream.select(channel='?[LMH]?'):
            for decimation_factor in decimate['decimation_factors']:
                decimate_SEIS(trace, decimation_factor, verbose=False)
                print(u'INFO: decimated data by factor %s' % decimation_factor)

            if decimate['final_location']:
                trace.stats.location = decimate['final_location']

            if decimate['final_bandcode']:
                trace.stats.channel = decimate['final_bandcode'] + trace.stats.channel[1:]


        # Output
        outfile = '.'.join(file.split('.')[:-1]) + '_processed.' + format_in.lower()
        stream.write(outfile)

        print(u'File out: %s' % outfile)
        print()

    print(u'Finished.')


# Ppol (P-wave Polarization)
class ppol():

    """
    If stream object is passed, note that ppol does not currently handle gaps/overlaps.
    In such case, only the first trace of each component is processed. This 
    is true also if `starttime` and/or `endtime` are passed (only the FIRST
    trace of each component corresponding to these times is processed).

    So, no merging / interpolation of any kind is implemented.


    results = PHI_2D,     PHI_err_2D, PHI_3D,     PHI_err_3D, \
              INC_2D,     INC_err_2D, INC_3D,     INC_err_3D, \
              SNR_HOR_2D, SNR_3D,     SNR_RZp_2D, SNR_RZp_3D, \
              POL_HOR_2D, POL_3D,     POL_RZp_2D, POL_RZp_3D, \
              eig_vec_2D, eig_val_2D, eig_vec_3D, eig_val_3D, \
              data_R_2D,  data_T_2D,  data_R_3D,  data_T_3D
    """

    def __init__(self, comp_N=None, comp_E=None, comp_Z=None, stream=None, starttime=None, endtime=None, demean=True, **kwargs):

        """
        `kwargs` are passed to `calc` and `plot `methods.
        """

        if starttime is not None:
            self.starttime = UTCDateTime(starttime)
        else:
            self.starttime = starttime
        if endtime is not None:
            self.endtime = UTCDateTime(endtime)
        else:
            self.endtime = endtime
        
        self.trace_N, self.trace_E, self.trace_Z = self._assign_traces(comp_N=comp_N, comp_E=comp_E, comp_Z=comp_Z, stream=stream)    # returns a copy of the input data
        self.traces_input = [self.trace_N.copy(), self.trace_E.copy(), self.trace_Z.copy()]  # have effect only for plotting and only if start- & endtime are given, so one can zoom out from the ppol data! ;)

        if self.starttime is not None:
            self.trace_N.trim(starttime=self.starttime)
            self.trace_E.trim(starttime=self.starttime)
            self.trace_Z.trim(starttime=self.starttime)
        if self.endtime is not None:
            self.trace_N.trim(endtime=self.endtime)
            self.trace_E.trim(endtime=self.endtime)
            self.trace_Z.trim(endtime=self.endtime)            

        if demean:
            mean_N = np.mean(self.trace_N.data)
            mean_E = np.mean(self.trace_E.data)
            mean_Z = np.mean(self.trace_Z.data)                # if Z is empty trace object, mean will be np.nan

            self.trace_N.data = self.trace_N.data - mean_N
            self.trace_E.data = self.trace_E.data - mean_E
            self.trace_Z.data = self.trace_Z.data - mean_Z     # if Z is empty trace object, this operation leaves the it still empty

            # for plotting, these data will be demeaned the same way as the data where the ppol calc is run on. What is plotted will hence correspond to mathematical calculation.
            self.traces_input[0].data = self.traces_input[0].data - mean_N
            self.traces_input[1].data = self.traces_input[1].data - mean_E
            self.traces_input[2].data = self.traces_input[2].data - mean_Z     # if Z is empty trace object, this operation leaves the it still empty

            self.demeaned = True
        else:
            self.demeaned = False

        self._sanity_check()
        self.results = self.calc(**kwargs)
    def __str__(self, tag='', print_3D=True):

        """
        When calling print(...)
        """

        if self.fix_angles == 'AMP':
            info = '(in this mode, radial-transverse data are ambiguous but BAZs & INCs are correct)'
        else:
            info = ''

        if print_3D:
            string =    u'\n' \
                      + u'PPOL CALCULATIONS: %s   \n' % tag \
                      + u'   %11s : %s %s         \n' % ('Fix angles', self.fix_angles, info) \
                      + u'   %11s : %s            \n' % ('demeaned',   self.demeaned)   \
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
                      + u'   %11s : %s            \n' % ('EigVecs_2D', self.results[16][0]) \
                      + u'   %11s : %s            \n' % (''          , self.results[16][1]) \
                      + u'   %11s : %s            \n' % ('EigVals_2D', self.results[17]) \
                      + u'   %11s : %s            \n' % ('EigVecs_3D', self.results[18][0]) \
                      + u'   %11s : %s            \n' % (''          , self.results[18][1]) \
                      + u'   %11s : %s            \n' % (''          , self.results[18][2]) \
                      + u'   %11s : %s              ' % ('EigVals_3D', self.results[19])
        else:
            string =    u'\n' \
                      + u'PPOL CALCULATIONS: %s   \n' % tag \
                      + u'   %11s : %s            \n' % ('Fix angles', self.fix_angles) \
                      + u'   %11s : %s            \n' % ('demeaned',   self.demeaned) \
                      + u'   %11s : %s            \n' % ('NPTS',       len(self.trace_N.data)) \
                      + u'   %11s : %-5.1f ± %3.1f\n' % ('BAZ_2D',     self.results[0], self.results[1]) \
                      + u'   %11s : %-10.1g       \n' % ('SNR_HOR_2D', self.results[8]) \
                      + u'   %11s : %-10.3f       \n' % ('POL_HOR_2D', self.results[12]) \
                      + u'   %11s : %s            \n' % ('EigVecs_2D', self.results[16][0]) \
                      + u'   %11s : %s            \n' % (''          , self.results[16][1]) \
                      + u'   %11s : %s              ' % ('EigVals_2D', self.results[17])

        return string
    def _assign_traces(self, comp_N=None, comp_E=None, comp_Z=None, stream=None):

        """
        """

        if stream is not None:
            stream_N = stream.select(component='N') or stream.select(component='T') or stream.select(component='U') or stream.select(component='1') or stream.select(component='T')
            stream_E = stream.select(component='E') or stream.select(component='L') or stream.select(component='V') or stream.select(component='2') or stream.select(component='R')
            stream_Z = stream.select(component='Z') or stream.select(component='Q') or stream.select(component='W') or stream.select(component='3')

            if not stream_N or not stream_E:
                print(u'WARNING: No idea how to perform polarization analysis on components: %s' % ', '.join( [tr.id for tr in stream] ))
                sys.exit()

            if self.starttime is not None:
                for trace in stream_N:
                    if trace.slice(starttime=self.starttime):
                        trace_N = trace
                        break
                else:
                    print(u'WARNING: No trace found in stream that corresponds to passed starttime %s.' % self.starttime)
                    sys.exit()
                for trace in stream_E:
                    if trace.slice(starttime=self.starttime):
                        trace_E = trace
                        break
                else:
                    print(u'WARNING: No trace found in stream that corresponds to passed starttime %s.' % self.starttime)
                    sys.exit()                        
                if stream_Z is not None:
                    for trace in stream_Z:
                        if trace.slice(starttime=self.starttime):
                            trace_Z = trace
                            break
                    else:
                        print(u'WARNING: No trace found in stream that corresponds to passed starttime %s.' % self.starttime)
                        sys.exit() 
                else:
                    trace_Z = Trace()

            elif self.endtime is not None:
                for trace in stream_N:
                    if trace.slice(starttime=self.endtime):
                        trace_N = trace
                        break
                else:
                    print(u'WARNING: No trace found in stream that corresponds to passed endtime %s.' % self.endtime)
                    sys.exit()                        
                for trace in stream_E:
                    if trace.slice(starttime=self.endtime):
                        trace_E = trace
                        break
                else:
                    print(u'WARNING: No trace found in stream that corresponds to passed endtime %s.' % self.endtime)
                    sys.exit()                          
                if stream_Z is not None:
                    for trace in stream_Z:
                        if trace.slice(starttime=self.endtime):
                            trace_Z = trace
                            break
                    else:
                        print(u'WARNING: No trace found in stream that corresponds to passed endtime %s.' % self.endtime)
                        sys.exit()                              
                else:
                    trace_Z = Trace()

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
                trace_N = Trace(data=np.asarray(comp_N), header={'channel':'N'})        # sampling rate of 1 Hz by default
            else:
                print(u'`comp_N` must be either a ObsPy `Trace` object or a `list`-like object containing your waveform data.')
                sys.exit()

            if isinstance(comp_E, Trace):
                trace_E = comp_E
            elif isinstance(comp_E, (list, tuple, np.ndarray)):
                trace_E = Trace(data=np.asarray(comp_E), header={'channel':'E'})        # sampling rate of 1 Hz by default
            else:
                print(u'`comp_E` must be either a ObsPy `Trace` object or a `list`-like object containing your waveform data.')
                sys.exit()

            if comp_Z is not None:
                if isinstance(comp_Z, Trace):
                    trace_Z = comp_Z
                elif isinstance(comp_Z, (list, tuple, np.ndarray)):
                    trace_Z = Trace(data=np.asarray(comp_Z), header={'channel':'Z'})    # sampling rate of 1 Hz by default
                else:
                    print(u'`comp_Z` must be either a ObsPy `Trace` object or a `list`-like object containing your waveform data.')
                    sys.exit()                    
            else:
                trace_Z = Trace()           # trace object with empty data array


        else:
            print(u'You must either specify an ObsPy `stream` object containing at least two components')
            print(u'or `comp_N` and `comp_E` (both either as ObsPy `Trace` objects or lists).')
            sys.exit()


        return trace_N.copy() , trace_E.copy(), trace_Z.copy()
    def _sanity_check(self):

        
        len_N   = self.trace_N.stats.npts
        len_E   = self.trace_E.stats.npts
        len_Z   = self.trace_Z.stats.npts

        start_N = self.trace_N.stats.starttime
        start_E = self.trace_E.stats.starttime
        start_Z = self.trace_Z.stats.starttime


        # N or E have no data
        if len_E==0 or len_N==0:
            print(u'ERROR: One or multiple components have no data. Cannot perform ppol analysis.')
            sys.exit()


        # Comps have not same data length
        if len_Z != 0:
            if len_N!=len_E or len_N!=len_Z or len_E!=len_Z:
                print(u'ERROR: Data lengths of components are not equal. Cannot perform ppol analysis.')
                sys.exit()

        else:
            if len_N!=len_E:
                print(u'ERROR: Data lengths of components are not equal. Cannot perform ppol analysis.')
                sys.exit()


        # Comps do not start at same time - Warning only
        if len_Z != 0:
            if start_N!=start_E or start_N!=start_Z or start_E!=start_Z:
                print(u'WARNING: Data do not start at the same time. Ppol analysis is performed nevertheless.')
                print('  '+self.trace_N.__str__())
                print('  '+self.trace_E.__str__())
                print('  '+self.trace_Z.__str__())

        else:
            if start_N!=start_E:
                print(u'WARNING: Data do not start at the same time. Ppol analysis is performed nevertheless.')
                print('  '+self.trace_N.__str__())
                print('  '+self.trace_E.__str__())
    def calc(self, bias=False, fix_angles='EQ', Xoffset_samples_for_amplitude=0, **kwargs):

        """
        DO INDIVIDUAL PPOL MEASUREMENTS

        Take data arrays and perform 2-D and 3-D principle component
        analysis (2-D and 3-D PCA). Components must be of equal length
        and correspond to same start and thus end times.


        Useful for seismology, for example.

               ^  
               | data_1
               |
               |
               |
               o---------> data_2
             data_Z 
         (pointing to you, left-hand rule)


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
        self.fix_angles = fix_angles.upper()
        Xoffset_samples_for_amplitude = int(Xoffset_samples_for_amplitude)

        data_1 = self.trace_N.data
        data_2 = self.trace_E.data
        data_Z = self.trace_Z.data



        ### IF NO Z-COMPONENT GIVEN
        if data_Z is None:
     

            ### 2-D phi, horizontal plane
            covariance_matrix_2D_hori      = covariance_matrix(data_1, data_2, bias=bias)                                  # does be default no demeaning internally!
            eig_val_2D, eig_vec_2D         = np.linalg.eig(covariance_matrix_2D_hori)
            index_array_descending_eig_val = np.abs(eig_val_2D.argsort()[::-1])
            eig_val_2D                     = eig_val_2D[index_array_descending_eig_val]                                    # Eig-values descending
            eig_vec_2D                     = eig_vec_2D[:,index_array_descending_eig_val]                                  # Eig-vectors sorted acc. to E-values
            eig_vec_2D_1                   = eig_vec_2D[:,0]
            
            # Derived
            PHI_2D                         = (np.arctan2( eig_vec_2D_1[1].real, eig_vec_2D_1[0].real ) * 180/np.pi )%360
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
            eig_val_3D                     = np.nan*np.ones(3)                                                            # Eig-values descending
            eig_vec_3D                     = np.nan*np.ones(3)                                     


            ### AMBIGUITY 180°
            if self.fix_angles=='EQ':                # Nothing we can do solving the 180° ambiguity in 2D (this code doesn't make use of first motion information)
                pass


            elif self.fix_angles=='AMP':
                data_1_offest = data_1[Xoffset_samples_for_amplitude:]
                data_2_offest = data_2[Xoffset_samples_for_amplitude:]
                
                amp_1         = data_1_offest[np.argmax(np.abs(data_1_offest))]
                amp_2         = data_2_offest[np.argmax(np.abs(data_2_offest))]
                
                PHI_2D_OLD    = PHI_2D

                # 2-D PHI
                if abs(amp_1)>=abs(amp_2):              # `data_1` more significant than `data_2` 
                    if amp_1>=0:                        # `data_1` positive
                        if PHI_2D>=90 and PHI_2D<180 or PHI_2D>=270 and PHI_2D<360:
                            PHI_2D = PHI_2D%180 + 180
                        else:
                            PHI_2D = PHI_2D%180
                    else:                               # `data_1` negative
                        if PHI_2D>=90 and PHI_2D<180 or PHI_2D>=270 and PHI_2D<360:
                            PHI_2D = PHI_2D%180
                        else:
                            PHI_2D = PHI_2D%180 + 180
                else:                                   # `data_2` more significant than `data_1` 
                    if amp_2>=0:                        # `data_2` positive
                        PHI_2D = PHI_2D%180
                    else:                               # `data_2` negative
                        PHI_2D = PHI_2D%180 + 180        

                # correct radial and transverse data
                data_R_2D, data_T_2D = rotate.rotate_ne_rt(data_1, data_2, PHI_2D)


            else:
                pass



        ### IF Z-COMPONENT GIVEN
        else:


            ### 2-D phi, horizontal plane
            covariance_matrix_2D_hori      = covariance_matrix(data_1, data_2, bias=bias)                               # does be default no demeaning internally!
            eig_val_2D, eig_vec_2D         = np.linalg.eig(covariance_matrix_2D_hori)
            eig_val_2D                     = np.abs(eig_val_2D)                                                         # Veeery mall eig-values may be negative
            index_array_descending_eig_val = eig_val_2D.argsort()[::-1]
            eig_val_2D                     = eig_val_2D[index_array_descending_eig_val]                                 # Eig-values descending
            eig_vec_2D                     = eig_vec_2D[:,index_array_descending_eig_val]                               # Eig-vectors sorted acc. to E-values
            eig_vec_2D_1                   = eig_vec_2D[:,0]
            
            # Derived
            PHI_2D                         = (np.arctan2( eig_vec_2D_1[1].real, eig_vec_2D_1[0].real ) * 180/np.pi ) % 360
            PHI_err_2D                     = np.arctan( np.sqrt( eig_val_2D[1]/eig_val_2D[0] )) * 180/np.pi             # Reymond (2010)
            SNR_HOR_2D                     = (eig_val_2D[0] - eig_val_2D[1]) / eig_val_2D[1]                            # De Meersman et al. (2006)
            POL_HOR_2D                     = 1 - eig_val_2D[1]/eig_val_2D[0]                                            # rectilinearity in horizontal plane (1 for linearised, 0 for circular polarisation). Jurkevics (1988)
            data_R_2D, data_T_2D           = rotate.rotate_ne_rt(data_1, data_2, PHI_2D)


            ### 2-D phi, radial-vertical plane
            data_R_2D, data_T_2D             = rotate.rotate_ne_rt(data_1, data_2, PHI_2D)
            covariance_matrix_2D_radZ        = covariance_matrix(data_Z, data_R_2D, bias=bias)                          # does be default no demeaning internally!
            eig_val_2D_radZ, eig_vec_2D_radZ = np.linalg.eig(covariance_matrix_2D_radZ)
            eig_val_2D_radZ                  = np.abs(eig_val_2D_radZ)                                                  # Veeery mall eig-values may be negative
            index_array_descending_eig_val   = np.argsort( eig_val_2D_radZ )[::-1]
            eig_val_2D_radZ                  = eig_val_2D_radZ[index_array_descending_eig_val]                          # Eig-values descending
            eig_vec_2D_radZ                  = eig_vec_2D_radZ[:,index_array_descending_eig_val]                        # Eig-vectors sorted acc. to E-values
            eig_vec_2D_radZ_1                = eig_vec_2D_radZ[:,0]

            # derived
            INC_2D                           = np.arctan( eig_vec_2D_radZ_1[1] / eig_vec_2D_radZ_1[0] ) * 180/np.pi
            INC_err_2D                       = np.arctan( np.sqrt( eig_val_2D_radZ[1]/eig_val_2D_radZ[0] )) * 180/np.pi
            SNR_RZp_2D                       = (eig_val_2D_radZ[0] - eig_val_2D_radZ[1]) / eig_val_2D_radZ[1]           # De Meersman et al. (2006)
            POL_RZp_2D                       = 1 - eig_val_2D_radZ[1]/eig_val_2D_radZ[0]                                # rectilinearity in radial-vertical plane (1 for linearised, 0 for circular polarisation). Jurkevics (1988)


            ### 3-D
            covariance_matrix_3D             = covariance_matrix(data_1, data_2, data_Z, bias=bias)                     # does be default no demeaning internally!
            eig_val_3D, eig_vec_3D           = np.linalg.eig(covariance_matrix_3D)
            eig_val_3D                       = np.abs(eig_val_3D)                                                       # Veeery mall eig-values may be negative
            index_array_descending_eig_val   = np.argsort( eig_val_3D )[::-1]
            eig_val_3D                       = eig_val_3D[index_array_descending_eig_val]                               # Eig-values descending
            eig_vec_3D                       = eig_vec_3D[:,index_array_descending_eig_val]                             # Eig-vectors sorted acc. to E-values
            eig_vec_3D_1                     = eig_vec_3D[:,0]

            # derived
            PHI_3D                           = (np.arctan2( eig_vec_3D_1[1].real, eig_vec_3D_1[0].real ) * 180/np.pi) % 360
            PHI_err_3D                       = np.arctan( np.sqrt( eig_val_3D[2]/(eig_val_3D[1]+eig_val_3D[0]) )) * 180/np.pi
            SNR_3D                           = abs((eig_val_3D[0] - (eig_val_3D[1] + eig_val_3D[2]))) / (eig_val_3D[1] + eig_val_3D[2])  # De Meersman et al. (2006)        
            POL_3D                           = 1 - ( eig_val_3D[1]+eig_val_3D[2] )/( 2*eig_val_3D[0] )                  # rectilinearity in 3-D. (1 for linearised, 0 for circular polarisation). Jurkevics (1988)
            data_R_3D, data_T_3D             = rotate.rotate_ne_rt(data_1, data_2, PHI_3D)


            ### 3-D phi, radial & Z-data plane
            data_R_3D, data_T_3D             = rotate.rotate_ne_rt(data_1, data_2, PHI_3D)
            covariance_matrix_3D_radZ        = covariance_matrix(data_Z, data_R_3D, bias=bias)                          # does be default no demeaning internally!
            eig_val_3D_radZ, eig_vec_3D_radZ = np.linalg.eig(covariance_matrix_3D_radZ)
            eig_val_3D_radZ                  = np.abs(eig_val_3D_radZ)                                                  # Veeery mall eig-values may be negative
            index_array_descending_eig_val   = np.argsort( eig_val_3D_radZ )[::-1]
            eig_val_3D_radZ                  = eig_val_3D_radZ[index_array_descending_eig_val]                          # Eig-values descending
            eig_vec_3D_radZ                  = eig_vec_3D_radZ[:,index_array_descending_eig_val]                        # Eig-vectors sorted acc. to E-values
            eig_vec_3D_radZ_1                = eig_vec_3D_radZ[:,0]

            # derived
            INC_3D                           = np.arctan( eig_vec_3D_radZ_1[1] / eig_vec_3D_radZ_1[0] ) * 180/np.pi
            INC_err_3D                       = np.arctan( np.sqrt( eig_val_3D_radZ[1]/eig_val_3D_radZ[0] )) * 180/np.pi
            SNR_RZp_3D                       = (eig_val_3D_radZ[0] - eig_val_3D_radZ[1]) / eig_val_3D_radZ[1]           # De Meersman et al. (2006)
            POL_RZp_3D                       = 1 - eig_val_3D_radZ[1]/eig_val_3D_radZ[0]                                # rectilinearity in radial-vertical plane (1 for linearised, 0 for circular polarisation). Jurkevics (1988)


            ### AMBIGUITY 180°
            if self.fix_angles=='EQ':    # Correct baz must deliver incidence angle>0, therefore can solve ambiguity
                if INC_2D < 0:
                    PHI_2D               = (PHI_2D+180)%360
                    INC_2D               = abs(INC_2D)
                    data_R_2D, data_T_2D = rotate_2D(data_R_2D, data_T_2D, 180)
                if INC_3D < 0:
                    PHI_3D               = (PHI_3D+180)%360
                    INC_3D               = abs(INC_3D)
                    data_R_3D, data_T_3D = rotate_2D(data_R_3D, data_T_3D, 180)


            elif self.fix_angles=='AMP':
                data_1_offest = data_1[Xoffset_samples_for_amplitude:]
                data_2_offest = data_2[Xoffset_samples_for_amplitude:]
                data_Z_offest = data_Z[Xoffset_samples_for_amplitude:]
                
                amp_1         = data_1_offest[np.argmax(np.abs(data_1_offest))]
                amp_2         = data_2_offest[np.argmax(np.abs(data_2_offest))]
                amp_Z         = data_Z_offest[np.argmax(np.abs(data_Z_offest))]
                
                PHI_2D_OLD    = PHI_2D
                PHI_3D_OLD    = PHI_3D

                # 2-D PHI
                if abs(amp_1)>=abs(amp_2):              # `data_1` more significant than `data_2` 
                    if amp_1>=0:                        # `data_1` positive
                        if PHI_2D>=90 and PHI_2D<180 or PHI_2D>=270 and PHI_2D<360:
                            PHI_2D = PHI_2D%180 + 180
                        else:
                            PHI_2D = PHI_2D%180
                    else:                               # `data_1` negative
                        if PHI_2D>=90 and PHI_2D<180 or PHI_2D>=270 and PHI_2D<360:
                            PHI_2D = PHI_2D%180
                        else:
                            PHI_2D = PHI_2D%180 + 180
                else:                                   # `data_2` more significant than `data_1` 
                    if amp_2>=0:                        # `data_2` positive
                        PHI_2D = PHI_2D%180
                    else:                               # `data_2` negative
                        PHI_2D = PHI_2D%180 + 180                  
    
                # 3-D PHI
                if abs(amp_1)>=abs(amp_2):              # `data_1` more significant than `data_2` 
                    if amp_1>=0:                        # `data_1` positive
                        if PHI_3D>=90 and PHI_3D<180 or PHI_3D>=270 and PHI_3D<360:
                            PHI_3D = PHI_3D%180 + 180
                        else:
                            PHI_3D = PHI_3D%180
                    else:                               # `data_1` negative
                        if PHI_3D>=90 and PHI_3D<180 or PHI_3D>=270 and PHI_3D<360:
                            PHI_3D = PHI_3D%180
                        else:
                            PHI_3D = PHI_3D%180 + 180
                else:                                   # `data_2` more significant than `data_1` 
                    if amp_2>=0:                        # `data_2` positive
                        PHI_3D = PHI_3D%180
                    else:                               # `data_2` negative
                        PHI_3D = PHI_3D%180 + 180        

                # correct radial and transverse data
                data_R_2D, data_T_2D = rotate.rotate_ne_rt(data_1, data_2, PHI_2D)
                data_R_3D, data_T_3D = rotate.rotate_ne_rt(data_1, data_2, PHI_3D)

                # 2-D INC
                INC_2D_OLD = INC_2D
                if amp_Z>=0:
                    INC_2D = np.abs(INC_2D) % 180
                else:
                    INC_2D = -1*np.abs(INC_2D) % 180

                # 3-D INC
                if amp_Z>=0:
                    INC_3D = np.abs(INC_3D) % 180
                else:
                    INC_3D = -1*np.abs(INC_3D) % 180


            else:
                pass
        


        ### RESULTS
        results = PHI_2D,     PHI_err_2D, PHI_3D,     PHI_err_3D, \
                  INC_2D,     INC_err_2D, INC_3D,     INC_err_3D, \
                  SNR_HOR_2D, SNR_3D,     SNR_RZp_2D, SNR_RZp_3D, \
                  POL_HOR_2D, POL_3D,     POL_RZp_2D, POL_RZp_3D, \
                  eig_vec_2D, eig_val_2D, eig_vec_3D, eig_val_3D, \
                  data_R_2D,  data_T_2D,  data_R_3D,  data_T_3D



        ### ASSIGN / RETURN RESULTS
        self.results = results
        return self.results
    def plot(self, title='', verts=(), outfile=None, show=True, original_data=False, **kwargs):

        """
        PLOT INDIVIDUAL PPOL MEASUREMENT

        Either provide numpy lists or Obspy traces.
        Always plots 2-D angle
        Expects that all traces have same start time and same amount of samples!

        `original_data` is plotted as provided. No demeaning or alike is done.
        `kwargs` are passed to the quick_plot function.
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
        BAZ_3D     = self.results[2]
        BAZ_err_3D = self.results[3]
        INC_3D     = self.results[6]
        INC_err_3D = self.results[7]
        eigvecs    = self.results[18]   # eig_vec_3D
        data_R     = self.results[20]   # comp_R_2D


        ### PLOT WAVEFORMS
        # Figure instance
        if self.starttime is not None:
            verts += (self.starttime, )
        if self.endtime is not None:
            verts += (self.endtime, )

        if original_data:
            rows = 3
            fig  = plt.figure(figsize=(9,9))
            gs   = fig.add_gridspec(rows, 3)
            ax0  = fig.add_subplot(gs[0, :])
            ax0  = quick_plot(*original_data, verts=verts, ylabel='Amplitude', xlabel='Time', legend_loc='upper right', axis=ax0, xlim=[self.starttime, self.endtime])
            ax1  = fig.add_subplot(gs[rows-2, :], sharex=ax0)
            ax1  = quick_plot(*self.traces_input, verts=verts, ylabel='Amplitude', xlabel='Time', legend_loc='upper right', axis=ax1, xlim=[self.starttime, self.endtime], **kwargs)

        else:
            rows = 2
            fig  = plt.figure(figsize=(9,6))
            gs   = fig.add_gridspec(rows, 3)
            ax1  = fig.add_subplot(gs[rows-2, :])
            ax1  = quick_plot(*self.traces_input, verts=verts, ylabel='Amplitude', xlabel='Time', legend_loc='upper right', axis=ax1, xlim=[self.starttime, self.endtime], **kwargs)
       
        if title:
            fig.canvas.set_window_title('Ppol plot individual measurement: %s' % title)
        else:
            fig.canvas.set_window_title('Ppol plot individual measurement')
        fig.suptitle(title, fontsize=11)


        ### SMALL HACKS to make sure for small amplitudes everything works out (precisely: arrow head of angle indicators)
        factor     = 1
        factor_str = ''

        if data_Z is not None:
            BAZ         = BAZ_3D
            BAZ_err     = BAZ_err_3D
            BAZ_label   = 'BAZ_3D'

            if max( [max(abs(data_Z)), max(abs(data_1)), max(abs(data_2)), max(abs(data_R))] ) <= 1e-3:
                factor      = 1e9
                factor_str  = '%.0g' % (factor**-1)
                data_Z     *= factor
                data_1     *= factor
                data_2     *= factor
                data_R     *= factor

        else:
            BAZ         = BAZ_2D
            BAZ_err     = BAZ_err_2D
            BAZ_label   = 'BAZ_2D'

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
        yy      = 1/np.tan(BAZ*np.pi/180) * xx
        yy_be   = 1/np.tan((BAZ-BAZ_err)*np.pi/180) * xx
        yy_af   = 1/np.tan((BAZ+BAZ_err)*np.pi/180) * xx
        
        x       = maxi/2*np.sin( (BAZ-7)*np.pi/180 )
        y       = maxi/2*np.cos( (BAZ-7)*np.pi/180 )
        x2      = maxi/2*np.sin( (BAZ+2)*np.pi/180 )
        y2      = maxi/2*np.cos( (BAZ+2)*np.pi/180 )

        # Set-up
        ax2 = fig.add_subplot(gs[rows-1, 0])
        ax2.grid(ls='-.', lw=0.5)
        ax2.set(xlabel=label_2, ylabel=label_1)
        ax2.set_ylim([-maxi, maxi])
        ax2.set_xlim([-maxi, maxi])
        ax2.set_aspect('equal', adjustable='box')
        ax2.text(1, 0, '%s points' % len(data_1), ha='right', color='red', transform=ax2.transAxes, fontsize=6, bbox=dict(boxstyle='round,pad=0', fc='white', ec="white", lw=0))
        ax2.text(0, 1, factor_str, ha='center', color='black', transform=ax2.transAxes, fontsize=9, bbox=dict(boxstyle='round,pad=0', fc='white', ec="white", lw=0))
        
        ax2.spines['right'].set_color('none')
        ax2.spines['top'].set_color('none')
        ax2.spines['left'].set_color('none')
        ax2.spines['bottom'].set_color('none')  
        
        # Plot commands; data & Phi angle + erorrs
        if (BAZ+BAZ_err)//180 != BAZ//180.001 or BAZ == 0:
            ax2.fill_betweenx(yy_af, xx, facecolor='red', alpha=0.075)
            ax2.fill_betweenx(yy,    xx, facecolor='red', alpha=0.075)
        else:
            ax2.fill_between(xx, yy, yy_af, facecolor='red', alpha=0.075)       # error area of PHI
        if (BAZ-BAZ_err)//180 != BAZ//180.001 or BAZ == 0:
            ax2.fill_betweenx(yy_be, xx, facecolor='red', alpha=0.075)
            ax2.fill_betweenx(yy,    xx, facecolor='red', alpha=0.075)
        else:
            ax2.fill_between(xx, yy_be, yy, facecolor='red', alpha=0.075)       # error area of PHI
        ax2.plot([-maxi, maxi], [0, 0], 'k', lw=1)                              # centered horizontal line
        ax2.plot([0, 0], [-maxi, maxi], 'k', lw=1)                              # centered vertical line
        ax2.scatter(data_2, data_1, s=7, c=colours, zorder=3)                             # data
        ax2.plot( xx, yy, 'indianred', lw=1.5, zorder=4)                                  # PHI results
                
        # Angle arc + arrow head
        ax2.add_patch( Arc([0,0], maxi,  maxi, 90, -BAZ, 0, color='indianred', lw=1.5, zorder=5))
        a = FancyArrowPatch([x,y], [x2,y2], mutation_scale=20, lw=1.5, arrowstyle="-|>", color="indianred", zorder=6)
        ax2.add_artist(a)       

        # Legend
        scatter_proxy = mpl.lines.Line2D([0],[0], c="indianred", marker='>')
        ax2.legend([scatter_proxy], ['%s=(%.1f\u00B1%.1f)\u00b0' % (BAZ_label, BAZ, BAZ_err)], numpoints=1, loc='upper right', prop={'size': 8})


        ## plot axis (part) 2 and 3, if there's Z-data
        if data_Z is not None:

            # Variables needed for (beautiful) plotting
            maxi  = max( list(np.abs(data_R))+list(np.abs(data_Z)) ) * 1.05
            
            xx    = np.linspace(-maxi, maxi, 100)
            yy    = 1/np.tan( INC_3D*np.pi/180) * xx
            yy_be = 1/np.tan((INC_3D-INC_err_3D)*np.pi/180) * xx
            yy_af = 1/np.tan((INC_3D+INC_err_3D)*np.pi/180) * xx
            
            x     = maxi/2*np.sin( (INC_3D-7)*np.pi/180 )
            y     = maxi/2*np.cos( (INC_3D-7)*np.pi/180 )
            x2    = maxi/2*np.sin( (INC_3D+2)*np.pi/180 )
            y2    = maxi/2*np.cos( (INC_3D+2)*np.pi/180 )
            
            # Set-up
            ax3 = fig.add_subplot(gs[rows-1, 1])
            ax3.grid(ls='-.', lw=0.5, zorder=-1)
            ax3.set(xlabel=label_R, ylabel=label_Z)
            ax3.set_ylim([-maxi, maxi])
            ax3.set_xlim([-maxi, maxi])
            ax3.set_aspect('equal', adjustable='box')
            ax3.text(1, 0, '%s points' % len(data_2), horizontalalignment='right', color='red', transform=ax3.transAxes, fontsize=6, bbox=dict(boxstyle='round,pad=0', fc='white', ec="white", lw=0))
            ax3.text(0, 1, factor_str, ha='center', color='black', transform=ax3.transAxes, fontsize=9, bbox=dict(boxstyle='round,pad=0', fc='white', ec="white", lw=0))
        
            ax3.spines['right'].set_color('none')
            ax3.spines['top'].set_color('none')
            ax3.spines['left'].set_color('none')
            ax3.spines['bottom'].set_color('none')  
        
            # Plot commands; data & Phi angle + erorrs
            if (INC_3D+INC_err_3D)//180 != INC_3D//180.001 or INC_3D == 0:
                ax3.fill_betweenx(yy_af, xx, facecolor='red', alpha=0.075)
                ax3.fill_betweenx(yy,    xx, facecolor='red', alpha=0.075)
            else:
                ax3.fill_between(xx, yy, yy_af, facecolor='red', alpha=0.075)           # error area of INC
            if (INC_3D-INC_err_3D)//180 != INC_3D//180.001 or INC_3D == 0:
                ax3.fill_betweenx(yy_be, xx, facecolor='red', alpha=0.075)
                ax3.fill_betweenx(yy,    xx, facecolor='red', alpha=0.075)
            else:
                ax3.fill_between(xx, yy_be, yy, facecolor='red', alpha=0.075)           # error area of INC         
            ax3.plot([-maxi, maxi], [0, 0],  'k', lw=1)
            ax3.plot( [0, 0], [-maxi, maxi], 'k', lw=1)
            ax3.scatter(data_R, data_Z, s=5, c=colours, zorder=3)
            ax3.plot( xx, yy, 'indianred', lw=1.5, zorder=4)

            # Angle arc + arrow head
            ax3.add_patch( Arc([0,0], maxi,  maxi, 90, -INC_3D, 0, color='indianred', lw=1.5, zorder=5))
            a = FancyArrowPatch([x,y], [x2,y2], mutation_scale=20, lw=1.5, arrowstyle="-|>", color="indianred", zorder=6)
            ax3.add_artist(a)
            
            # Legend
            scatter_proxy = mpl.lines.Line2D([0],[0], c="indianred", marker='>')
            ax3.legend([scatter_proxy], ['INC_3D=(%.1f\u00B1%.1f)\u00b0' % (INC_3D, INC_err_3D)], loc='upper right', prop={'size': 8})



            ## 3-D plot
            # Variables needed for (beautiful) plotting
            mean_N   = np.average(data_1)
            mean_E   = np.average(data_2)
            mean_Z   = np.average(data_Z)
            maxi     = max( list(np.abs(data_Z))+list(np.abs(data_1))+list(np.abs(data_2)) ) * 1.05
            max_dist = np.sqrt(3*(maxi*2)**2)

            # Set-up
            ax4 = fig.add_subplot(gs[rows-1, 2], projection='3d')
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
def quick_plot(*y, 
    x            = None, 
    data_labels  = (), 
    ls           = '-',
    lw           = 1.5,
    lc           = None,
    win_title    = '', 
    title        = '',
    title_kwargs = {},
    xlabel       = 'Data points', 
    ylabel       = 'Amplitude',
    grid         = True,
    xscale       = None,
    yscale       = None,
    xinvert      = False, 
    yinvert      = False, 
    xlim         = None, 
    ylim         = None,
    date_format  = '%Y-%m-%d\n%H:%M:%S',
    verts        = None,
    vertsc       = None,
    horis        = None,
    horisc       = None,
    figsize      = (10,8),
    legend_loc   = 'best',
    axis         = False,
    outfile      = None, 
    show         = True, 
    keep_open    = False,
    **kwargs):

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


    ### FIXED VARIABLES
    fs_title  = 12
    fs_label  = 11
    fs_legend = 9



    ### Figure instance & settings
    if not axis:
        fig = plt.figure(num=title, figsize=figsize)
        fig.canvas.set_window_title(win_title)
        ax = fig.add_subplot(111)
    else:
        ax = axis

    ax.ticklabel_format(axis='y', style='sci', scilimits=(-3, 3))
    ax.set_title(title, fontsize=fs_title, **title_kwargs)
    ax.set_xlabel(xlabel, fontsize=fs_label)
    ax.set_ylabel(ylabel, fontsize=fs_label)
    if grid:
        ax.grid(ls='-.', lw=0.5)
    if xscale:
        ax.set_xscale(xscale)
    if yscale:
        ax.set_yscale(yscale)

    if xinvert:
        ax.invert_xaxis()
    if yinvert:
        ax.invert_yaxis()



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
        legend_loc = None
        labels     = [None for i in range(len(y))]



    ### Data plotting (empty data are not plotted)
    if isinstance(ls, str):
        ls = [ls]

    if isinstance(lw, (int, float)):
        lw = [lw]

    if isinstance(lc, str) or lc is None:
        lc = [lc]

    if x is None:
        x = [None for i in y]
    elif len(x)==len(y) and all(isinstance(i, (tuple, list, np.ndarray)) for i in x):
        pass
    else:
        x = [x for _ in range(len(y))]


    for j, data in enumerate(y):

        if len(ls)!=len(y):
            linestyle = ls[0]
        else:
            linestyle = ls[j]

        if len(lw)!=len(y):
            linewidth = lw[0]
        else:
            linewidth = lw[j]

        if len(lc)!=len(y):
            linecolor = lc[0]
        else:
            linecolor = lc[j]

        if isinstance(data, Trace):
            if x[j] is not None:
                xdata = x[j]
                if all(isinstance(ele, (UTCDateTime, datetime.datetime, str)) for ele in xdata):
                    xdata = [UTCDateTime(e).datetime for e in xdata]    # convert all to datetime.datetime objects
                    xdata = mdates.date2num(xdata)                      # convert all to matplotlib times (wouldn't need to, but for later it is need as ax.get_xlim() retrieves matplotlib times!)
                    ax.plot_date(xdata, data.data, ls=linestyle, lw=linewidth, c=linecolor, label=labels[j])
                    myFmt = mdates.DateFormatter(date_format)
                    ax.xaxis.set_major_formatter(myFmt)                # because we want dates in customised way, not matplotlib times
                else:   # normal array, relative times, matlab times, or POSIX timestamps
                    ax.plot(xdata, data, ls=linestyle, lw=linewidth, c=linecolor, label=labels[j])        
            else:
                xdata = data.times(type="matplotlib")
                ax.plot_date(xdata, data.data, ls=linestyle, marker=None, lw=linewidth, c=linecolor, label=labels[j])
                myFmt = mdates.DateFormatter(date_format)
                ax.xaxis.set_major_formatter(myFmt)

        elif isinstance(data, (Stream)):
            print(u'WARNING: Using plot() method of %s object and then return.' % type(data))
            plt.close()
            data.plot()     # use Stream, respectively, object's plot function
            return

        else:
            if x[j] is not None:
                xdata = x[j]
                if all(isinstance(ele, (UTCDateTime, datetime.datetime, str)) for ele in xdata):
                    xdata = [UTCDateTime(e).datetime for e in xdata]    # convert all to datetime.datetime objects
                    myFmt = mdates.DateFormatter(date_format)
                    ax.plot(xdata, data, ls=linestyle, lw=linewidth, c=linecolor, label=labels[j])
                    ax.xaxis.set_major_formatter(myFmt)                # because we want dates in customised way, not matplotlib times
                else:   # normal array, relative times, matlab times, or POSIX timestamps
                    ax.plot(xdata, data, ls=linestyle, lw=linewidth, c=linecolor, label=labels[j])
            else:
                xdata = np.arange(len(data))
                ax.plot(xdata, data, ls=linestyle, lw=linewidth, c=linecolor, label=labels[j])



    ### Limits of x- and y-axis
    if xlim is not None and not all(x is None for x in xlim):
        if all(isinstance(ele, (UTCDateTime, datetime.datetime, str)) for ele in xlim):
            xlim = [UTCDateTime(e).datetime for e in xlim]
            xlim = [mdates.date2num(e) for e in xlim]
        ax.set_xlim( xlim )
    

    if ylim is not None and not all(y is None for y in ylim):
        ax.set_ylim( ylim )
    else:           # make sure ylim is according to newly set xlim
        y_mins = []
        y_maxs = []
        for line in ax.lines:
            x = line.get_xdata()
            y = line.get_ydata()
            try:                                                                    # e.g. if one of the `y` is empty
                i = np.where( (x >= ax.get_xlim()[0]) & (x <= ax.get_xlim()[1]) )[0]    # all indexes of y_data according to xlim
                line_min = np.nanmin(y[i])                                          # get minimum y within all data according to xlim
                line_max = np.nanmax(y[i])                                          # get maximum y within all data according to xlim
                y_mins.append(line_min)
                y_maxs.append(line_max)
            except:
                continue

        try:
            y_min = np.nanmin(y_mins)
            y_max = np.nanmax(y_maxs)
            ylim = [y_min-0.025*np.abs(y_max-y_min), y_max+0.025*np.abs(y_max-y_min)] # give 2.5% margin that works both for pos. & neg. values 
            ax.set_ylim( ylim )
        except:
            pass
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    

    ### Vertical lines for data indications
    if isinstance(vertsc, str):
        colours = [vertsc]
    elif vertsc is None:
        colours = ['k']
    else:
        colours = vertsc
    
    if verts is not None: 

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
    if isinstance(horisc, str):
        colours = [horisc]
    elif horisc is None:
        colours = ['k']
    else:
        colours = horisc

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
    if legend_loc is not None:
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys(), loc=legend_loc, prop={'size': fs_legend})

    if axis:
        return ax

    plt.tight_layout()

    if outfile:
        plt.savefig(outfile)

    if show:
        plt.show()

    if not keep_open:
        plt.close()
def on_click(event, fig):
                axes = fig.axes
                axis = event.inaxes
                if axes is None or axis is None:
                    return                                                                        # case clicked not on axis but frame ..

                if event.dblclick:
                    if event.button == 1:                                                         # nothing right now
                        pass


                    elif event.button == 2:                                                       # middle-click (set ylim for each axis to maximum available ylim of all axes whilst centering data of each axis == equal scale)
                        xlim  = axis.get_xlim()
                        ylims = []

                        for axis in axes:

                            y_mins = []
                            y_maxs = []     

                            for line in axis.lines:
                                x = line.get_xdata()
                                if x[0]==x[-1]:                                                 # exclude vertical lines
                                    continue
                                y = line.get_ydata()
                                i = np.where( (x >= xlim[0]) & (x <= xlim[1]) )[0]              # all indexes of y_data according to xlim
                                if not i.any():
                                    continue                                                    # e.g. data gaps, that is, 2 lines within axes where one does not lie within chose xlims

                                line_min = y[i].min()                                           # get minimum y within all data according to xlim
                                line_max = y[i].max()                                           # get maximum y within all data according to xlim
                                y_mins.append(line_min)
                                y_maxs.append(line_max)
                            
                            try:        # e.g. empty axis by axis.twinx() or axis.twiny() if no data plotted in that axis ...
                                y_min = min(y_mins)
                                y_max = max(y_maxs)
                                ylims.append([y_min, y_max])
                            except:
                                continue

                        ylims        = np.array(ylims)
                        ylims_vals   = [ylim[1]-ylim[0] for ylim in ylims]
                        ylim_max_ind = np.argmax(ylims_vals)

                        for r, ylim in enumerate(ylims):
                            if r == ylim_max_ind:
                                pass
                            else:
                                ylim = [ylim[0]-ylims_vals[ylim_max_ind]/2+(ylim[1]-ylim[0])/2, ylim[1]+ylims_vals[ylim_max_ind]/2-(ylim[1]-ylim[0])/2]
                            ylim_set = [ylim[0]-0.025*np.abs(ylim[1]-ylim[0]), ylim[1]+0.025*np.abs(ylim[1]-ylim[0])] # give 2.5% margin that works both for pos. & neg. values
                            axes[r].set_ylim( ylim_set ) 
                        fig.canvas.draw()


                    elif event.button == 3:                                                     # right-click (for each axis individually, set ylim to respective data minimum and maximum + margin)
                        xlim   = axis.get_xlim()

                        for axis in axes:
                            y_mins = []
                            y_maxs = []
                            for line in axis.lines:
                                x = line.get_xdata()
                                if x[0]==x[-1]:                                                 # exclude vertical lines
                                    continue
                                y = line.get_ydata()
                                i = np.where( (x >= xlim[0]) & (x <= xlim[1]) )[0]              # all indexes of y_data according to xlim
                                if not i.any():
                                    continue                                                    # e.g. data gaps, that is, 2 lines within axes where one does not lie within chose xlims
                                line_min = y[i].min()                                           # get minimum y within all data according to xlim
                                line_max = y[i].max()                                           # get maximum y within all data according to xlim
                                y_mins.append(line_min)
                                y_maxs.append(line_max)

                            try:        # e.g. empty axis by axis.twinx() or axis.twiny() if no data plotted in that axis ...
                                y_min = min(y_mins)
                                y_max = max(y_maxs)
                                ylim  = [y_min-0.025*np.abs(y_max-y_min), y_max+0.025*np.abs(y_max-y_min)] # give 2.5% margin that works both for pos. & neg. values 
                                axis.set_ylim( ylim )
                            except:
                                continue
                        fig.canvas.draw()
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


### _ _ N A M E _ _ = = " _ _ M A I N _ _ "  
if __name__ == "__main__":
    pass