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


"""
Seisglitch script. This file is called when typing `seisglitch function path/to/config.yml`
into the terminal after installation of the package.

Handed arguments are parsed and the desired function is executed.
"""


#####  python modules import  #####
import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt


#####  obspy modules import  #####
from obspy import read


#####  seisglitch modules import  #####
import seisglitch
from seisglitch import detect, plot, remove, util

version = seisglitch.__version__
author  = seisglitch.__author__


### PARSING ARGUEMENTS
parser = argparse.ArgumentParser(
    description = """
                  SEISglitch Toolbox
                  """,
    prog        = 'seisglitch')

parser.add_argument('mode', 
        metavar = 'MODE',
        type    = str,
        help    = """
                  Choose which SEISglitch mode to run. 
                  Choices: 'detect', 'plot', 'remove', 'download', 'process', 'merge', or 'time'.
                  """,
        choices = ['detect', 'plot', 'remove', 'download', 'process', 'merge', 'time'])

parser.add_argument('config_file', 
        metavar = 'CONFIG_FILE',
        type    = str,
        help    = """
                  Choose config file for specific parameters.
                  """)
parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + version)

args   = vars( parser.parse_args() )
mode   = args['mode'].lower()
config = args['config_file']



### EXECUTE ACCORDING TO SPECIFIED PARAMETERS
params         = util.read_config(config)
waveform_files = params['waveform_files']
inventory_file = params['inventory_file']

if not inventory_file:
    inventory_file = 'IPGP'



if mode=='detect':

    detect_params = params['detect']['detector']
    detect.detect( *waveform_files, inventory_file=inventory_file, **detect_params )



elif mode=='remove':

    glitch_detector_files = params['remove']['glitch_detector_files']
    remove_params         = params['remove']['remover']

    remove.remove( *glitch_detector_files, waveform_files=waveform_files, inventory_file=inventory_file, **remove_params )



elif mode=='plot':
    
    glitch_detector_files = params['plot']['glitch_detector_files']
    plot_params           = params['plot']['plotter']

    plot.plot_glitch_remover(     *glitch_detector_files, **plot_params, **params['plot']['glitch_remove_plot'])
    plot.plot_glitch_detector(    *glitch_detector_files, **plot_params, **params['plot']['glitch_detector_plot'], waveform_files=waveform_files)
    plot.plot_glitch_overview(    *glitch_detector_files, **plot_params, **params['plot']['glitch_overview_plot'], waveform_files=waveform_files)
    plot.glitch_SOLoverLMST_plot( *glitch_detector_files, **plot_params, **params['plot']['glitch_SOLoverLMST_plot'])
    plot.plot_glitch_gutenberg(   *glitch_detector_files, **plot_params, **params['plot']['glitch_gutenberg_plot'] )
    plot.plot_glitch_ppol(                                **plot_params, **params['plot']['glitch_ppol_plot'],     waveform_files=waveform_files, inventory_file=inventory_file)

    if plot_params['show']:
        plt.show()
    plt.close('all')



elif mode=='download':

    kwargs = params['download']
    util.download_data(**kwargs)



elif mode=='process':

    process_params = params['process']
    util.process_data( *waveform_files, inventory_file=inventory_file, **process_params )



elif mode=='merge':

    glitch_detector_files = params['merge']['glitch_detector_files']
    if not glitch_detector_files:
        waveform_files        = params['waveform_files']        
        glitch_detector_files = [os.path.join( os.path.dirname(file), 'glitches_' + '.'.join(os.path.basename(file).split('.')[:-1]) + '.txt' ) for file in waveform_files]

    outfile       = params['merge']['outfile']
    multiples_out = params['merge']['multiples_out']

    util.merge_glitch_detector_files(glitch_detector_files, outfile, starttime_sort=True, multiples_out=multiples_out)



elif mode=='time':

    file    = params['time']['file']
    convert = params['time']['convert']
    print()

    if convert:
        print(u"Processing 'convert' option:")        

        time_converted = util.time_funnel(convert)
        if time_converted is None:
            print(u'Could not convert time: Please follow the following formats (h:m:s must not be fully stated):')
            print(u"For UTC:  'XXXX-XX-XXTxx:xx:xx.x'")
            print(u"For LMST: 'XXXMxx:xx:xx.x'")
            sys.exit()
        print(time_converted)

    elif file:
        print(u"Processing 'file' option:")        

        file_read = np.loadtxt(file, dtype=str)
        try:
            times = file_read[:,0]
        except IndexError:
            times = file_read
        times_converted = []

        for t, time in enumerate(times):
            time_converted = util.time_funnel(time)
            if time_converted is None:
                print(u"Time '%s' in input line %d is no valid. Skipped." % (time, t+1))
            times_converted.append( time_converted )

        file_read = np.c_[ file_read, times_converted ]
        np.savetxt(file, file_read, fmt='%25s')
        print(u'Finished. Conversions written to new column %s.' % len(file_read[1]))

    else:
        print(u'At least one option must be specified. Nothing done.')


else:
    pass