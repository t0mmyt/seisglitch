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
....
"""



#####  python modules import  #####
import sys
import argparse


#####  ppol package imports  #####
from ppol import run


#####  _ _ N A M E _ _ = = " _ _ M A I N _ _ "  #####
if __name__ == "__main__":


    ### PARSING ARGUEMENTS
    parser = argparse.ArgumentParser(
        description = """
                      PPOL PROGRAM
                      """)
    parser.add_argument('mode', 
            metavar = 'MODE',
            type    = str,
            help    = """
                      Choose which Ppol function to run. 
                      Choices: 'main' or 'eval'.
                      """,
            choices = ['main', 'eval'])
    parser.add_argument('config_file', 
            metavar = 'CONFIG_FILE',
               type = str,
               help = """
                      Path to config yaml file.
                      Needed for both `MODE main` and `MODE eval`.
                      """)
    parser.add_argument('-r','--ppol_results_file_all', 
               type = str,
               help = """
                      Path to results_all file created by `MODE main`.
                      Needed only for `MODE eval`.
                      """,
            default = None,
           required = False)

    args                  = vars( parser.parse_args() )
    mode                  = args['mode'].lower()
    config_file           = args['config_file']
    ppol_results_file_all = args['ppol_results_file_all']



    ### EXECUTE ACCORDING TO SPECIFIED PARAMETERS
    if mode=='main':
        
        run.ppol_main( config_file )


    elif mode=='eval':

        if not ppol_results_file_all:
            print(u"For Ppol mode `eval` need variable '--ppol_results_file_all' (short '-r').")
            print(u"Check Ppol help via '-h' or '--help'." )
            sys.exit()

        run.ppol_eval( config_file, ppol_results_file_all ) 