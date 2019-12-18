#!/usr/bin/env python
# -*- coding: utf-8 -*-


""""
Glitch detector for InSight's VBB and SP time series glitches.


#----------------------------------------------------------------------
#   Author:    John - Robert Scholz
#   Email:     john.robert.scholz@gmail.com
#   Date:      Dec 2019
#---------------------------------------------------------------------
"""


from setuptools import setup, find_packages

setup(name             = 'glitch_detector',
	  version          = '0.2',
	  description      = "Glitch detector for SEIS' VBB and SP seismic sensors.",
	  url              = 'https://pss-gitlab.math.univ-paris-diderot.fr/data-processing-wg/glitch-detector',
	  author           = 'John-Robert Scholz',
	  author_email     = 'john.robert.scholz@gmail.com',
	  license          = 'TBD',
	  packages         = find_packages(), 
	  install_requires = ['obspy'])
