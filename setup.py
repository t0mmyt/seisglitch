#!/usr/bin/env python
# -*- coding: utf-8 -*-


""""
Some python tools for the InSight mars mission.

:copyright:
    Martin van Driel (Martin@vanDriel.de), 2018
    Simon Stähler (mail@simonstaehler.com)
"""


from setuptools import setup, find_packages

setup(name         = 'glitch_detector_v1',
	  version      = '0.1',
	  description  = "Glitch detector for SEIS' VBB and SP seismic sensors.",
	  url          = 'https://pss-gitlab.math.univ-paris-diderot.fr/data-processing-wg/glitch-detector',
	  author       = 'John-Robert Scholz',
	  author_email = 'john.robert.scholz@gmail.com',
	  license      = 'None',
	  packages     = find_packages(), install_requires=['obspy'])
