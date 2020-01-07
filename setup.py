import os
import re
from numpy.distutils.core import setup


setup(
    name             = 'seisglitch',
    version          = '0.0.1',
    keywords         = ["InSight mission", "Planet Mars", "seismology", "glitch", "Obspy"],
    description      = "Detect, analyse, remove and evaluate removal for glitches in SEIS' (VBB & SP) time series data",
    author           = 'John-Robert Scholz',
    maintainer       = 'John-Robert Scholz',
    maintainer_email = 'john.robert.scholz@gmail.com',
    classifiers      = ['Development Status :: 3 - Alpha',
                       'License :: OSI Approved :: MIT License',
                       'Programming Language :: Python :: 3.6',
                       'Programming Language :: Python :: 3.7'],
    install_requires = ['obspy', 
                        'pandas', 
                        'pyyaml'],
    python_requires  =  '>=3.6',
    packages         = ['seisglitch'],
    scripts          = [os.path.join('Scripts',file) for file in os.listdir('Scripts/')],
    url              = 'https://gitlab.com/johnrobertscholz/ppol'
    )