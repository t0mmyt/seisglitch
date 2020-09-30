import os
import re
from numpy.distutils.core import setup


### Small helper function
def find_version_author(*paths):
    fname = os.path.join(os.path.dirname(__file__), *paths)
    with open(fname, encoding='utf-8') as fp:
        code = fp.read()

    match1 = re.search(r"^__version__\s*=\s*['\"]([^'\"]*)['\"]", code, re.M)
    if match1:
        version = match1.group(1)
    else:
        version = ''

    match2 = re.search(r"^__author__\s*=\s*['\"]([^'\"]*)['\"]", code, re.M)
    if match2:
        author = match2.group(1)
    else:
        author = ''

    return (version, author)



### Setup
version, author = find_version_author('seisglitch', '__init__.py')

setup(
    name             = 'seisglitch',
    version          = version,
    keywords         = ["glitches", "data disturbances", "seismology", "planet Mars", "InSight mission", "VBB seismometer", "SP seismometer"],
    description      = """
                       Toolbox to detect, remove, and plot glitches on SEIS' seismometers VBB and SP.
                       This code has been developed in the context of NASA's Discovery mission 'InSight'.
                       """,
    author           = 'John-Robert Scholz',
    maintainer       = 'John-Robert Scholz',
    maintainer_email = 'john.robert.scholz@gmail.com',
    license          = 'MIT',
    classifiers      = ['Development Status :: 3 - Alpha',
                        'Topic :: Scientific/Engineering :: Physics',
                        'Programming Language :: Python :: 3.7',
                        'License :: OSI Approved :: MIT License',
                        'Operating System :: OS Independent'],
    python_requires  =  '>=3.7',
    install_requires = ['numpy==1.16'],
    packages         = ['seisglitch'],
    scripts          = [os.path.join('Scripts',file) for file in os.listdir('Scripts/')],
    url              = 'https://pss-gitlab.math.univ-paris-diderot.fr/data-processing-wg/seisglitch'
    )


### Only displayed anyway if `pip install .` is run with "-v" option
if version:
    print(u"Successfully installed seisglitch %s." % version)
else:
    print(u"Successfully installed seisglitch (no version specified).")

if author:
    print(u"Author: %s" % author)