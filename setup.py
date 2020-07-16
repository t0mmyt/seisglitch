import os
import re
from numpy.distutils.core import setup


def find_version(*paths):
    fname = os.path.join(os.path.dirname(__file__), *paths)
    with open(fname) as fp:
        code = fp.read()
    match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", code, re.M)
    if match:
        return match.group(1)
    raise RuntimeError("Unable to find version string.")


setup(
    name             = 'seisglitch',
    version          = find_version('seisglitch', '__init__.py'),
    keywords         = ["glitches", "seismology", "planet mars"],
    description      = """
                       Toolbox to detect, plot, and remove glitches on SEIS' seismometers VBB and SP.
                       This code has been developed in the context of NASA's discovery mission 'InSight'.
                       """,
    author           = 'John-Robert Scholz',
    maintainer       = 'John-Robert Scholz',
    maintainer_email = 'john.robert.scholz@gmail.com',
    classifiers      = ['Development Status :: 3 - Alpha',
                       'License :: OSI Approved :: MIT License',
                       'Programming Language :: Python :: 3.7'],
    python_requires  =  '>=3.7',
    install_requires = ['numpy==1.16'],
    packages         = ['seisglitch'],
    scripts          = [os.path.join('Scripts',file) for file in os.listdir('Scripts/')],
    url              = 'https://pss-gitlab.math.univ-paris-diderot.fr/data-processing-wg/seisglitch'
    )