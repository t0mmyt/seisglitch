.. _installation:

Installation
============

We recommend using ``conda`` to install Python packages. 
It comes shipped with Anaconda_ (or Miniconda_ for a smaller footprint), a Python 
distribution that allows to quickly configure a specific Python environment (e.g. specific package versions etc..). 
If you are interested in the whole Python toolbox (corresponding to method 'MPS', see `Scholz et al`_., 2020_), carry out step 1-4.
If you are aiming for the MATLAB files only (corresponding to methods 'ISAE', 'UCLA', or 'IPGP') or the deglitched data
provided along with this package for a selection of seismic events, you just need to do step 3. 


1. Install Anaconda_ or Miniconda_.
2. Open a terminal (Windows: Anaconda powershell), create a new environment and activate it (choose an environment name you prefer, here `seisglitch` is chosen)::

    conda create -n seisglitch python=3.7 numpy=1.16 obspy pyyaml -c conda-forge
    conda activate seisglitch


3. Go to to your prefered parent folder (no need to create an extra package folder, this will be done anyway), 
download ``seisglitch`` (on Windows, you may need to install "git" first)::

    cd your/prefered/folder
    git clone https://pss-gitlab.math.univ-paris-diderot.fr/data-processing-wg/seisglitch.git


4. Go in the seisglitch folder and install it via ``pip`` (on Windows, you may need to install "pip" first)::

    cd seisglitch
    pip install .

Good news, everything is installed and ready to use!
Remember to activate the seisglitch environment each time you want to use it (e.g. "conda activate seisglitch").
The package should be available from all paths in your system. 
You can check the seisglitch version you have installed via::

    seisglitch -v

**The most recent version is 1.0.0.**

----

If you'd like to introduce your own changes to the code, you may want to replace step 4 by:

4. Install seisglitch::

	pip install -e .

thus creating symbolic links to the module files, avoiding having to re-compile after each code alteration.

.. _Anaconda: https://docs.anaconda.com/anaconda/install/
.. _Miniconda: https://docs.conda.io/en/latest/miniconda.html
.. _Scholz et al: https://www.essoar.org/doi/10.1002/essoar.10503314.2
.. _2020: https://www.essoar.org/doi/10.1002/essoar.10503314.2

----

Update / Upgrade
================

If you already have ``seisglitch`` installed and just would like to get the newest version, simply:


Open a terminal (Windows: Anaconda powershell) and type::

    conda activate seisglitch
    cd your/seisglitch/folder  # see Installation
    git pull
    pip install . (or pip install -e .)

If you are aiming only for the potentially updated MATLAB files or updated deglitched data, you just need to do::

    cd your/seisglitch/folder  # see Installation
    git pull