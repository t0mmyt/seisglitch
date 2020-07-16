.. _installation:

Installation
============

I recommend using ``conda`` to install Python packages. 
It comes shipped with Anaconda_ (or Miniconda_ for a smaller footprint), a Python 
distribution that allows to quickly configure a specific Python environment (e.g. specific package versions etc..). 


1. Install Anaconda_ or Miniconda_.
2. Open a terminal (Windows: Anaconda or Windows powershell), create a new environment and activate it (choose an environment name you prefer, here ``seisglitch`` is chosen)::

    conda create -n seisglitch python=3.7 numpy=1.16 obspy pyyaml -c conda-forge
    conda activate seisglitch


3. Go to to your prefered folder, download ``seisglitch`` (on Windows, you may need to download "Git" first)::

    cd your/prefered/folder
    git clone https://pss-gitlab.math.univ-paris-diderot.fr/data-processing-wg/seisglitch.git
    cd seisglitch


4. Install ``seisglitch``::

    pip install .

Good news, everything is installed and ready to use!
Remember to activate the ``seisglitch`` environment each time you want to use it (e.g. "conda activate seisglitch").
The package should be available from all paths in your system.

| 

----

| 

If you'd like to introduce your own changes to the code, you may want to replace step 4 by:

4. Install ``seisglitch``::

	pip install -e .

thus creating symbolic links to the module files, avoiding having to re-compile after each code alteration.

.. _Anaconda: https://docs.anaconda.com/anaconda/install/
.. _Miniconda: https://docs.conda.io/en/latest/miniconda.html
