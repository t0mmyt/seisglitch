.. _installation:

Installation
============

I recommend using ``conda`` to install Python packages. 
It comes shipped with Anaconda_ (or Miniconda_ for a smaller footprint), a Python 
distribution that allows to quickly configure a specific Python environment (e.g. specific package versions etc..). 


1. Install Anaconda_ or Miniconda_.
2. Create a new environment and activate it (choose an environment name you prefer, here ``seisglitch`` is chosen)::

    cd your/prefered/folder
    conda create -n seisglitch python=3.7 numpy=1.16 obspy pyyaml -c conda-forge
    conda activate seisglitch


3. Go to your prefered folder and download the repository::

    git clone pss-gitlab.math.univ-paris-diderot.fr:data-processing-wg/seisglitch.git
    cd seisglitch


4. Install ``seisglitch``::

    pip install .

Good news, everything is installed and ready to use!
Remember to activate the ``seisglitch`` environment each time you want to use it.
The package should be available from all paths.

| 

----

| 

If you'd like to introduce your own changes to the code, you may want to replace step 4 by:

4. Install ``seisglitch``::

	pip install -e .

thus creating symbolic links to the module files, avoiding having to re-compile after each code alteration.

.. _Anaconda: https://www.anaconda.com/
.. _Miniconda: https://docs.conda.io/en/latest/miniconda.html
