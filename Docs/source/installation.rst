.. _installation:

Installation
============

I recommend using ``conda`` to install Python packages. 
It comes shipped with Anaconda_ (or Miniconda_ for a smaller footprint), a Python 
distribution that allows to quickly configure a specific Python environment (e.g. specific package versions etc..). 


1. Install Anaconda_ or Miniconda_.
2. Create a new environment and activate it (choose an environment name you prefer, here ``ppol`` is chosen)::

    conda create -n seisglitch python=3.7 obspy pandas pyyaml -c conda-forge
    conda activate seisglitch


3. Go to your prefered folder and download the repository::

    git clone https://gitlab.com/johnrobertscholz/seisglitch.git  
    cd seisglitch


4. Install ``seisglitch``::

    pip install .

Good news, everything is installed and ready to use!
Remember to activate the ``seisglitch`` environment each time you want to use it.

| 

----

| 

If you'd like to introduce your own changes to the code, you may want to replace step 4 by:

4. Install ``seisglitch``::

	pip install -e .

thus creating symbolic links to the module files, avoiding having to re-compile after each code alteration.

.. _Anaconda: https://www.anaconda.com/
.. _Miniconda: https://docs.conda.io/en/latest/miniconda.html
