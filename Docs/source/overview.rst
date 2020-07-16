.. _overview:

Overview
========

Once ``seisglitch`` is :ref:`installed <installation>`, all you have to do is:

* data preparation
* setup of ``config.yml`` file that came with the download
* execute the intended ``seisglitch`` function


.. _data_prep:

Data preparation
^^^^^^^^^^^^^^^^

Your waveform data should all comply with the following conditions:

* all waveform files must contain all three seismic components (i.e., "U", "V", "W" of either the VBB or SP seismometer)
* all waveform files must be readable as seismological files (e.g. MSEED, SAC, ...)




Setup config.yml
^^^^^^^^^^^^^^^^

This file – that you downloaded along with the ``seisglitch`` package – allows you to specify all needed options. 
Go ahead, open the config file and enter your settings – they should be straightforward.




Execute seisglitch functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There are three main ``seisglitch`` functions that you most likely want to use. 
Before you execute them, remember to activate your environment beforehand (see :ref:`installation`)!
Then, type one of the following commands into your terminal:
::

    seisglitch detect path/to/config.yml
    seisglitch plot path/to/config.yml
    seisglitch remove path/to/config.yml

You can see right away, you can have multiple ``config.yml`` files for multiple setups. 
Personally, I had one for the VBB and one for the SP seismometer.

There are furthermore additional ``seisglitch`` functions that shall ease data handling for your convenience.
Access them from your terminal like so:
::

    seisglitch download path/to/config.yml
    seisglitch decimate path/to/config.yml
    seisglitch merge path/to/config.yml
    seisglitch time path/to/config.yml

Each of these seven functions has a dedicated section in the ``config.yml``. 
Typically, each function only accesses the accordant section in the ``config.yml`` file, however,
there is one small exception; the 'merge' function will look for the waveform_files of the 'detect' function if
no glitch detector files to merge are specified. :)