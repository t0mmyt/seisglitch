.. _overview:

Overview
========

Once ``seisglitch`` is :ref:`installed <installation>`, all you have to do is:

* prepare the data
* setup of ``config.yml`` file that came with the download
* execute the intended ``seisglitch`` function


.. _data_prep:

Data preparation
^^^^^^^^^^^^^^^^

Your waveform data should all comply with the following conditions:

* each waveform file must contain all three seismic components (i.e., "U", "V", "W" of either the VBB or SP seismometer). That is, if you have individual files for each component you must merge them into **one** file before. For data retrieval, see e.g. :ref:`download`.
* each waveform file must be readable as seismological file (e.g. MSEED, SAC, ...)




Setup config.yml
^^^^^^^^^^^^^^^^

This file – that you downloaded along with the ``seisglitch`` package – allows you to specify all needed options. 
Go ahead, open the ``config.yml`` and enter your settings. They should be straightforward.




Execute seisglitch functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There are three main ``seisglitch`` functions that you most likely want to use. 
Before you execute them, remember to activate your environment beforehand (see :ref:`installation`)!
Then, type one of the following commands into your terminal, depending on what you want to do:
::

    seisglitch detect path/to/config.yml
    seisglitch remove path/to/config.yml
    seisglitch plot path/to/config.yml

You can see right away, you can have multiple ``config.yml`` files accounting for multiple setups. 
Personally, I had one for the VBB and one for the SP seismometer.

There are furthermore additional ``seisglitch`` functions that shall ease data handling for your convenience.
Access them from your terminal like so:
::

    seisglitch download path/to/config.yml
    seisglitch process path/to/config.yml
    seisglitch merge path/to/config.yml
    seisglitch time path/to/config.yml

Each of these seven functions has a dedicated section in the ``config.yml``.