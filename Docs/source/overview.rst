.. _overview:

Overview
========

Once ``seisglitch`` is :ref:`installed <installation>`, all you have to do is:

* prepare your data
* specify your settings in the ``config.yml`` file, and 
* execute the intended ``seisglitch`` module.


----



Data preparation
^^^^^^^^^^^^^^^^

Your waveform data should all comply with the following conditions:

* all waveform files must contain all three seismic components (i.e., "U", "V", "W" of either VBB or SP)
* all waveform files must contain no other components
* all waveform files must be readable as seismological files (e.g. MSEED, SAC, ...)
* all waveform files must contain no gaps (see: seisglitch.util.pierce_stream)




Setup config.yml
^^^^^^^^^^^^^^^^

This file – that you downloaded along with the ``seisglitch`` package – allows you to specify all needed options. 
Go ahead, open the config file and enter your settings – they should be straightforward. 
For a detailed parameter discussion, see :ref:`config_file`.



----



**Code execution**
^^^^^^^^^^^^^^^^^^

In general, ``seisglitch`` is intended to be run from as terminal command.
Don't forget to activate your correct environment beforhand, see :ref:`installation`.
::