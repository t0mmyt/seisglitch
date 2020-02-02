.. _detect:

seisglitch.detect
=================

Activate your correct Python :ref:`environment <installation>`, fill the options
in the ``config.yml`` under the section ``detect``, then run from terminal:
::

    seisglitch detect path/to/config.yml

The glitch detector is now running according to your specifications.

Generally, the glitch detector performs the following steps internally:

1. cut RAW waveforms into 12 hours parts and deconvolute instrument response to acceleration.
2. filter acceleration data according to specified filter ``ACCfilter`` (glitches are steps in acceleration).
3. calculate derivative of filtered accerleration data (glitches will be close to a delta-function).
4. detect delta-functions via threshold that equals moving average absolute peak height * ``average_peak_height_times``.
5. once a glitch is triggered, verify if the same glitch is found or not on other components (parameter: ``glitch_min_dist``).
6. for each detected glitch, check its polarization that must be equal or larger to ``glitch_min_polarization``, making sure it's indeed a glitch (glitches have a high, linear polarization).


In the following, the glitch detector's parameters (that you find in the ``config.yml``) are explained in more detail:


* ``waveform_files``: all RAW waveform files (components U, V, W) where the glitch detector shall run on. Note, for each file, the glitch detector creates one output file. (For a helping tool to create the needed files as :ref:`specified <data_prep>`, see: ``seisglitch.util.merge_glitch_detector_files``).
* ``inventory_file``: Need for removal of instrument response. Specify either path to (downloaded) file or online data centre, e.g. 'IRIS' or 'IPGP'.
* ``taper_length_per_side``: length of taper legth per side in fraction of per cent, i.e. 0.05 means 5% for each side.
* ``pre_filt``: pre-filter used to deconvolute instrument response.
* ``water_level``: water level used to deconvolute instrument response.
* ``ACCfilter``: dictionary containing the type of filter_, its options, and a discriptive string used for output.
* ``window_length_minutes``: length of moving window to determine peaks in derivative of fitered accerleration data.
* ``average_peak_height_times``: the average absolute peak height in the moving window times this parameter is the threshold for glitch detection.
* ``show_triggering``: If `True`, a plot is shown in between that shows all glitch triggerings. Note the terminal infomration when doing so.
* ``glitch_min_dist``: In seconds. Defines the minimum distance in which a new glitch cannot be detected once a glitch has been detected.
* ``glitch_length``: Fixed, in seconds.
* ``glitch_min_polarization``: Glitches have a high, linear polarization (1=linear, 0=circular). Therefore, this condition throws out glitch triggerings that may not be glitches. 
* ``plot_individual``: If `True`, for each glitch (!) four plots are stored in the current working directory showing the polarization analysis. One plot for gain corrected RAW waveforms, and for each waveforms corrected to displacement, velocity and accerleration. Note, this option siginficantly slows down the detection.

.. _filter: https://docs.obspy.org/packages/autogen/obspy.core.stream.Stream.filter.html