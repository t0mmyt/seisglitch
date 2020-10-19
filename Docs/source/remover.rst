.. _remove:

seisglitch remove
=================


The remove function allows to remove the glitches from the raw waveform files. 
As input it needs a glitch list as produced by the glitch :ref:`detector module <detect>`. 
Multiple output files can also be passed (one file, multiple files, or multiple files :ref:`merged <merge>` into one).
To run the remover, activate your correct Python :ref:`environment <installation>`, fill the ``remove`` options
in the ``config.yml``, then run from terminal:
::

    seisglitch remove path/to/config.yml

The glitch remover is now running according to your specifications. Note that each listed glitch is attempted to be removed from all three components,
regardless on which component this glitch was originally detected. This mitigates propagating missing glitch detections into the corrected waveform files.
Corrected data can take up more storage, as glitch fitting for removal requires to introduce float raw values.
The general working principle of the glitch removal is described in `Scholz et al`_, (2020_), Section Glitch Removal (MPS method).

In the following, the glitch detector's parameters (that you find in the ``config.yml``) are detailed:


* ``glitch_window_leftright``: in seconds, extension of glitch window left *and* right, i.e., prepended to the glitch start and appended to the glitch end, both given in the glitch detector files (the length of the modeled glitch is detected glitch end minus glitch start time). The resulting window is the window in which, for each sample, fits between the modeled glitch and the data are attempted. Whilst the glitch detection mostly delivers accurate glitch starts, meaning this parameter can be chosen smaller, it can help to increase it as the glitch onsets may not have been detected cleanly in all instances. 
* ``glitch_prepend_zeros``: in seconds, length of zeros prepended to glitch model. This can help improving the fits especially for low signal-to-noise ratios. This time is added to the data window of `glitch_window_leftright`. A good value is 1 second.
* ``glitch_interpolation_samples``: integer. This is the option that can help removing glitches that do not well fit the acceleration step-model. Set to e.g. "100". Glitches are then fitted with acceleration steps with maximum 100 samples rise time (in steps of 10 samples). If the parameter is "0" (true acceleration step of zero rise-time), the algorithm is *significantly* faster but some glitches onsets are not well removed.
* ``glitch_subsample_factor``: integer, determines how many times between two samples a glitch shall be modeled (true for all sampling periods). For glitches 1 is typically OK. Higher values slow computations.
* ``spike_fit``: If `True`, attempt to fit spike after glitch fit. If `True`, removal takes longer. If `False`, no spikes are attempted to be removed.
* ``spike_subsample_factor``: integer, determines how many times between two samples a glitch spike shall be modeled (true for all sampling periods). For spikes 5 or higher may significantly improve their fits.
* ``spike_fit_samples_leftright``: samples left *and* right around fitted glitch onset where it is attempted to fit the modeled spike. Larger is slower. Default is 7 samples.
* ``var_reduction_spike``: in %. Minimum spike variance reduction to be achieved if spike fit shall be removed. 2% is default.
* ``var_reduction_total``: in %. Minimum total variance reduction to be achieved if fit (either only glitch or glitch+spike) shall be removed. 80% is default.
* ``show_fit``: Attention, if `True`, an interactive plot is shown for each attempted fit on each component!
* ``store_glitches``: If `True`, also the corrections that were subtracted from the data are saved to file, i.e., "glitch(+spike) times series".
* ``plot_removal_statistic``: If `True`, two interactive plots are shown summarizing overall statistics of glitch(+spike) removal. These plots are not perfected though.


.. _Scholz et al: https://doi.org/10.1029/2020EA001317
.. _2020: https://doi.org/10.1029/2020EA001317