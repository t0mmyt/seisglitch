.. _remove:

seisglitch remove
=================

The remove function allows you to remove the glitches from the waveform files. 
As input the functions needs a glitch list as produced by the glitch detector module. 
Multiple glitch lists can also be passed.
To run the remover, activate your correct Python :ref:`environment <installation>`, fill the options
in the ``config.yml`` under the section ``remove``, then run from terminal:
::

    seisglitch remove path/to/config.yml

The glitch remover is now running according to your specifications.
The general working principle of the glitch removal is described in our paper, Section Glitch Removal (MPS).

In the following, the glitch detector's parameters (that you find in the ``config.yml``) are explained 
a bit more:

Cleaned data can take up more staorage, as fitting required to set important as glitch removal introduces float RAW values


* ``glitch_window_leftright``: in seconds, data window appended left and right to window of detected glitch (i.e., start and end time of glitches from the glitch detector files), buildung the overall window glitch fits are attempted. Whilst the glitch detection mostly delivers accurate glitch starts, meaning this parameter can be chosen smaller, for some glitches it can help increasing this parameter as the glitch onsets may not have been detected cleanly. 
* ``glitch_prepend_zeros``: in seconds, length of zeros prepended to glitch model. This can help improving fits. This time is added to the data window of 'glitch_window_leftright'. A good first try is 1 second.
* ``glitch_interpolation_samples``: integer. This is the option that can help removing glitches that do not well fit the acceleration step-model. To do so, set to e.g. 100. Glitches are then fitted with accleration steps with maximum 100 samples rise time (in steps of 10 samples). If not 0, algorithm is significantly slower! 
* ``glitch_subsample_factor``: integer, determines how many times beween two samples a glitch shall be modeled. For glitches 1 is typically OK. Default is one.
* ``spike_fit``: If True, attempt to fit spike after glitch fit. If True, removal takes longer.
* ``spike_subsample_factor``: integer, determines how many times beween two samples a glitch spike shall be modeled. For spikes 5 or higher may significantly improve their fits but is slower.
* ``spike_fit_samples_leftright``: samples left and right around fitted glitch onset where it is attempted to fit spike. Larger is slower. Default is 7.
* ``var_reduction_spike``: in %. Minimum spike variance reduction to be achieved if spike fit shall be removed. 2% is default.
* ``var_reduction_total``: in %. Minimum total variance reduction to be achieved if fit (either only glitch or glitch+spike) shall be removed. 80% is default.
* ``show_fit``: Attention, if True, an interactive plot is shown for each attempted fit on each component!
* ``store_glitches``: If True, also the corrections that were subtracted from the data are saved to file.
* ``plot_removal_statistic``: If True, two interactive plots are shown summarising overall statistics of deglitching. These plots are not perfectionized though.