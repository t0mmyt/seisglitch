.. _remove:

seisglitch remove
=================


The remove function allows to remove the glitches from the raw waveform files. 
As input the functions needs a glitch list as produced by the glitch :ref:`detector module <detect>`. 
Multiple output files of the detect function can also be passed (one file, multiple files, or multiple files :ref:`merged <merge>` into one)
To run the remover, activate your correct Python :ref:`environment <installation>`, fill the ``remove`` options
in the ``config.yml``, then run from terminal:
::

    seisglitch remove path/to/config.yml

The glitch remover is now running according to your specifications.
Cleaned data can take up more storage, as glitch fitting for subsequent removal required to introduce float RAW values.
The general working principle of the glitch removal is described in `Scholz et al`_, (2020_), Section Glitch Removal (MPS method).

In the following, the glitch detector's parameters (that you find in the ``config.yml``) are detailed:


* ``glitch_window_leftright``: in seconds, data window appended left and right to window of detected glitch (i.e., start and end time of glitches from the glitch detector files), building the overall window glitch fits are attempted. Whilst the glitch detection mostly delivers accurate glitch starts, meaning this parameter can be chosen smaller, for some glitches it can help increasing this parameter as the glitch onsets may not have been detected cleanly. 
* ``glitch_prepend_zeros``: in seconds, length of zeros prepended to glitch model. This can help improving fits. This time is added to the data window of 'glitch_window_leftright'. A good first try is 1 second.
* ``glitch_interpolation_samples``: integer. This is the option that can help removing glitches that do not well fit the acceleration step-model. To do so, set to e.g. 100. Glitches are then fitted with acceleration steps with maximum 100 samples rise time (in steps of 10 samples). If not 0, algorithm is significantly slower! 
* ``glitch_subsample_factor``: integer, determines how many times between two samples a glitch shall be modeled. For glitches 1 is typically OK. Default is one.
* ``spike_fit``: If True, attempt to fit spike after glitch fit. If True, removal takes longer.
* ``spike_subsample_factor``: integer, determines how many times between two samples a glitch spike shall be modeled. For spikes 5 or higher may significantly improve their fits but is slower.
* ``spike_fit_samples_leftright``: samples left and right around fitted glitch onset where it is attempted to fit spike. Larger is slower. Default is 7.
* ``var_reduction_spike``: in %. Minimum spike variance reduction to be achieved if spike fit shall be removed. 2% is default.
* ``var_reduction_total``: in %. Minimum total variance reduction to be achieved if fit (either only glitch or glitch+spike) shall be removed. 80% is default.
* ``show_fit``: Attention, if True, an interactive plot is shown for each attempted fit on each component!
* ``store_glitches``: If True, also the corrections that were subtracted from the data are saved to file.
* ``plot_removal_statistic``: If True, two interactive plots are shown summarizing overall statistics of deglitching. These plots are not perfected though.


.. _Scholz et al: https://www.essoar.org/doi/10.1002/essoar.10503314.2
.. _2020: https://www.essoar.org/doi/10.1002/essoar.10503314.2