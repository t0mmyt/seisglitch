.. _detect:

seisglitch detect
=================

The detect function allows detect the glitches on the waveform files. 
This is a necessary step if you want remove glitches (and spikes), as the removal function
takes the output of the glitch detection as input. 
The glitch detector is also responsible for the calculation of glitch azimuths and incidence angles that
are needed for certain plots in the plot function.


To use the glitch detector, activate your correct Python :ref:`environment <installation>`, fill the options
in the ``config.yml`` under the section ``detect``, then run from terminal:
::

    seisglitch detect path/to/config.yml

The glitch detector is now running according to your specifications.
The general working principle of the glitch detector is described in our paper, Section Glitch Detection (MPS).

In the following, the glitch detector's parameters (that you find in the ``config.yml``) are explained 
a bit more:


* ``taper_length_per_side``: in %, i.e. 0.05 means 5% for each side. On tapered data no glitch detection is done.
* ``pre_filt``: used to deconvolute instrument response.
* ``water_level``: used to deconvolute instrument response. If both specified (pre_filt as well), both are performed.
* ``ACCfilter``: dictionary specifying the type of filter_ used for the acceleration data.
* ``threshold``: used to trigger glitches on derivative of filtered acceleration data. Unit is therefore m/s**3. Is applied to positive and negative values.
* ``plot_triggering``: If `True`, a plot is shown in between that shows all glitch triggers that passed the threshold condition. These candidates are not yet checked against their polarization.
* ``glitch_min_length``: In seconds. Defines the minimum length in which a new glitch cannot be detected once a glitch has been detected. 5 seconds is fine,
* ``glitch_length``: Fixed, in seconds. For VBB 25 seconds is fine, for SP 50 seconds is fine.
* ``glitch_min_polarization``: Glitches have a high, linear polarization (1=linear, 0=circular). Therefore, this condition throws out glitch candidates that have too little linear polarization. A value of 0.9 may be a good start, higher is stricter.

----

There are only two parameters that really affect the sensitivity of the glitch detection:
``threshold`` and ``glitch_min_polarization``.

- ``threshold``: the lower you choose it, the more often it triggers on the derivative of the filtered acceleration data (remember, glitches are steps in acceleration, meaning their derivative should be like a delta-impulse whilst all other signals should not be delta-like). Obviously, at some stage, these triggers do not represent glitches anymore but just signal / seismic noise. The threshold therefore should be chosen high enough to not start triggering noise, but low enough so you detect the right amount of glitches (for smaller glitches, there are hundreds per Martian day in the VBB data). Keep further in mind that seismic amplitudes vary significantly during a Martian day (weather influence, amplitudes may change by a factor of 100 and more) so this complicates things. To circumvent both effects and really detect glitches only, each candidate is checked for its linear polarization (should be high for glitches). That is, the following parameter also has influence:

- ``glitch_min_polarization``: can be between 0 and 1, where 1 means full linear polarization. The lower you choose it (e.g., 0.9), the more candidates will pass and be declared as glitch. In combination with the parameter ``threshold`` this should enable to minimize potential false-positives, however, some certainly remain. Note the polarization analysis is performed on the gain corrected raw data rotated into the ZNE-system.

In our experience, for the VBB seismometer, ``threshold`` may range from 2e-09 to 0.5e-9 m/s**3, i.e., the latter leading to more triggers 
(if the given filter is used: acausal, 3rd order, Butterworth 0.001-0.1 Hz). 
For SP, this parameter should generally be chosen a bit higher than for VBB, e.g. 5e-9 m/s**3.
Smaller thresholds take longer calculation times as more glitch candidates will be processed (e.g. checked for their polarization).
However, typically the whole glitch detection may take up only 2 minutes per 24 hours of data. The reason is internally all input data are decimated to 2 (or 2.5) SPS data prior to detection. 
Only the polarization analysis of glitch candidates is performed on the raw, unfiltered input data in the sampling rate provided! 
For the parameter ``glitch_min_polarization``, values may range 0.90 to 0.98, i.e. coarse to strict settings. That should work OK for both VBB and SP.


.. _filter: https://docs.obspy.org/packages/autogen/obspy.core.stream.Stream.filter.html
