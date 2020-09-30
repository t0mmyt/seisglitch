.. _detect:

seisglitch detect
=================

The detect function allows detect the glitches in the VBB and SP waveform files. 
This is a necessary step if you want remove glitches (and spikes), as the removal function
takes the output of the glitch detection as input. 
The glitch detector is also responsible for the calculation of the glitch azimuths and incidence angles 
needed for certain plots in the :ref:`plot function <plot>`.

**As a rule of thumb; whole glitch detection may take around 2 minutes per 24 hours of 3-D data. DOWNSAMPLING DATA**

To use the glitch detector, activate your correct Python :ref:`environment <installation>`, fill the ``detect`` options
in the ``config.yml``, then run from terminal:
::

    seisglitch detect path/to/config.yml

The glitch detector is now running according to your specifications.
The general working principle of the glitch detector is described in `Scholz et al`_, (2020_), Section Glitch Detection (MPS method).

In the following, the glitch detector's parameters (that you find in the ``config.yml``) are detailed:


* ``taper_length_per_side``: in %, i.e., 0.05 means 5% for each side. On tapered data bits no glitch detection is done.
* ``pre_filt``: used to deconvolute instrument response.
* ``water_level``: used to deconvolute instrument response. If both specified (`pre_filt` as well), both are performed.
* ``ACCfilter``: dictionary specifying the type of filter_ used for the acceleration data. There is no obvious reason to change the default (bandpass: 0.001 - 0.1 Hz).
* ``threshold``: used to trigger glitches on the derivative of the filtered acceleration data (filtered jerk). Unit is therefore ms:sup:`-3`. Applied to positive and negative filtered jerk.
* ``plot_triggering``: If `True`, a plot is shown in between that shows all glitch triggers that passed the threshold condition. These candidates are not yet checked against their polarization.
* ``glitch_min_length``: In seconds. Defines the minimum length in which a new glitch cannot be detected once a glitch has been declared. 5 seconds is fine. This parameter allows to detect "poly-glitches".
* ``glitch_length``: Fixed, in seconds. For VBB 25 seconds is fine, for SP 50 seconds is fine.
* ``glitch_min_polarization``: Glitches have a high, linear polarization (circular=0 .. 1=linear). Therefore, this condition throws out glitch candidates that have too little linear polarization. A value of 0.9 may be a good start, higher is stricter and declares less glitches.

----

There are only two parameters that really affect the sensitivity of the glitch detection:
``threshold`` and ``glitch_min_polarization``.

* ``threshold``: the lower you choose it, the more often it triggers on the filtered jerk (remember, glitches are steps in acceleration, meaning their derivative should be like a delta-impulse whilst all other signals should not be delta-like). Obviously, at some stage, these triggers do not represent glitches anymore but signal / seismic noise. The threshold therefore should be chosen high enough to not trigger noise but low enough to detect small glitches, too (for smaller glitches, there are hundreds per Martian day in the VBB data). Keep further in mind that seismic amplitudes vary significantly during a Martian day (weather influence, amplitudes may change by a factor of 100 and more) so this complicates glitch detections. To circumvent, each glitch candidate is checked for its linear polarization (should be high for glitches). That is, the following parameter also has influence:
* ``glitch_min_polarization``: can be between 0 and 1, where 1 is full linear polarization and 0 circular polarization. Glitches, given they represent acceleration steps seen by the sensors, have a high, linear polarization. The lower you choose this parameter (e.g., 0.9), the more glitch candidates will pass and indeed be declared as glitch. In combination with ``threshold`` this should minimize potential false-positives, however, some certainly remain. Note the polarization analysis is performed on the original, gain corrected raw data rotated into the ZNE-system, i.e., no downsampling of the input data is done for this analysis. Table 1 summarizes sensible ranges for two parameters ``threshold`` and ``glitch_min_polarization``.


.. list-table:: Table 1: Sensible ranges for the two most important parameters influencing the performance of the glitch detection. `Strict` means less detections, `lax` means more detections. Stricter settings typically result in slightly faster run times.
   :widths: 25 25 50 50
   :header-rows: 1

   * - 
     - 
     - ``threshold``
     - ``glitch_min_polarization``
   * - **VBB**
     - strict
     - 2.0 x 10 :sup:`-9`
     - 0.98
   * - 
     - lax
     - 0.5 x 10 :sup:`-9`
     - 0.85 x 10 :sup:`-9`
   * - **SP**
     - strict
     - 
     - 
   * - 
     - lax
     - 
     - 


.. _filter: https://docs.obspy.org/packages/autogen/obspy.core.stream.Stream.filter.html
.. _Scholz et al: https://www.essoar.org/doi/10.1002/essoar.10503314.2
.. _2020: https://www.essoar.org/doi/10.1002/essoar.10503314.2