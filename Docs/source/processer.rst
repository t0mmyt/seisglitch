.. _process:

seisglitch process
==================


The process function allows to apply some automated, basic processing to the targeted waveforms.
Currently implemented is:

* removing gain from seismic raw data
* removing instrument response from seismic raw data
* rotate non-orthogonal UVW components into ZNE-system
* decimate data with FIR-coefficients

After having entered the `process` options in the ``config.yml``,
run the process function from terminal like so:
::

    seisglitch process path/to/config.yml