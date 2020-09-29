.. _process:

seisglitch process
==================


The process function allows to apply some automated, basic processing to the specified waveforms.
Currently implemented is:

* removing gain from seismic traces
* removing instrument response of seismic traces
* rotate non-orthogonal UVW components into ZNE-system
* decimate data with FIR-coefficients

The options in the ``config.yml`` are straightforward. More informatione are given there.