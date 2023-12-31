.. _download:

seisglitch download
===================


The download function allows to download data depending on the channel IDs. 
The options in the ``config.yml`` are straightforward. This function should
also work to retrieve seismic data not related to the InSight mission.
More information are given in the config.yml file.

After having entered the `download` options in the ``config.yml``,
run the download function from terminal like so:
::

    seisglitch download path/to/config.yml


Note: All details on channel names and the SEIS instrument itself can be found in:
::

    Lognonné, P., W. B. Banerdt, D. Giardini, W. T. Pike, U. Christensen, P. Laudet, 
    S. de Raucourt, et al. “SEIS: Insight’s Seismic Experiment for Internal Structure 
    of Mars.” Space Science Reviews 215, no. 1 (January 2019). 
    https://doi.org/10.1007/s11214-018-0574-6.