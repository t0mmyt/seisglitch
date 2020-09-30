.. _plot:

seisglitch plot
===============


The plot function allows to plot the detected glitches (see :ref:`detect`) in different ways.
Most of the plots available correspond to the plots shown in the paper and its Supplementary Information 2.
The idea is to specify the output files of the detect function (one file, multiple files, or multiple files :ref:`merged <merge>` inti one),
subsequently select those glitches you are interested in, and finally specify which plots shall be generated (see `run` option for each plot).
More information are given in the config.yml file.

After having entered the `plot` options in the ``config.yml``,
run the plot function from terminal like so:
::

    seisglitch plot path/to/config.yml


The following plots show an example with 2019 list of glitches detected on the VBB seismometer
(this list is not distributed as part of the package, however, a similar list is published as Supplementarty Information 1).
The plots were produced by running the ``glitch_overview_plot`` option.


.. image:: overview1.png
    :alt: glitch azimuth and incidence angles

.. image:: overview2.png
    :alt: glitch sphre and histogram