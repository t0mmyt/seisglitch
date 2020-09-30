.. _plot:

seisglitch plot
===============


The plot function allows to plot the :ref:`detected glitches <detect>` in different ways.
Most of the plots available correspond to the plots shown in the `Scholz et al`_, (2020_) and its Supplementary Information 2.
The idea is to specify the output files of the detect function (one file, multiple files, or multiple files :ref:`merged <merge>` into one),
subsequently select those glitches you are interested in (e.g. horizontal polarizations only), and finally specify 
which plots shall be generated (see `run` option for each plot).
More information are given in the config.yml file.

After having entered the `plot` options in the ``config.yml``,
run the plot function from terminal like so:
::

    seisglitch plot path/to/config.yml


The following plots show an example with list of glitches on the VBB seismometer detected for 2019
(this list is not distributed as part of the package, however, a similar list is published as Supplementarty Information 1).
The following two plots were produced by running the ``glitch_overview_plot`` option.


.. image:: _static/overview1.png
    :alt: glitch azimuth and incidence angles

.. image:: _static/overview2.png
    :alt: glitch sphre and histogram


.. _Scholz et al: https://www.essoar.org/doi/10.1002/essoar.10503314.2
.. _2020: https://www.essoar.org/doi/10.1002/essoar.10503314.2