.. _merge:

seisglitch merge
================


The merge function allows to merge multiple output files from the :ref:`detection function <detect>`.
The idea is that when you have run the glitch detector multiple times, e.g. for different data periods,
you can merge the results all into one file using this function. 
This can be handy especially for the :ref:`plot function <plot>`.

After having entered the `merge` options in the ``config.yml``,
run from the merge function from terminal like so:
::

    seisglitch merge path/to/config.yml