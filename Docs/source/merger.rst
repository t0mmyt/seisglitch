.. _merge:

seisglitch merge
================


The merge function allows to merge multiple output files from the :ref:`detection function <detect>`.
The idea is when you run the glitch detector multiple times, e.g. for different data periods,
you can merge the results all into one file using this function. 
The merge function also allows to remove duplicated detections e.g. caused by merging glitch detector files whose analysis periods overlap.
Merging glitch detector files can be handy especially for the :ref:`plot function <plot>` however it is not necessary.

After having entered the `merge` options in the ``config.yml``,
run the merge function from terminal like so:
::

    seisglitch merge path/to/config.yml