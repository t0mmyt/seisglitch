.. _time:

seisglitch time
===============


The time function allows you to convert times both ways between UTC (earth time) and InSight's local mean solar time (LMST).
Note that InSight landed on Mars on November 26, 2018 (Sol 0, a sol is a Martian day with around 24h 40m).

The time formats should follow:

 | UTC: 2019-05-23T02:23:16
 | LMST: 173M02:58:53 (where 173 is Sol 173 with respect to the InSight landing)

In the ``config.yml`` you have two options:

  1. file - if you want to convert many times you can pass a text file with your times to be converted listed in a column, or 
  2. convert - if you just want to convert one time you can pass it with the result immediately printed. 

For specifics regarding each option see the config.yml. On the other hand, it is also 
possible to access the time conversion tool from within the Python interpreter (correct :ref:`environment <installation>`)):

.. code:: python

    from seisglitch.util import marstime

    time_instance = marstime(UTC='2019-05-23T02:23:16') # pass a string or time object
    LMST_string   = time_instance.LMST_string           # returns string
    print(LMST_string)

or similarly to convert LMST to UTC:

.. code:: python

    from seisglitch.util import marstime

    time_instance = marstime(LMST='173M02:58:53')       # pass a string
    UTC_string    = time_instance.UTC_string            # returns string
    UTC_time      = time_instance.UTC_time              # returns UTCDateTime object
    print(UTC_time)

If you would like to plot the current time, you could simply do:

.. code:: python

    from seisglitch.util import marstime

    time_instance = marstime()       # if no time passed, take time of now
    print(time_instance)             # prints UTC and LMST times