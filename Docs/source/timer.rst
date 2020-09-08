.. _time:

seisglitch time
===============


The time function allows you to convert times both ways between UTC (earth time) and InSight's local mean solar time (LMST).
Note that InSight landed on Mars on November 26, 2018 (Sol 0, a sol is a Martian day with around 24h 40m).

The time formats should follow:
| UTC: 2019-05-23T02:23:16
| LMST: 173M02:58:53 (where 173 is Sol 173 with respect to the InSight landing)

In the ``config.yml`` you have two options; if you want to convert many times you can pass a text file, or if you just want to convert one time 
you can simply pass one with the result immediately printed. You can also access the time conversion tool from within the Python interpreter:

.. code:: python

    from seisglitch.util import marstime

    time_instance = marstime(UTC='2019-05-23T02:23:16') # you can pass a string or time object
    LMST_string   = time_instance.LMST_string           # returns a String object
    print(LMST_string)

or similarly to convert LMST to UTC:

.. code:: python

    from seisglitch.util import marstime

    time_instance = marstime(LMST='173M02:58:53')       # you can pass a string
    UTC_string    = time_instance.UTC_string            # returns a String object
    UTC_time      = time_instance.UTC_time              # returns a UTCDateTime object
    print(UTC_time)