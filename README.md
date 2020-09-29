# SEISglitch


**SEISglitch** is a Python package based on the [ObsPy](https://github.com/obspy/obspy/wiki) library.
Its purpose is to detect "glitches" in the seismic raw data of InSight's seismometers
VBB and SP, plot the detected glitches in different ways, and
remove them if possible. The SEISglitch package corresponds to the 'MPS' method described in the
peer-reviewed paper.

Full installation guide and documentation can be found at:  
https://seisglitch.readthedocs.io/en/latest/

If you are aiming for the deglitched data only, you can just go here:  
https://pss-gitlab.math.univ-paris-diderot.fr/data-processing-wg/seisglitch/tree/master/DEGLITCHED_DATA

If you are aiming for the MATLAB codes of the other methods ('ISAE', 'UCLA', or 'IPGP', if present) only, you can just go here:  
https://pss-gitlab.math.univ-paris-diderot.fr/data-processing-wg/seisglitch/tree/master/MATLAB_ALTERNATIVES