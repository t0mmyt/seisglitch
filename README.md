# SEISglitch


*SEISglitch* is a seismological Python package based on the [ObsPy](https://github.com/obspy/obspy/wiki) library.
Its purpose is to detect "glitches" - also referred to as "long-period data disturbances" - in the seismic raw data of [InSight's](https://mars.nasa.gov/insight/) seismometers
VBB and SP, plot the detected glitches in different ways, and
remove them if possible. The SEISglitch package corresponds to the 'MPS' method described in
[Scholz et al., 2020](https://www.essoar.org/doi/10.1002/essoar.10503314.2).

Full installation guide and documentation can be found at:  
https://seisglitch.readthedocs.io/en/latest/

If you are only aiming for the deglitched data provided for a selection of seismic events, go here:  
https://pss-gitlab.math.univ-paris-diderot.fr/data-processing-wg/seisglitch/tree/master/DEGLITCHED_DATA

If you are only aiming for the MATLAB codes of the other methods ('ISAE', 'UCLA', or 'IPGP', if present), go here:  
https://pss-gitlab.math.univ-paris-diderot.fr/data-processing-wg/seisglitch/tree/master/MATLAB_ALTERNATIVES