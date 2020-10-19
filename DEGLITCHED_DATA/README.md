# Deglitched (corrected) Data


As outlined in Scholz et al. (2020, see below), four groups within the [InSight](https://mars.nasa.gov/insight/) project 
have developed code for the glitch and spike removal (=deglitching).
Here, we provide a selection of deglitched seismic events from Mars that the interested user may use or compare to her/his corrections.
Future updates of these files will have no official announcements, however, the code versions used for deglitching are 
indicated in the file names of the deglitched data. Version numbers differ across the four groups and have no temporal
implication, meaning version 4 of one group is not necessarily indicating more recent updates compared to version 1 of another group.
On the other hand, version numbers within a group allow you to conclude temporality, i.e., higher version numbers are more recent.


The four groups and their lead analysts are:
- MPS (John-Robert Scholz & Rudolf Widmer-Schnidrig, SEISglitch toolbox)
- ISAE (Baptiste Pinot & Raphaël Garcia, MATLAB)
- UCLA (Paul Davis, MATLAB)
- IPGP (Philippe Lognonné, MATLAB)


Note further that the timing of the waveform files (start & end times) may slightly differ across the different groups.
This is because some groups may have cut the waveforms differently and/or merged overlaps prior to processing whilst others did not. 
These differences are typically not significant, especially for total data lengths of less than 1 day (about 1 sol).

For more information on the Mars seismic catalogue curated by the Marsquake Service (MQS), see:  
https://www.seis-insight.eu/en/science/seis-products/mqs-catalogs

Please remember citing if you used the corrected data provided here:

    Scholz, J.‐R., Widmer‐Schnidrig, R., Davis, P., Lognonné, P., Pinot, B., Garcia, R. F., 
    et al. (2020). Detection, analysis and removal of glitches from InSight’s seismic 
    data from Mars. Earth and Space Science, 7, e2020EA001317‐T. 
    https://doi.org/10.1029/2020EA001317

and:

    Lognonné, P., W. B. Banerdt, W. T. Pike, D. Giardini, U. Christensen, R. F. Garcia, 
    T. Kawamura, et al. “Constraints on the Shallow Elastic and Anelastic Structure of
    Mars from InSight Seismic Data.” Nature Geoscience 13, no. 3 (March 2020), 213–20. 
    https://doi.org/10.1038/s41561-020-0536-y.

You can also help the Obspy developers by [citing](https://github.com/obspy/obspy/wiki#acknowledging) their work that the SEISglitch package heavily relies on.