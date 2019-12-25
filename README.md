# orbit-determination

The OD Code:

Takes in information from an input file in the following format:
YYYY MM DD HH:MM:SS HH MM SS.SS DEG ARCMIN ARCSEC SUNX SUNY SUNZ

One line should be used per observing session. Dates should be in UT and the sun vector should be in the equatorial coordinate system from Etscorn to the sun.

If a negative right ascension or declination is used, include a negative sign with all "sub-values" of RA and DEC (minutes and seconds).

The OD Code tests all feasible combinations from the input file and then prints out all possible outputs. It asks the user to pick one of these outputs to save to a result file which the Monte Carlo program reads from.

A sample result file has already been provided in the folder (1992JE Results) which includes a section for the user to manually input astrometric uncertainties.

Two input files have been provided (and one sample result file for the Monte Carlo). "1992JE Input" uses team 5's observations of 1992 JE which has 3 observations while OD Test Input has 5 observations that were randomly assembled from JPL data.

All parameters that would be of interest to change are located at the top of the OD and Monte Carlo codes. Many functions needed for both are included in the ODFunctions file.

The Monte Carlo file will save histograms of the plots of the orbital elements to the folder. OrbitViz is a visualization of the orbit of 1992JE. The standard deviations and medians calculated from this are not saved to a file.
