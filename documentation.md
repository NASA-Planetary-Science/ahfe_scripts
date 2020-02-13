# Apollo Heat Flow Experiment Concatenation Scripts

## introduction

These scripts produce ASCII tables containing corrected, reduced, and
concatenated versions of all available calibrated data from the Apollo 15 and
17 Heat Flow Experiments. This file serves as minimal documentation for the
scripts as such. For more detailed background on the AHFE and its data, please
refer to documentation included in PDS4 bundle
urn:nasa:pds:a15_17_hfe_concatenated (not yet available in the PDS but
permanently mirrored at https://github.com/MillionConcepts/ahfe).

These scripts are useless without their highly-specific source data, which are
also included in this repository. Sources for these data are listed at the
bottom of this document.

They are Python 3 scripts. Dependencies beyond the Python standard library are
``numpy,`` ``astropy,`` and ``pandas.`` They were prepared using the Anaconda
distribution of Python, and are most safely used in Python environments
prepared by that distribution.

## usage 

Running hfe_cleaner.py in the root directory of this repository will produce a
directory named 'data_local' with subdirectories named 'clean,' 'depth,' and
'split' that contain differently-formatted versions of the AHFE data (see the
'formats' section below for details).

hfe_corrections.py contains flags and corrections for individual source files:
locations of impulsive outliers, orders of magnitude for numerically damaged
data, and so on. 

hfe_utils.py contains various utilities for converting and organizing AHFE
data.

## format

These scripts output data in three 'flavors' designed for different use cases.
These sets may also be viewed as templates for the assembly of various
versions of the AHFE data. For instance, users who wish to produce a set like
``depth`` but which retains impulsive outliers might wish to examine the
``forbidden_flags`` variable in hfe_utils.py.

### ``clean``

This set follows the basic table and number format of the 1971-1974 Langseth
et al. PDS-archived data for maximum compatibility with existing AHFE
workflows. It provides time, T, and dT for files 1, 2, 4, and 5; and time,
reference thermometer temperature, thermocouple temperatures, and heater state
for file 3. It also adds one or two additional columns per file: ``flags`` and
sometimes also ``dT_corr.`` ``flags`` is a bitmask marking special features of
some data points, such as the presence of a missing data constant or an
artifact. See section 'Bitmask Values' below for an index to these values.
``dT_corr`` contains dT values after correction for numerical errors that
appear in some of the Lamont data. Please refer to the bundle documentation
for descriptions of these errors and artifacts.

It takes data from the sets produced by Nagihara et al. and converts them to
the table and number format of the Langseth et al. data, organizing them by
file but keeping 1975, 1976, and 1977 data separate from one another and from
the Langseth et al. data. It also converts time units in these files from DOY
UTC to mission epoch  (milliseconds from day before beginning of 1971 or
1972). It deviates from the number format of the Langseth et al. data by
adding several additional digits to the ``Time`` field in order to retain the
millisecond precision of the Nagihara et al. sets. It generates only one dT
column from the Nagihara et al. PDS sets: at each point, it selects the dT
measurement that the Nagihara team used to calculate per-thermometer
temperature values from average temperature. It converts data from the *JGR
Planets* release to T and dT by, respectively, averaging and computing the
difference of the explicitly-given thermometer temperatures.

### ``split``

This set presents the data in a more user-friendly fashion, with the following
features:

* retains the distinction between files 1, 2, 3, 4, and 5, but concatenates
	data from all sources. For instance, a17p1f1 contains data from 1972-1977.
* presents time as YMD UTC rather than mission epoch.  substitutes the
	``dT_corr`` field for ``dT.``  
* discards most flagged points from ``clean,``
	pruning impulsive outliers and missing data to reduce data cleaning time 
* discards the mangled ``HTR`` field from file 3. 
* gives explicit temperature values at each differential thermometer rather 
	than T and dT (it takes explicit temperatures from Nagihara et al.'s data
	and computes them from the Langseth et al. data).

### ``depth``

This set is basically ``split`` further reduced for ease of use in some
workflows, with the following features:

* file and probe distinction is discarded: temperature, time, and depth are
* given in one table per mission this set *only* includes depth values for
* sensors below the surface and discards values for TC4 to avoid ambiguity
* with topmost gradient thermometer

## ``flags`` Bitmask Values

* 1: missing data 

*only used in files 1,2,4,5:*

* 10: ambiguous bitflip correction 
* 100: outlier in dT 
* 1000: outlier in T

*only used in file 3:*

* 10000: outlier in HTR (placeholder; not used due to mangled HTR values)
* 100000: outlier in TREF 
* 1000000: outlier in TC1 
* 10000000: outlier in TC2
* 100000000: outlier in TC3 
* 1000000000: outlier in TC4 
* 10000000000: time discrepancy

*only used in a15p2f1:*

* 100000000000: cycle change spikes 
* 1000000000000: presumed sensitivity error
* 10000000000000: bit-flipped eclipse events

## sources

**Apollo 15 Data from 1971-1974**

M. G. S. Langseth, S. J. Keihm, J. L. Chute, H. K. Hills, and D. R. Williams.
APOLLO 15 HEAT FLOW THERMAL CONDUCTIVITY RDR SUBSAMPLED V1.0,
15A-L-HFE-3-THERMAL-CONDUCTIVITY-V1.0, 2014. NSSDC ID PSPG-00093, NSSDC ID of
older tapes 71-063C-06A.

**Apollo 17 Data from 1972-1974**

M. G. S. Langseth, S. J. Keihm, J. L. Chute, H. K. Hills, and D. R. Williams. 
APOLLO 17 HEAT FLOW THERMAL CONDUCTIVITY RDR SUBSAMPLED V1.0,
A17A-L-HFE-3-THERMAL-CONDUCTIVITY-V1.0, 2014b. NSSDC ID PSPG-00022, NSSDC ID
of older tapes 72-096C-01A

**Apollo 15 Data from 1975**

S. Nagihara and Y. Nakamura. Apollo 15 ALSEP ARCSAV Heat Flow Experiment
Calibrated Gradient Bridge Temperatures Collection (1975-092 to 1975-181),
2019. URL
http://pds-geosciences.wustl.edu/lunar/urn-nasa-pds-a15hfe_calibrated_arcsav/data/collection.xml.

**Apollo 17 Data from 1975**

S. Nagihara and Y. Nakamura. Apollo 17 ALSEP ARCSAV Heat Flow Experiment
Calibrated Gradient Bridge Temperatures Collection (1975-092 to 1975-181),
2019. URL
http://pds-geosciences.wustl.edu/lunar/urn-nasa-pds-a17hfe_calibrated_arcsav/data/collection.xml.

**Apollo 17 Data from 1976-1977**

Ancillary material released along with: S. Nagihara, W.S. Kiefer, P.T. Taylor,
and Y. Nakamura. Examination of the long-term subsurface warming observed at
the Apollo 15 and 17 sites utilizing the newly restored heat flow experiment
data from 1975 to 1977. J Geophysical Research: Planets, April 2018. doi:
10.1029/2018JE005579
