# CO-OPS-Tidal-Analysis-Datum-Calculator



# Overview

This is code for the NOAA Tidal Analysis Datums Calculator tool. 

This tool computes tidal datums from any water level data time series with a consistent time step. The Tidal Analysis and Datums Calculator (TAD) uses a Butterworth digital filter to remove high frequency (> 4 cycles/day) water level variability in order to identify tidal high and low waters from observed water level data. For official datums, CO-OPS uses a procedure called Curve Fit Manual Verification (CFMV) approach to identify tidal high and low waters that requires human oversight. The accuracy of this tool relative to the CO-OPS approach is described in detail in the following technical document:

Licate LA, Huang L and Dusek G (2017) A Comparison of Datums Derived from CO-OPS Verified Data Products and Tidal Analysis Datum Calculator. NOAA Technical Report NOS CO-OPS 085. https://access.co-ops.nos.noaa.gov/datumcalc/docs/TechnicalReport.pdf

# Configuration file (config.cfg)

The online version of this tool has a GUI that allows users to select files, choose parameters and setup options. Since this is a standalone bit of code, included with the python files is a config.cfg file that can be modified in a text editor to set the required parameters for the code to run. Those parameters are: filename to be loaded, control station ID (if desired), Method, Time_Zone, ata units, your station’s latitude, your station’s longitude. 

Typically a user would simply need to set the filename, units, time zone, lat/lon, and choose a control station (if desired). The rest of the parameters can be determined by the code. 

# Code Description

The primary code is the SDC.py file. To run the code, open your preferred Python-enabled command prompt and enter:  

<code> python SDC.py </code>


## SDC.py

This is the main function that runs all of the following functions to calculate datums of a given water level station using a .csv file input (example station water level data included). 

This function looks to the config.cfg file for the user's file name, latitude and longitude to define the time zone, the data unit of the input file, and, if desired, the specified 7-digit control station ID to be used as a reference. The file and data format should follow the standards stated in the NOAA Tidal Analysis Datums Calculator User's Guide: https://access.co-ops.nos.noaa.gov/datumcalc/docs/UserGuide.pdf.

The final output includes .csv files containing the High Low tides, tide plots image files, and a main log text file called SDC.out. 

## filter_defs.py

This function defines a low-pass Butterworth filter to preserve the tidal energy and remove the meteorological effects of the water level.

## control_data.py

This function retrieves, prepares and returns Monthly Means, High Lows, Accepted Datums of the control station (if selected) from CO-OPS' API (https://tidesandcurrents.noaa.gov/api-helper/url-generator.html) and determines the datum computation method to be used.

## tides.py

This function identifies, tabulates and flags tides within the series and designates them as Lower Low, Higher Low, Lower High, or Higher High.  Tide picking algorithms for both semi-diurnal and diurnal tide signals are used and depend on the nature of the time series. 

For additional information, contact:
Jerry Hovis (Gerald.Hovis@noaa.gov),
NOAA Center for Operational Oceanographic Products and Services

## Example.csv

This is an example data file that is formatted correctly to run with the code. 

# NOAA Open Source Disclaimer

This repository is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project code is provided on an 'as is' basis and the user assumes responsibility for its use. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.

# License

Software code created by U.S. Government employees is not subject to copyright in the United States (17 U.S.C. �105). The United States/Department of Commerce reserves all rights to seek and obtain copyright protection in countries other than the United States for Software authored in its entirety by the Department of Commerce. To this end, the Department of Commerce hereby grants to Recipient a royalty-free, nonexclusive license to use, copy, and create derivative works of the Software outside of the United States.
