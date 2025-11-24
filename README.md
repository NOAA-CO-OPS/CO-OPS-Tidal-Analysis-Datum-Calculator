# CO-OPS Tidal Analysis Datum Calculator (TADC)


# Overview

This is code for the NOAA Tidal Analysis Datum Calculator tool. 

This tool computes tidal datums from a water level data time series. The TADC uses a Butterworth digital filter to remove high frequency (> 4 cycles/day) water level variability in order to identify tidal high and low waters from observed water level data. For official datums, CO-OPS uses a procedure called Curve Fit Manual Verification (CFMV) approach to identify tidal high and low waters that requires human oversight. The accuracy of this tool relative to the CO-OPS approach is described in detail in the following technical document:

Licate LA, Huang L and Dusek G (2017) A Comparison of Datums Derived from CO-OPS Verified Data Products and Tidal Analysis Datum Calculator. [NOAA Technical Report NOS CO-OPS 085](https://access.co-ops.nos.noaa.gov/datumcalc/docs/TechnicalReport.pdf).

For additional information, contact: NOAA's [Center for Operational Oceanographic Products and Services](https://tidesandcurrents.noaa.gov/)
    User Services: tide.predictions@noaa.gov


# Installation
```bash
pip install git+https://github.com/NOAA-CO-OPS/CO-OPS-Tidal-Analysis-Datum-Calculator
```
or
```bash
git clone https://github.com/NOAA-CO-OPS/dev_CO-OPS-Tidal-Analysis-Datum-Calculator
cd dev_CO-OPS-Tidal-Analysis-Datum-Calculator
conda env create -f environment.yml
conda activate tadc
```


# Usage 
You can run the TADC from the command line. Results will print to the terminal window:
```bash
python -m tadc.run --fname data\test_timeseries_simple.csv --Subordinate_Lat 37 --Subordinate_Lon -79
```
> [!NOTE]
> Use `--outfile_save_dir C:\your\dir` and `--make_plots True` to save results and monthly diagnostic plots.

Or, you can use TADC within other Python programs like so:
```Python
import tadc
out = tadc.run(fname='data/test_timeseries_simple.csv', Subordinate_Lat=37, Subordinate_Lon=-79)
summary = out.readme
high_lows = out.high_lows
monthly_plots = out.plots
sub_monthly_means = out.subordinate_monthly_means
datums = out.datums
```
Note that you can alternatively pass your timeseries directly to tadc as a `Pandas` `DataFrame`:
```Python
import pandas as pd
import tadc
timeseries = pd.read_csv('data/test_timeseries_simple.csv')
out = tadc.run(data=timeseries, Subordinate_Lat=37.8, Subordinate_Lon=-79.6)
```
> [!NOTE]
> One and only one of either `fname` or `data` must be specified when running `tadc()`


# Detailed Description of Modules

## run.py

This is the main module that imports all of the following modules to calculate datums of a given water level station using a .csv file input (example station water level data included) or data as a Pandas DataFrame. 

The file and data format should follow the standards stated in the NOAA Tidal Analysis Datums Calculator [User's Guide](https://access.co-ops.nos.noaa.gov/datumcalc/docs/UserGuide.pdf).

## qa.py and qc.py

These module defines some basic QA/QC assurances/checks for the input data. If a check is failed, a known `Error` will be generated.

## filter_defs.py

This module defines a low-pass Butterworth filter to preserve the tidal energy and remove the meteorological effects of the water level.

## control_data.py

This module retrieves, prepares and returns Monthly Means, High Lows, Accepted Datums of the control station (if selected) from [CO-OPS' API](https://tidesandcurrents.noaa.gov/api-helper/url-generator.html) and determines the datum computation method to be used.

## tides.py

This module identifies, tabulates and flags tides within the series and designates them as Lower Low, Higher Low, Lower High, or Higher High.  Tide picking algorithms for both semi-diurnal and diurnal tide signals are used and depend on the nature of the time series. 


# NOAA Open Source Disclaimer

This repository is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project code is provided on an 'as is' basis and the user assumes responsibility for its use. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.


# License

Software code created by U.S. Government employees is not subject to copyright in the United States (17 U.S.C. �105). The United States/Department of Commerce reserves all rights to seek and obtain copyright protection in countries other than the United States for Software authored in its entirety by the Department of Commerce. To this end, the Department of Commerce hereby grants to Recipient a royalty-free, nonexclusive license to use, copy, and create derivative works of the Software outside of the United States.
