# CO-OPS Tidal Analysis Datum Calculator (TADC)
[![PyPI version](https://img.shields.io/pypi/v/tadc.svg?color=blue)](https://pypi.org/project/tadc/)
[![Conda Forge version](https://img.shields.io/conda/vn/conda-forge/tadc.svg?color=blue)](https://anaconda.org/conda-forge/tadc)
[![Web App](https://img.shields.io/badge/Web_App_GUI-005b96.svg)](https://access.co-ops.nos.noaa.gov/datumcalc/)

# Overview

This is code for the NOAA Tidal Analysis Datum Calculator tool. 

A [web-based GUI](https://access.co-ops.nos.noaa.gov/datumcalc/) for single-file computations is also available.

This tool computes tidal datums from a water level data time series. The TADC uses a Butterworth digital filter to remove high frequency (> 4 cycles/day) water level variability in order to identify tidal high and low waters from observed water level data. For official datums, CO-OPS uses a procedure called Curve Fit Manual Verification (CFMV) approach to identify tidal high and low waters that requires human oversight. The accuracy of this tool relative to the CO-OPS approach is described in detail in the following technical document:

Licate LA, Huang L and Dusek G (2017) A Comparison of Datums Derived from CO-OPS Verified Data Products and Tidal Analysis Datum Calculator. [NOAA Technical Report NOS CO-OPS 085](https://access.co-ops.nos.noaa.gov/datumcalc/docs/TechnicalReport.pdf).

For additional information, contact: NOAA's [Center for Operational Oceanographic Products and Services](https://tidesandcurrents.noaa.gov/)
    User Services: tide.predictions@noaa.gov


# Installation
```bash
pip install tadc
```
or
```bash
conda install tadc -c conda-forge
```
or
```bash
git clone https://github.com/NOAA-CO-OPS/CO-OPS-Tidal-Analysis-Datum-Calculator
cd CO-OPS-Tidal-Analysis-Datum-Calculator
conda env create -f environment.yml
conda activate tadc
```

# Required Data Format
The input `.csv` must have only two columns and one header row, with the date/time column being on the left and water levels on the right. Most any date/time format can be used. 


# Usage 
You can run the TADC from the command line. You must clone the repo to do this. Results will print to the terminal window:
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

The object returned by `tadc.run()` also supports some further analyses on the data, in particular `daily_max_analysis`, `daily_min_analysis`, and `inundation_analysis`:
```Python
out = tadc.run(data=timeseries, Subordinate_Lat=37.8, Subordinate_Lon=-79.6)
daily_max_obj = out.daily_max_analysis(datum='MHHW')  # Do an analysis of daily maximums on MHHW #
daily_maxs = daily_max_obj.daily_maxs  # returns the daily max timeseries #
daily_max_obj.plot()  # creates a timeseries plot of the daily maxs #
daily_min_obj = out.daily_min_analysis(datum='MHHW')  # Do an analysis of daily minimums on MHHW #
# daily_min_obj has the same sub-methods as daily_max_obj #
inundation_analysis_obj = out.inundation_analysis(threshold=0.5, threshold_datum='MHHW') # Do an analysis of inundations above MHHW+0.5 m #
inundation_analysis_obj.inundations()  # returns the intervals of threshold exceedance #
inundation_analysis_obj.plot()  # creates a plot of inundation time vs. depth #
```

# Parameters:
Whether running from the command line or in Python, the following input arguments can be passed:

| Argument | Type  | Default | Description |
| :--- | :--- | :--- | :--- |
| `--fname` | `file path` | `None` | Required unless `--data` is specified.<br>The full path of the timeseries data to be analyzed.
| `--data` | `pandas DataFrame` | `None` | Required unless `--fname` is specified.<br>The timeseries data to be analyzed as a Pandas DataFrame object. Only used in Python workflows.
| `--Subordinate_Lat` | `float` | `None` | The latitude (in decimal degrees) of the station for which datums are being computed.
| `--Subordinate_Lon` | `float` | `None` | The longitude (in decimal degrees) of the station for which datums are being computed.
| `--Units` | {`'Meters'`,<br>`'Centimeters'`,<br>`'Feet'`,<br>`'Inches'`} | `'Meters'` | The units of the input data. All calculations will be done on these units.<br>**Warning** Must match the units of the input data.
| `--Time_Zone` | `str` | `'GMT'` | The time zone of the input data. All calculations will be done in this timezone.<br>**Warning** Must match the time zone of the input data.
| `--resample_minutes` | `int` | `None` | Resample the input data to this sampling rate (in minutes).<br> **Note** If your data are unevenly spaced in time, resampling will be done automatically at a detected interval which can be modified by setting this parameter.
| `--Control_Station_ID` | `int` | `None` | The seven-digit station ID of the NWLON control station being used, if using a control station.
| `--epoch_start_year` | `int` | `1983` | The beginning of the 19 year datum epoch on which to compute datums. To use, must also pass `--Control_Station_ID`.
| `--Method_Option` | {`'AUTO'`,<br>`'FRED'`,<br>`'TBYT'`,<br>`'MMSC'`} | `'AUTO'` | The datum adjustment method to use when using a control station.<br>`'AUTO'`: Auto select. If using a control station and <1 month of data 'TBYT' will be used, 'MMSC' otherwise. If not using a control station 'FRED' will be used.<br>`'FRED'`: First reduction method. No adjustment from control station, simple average of highs and lows.<br>`'TBYT'`: Tide by tide analysis. Compares simultaneous high/low waters at the input and control stations.<br>`'MMSC'`: Monthly means simultaneous comparison. Compares monthly means at the input and control stations.
| `--Pick_Method` | {`'PolyFit'`,<br>`'Direct'`} | `'PolyFit'` | The method to pick high and low tide times and values.<br>`'PolyFit'`: Fit a polynomial to the observations and extract highs and lows from the fit.<br>`'Direct'`: Extract highs and lows directly from the input timeseries.
| `--outfile_save_dir` | `directory path` | `None` | The directory into which to save output files, if run from the command line. No effect if run in a Python workflow.<br>**Note** A timestamped folder will be created containing the output files.
| `--make_plots` | `bool` | `False` | Whether or not to create plots, one per month, showing the extracted high and low tides.<br>**Warning** To save plots when running from the command line, `--outfile_save_dir` must also be specified. If running in a Python workflow, plots will be created as Python objects.
| `--daily_max_analysis` | {`'MHHW'`,<br>`'MHW'`,<br>`'MSL'`,<br>`'MLW'`,<br>`'MLLW'`}  | `None` |  Complete an analysis of daily maxima on this datum.
| `--daily_min_analysis` | {`'MHHW'`,<br>`'MHW'`,<br>`'MSL'`,<br>`'MLW'`,<br>`'MLLW'`}  | `None` |  Complete an analysis of daily minima on this datum.
| `--inundation_analysis` | {`'MHHW'`,<br>`'MHW'`,<br>`'MSL'`,<br>`'MLW'`,<br>`'MLLW'`}  | `None` |  Complete an analysis of inundations on this datum.
| `--inundation_threshold` | `int` | `None` |  Complete an analysis of inundations above this threshold.


# NOAA Open Source Disclaimer

This repository is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project code is provided on an 'as is' basis and the user assumes responsibility for its use. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.


# License

Software code created by U.S. Government employees is not subject to copyright in the United States (17 U.S.C. §105). The United States/Department of Commerce reserves all rights to seek and obtain copyright protection in countries other than the United States for Software authored in its entirety by the Department of Commerce. To this end, the Department of Commerce hereby grants to Recipient a royalty-free, nonexclusive license to use, copy, and create derivative works of the Software outside of the United States.
