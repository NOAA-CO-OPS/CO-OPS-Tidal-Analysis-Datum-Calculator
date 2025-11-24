#!/usr/bin/python
"""
NOAA/CO-OPS Simple Datum Calculator Engine

@author: George Story

Convert to Python 3  February 2020

"""

import argparse
import configparser
from datetime import datetime, date, time, timedelta
from dateutil.parser import parse
import matplotlib
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
from numpy import diff, sign, mean
import os
import pandas as pd
from scipy.signal import butter, filtfilt
import sys
from time import gmtime, strftime
import warnings

from . import control_data as cd
from . import daily_max_analysis as dma
from . import filter_defs as fd
from . import inundation_analysis as ia
from . import qa
from . import qc
from . import tides as tf


def Get_Method(xtimes, OutFile):
    #This function chooses method to use
    #If there is at least 1 complete month of data MMSC else TBYT
    i = 0
    while xtimes[i].day > 1 and i < (len(xtimes)-1):
    #Find the first month start
        i = i+1

    if i == len(xtimes)-1:
        OutFile = SDC_Print(['Could not find month start'], OutFile)
        return 'TBYT', OutFile
    ld = tf.Last_Day_In_Month(xtimes[i].year,xtimes[i].month)

    #See if the data continues to the end of the month
    while xtimes[i].day < ld and i < (len(xtimes)-1):
        i = i + 1

    if xtimes[i].day == ld:
        return "MMSC", OutFile
    else:
        return "TBYT", OutFile
            

def Fill_Gaps(x,y, OutFile):
    #This function fills gaps   
    ts = pd.DataFrame({'time':x,'val':y})
    ts.loc[ts['val']==-99999.99,'val'] = np.nan
    # Interpolate through all the gaps with a 2nd order polynomial #
    ts['val'] = ts['val'].interpolate(method='polynomial', order=2)
    # Go back and get rid of the interpolations through gaps longer than 3 hr #
    interval = ts['time'].diff().iloc[1].seconds
    max_gap_size = (3 * 60 * 60) / interval
    s = pd.Series(pd.Series(y)[pd.Series(y)==-99999.99].index.values)
    group_ids = (s.diff() != 1).cumsum()
    grouped = s.groupby(group_ids)
    groups_list = grouped.apply(list)
    groups_list_nofill = groups_list[grouped.size() > max_gap_size]
    for i in range(len(groups_list_nofill)):
        ts.loc[groups_list_nofill.iloc[i], 'val'] = -99999.99
    x = ts['time'].values
    y = ts['val'].values
    remaining_gaps = len(groups_list_nofill)
    OutFile = SDC_Print([str(len(groups_list)) + ' gaps identified. ' + str(len(groups_list) - remaining_gaps) + ' gaps filled, the remaining are too long to be filled.'], OutFile)
    return y, remaining_gaps, OutFile


def Longest_Segment(x,y, OutFile):
    #This function trims dataset to just the longest continuous segment
    if len(x) != len(y):
        OutFile = SDC_Print(['*** Error x and y are different lengths'], OutFile)
        raise RuntimeError('x and y are different lengths')
    TS = []
    #Create an empty list and build (append) a list of the segments
    inseg = False
    for i in range(len(x)-1):
        if inseg:
            #In a segment, end it if a missing point is found
            if y[i] == -99999.99:
                inseg = False
                TS.append((start_loc, i-1))
        else:
            #Not in a segment - start one if there is data here
            if y[i] != -99999.99:
                inseg = True
                start_loc = i
                
    if inseg:
        TS.append((start_loc, len(x)-1))    
    #find the longest string
    slen = timedelta(hours=0)
    sloc = None
    for i in range(len(TS)):
         if (x[TS[i][1]] - x[TS[i][0]]) > slen:
            slen = x[TS[i][1]] - x[TS[i][0]]
            sloc = i
    #slice out the longest continuous segment
    x2 = x[TS[sloc][0]:TS[sloc][1]+1]
    y2 = y[TS[sloc][0]:TS[sloc][1]+1]
    return(x2,y2, OutFile)


def SDC_Print(PLines, OutFile):
    #This is a print function 
    OutLine = ''
    for i in range(len(PLines)):
        OutLine = OutLine + str(PLines[i]) + ' '
    OutFile = OutFile + OutLine + '\n'
    return OutFile


def check_required_inputs(fname, data):
    if fname is None and data is None:
        raise ValueError('Either `fname` or `data` must be specified.')
    elif fname is not None and data is not None:
        raise ValueError('Only one of `fname` or `data` can be specified.')
    else:
        if fname is not None:
            if not isinstance(fname,str):
                raise ValueError('Input arg `fname` must be a string')
        elif data is not None:
            if not isinstance(data,pd.core.frame.DataFrame):
                raise ValueError('Input arg `data` must be a two-column pandas DataFrame.')


class Out:
    def __init__(self, data, readme, plots, high_lows, subordinate_monthly_means, datums):
        self.data = data
        self.readme = readme
        self.plots = plots
        self.high_lows = high_lows
        self.subordinate_monthly_means = subordinate_monthly_means
        self.datums = datums
        self.datums = {key: float(value) if isinstance(value, np.float64) else value for key, value in self.datums.items()}

    def inundation_analysis(self, threshold, threshold_datum):
        out_ia = ia.run(threshold, threshold_datum, self.data, self.datums, self.high_lows)
        return out_ia

    def daily_max_analysis(self, datum):
        out_dma = dma.run(datum, self.data, self.datums)
        return out_dma


def run(*, fname=None, data=None, Pick_Method='PolyFit', Control_Station_ID=None, Method_Option='AUTO',
         Time_Zone='GMT', Units='Meters', Subordinate_Lat=None, Subordinate_Lon=None, outfile_save_dir='None', make_plots=False):

    #Init the output readme
    OutFile = ''

    if Control_Station_ID == None or len(str(Control_Station_ID)) < 7:
        Method_Option = 'FRED'
    else:
        Control_Station_ID = str(Control_Station_ID)

    #when the script fails to execute as a result of incorrect parameter submission, reason will be indicated in GUI as well as the logfile 
    if fname == '':
        OutFile = SDC_Print(['No filename on command line. Arguments:'], OutFile)
        OutFile = SDC_Print(['filename = filespec of input to use'], OutFile)
        OutFile = SDC_Print(['control  = 7 character control station ID (or None)'], OutFile)
        OutFile = SDC_Print(['method  = AUTO, TBYT or FRED'], OutFile)
        OutFile = SDC_Print(['time_zone = some tring with the gmt offset at the end e.g. UST5'], OutFile)
        OutFile = SDC_Print(['units = Meters, Centimeters, Millimeters, Feet, Inches'], OutFile)
        OutFile = SDC_Print(['lat  = latitude of station'], OutFile)
        OutFile = SDC_Print(['lon  = longitude of station'], OutFile)
     
    OutFile = SDC_Print (['Run Time: ', strftime("%Y-%m-%d %H:%M:%S", gmtime())], OutFile)

    # Run the QA/QC to ensure input data is ok and then initiate the data #
    check_required_inputs(fname,data)
    if fname is not None:
        end_of_path = fname.rfind('/')
        if end_of_path > -1:
            path = fname[0:end_of_path+1]
        else:
            path = ''
        OutFile = SDC_Print(['Using ', fname[end_of_path+1:]], OutFile)
        ts_qa = qa.run(pd.read_csv(fname))
        qc.run(ts_qa, Control_Station_ID, Subordinate_Lat, Subordinate_Lon)      
    else:
        OutFile = SDC_Print(['Using user input timeseries'], OutFile)
        ts_qa = qa.run(data)
        qc.run(ts_qa, Control_Station_ID, Subordinate_Lat, Subordinate_Lon)
        
    #Get time offset if subordinate is not gmt
    OutFile = SDC_Print(['Time Zone = ' + Time_Zone], OutFile)
    hrstr = ''
    n = len(Time_Zone)-1
    while Time_Zone[n].isdigit():
        hrstr = hrstr + Time_Zone[n]
        n = n-1
    if len(hrstr) > 0:
        gmt_offset = int(hrstr)
    else:
        gmt_offset = 0

    #Get Date-Times and Water Levels from csv file using Pandas
    data = ts_qa
    data[data.columns[0]] = pd.to_datetime(data[data.columns[0]])
    MissingPoints = len(data[data.columns[1]][data[data.columns[1]].isna()])
    data = data.fillna(value=-99999.99)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', category=FutureWarning)
        x = np.array(data[data.columns[0]].dt.to_pydatetime())
    y = np.array(data[data.columns[1]])

    OutFile = SDC_Print([len(x), ' data points loaded.'], OutFile)

    #Determine interval and check for consistency
    Interval = x[1] - x[0]
    OutFile = SDC_Print(['Interval: ', Interval], OutFile)
    i = 1
    ni = len(x)
    #for i in range(len(x)-1):
    while i < ni-3:
        if (x[i+1] - x[i]) > Interval and ((x[i+1] - int(x[i]).seconds % Interval.seconds) == 0):     
            OutFile = SDC_Print(['Interval of ', x[i+1] - x[i], ' at ', x[i]], OutFile)
            nmissing = ((x[i+1] - x[i]).seconds / Interval.seconds) - 1
            misstime = x[i+1] - Interval
            for ii in range(int(nmissing)):
                x = np.insert(x, i+1, misstime)
                y = np.insert(y, i+1, -99999.99)

                misstime = misstime - Interval
                MissingPoints = MissingPoints + 1
                ni = ni+1
            i = i + nmissing
        elif(((x[i+1] - x[i]) < Interval) or ((x[i+1] - x[i]).seconds % Interval.seconds != 0)):
            OutFile = SDC_Print(['Interval of ', x[i+1] - x[i], ' at ', x[i]], OutFile)
            OutFile = SDC_Print(['***Error*** Time interval is inconsistent.'], OutFile)
            raise RuntimeError('Time interval is inconsistent.')
        else:
            i = i+1

    ngaps = 0
    if MissingPoints > 0:
        OutFile = SDC_Print([MissingPoints, ' missing data points'], OutFile)
        y, remaining_gaps, OutFile = Fill_Gaps(x,y, OutFile)
    else:
        remaining_gaps = 0

    #It there are gaps that can not be filled (>3 hrs) then 
    #just take longest string withou breaks
    if (remaining_gaps > 0):
        x, y, OutFile = Longest_Segment(x,y, OutFile)
        OutFile = SDC_Print([' '], OutFile)
        OutFile = SDC_Print(['Dataset trimmed to longest continuous segment.'], OutFile)
        OutFile = SDC_Print(['New Start:', x[0]], OutFile)
        OutFile = SDC_Print(['New End:  ', x[len(x)-1]], OutFile)
        warnings.warn(('Dataset trimmed to longest continuous segment with gaps <= 3 hours: ' + datetime.strftime(x[0],'%Y-%m-%d %H:%M') +
                       ' - ' + datetime.strftime(x[-1],'%Y-%m-%d %H:%M')))

    if (x[len(x)-1] - x[0]) < timedelta(days=14):
        OutFile = SDC_Print(['***Error*** Not enough data for analysis. 2 weeks minimum'], OutFile)
        raise RuntimeError('Not enough data for analysis. 2 weeks minimum')

    #Check for input units.  Get conversion factor from meters
    OutFile = SDC_Print([Units.upper()], OutFile)
    if Units.upper() == 'METERS':
        CFactor = 1.0
        fmt = "%.3f"
    elif Units.upper() == 'CENTIMETERS':
        CFactor = 100.0
        fmt = "%.1f"
    elif Units.upper() == 'MILLIMETERS':
        CFactor = 1000.0
        fmt = "%.0f"
    elif Units.upper() == 'FEET':
        CFactor = 3.28084
        fmt = "%.2f"
    elif Units.upper() == 'INCHES':
        CFactor = 39.3701
        fmt = "%.1f"
    else:
        OutFile = SDC_Print(['***Error*** Input units of', Units, 'not defined.'], OutFile)
        raise RuntimeError('Input units of ' + Units + ' not defined.')

    OutFile = SDC_Print([''], OutFile)
    OutFile = SDC_Print(['All calculations and results are in', Units], OutFile)

    #Determine calulation method MMSC or TBYT
    #Use MMSC if there is a complete calendar month of data and a control station is selected
    Calc_Method, OutFile = Get_Method(x, OutFile)
    if Method_Option == 'TBYT':
        Calc_Method = 'TBYT'
    if Method_Option == 'FRED':
        Calc_Method = 'FRED'

    if Method_Option == 'AUTO':
        if Calc_Method == 'TBYT':
            OutFile = SDC_Print(['Less than one month of data loaded.  Use Tide-By-Type comparison.'], OutFile)
        else:
            OutFile = SDC_Print(['At least one complete month loaded.  Use Monthly Means comparison.'], OutFile)

    if Calc_Method == 'FRED':
        #Get Sub-Method from subordinate station location
        try:
            flon = float(Subordinate_Lon)
        except:
            OutFile = SDC_Print(['***Error***  Station Longitude is not a number:', Subordinate_Lon], OutFile)
            raise RuntimeError('Station Longitude is not a number:' + str(Subordinate_Lon))
        if flon < -100.:
            Sub_Method = 'Standard'
        else:
            Sub_Method = 'Modified'
    else:
        #Get Sub-Method from Control Station Location
        Sub_Method = cd.Get_SubMethod(Control_Station_ID)
        if (Sub_Method not in ['Standard', 'Modified']):
            OutFile = SDC_Print([Sub_Method], OutFile)
    OutFile = SDC_Print([''], OutFile)
    if Sub_Method == 'Standard':
        OutFile = SDC_Print(['West coast/Pacific station:\n   Using Standard Range Ratio Method'], OutFile)
    else:
        OutFile = SDC_Print(['Gulf/East coast station:\n   Using Modified Range Ratio Method'], OutFile)
    OutFile = SDC_Print([''], OutFile)


    """Set up Filter parameters."""
    fs= 86400 / Interval.seconds 
    order = 6
    #cutoff = 5.0  #desired cutoff frequency of the filter, per day
    cutoff = 4.0
    OutFile = SDC_Print(['Sampling Rate: ', fs, ' per day.   Using cutoff frequency of ', cutoff , ' per day'], OutFile)
    #Get the filter coefficients 
    b, a = fd.butter_lowpass(cutoff, fs, order)

    #Filter the data, and plot both the original and filtered signals.
    filt = fd.butter_lowpass_filter(y, cutoff, fs, order)

    #find inflection points (tides) in fitered signal
    highs = (diff(sign(diff(filt))) < 0).nonzero()[0] + 1 #local max
    lows  = (diff(sign(diff(filt))) > 0).nonzero()[0] + 1 #local min

    #check potential tides for spacing in time and height
    highs_mask, lows_mask = tf.Check_Tides(x, y, highs, lows, CFactor) 

    # Delete bad tides
    highs = highs[highs_mask]
    lows = lows[lows_mask]

    #Check Tide Order
    CTO = tf.Check_Tide_Order(x, highs, lows)
    if (CTO < 0):
        OutFile = SDC_Print(["***Warning*** - Tides are out of order"], OutFile)

    high_values = []
    high_dts = []
    low_values = []
    low_dts = []
    if Pick_Method == 'PolyFit':
        #Use a polynomial curve fit to select the extreme
        for i in range(len(highs)):
            high_dt, high_val = tf.Local_Max_Fit(x, y, highs[i])
            high_values.append(high_val)
            high_dts.append(high_dt)
        for i in range(len(lows)):
            low_dt, low_val = tf.Local_Min_Fit(x, y, lows[i])
            low_values.append(low_val)
            low_dts.append(low_dt)
    else:
        #Just pick the highest/lowest point within specified window (+- 30 minutes)
        for i in range(len(highs)):
            high_dt, high_val = tf.Local_Max(x, y, highs[i], timedelta(minutes=30))
            high_values.append(high_val)
            high_dts.append(high_dt)
        for i in range(len(lows)):
            low_dt, low_val = tf.Local_Min(x, y, lows[i], timedelta(minutes=30))
            low_values.append(low_val)
            low_dts.append(low_dt)

    ntides = len(highs) + len(lows)
    t = (x[len(x)-1]-x[0])

    #Highest and lowest water levels
    LWL = np.amin(y)
    HWL = np.amax(y)

    Mean_Level = np.mean(y, dtype=np.float64)
    OutFile = SDC_Print(['Data Start: ' , str(x[0])], OutFile)
    OutFile = SDC_Print(['Data End  : ' , str(x[len(x)-1])], OutFile)
    OutFile = SDC_Print(['Mean Water Level: ', fmt % Mean_Level], OutFile)
    OutFile = SDC_Print(['Highest Water Level: ', fmt % HWL], OutFile)
    OutFile = SDC_Print(['Lowest Water Level: ', fmt % LWL], OutFile)
    OutFile = SDC_Print(['Duration: ', t], OutFile)
    OutFile = SDC_Print(['High Tides Found: ', len(highs)], OutFile)
    OutFile = SDC_Print(['Low Tides Found : ', len(lows)], OutFile)
    ndays = float(t.days) + (t.seconds + 360.)/86400.0
    tpd = float(ntides)/float(ndays)
    OutFile = SDC_Print(['Tides per day: {0:.1f}'.format(tpd)], OutFile)
    if (tpd > 3.0):
        OutFile = SDC_Print([('Semi-Diurnal - Using EXHL')], OutFile)
        high_types, low_types = tf.EXHL(high_values, low_values)
    else:
        OutFile = SDC_Print([('Diurnal Using DIUR')], OutFile)
        high_types, low_types = tf.DIUR(high_dts, high_values, low_dts, low_values, x[0])

    #Summarize Extremes
    OutFile = SDC_Print([high_types.count('H'), ' Highs'], OutFile)
    OutFile = SDC_Print([high_types.count('HH'), ' Higher Highs'], OutFile)
    OutFile = SDC_Print([low_types.count('L'), ' Lows'], OutFile)
    OutFile = SDC_Print([low_types.count('LL'), ' Lower Lows'], OutFile)
    OutFile = SDC_Print([' '], OutFile)

    #Store the highs and lows in an object
    li=0
    hi=0
    rows = []
    while ((hi < len(highs)) and (li < len(lows))):
        if (low_dts[li] < high_dts[hi]):
            row = {'time':['{0:%Y-%m-%d %H:%M}'.format(low_dts[li])], 'value':low_values[li], 'tide type':low_types[li]}
            rows.append(row)
            li = li + 1
        else:
            row = {'time':['{0:%Y-%m-%d %H:%M}'.format(high_dts[hi])], 'value':high_values[hi], 'tide type':high_types[hi]}
            rows.append(row)
            hi = hi + 1
    else:
        if li >= (len(lows)):
            while (hi < len(highs)):
                row = {'time':['{0:%Y-%m-%d %H:%M}'.format(high_dts[hi])], 'value':high_values[hi], 'tide type':high_types[hi]}
                rows.append(row)
                hi = hi + 1
        if hi >= (len(highs)):
            while (li < len(lows)):
                row = {'time':['{0:%Y-%m-%d %H:%M}'.format(low_dts[li])], 'value':low_values[li], 'tide type':low_types[li]}
                rows.append(row)
                li = li + 1
    high_lows = pd.DataFrame(rows)
    high_lows['time'] = [datetime.strptime(high_lows.iloc[i]['time'][0],'%Y-%m-%d %H:%M') for i in range(len(high_lows))]  # Convert str times to datetime objects #
    
    #Generate Plots of wl and high-lows for each month and save as object
    pn = 1
    m1 = x[0].year * 12 + x[0].month - 1
    m2 = x[len(x)-1].year * 12 + x[len(x)-1].month - 1
    nyears = x[len(x)-1].year-x[0].year + 1
    times = []
    figs = []
    for m in range(m1,m2+1):
        yr = m//12
        mn = m+1-yr*12
        p1 , p2 = tf.first_last_in_month(x, mn, yr) 

        MHighs = []
        MHHighs = []
        MHTimes = []
        MHHTimes = []
        MLows = []
        MLLows = []
        MLTimes = []
        MLLTimes = []
        for ii in range(len(high_dts)):
            if high_dts[ii] >= x[p1] and high_dts[ii]<= x[p2]:
                if high_types[ii] == 'HH':
                    MHHTimes.append(high_dts[ii])
                    MHHighs.append(high_values[ii])
                else:
                    MHTimes.append(high_dts[ii])
                    MHighs.append(high_values[ii])
        for ii in range(len(low_dts)):
            if low_dts[ii] >= x[p1] and low_dts[ii]<= x[p2]:
                if low_types[ii] == 'LL':
                    MLLTimes.append(low_dts[ii])
                    MLLows.append(low_values[ii])
                else:
                    MLTimes.append(low_dts[ii])
                    MLows.append(low_values[ii])

        if make_plots:
            fig,ax = plt.subplots(1)
            ax.plot(x[p1:p2], y[p1:p2], 'b-', label='Water Level')
            ax.plot(MHHTimes, MHHighs, label = 'Higher Highs', marker='D', markersize=3, linestyle='None', color='r')
            ax.plot(MHTimes, MHighs, label = 'Highs', marker='o', markersize=3, linestyle='None', color='m')
            ax.plot(MLLTimes, MLLows, label = 'Lower Lows', marker='D', markersize=3, linestyle='None', color='r')
            ax.plot(MLTimes, MLows, label = 'Lows', marker='o', markersize=3, linestyle='None', color='m')

            ax.set_ylabel(Units)
            ax.grid('on')

            majorLocator = matplotlib.ticker.MultipleLocator(5)
            minorLocator = matplotlib.ticker.MultipleLocator(1)
            xax = ax.get_xaxis() 
            xax.set_major_locator(majorLocator)
            #format major xtick label
            xax.set_major_formatter(mdates.DateFormatter('%m/%d/%y'))
            times.append('Month' + str(pn))
            figs.append(fig)
            plt.close(fig)
        pn = pn + 1

    if make_plots:
        out_plots = pd.DataFrame({'time':times,'plot':figs})
        OutFile = SDC_Print([pn-1, ' Monthly plots generated\n'], OutFile)
    else:
        out_plots = pd.DataFrame({})

    # Calculate Datums by First Reduction
    #High Means and Highest Tide 
    MHHW = 0.0
    MHW = 0.0
    HWL = -99999.99
    nhighs = 0
    nhhighs = 0
    for i in range(len(highs)):
        if (high_types[i] == 'HH'):
            MHHW = MHHW + high_values[i]
            nhhighs = nhhighs + 1
            MHW = MHW + high_values[i]
            nhighs = nhighs + 1
        if (high_types[i] == 'H'):
            MHW = MHW + high_values[i]
            nhighs = nhighs + 1
        if (high_types[i] != 'H' and high_types[i] != 'HH'):
            SDC_Print(["Bad high type", high_types[i]])

        if high_values[i] > HWL:
            HWL = high_values[i]
            HWL_DT = high_dts[i]

    MHHW = MHHW / nhhighs
    MHW  = MHW  / nhighs
        
    #Low Means and Lowest Tide 
    MLLW = 0.0
    MLW = 0.0
    LWL = 99999.99
    nlows = 0
    nllows = 0
    for i in range(len(lows)):
        if (low_types[i] == 'LL'):
            MLLW = MLLW + low_values[i]
            nllows = nllows + 1
            MLW = MLW + low_values[i]
            nlows = nlows + 1
        if (low_types[i] == 'L'):
            MLW = MLW + low_values[i]
            nlows = nlows + 1
        if (low_types[i] != 'L' and low_types[i] != 'LL'):
            SDC_Print(["Bad low type", low_types[i]])

        if low_values[i] < LWL:
            LWL = low_values[i]
            LWL_DT = low_dts[i]

    MLLW = MLLW / nllows
    MLW  = MLW  / nlows

    if Calc_Method == 'FRED':
        OutFile = SDC_Print([' '], OutFile)
        OutFile = SDC_Print([' TIDAL Datums by Arithmetic Mean of Your Data (First Reduction):'], OutFile)
        OutFile = SDC_Print(['HWL  = ', fmt % HWL, '  (' + HWL_DT.strftime("%Y/%m/%d %H:%M") + ')'], OutFile)
        OutFile = SDC_Print(['MHHW = ', fmt % MHHW], OutFile)
        OutFile = SDC_Print(['MHW  = ', fmt % MHW], OutFile)
        OutFile = SDC_Print(['DTL  = ', fmt % (0.5 * (MHHW + MLLW))], OutFile)
        OutFile = SDC_Print(['MTL  = ', fmt % (0.5 * (MHW + MLW))], OutFile)
        OutFile = SDC_Print(['MSL  = ', fmt % Mean_Level], OutFile)
        OutFile = SDC_Print(['MLW  = ', fmt % MLW], OutFile)
        OutFile = SDC_Print(['MLLW = ', fmt % MLLW], OutFile)
        OutFile = SDC_Print(['DHQ  = ', fmt % (MHHW - MHW)], OutFile)
        OutFile = SDC_Print(['DLQ  = ', fmt % (MLW - MLLW)], OutFile)
        OutFile = SDC_Print(['MN   = ', fmt % (MHW - MLW)], OutFile)
        OutFile = SDC_Print(['GT   = ', fmt % (MHHW - MLLW)], OutFile)
        OutFile = SDC_Print(['LWL  = ', fmt % LWL , '  (' + LWL_DT.strftime("%Y/%m/%d %H:%M") + ')'], OutFile)
        OutFile = SDC_Print([' '], OutFile)

        datums = {'HWL':HWL,
                  'MHHW':MHHW,
                  'MHW':MHW,
                  'DTL':0.5 * (MHHW + MLLW),
                  'MTL':0.5 * (MHW + MLW),
                  'MSL':Mean_Level,
                  'MLW':MLW,
                  'MLLW':MLLW,
                  'DHQ':MHHW - MHW,
                  'DLQ':MLW - MLLW,
                  'MN':MHW - MLW,
                  'GT':MHHW - MLLW,
                  'LWL':LWL}

    if Calc_Method == 'MMSC' or Calc_Method == 'TBYT':
    #Get Accepted Datums for Control Station
        Control_Acc_Datums = cd.Get_Accepted_Datums(Control_Station_ID, CFactor)
        if (Control_Acc_Datums[0] == None or Control_Acc_Datums[1] == None or 
           Control_Acc_Datums[5] == None or Control_Acc_Datums[6] == None):
            OutFile = SDC_Print(['***Error*** Problem retrieving Accepted Datums for station ', Control_Station_ID], OutFile)
            raise RuntimeError('Problem retrieving Accepted Datums for station ' + str(Control_Station_ID))
        OutFile = SDC_Print(['Control Datums for: ' , Control_Station_ID], OutFile)
        OutFile = SDC_Print(['\nMHHW,  MHW,  DTL,  MTL,  MSL,  MLW,  MLLW'], OutFile)
        MeanString = ''
        for di in range(0,7):
            MeanString = MeanString + fmt % Control_Acc_Datums[di] + ' '
        OutFile = SDC_Print([MeanString], OutFile)

        OutFile = SDC_Print(['GT,  MN,  DHQ,  DLQ,  NAVD,  LWI,  HWI'], OutFile)
        MeanString = ''
        for di in range(7,14):
            if (Control_Acc_Datums[di] == None):
                MeanString = MeanString + 'Null '
            else:
                MeanString = MeanString + fmt % Control_Acc_Datums[di] + ' '
        OutFile = SDC_Print([MeanString], OutFile)

        Control_DHQ = Control_Acc_Datums[9]
        Control_DLQ = Control_Acc_Datums[10]
        Control_MLW = Control_Acc_Datums[5]
        Control_MLLW = Control_Acc_Datums[6]
        Control_MHW = Control_Acc_Datums[1]
        Control_MHHW = Control_Acc_Datums[0]
        Control_MTL = Control_Acc_Datums[3]
        Control_MSL = Control_Acc_Datums[4]
        Control_DTL = Control_Acc_Datums[2] 
        Control_MN = Control_Acc_Datums[1] - Control_Acc_Datums[5]
        Control_GT = Control_Acc_Datums[0] - Control_Acc_Datums[6]

    #if Calc_Method == 'MMSC', Loop Month-by-month and calculate means and Store in MM_Subordinate list
    MM_Subordinate = []
    OutFile = SDC_Print(['\nSUBORDINATE MONTHLY MEANS:'], OutFile)

    # Calc subordinate monthly means #
    if x[0].year == x[-1].year and x[0].month == x[-1].month and (x[0].day > 1 or x[-1].day < tf.Last_Day_In_Month(x[-1].year,x[-1].month)):  # If the dataset is < 1 month long, cannot calculate subordinate monthly means
        subordinate_monthly_means = pd.DataFrame({})
    else:
        #Find start of first full month
        p1 = 0
        while x[p1].day > 1:
            p1 = p1+1
        start_month, start_year = x[p1].month, x[p1].year

        m1 = x[p1].year * 12 + x[p1].month - 1

        #Find the end of the last full month
        p2 = len(x) - 1
        while x[p2].day < tf.Last_Day_In_Month(x[p2].year,x[p2].month):
            p2 = p2-1
        end_month, end_year = x[p2].month, x[p2].year

        m2 = x[p2].year * 12 + x[p2].month - 1
        nyears = x[p2].year-x[p1].year + 1
        rows = []
        for m in range(m1,m2+1):
            yr = m//12
            mn = m+1-yr*12
            p1 , p2 = tf.first_last_in_month(x, mn, yr) 

            OutFile = SDC_Print([mn, '/',  yr, ':'], OutFile)
            #Calculate MSL from input heights
            MSL = 0.0
            npts = 0
            for i in range(p1,p2+1):
                 MSL = MSL + y[i]
                 npts = npts + 1
            MSL = MSL/npts

            #High Means and Highest Tide 
            MHHW = 0.0
            MHW = 0.0
            Mon_HWL = -99999.99
            nhighs = 0
            nhhighs = 0
            for i in range(len(highs)):
                if ((high_dts[i] >= x[p1]) and (high_dts[i] <= x[p2])):
                    if (high_types[i] == 'HH'):
                        MHHW = MHHW + high_values[i]
                        nhhighs = nhhighs + 1
                        MHW = MHW + high_values[i]
                        nhighs = nhighs + 1
                    if (high_types[i] == 'H'):
                        MHW = MHW + high_values[i]
                        nhighs = nhighs + 1
                    if (high_types[i] != 'H' and high_types[i] != 'HH'):
                        OutFile = SDC_Print(["Bad high type", high_types[i]], OutFile)

                    if high_values[i] > Mon_HWL:
                        Mon_HWL = high_values[i]

            MHHW = MHHW / nhhighs
            MHW  = MHW  / nhighs

            #Low Means and Lowest Tide 
            MLLW = 0.0
            MLW = 0.0
            Mon_LWL = 99999.99
            nlows = 0
            nllows = 0
            for i in range(len(lows)):
                if ((low_dts[i] >= x[p1]) and (low_dts[i] <= x[p2])):
                    if (low_types[i] == 'LL'):
                        MLLW = MLLW + low_values[i]
                        nllows = nllows + 1
                        MLW = MLW + low_values[i]
                        nlows = nlows + 1
                    if (low_types[i] == 'L'):
                        MLW = MLW + low_values[i]
                        nlows = nlows + 1
                    if (low_types[i] != 'L' and low_types[i] != 'LL'):
                        OutFile = SDC_Print(["Bad low type", low_types[i]], OutFile)

                    if low_values[i] < Mon_LWL:
                        Mon_LWL = low_values[i]

            MLLW = MLLW / nllows
            MLW  = MLW  / nlows

            OutFile = SDC_Print(['HWL  = ', fmt % Mon_HWL], OutFile)
            OutFile = SDC_Print(['MHHW = ', fmt % MHHW], OutFile)
            OutFile = SDC_Print(['MHW  = ', fmt % MHW], OutFile)
            OutFile = SDC_Print(['MSL  = ', fmt % MSL], OutFile)
            OutFile = SDC_Print(['MLW  = ', fmt % MLW], OutFile)
            OutFile = SDC_Print(['MLLW = ', fmt % MLLW], OutFile)
            OutFile = SDC_Print(['LWL  = ', fmt % Mon_LWL], OutFile)
            
            MM_Subordinate.append([Mon_HWL,MHHW,MHW, MSL, MLW, MLLW, Mon_LWL])

            row = {'time':datetime(yr,mn,1),
                   'HWL':Mon_HWL,
                   'MHHW':MHHW,
                   'MHW':MHW,
                   'MSL':MSL,
                   'MLW':MLW,
                   'MLLW':MLLW,
                   'LWL':Mon_LWL}
            rows.append(row)
        subordinate_monthly_means = pd.DataFrame(rows)


        ##################################################################
        #Calculate Means by MonthlyMeansSimultaneousComparison          #
        ##################################################################
    if Calc_Method == 'MMSC':
        OutFile = SDC_Print([' '], OutFile)
        OutFile = SDC_Print([' TIDAL DATUMS BY Monthly Means Simultaneous Comparison:'], OutFile)
        OutFile = SDC_Print([' '], OutFile)

        #Get Means for Control Station
        MM_Control = cd.Get_Monthly_Means(Control_Station_ID, start_month, start_year, end_month, end_year, CFactor)
        if len(MM_Control) == 0:
            OutFile = SDC_Print(['***Error*** No Monthly Means Returned for Control station: ', Control_Station_ID], OutFile)
            OutFile = SDC_Print(['Can not continue.'], OutFile)
            raise RuntimeError('No Monthly Means Returned for Control station: ' + str(Control_Station_ID) + '. Can not continue.')
        OutFile = SDC_Print([len(MM_Control), 'Months of control station means retrieved.'], OutFile)

        #Check That means tables are same size
        if len(MM_Subordinate) != len(MM_Control):
            OutFile = SDC_Print(['***Error*** Monthly means tables are different lengths!'], OutFile)
            OutFile = SDC_Print(['Sub Len= ', len(MM_Subordinate), '   Control Len = ', len(MM_Control)], OutFile)
            OutFile = SDC_Print([len(MM_Control)-len(MM_Subordinate), ' Missing months.'], OutFile)
            OutFile = SDC_Print(['Can not continue.'], OutFile)
            raise RuntimeError(('Monthly means tables are different lengths! ' +
                                'Sub Len= ' + str(len(MM_Subordinate)) + '   Control Len = ' + str(len(MM_Control)) + ' ' +
                                str(len(MM_Control)-len(MM_Subordinate)) + ' Missing months. ' +
                                'Can not continue.'))

        #Pass through control means and delete missing months from the analysis
        myear = start_year
        mmonth = start_month
        MM_S = []
        MM_C = []
        for i in range(len(MM_Control)):
            if not (-99999.99 in MM_Control[i]):
                MM_S.append(MM_Subordinate[i])
                MM_C.append(MM_Control[i])
            else:
                OutFile = SDC_Print(['Missing control means: ' + str(mmonth) + '/' + str(myear) + ' EXCLUDED from the analysis.'], OutFile)
            if mmonth == 12:
                mmonth = 1
                myear = myear + 1
            else:
                mmonth = mmonth + 1
        if (len(MM_C) < 1):
            OutFile = SDC_Print(['***Error*** No Monthly Means Returned for Control station: ', Control_Station_ID], OutFile)
            OutFile = SDC_Print(['Can not continue.'], OutFile)
            raise RuntimeError('No Monthly Means Returned for Control station: ' + str(Control_Station_ID) + '. Can not continue.')

        MM_Subordinate = MM_S
        MM_Control = MM_C

        myear = start_year
        mmonth = start_month
        for i in range(len(MM_Control)):

            if mmonth == 12:
                mmonth = 1
                myear = myear + 1
            else:
                mmonth = mmonth + 1

        nmonths = len(MM_Control)

        #Calculate Mean differences of monthly MSL, MTL, DTL, MN and GT, DHQ, DLQ
        Mean_Diff_MSL = 0.0
        Mean_Diff_MTL = 0.0
        Mean_Diff_DTL = 0.0
        Mean_Ratio_MN = 0.0
        Mean_Ratio_GT = 0.0
        Mean_Ratio_DHQ = 0.0
        Mean_Ratio_DLQ = 0.0
        Mean_Diff_MHW = 0.0
        Mean_Diff_MHHW = 0.0
        Mean_Diff_MLW = 0.0
        Mean_Diff_MLLW = 0.0

        for i in range(len(MM_Subordinate)):
            SMEANS = MM_Subordinate[i]
            CMEANS = MM_Control[i]
            Mean_Diff_MSL = Mean_Diff_MSL + SMEANS[3] - CMEANS[3]
            Sub_DTL = 0.5 * (SMEANS[1] + SMEANS[5])
            Con_DTL = 0.5 * (CMEANS[1] + CMEANS[5])
            Mean_Diff_DTL = Mean_Diff_DTL + Sub_DTL - Con_DTL
            Sub_MTL = 0.5 * (SMEANS[2] + SMEANS[4])
            Con_MTL = 0.5 * (CMEANS[2] + CMEANS[4])
            Mean_Diff_MTL = Mean_Diff_MTL + Sub_MTL - Con_MTL
            Sub_MN = SMEANS[2] - SMEANS[4]
            Con_MN = CMEANS[2] - CMEANS[4]
            Mean_Ratio_MN = Mean_Ratio_MN + (Sub_MN / Con_MN)
            Sub_GT = SMEANS[1] - SMEANS[5]
            Con_GT = CMEANS[1] - CMEANS[5]
            Mean_Ratio_GT = Mean_Ratio_GT + (Sub_GT / Con_GT)
            Sub_DHQ = SMEANS[1] - SMEANS[2]
            Con_DHQ = CMEANS[1] - CMEANS[2]
            Sub_DLQ = SMEANS[4] - SMEANS[5]
            Con_DLQ = CMEANS[4] - CMEANS[5]
            Mean_Diff_MHHW = Mean_Diff_MHHW + SMEANS[1] - CMEANS[1]
            Mean_Diff_MHW = Mean_Diff_MHW + SMEANS[2] - CMEANS[2]
            Mean_Diff_MLW = Mean_Diff_MLW + SMEANS[4] - CMEANS[4]
            Mean_Diff_MLLW = Mean_Diff_MLLW + SMEANS[5] - CMEANS[5]
            if Sub_Method == 'Standard':
                Mean_Ratio_DHQ = Mean_Ratio_DHQ + (Sub_DHQ / Con_DHQ)
                Mean_Ratio_DLQ = Mean_Ratio_DLQ + (Sub_DLQ / Con_DLQ)

        Mean_Diff_MSL = Mean_Diff_MSL / nmonths
        Mean_Diff_MTL = Mean_Diff_MTL / nmonths
        Mean_Diff_DTL = Mean_Diff_DTL / nmonths
        Mean_Ratio_MN = Mean_Ratio_MN / nmonths
        Mean_Ratio_GT = Mean_Ratio_GT / nmonths
        Mean_Diff_MHHW = Mean_Diff_MHHW / nmonths
        Mean_Diff_MHW = Mean_Diff_MHW / nmonths
        Mean_Diff_MLW = Mean_Diff_MLW / nmonths
        Mean_Diff_MLLW = Mean_Diff_MLLW / nmonths
        if Sub_Method == 'Standard':
            Mean_Ratio_DHQ = Mean_Ratio_DHQ / nmonths
            Mean_Ratio_DLQ = Mean_Ratio_DLQ / nmonths
        
        OutFile = SDC_Print([nmonths, ' months in the analysis\n'], OutFile)
        OutFile = SDC_Print(['Mean_Diff_MSL  = ', fmt % Mean_Diff_MSL ], OutFile)
        OutFile = SDC_Print(['Mean Diff MTL  = ', fmt % Mean_Diff_MTL ], OutFile)
        OutFile = SDC_Print(['Mean_Diff_DTL  = ', fmt % Mean_Diff_DTL ], OutFile)
        OutFile = SDC_Print(['Mean_Ratio_MN  = ', fmt % Mean_Ratio_MN ], OutFile)
        OutFile = SDC_Print(['Mean Ratio GT  = ', fmt % Mean_Ratio_GT], OutFile)
        OutFile = SDC_Print(['Mean_Diff_MHHW = ', fmt % Mean_Diff_MHHW], OutFile)
        OutFile = SDC_Print(['Mean_Diff_MHW  = ', fmt % Mean_Diff_MHW], OutFile)
        OutFile = SDC_Print(['Mean_Diff_MLW  = ', fmt % Mean_Diff_MLW], OutFile)
        OutFile = SDC_Print(['Mean_Diff_MLLW = ', fmt % Mean_Diff_MLLW], OutFile)
        if Sub_Method == 'Standard':
            OutFile = SDC_Print(['Mean Ratio DHQ = ', fmt % Mean_Ratio_DHQ], OutFile)
            OutFile = SDC_Print(['Mean Ratio DLQ = ', fmt % Mean_Ratio_DLQ], OutFile)

        #Add the Mean differences to correct the Control Datums to this station
        Subordinate_MSL = Control_MSL + Mean_Diff_MSL
        Subordinate_MTL = Control_MTL + Mean_Diff_MTL
        Subordinate_DTL = Control_DTL + Mean_Diff_DTL
        
        #Multlpy the Mean Control Datums by mean ratios to correct to this station
        Subordinate_MN = Control_MN * Mean_Ratio_MN
        Subordinate_GT = Control_GT * Mean_Ratio_GT
        Subordinate_MHHW = Control_Acc_Datums[0] + Mean_Diff_MHHW
        Subordinate_MHW = Control_Acc_Datums[1] + Mean_Diff_MHW
        Subordinate_MLLW = Control_Acc_Datums[6] + Mean_Diff_MLLW
        Subordinate_MLW = Control_Acc_Datums[5] + Mean_Diff_MLW
        if Sub_Method == 'Standard':
            Subordinate_DHQ = Control_Acc_Datums[9] * Mean_Ratio_DHQ
            Subordinate_DLQ = Control_Acc_Datums[10] * Mean_Ratio_DLQ

        OutFile = SDC_Print(['\n Corrected values for MN, GT, MTL, DTL'], OutFile)
        OutFile = SDC_Print([fmt % Subordinate_MN, fmt % Subordinate_GT, fmt % Subordinate_MTL, fmt % Subordinate_DTL], OutFile)
        if Sub_Method == 'Standard':
            OutFile = SDC_Print([' Corrected values for DHQ, DLQ'], OutFile)
            OutFile = SDC_Print([fmt % Subordinate_DHQ, fmt % Subordinate_DLQ], OutFile)
        OutFile = SDC_Print([' Corrected values for MHHW, MHW, MLW, MLLW'], OutFile)
        OutFile = SDC_Print([fmt % Subordinate_MHHW, fmt % Subordinate_MHW, 
                  fmt % Subordinate_MLW, fmt % Subordinate_MLLW], OutFile)

        #Save MHHW, MHW, MLW and MLLW for Direct computation
        Direct_MHHW = Subordinate_MHHW
        Direct_MHW = Subordinate_MHW
        Direct_MLW = Subordinate_MLW
        Direct_MLLW = Subordinate_MLLW

        if Sub_Method == 'Modified':
            #Calculate remaining subordinate datums with Modified Range Ratio Method
            Subordinate_MLW = Subordinate_MTL - (0.5 * Subordinate_MN)
            Subordinate_MHW = Subordinate_MLW + Subordinate_MN
            Subordinate_MLLW = Subordinate_DTL - (0.5 * Subordinate_GT)
            Subordinate_MHHW = Subordinate_MLLW + Subordinate_GT

            OutFile = SDC_Print(['\nDatums by Monthly Means Simultaneous Comparison (MMSC):'], OutFile)
            OutFile = SDC_Print(['HWL  = ', fmt % HWL, '  (' + HWL_DT.strftime("%Y/%m/%d %H:%M") + ')'], OutFile)
            OutFile = SDC_Print(['MHHW = ', fmt % Subordinate_MHHW], OutFile)
            OutFile = SDC_Print(['MHW  = ', fmt % Subordinate_MHW], OutFile)
            OutFile = SDC_Print(['DTL  = ', fmt % Subordinate_DTL], OutFile)
            OutFile = SDC_Print(['MTL  = ', fmt % Subordinate_MTL], OutFile)
            OutFile = SDC_Print(['MSL  = ', fmt % Subordinate_MSL], OutFile)
            OutFile = SDC_Print(['MLW  = ', fmt % Subordinate_MLW], OutFile)
            OutFile = SDC_Print(['MLLW = ', fmt % Subordinate_MLLW], OutFile)
            OutFile = SDC_Print(['DHQ  = ', fmt % (Subordinate_MHHW - Subordinate_MHW)], OutFile)
            OutFile = SDC_Print(['DLQ  = ', fmt % (Subordinate_MLW - Subordinate_MLLW)], OutFile)
            OutFile = SDC_Print(['GT   = ', fmt % Subordinate_GT], OutFile)
            OutFile = SDC_Print(['MN   = ', fmt % Subordinate_MN], OutFile)
            OutFile = SDC_Print(['LWL  = ', fmt % LWL, '  (' + LWL_DT.strftime("%Y/%m/%d %H:%M") + ')'], OutFile)

            datums = {'HWL':HWL,
                      'MHHW':Subordinate_MHHW,
                      'MHW':Subordinate_MHW,
                      'DTL':Subordinate_DTL,
                      'MTL':Subordinate_MTL,
                      'MSL':Subordinate_MSL,
                      'MLW':Subordinate_MLW,
                      'MLLW':Subordinate_MLLW,
                      'DHQ':Subordinate_MHHW - Subordinate_MHW,
                      'DLQ':Subordinate_MLW - Subordinate_MLLW,
                      'MN':Subordinate_MN,
                      'GT':Subordinate_GT,
                      'LWL':LWL}

     
        if Sub_Method == 'Standard':
            #Calculate remaining subordinate datums with Standard Method
            Subordinate_MLW = Subordinate_MTL - (0.5 * Subordinate_MN)
            Subordinate_MHW = Subordinate_MLW + Subordinate_MN
            Subordinate_MLLW = Subordinate_MLW - Subordinate_DLQ
            Subordinate_MHHW = Subordinate_MHW + Subordinate_DHQ

            OutFile = SDC_Print(['\nDatums by Monthly Means Simultaneous Comparison (MMSC):'], OutFile)
            OutFile = SDC_Print(['HWL  = ', fmt % HWL, '  (' + HWL_DT.strftime("%Y/%m/%d %H:%M") + ')'], OutFile)
            OutFile = SDC_Print(['MHHW = ', fmt % Subordinate_MHHW], OutFile)
            OutFile = SDC_Print(['MHW  = ', fmt % Subordinate_MHW], OutFile)
            OutFile = SDC_Print(['DTL  = ', fmt % Subordinate_DTL], OutFile)
            OutFile = SDC_Print(['MTL  = ', fmt % Subordinate_MTL], OutFile)
            OutFile = SDC_Print(['MSL  = ', fmt % Subordinate_MSL], OutFile)
            OutFile = SDC_Print(['MLW  = ', fmt % Subordinate_MLW], OutFile)
            OutFile = SDC_Print(['MLLW = ', fmt % Subordinate_MLLW], OutFile)
            OutFile = SDC_Print(['DHQ  = ', fmt % (Subordinate_MHHW - Subordinate_MHW)], OutFile)
            OutFile = SDC_Print(['DLQ  = ', fmt % (Subordinate_MLW - Subordinate_MLLW)], OutFile)
            OutFile = SDC_Print(['GT   = ', fmt % Subordinate_GT], OutFile)
            OutFile = SDC_Print(['MN   = ', fmt % Subordinate_MN], OutFile)
            OutFile = SDC_Print(['LWL  = ', fmt % LWL, '  (' + LWL_DT.strftime("%Y/%m/%d %H:%M") + ')'], OutFile)

            datums = {'HWL':HWL,
                      'MHHW':Subordinate_MHHW,
                      'MHW':Subordinate_MHW,
                      'DTL':Subordinate_DTL,
                      'MTL':Subordinate_MTL,
                      'MSL':Subordinate_MSL,
                      'MLW':Subordinate_MLW,
                      'MLLW':Subordinate_MLLW,
                      'DHQ':Subordinate_MHHW - Subordinate_MHW,
                      'DLQ':Subordinate_MLW - Subordinate_MLLW,
                      'MN':Subordinate_MN,
                      'GT':Subordinate_GT,
                      'LWL':LWL}

        if Sub_Method == 'Direct':
        # Datums with Direct Method
            OutFile = SDC_Print(['\nDatums by Monthly Means Simultaneous Comparison (MMSC):'], OutFile)
            OutFile = SDC_Print(['HWL  = ', fmt % HWL, '  (' + HWL_DT.strftime("%Y/%m/%d %H:%M") + ')'], OutFile)
            OutFile = SDC_Print(['MHHW = ', fmt % Direct_MHHW], OutFile)
            OutFile = SDC_Print(['MHW  = ', fmt % Direct_MHW], OutFile)
            OutFile = SDC_Print(['DTL  = ', fmt % (0.5 * (Direct_MHHW + Direct_MLLW))], OutFile)
            OutFile = SDC_Print(['MTL  = ', fmt % (0.5 * (Direct_MHW + Direct_MLW))], OutFile)
            OutFile = SDC_Print(['MSL  = ', fmt % Subordinate_MSL], OutFile)
            OutFile = SDC_Print(['MLW  = ', fmt % Direct_MLW], OutFile)
            OutFile = SDC_Print(['MLLW = ', fmt % Direct_MLLW], OutFile)
            OutFile = SDC_Print(['DHQ  = ', fmt % (Direct_MHHW - Direct_MHW)], OutFile)
            OutFile = SDC_Print(['DLQ  = ', fmt % (Direct_MLW - Direct_MLLW)], OutFile)
            OutFile = SDC_Print(['GT   = ', fmt % (Direct_MHHW - Direct_MLLW)], OutFile)
            OutFile = SDC_Print(['MN   = ', fmt % (Direct_MHW - Direct_MLW)], OutFile)
            OutFile = SDC_Print(['LWL  = ', fmt % LWL, '  (' + LWL_DT.strftime("%Y/%m/%d %H:%M") + ')'], OutFile)

            datums = {'HWL':HWL,
                      'MHHW':Direct_MHHW,
                      'MHW':Direct_MHW,
                      'DTL':0.5 * (Direct_MHHW + Direct_MLLW),
                      'MTL':0.5 * (Direct_MHW + Direct_MLW),
                      'MSL':Subordinate_MSL,
                      'MLW':Direct_MLW,
                      'MLLW':Direct_MLLW,
                      'DHQ':Direct_MHHW - Direct_MHW,
                      'DLQ':Direct_MLW - Direct_MLLW,
                      'MN':Direct_MHW - Direct_MLW,
                      'GT':Direct_MHHW - Direct_MLLW,
                      'LWL':LWL}


    if Calc_Method == 'TBYT':
        ######################################################################################
        #Calculate Datums using Tide-By-Tide                                                #
        ######################################################################################
        OutFile = SDC_Print([' '], OutFile)
        OutFile = SDC_Print(['Calculating Datums by TBYT:'], OutFile)

        #Get the Control Station Tides
        HL_Control = cd.Get_High_Lows(Control_Station_ID, x[0], x[len(x)-1], gmt_offset, CFactor)
        if len(HL_Control) == 0:
            OutFile = SDC_Print(['***Error*** No tides for control station: ', Control_Station_ID], OutFile)
            OutFile = SDC_Print(['Can not continue.'], OutFile)
            raise RuntimeError('No tides for control station: ' + str(Control_Station_ID) + '. Can not continue.')

        #Get the Subordinate Tides
        HL_Subordinate = []
        for i in range(len(high_lows)):
            HL_Subordinate.append([high_lows.iloc[i]['time'],high_lows.iloc[i]['value'],high_lows.iloc[i]['tide type']])

        #Calculate an expected time difference
        ETD = tf.Calc_Expected_Diff(HL_Subordinate, HL_Control)

        #Limit Magnitude of diff to 3 hrs
        if (abs(ETD) > 180.):
            if ETD < 0:
                ETD = -180.
            else:
                ETD = 180.

        OutFile = SDC_Print([' '], OutFile)
        OutFile = SDC_Print(['Estimated Time Difference: ', ETD, ' Minutes'], OutFile)
        OutFile = SDC_Print([' '], OutFile)

        Match_Window = 1  #hours for non-diurnal stations else 4
        if (tpd > 3.0):
            Match_Window = 1  #for non-diurnal stations else 4
        else:
            Match_Window = 4

        #Pair the tides
        Pairs = []
        for i in range(len(HL_Subordinate)):
            dt1 = HL_Subordinate[i][0]
            type1 = HL_Subordinate[i][2]
            val1 = HL_Subordinate[i][1]

            Win_Start = dt1 + timedelta(minutes=ETD) - timedelta(hours=Match_Window)
            Win_End = Win_Start + timedelta(hours=12)
            for j in range(len(HL_Control)):
                dt2 = HL_Control[j][0]
                type2 = HL_Control[j][2]
                val2 = HL_Control[j][1]
                if (dt2 >= Win_Start) and (dt2 <= Win_End) and type2[0] == type1[0]:
                    tdiff = (dt2-dt1).days*1440 + (dt2-dt1).seconds/60
                    Pairs.append([dt1, val1, dt2, val2, type2, tdiff])

                    break

        #Find means and mean differences for MHHW, MHW, MLW, MLLW
        MHHa = 0.0
        MLHa = 0.0
        MHLa = 0.0
        MLLa = 0.0
        DeltaHH = 0.0
        DeltaH  = 0.0
        DeltaL = 0.0
        DeltaLL = 0.0
        NHH = 0
        NH = 0
        NL = 0
        NLL = 0
        for i in range(len(Pairs)):
            PairType = Pairs[i][4]
            PairDiff = (Pairs[i][1] - Pairs[i][3])
            if PairType == 'HH':
                MHHa = MHHa + Pairs[i][1]
                DeltaHH = DeltaHH + PairDiff
                NHH = NHH + 1
            elif PairType == 'H':
                MLHa = MLHa + Pairs[i][1]
                DeltaH = DeltaH + PairDiff
                NH = NH + 1
            elif PairType == 'L':
                MHLa = MHLa + Pairs[i][1]
                DeltaL = DeltaL + PairDiff
                NL = NL + 1
            elif PairType == 'LL':
                MLLa = MLLa + Pairs[i][1]
                DeltaLL = DeltaLL + PairDiff
                NLL = NLL + 1

        if NHH == 0:
            OutFile = SDC_Print(['No Higher Highs were paired.'], OutFile)
        if NH == 0:
            OutFile = SDC_Print(['No Highs were paired.'], OutFile)
        if NL == 0:
            OutFile = SDC_Print(['No Lows were paired.'], OutFile)
        if NLL == 0:
            OutFile = SDC_Print(['No Lower-lows were paired.'], OutFile)

        if NHH == 0 or NH == 0 or NL == 0 or NLL == 0:
            OutFile = SDC_Print(['***Error*** Fatal issue. Exiting Analysis.'], OutFile)
            raise RuntimeError('Fatal issue. Exiting Analysis.')

        DeltaHH = DeltaHH / NHH
        DeltaH = DeltaH / NH
        DeltaL = DeltaL / NL
        DeltaLL = DeltaLL / NLL
        MHHa = MHHa / NHH
        MLHa = MLHa / NH
        MHLa = MHLa / NL
        MLLa = MLLa / NLL

        DHQa = 0.5 * (MHHa - MLHa)
        DLQa = 0.5 * (MHLa - MLLa)
        MHWa = 0.5 * (MHHa + MLHa)
        MLWa = 0.5 * (MHLa + MLLa)
        GTa  = MHHa - MLLa
        MNa  = MHWa - MLWa
        DTLa = 0.5 * (MHHa + MLLa)
        MTLa = 0.5 * (MHWa + MLWa)

        DeltaDHQ = 0.5 * (DeltaHH - DeltaH)
        DeltaDLQ = 0.5 * (DeltaL - DeltaLL)
        DeltaMHW = 0.5 * (DeltaHH + DeltaH)
        DeltaMLW = 0.5 * (DeltaL + DeltaLL)
        DeltaMN  = DeltaMHW - DeltaMLW
        DeltaDTL = 0.5 * (DeltaHH + DeltaLL)
        DeltaMTL = 0.5 * (DeltaMHW + DeltaMLW)
        DeltaGT  = DeltaHH - DeltaLL
        RatioMN  = MNa / (MNa - DeltaMN)
        RatioGT  = GTa / (GTa - DeltaGT)

        MTL  = Control_MTL + DeltaMTL
        DTL  = Control_DTL + DeltaDTL
        DHQ  = Control_DHQ + DeltaDHQ
        DLQ  = Control_DLQ + DeltaDLQ
        MN   = Control_MN * RatioMN
        GT   = Control_GT * RatioGT
     
        if Sub_Method == 'Modified':
            #Modified Range Ratio Method
            MLW  = MTL - 0.5 * MN
            MHW  = MLW + MN
            MLLW = DTL - 0.5 * GT
            MHHW = MLLW + GT

            OutFile = SDC_Print(['\nDatums by Tide By Tide Comparison (TBYT):'], OutFile)
            OutFile = SDC_Print(['HWL  = ', fmt % HWL, '  (' + HWL_DT.strftime("%Y/%m/%d %H:%M") + ')'], OutFile)
            OutFile = SDC_Print(['MHHW = ', fmt % MHHW], OutFile)
            OutFile = SDC_Print(['MHW  = ', fmt % MHW], OutFile)
            OutFile = SDC_Print(['DTL  = ', fmt % DTL], OutFile)
            OutFile = SDC_Print(['MTL  = ', fmt % MTL], OutFile)
            OutFile = SDC_Print(['MSL  = ', fmt % Mean_Level], OutFile)
            OutFile = SDC_Print(['MLW  = ', fmt % MLW], OutFile)
            OutFile = SDC_Print(['MLLW = ', fmt % MLLW], OutFile)
            OutFile = SDC_Print(['DHQ  = ', fmt % (MHHW - MHW)], OutFile)
            OutFile = SDC_Print(['DLQ  = ', fmt % (MLW - MLLW)], OutFile)
            OutFile = SDC_Print(['GT   = ', fmt % GT], OutFile)
            OutFile = SDC_Print(['MN   = ', fmt % MN], OutFile)
            OutFile = SDC_Print(['LWL  = ', fmt % LWL, '  (' + LWL_DT.strftime("%Y/%m/%d %H:%M") + ')'], OutFile)

            datums = {'HWL':HWL,
                      'MHHW':MHHW,
                      'MHW':MHW,
                      'DTL':DTL,
                      'MTL':MTL,
                      'MSL':Mean_Level,
                      'MLW':MLW,
                      'MLLW':MLLW,
                      'DHQ':MHHW - MHW,
                      'DLQ':MLW - MLLW,
                      'MN':MN,
                      'GT':GT,
                      'LWL':LWL}

        if Sub_Method == 'Standard':
            #Standard Method is selected, calculate subordinate station datums as follows  
            MLW  = MTL - 0.5 * MN
            MHW  = MLW + MN
            MLLW = MLW - DLQ
            MHHW = MHW + DHQ
            OutFile = SDC_Print(['\nDatums by Tide By Tide Comparison (TBYT):'], OutFile)
            OutFile = SDC_Print(['HWL  = ', fmt % HWL, '  (' + HWL_DT.strftime("%Y/%m/%d %H:%M") + ')'], OutFile)
            OutFile = SDC_Print(['MHHW = ', fmt % MHHW], OutFile)
            OutFile = SDC_Print(['MHW  = ', fmt % MHW], OutFile)
            OutFile = SDC_Print(['DTL  = ', fmt % DTL], OutFile)
            OutFile = SDC_Print(['MTL  = ', fmt % MTL], OutFile)
            OutFile = SDC_Print(['MSL  = ', fmt % Mean_Level], OutFile)
            OutFile = SDC_Print(['MLW  = ', fmt % MLW], OutFile)
            OutFile = SDC_Print(['MLLW = ', fmt % MLLW], OutFile)
            OutFile = SDC_Print(['DHQ  = ', fmt % (MHHW - MHW)], OutFile)
            OutFile = SDC_Print(['DLQ  = ', fmt % (MLW - MLLW)], OutFile)
            OutFile = SDC_Print(['GT   = ', fmt % GT], OutFile)
            OutFile = SDC_Print(['MN   = ', fmt % MN], OutFile)
            OutFile = SDC_Print(['LWL  = ', fmt % LWL, '  (' + LWL_DT.strftime("%Y/%m/%d %H:%M") + ')'], OutFile)

            datums = {'HWL':HWL,
                      'MHHW':MHHW,
                      'MHW':MHW,
                      'DTL':DTL,
                      'MTL':MTL,
                      'MSL':Mean_Level,
                      'MLW':MLW,
                      'MLLW':MLLW,
                      'DHQ':MHHW - MHW,
                      'DLQ':MLW - MLLW,
                      'MN':MN,
                      'GT':GT,
                      'LWL':LWL}

        if Sub_Method == 'Direct':
            #Direct Method is selected, calculate subordinate station datums as follows 
            MLW = Control_MLW + DeltaMLW
            MLLW = Control_MLLW + DeltaLL
            MHW = Control_MHW + DeltaMHW
            MHHW = Control_MHHW + DeltaHH
            OutFile = SDC_Print(['\nDatums by Tide By Tide Comparison (TBYT):'], OutFile)
            OutFile = SDC_Print(['HWL  = ', fmt % HWL, '  (' + HWL_DT.strftime("%Y/%m/%d %H:%M") + ')'], OutFile)
            OutFile = SDC_Print(['MHHW = ', fmt % MHHW], OutFile)
            OutFile = SDC_Print(['MHW  = ', fmt % MHW], OutFile)
            OutFile = SDC_Print(['DTL  = ', fmt % DTL], OutFile)
            OutFile = SDC_Print(['MTL  = ', fmt % MTL], OutFile)
            OutFile = SDC_Print(['MSL  = ', fmt % Mean_Level], OutFile)
            OutFile = SDC_Print(['MLW  = ', fmt % MLW], OutFile)
            OutFile = SDC_Print(['MLLW = ', fmt % MLLW], OutFile)
            OutFile = SDC_Print(['DHQ  = ', fmt % (MHHW - MHW)], OutFile)
            OutFile = SDC_Print(['DLQ  = ', fmt % (MLW - MLLW)], OutFile)
            OutFile = SDC_Print(['GT   = ', fmt % GT], OutFile)
            OutFile = SDC_Print(['MN   = ', fmt % MN], OutFile)
            OutFile = SDC_Print(['LWL  = ', fmt % LWL, '  (' + LWL_DT.strftime("%Y/%m/%d %H:%M") + ')'], OutFile)

            datums = {'HWL':HWL,
                      'MHHW':MHHW,
                      'MHW':MHW,
                      'DTL':DTL,
                      'MTL':MTL,
                      'MSL':Mean_Level,
                      'MLW':MLW,
                      'MLLW':MLLW,
                      'DHQ':MHHW - MHW,
                      'DLQ':MLW - MLLW,
                      'MN':MN,
                      'GT':GT,
                      'LWL':LWL}

    OutFile = SDC_Print(['\n', Units], OutFile)
    OutFile = SDC_Print(['\nThat is all.'], OutFile)
    
    
    out = Out(pd.DataFrame({'time':x,'val':y}), OutFile, out_plots, high_lows, subordinate_monthly_means, datums)
        
    return out


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run TADC from the command line')
    parser.add_argument('--fname',
                        type=str,
                        default=None,
                        help='The filespec of the timeseries data to be analyzed.')
    parser.add_argument('--data',
                        type=pd.core.frame.DataFrame,
                        default=None,
                        help='The timeseries data to be analyzed as a Pandas DataFrame object.')    
    parser.add_argument('--Pick_Method',
                        type=str,
                        default='PolyFit',
                        help="The method to use to pick high/low tides. (Default: 'PolyFit')")
    parser.add_argument('--Control_Station_ID',
                        type=int,
                        default=None,
                        help=("ID of NWLON control station to use.\n"
                              "Choose a control station or leave blank to default to First Reduction Datum (FRED) method.\n"
                              "If there is no control station Subordinate_Lon (the longitude of short-term station) is required for picking tide type.\n"
                              "(Default: 'None')"))
    parser.add_argument('--Method_Option',
                        type=str,
                        default='AUTO',
                        help=("Method to use ('AUTO', 'TBYT' or 'FRED').\nLess than one month of data, use Tide-By-Type (TBYT) comparison.\n"
                              "At least one complete month data, use Monthly Means comparison (AUTO).\n"
                              "In FRED method, datums are computed by averaging values over the observation time period.\n"
                              "If no control station ID is provided, the code chooses FRED method by default.\n"
                              "(Default: AUTO)"))
    parser.add_argument('--Time_Zone',
                        type=str,
                        default='GMT',
                        help=("The time zone.\n"
                              "Make sure the selected time zone is consistent with the data file uploaded.\n"
                              "The purpose of defining the time zone is to make sure appropriate time zone adjustment are applied to the control station.\n"
                              "String with the gmt offset at the end e.g. UST5.\n"
                              "(Default: 'GMT')"))
    parser.add_argument('--Units',
                        type=str,
                        default='Meters',
                        help=("The units.\n"
                              "The purpose of defining the unit is to make sure appropriate data unit conversion are applied to the control station.\n"
                              "Define units ('Meters', 'Centimeters', 'Millimeters', 'Feet', or 'Inches').\n"
                              "(Default: 'Meters')"))
    parser.add_argument('--Subordinate_Lat',
                        type=float,
                        default=None,
                        help=("The latitude of the subordinate station.\n"
                              "If no control station provided, subordinate station latitude is required.\n"
                              "Set this parameter to None if you are using control station.\n"
                              "(Default: 'None')"))
    parser.add_argument('--Subordinate_Lon',
                        type=float,
                        default=None,
                        help=("The longitude of the subordinate station.\n"
                              "Choosing Control_Station=None will allow users to compute tidal datum by averaging values of each\n"
                              "tide parameter over the observation time period using FRED.\n"
                              "If no control station provided, subordinate station longitude is required.\n"
                              "Set this parameter to None, if you are using control station.\n"
                              "(Default: 'None')"))
    parser.add_argument('--outfile_save_dir',
                        type=str,
                        default='None',
                        help=("The directory into which to save output files, if run from the command line.\n"
                              "If running from within a Python script, this has no effect as output files are not created.\n"
                              "If left as None, output files will not be saved.\n"
                              "(Default: None)"))
    parser.add_argument('--make_plots',
                        type=bool,
                        default=False,
                        help=("Whether or not to generate monthly plots. It is somewhat slow to generate the plots, \n"
                              "particularly for large input datasets."
                              "(Default: False)"))
    try:
        args = parser.parse_args()
    except SystemExit:
        sys.exit(1)

    out = run(fname=args.fname,
              Pick_Method=args.Pick_Method,
              Control_Station_ID=args.Control_Station_ID,
              Method_Option=args.Method_Option,
              Time_Zone=args.Time_Zone,
              Units=args.Units,
              Subordinate_Lat=args.Subordinate_Lat,
              Subordinate_Lon=args.Subordinate_Lon,
              outfile_save_dir=args.outfile_save_dir,
              make_plots=args.make_plots)
    print(out.readme)

    # Save files when run as a script #
    if args.outfile_save_dir != 'None':
        out_dir = os.path.join(args.outfile_save_dir, 'outfiles_' + datetime.now().strftime("%Y-%m-%d-%H%M%S"))
        os.mkdir(out_dir)
        with open(os.path.join(out_dir,'tadc.out'),'w') as f:
            f.write(out.readme)
        out.high_lows.to_csv(os.path.join(out_dir,'High_Lows.csv'),index=False)
        if args.make_plots:
            for i in range(len(out.plots)):
                out.plots.iloc[i]['plot'].savefig(os.path.join(out_dir,out.plots.iloc[i]['time']+'.png'),dpi=300)
    
        

              

