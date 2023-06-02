#!/usr/bin/python
"""
NOAA/CO-OPS Simple Datum Calculator Engine

@author: George Story

Convert to Python 3  February 2020

"""

import numpy as np
from numpy import diff, sign, mean
from scipy.signal import butter, filtfilt
import matplotlib
import matplotlib.dates as mdates
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from datetime import datetime, date, time, timedelta
import filter_defs as fd
import tides as tf
import control_data as cd
from dateutil.parser import parse
import sys
import configparser
from time import gmtime, strftime


def Get_Method(xtimes):
    #This function chooses method to use
    #If there is at least 1 complete month of data MMSC else TBYT
    i = 0
    while xtimes[i].day > 1 and i < (len(xtimes)-1):
    #Find the first month start
        i = i+1

    if i == len(xtimes)-1:
        print(' Could not find month start')
        return 'TBYT'
    ld = tf.Last_Day_In_Month(xtimes[i].year,xtimes[i].month)

    #See if the data continues to the end of the month
    while xtimes[i].day < ld and i < (len(xtimes)-1):
        i = i + 1

    if xtimes[i].day == ld:
        return("MMSC")
    else:
        return("TBYT")

#################################################################### 

def Fill_Gaps(x,y):
    #This function fills gaps
    global remaining_gaps
    loc = 0
    while loc < len(y)-1:
        #Find the start of the next gap
        while y[loc] != -99999.99 and loc < len(y)-1:
            loc = loc + 1
        if loc == len(y) -1:
            return(y)
        gap_start = loc
        
        #Find the end
        while y[loc] == -99999.99 and loc < len(y)-1:
            loc = loc + 1
        gap_end = loc - 1
        if x[gap_end] - x[gap_start] > timedelta(hours=3):
            remaining_gaps = remaining_gaps + 1
            print('Not Filling gap from:',x[gap_start], ' to ', x[gap_end] )
        else:
            print('Filling gap from:',x[gap_start], ' to ', x[gap_end] )
            gap_width = gap_end - gap_start + 1
            
            #Fill the gap    
            #Collect known values on either side of gap
            known_ys = []
            known_xs = []
            start_knowns = max(gap_start - max(gap_width,3), 0)
            end_knowns = min(gap_end + max(gap_width,3), len(x)-1)
            for l in range(start_knowns, end_knowns+1):
                if y[l] != -99999.99:
                    known_xs.append(l)
                    known_ys.append(y[l])
            
            #Calculate polynomial
            z = np.polyfit(known_xs, known_ys, 2)
            f = np.poly1d(z)

            #Calculate new values
            unknown_xs = []
            for i in range(gap_start, gap_end+1):
                unknown_xs.append(i)
            new_ys = f(unknown_xs)
            for i in range(gap_start, gap_end+1):
               y[i] = new_ys[i-gap_start]
    return y

#################################################################### 


def Longest_Segment(x,y):
    #This function trims dataset to just the longest continuous segment
    if len(x) != len(y):
        SDC_Print(['*** Error x and y are different lengths'])
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
    return(x2,y2)

#################################################################### 
def SDC_Print(PLines):
#This is a print function 

    OutLine = ''
    for i in range(len(PLines)):
        OutLine = OutLine + str(PLines[i]) + ' '
    OutFile.write(OutLine + '\n')
    print(OutLine)
    return()
  
#################################################################### 
#################################################################### 

"""This is the main routine  """

####################################################################################################
############### The following blocks are only valid for running the code without the web-based GUI.#
############### As is stated below, a configuration file named config.cfg is where #################
############### users can customize the inputs for the calculation. ################################
#################################################################### ###############################
#################################################################### ###############################
CONFIG_FILE = "config.cfg"

""" 
Read a configuration file section (denoted as [section]) and
store in a hash {}. Returns empty hash if no parameters found
Args:
  config_file (str): configuration file to read
  section (str): name of the section to read

Returns:
  List with configuration parameters from requested section

Raises:
  None
"""
def read_config(config_file, section):
    params = {}
    try:
        config = configparser.ConfigParser()
        with open(config_file) as f:
            #config.readfp(f)
            config.read_file(f)
            options = config.options(section)
            for option in options:            
                try:
                    params[option] = config.get(section, option)
                    if params[option] == -1:
                        print("Could not read option: %s" % option)                    
                except:
                    print("Exception reading option %s!" % option)
                    params[option] = None
    except configparser.NoSectionError as nse:
        print("No section %s found reading %s: %s", section, config_file, nse)
    except IOError as ioe:
        print("Config file not found: %s: %s", config_file, ioe)

    return params


# Obtain parameters info from config file #
dcalc_params = read_config(CONFIG_FILE, "par")
Pick_Method = dcalc_params['pick_method']
fname = dcalc_params['fname']
Control_Station_ID = dcalc_params['control_station']
Method_Option = dcalc_params['method_option']
Units = dcalc_params['units']
Time_Zone = dcalc_params['time_zone']
Subordinate_Lon = dcalc_params['subordinate_lon']
Subordinate_Lat = dcalc_params['subordinate_lat']


"""Process the command line arguments
SDC.py  filename  control  method  time_zone  units  datum  lat  lon
Where:
    filename = filespec of input to use
    control  = 7 character control station ID (or 'None')
    method  = "AUTO', 'TBYT' or 'FRED'
    time_zone = some tring with the gmt offset at the end e.g. UST5
    units = 'Meters, Centimeters, Millimeters, Feet, Inches
lat  = lattitude of station
lon  = longitude of station"""
remaining_gaps = 0

if len(sys.argv) > 1:
   fname = sys.argv[1]
if len(sys.argv) > 2:
    Control_Station_ID = sys.argv[2]
if len(sys.argv) > 3:
    Method_Option = sys.argv[3]
if len(sys.argv) > 4:
    Time_Zone = sys.argv[4]
if len(sys.argv) > 5:
    Units = sys.argv[5]
if len(sys.argv) > 6:
    Subordinate_Lat = sys.argv[6]
if len(sys.argv) > 7:
    Subordinate_Lon = sys.argv[7]

#get directory of the input filespec to use for the output files.
end_of_path = fname.rfind('/')
if end_of_path > -1:
    path = fname[0:end_of_path+1]
else:
    path = ''

#Open the output file
OutFile = open(path + 'SDC.out', 'w')

if Control_Station_ID == 'None' or len(Control_Station_ID) < 7:
    Method_Option = 'FRED'

#when the script fails to execute as a result of incorrect parameter submission, reason will be indicated in GUI as well as the logfile 
if fname == '':
    SDC_Print(['No filename on command line. Arguments:'])
    SDC_Print(['filename = filespec of input to use'])
    SDC_Print(['control  = 7 character control station ID (or None)'])
    SDC_Print(['method  = AUTO, TBYT or FRED'])
    SDC_Print(['time_zone = some tring with the gmt offset at the end e.g. UST5'])
    SDC_Print(['units = Meters, Centimeters, Millimeters, Feet, Inches'])
    SDC_Print(['lat  = latitude of station'])
    SDC_Print(['lon  = longitude of station'])

    OutFile.close
    exit(-1)
 
SDC_Print (['Run Time: ', strftime("%Y-%m-%d %H:%M:%S", gmtime())])
SDC_Print(['Using ', fname[end_of_path+1:]])
f = open(fname, 'r')

#Get time offset if subordinate is not gmt
print('Time Zone = ', Time_Zone)
hrstr = ''
n = len(Time_Zone)-1
while Time_Zone[n].isdigit():
    hrstr = hrstr + Time_Zone[n]
    n = n-1
if len(hrstr) > 0:
    gmt_offset = int(hrstr)
else:
    gmt_offset = 0

#Get Date-Times and Water Levels from csv file
dt = []
wl = []
MissingPoints = 0
lineno=1
for line in f:
    if len(line) > 12:
        try:
            comma=line.index(',')
            try:
                thedt = parse(line[0:comma])
                dt.append(thedt)
                try:
                    field_end = line[comma+1:].find(',')
                    if field_end >= 0:
                        field_end = comma + 1 + field_end
                    else:
                        field_end = len(line) - 1
                    if (line[comma+1:min(field_end, comma+1+3)].upper() == 'NAN'):
                        wl.append(-99999.99)
                        MissingPoints = MissingPoints + 1
                    else:
                        wl.append(float(line[comma+1:field_end]))
                except ValueError:
                    wl.append(-99999.99)
                    MissingPoints = MissingPoints + 1
            except ValueError:
                pass
        except ValueError:
            pass
    lineno=lineno+1

f.close

#Convert to Numpy Arrays
x = np.array(dt)
y = np.array(wl)

SDC_Print([len(x), ' data points loaded.'])

#Determine interval and check for consistency
Interval = x[1] - x[0]
SDC_Print(['Interval: ', Interval])
i = 1
ni = len(x)
#for i in range(len(x)-1):
while i < ni-3:
    if (x[i+1] - x[i]) > Interval and ((x[i+1] - int(x[i]).seconds % Interval.seconds) == 0):     
        SDC_Print(['Interval of ', x[i+1] - x[i], ' at ', x[i]])
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
        SDC_Print(['Interval of ', x[i+1] - x[i], ' at ', x[i]])
        SDC_Print(['***Error*** Time interval is inconsistent.'])
        exit(-1)
    else:
        i = i+1

ngaps = 0
if MissingPoints > 0:
    SDC_Print([MissingPoints, ' missing data points'])
    y = Fill_Gaps(x,y)

#It there are gaps that can not be filled (>3 hrs) then 
#just take longest string withou breaks
if (remaining_gaps > 0):
    x, y = Longest_Segment(x,y)
    SDC_Print([' '])
    SDC_Print(['Dataset trimmed to longest continuous segment.'])
    SDC_Print(['New Start:', x[0]])
    SDC_Print(['New End:  ', x[len(x)-1]])

if (x[len(x)-1] - x[0]) < timedelta(days=14):
    SDC_Print(['***Error*** Not enough data for analysis. 2 weeks minimum'])
    exit(-1)

#Check for input units.  Get conversion factor from meters
print(Units.upper())
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
    SDC_Print(['***Error*** Input units of', Units, 'not defined.'])
    OutFile.close
    exit(-1)

SDC_Print([''])
SDC_Print(['All calculations and results are in', Units])

#Determine calulation method MMSC or TBYT
#Use MMSC if there is a complete calendar month of data and a control station is selected
Calc_Method = Get_Method(x)
if Method_Option == 'TBYT':
    Calc_Method = 'TBYT'
if Method_Option == 'FRED':
    Calc_Method = 'FRED'

if Method_Option == 'AUTO':
    if Calc_Method == 'TBYT':
        SDC_Print(['Less than one month of data loaded.  Use Tide-By-Type comparison.'])
    else:
        SDC_Print(['At least one complete month loaded.  Use Monthly Means comparison.'])

if Calc_Method == 'FRED':
    #Get Sub-Method from subordinate station location
    try:
        flon = float(Subordinate_Lon)
    except:
        SDC_Print(['***Error***  Station Longitude is not a number:', Subordinate_Lon])
        exit(-1)
    if flon < -100.:
        Sub_Method = 'Standard'
    else:
        Sub_Method = 'Modified'
else:
    #Get Sub-Method from Control Station Location
    Sub_Method = cd.Get_SubMethod(Control_Station_ID)
    if (Sub_Method not in ['Standard', 'Modified']):
        SDC_Print([Sub_Method])
        exit(-1)
SDC_Print([''])
if Sub_Method == 'Standard':
    SDC_Print(['West coast/Pacific station:\n   Using Standard Range Ratio Method'])
else:
    SDC_Print(['Gulf/East coast station:\n   Using Modified Range Ratio Method'])
SDC_Print([''])


"""Set up Filter parameters."""
fs= 86400 / Interval.seconds 
order = 6
#cutoff = 5.0  #desired cutoff frequency of the filter, per day
cutoff = 4.0
SDC_Print(['Sampling Rate: ', fs, ' per day.   Using cutoff frequency of ', cutoff , ' per day'])
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
    SDC_Print(["***Warning*** - Tides are out of order"])
    OutFile.close

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
SDC_Print(['Data Start: ' , str(x[0])])
SDC_Print(['Data End  : ' , str(x[len(x)-1])])
SDC_Print(['Mean Water Level: ', fmt % Mean_Level])
SDC_Print(['Highest Water Level: ', fmt % HWL])
SDC_Print(['Lowest Water Level: ', fmt % LWL])
SDC_Print(['Duration: ', t])
SDC_Print(['High Tides Found: ', len(highs)])
SDC_Print(['Low Tides Found : ', len(lows)])
ndays = float(t.days) + (t.seconds + 360.)/86400.0
tpd = float(ntides)/float(ndays)
SDC_Print(['Tides per day: {0:.1f}'.format(tpd)])
if (tpd > 3.0):
    SDC_Print([('Semi-Diurnal - Using EXHL')])
    high_types, low_types = tf.EXHL(high_values, low_values)
else:
    SDC_Print([('Diurnal Using DIUR')])
    high_types, low_types = tf.DIUR(high_dts, high_values, low_dts, low_values, x[0])

#Summarize Extremes
SDC_Print([high_types.count('H'), ' Highs'])
SDC_Print([high_types.count('HH'), ' Higher Highs'])
SDC_Print([low_types.count('L'), ' Lows'])
SDC_Print([low_types.count('LL'), ' Lower Lows'])
SDC_Print([' '])

#Store the highs and lows in a file
f = open(path +  'High-Lows.csv', 'w')
li=0
hi=0
while ((hi < len(highs)) and (li < len(lows))):
    if (low_dts[li] < high_dts[hi]):
        f.write('{0:%Y-%m-%d %H:%M}, {1: f}, {2:s}\n'.format(low_dts[li], low_values[li], low_types[li]))
        li = li + 1
    else:
        f.write('{0:%Y-%m-%d %H:%M}, {1: f}, {2:s}\n'.format(high_dts[hi], high_values[hi], high_types[hi]))
        hi = hi + 1
else:
    if li >= (len(lows)):
        while (hi < len(highs)):
            f.write('{0:%Y-%m-%d %H:%M}, {1: f}, {2:s}\n'.format(high_dts[hi], high_values[hi], high_types[hi]))
            hi = hi + 1
    if hi >= (len(highs)):
        while (li < len(lows)):
            f.write('{0:%Y-%m-%d %H:%M}, {1: f}, {2:s}\n'.format(low_dts[li], low_values[li], low_types[li]))
            li = li + 1
f.close

#Generate Plots of wl and high-lows for each month
pn = 1
m1 = x[0].year * 12 + x[0].month - 1
m2 = x[len(x)-1].year * 12 + x[len(x)-1].month - 1
nyears = x[len(x)-1].year-x[0].year + 1

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

    plt.clf()
    plt.plot(x[p1:p2], y[p1:p2], 'b-', label='Water Level')
    plt.plot(MHHTimes, MHHighs, label = 'Higher Highs', marker='D', markersize=3, linestyle='None', color='r')
    plt.plot(MHTimes, MHighs, label = 'Highs', marker='o', markersize=3, linestyle='None', color='m')
    plt.plot(MLLTimes, MLLows, label = 'Lower Lows', marker='D', markersize=3, linestyle='None', color='r')
    plt.plot(MLTimes, MLows, label = 'Lows', marker='o', markersize=3, linestyle='None', color='m')

    plt.ylabel(Units)
    plt.grid()

    majorLocator = matplotlib.ticker.MultipleLocator(5)
    minorLocator = matplotlib.ticker.MultipleLocator(1)
    xax = plt.gca().get_xaxis() 
    xax.set_major_locator(majorLocator)
    #format major xtick label
    xax.set_major_formatter(mdates.DateFormatter('%m/%d/%y'))
    plt.savefig(path + 'Month' + str(pn))
    pn = pn + 1

SDC_Print([pn-1, ' Monthly plots generated\n'])

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
    SDC_Print([' '])
    SDC_Print([' TIDAL Datums by Arithmetic Mean of Your Data (First Reduction):'])
    SDC_Print(['HWL  = ', fmt % HWL, '  (' + HWL_DT.strftime("%Y/%m/%d %H:%M") + ')'])
    SDC_Print(['MHHW = ', fmt % MHHW])
    SDC_Print(['MHW  = ', fmt % MHW])
    SDC_Print(['DTL  = ', fmt % (0.5 * (MHHW + MLLW))])
    SDC_Print(['MTL  = ', fmt % (0.5 * (MHW + MLW))])
    SDC_Print(['MSL  = ', fmt % Mean_Level])
    SDC_Print(['MLW  = ', fmt % MLW])
    SDC_Print(['MLLW = ', fmt % MLLW])
    SDC_Print(['DHQ  = ', fmt % (MHHW - MHW)])
    SDC_Print(['DLQ  = ', fmt % (MLW - MLLW)])
    SDC_Print(['MN   = ', fmt % (MHW - MLW)])
    SDC_Print(['GT   = ', fmt % (MHHW - MLLW)])
    SDC_Print(['LWL  = ', fmt % LWL , '  (' + LWL_DT.strftime("%Y/%m/%d %H:%M") + ')'])
    SDC_Print([' '])


if Calc_Method == 'MMSC' or Calc_Method == 'TBYT':
#Get Accepted Datums for Control Station
    Control_Acc_Datums = cd.Get_Accepted_Datums(Control_Station_ID, CFactor)
    if (Control_Acc_Datums[0] == None or Control_Acc_Datums[1] == None or 
       Control_Acc_Datums[5] == None or Control_Acc_Datums[6] == None):
        SDC_Print(['***Error*** Problem retrieving Accepted Datums for station ', Control_Station_ID])
        OutFile.close
        exit(-1)
    SDC_Print(['Control Datums for: ' , Control_Station_ID])
    SDC_Print(['\nMHHW,  MHW,  DTL,  MTL,  MSL,  MLW,  MLLW'])
    MeanString = ''
    for di in range(0,7):
        MeanString = MeanString + fmt % Control_Acc_Datums[di] + ' '
    SDC_Print([MeanString])

    SDC_Print(['GT,  MN,  DHQ,  DLQ,  NAVD,  LWI,  HWI'])
    MeanString = ''
    for di in range(7,14):
        if (Control_Acc_Datums[di] == None):
            MeanString = MeanString + 'Null '
        else:
            MeanString = MeanString + fmt % Control_Acc_Datums[di] + ' '
    SDC_Print([MeanString])

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
SDC_Print(['\nSUBORDINATE MONTHLY MEANS:'])
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
for m in range(m1,m2+1):
    yr = m//12
    mn = m+1-yr*12
    p1 , p2 = tf.first_last_in_month(x, mn, yr) 

    SDC_Print([mn, '/',  yr, ':'])
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
                SDC_Print(["Bad high type", high_types[i]])

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
                SDC_Print(["Bad low type", low_types[i]])

            if low_values[i] < Mon_LWL:
                Mon_LWL = low_values[i]

    MLLW = MLLW / nllows
    MLW  = MLW  / nlows

    SDC_Print(['HWL  = ', fmt % Mon_HWL])
    SDC_Print(['MHHW = ', fmt % MHHW])
    SDC_Print(['MHW  = ', fmt % MHW])
    SDC_Print(['MSL  = ', fmt % MSL])
    SDC_Print(['MLW  = ', fmt % MLW])
    SDC_Print(['MLLW = ', fmt % MLLW])
    SDC_Print(['LWL  = ', fmt % Mon_LWL])
    
    MM_Subordinate.append([Mon_HWL,MHHW,MHW, MSL, MLW, MLLW, Mon_LWL])


    ##################################################################
    #Calculate Means by MonthlyMeansSimultaneousComparison          #
    ##################################################################
if Calc_Method == 'MMSC':
    SDC_Print([' '])
    SDC_Print([' TIDAL DATUMS BY Monthly Means Simultaneous Comparison:'])
    SDC_Print([' '])

    #Get Means for Control Station
    MM_Control = cd.Get_Monthly_Means(Control_Station_ID, start_month, start_year, end_month, end_year, CFactor)
    if len(MM_Control) == 0:
        SDC_Print(['***Error*** No Monthly Means Returned for Control station: ', Control_Station_ID])
        SDC_Print(['Can not continue.'])
        OutFile.close
        exit(-1)
    SDC_Print([len(MM_Control), 'Months of control station means retrieved.'])

    #Check That means tables are same size
    if len(MM_Subordinate) != len(MM_Control):
        SDC_Print(['***Error*** Monthly means tables are different lengths!'])
        SDC_Print(['Sub Len= ', len(MM_Subordinate), '   Control Len = ', len(MM_Control)])
        SDC_Print([len(MM_Control)-len(MM_Subordinate), ' Missing months.'])
        SDC_Print(['Can not continue.'])
        OutFile.close
        exit(-1)

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
            SDC_Print(['Missing control means: ' + str(mmonth) + '/' + str(myear) + ' EXCLUDED from the analysis.'])
        if mmonth == 12:
            mmonth = 1
            myear = myear + 1
        else:
            mmonth = mmonth + 1
    if (len(MM_C) < 1):
        SDC_Print(['***Error*** No Monthly Means Returned for Control station: ', Control_Station_ID])
        SDC_Print(['Can not continue.'])
        OutFile.close
        exit(-1)
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
    
    SDC_Print([nmonths, ' months in the analysis\n'])
    SDC_Print(['Mean_Diff_MSL  = ', fmt % Mean_Diff_MSL ])
    SDC_Print(['Mean Diff MTL  = ', fmt % Mean_Diff_MTL ])
    SDC_Print(['Mean_Diff_DTL  = ', fmt % Mean_Diff_DTL ])
    SDC_Print(['Mean_Ratio_MN  = ', fmt % Mean_Ratio_MN ])
    SDC_Print(['Mean Ratio GT  = ', fmt % Mean_Ratio_GT])
    SDC_Print(['Mean_Diff_MHHW = ', fmt % Mean_Diff_MHHW])
    SDC_Print(['Mean_Diff_MHW  = ', fmt % Mean_Diff_MHW])
    SDC_Print(['Mean_Diff_MLW  = ', fmt % Mean_Diff_MLW])
    SDC_Print(['Mean_Diff_MLLW = ', fmt % Mean_Diff_MLLW])
    if Sub_Method == 'Standard':
        SDC_Print(['Mean Ratio DHQ = ', fmt % Mean_Ratio_DHQ])
        SDC_Print(['Mean Ratio DLQ = ', fmt % Mean_Ratio_DLQ])

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

    SDC_Print(['\n Corrected values for MN, GT, MTL, DTL'])
    SDC_Print([fmt % Subordinate_MN, fmt % Subordinate_GT, fmt % Subordinate_MTL, fmt % Subordinate_DTL])
    if Sub_Method == 'Standard':
        SDC_Print([' Corrected values for DHQ, DLQ'])
        SDC_Print([fmt % Subordinate_DHQ, fmt % Subordinate_DLQ])
    SDC_Print([' Corrected values for MHHW, MHW, MLW, MLLW'])
    SDC_Print([fmt % Subordinate_MHHW, fmt % Subordinate_MHW, 
              fmt % Subordinate_MLW, fmt % Subordinate_MLLW])

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

        SDC_Print(['\nDatums by Monthly Means Simultaneous Comparison (MMSC):'])
        SDC_Print(['HWL  = ', fmt % HWL, '  (' + HWL_DT.strftime("%Y/%m/%d %H:%M") + ')'])
        SDC_Print(['MHHW = ', fmt % Subordinate_MHHW])
        SDC_Print(['MHW  = ', fmt % Subordinate_MHW])
        SDC_Print(['DTL  = ', fmt % Subordinate_DTL])
        SDC_Print(['MTL  = ', fmt % Subordinate_MTL])
        SDC_Print(['MSL  = ', fmt % Subordinate_MSL])
        SDC_Print(['MLW  = ', fmt % Subordinate_MLW])
        SDC_Print(['MLLW = ', fmt % Subordinate_MLLW])
        SDC_Print(['DHQ  = ', fmt % (Subordinate_MHHW - Subordinate_MHW)])
        SDC_Print(['DLQ  = ', fmt % (Subordinate_MLW - Subordinate_MLLW)])
        SDC_Print(['GT   = ', fmt % Subordinate_GT])
        SDC_Print(['MN   = ', fmt % Subordinate_MN])
        SDC_Print(['LWL  = ', fmt % LWL, '  (' + LWL_DT.strftime("%Y/%m/%d %H:%M") + ')'])

 
    if Sub_Method == 'Standard':
        #Calculate remaining subordinate datums with Standard Method
        Subordinate_MLW = Subordinate_MTL - (0.5 * Subordinate_MN)
        Subordinate_MHW = Subordinate_MLW + Subordinate_MN
        Subordinate_MLLW = Subordinate_MLW - Subordinate_DLQ
        Subordinate_MHHW = Subordinate_MHW + Subordinate_DHQ

        SDC_Print(['\nDatums by Monthly Means Simultaneous Comparison (MMSC):'])
        SDC_Print(['HWL  = ', fmt % HWL, '  (' + HWL_DT.strftime("%Y/%m/%d %H:%M") + ')'])
        SDC_Print(['MHHW = ', fmt % Subordinate_MHHW])
        SDC_Print(['MHW  = ', fmt % Subordinate_MHW])
        SDC_Print(['DTL  = ', fmt % Subordinate_DTL])
        SDC_Print(['MTL  = ', fmt % Subordinate_MTL])
        SDC_Print(['MSL  = ', fmt % Subordinate_MSL])
        SDC_Print(['MLW  = ', fmt % Subordinate_MLW])
        SDC_Print(['MLLW = ', fmt % Subordinate_MLLW])
        SDC_Print(['DHQ  = ', fmt % (Subordinate_MHHW - Subordinate_MHW)])
        SDC_Print(['DLQ  = ', fmt % (Subordinate_MLW - Subordinate_MLLW)])
        SDC_Print(['GT   = ', fmt % Subordinate_GT])
        SDC_Print(['MN   = ', fmt % Subordinate_MN])
        SDC_Print(['LWL  = ', fmt % LWL, '  (' + LWL_DT.strftime("%Y/%m/%d %H:%M") + ')'])

    if Sub_Method == 'Direct':
    # Datums with Direct Method
        SDC_Print(['\nDatums by Monthly Means Simultaneous Comparison (MMSC):'])
        SDC_Print(['HWL  = ', fmt % HWL, '  (' + HWL_DT.strftime("%Y/%m/%d %H:%M") + ')'])
        SDC_Print(['MHHW = ', fmt % Direct_MHHW])
        SDC_Print(['MHW  = ', fmt % Direct_MHW])
        SDC_Print(['DTL  = ', fmt % (0.5 * (Direct_MHHW + Direct_MLLW))])
        SDC_Print(['MTL  = ', fmt % (0.5 * (Direct_MHW + Direct_MLW))])
        SDC_Print(['MSL  = ', fmt % Subordinate_MSL])
        SDC_Print(['MLW  = ', fmt % Direct_MLW])
        SDC_Print(['MLLW = ', fmt % Direct_MLLW])
        SDC_Print(['DHQ  = ', fmt % (Direct_MHHW - Direct_MHW)])
        SDC_Print(['DLQ  = ', fmt % (Direct_MLW - Direct_MLLW)])
        SDC_Print(['GT   = ', fmt % (Direct_MHHW - Direct_MLLW)])
        SDC_Print(['MN   = ', fmt % (Direct_MHW - Direct_MLW)])
        SDC_Print(['LWL  = ', fmt % LWL, '  (' + LWL_DT.strftime("%Y/%m/%d %H:%M") + ')'])


if Calc_Method == 'TBYT':
    ######################################################################################
    #Calculate Datums using Tide-By-Tide                                                #
    ######################################################################################
    SDC_Print([' '])
    SDC_Print(['Calculating Datums by TBYT:'])

    #Get the Control Station Tides
    HL_Control = cd.Get_High_Lows(Control_Station_ID, x[0], x[len(x)-1], gmt_offset, CFactor)
    if len(HL_Control) == 0:
        SDC_Print(['***Error*** No tides for control station: ', Control_Station_ID])
        SDC_Print(['Can not continue.'])
        OutFile.close
        exit(-1)

    #Get the Subordinate Tides
    HL_Subordinate = []
    f = open(path + 'High-Lows.csv', 'r')
    for line in f:
        Vals = line.split(',')
        HL_Subordinate.append([datetime.strptime(Vals[0], "%Y-%m-%d %H:%M"), float(Vals[1]), Vals[2].strip()])
    f.close

    #Calculate an expected time difference
    ETD = tf.Calc_Expected_Diff(HL_Subordinate, HL_Control)

    #Limit Magnitude of diff to 3 hrs
    if (abs(ETD) > 180.):
        if ETD < 0:
            ETD = -180.
        else:
            ETD = 180.

    SDC_Print([' '])
    SDC_Print(['Estimated Time Difference: ', ETD, ' Minutes'])
    SDC_Print([' '])

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
        SDC_Print(['No Higher Highs were paired.'])
    if NH == 0:
        SDC_Print(['No Highs were paired.'])
    if NL == 0:
        SDC_Print(['No Lows were paired.'])
    if NLL == 0:
        SDC_Print(['No Lower-lows were paired.'])

    if NHH == 0 or NH == 0 or NL == 0 or NLL == 0:
        SDC_Print(['***Error*** Fatal issue. Exiting Analysis.'])
        exit(-1)

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

        SDC_Print(['\nDatums by Tide By Tide Comparison (TBYT):'])
        SDC_Print(['HWL  = ', fmt % HWL, '  (' + HWL_DT.strftime("%Y/%m/%d %H:%M") + ')'])
        SDC_Print(['MHHW = ', fmt % MHHW])
        SDC_Print(['MHW  = ', fmt % MHW])
        SDC_Print(['DTL  = ', fmt % DTL])
        SDC_Print(['MTL  = ', fmt % MTL])
        SDC_Print(['MSL  = ', fmt % Mean_Level])
        SDC_Print(['MLW  = ', fmt % MLW])
        SDC_Print(['MLLW = ', fmt % MLLW])
        SDC_Print(['DHQ  = ', fmt % (MHHW - MHW)])
        SDC_Print(['DLQ  = ', fmt % (MLW - MLLW)])
        SDC_Print(['GT   = ', fmt % GT])
        SDC_Print(['MN   = ', fmt % MN])
        SDC_Print(['LWL  = ', fmt % LWL, '  (' + LWL_DT.strftime("%Y/%m/%d %H:%M") + ')'])

    if Sub_Method == 'Standard':
        #Standard Method is selected, calculate subordinate station datums as follows  
        MLW  = MTL - 0.5 * MN
        MHW  = MLW + MN
        MLLW = MLW - DLQ
        MHHW = MHW + DHQ
        SDC_Print(['\nDatums by Tide By Tide Comparison (TBYT):'])
        SDC_Print(['HWL  = ', fmt % HWL, '  (' + HWL_DT.strftime("%Y/%m/%d %H:%M") + ')'])
        SDC_Print(['MHHW = ', fmt % MHHW])
        SDC_Print(['MHW  = ', fmt % MHW])
        SDC_Print(['DTL  = ', fmt % DTL])
        SDC_Print(['MTL  = ', fmt % MTL])
        SDC_Print(['MSL  = ', fmt % Mean_Level])
        SDC_Print(['MLW  = ', fmt % MLW])
        SDC_Print(['MLLW = ', fmt % MLLW])
        SDC_Print(['DHQ  = ', fmt % (MHHW - MHW)])
        SDC_Print(['DLQ  = ', fmt % (MLW - MLLW)])
        SDC_Print(['GT   = ', fmt % GT])
        SDC_Print(['MN   = ', fmt % MN])
        SDC_Print(['LWL  = ', fmt % LWL, '  (' + LWL_DT.strftime("%Y/%m/%d %H:%M") + ')'])

    if Sub_Method == 'Direct':
        #Direct Method is selected, calculate subordinate station datums as follows 
        MLW = Control_MLW + DeltaMLW
        MLLW = Control_MLLW + DeltaLL
        MHW = Control_MHW + DeltaMHW
        MHHW = Control_MHHW + DeltaHH
        SDC_Print(['\nDatums by Tide By Tide Comparison (TBYT):'])
        SDC_Print(['HWL  = ', fmt % HWL, '  (' + HWL_DT.strftime("%Y/%m/%d %H:%M") + ')'])
        SDC_Print(['MHHW = ', fmt % MHHW])
        SDC_Print(['MHW  = ', fmt % MHW])
        SDC_Print(['DTL  = ', fmt % DTL])
        SDC_Print(['MTL  = ', fmt % MTL])
        SDC_Print(['MSL  = ', fmt % Mean_Level])
        SDC_Print(['MLW  = ', fmt % MLW])
        SDC_Print(['MLLW = ', fmt % MLLW])
        SDC_Print(['DHQ  = ', fmt % (MHHW - MHW)])
        SDC_Print(['DLQ  = ', fmt % (MLW - MLLW)])
        SDC_Print(['GT   = ', fmt % GT])
        SDC_Print(['MN   = ', fmt % MN])
        SDC_Print(['LWL  = ', fmt % LWL, '  (' + LWL_DT.strftime("%Y/%m/%d %H:%M") + ')'])

SDC_Print(['\n', Units])
SDC_Print(['\nThat is all.'])
OutFile.close
exit(1)

