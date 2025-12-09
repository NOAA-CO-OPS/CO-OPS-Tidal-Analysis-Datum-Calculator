"""Tide analysis functions for CO-OPS Datums Calculator"""

from datetime import datetime, date, time, timedelta
import logging
import numpy as np
logger = logging.getLogger(__name__)


def first_last_in_month(dates, month, year):
    #This function return the indexes of the first and last dates within month
    p = 0
    m = dates[p].month
    y = dates[p].year
    #find first point in month
    while not (y==year and m==month):
        p = p+1
        m = dates[p].month
        y = dates[p].year
    p1=p
    while (y==year and m==month and p < (len(dates)-1)):
        p = p+1
        m = dates[p].month
        y = dates[p].year
    if (p==(len(dates)-1)):
        p2=len(dates)-1
    else:
        p2= p-1
    return p1, p2


def Check_Tide_Order(dt, h, l):
    #This function checks that tides are in High-Low order and 
    #merges the highs and lows into a single time-ordered list
    hi = 0
    li = 0
    tides = []
    tide_types = []
    tide_indexes = []
    while ((hi<len(h)) and (li<len(l))):
        if (dt[h[hi]] <= dt[l[li]]):
            tides.append(h[hi])
            tide_types.append('H')
            tide_indexes.append(hi)
            hi = hi + 1
        else:
            tides.append(l[li])
            tide_types.append('L')
            tide_indexes.append(li)
            li = li + 1
    while hi < len(h):
        tides.append(h[hi])
        tide_types.append('H')
        tide_indexes.append(hi)
        hi = hi + 1
    while li < len(l):
        tides.append(l[li])
        tide_types.append('L')
        tide_indexes.append(li)
        li = li + 1

    ttype = tide_types[0]
    for i in range(1,len(tides)-1):
        if tide_types[i] == ttype:
            logger.warning('Tides are out of order at:', dt[tides[i]])
            return -1
        ttype = tide_types[i]
        i = i+1

    return 1


def EXHL(hvals, lvals):
    #This function chooses and flags highs and lows for semidiurnal/mixed tide types
    htypes = []
    for i in range(len(hvals)):
        htypes.append('H')
    for i in range(0, len(hvals)-1,2):
        if (hvals[i] >= hvals[i+1]):
            htypes[i] = 'HH'
        else:
            htypes[i+1] = 'HH'

    ltypes = []
    for i in range(len(lvals)):
        ltypes.append('L')
    for i in range(0, len(lvals)-1,2):
        if (lvals[i] <= lvals[i+1]):
            ltypes[i] = 'LL'
        else:
            ltypes[i+1] = 'LL'

    return htypes, ltypes


def Highest(h_dts, h_vals, t1, t2):
    #This function returns the index of the highest value between t1 and t2
    mxindex = -1
    mxval = -99999.99
    for i in range(len(h_dts)):
        if ((h_dts[i] >= t1) and (h_dts[i] <= t2)):
            if (h_vals[i] > mxval):
               mxval = h_vals[i]
               mxindex = i
    return mxindex


def Lowest(l_dts, l_vals, t1, t2):
    #This function returns the index of the lowest value between t1 and t2
    mxindex = -1
    minval = 99999.99
    for i in range(len(l_dts)):
        if ((l_dts[i] >= t1) and (l_dts[i] <= t2)):
            if (l_vals[i] < minval):
               mxval = l_vals[i]
               mxindex = i
    return mxindex


def Nearest_Tide(t_dts, dt):
    #This function picks and returns the index of the tide nearest
    t = t_dts[0]
    i=0
    while (t_dts[i] < dt and i<(len(t_dts)-1)):
        i = i+1
    if ((dt-t_dts[i-1]) < (t_dts[i]-dt)):
        i = i-1
    return i

    
def DIUR(h_dts, h_vals, l_dts, l_vals, t0):
    #This function chooses and flags highs and lows for diurnal tide types
    """Note:  DIUR selects:  point1 tide 1  point1  check for extreme tide
                                timeA    timeC   between point1 and timeC
                                #point2         tide 2
                                #timeB
                                #point3
            Point1 is Highest/Lowest tide in first 13/25 hours of month (or 1st tide).
            Point2 is tide nearest 25 hours beyond point1.
            Point3 is tide nearest 25 hours beyond point2.
            
            TimeA is the midpoint between point1 and point2.
            TimeB is the midpoint between point2 and point3.

            Tide 1 is point1.
            Tide 2 is the more extreme point1 or tide between timeA  and timeB.
            TimeC is the midpoint between tide1 and tide2.

            Tide 1 is the more extreme point1 or tide between point1 and timeC.

            Flag the extreme tide between tide1 and timeC.

            Subsequent tides are found by shifting point2 to point1 (tide2
            to tide1), and then repeating the process from Step 2 (below). """

    htypes = []
    for i in range(len(h_dts)):
        htypes.append('H')
#Highest tide in first 13 hours, else first tide
    tide = Highest(h_dts, h_vals, t0, t0 + timedelta(hours=13))
    if (tide < 0):
        tide = 0
    Point1 = h_dts[tide]
    Tide1 = tide
    while (Tide1 < (len(h_dts)-1)):
        tide = Nearest_Tide(h_dts, Point1+timedelta(hours=25))
        if h_dts[tide] == Point1:
            tide = tide+1
        Point2 = h_dts[tide]
        tide = Nearest_Tide(h_dts, Point2+timedelta(hours=25))
        if h_dts[tide] == Point2:
            tide = min(tide+1, len(h_dts)-1)
        Point3 = h_dts[tide]
        TimeA = Point1 + (Point2-Point1)/2
        TimeB = Point2 + (Point3-Point2)/2
        Tide2 = Highest(h_dts, h_vals, TimeA, TimeB)
        TimeC = h_dts[Tide1] + (h_dts[Tide2] - h_dts[Tide1])/2
        Tide1 = Highest(h_dts, h_vals, Point1, TimeC)
        extide =  Highest(h_dts, h_vals, h_dts[Tide1], TimeC)
        htypes[extide] = 'HH'
        Point1 = Point2
        Tide1 = Tide2


    ltypes = []
    for i in range(len(l_dts)):
        ltypes.append('L')
        #Lowest tide in first 13 hours, else first tide
    tide = Lowest(l_dts, l_vals, t0, t0 + timedelta(hours=13))
    if (tide < 0):
        tide = 0
    Point1 = l_dts[tide]
    Tide1 = tide
    while (Tide1 < (len(l_dts)-1)):
        tide = Nearest_Tide(l_dts, Point1+timedelta(hours=25))
        if l_dts[tide] == Point1:
            tide = tide+1
        Point2 = l_dts[tide]
        tide = Nearest_Tide(l_dts, Point2+timedelta(hours=25))
        if l_dts[tide] == Point2:
            tide = min(tide+1, len(l_dts)-1)
        Point3 = l_dts[tide]
        TimeA = Point1 + (Point2-Point1)/2
        TimeB = Point2 + (Point3-Point2)/2
        Tide2 = Lowest(l_dts, l_vals, TimeA, TimeB)
        TimeC = l_dts[Tide1] + (l_dts[Tide2] - l_dts[Tide1])/2
        Tide1 = Lowest(l_dts, l_vals, Point1, TimeC)
        extide =  Lowest(l_dts, l_vals, l_dts[Tide1], TimeC)
        ltypes[extide] = 'LL'
        Point1 = Point2
        Tide1 = Tide2

    return htypes, ltypes


def Check_Tides(dt, wl, h, l, Units_Factor):
    #This function Checks tides for minimum time and height between neighbors,
    #selects tides satisfying criteria and  
    #merges the highs and lows into a single time-ordered list

    Min_Height_Diff = 0.03 * Units_Factor
    hi = 0
    li = 0
    tides = []
    tide_types = []
    tide_indexes = []
    while ((hi<len(h)) and (li<len(l))):
        if (dt[h[hi]] <= dt[l[li]]):
            tides.append(h[hi])
            tide_types.append('H')
            tide_indexes.append(hi)
            hi = hi + 1
        else:
            tides.append(l[li])
            tide_types.append('L')
            tide_indexes.append(li)
            li = li + 1
    while hi < len(h):
        tides.append(h[hi])
        tide_types.append('H')
        tide_indexes.append(hi)
        hi = hi + 1
    while li < len(l):
        tides.append(l[li])
        tide_types.append('L')
        tide_indexes.append(li)
        li = li + 1

    hi_mask = np.ones(len(h), dtype=bool)
    lo_mask = np.ones(len(l), dtype=bool)

    t1 = 0
    t2 = 1
    aredeletedtides = 0
    while (t2<len(tides)):
        #/* Walk through the tides and mark offending pairs for deletion */

        if(((dt[tides[t2]] - dt[tides[t1]]) < timedelta(hours=2)) and (tide_types[t1] == tide_types[t2])):
            #check if time differnce is less than threshold and if tide types are the same, delete second tide 
            if tide_types[t2] == "H":
                hi_mask[tide_indexes[t2]] = False
            else:
                lo_mask[tide_indexes[t2]] = False
            t2=t2+1
            aredeletedtides = 1
        elif (((dt[tides[t2]] - dt[tides[t1]]) < timedelta(hours=2)) or (abs(wl[tides[t2]]-wl[tides[t1]]) < Min_Height_Diff)) and \
                   (tide_types[t1] != tide_types[t2]):
            #check if tide height/time difference is less than threshold and if tide type is not the same, delete both tides 
            if tide_types[t1] == "H":
                hi_mask[tide_indexes[t1]] = False
            else:
                lo_mask[tide_indexes[t1]] = False
            if tide_types[t2] == "H":
                hi_mask[tide_indexes[t2]] = False
            else:
                lo_mask[tide_indexes[t2]] = False
            t1=t2+1
            t2=t1+1
        else:
            t1=t2
            t2=t1+1
    return hi_mask, lo_mask


def Last_Day_In_Month(y,m):
    #This function determines last day of a month 
    if (m == 1 or m == 3 or m == 5 or m == 7 or m == 8 or m == 10 or m == 12):
        Last_Day = 31
    elif (m == 2):
        if (((y % 4 == 0) and (y % 100 != 0)) or (y % 400) == 0):
            Last_Day = 29
        else:
            Last_Day = 28
    else:
        Last_Day = 30

    return Last_Day


def Local_Max(dt, wl, h, t):
    #This function checks wl for max value in +-t around time of h
    win_center = dt[h]
    win_start = win_center - t
    win_end = win_center + t
    max_val = wl[h]
    max_loc = h
    loc = h
    while dt[loc] >= win_start:
        if wl[loc] > max_val:
            max_val = wl[loc]
            max_loc = loc
        loc = loc - 1
    loc = h
    while dt[loc] <= win_end:
        if wl[loc] > max_val:
            max_val = wl[loc]
            max_loc = loc
        loc = loc + 1
    return dt[max_loc], wl[max_loc]


def Local_Min(dt, wl, h, t):
    #This function checks wl for min value in +-t around time of h
    win_center = dt[h]
    win_start = win_center - t
    win_end = win_center + t
    min_val = wl[h]
    min_loc = h
    loc = h
    while dt[loc] >= win_start:
        if wl[loc] < min_val:
            min_val = wl[loc]
            min_loc = loc
        loc = loc - 1
    loc = h
    while dt[loc] <= win_end:
        if wl[loc] < min_val:
            min_val = wl[loc]
            min_loc = loc
        loc = loc + 1
    return dt[min_loc], wl[min_loc]


def Calc_Expected_Diff(HL_Sub, HL_Con):
    #This function calculates the expected time difference between control and subrdinate tide occurences 

    Diff = 0
    start = 0
    Control_High_Dates = []
    Control_High_Indexes = []
    Control_Low_Dates = []
    Control_Low_Indexes = []
    for i in range(len(HL_Con)):
        if (HL_Con[i][2][0] == 'H'):
            Control_High_Dates.append(HL_Con[i][0]) 
            Control_High_Indexes.append(i)
        else:
            Control_Low_Dates.append(HL_Con[i][0])
            Control_Low_Indexes.append(i)

    Pairs = []

    for i in range(len(HL_Sub)):
        dt1 = HL_Sub[i][0]
        type1 = HL_Sub[i][2]

        if type1[0] == 'H':
            nt = Control_High_Indexes[Nearest_Tide(Control_High_Dates, dt1)]
        else:
            nt = Control_Low_Indexes[Nearest_Tide(Control_Low_Dates, dt1)]
        dt2 = HL_Con[nt][0]
        type2 = HL_Con[nt][2]
        tdiff = (dt2-dt1).days*1440 + (dt2-dt1).seconds/60
        if abs(tdiff) < 745:  #745 min = 12.42 hours)
            Pairs.append([dt1, dt2, type2, tdiff])
          
    #Calculate mean high and low differences
    MeanHDiff = 0.0
    NHighs = 0
    MeanLDiff = 0.0
    NLows = 0
    for i in range(len(Pairs)):
        if Pairs[i][2][0] == "H":
            MeanHDiff = MeanHDiff + Pairs[i][3]
            NHighs = NHighs + 1
        else:
            MeanLDiff = MeanLDiff + Pairs[i][3]
            NLows = NLows + 1
    MeanHDiff = MeanHDiff / NHighs
    MeanLDiff = MeanLDiff / NLows
    #Calculate mean of diffs above the mean and below the mean
    MeanHDiffAbove = 0.0
    MeanHDiffBelow = 0.0
    MeanLDiffAbove = 0.0
    MeanLDiffBelow = 0.0
    NHighsAbove = 0
    NHighsBelow = 0
    NLowsAbove = 0
    NLowsBelow = 0
    for i in range(len(Pairs)):
        if Pairs[i][2][0] == "H":
            if abs(Pairs[i][3]) > abs(MeanHDiff):
                MeanHDiffAbove = MeanHDiffAbove + Pairs[i][3]
                NHighsAbove = NHighsAbove + 1
            if abs(Pairs[i][3]) < abs(MeanHDiff):
                MeanHDiffBelow = MeanHDiffBelow + Pairs[i][3]
                NHighsBelow = NHighsBelow + 1
        else:
            if abs(Pairs[i][3]) > abs(MeanLDiff):
                MeanLDiffAbove = MeanLDiffAbove + Pairs[i][3]
                NLowsAbove = NLowsAbove + 1
            if abs(Pairs[i][3]) < abs(MeanLDiff):
                MeanLDiffBelow = MeanLDiffBelow + Pairs[i][3]
                NLowsBelow = NLowsBelow + 1
    if NHighsAbove==0:
        SDC_Print(['Error. No Highs above mean.'])
    else:
        MeanHDiffAbove = MeanHDiffAbove / NHighsAbove

    if NLowsAbove ==0:
        logger.warning('Error. No Lows above mean.')
    else:
        MeanLDiffAbove = MeanLDiffAbove / NLowsAbove

    if NHighsBelow == 0:
        logger.warning('Error. No Highs below mean.')
    else:
        MeanHDiffBelow = MeanHDiffBelow / NHighsBelow

    if NLowsBelow == 0:
        logger.warning('Error. No Lows below mean.')
    else:
        MeanLDiffBelow = MeanLDiffBelow / NLowsBelow

    if NHighsAbove == 0 or NLowsAbove == 0 or NHighsBelow == 0 or NLowsBelow == 0:
        logger.warning('***Error*** Fatal issue. Exiting Analysis.')
        exit(-1)

    if NHighsAbove > NHighsBelow:
        Diff = MeanHDiffAbove
    else:
        Diff = MeanHDiffBelow
    if NLowsAbove > NLowsBelow:
        Diff = (Diff + MeanLDiffAbove) / 2
    else:
        Diff = (Diff + MeanLDiffBelow) / 2

    return int(Diff)


def Local_Max_Fit(dt, wl, h):
    #This function fits the points around h to 3rd degree polynomial
    #and return the max of that function evaluated at the data's interval
    #Use a variable window width from 2.5 hrs to 6.5 hrs based on the range of the data

    #start with a 6.5 hr window and find range of included data
    win_center = dt[h]
    win_start = win_center - timedelta(hours=3) - timedelta(minutes=15)
    win_end = win_center + timedelta(hours=3) + timedelta(minutes=15)
    max_val = wl[h]
    min_val = wl[h]
    max_loc = h
    min_loc = h
    loc = h
    while dt[loc] >= win_start:
        if wl[loc] > max_val:
            max_val = wl[loc]
            max_loc = loc
        if wl[loc] < min_val:
            min_val = wl[loc]
            min_loc = loc
        loc = loc - 1
        if loc == 0:
            break
    loc = h
    while dt[loc] <= win_end:
        if wl[loc] > max_val:
            max_val = wl[loc]
            max_loc = loc
        if wl[loc] < min_val:
            min_val = wl[loc]
            min_loc = loc
        loc = loc + 1
        if loc == len(dt) - 1:
            break

    WL_Range = max_val-min_val
    #Bound the range between 0.5 and 1.0
    if WL_Range > 1.0:
        WL_Range = 1.0
    if WL_Range < 0.5:
        WL_Range = 1.0
    #/* Scale the window width between 25 and 65 slots based on the data range (only for small interval data) */
    if ((dt[h+1] - dt[h]) < timedelta(minutes=15)):

        win_width = timedelta(hours=6) + timedelta(minutes=30) - timedelta(minutes= ((WL_Range - 0.5) * 480.0))
        win_start = win_center - win_width/2
        win_end = win_center + win_width/2

    #collect the known points withn the window
    loc = h
    #find first value in window
    while dt[loc] >= win_start:
        first_in_window = loc
        loc = loc - 1
        if loc == 0:
            break

    #find the last value
    loc = h
    while dt[loc] <= win_end:
        last_in_window = loc 
        loc = loc + 1
        if loc == len(dt) - 1:
            break

    known_dates = []
    known_xs = []
    known_ys = []

    for l in range(first_in_window, last_in_window+1):
        known_dates.append(dt[l])
        known_xs.append((dt[l]-win_start).days*1440.0  + (dt[l]-win_start).seconds/60.0)
        known_ys.append(wl[l])
    
    #calculate polynomial
    z = np.polyfit(known_xs, known_ys, 3)
    f = np.poly1d(z)

    #calculate new x's and y's
    new_ys = f(known_xs)
    max = -999999.99
    for l in range(len(new_ys)):
        if new_ys[l] > max:
            max = new_ys[l]
            max_loc = first_in_window + l

    return dt[max_loc], max


def Local_Min_Fit(dt, wl, l):
    #This function fits the points around h to 3rd degree polynomial
    #and return the minimum of that function evaluated at the data's interval
    #Use a variable window width from 2.5 hrs to 6.5 hrs based on the range of the data

    #start with a 6.5 hr window and find range of included data
    win_center = dt[l]
    win_start = win_center - timedelta(hours=3) - timedelta(minutes=15)
    win_end = win_center + timedelta(hours=3) + timedelta(minutes=15)
    max_val = wl[l]
    min_val = wl[l]
    max_loc = l
    min_loc = l
    loc = l
    while dt[loc] >= win_start:
        if wl[loc] > max_val:
            max_val = wl[loc]
            max_loc = loc
        if wl[loc] < min_val:
            min_val = wl[loc]
            min_loc = loc
        loc = loc - 1
        if loc == 0:
            break
    loc = l
    while dt[loc] <= win_end:
        if wl[loc] > max_val:
            max_val = wl[loc]
            max_loc = loc
        if wl[loc] < min_val:
            min_val = wl[loc]
            min_loc = loc
        loc = loc + 1
        if loc == (len(dt) - 1):
            break

    WL_Range = max_val-min_val
    #Bound the range between 0.5 and 1.0
    if WL_Range > 1.0:
        WL_Range = 1.0
    if WL_Range < 0.5:
        WL_Range = 1.0
    #/* Scale the window width between 25 and 65 slots based on the data range (only for small interval data) */
    if ((dt[l+1] - dt[l]) < timedelta(minutes=15)):
        win_width = timedelta(hours=6) + timedelta(minutes=30) - timedelta(minutes= ((WL_Range - 0.5) * 480.0))
        win_start = win_center - win_width/2
        win_end = win_center + win_width/2

    #collect the known points withn the window
    loc = l
    #find first value in winow
    while dt[loc] >= win_start:
        first_in_window = loc
        loc = loc - 1
        if loc == 0:
            break

    #find the last value
    loc = l
    while dt[loc] <= win_end:
        last_in_window = loc 
        loc = loc + 1
        if loc == len(dt) - 1:
            break

    known_dates = []
    known_xs = []
    known_ys = []

    for ll in range(first_in_window, last_in_window+1):
        known_dates.append(dt[ll])
        known_xs.append((dt[ll]-win_start).days*1440.0  + (dt[ll]-win_start).seconds/60.0)
        known_ys.append(wl[ll])

    #calculate polynomial
    z = np.polyfit(known_xs, known_ys, 3)
    f = np.poly1d(z)

    #calculate new x's and y's
    new_ys = f(known_xs)
    min = 999999.99
    for ll in range(len(new_ys)):
        if new_ys[ll] < min:
            min = new_ys[ll]
            min_loc = first_in_window + ll

    return dt[min_loc], min

