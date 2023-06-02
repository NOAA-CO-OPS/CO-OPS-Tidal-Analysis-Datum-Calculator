"""CO-OPS API data retrieval functions for CO-OPS Datums Calculator"""

import urllib.request
import tides as tf
from datetime import datetime, date, time, timedelta

def Get_Monthly_Means(Control_Station_ID, Begin_Month, Begin_Year, End_Month, End_Year, Conversion):
    #This function retrieves the control station's monthly means using CO-OPS data api 
    
    MM = []
    #Creat empty list to append Monthly Means
    Months_Loaded = 0
    Months_Needed = (End_Year*12 + End_Month) - (Begin_Year*12 + Begin_Month) + 1
    #api for monthly means is limited to 10 years (120 months), so break the query into chunks if neccessary
    
    QBegYear = Begin_Year
    QBegMonth = Begin_Month
    if (Months_Needed) > 120:
        QEndYear = QBegYear
        print(QEndYear)
        QEndMonth = QBegMonth
        print( QEndMonth)
        for nn in range(118):
            if QEndMonth == 12:
                QEndYear = QEndYear + 1
                QEndMonth = 1
            else:
                QEndMonth = QEndMonth + 1
    else:
        QEndYear = End_Year
        QEndMonth = End_Month

    theYear = QBegYear
    print(theYear)
    theMonth = QBegMonth
    print(theMonth)
    while (Months_Loaded < Months_Needed):
        print ('From', QBegMonth, '/', QBegYear, '  to  ', QEndMonth, '/', QEndYear)
        Begin_Date = str(QBegYear) + str(QBegMonth).zfill(2) + '01%2000:00'
        End_Date = str(QEndYear) + str(QEndMonth).zfill(2) + str(tf.Last_Day_In_Month(QEndYear, QEndMonth)).zfill(2) + '%2023:54'
        url = 'https://tidesandcurrents.noaa.gov/api/datagetter?begin_date=' + Begin_Date + '&end_date=' + End_Date
        url = url + '&station=' + Control_Station_ID
        url = url + '&product=monthly_mean&datum=stnd&units=metric&time_zone=gmt&application=SDC&format=csv '
        req = urllib.request.Request(url)
        uh = urllib.request.urlopen(req)
        #Year, Month,  Highest, MHHW, MHW, MSL, MTL, MLW, MLLW, DTL, 
        #GT, MN, DHQ, DLQ, HWI, LWI, Lowest, Inferred 
        i = 0
        for line in uh:
            line = line.decode('utf-8')
            if i > 0:  #skip the header
                if not ('Error' in line): 
                    AM = line.split(',')
                    while not(int(AM[0]) == theYear and int(AM[1]) == theMonth):
                        MM.append([-99999.99, -99999.99, -99999.99, -99999.99, -99999.99, -99999.99])
                        if theMonth < 12:
                            theMonth = theMonth+1
                        else:
                            theYear = theYear+1
                            theMonth = 1  
                        Months_Loaded = Months_Loaded + 1  
                    #HWL,MHHW,MHW, MSL, MLW, MLLW, LWL
                    MM.append([float(AM[2])*Conversion, float(AM[3])*Conversion, float(AM[4])*Conversion, 
                               float(AM[5])*Conversion, float(AM[7])*Conversion, float(AM[8])*Conversion, 
                               float(AM[16])*Conversion])
                    if theMonth < 12:
                        theMonth = theMonth+1
                    else:
                        theYear = theYear+1
                        theMonth = 1    
                    Months_Loaded = Months_Loaded + 1
                else:
                    while (Months_Loaded < Months_Needed):
                        MM.append([-99999.99, -99999.99, -99999.99, -99999.99, -99999.99, -99999.99])
                        Months_Loaded = Months_Loaded + 1
            i = i+1
        if QEndMonth == 12:
            QBegMonth = 1
            QBegYear = QEndYear + 1
        else:
            QBegYear = QEndYear
            QBegMonth = QEndMonth + 1

        if ((End_Year*12 + End_Month) - (QBegYear*12 + QBegMonth) + 1) > 120:
            QEndYear = QBegYear
            QEndMonth = QBegMonth
            for nn in range(118):
                if QEndMonth == 12:
                    QEndYear = QEndYear + 1
                    QEndMonth = 1
                else:
                    QEndMonth = QEndMonth + 1
        else:
            QEndYear = End_Year
            QEndMonth = End_Month

    return MM

####################################################################

def Get_High_Lows(Control_Station_ID, Start_DT, End_DT, gmt_offset, Conversion):
    #This function retrieves control station high and low tides using CO-OPS data api

    print ('shifting gmt control data by', gmt_offset, 'hours')
    #if subordinate (short-term) station time is not in gmt, get time offset

    Start_DT = Start_DT + timedelta(hours=gmt_offset)
    End_DT = End_DT + timedelta(hours=gmt_offset)

    HL = []
    #Creat empty list to append Highs and Lows
    QBDT = Start_DT
    if (Start_DT + timedelta(weeks=52) - timedelta(minutes=1)) <= End_DT:
        #api is limited to 1 year at a time, so break up query into parts smaller thay 1 yr
        QEDT = Start_DT + timedelta(weeks=52) - timedelta(minutes=1)
    else:
        QEDT = End_DT
    while QEDT <= End_DT:
        Begin_Date = (str(QBDT.year) + str(QBDT.month).zfill(2) + str(QBDT.day).zfill(2) + '%20' + 
                      str(QBDT.hour).zfill(2) + ':' + str(QBDT.minute).zfill(2))
        End_Date = (str(QEDT.year) + str(QEDT.month).zfill(2) + str(QEDT.day).zfill(2) + '%20' +
                    str(QEDT.hour).zfill(2) + ':' + str(QEDT.minute).zfill(2))
        url = 'https://tidesandcurrents.noaa.gov/api/datagetter?begin_date=' + Begin_Date + '&end_date=' + End_Date
        url = url + '&station=' + Control_Station_ID
        url = url + '&product=high_low&datum=stnd&units=metric&time_zone=gmt&application=SDC&format=csv '
        req = urllib.request.Request(url)
        uh = urllib.request.urlopen(req)

        i = 0
        for line in uh:
            line = line.decode('utf-8')
            if i > 0:
                if not ('Error' in line) and len(line) > 10: 
                    AM = line.split(',')
                    HL.append([datetime.strptime(AM[0], "%Y-%m-%d %H:%M") - timedelta(hours=gmt_offset), 
                               float(AM[1])*Conversion, AM[2].strip()])
            i = i+1
        if QEDT == End_DT:
            QEDT = End_DT + timedelta(minutes=1)
        else:
            QBDT = QEDT + timedelta(minutes=1)
            if (QBDT + timedelta(weeks=52) - timedelta(minutes=1)) < End_DT:
                QEDT = QBDT + timedelta(weeks=52) - timedelta(minutes=1)
            else:
                QEDT = End_DT
    return HL

####################################################################

def Get_Accepted_Datums(Station_ID, Conversion):
    #This function retrieves the accepted control station datums using CO-OPS data api
    MHHW = None
    MHW = None
    DTL = None
    MTL = None
    MSL = None
    MLW = None
    MLLW = None
    GT = None
    MN = None
    DHQ = None
    DLQ = None
    NAVD = None
    LWI = None
    HWI = None
    HWL = None
    url = 'https://tidesandcurrents.noaa.gov/api/datagetter?&station=' + Station_ID
    url = url + '&product=datums&datum=stnd&units=metric&application=SDC&format=csv '
    req = urllib.request.Request(url)
    uh = urllib.request.urlopen(req)
    #Datum Name, Value (e.g. HWL, 6.444)
    i = 0
    for line in uh:
        line = line.decode('utf-8')
        if i > 0: 
            D = line.split(',')
            if D[0] == 'HWL':
                HWL = float(D[1])*Conversion
            elif D[0] == 'MHHW':
                MHHW = float(D[1])*Conversion
            elif D[0] == 'MHW':
                MHW = float(D[1])*Conversion
            elif D[0] == 'DTL':
                DTL = float(D[1])*Conversion
            elif D[0] == 'MTL':
                MTL = float(D[1])*Conversion
            elif D[0] == 'MSL':
                MSL = float(D[1])*Conversion
            elif D[0] == 'MLW':
                MLW = float(D[1])*Conversion
            elif D[0] == 'MLLW':
                MLLW = float(D[1])*Conversion
            elif D[0] == 'GT':
                GT = float(D[1])*Conversion
            elif D[0] == 'MN':
                MN = float(D[1])*Conversion
            elif D[0] == 'DHQ':
                DHQ = float(D[1])*Conversion
            elif D[0] == 'DLQ':
                DLQ = float(D[1])*Conversion
            elif D[0] == 'NAVD':
                NAVD = float(D[1])*Conversion
            elif D[0] == 'LWI':
                LWI = float(D[1])
            elif D[0] == 'HWI':
                HWI = float(D[1])
            else:
                print ('Unknown Datum Type: ', D[0])
        i = i+1
    SD = [MHHW, MHW, DTL, MTL, MSL, MLW, MLLW, GT, MN, DHQ, DLQ, NAVD, LWI, HWI] 
    return SD

####################################################################

def Get_SubMethod(Station_ID):
    #This function checks if the control station is a West coast/Pacific or  East Coast/Gulf Coast/Caribbean Island station
    #for choosing datum computation method
    import json
    url = 'https://api.tidesandcurrents.noaa.gov/mdapi/prod/webapi/stations/' + Station_ID + '.json'
    req = urllib.request.Request(url)
    uh = urllib.request.urlopen(req)
    d = json.load(uh)
    try:
        lat = d['stations'][0]['lat']
        lon = d['stations'][0]['lng']
    except:
        return("No location info for station ID " + Station_ID)
    if lon < -100:
        return('Standard')
    else:
        return('Modified')


