"""CO-OPS API data retrieval functions for CO-OPS Datums Calculator"""

from datetime import datetime, date, time, timedelta
import numpy as np
import pandas as pd
import requests

from . import tides as tf


def Get_Monthly_Means(Control_Station_ID, Begin_Month, Begin_Year, End_Month, End_Year, Conversion):
    #This function retrieves the control station's monthly means using CO-OPS data api 
    end_days = tf.Last_Day_In_Month(int(End_Year),int(End_Month))
    if int(Begin_Month) < 10:
        sb = '0'
    else:
        sb = ''
    if int(End_Month) < 10:
        se = '0'
    else:
        se = ''    
    url1 = 'https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?'
    url2 = 'begin_date=' + str(Begin_Year) + sb + str(Begin_Month) + '01' + '&end_date=' + str(End_Year) + se + str(End_Month) + str(end_days) + '&station=' + str(Control_Station_ID)
    url3 = '&product=monthly_mean&datum=stnd&units=metric&time_zone=gmt&application=TADC&format=json'
    r = requests.get(url1 + url2 + url3)
    MM = pd.DataFrame(r.json()['data'])
    for c in ['highest','MHHW','MHW','MSL','MLW','MLLW','lowest']:
        MM[c]  = MM[c].astype(float) * Conversion
    MM_lists = [MM[['highest','MHHW','MHW','MSL','MLW','MLLW','lowest']].iloc[i].values.tolist() for i in range(len(MM))]  # Convert to the list of lists format needed by run.py #
    return MM_lists


def Get_High_Lows(Control_Station_ID, Start_DT, End_DT, gmt_offset, Conversion):
    #This function retrieves control station high and low tides using CO-OPS data api

    #if subordinate (short-term) station time is not in gmt, get time offset
    Start_DT += timedelta(hours=gmt_offset)
    End_DT += timedelta(hours=gmt_offset)

    if End_DT - Start_DT > timedelta(days=365):
        chunks = pd.date_range(Start_DT, End_DT, periods=int(np.ceil((End_DT - Start_DT).days/365))+1)
    else:
        chunks = (Start_DT, End_DT)

    hl_chunks = []
    for i in range(len(chunks)-1):
        start_dt = chunks[i]
        end_dt = chunks[i+1]

        start_datestr = datetime.strftime(start_dt,'%Y%m%d')
        end_datestr = datetime.strftime(end_dt,'%Y%m%d')

        url1 = 'https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?'
        url2 = 'begin_date=' + start_datestr + '&end_date=' + end_datestr + '&station=' + str(Control_Station_ID)
        url3 = '&product=High_low&datum=stnd&units=metric&time_zone=gmt&application=TADC&format=json'    
        r = requests.get(url1 + url2 + url3)

        hl_chunks.append(pd.DataFrame(r.json()['data']))
    HL = pd.concat(hl_chunks,ignore_index=True)
    HL['t'] = pd.to_datetime(HL['t']) - timedelta(hours=gmt_offset)
    HL['v']  = HL['v'].astype(float) * Conversion
    HL['ty'] = [HL['ty'].iloc[i].replace(' ','') for i in range(len(HL))]
    HL_lists = [HL[['t','v','ty']].iloc[i].values.tolist() for i in range(len(HL))]  # Convert to the list of lists format needed by run.py #
    return HL_lists


def Get_Accepted_Datums(Station_ID, Conversion):
    #This function retrieves the accepted control station datums using CO-OPS metadata api
    url = 'https://api.tidesandcurrents.noaa.gov/mdapi/prod/webapi/stations/' + str(Station_ID) + '/datums.json?units=metric'
    r = requests.get(url)
    datums = pd.DataFrame(r.json()['datums'])
    SD = []
    for datum in ['MHHW','MHW','DTL','MTL','MSL','MLW','MLLW','GT','MN','DHQ','DLQ','NAVD88','LWI','HWI']:
        try:
            val = datums.loc[datums['name'] == datum,'value'].values[0]
        except IndexError:
            SD.append(np.nan)
        else:
            SD.append(val)
    return SD


def Get_SubMethod(Station_ID):
    #This function checks if the control station is a West coast/Pacific or  East Coast/Gulf Coast/Caribbean Island station
    #for choosing datum computation method
    url = 'https://api.tidesandcurrents.noaa.gov/mdapi/prod/webapi/stations/' + str(Station_ID) + '.json?units=metric'
    r = requests.get(url)
    lon = r.json()['stations'][0]['lng']
    if lon < -100:
        return('Standard')
    else:
        return('Modified')


