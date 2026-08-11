"""CO-OPS API data retrieval functions for CO-OPS Datums Calculator"""

from datetime import datetime, date, time, timedelta
import logging
import numpy as np
import pandas as pd
import requests
logger = logging.getLogger(__name__)


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
        MM[c]  = MM[c].replace('','nan').astype(float) * Conversion
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
        if 'data' in r.json():
            hl_chunks.append(pd.DataFrame(r.json()['data']))
        else:
            fill_t = pd.date_range(start_dt,end_dt,freq='6h')
            fill = pd.DataFrame({'t':fill_t,'v':np.nan,'ty':'unknown','f':'unknown'})
            hl_chunks.append(fill)
    HL = pd.concat(hl_chunks,ignore_index=True)
    HL['t'] = pd.to_datetime(HL['t']) - timedelta(hours=gmt_offset)
    HL['v']  = HL['v'].astype(float) * Conversion
    HL['ty'] = [HL['ty'].iloc[i].replace(' ','') for i in range(len(HL))]
    HL_lists = [HL[['t','v','ty']].iloc[i].values.tolist() for i in range(len(HL))]  # Convert to the list of lists format needed by run.py #
    return HL_lists


def Get_Accepted_Datums(Station_ID, epoch_start_year, gmt_offset, Conversion):
    #This function retrieves the accepted control station datums using CO-OPS metadata api
    if epoch_start_year == 1983:
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
                if datum not in ['LWI','HWI']:
                    SD.append(val * Conversion)
                else:
                    SD.append(val)
    else:
        if epoch_start_year == 2002:
            logger.warning(('WARNING: You have requested to compute datums for your data relative to the 2002-2020 National Tidal Datum Epoch, ' +
                            'which has not yet been released by NOAA. Preliminary, unofficial datums at the selected control station have been computed ' +
                            'and used to adjust your data. These control datums may be different from those that are planned to be officially released ' +
                            'in early 2029 and should be used for planning purposes only. The chosen control station should have similar tidal characteristics ' +
                            'to the subordinate station. More details on choosing a suitable control station can be found here: ' +
                            'https://access.co-ops.nos.noaa.gov/datumcalc/docs/FAQs.pdf.'))
        else:
            logger.warning(('WARNING: You have requested to compute datums for your data relative to a custom 19 yr Datum Epoch. Preliminary, unofficial datums ' +
                            'at the selected control station will be computed and used to adjust your data. These control datums should be used for planning purposes only. ' +
                            'The chosen control station should have similar tidal characteristics to the subordinate station. More details on choosing a suitable control ' +
                            'station can be found here: https://access.co-ops.nos.noaa.gov/datumcalc/docs/FAQs.pdf.'))            
        mm = Get_Monthly_Means(Station_ID,
                               1,
                               epoch_start_year,
                               12,
                               epoch_start_year + 18,
                               Conversion)
        MM = pd.DataFrame(mm,columns=['highest','MHHW','MHW','MSL','MLW','MLLW','lowest'])
        if len(MM) >= 120: # If at least 10 years of data, do the calculation
            if len(MM) < 192: # But if less than 16 years of data, throw a warning #
                logger.warning(('WARNING: Control station is missing more than 3 yr of data for the selected 19 year epoch. ' +
                                'Control datums may be unreliable. Consider choosing a different control station.'))
            MHHW = MM['MHHW'].mean()
            MHW = MM['MHW'].mean()
            MSL = MM['MSL'].mean()
            MLW = MM['MLW'].mean()
            MLLW = MM['MLLW'].mean()
            MTL = 0.5 * (MHW + MLW)
            DTL = 0.5 * (MHHW + MLLW)
            GT = MHHW - MLLW
            MN = MHW - MLW
            DHQ = MHHW - MHW
            DLQ = MLW - MLLW
            NAVD88 = np.nan
            LWI = np.nan
            HWI = np.nan
            SD = [MHHW*Conversion,MHW*Conversion,DTL*Conversion,MTL*Conversion,MSL*Conversion,
                  MLW*Conversion,MLLW*Conversion,GT*Conversion,MN*Conversion,DHQ*Conversion,
                  DLQ*Conversion,NAVD88*Conversion,LWI,HWI]              
        else:
            raise RuntimeError('Control station is missing more than 9 yr of data for the selected 19 year epoch. Please select a different control station.')
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
