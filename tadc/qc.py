import datetime
import numpy as np
import pandas as pd
import requests
from scipy.signal import periodogram
import warnings


class Tests:
    def __init__(self, ts, control_station_id, subordinate_lat, subordinate_lon):
        self.ts = ts
        self.control_station_id = control_station_id
        self.subordinate_lat = subordinate_lat
        self.subordinate_lon = subordinate_lon
        
    def check_csv_format(self):
        if len(self.ts.columns) != 2:
            raise RuntimeError("Input csv file must contain two columns: time and water level")
        else:
            self.ts = self.ts.rename(columns={self.ts.columns[0]:'time',self.ts.columns[1]:'val'})

    def check_date_format(self):
        try:
            self.ts['time'] = pd.to_datetime(self.ts['time'])
        except ValueError:
            raise ValueError("Timestamps could not be interpreted.")

    def check_control_station_distance(self):
        if self.control_station_id != None:
            r = requests.get('https://api.tidesandcurrents.noaa.gov/mdapi/prod/webapi/stations/' + str(self.control_station_id) + '.json?units=english')
            lat_control = r.json()['stations'][0]['lat']
            lon_control = r.json()['stations'][0]['lng']
            d = self._haversine(self.subordinate_lat,self.subordinate_lon,lat_control,lon_control)
            if d > 10:
                warnings.warn('Control station is ' + str(round(d,2)) + ' km from subordinate station.')

    def check_rough_tidal_amplitude(self):
        fs = 1/self.ts['time'].diff()[1].seconds  # Sampling frequency #
        y_demean = self.ts['val'] - self.ts['val'].mean()  # De-mean the timeseries #
        y_demean[y_demean.isna()] = 0
        f, Pxx = periodogram(y_demean,fs)  # Get the PSD curve. Units are m^2/HZ #
        cpd = f*60*60*24  # Convert from HZ to cycles per day #
        var = np.trapz(Pxx, x=cpd)/60/60/24  # Total variance is the integral under the PSD curve assuming the input data is demeaned #
        var48 = np.trapz(Pxx[cpd<=2], x=cpd[cpd<=2])/60/60/24  # Variance at periods less than or equal to 2 days #
        if var48 < 0.02 or var48/var < 0.75:  # These are thresholds that seem to be reasonable - may need tweaking.
            raise RuntimeError('Input data do not appear to have a strong enough tidal signal for reliable analysis.')

    @staticmethod
    def _haversine(lat1, lon1, lat2, lon2):
        lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
        R = 6371 
        dlat = lat2 - lat1
        dlon = lon2 - lon1
        a = np.sin(dlat/2.0)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2.0)**2
        c = 2 * np.arcsin(np.sqrt(a))
        km = R * c
        return km   
                

def run(ts, control_station_id, subordinate_lat, subordinate_lon):
    tests = Tests(ts, control_station_id, subordinate_lat, subordinate_lon)
    tests.check_csv_format()
    tests.check_date_format()
    tests.check_control_station_distance()
    tests.check_rough_tidal_amplitude()
