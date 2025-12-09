import datetime
import logging
import numpy as np
import pandas as pd
import requests
from scipy.signal import periodogram
logger = logging.getLogger(__name__)


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
                logger.warn('WARNING: Control station is ' + str(round(d,2)) + ' km from subordinate station.')

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
