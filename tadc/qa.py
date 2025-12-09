import datetime
import logging
import numpy as np
import pandas as pd
import requests
logger = logging.getLogger(__name__)


class Assurances:
    def __init__(self, ts, resample_minutes):
        self.ts = ts
        self.resample_minutes = resample_minutes
        
    def assure_no_unreadable_values(self):
        ts = self.ts.rename(columns={self.ts.columns[0]:'time',self.ts.columns[1]:'val'})
        try:
            ts['val'].astype(float)
        except Exception as error:
            bad_val = error.args[0].split(':')[1].replace(' ','').replace("'","")
            ts.loc[ts['val']==bad_val,'val'] = 'NaN'
            ts['val'] = ts['val'].astype(float)
            logger.warning('WARNING: Unreadable value found in data: ' + bad_val + '. Replacing all occurrences with NaNs.')            
        self.ts = ts
        
    def assure_even_temporal_spacing(self):
        def resamp(ts, interval_want, interval_mean):
            if interval_want < interval_mean:
                logger.warning('WARNING: Input resampling rate is higher than data sampling rate. This may result in unstable behavior. Consider resampling to a lower rate.')
            ti = pd.date_range(ts['time'].iloc[0],ts['time'].iloc[-1],freq=interval_want)
            ts_interp = ts.set_index('time').reindex(ti,method='nearest',tolerance=interval_mean).reset_index().rename(columns={'index':'time'})  # Reinterpolate being careful to preserve gaps #
            return ts_interp         
        ts = self.ts
        ts['time'] = pd.to_datetime(ts['time'])         
        ts = ts.groupby('time').first().reset_index()  # Remove duplicates #
        time_diffs_all = ts['time'].diff()
        interval_mean = ts['time'].diff().mean().floor('min')
        if len(time_diffs_all.unique().dropna()) > 1:
            if self.resample_minutes is None:
                interval_want = interval_mean
                logger.warning('WARNING: Input timeseries has uneven temporal spacing. Re-interpolating to a spacing of ' + str(round(interval_want.seconds/60, 2)) + ' minutes per sample.')

            else:
                interval_want = pd.Timedelta(minutes = self.resample_minutes)
                logger.warning('WARNING: Input timeseries has uneven temporal spacing. Re-interpolating to a spacing of ' + str(round(interval_want.seconds/60, 2)) + ' minutes per sample.')
            self.ts = resamp(ts, interval_want, interval_mean)
        else:
            if self.resample_minutes is not None:
                interval_want = pd.Timedelta(minutes = self.resample_minutes)
                logger.warning('WARNING: Re-interpolating to a spacing of ' + str(round(interval_want.seconds/60, 2)) + ' minutes per sample.')
                self.ts = resamp(ts, interval_want, interval_mean)
            else:
                self.ts = ts
                
    def assure_flatlines_are_gaps(self):
        is_flatline = self.ts['val'].diff().abs() < 0.001
        consecutive_groups = is_flatline.ne(is_flatline.shift()).cumsum()
        group_sizes = consecutive_groups.groupby(consecutive_groups).transform('size')
        is_long_flatline = is_flatline & (group_sizes > 1)
        is_end_of_run = is_long_flatline & is_long_flatline.ne(is_long_flatline.shift(-1))
        is_long_flatline_final = is_long_flatline & ~is_end_of_run
        if is_long_flatline_final.sum() > 0:
            self.ts.loc[is_long_flatline_final,'val'] = np.nan
            logger.warning('WARNING: Flatlines detected. Treating flatlines as missing data.')
            
        
def run(ts, resample_minutes):
    assurances = Assurances(ts, resample_minutes)
    assurances.assure_no_unreadable_values()
    assurances.assure_even_temporal_spacing()
    assurances.assure_flatlines_are_gaps()
    return assurances.ts
