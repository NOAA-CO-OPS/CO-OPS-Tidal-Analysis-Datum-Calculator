import datetime
import numpy as np
import pandas as pd
import requests
import warnings


class Assurances:
    def __init__(self, ts):
        self.ts = ts
        
    def assure_even_temporal_spacing(self):
        ts = self.ts.rename(columns={self.ts.columns[0]:'time',self.ts.columns[1]:'val'})
        ts['time'] = pd.to_datetime(ts['time'])
        ts = ts.groupby('time').first().reset_index()  # Remove duplicates #
        time_diffs_all = ts['time'].diff()
        if len(time_diffs_all.unique().dropna()) > 1:
            interval_want = ts['time'].diff().mean().floor('min')
            ti = pd.date_range(ts['time'].iloc[0],ts['time'].iloc[-1],freq=interval_want)
            ts_interp = ts.set_index('time').reindex(ti,method='nearest',tolerance=interval_want).reset_index().rename(columns={'index':'time'})  # Reinterpolate being careful to preserve gaps #
            self.ts = ts_interp
            warnings.warn('Input timeseries has uneven temporal spacing. Re-interpolating to a spacing of ' + str(round(interval_want.seconds/60, 2)) + ' minutes per sample.')

    def assure_flatlines_are_gaps(self):
        is_flatline = self.ts['val'].diff().abs() < 0.001
        consecutive_groups = is_flatline.ne(is_flatline.shift()).cumsum()
        group_sizes = consecutive_groups.groupby(consecutive_groups).transform('size')
        is_long_flatline = is_flatline & (group_sizes > 1)
        is_end_of_run = is_long_flatline & is_long_flatline.ne(is_long_flatline.shift(-1))
        is_long_flatline_final = is_long_flatline & ~is_end_of_run
        if is_long_flatline_final.sum() > 0:
            self.ts.loc[is_long_flatline_final,'val'] = np.nan
            warnings.warn('Flatlines detected. Treating flatlines as missing data.')
            
        
def run(ts):
    assurances = Assurances(ts)
    assurances.assure_even_temporal_spacing()
    assurances.assure_flatlines_are_gaps()
    return assurances.ts
