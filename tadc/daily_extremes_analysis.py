import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


class Out:
    def __init__(self, daily_extremes, extremes_type, datum, units):
        self.__daily_extremes = daily_extremes
        if extremes_type == 'max':
            self.daily_maxs = daily_extremes
            self.daily_mins = None
        elif extremes_type == 'min':
            self.daily_mins = daily_extremes
            self.daily_maxs = None
        self.datum = datum
        self.units = units
            
    def percentile(self, prctile):
        return self.__daily_extremes['elevation'].quantile(prctile/100)

    def plot(self, prctile=None):
        fig,ax = plt.subplots(1,figsize=(9,5))
        ax.tick_params(axis='both',labelsize=8)
        ax.grid('on',linestyle='--')
        ax.plot(self.__daily_extremes['time'],self.__daily_extremes['elevation'],'-o',label='Daily max',zorder=2)
        ax.set_ylabel('Elevation ('+self.units+', '+self.datum+')',fontsize=8)
        if prctile != None:
            prctile_elev = self.percentile(prctile)
            ax.set_xlim(ax.get_xlim())
            ax.plot(ax.get_xlim(),[prctile_elev,prctile_elev],'k--',label=str(prctile)+' percentile',zorder=3)
            ax.legend(fontsize=8)
        if self.daily_maxs is not None:
            ax.set_title('Daily Maximum Water Levels',fontsize=8)
        else:
            ax.set_title('Daily Minimum Water Levels',fontsize=8)
        total_dt = self.__daily_extremes['time'].iloc[-1] - self.__daily_extremes['time'].iloc[0]
        ticks = pd.date_range(self.__daily_extremes['time'].iloc[0],
                              self.__daily_extremes['time'].iloc[-1],
                              freq=total_dt/8)
        ax.set_xlim(self.__daily_extremes['time'].iloc[0] - (total_dt/8/4),
                    self.__daily_extremes['time'].iloc[-1] + (total_dt/8/4))
        ax.set_xticks(ticks)
        fig.autofmt_xdate()
        fig.show()
        return fig

       
def run(extremes_type, datum, data, datums, units):
    # Get timestamps into a usable format #
    data = data.rename(columns={data.columns[0]:'time',data.columns[1]:'val'})
    data['time'] = pd.to_datetime(data['time'])
    data = data.replace(-99999.99, np.nan)

    # Put the data onto the threshold datum and onto MHHW #
    data_dwant = pd.DataFrame({'time':data['time'],'val':data['val']-datums[datum]})

    # Calc daily maxes #
    data_dwant = data_dwant.set_index('time')
    interval_hrs = (data_dwant.index[1] - data_dwant.index[0]).seconds/3600
    n = data_dwant.groupby(data_dwant.index.date)['val'].size()
    per_complete = n / (24 / interval_hrs) * 100
    if extremes_type == 'max':
        dmi = data_dwant.groupby(data_dwant.index.date)['val'].idxmax()
    elif extremes_type == 'min':
        dmi = data_dwant.groupby(data_dwant.index.date)['val'].idxmin()
    dm = data_dwant.loc[dmi].reset_index()
    dm = dm.rename(columns={'time':'time','val':'elevation'})
    dm['completeness'] = per_complete.values

    return Out(dm, extremes_type, datum, units)
