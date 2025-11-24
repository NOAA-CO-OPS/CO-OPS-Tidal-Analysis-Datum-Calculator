import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


class Out:
    def __init__(self, daily_maxs):
        self.daily_maxs = daily_maxs

    def percentile(self, prctile):
        return self.daily_maxs['elevation'].quantile(prctile/100)

    def plot(self, prctile=None):
        fig,ax = plt.subplots(1,figsize=(9,5))
        ax.tick_params(axis='both',labelsize=8)
        ax.grid('on',linestyle='--')
        ax.plot(self.daily_maxs['time'],self.daily_maxs['elevation'],'-o',label='Daily max',zorder=2)
        ax.set_ylabel('Elevation (m)',fontsize=8)
        if prctile != None:
            prctile_elev = self.percentile(prctile)
            ax.set_xlim(ax.get_xlim())
            ax.plot(ax.get_xlim(),[prctile_elev,prctile_elev],'k--',label=str(prctile)+' percentile',zorder=3)
            ax.legend(fontsize=8)
        fig.show()
        return fig

        
def run(datum, data, datums):
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
    dmi = data_dwant.groupby(data_dwant.index.date)['val'].idxmax()
    dm = data_dwant.loc[dmi].reset_index()
    dm = dm.rename(columns={'time':'time','val':'elevation'})
    dm['completeness'] = per_complete.values

    return Out(dm)
