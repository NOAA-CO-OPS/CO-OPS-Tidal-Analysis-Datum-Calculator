from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import pandas as pd


class Out:
    def __init__(self, inundations, datum, units, fun_inps):
        self.inundations = inundations
        self.datum = datum
        self.units = units
        self.__fun_inps = fun_inps

    def plot(self):
        fig_picks,ax = plt.subplots(1,figsize=(9,5))
        ax.plot(self.__fun_inps['data']['time'],self.__fun_inps['data']['val'] - self.__fun_inps['datums'][self.__fun_inps['threshold_datum']],zorder=2,label='Data')
        ax.plot(ax.get_xlim(),[self.__fun_inps['threshold'],self.__fun_inps['threshold']],'k--',zorder=3,label='Threshold')
        ax.set_ylabel('Elevation ('+self.units+', '+self.datum+')',fontsize=8)
        ax.grid('on',linestyle='--')
        ax.tick_params(axis='both',labelsize=8)
        ax.plot(self.inundations['Period Start'],np.tile(self.__fun_inps['threshold'],len(self.inundations)),'ks',markerfacecolor='r',zorder=2,label='Threshold crossing')
        ax.plot(self.inundations['Period End'],np.tile(self.__fun_inps['threshold'],len(self.inundations)),'ks',markerfacecolor='r',zorder=2)
        ax.legend(fontsize=8)
        total_dt = self.__fun_inps['data']['time'].iloc[-1] - self.__fun_inps['data']['time'].iloc[0]
        ticks = pd.date_range(self.__fun_inps['data']['time'].iloc[0],
                              self.__fun_inps['data']['time'].iloc[-1],
                              freq=total_dt/8)
        ax.set_xlim(self.__fun_inps['data']['time'].iloc[0] - (total_dt/8/4),
                    self.__fun_inps['data']['time'].iloc[-1] + (total_dt/8/4))
        ax.set_xticks(ticks)
        fig_picks.autofmt_xdate()
        ax.set_title(('Inundation History\n' + 
              'Threshold = ' + str(self.__fun_inps['threshold']) + ' '+self.units+' above ' + self.__fun_inps['threshold_datum'] + '\n' +
              'Time range =  ' + datetime.strftime(self.__fun_inps['data']['time'].iloc[0],'%Y-%m-%d') + ' to ' +  datetime.strftime(self.__fun_inps['data']['time'].iloc[-1],'%Y-%m-%d') + '\n' +
              'Results: ' + str(len(self.inundations)) + ' Inundations. Total Duration = ' + str(round(self.inundations['Duration (hours)'].sum(),2)) + ' hours '+
              '(' + str(round(self.inundations['Duration (hours)'].sum()/((self.__fun_inps['data']['time'].iloc[-1]- self.__fun_inps['data']['time'].iloc[0]).total_seconds()/60/60)*100,2)) + '%)'),
              fontsize=8, fontweight='normal', loc='left', ha='left')
        fig_picks.show()

        
        fig_d_vs_h,axx = plt.subplots(1)
        axx.plot(self.inundations['Duration (hours)'],self.inundations['Maximum Elevation Above Threshold'],'ko',markerfacecolor='gray',zorder=2)
        axx.grid('on',linestyle='--')
        axx.set_xlabel('Duration of Inundation (Hours)',fontsize=8)
        axx.set_ylabel('Maximum Elevation ('+self.units+') Above Threshold',fontsize=8)
        axx.tick_params(axis='both',labelsize=8)
        axx.set_ylim(0,axx.get_ylim()[-1])
        axx.set_title(('Maximum Elevation vs. Duration of Inundation\n' + 
                      'Threshold = ' + str(self.__fun_inps['threshold']) + ' '+self.units+' above ' + self.__fun_inps['threshold_datum'] + '\n' +
                      'Time range =  ' + datetime.strftime(self.__fun_inps['data']['time'].iloc[0],'%Y-%m-%d') + ' to ' +  datetime.strftime(self.__fun_inps['data']['time'].iloc[-1],'%Y-%m-%d') + '\n' +
                      'Results: ' + str(len(self.inundations)) + ' Inundations. Total Duration = ' + str(round(self.inundations['Duration (hours)'].sum(),2)) + ' hours '+
                      '(' + str(round(self.inundations['Duration (hours)'].sum()/((self.__fun_inps['data']['time'].iloc[-1]- self.__fun_inps['data']['time'].iloc[0]).total_seconds()/60/60)*100,2)) + '%)'),
                      fontsize=8, fontweight='normal', loc='left', ha='left')
        fig_d_vs_h.show()

        if self.units == 'Feet' or self.units == 'feet' or self.units == 'ft':
            bs = 0.1
        elif self.units == 'Inches' or self.units == 'inches' or self.units == 'in':
            bs = 1
        elif self.units == 'Meters' or self.units == 'meters' or self.units == 'm':
            bs = 0.05
        elif self.units == 'Centimeters' or self.units == 'centimeters' or self.units == 'cm':
            bs = 5
        elif self.units == 'Millimeters' or self.units == 'millimeters' or self.units == 'mm':
            bs = 50
        bins = np.arange(0,self.inundations['Maximum Elevation Above Threshold'].max()+bs,bs)
        labels = [
            f"{start if i == 0 else start + 0.001:.3f} - {end:.3f}" 
            for i, (start, end) in enumerate(zip(bins[:-1], bins[1:]))
        ]
        bin_vals = pd.cut(self.inundations['Maximum Elevation Above Threshold'],bins,right=True,labels=labels)
        bin_counts = bin_vals.value_counts(sort=False)
        reverse_cumsum = bin_counts.sort_index(ascending=False).cumsum().sort_index()
        reverse_percentages = (reverse_cumsum / bin_counts.sum()) * 100
        fig_f_of_elev,axxx = plt.subplots(1,figsize=(9,5))
        axxx2 = axxx.twinx()
        fig_f_of_elev.subplots_adjust(bottom=0.2)
        axxx.tick_params(axis='x',rotation=25)
        axxx.tick_params(axis='both',labelsize=8)
        l1 = bin_counts.plot.bar(ax=axxx, width=0.9, edgecolor='black', color='royalblue',
                            rot=45, label='Frequency',zorder=2)
        axxx.set_xlabel('Inundation Exceedance ('+self.units+')',fontsize=8)
        axxx.set_ylabel('Frequency',fontsize=8)
        axxx.set_ylim(0,bin_counts.max()+0.5)
        l2 = reverse_percentages.plot.line(ax=axxx2, color='red', marker='o', linewidth=2,
                                      label='Percentage of Inundations',zorder=2)
        axxx2.set_ylim(0,105)
        axxx2.set_ylabel('Percentage of Inundation',fontsize=8)
        axxx2.tick_params(axis='both',labelsize=8)
        axxx.grid('on',linestyle='--')
        axxx.set_title(('Frequency of Elevations\n' + 
              'Threshold = ' + str(self.__fun_inps['threshold']) + ' '+self.units+' above ' + self.__fun_inps['threshold_datum'] + '\n' +
              'Time range =  ' + datetime.strftime(self.__fun_inps['data']['time'].iloc[0],'%Y-%m-%d') + ' to ' +  datetime.strftime(self.__fun_inps['data']['time'].iloc[-1],'%Y-%m-%d') + '\n' +
              'Results: ' + str(len(self.inundations)) + ' Inundations. Total Duration = ' + str(round(self.inundations['Duration (hours)'].sum(),2)) + ' hours '+
              '(' + str(round(self.inundations['Duration (hours)'].sum()/((self.__fun_inps['data']['time'].iloc[-1]- self.__fun_inps['data']['time'].iloc[0]).total_seconds()/60/60)*100,2)) + '%)'),
              fontsize=8, fontweight='normal', loc='left', ha='left')
        handles1, labels1 = axxx.get_legend_handles_labels()
        handles2, labels2 = axxx2.get_legend_handles_labels()
        axxx.legend(handles1 + handles2, labels1 + labels2, loc='upper right',fontsize=8)
        fig_f_of_elev.show()

        bs = 1
        bins = np.arange(0,self.inundations['Duration (hours)'].max()+bs,bs)
        labels = [
            f"{start if i == 0 else start + 0.001:.3f} - {end:.3f}" 
            for i, (start, end) in enumerate(zip(bins[:-1], bins[1:]))
        ]
        bin_vals = pd.cut(self.inundations['Duration (hours)'],bins,right=True,labels=labels)
        bin_counts = bin_vals.value_counts(sort=False)
        reverse_cumsum = bin_counts.sort_index(ascending=False).cumsum().sort_index()
        reverse_percentages = (reverse_cumsum / bin_counts.sum()) * 100
        fig_f_of_dur,axxx = plt.subplots(1,figsize=(9,5))
        axxx2 = axxx.twinx()
        fig_f_of_dur.subplots_adjust(bottom=0.2)
        axxx.tick_params(axis='x',rotation=25)
        axxx.tick_params(axis='both',labelsize=8)
        l1 = bin_counts.plot.bar(ax=axxx, width=0.9, edgecolor='black', color='royalblue',
                            rot=45, label='Frequency',zorder=2)
        axxx.set_xlabel('Inundation Exceedance (Hours)',fontsize=8)
        axxx.set_ylabel('Frequency',fontsize=8)
        axxx.set_ylim(0,bin_counts.max()+0.5)
        l2 = reverse_percentages.plot.line(ax=axxx2, color='red', marker='o', linewidth=2,
                                      label='Percentage of Inundations',zorder=2)
        axxx2.set_ylim(0,105)
        axxx2.set_ylabel('Percentage of Inundation',fontsize=8)
        axxx2.tick_params(axis='both',labelsize=8)
        axxx.grid('on',linestyle='--')
        axxx.set_title(('Frequency of Durations\n' + 
              'Threshold = ' + str(self.__fun_inps['threshold']) + ' '+self.units+' above ' + self.__fun_inps['threshold_datum'] + '\n' +
              'Time range =  ' + datetime.strftime(self.__fun_inps['data']['time'].iloc[0],'%Y-%m-%d') + ' to ' +  datetime.strftime(self.__fun_inps['data']['time'].iloc[-1],'%Y-%m-%d') + '\n' +
              'Results: ' + str(len(self.inundations)) + ' Inundations. Total Duration = ' + str(round(self.inundations['Duration (hours)'].sum(),2)) + ' hours '+
              '(' + str(round(self.inundations['Duration (hours)'].sum()/((self.__fun_inps['data']['time'].iloc[-1]- self.__fun_inps['data']['time'].iloc[0]).total_seconds()/60/60)*100,2)) + '%)'),
              fontsize=8, fontweight='normal', loc='left', ha='left')
        handles1, labels1 = axxx.get_legend_handles_labels()
        handles2, labels2 = axxx2.get_legend_handles_labels()
        axxx.legend(handles1 + handles2, labels1 + labels2, loc='upper right',fontsize=8)
        fig_f_of_dur.show()        
        return [fig_picks, fig_d_vs_h, fig_f_of_elev, fig_f_of_dur]


def run(threshold, threshold_datum, data, datums, high_lows, units):
    # Get timestamps into a usable format #
    data = data.rename(columns={data.columns[0]:'time',data.columns[1]:'val'})
    data['time'] = pd.to_datetime(data['time'])
    data = data.replace(-99999.99, np.nan)

    # Put the data onto the threshold datum and onto MHHW #
    data_mhhw = pd.DataFrame({'time':data['time'],'val':data['val']-datums['MHHW']})
    data_dwant = pd.DataFrame({'time':data['time'],'val':data['val']-datums[threshold_datum]})

    # Separate threshold exceedances into temporally separate groups #
    up_crosses = np.where((data_dwant['val']>threshold) & (data_dwant['val'].shift(1)<=threshold))[0]
    down_crosses = np.where((data_dwant['val']<=threshold) & (data_dwant['val'].shift(1)>threshold))[0]
    if len(up_crosses) == 0 and len(down_crosses) == 0:
        return pd.DataFrame()
    else:
        exceedance_groups = []
        i_up = 0
        i_down = 0
        while i_up < min(len(up_crosses)-1,len(down_crosses)) and i_down < min(len(up_crosses)-1,len(down_crosses)):
            group_start = up_crosses[i_up]
            if up_crosses[i_up] < down_crosses[i_down] < up_crosses[i_up+1]:
                group_end = down_crosses[i_down]
                group = data_dwant.iloc[group_start:group_end]
                exceedance_groups.append(group)
            else:
                i_up -= 1
            i_up += 1
            i_down += 1
        group_start_final = up_crosses[np.argmin(np.abs(down_crosses[-1]-up_crosses))]
        group_end_final = down_crosses[-1]
        group_final = data_dwant.iloc[group_start_final:group_end_final]
        exceedance_groups.append(group_final)

        #if np.where(up_crosses)[0][0]<np.where(down_crosses)[0][0]:
        #    exceedance_groups = [data_dwant.iloc[np.where(up_crosses)[0][i]:np.where(down_crosses)[0][i]] for i in range(len(np.where(up_crosses)[0]))]
          
        # For each exceedance, interpolate to before first and after last points to find precise exceedance time, and make a nice DataFrame with results #  
        c = -1
        for group in exceedance_groups:
            c += 1
            up_cross_df_i = data_dwant.iloc[group.index[0]-1:group.index[0]+1].resample('1min',kind='timestamp',on='time').mean().interpolate().reset_index()
            up_cross_time = up_cross_df_i.iloc[(up_cross_df_i['val'] - threshold).abs().argmin()]['time']        
            down_cross_df_i = data_dwant.iloc[group.index[-1]:group.index[-1]+2].resample('1min',kind='timestamp',on='time').mean().interpolate().reset_index()
            down_cross_time = down_cross_df_i.iloc[(down_cross_df_i['val'] - threshold).abs().argmin()]['time']
            peak_time = group['time'].iloc[group['val'].argmax()]
            try:
                tide_type = high_lows[high_lows['time'] == peak_time]['tide type'].values[0]
            except IndexError:
                tide_type = 'Unknown'
            row = pd.DataFrame({'Peak Date/Time':group['time'].iloc[group['val'].argmax()],
                                'Period Start':[up_cross_time],
                                'Period End':down_cross_time,
                                'Duration (hours)':(down_cross_time - up_cross_time).total_seconds()/60/60,
                                'Maximum Elevation Above Threshold':group['val'].max() - threshold,
                                'Maximum Elevation (MHHW)':data_mhhw['val'][data_mhhw['time'] == group['time'].iloc[group['val'].argmax()]].values[0],
                                'Tide Type':tide_type})
            if c == 0:
                inundations = row
            else:
                inundations = pd.concat([inundations,row],ignore_index=True) 

        return Out(inundations, threshold_datum, units, {'threshold' : threshold , 'threshold_datum' : threshold_datum , 'data' : data , 'datums' : datums})
    

    
    
