"""
RainEvent_Toolbox.py
"""

#import the python packages which these functions use
import pandas as pd
import xarray as xr
import numpy as np
import datetime as dt
import scipy
import copy
from matplotlib import pyplot as plt

def find_rain_events(dataset, min_duration, min_separation, threshold, noise_floor, front_pad, end_pad):
    '''
    This function splits an xarray dataset into individual datasets for each precipitation event found in the timeseries, based on the input parameters.

    Inputs:
        dataset        - xarray dataset containing a precipitation field labeled 'P'
        min_duration   - minimum consecutive timesteps of rainfall required to be considered an event
        min_separation - minimum consecutive timesteps of no rainfall required to be considered an event end
        threshold      - minimum peak rainfall rate - if the event doesn't get past threshold, it is thrown out
        noise_floor    - maximum value to treat as a zero (prevents very long tails on events where precipitation was juuuust above zero)
        front_pad      - how many minutes prior to rain onset to include in the rain event
        end_pad        - how many minutes after rain end to include in the rain event

    Outputs:
        rain_event_list - a list of xarray datasets, each of which contains one rain event and has the same format as 'dataset' and contains the following additional metadata:
            Rain Event # - monotonic event counter for ease of referring to events
            Rain Onset   - time of first rainfall in event
            Rain End     - time of first zero after event
            Peak Time    - time of peak rainfall
            Peak Rate    - rate of peak rainfall
    Plots:
        This function plots all of the events in a giant row for human review, with the onset and offset times marked by green and red circles
    '''

    rain_rate = dataset.P
    rain_event_list = []
    rain_event_counter = 0
    timestep = 0
    while timestep < len(rain_rate):
        #if the current rain rate isn't zero:
        if (rain_rate.isel(index=timestep).values > noise_floor):
            #then we've got ourselves the start of a rain event!
            rain_start = timestep-1 #preceding zero value is the start, for consistency with Deb's approach
            raining = True
            #enter a new loop that finds the end index of this rain event
            while raining:
                #if the current rain rate is zero...AND...it's going to be zero from now until 'min_separation' from now:
                if (rain_rate.isel(index=timestep).values <= noise_floor and
                    np.mean(rain_rate.isel(index=slice(timestep,timestep+min_separation)).values) <= noise_floor):
                        #then we've located the end of the rain event
                        rain_end = timestep
                        raining = False
                        #create new dataarray containing solely this rain event, plus some extra timesteps on either side
                        rainevent = dataset.isel(index=slice(rain_start-front_pad,rain_end+end_pad))
                        #calculate peak rainfall rate - if you encounter any issues, set the rate to zero
                        try:
                            peak_rate = dataset.isel(index=slice(rain_start,rain_end)).P.max().values.item()
                        except ValueError:
                            peak_rate = 0
                        #check if the peak of the rain event exceeds 'threshold' and the rain event is long enough
                        if (peak_rate >= threshold and (rain_end - rain_start) >= min_duration):
                            #if so, fill out metadata and add this to the master list of rain events
                            rain_event_counter = rain_event_counter + 1
                            rainevent.attrs['Rain Event #'] = rain_event_counter
                            rainevent.attrs['Rain Onset'] = rain_rate.isel(index=rain_start).index.values
                            rainevent.attrs['Rain End'] = rain_rate.isel(index=rain_end).index.values
                            rainevent.attrs['Peak Rate'] = peak_rate
                            rainevent.attrs['Peak Time'] = rainevent.P.where(rainevent.P==peak_rate,drop=True).index.values[0]
                            rain_event_list.append(rainevent)
                            #(if not, just discard and move on)
                #as long as it keeps raining, just increment the timestep and start over
                else:
                    timestep = timestep + 1
            #[spits us out at the timestep after the event ends]
        #if we haven't found a rain event, just increment the timestep and start over
        else:
            timestep = timestep + 1

    #plot all detected rain events for human review (comment out for speed)
    '''
    figlength = 5*len(rain_event_list)
    fig, axx = plt.subplots(nrows=1, ncols=len(rain_event_list),facecolor='w',figsize=(figlength,3))

    for event_num in np.arange(0,len(rain_event_list)):
        rain_event_list[event_num].P.plot.line('-o',ax=axx[event_num],markersize=3,fillstyle=None)
        axx[event_num].set_title(f'Rain Event # {event_num+1}')
        start = rain_event_list[event_num].attrs['Rain Onset']
        axx[event_num].plot(start,rain_event_list[event_num].P.sel(index=start),'.g',markersize=12,fillstyle=None)
        end = rain_event_list[event_num].attrs['Rain End']
        axx[event_num].plot(end,rain_event_list[event_num].P.sel(index=end),'.r',markersize=12,fillstyle=None)
    plt.tight_layout()
    '''

    return rain_event_list

###########################################################################################

def sst_rain_response(rain_event_list, sst, pre_onset_averaging):
    '''
    This function takes a list of rain events and a dataset of sst and slices out the sst for each rain event, calculating the max sst deviation over the rain event along the way.
    It then plots the rain rate & sst for each event, plus a second subplot of heat fluxes and wind speeds.
    The plot includes markers for rain onset (green), peak (blue), and end (red) times, as well as a purple cross at the max SST deviation and a black line with error bars, located at the start of the pre-onset averaging period and indicating the mean+std of SST prior to rain onset.

    Inputs:
        rain_event_list     - a list of xarray datasets, each of which contains one rain event
        sst                 - an xarray dataarray of sst values that (hopefully) covers the time period of all the rain events
        pre_onset_averaging - how many minutes before rain onset to use as the average SST from which to calculate SST deviation

    Outputs:
        sst_event_list      - a list of xarray datasets containing sst data for the same rain events as in rain_event_list. It contains:
            data:
                sst              - timeseries of sst during the event
                δsst             - timeseries of sst deviation from pre-onset mean
            metadata:
                pre-onset mean   - mean of sst for 'pre_onset_averaging' minutes prior to rain onset
                pre-onset std    - std of sst for 'pre_onset_averaging' minutes prior to rain onset
                Max δsst         - maximum deviation of sst from pre-onset mean during rain event
                Time of max δsst - time of maximum deviation
                Time of sst recovery - time the sst returns to within a std of the pre-onset mean
                Minutes to sst recovery - number of minutes between max δsst and sst recovery
                Rain Event #     - monotonic event counter, taken directly from the rain_event_list

    Plots:
        This function plots all of the events in a giant row for human review, with the rain rate & sst for each event in the upper panel, plus a lower panel of heat fluxes and wind speeds during the event.
        The plots include markers for rain onset (green), peak (blue), and end (red) times, as well as a purple cross at the max SST deviation and a black line with error bars, located at the start of the pre-onset averaging period and indicating the mean+std of SST prior to rain onset.
    '''

    sst_event_list = []
    #initialize figure to plot into based on how many rain events we have
    figlength = 5*len(rain_event_list)
    fig, axx = plt.subplots(nrows=2, ncols=len(rain_event_list),facecolor='w',figsize=(figlength,6))

    #cycle through each rain event
    for event_num in np.arange(0,len(rain_event_list)):
        #extract start and end times
        start = pd.to_datetime(rain_event_list[event_num].attrs['Rain Onset'])
        end = pd.to_datetime(rain_event_list[event_num].attrs['Rain End'])
        first = rain_event_list[event_num].index[0]
        last = rain_event_list[event_num].index[-1]

        #-----CREATE SST DATASET------------
        #slice out sst from this time and make new dataset
        sst_event = xr.Dataset()
        sst_event['sst'] = sst.sel(index=slice(first,last))

        #calculate mean sst of previous 'pre_onset_averaging' minutes
        pre_onset = start - dt.timedelta(minutes = pre_onset_averaging)
        sst_event.attrs['pre-onset time'] = pre_onset
        sst_event.attrs['pre-onset mean'] = sst_event.sst.sel(index=slice(pre_onset,start)).mean().item()
        sst_event.attrs['pre-onset std'] = sst_event.sst.sel(index=slice(pre_onset,start)).std().item()

        #calculate max sst deviation from pre-onset mean, time of deviation, and time of return to mean
        sst_event['δsst'] = sst_event.sst - sst_event.attrs['pre-onset mean']
        try:
            #(if there is no data in this slice of sst, these lines will throw an error. But I want it to just pass a dataset full of nans instead.)
            sst_event.attrs['Max δsst'] = sst_event.δsst.sel(index=slice(start,end)).min().values #assumes that a rain event leads to a reduction in SST
            sst_event.attrs['Time of max δsst'] = sst_event.δsst.where(sst_event.δsst==sst_event.attrs['Max δsst'], drop=True).index.values[0]
            
        except IndexError:
            sst_event.attrs['Max δsst'] = np.NaN
            sst_event.attrs['Time of max δsst'] = np.datetime64("NaT")
            sst_event.attrs['Time of sst recovery'] = np.datetime64("NaT")
            sst_event.attrs['Minutes to sst recovery'] = np.NaN
            
        try: 
            #find first time after the max sst deviation that it returns to within 'pre-onset std' of zero, use 5min rolling mean to avoid catching a single outlier too early
            sst_event.attrs['Time of sst recovery'] = sst_event.δsst.rolling(index=5,center=True).mean().where(sst_event.δsst.sel(index=slice(sst_event.attrs['Time of max δsst'],sst_event.index[-1])) > -sst_event.attrs['pre-onset std'], drop=True).index.values[0]
            sst_event.attrs['Minutes to sst recovery'] = (sst_event.attrs['Time of sst recovery'] - sst_event.attrs['Time of max δsst']).astype('timedelta64[m]').astype('float')
        except IndexError:
            sst_event.attrs['Time of sst recovery'] = np.datetime64("NaT")
            sst_event.attrs['Minutes to sst recovery'] = np.NaN

        #add in rain event number to event metadata, and add event to list
        sst_event.attrs['Rain Event #'] = rain_event_list[event_num].attrs['Rain Event #']
        sst_event_list.append(sst_event)

        #----TOP PLOT: Precip & SST-------
        #plot precipitation rate
        rain_event_list[event_num].P.plot.line('-o',ax=axx[0,event_num],markersize=3,fillstyle=None)
        axx[0,event_num].set_ylabel('Rain Rate [mm/hr]',color='C0')
        #plot rainfall start and end times
        axx[0,event_num].plot(start,rain_event_list[event_num].P.sel(index=start),'.g',markersize=12,fillstyle=None)
        axx[0,event_num].plot(end,rain_event_list[event_num].P.sel(index=end),'.r',markersize=12,fillstyle=None)
        #plot rainfall peak
        axx[0,event_num].plot(rain_event_list[event_num].attrs['Peak Time'],rain_event_list[event_num].attrs['Peak Rate'],'.b',markersize=12,fillstyle=None)

        #plot SST
        ax2 = axx[0,event_num].twinx()
        sst_event.sst.plot.line('C1',ax=ax2,fillstyle=None)
        ax2.set_ylabel('SST [$^\circ$C]',color='C1')
        #plot max δSST and pre-onset mean + std
        try:
            ax2.plot(sst_event.attrs['Time of max δsst'],sst_event.sst.sel(index=sst_event.attrs['Time of max δsst']),'x',color='darkmagenta',markersize=12)
            ax2.errorbar(pre_onset,sst_event.attrs['pre-onset mean'],yerr=sst_event.attrs['pre-onset std'],fmt='+k',ecolor='k',capsize=10)
        except KeyError:
            None
            
        try:
            ax2.plot(sst_event.attrs['Time of sst recovery'],sst_event.sst.sel(index=sst_event.attrs['Time of sst recovery']),'^',color='purple',markersize=12)
        except KeyError:
            None
            
        #title
        axx[0,event_num].set_title(f'Rain Event # {event_num+1}')

        #-----BOTTOM PLOT: MJO Phase----
        rain_event_list[event_num].mjo_phase.plot.line('-',color='purple',ax=axx[1,event_num],markersize=3,fillstyle=None)
        axx[1,event_num].set_ylabel('MJO Phase', color='purple')
        axx[1,event_num].set_ylim([0,9])
        axx[1,event_num].grid()
        
        #-----BOTTOM PLOT: Winds & Heat Fluxes----
        #rain_event_list[event_num].lhf.plot.line('-o',color='purple',ax=axx[1,event_num],markersize=3,fillstyle=None)
        #rain_event_list[event_num].shf.plot.line('-o',color='firebrick',ax=axx[1,event_num],markersize=3,fillstyle=None)
        #axx[1,event_num].set_ylabel('Latent Heat Flux [W/m^2]', color='purple')
        #ax1 = axx[1,event_num].twinx()
        #rain_event_list[event_num].U10.plot.line('-o',color='slategray',ax=ax1,markersize=3,fillstyle=None)
        #ax1.set_ylabel('Wind Speed [m/s]',color='slategray')

    plt.tight_layout()

    return sst_event_list

##################################################################################

def plot_rain_events(rain_event_list, sst_event_list, rain_ylims, δsst_ylims, rhf_ylims, wind_ylims):
    '''
    This function is a plotter for once you have developed lists of rain events for both revelle data and sst data.
    The key difference from the plots made by the other functions is that these are on standardized y-axes.
    At the moment, it plots rain rate + δsst on the top plot, and rain heat flux + wind speed on the bottom plot.
    It allows you to specify the universal y limits for each variable as an input (2-element list).
    '''
    #initialize figure to plot into based on how many rain events we have
    figlength = 5*len(rain_event_list)
    fig, axx = plt.subplots(nrows=2, ncols=len(rain_event_list),facecolor='w',figsize=(figlength,6))

    for event_num in np.arange(0,len(rain_event_list)):

        start = pd.to_datetime(rain_event_list[event_num].attrs['Rain Onset'])
        end = pd.to_datetime(rain_event_list[event_num].attrs['Rain End'])

        #----TOP PLOT: Precip & δSST-------
        #plot precipitation rate
        rain_event_list[event_num].P.plot.line('-o',ax=axx[0,event_num],markersize=3,fillstyle=None)
        axx[0,event_num].set_ylabel('Rain Rate [mm/hr]',color='C0')
        #plot rainfall start and end times
        axx[0,event_num].plot(start,rain_event_list[event_num].P.sel(index=start),'.g',markersize=12,fillstyle=None)
        axx[0,event_num].plot(end,rain_event_list[event_num].P.sel(index=end),'.r',markersize=12,fillstyle=None)
        #plot rainfall peak
        axx[0,event_num].plot(rain_event_list[event_num].attrs['Peak Time'],rain_event_list[event_num].attrs['Peak Rate'],'.b',markersize=12,fillstyle=None)
        #set ylim
        axx[0,event_num].set_ylim(rain_ylims)

        #plot SST
        ax2 = axx[0,event_num].twinx()
        sst_event_list[event_num].δsst.plot.line('C1',ax=ax2,fillstyle=None)
        ax2.set_ylabel('$\delta$SST [$^\circ$C]',color='C1')
        #plot max δSST and pre-onset mean + std
        try:
            ax2.plot(sst_event_list[event_num].attrs['Time of max δsst'],sst_event_list[event_num].δsst.sel(index=sst_event_list[event_num].attrs['Time of max δsst']),'x',color='darkmagenta',markersize=12)
            ax2.errorbar(sst_event_list[event_num].attrs['pre-onset time'],0,yerr=sst_event_list[event_num].attrs['pre-onset std'],fmt='+k',capsize=10)
            ax2.plot(sst_event_list[event_num].attrs['Time of sst recovery'],sst_event_list[event_num].δsst.sel(index=sst_event_list[event_num].attrs['Time of sst recovery']),'^',color='purple',markersize=12)
        except KeyError:
            None
        #set ylim
        ax2.set_ylim(δsst_ylims)
        #title
        axx[0,event_num].set_title(f'Rain Event # {event_num+1}')

        #-----BOTTOM PLOT: Winds & Sensible Heat Flux due to Rain----
        rain_event_list[event_num].rhf.plot.line('-o',color='purple',ax=axx[1,event_num],markersize=3,fillstyle=None)
        axx[1,event_num].set_ylabel('Sensible Heat Flux Due to Rain [W/m^2]', color='purple')
        axx[1,event_num].set_ylim(rhf_ylims)
        ax1 = axx[1,event_num].twinx()
        rain_event_list[event_num].U10.plot.line('-o',color='slategray',ax=ax1,markersize=3,fillstyle=None)
        ax1.set_ylabel('Wind Speed [m/s]',color='slategray')
        ax1.set_ylim(wind_ylims)

    plt.tight_layout()

################################################################################

def extract_event_characteristics(rain_event_list, sst_event_list):
    '''
    This function takes the lists of rain events and corresponding sst events and extracts the metadata fields for each event, returning a new xarray dataarray with event number as the coordinate.

    Inputs:
        rain_event_list     - a list of xarray datasets, each of which contains one rain event
        sst_event_list      - a list of xarray datasets, each of which contains the sst data corresponding to the rain events in rain_event_list

    Outputs:
        rain_events_summary - an xarray dataset with event number as the shared coordinate, and the following variables:
            δsst_max        - Maximum SST deviation from mean SST prior to rain onset (C)
            t_δsst_max      - Time of Max SST deviation, in minutes after rain onset
            rain_max        - Maximum rainfall rate (mm/hr)
            t_rain_max      - Time of maximum rainfall rate, in minutes after rain onset
            L_rain          - Length of rain event, in minutes
            rain_total      - Total rainfall (mm)
    '''
    #initialize empty lists
    δsst_max = []
    t_δsst_max = []
    rain_max = []
    t_rain_max = []
    L_rain = []
    rain_total = []

    #loop through each event and extract characteristics from metadata
    for event_num in np.arange(0,len(rain_event_list)):
        #max δSST
        try:
            δsst_max.append(sst_event_list[event_num].attrs['Max δsst'].item())
        except AttributeError:
            δsst_max.append(np.nan)
        if np.isnat(sst_event_list[event_num].attrs['Time of max δsst']) == False:
            t_δsst_max.append((sst_event_list[event_num].attrs['Time of max δsst'] - rain_event_list[event_num].attrs['Rain Onset']).astype('timedelta64[m]').astype('float'))
        else:
            t_δsst_max.append(np.nan)
        #max rain
        rain_max.append(rain_event_list[event_num].attrs['Peak Rate'])
        #time of max rain
        t_rain_max.append((rain_event_list[event_num].attrs['Peak Time'] - rain_event_list[event_num].attrs['Rain Onset']).astype('timedelta64[m]').astype('float'))
        #length of rain event
        L_rain.append((rain_event_list[event_num].attrs['Rain End'] - rain_event_list[event_num].attrs['Rain Onset']).astype('timedelta64[m]').astype('float'))
        #cumulative rain
        rain_total.append(rain_event_list[event_num].Precip[-1].values - rain_event_list[event_num].Precip[0].values)

    #create xarray dataarray of metadata for each event
    rain_events_summary = xr.Dataset(data_vars = {'δsst_max':('event_num',δsst_max),
                                                  't_δsst_max':('event_num',t_δsst_max),
                                                  'rain_max':('event_num',rain_max),
                                                  't_rain_max':('event_num',t_rain_max),
                                                  'L_rain':('event_num',L_rain),
                                                  'rain_total':('event_num',rain_total)
                                                 },coords={'event_num':np.arange(1,len(sst_event_list)+1)})
    rain_events_summary.attrs['experiment'] = rain_event_list[0].attrs['experiment']
    rain_events_summary.δsst_max.attrs['long_name'] = 'Max δSST'
    rain_events_summary.δsst_max.attrs['units'] = 'C'
    rain_events_summary.t_δsst_max.attrs['long_name'] = 'Time of Max δSST'
    rain_events_summary.t_δsst_max.attrs['units'] = 'minutes after rain onset'
    rain_events_summary.rain_max.attrs['long_name'] = 'Max Rain Rate'
    rain_events_summary.rain_max.attrs['units'] = 'mm/hr'
    rain_events_summary.t_rain_max.attrs['long_name'] = 'Time of Max Rain Rate'
    rain_events_summary.t_rain_max.attrs['units'] = 'minutes after rain onset'
    rain_events_summary.L_rain.attrs['long_name'] = 'Length of Rain Event'
    rain_events_summary.L_rain.attrs['units'] = 'minutes'
    rain_events_summary.L_rain.attrs['long_name'] = 'Total Rainfall'
    rain_events_summary.L_rain.attrs['units'] = 'mm'

    return rain_events_summary

###########################################################################################

def plot_histograms(rain_events_summary, nbins):
    '''
    Simple function for plotting histograms of the 6 major variables contained in rain_events_summary.
    '''

    fig,axx = plt.subplots(nrows=2,ncols=3,facecolor='w',figsize=(8,5))

    rain_events_summary.δsst_max.plot.hist(ax=axx[1,0],bins=nbins);
    axx[1,0].set_title('Max δSST')
    axx[1,0].set_xlabel('$^\circ$C')
    rain_events_summary.t_δsst_max.plot.hist(ax=axx[1,1],bins=nbins);
    axx[1,1].set_title('Time of Max δSST')
    axx[1,1].set_xlabel('minutes after rain onset')
    rain_events_summary.rain_max.plot.hist(ax=axx[0,0],bins=nbins);
    axx[0,0].set_title('Max Rain Rate')
    axx[0,0].set_xlabel('mm/hr')
    rain_events_summary.t_rain_max.plot.hist(ax=axx[0,1],bins=nbins);
    axx[0,1].set_title('Time of Max Rain Rate')
    axx[0,1].set_xlabel('minutes after rain onset')
    rain_events_summary.L_rain.plot.hist(ax=axx[0,2],bins=nbins);
    axx[0,2].set_title('Length of Rain Event')
    axx[0,2].set_xlabel('minutes')
    rain_events_summary.rain_total.plot.hist(ax=axx[1,2],bins=nbins)
    axx[1,2].set_title('Total Rainfall')
    axx[1,2].set_xlabel('mm')

    plt.tight_layout()

##############################################################################################

def scatterplot_with_linreg(axis, x, xlabel, y, ylabel):
    #run linear regression and plot
    coef = np.polyfit(x,y,1)
    slope, intercept = coef[0], coef[1]
    r_value = np.corrcoef(x, y)[0, 1]
    
    poly1d_fn = np.poly1d(coef) # poly1d_fn is now a function which takes in x and returns an estimate for y
    axis.plot(x, poly1d_fn(x), 'C1', linewidth=1, zorder=2)
    axis.annotate('$R^2$ = %0.2f' % r_value**2, xy=(0.05, 0.89), xycoords='axes fraction')

    #scatterplot and label
    axis.scatter(x, y, s=5, zorder=1)
    axis.set_xlabel(xlabel)
    axis.set_ylabel(ylabel)

##############################################################################################

def plot_relationships(rain_events_summary):
    '''
    Simple function for plotting linear relationships between the 6 major variables in rain_event_summary
    '''

    idx = np.isfinite(rain_events_summary.δsst_max.values)

    δsst_max = rain_events_summary.δsst_max[idx]
    t_δsst_max = rain_events_summary.t_δsst_max[idx]
    rain_max = rain_events_summary.rain_max[idx]
    t_rain_max = rain_events_summary.t_rain_max[idx]
    L_rain = rain_events_summary.L_rain[idx]
    rain_total = rain_events_summary.rain_total[idx]

    fig,axx = plt.subplots(nrows=3,ncols=5,facecolor='w',figsize=(14,8))
    #row 1: Max δSST vs. other variables
    scatterplot_with_linreg(axx[0,0], δsst_max, 'Max δSST', rain_total, 'Total Rainfall [mm]')
    scatterplot_with_linreg(axx[0,1], δsst_max, 'Max δSST', L_rain, 'Length of Rain Event [minutes]')
    scatterplot_with_linreg(axx[0,2], δsst_max, 'Max δSST', rain_max, 'Max Rain Rate [mm/hr]')
    scatterplot_with_linreg(axx[0,3], δsst_max, 'Max δSST', t_rain_max, 'Minutes to Max Rain')
    scatterplot_with_linreg(axx[0,4], δsst_max, 'Max δSST', t_δsst_max, 'Minutes to Max δSST')

    scatterplot_with_linreg(axx[1,0], t_δsst_max, 'Minutes to Max δSST', rain_total, 'Total Rainfall [mm]')
    scatterplot_with_linreg(axx[1,1], t_δsst_max, 'Minutes to Max δSST', L_rain, 'Length of Rain Event [minutes]')
    scatterplot_with_linreg(axx[1,2], t_δsst_max, 'Minutes to Max δSST', rain_max, 'Max Rain Rate [mm/hr]')
    scatterplot_with_linreg(axx[1,3], t_δsst_max, 'Minutes to Max δSST', t_rain_max, 'Minutes to Max Rain')
    scatterplot_with_linreg(axx[1,4], t_δsst_max, 'Minutes to Max δSST', δsst_max, 'Max δSST')

    scatterplot_with_linreg(axx[2,0], t_rain_max, 'Minutes to Max Rain', rain_total, 'Total Rainfall [mm]')
    scatterplot_with_linreg(axx[2,1], t_rain_max, 'Minutes to Max Rain', L_rain, 'Length of Rain Event [minutes]')
    scatterplot_with_linreg(axx[2,2], t_rain_max, 'Minutes to Max Rain', rain_max, 'Max Rain Rate [mm/hr]')
    scatterplot_with_linreg(axx[2,3], t_rain_max, 'Minutes to Max Rain', t_δsst_max, 'Minutes to Max δSST')
    scatterplot_with_linreg(axx[2,4], t_rain_max, 'Minutes to Max Rain', δsst_max, 'Max δSST')

    plt.tight_layout()

##########################################################################################

def calculate_deltas(rain_event_list, param_list, pre_onset_averaging):
    '''
    Function for calculating the deviation from pre-onset mean of any parameter contained in the rain events in rain_event_list
    '''

    #cycle through rain events
    for event_num in np.arange(0,len(rain_event_list)):
        #extract pre-onset averaging time
        start = pd.to_datetime(rain_event_list[event_num].attrs['Rain Onset'])
        pre_onset = start - dt.timedelta(minutes = pre_onset_averaging)
        rain_event_list[event_num].attrs['pre-onset time'] = pre_onset

        #cycle through list of parameters to calculate deltas for
        for param in param_list:
            try:
                rain_event_list[event_num][param].attrs['pre-onset mean'] = rain_event_list[event_num][param].sel(index=slice(pre_onset,start)).mean().item()
                rain_event_list[event_num][param].attrs['pre-onset std'] = rain_event_list[event_num][param].sel(index=slice(pre_onset,start)).std().item()
                rain_event_list[event_num][f'δ{param}'] = rain_event_list[event_num][param] - rain_event_list[event_num][param].attrs['pre-onset mean']
            except KeyError:
                print(f'{param} is not a valid variable name in rain_event_list')
                param_list.remove(param)

    return rain_event_list

##################################################################################################

def extract_composite_event(rain_event_list, sst_event_list, param_list, start, stop, spacing, plotflag, title):
    '''
    This function takes all of the rain events and creates a single composite event by normalizing the timescale.
    Timescale normalization is calculated as (current time - rain onset time)[minutes]/(minutes to max δSST).
    The function then plots the composite event for 4 variables (air temp, sst, bulk temp, bulk salinity).

    Inputs:
        rain_event_list  - list of xarray datasets for each rain event
        sst_event_list   - list of xarray datasets containing sst information for each rain event
        param_list       - list of data variable names for which there is a δ{param} to be found in rain_event_list (should be same as the param_list input to calculate_deltas)
        start            - desired start time in normalized coordinates (-2 recommended)
        stop             - desired end time in normalized coordinates (8-12 recommended)
        spacing          - step size ("bin size") for normalized time coordinates

    Outputs:
        composite_event  - xarray dataset with two dimensions - event number, and normalized timescale

    Also puts out a plot of the composite events in air temp, sst, bulk temp, and bulk salinity.
    '''

    #initialize array of resampling coordinates
    resample_coords = np.linspace(start, stop, num=int((stop-start)/spacing),endpoint=False)

    #create new lists of datasets to do time-normalization on
    sst_normed_list = copy.deepcopy(sst_event_list)
    rain_normed_list = copy.deepcopy(rain_event_list)

    #initialize master arrays (shaped by event # vs. normalized timescale)
    #sst
    normed_δsst = np.empty(shape=(len(rain_event_list),len(resample_coords)))
    #rain_event variables
    array_list = []
    for param in param_list:
        array_list.append(np.empty(shape=(len(rain_event_list),len(resample_coords))))

    #cycle through rain events
    for event_num in np.arange(0,len(rain_event_list)):
        #find times of rain onset and max sst response
        onset = rain_event_list[event_num].attrs['Rain Onset']
        t_max_δsst = (sst_event_list[event_num].attrs['Time of max δsst'] - onset).astype('timedelta64[m]').astype('float64')

        #normalized time = (current time - rain onset time)[minutes]/(minutes to max δSST)
        rain_normtime = (rain_event_list[event_num].index - onset).values.astype('timedelta64[s]').astype('float64')/60/t_max_δsst
        sst_normtime  = ((sst_event_list[event_num].index - onset).values.astype('timedelta64[s]').astype('float64')/60/t_max_δsst)

        #replace index with normalized time variable in the new list
        sst_normed_list[event_num]['index'] = sst_normtime
        rain_normed_list[event_num]['index'] = rain_normtime

        #interpolate data onto resampling coordinates
        sst_interp = sst_normed_list[event_num].interp(index=resample_coords)
        rain_interp = rain_normed_list[event_num].interp(index=resample_coords)

        #insert data into row of master array
        normed_δsst[event_num,:] = sst_interp.δsst.values
        for idx in np.arange(0,len(param_list)):
            param = param_list[idx]
            array_list[idx][event_num,:] = rain_interp[f'δ{param}'].values

    #create and fill new dataset containing normalized events
    composite_event = xr.Dataset(coords={'index':resample_coords,
                                         'event':np.arange(1,len(rain_event_list)+1)})
    #sst
    composite_event = composite_event.assign({'δsst':(['event','index'],normed_δsst)})
    #rain event variables
    for idx in np.arange(0,len(param_list)):
        param = param_list[idx]
        composite_event = composite_event.assign({f'δ{param}':(['event','index'],array_list[idx])})

    #count the number of data points at each timestep to size plots
    sizes = np.empty_like(resample_coords)
    sizes_DC = np.empty_like(resample_coords)
    for idx in np.arange(0,len(resample_coords)):
        #count how many non-nan entries there are at this timestep
        sizes[idx] = composite_event.δsst.sel(index=resample_coords[idx]).count().values.item()
        sizes_DC[idx] = composite_event.δHsC.sel(index=resample_coords[idx]).count().values.item()
        
    #takes means and stds across all events
    means = composite_event.mean(dim='event',skipna=True)
    stds = composite_event.std(dim='event',skipna=True)
    
    means_DC = composite_event.where(~np.isnan(composite_event.δHsC)).mean(dim='event',skipna=True)
    stds_DC = composite_event.where(~np.isnan(composite_event.δHsC)).mean(dim='event',skipna=True)

    if plotflag == 1:
        #plot
        fig,axx= plt.subplots(nrows=4,ncols=2,figsize=(16,12),facecolor='w')
        ticks = np.arange(start,stop+1,1)
        ticklabels = ticks.astype('str').tolist()
        ticklabels[-start] = 'Onset'
        ticklabels[-start + 1] = r'$ \delta SST_{max} $'
        plt.suptitle(title, y=0.98, fontsize=18)

        temp_ylims = [-1,0.2]
        panel_x = 0.995
        panel_y = 0.97
        panel_font = 14
        title_font = 16

        #----------LEFT SIDE----------------
        #Subplot 0,0: Rain Rate
        axx[0,0].hlines(y=0,xmin=start-1,xmax=stop+1,color='k',linewidth = 0.5,zorder=0)
        axx[0,0].grid(alpha=0.2)
        axx[0,0].errorbar(x=resample_coords, y=means.δPrecip, yerr=stds.δPrecip, fmt='-k', ecolor='k', linewidth=0.5, capsize=10*spacing, zorder=100)
        sc = axx[0,0].scatter(x=resample_coords, y=means.δPrecip, s=sizes, c='C0',zorder=200)
        axx[0,0].set_xlim([resample_coords[0],resample_coords[-1]])
        axx[0,0].set_ylim([-12,31])
        axx[0,0].set_xticks(ticks)
        axx[0,0].set_xticklabels(ticklabels, fontsize=panel_font)
        #axx[0,0].set_xlabel(r'Normalized Time $ (t-t_{onset})/t_{δSST_{max}}$', fontsize=panel_font)
        axx[0,0].set_ylabel('Rain Rate $(mm/hr)$', fontsize=panel_font)
        axx[0,0].set_title('Rain Rate',fontsize=title_font)
        axx[0,0].legend(*sc.legend_elements("sizes", num=6),title='# of Events',loc='upper left', framealpha=0.9).set_zorder(1000)
        axx[0,0].text(panel_x, panel_y, '(a)', fontsize=panel_font, horizontalalignment='right', verticalalignment='top', transform=axx[0,0].transAxes, zorder=300)
        axx[0,0].yaxis.set_tick_params(labelsize=panel_font)
        

        
        #Subplot 1,0: δSST
        axx[1,0].hlines(y=0,xmin=start-1,xmax=stop+1,color='k',linewidth = 0.5,zorder=0)
        axx[1,0].grid(alpha=0.2)
        axx[1,0].errorbar(x=resample_coords, y=means.δsst, yerr=stds.δsst, fmt='-k', ecolor='k', linewidth=0.5, capsize=10*spacing, zorder=100)
        h1 = axx[1,0].scatter(x=resample_coords, y=means.δsst, s=sizes, c='red',zorder=200)
        h2 = axx[1,0].scatter(x=resample_coords, y=means.δSST, s=sizes*0.8, marker='v', facecolors='none', edgecolors='darkred' ,zorder=150) #plots COARE skin-corrected SST
        h3 = axx[1,0].scatter(x=resample_coords, y=means.δt_int_vF, s=sizes, marker='o', facecolors='none', edgecolors='darkmagenta' ,zorder=150) #plots Bellenger model (Original)
        h4 = axx[1,0].scatter(x=resample_coords, y=means.δt_int, s=sizes+2, marker='*', facecolors='none', edgecolors='navy' ,zorder=200) #plots Bellenger model (Modified)

        axx[1,0].set_xlim([resample_coords[0],resample_coords[-1]])
        axx[1,0].set_xticks(ticks)
        axx[1,0].set_xticklabels(ticklabels, fontsize=panel_font)
        #axx[1,0].set_xlabel(r'Normalized Time $ (t-t_{onset})/t_{δSST_{max}}$', fontsize=panel_font)
        axx[1,0].set_ylabel('δSST ($^\circ C$)', fontsize=panel_font)
        axx[1,0].set_ylim(temp_ylims)
        axx[1,0].set_title('Skin SST',fontsize=title_font)
        axx[1,0].legend([h1,h2,h3,h4],['Observed by KT-15 (10μm)', 'COARE3.5 Model', 'Original Bellenger Model', 'Modified Bellenger Model'])
        #axx[1,0].legend([h1,h2],['Observed by KT-15 (10μm)', 'COARE'])
        axx[1,0].text(panel_x, panel_y, '(b)', fontsize=panel_font, horizontalalignment='right', verticalalignment='top', transform=axx[1,0].transAxes, zorder=300)
        axx[1,0].yaxis.set_tick_params(labelsize=panel_font)

        #Subplot 2,0: δT_bulk (sea snake)
        axx[2,0].hlines(y=0,xmin=start-1,xmax=stop+1,color='k',linewidth = 0.5,zorder=0)
        axx[2,0].grid(alpha=0.2)
        axx[2,0].errorbar(x=resample_coords, y=means.δTsea, yerr=stds.δTsea, fmt='-k', ecolor='k', linewidth=0.5, capsize=10*spacing, zorder=100)
        h3 = axx[2,0].scatter(x=resample_coords, y=means.δTsea, s=sizes, c='darkred',zorder=200)
        h4 = axx[2,0].scatter(x=resample_coords, y=means.δTseaTSG, s=sizes, c='k',marker='1',zorder=100)
        h5 = axx[2,0].scatter(x=resample_coords, y=means.δt_int_vF, s=sizes, marker='o', facecolors='none', edgecolors='darkmagenta' ,zorder=200) #plots Bellenger model SST (Original)
        h6 = axx[2,0].scatter(x=resample_coords, y=means.δt_int, s=sizes+5, marker='*', facecolors='none', edgecolors='navy' ,zorder=200) #plots Bellenger model SST
        axx[2,0].legend([h3,h4,h5,h6],['Sea Snake (5cm)', 'Ship TSG (5m)', 'Original Bellenger Model (0m)', 'Modified Bellenger Model (0m)'])

        axx[2,0].set_xlim([resample_coords[0],resample_coords[-1]])
        axx[2,0].set_xticks(ticks)
        axx[2,0].set_xticklabels(ticklabels, fontsize=panel_font)
        #axx[2,0].set_xlabel(r'Normalized Time $ (t-t_{onset})/t_{δSST_{max}}$', fontsize=panel_font)
        axx[2,0].set_ylabel('$\delta T_{bulk} (^\circ C)$', fontsize=panel_font)
        #axx[2,0].set_ylim(temp_ylims)
        axx[2,0].set_title('Bulk Sea Temperature',fontsize=title_font)
        axx[2,0].text(panel_x, panel_y, '(c)', fontsize=panel_font, horizontalalignment='right', verticalalignment='top', transform=axx[2,0].transAxes, zorder=300)
        axx[2,0].yaxis.set_tick_params(labelsize=panel_font)
        
         #Subplot 3,0: Sensible Fluxes
        axx[3,0].hlines(y=0,xmin=start-1,xmax=stop+1,color='k',linewidth = 0.5,zorder=0)
        axx[3,0].grid(alpha=0.2)
        axx[3,0].errorbar(x=resample_coords, y=means_DC.δHsC, yerr=stds.δHsC, fmt='-g', ecolor='g', linewidth=0.5, capsize=10*spacing, zorder=300)
        h3 = axx[3,0].scatter(x=resample_coords, y=means_DC.δHsC, s=sizes_DC, c='forestgreen',zorder=400)
        axx[3,0].errorbar(x=resample_coords, y=means_DC.δshf, yerr=stds.δshf, fmt='-m', ecolor='m', linewidth=0.5, capsize=10*spacing, zorder=100)
        h4 = axx[3,0].scatter(x=resample_coords, y=means_DC.δshf, s=sizes_DC, c='purple',zorder=200)
        axx[3,0].errorbar(x=resample_coords, y=means_DC.δrhf, yerr=stds.δrhf, fmt='-k', ecolor='k', linewidth=0.5, capsize=10*spacing, zorder=100)
        h5 = axx[3,0].scatter(x=resample_coords, y=means_DC.δrhf, s=sizes_DC, c='darkblue',zorder=200)
        axx[3,0].set_xlim([resample_coords[0],resample_coords[-1]])
        axx[3,0].set_xticks(ticks)
        axx[3,0].set_xticklabels(ticklabels, fontsize=panel_font)
        axx[3,0].legend([h3,h4,h5],['DC','COARE','Gosnell Rain'])
        #axx[3,0].legend([h4],['COARE'],loc='upper left')
        axx[3,0].set_xlabel(r'Normalized Time $ (t-t_{onset})/t_{δSST_{max}}$', fontsize=panel_font)
        axx[3,0].set_ylabel('$\delta$ Flux $(W/m^2)$', fontsize=panel_font)
        axx[3,0].set_title('Sensible Fluxes',fontsize=title_font)
        axx[3,0].text(panel_x, panel_y, '(d)', fontsize=panel_font, horizontalalignment='right', verticalalignment='top', transform=axx[3,0].transAxes, zorder=300)
        axx[3,0].yaxis.set_tick_params(labelsize=panel_font)
        
        
        #--------------RIGHT SIDE---------------
        
        #Subplot 0,1: Air Temp
        axx[0,1].hlines(y=0,xmin=start-1,xmax=stop+1,color='k',linewidth = 0.5,zorder=0)
        axx[0,1].grid(alpha=0.2)
        axx[0,1].errorbar(x=resample_coords, y=means.δT02, yerr=stds.δT02, fmt='-k', ecolor='k', linewidth=0.5, capsize=10*spacing, zorder=100)
        axx[0,1].scatter(x=resample_coords, y=means.δT02, s=sizes, c='tomato',zorder=200)
        axx[0,1].set_xlim([resample_coords[0],resample_coords[-1]])
        axx[0,1].set_ylim([-2.3,0.3])
        axx[0,1].set_xticks(ticks)
        axx[0,1].set_xticklabels(ticklabels, fontsize=panel_font)
        #axx[0,1].set_xlabel(r'Normalized Time $ (t-t_{onset})/t_{δSST_{max}}$', fontsize=panel_font)
        axx[0,1].set_ylabel('$\delta T_{air} (^\circ C)$', fontsize=panel_font)
        axx[0,1].set_title('Air Temperature',fontsize=title_font)
        axx[0,1].text(panel_x, panel_y, '(e)', fontsize=panel_font, horizontalalignment='right', verticalalignment='top', transform=axx[0,1].transAxes, zorder=300)
        
        #Subplot 1,1: Wind Speed
        axx[1,1].hlines(y=0,xmin=start-1,xmax=stop+1,color='k',linewidth = 0.5,zorder=0)
        axx[1,1].grid(alpha=0.2)
        axx[1,1].errorbar(x=resample_coords, y=means.δwspd, yerr=stds.δwspd, fmt='-k', ecolor='k', linewidth=0.5, capsize=10*spacing, zorder=100)
        axx[1,1].scatter(x=resample_coords, y=means.δwspd, s=sizes, c='darkslategray',zorder=200)
        axx[1,1].set_xlim([resample_coords[0],resample_coords[-1]])
        axx[1,1].set_xticks(ticks)
        axx[1,1].set_xticklabels(ticklabels, fontsize=panel_font)
        #axx[1,1].set_xlabel(r'Normalized Time $ (t-t_{onset})/t_{δSST_{max}}$', fontsize=panel_font)
        axx[1,1].set_ylabel('$\delta U (m/s)$', fontsize=panel_font)
        axx[1,1].set_title('Wind Speed',fontsize=title_font)
        axx[1,1].text(panel_x, panel_y, '(f)', fontsize=panel_font, horizontalalignment='right', verticalalignment='top', transform=axx[1,1].transAxes, zorder=300)
        axx[1,1].yaxis.set_tick_params(labelsize=panel_font)
        
        #Subplot 2,1: δQ
        axx[2,1].hlines(y=0,xmin=start-1,xmax=stop+1,color='k',linewidth = 0.5,zorder=0)
        axx[2,1].grid(alpha=0.2)
        axx[2,1].errorbar(x=resample_coords, y=means.δQ02, yerr=stds.δQ02, fmt='-k', ecolor='k', linewidth=0.5, capsize=10*spacing, zorder=100)
        axx[2,1].scatter(x=resample_coords, y=means.δQ02, s=sizes, c='darkslateblue',zorder=200)
        axx[2,1].set_xlim([resample_coords[0],resample_coords[-1]])
        axx[2,1].set_xticks(ticks)
        axx[2,1].set_xticklabels(ticklabels, fontsize=panel_font)
        #axx[2,1].set_xlabel(r'Normalized Time $ (t-t_{onset})/t_{δSST_{max}}$', fontsize=panel_font)
        axx[2,1].set_ylabel('$\delta Q_{2m} (g/kg)$', fontsize=panel_font)
        axx[2,1].set_title('Specific Humidity',fontsize=title_font)
        axx[2,1].text(panel_x, panel_y, '(g)', fontsize=panel_font, horizontalalignment='right', verticalalignment='top', transform=axx[2,1].transAxes, zorder=300)
        axx[2,1].yaxis.set_tick_params(labelsize=panel_font)

        
        #Subplot 3,1: Latent fluxes
        #axx[3,1].plot(resample_coords, np.zeros(len(resample_coords)),linewidth = 0.5, zorder=0)
        #axx[3,1].errorbar(x=resample_coords, y=means_DC.δHlC, yerr=stds_DC.δHlC, fmt='-g', ecolor='g', linewidth=0.5, capsize=10*spacing, zorder=100)
        #h1 = axx[3,1].scatter(x=resample_coords, y=means_DC.δHlC, s=sizes_DC, c='forestgreen',zorder=200)
        #axx[3,1].errorbar(x=resample_coords, y=means_DC.δlhf, yerr=stds.δlhf, fmt='-m', ecolor='m', linewidth=0.5, capsize=10*spacing, zorder=100)
        #h2 = axx[3,1].scatter(x=resample_coords, y=means_DC.δlhf, s=sizes_DC, c='purple',zorder=200)
        #axx[3,1].set_xlim([resample_coords[0],resample_coords[-1]])
        #axx[2,0].set_ylim([-80,300])
        #axx[3,1].set_xticks(ticks)
        #axx[3,1].set_xticklabels(ticklabels)
        #axx[3,1].legend([h1,h2],['DC','COARE'])
        #axx[3,1].set_xlabel(r'Normalized Time $ (t-t_{onset})/t_{δSST_{max}}$')
        #axx[3,1].set_ylabel('$\delta$ Flux $(W/m^2)$')
        #axx[3,1].set_title('Latent Fluxes')
        
        
        #Subplot 3,1: Bulk Salinity
        axx[3,1].hlines(y=0,xmin=start-1,xmax=stop+1,color='k',linewidth = 0.5,zorder=0)
        axx[3,1].grid(alpha=0.2)
        axx[3,1].errorbar(x=resample_coords, y=means.δSalTSG, yerr=stds.δSalTSG, fmt='-k', ecolor='k', linewidth=0.5, capsize=10*spacing, zorder=10)
        h1 = axx[3,1].scatter(x=resample_coords, y=means.δSalTSG, s=sizes, c='blue',zorder=20)
        axx[3,1].errorbar(x=resample_coords, y=means.δs_int, yerr=stds.δs_int, fmt='-k', ecolor='k', linewidth=0.5, capsize=10*spacing, zorder=10)
        h2 = axx[3,1].scatter(x=resample_coords, y=means.δs_int, s=sizes, marker='*', facecolor='None', c='k')
        axx[3,1].set_xlim([resample_coords[0],resample_coords[-1]])
        #axx[3,1].set_ylim([-0.18,0.28])
        axx[3,1].set_xticks(ticks)
        axx[3,1].set_xticklabels(ticklabels, fontsize=panel_font)
        axx[3,1].set_xlabel(r'Normalized Time $ (t-t_{onset})/t_{δSST_{max}}$', fontsize=panel_font)
        axx[3,1].set_ylabel('$\delta S_{bulk} (psu)$', fontsize=panel_font)
        axx[3,1].set_title('Salinity',fontsize=title_font)
        axx[3,1].legend([h1,h2],['Ship TSG (5m)','Bellenger Model (0m)'])
        axx[3,1].text(panel_x, panel_y, '(h)', fontsize=panel_font, horizontalalignment='right', verticalalignment='top', transform=axx[3,1].transAxes, zorder=300)
        axx[3,1].yaxis.set_tick_params(labelsize=panel_font)


        #Subplot 2,0: Sensible Heat Flux Due to Rain
        #axx[2,0].plot(resample_coords, np.zeros(len(resample_coords)),linewidth = 0.5, zorder=0)
        #axx[2,0].errorbar(x=resample_coords, y=means.δrhf, yerr=stds.δrhf, fmt='-k', ecolor='k', linewidth=0.5, capsize=10*spacing, zorder=100)
        #axx[2,0].scatter(x=resample_coords, y=means.δrhf, s=sizes, c='forestgreen',zorder=200)
        #axx[2,0].set_xlim([resample_coords[0],resample_coords[-1]])
        #axx[2,0].set_ylim([-80,300])
        #axx[2,0].set_xticks(ticks)
        #axx[2,0].set_xticklabels(ticklabels)
        #axx[2,0].set_xlabel(r'Normalized Time $ (t-t_{onset})/t_{δSST_{max}}$')
        #axx[2,0].set_ylabel('Rain Heat Flux $(W/m^2)$')
        #axx[2,0].set_title('Sensible Heat Flux Due to Rain (COARE)')
        

        plt.tight_layout()

    return composite_event


def extract_composite_event_SS(rain_event_list, sst_event_list, param_list, start, stop, spacing, plotflag, title):
    '''
    **Same function as extract_composite_event but plotting the sea snake-initialized model results rather than the TSG-initialized model results**
    This function takes all of the rain events and creates a single composite event by normalizing the timescale.
    Timescale normalization is calculated as (current time - rain onset time)[minutes]/(minutes to max δSST).
    The function then plots the composite event for 4 variables (air temp, sst, bulk temp, bulk salinity).

    Inputs:
        rain_event_list  - list of xarray datasets for each rain event
        sst_event_list   - list of xarray datasets containing sst information for each rain event
        param_list       - list of data variable names for which there is a δ{param} to be found in rain_event_list (should be same as the param_list input to calculate_deltas)
        start            - desired start time in normalized coordinates (-2 recommended)
        stop             - desired end time in normalized coordinates (8-12 recommended)
        spacing          - step size ("bin size") for normalized time coordinates

    Outputs:
        composite_event  - xarray dataset with two dimensions - event number, and normalized timescale

    Also puts out a plot of the composite events in air temp, sst, bulk temp, and bulk salinity.
    '''

    #initialize array of resampling coordinates
    resample_coords = np.linspace(start, stop, num=int((stop-start)/spacing),endpoint=False)

    #create new lists of datasets to do time-normalization on
    sst_normed_list = copy.deepcopy(sst_event_list)
    rain_normed_list = copy.deepcopy(rain_event_list)

    #initialize master arrays (shaped by event # vs. normalized timescale)
    #sst
    normed_δsst = np.empty(shape=(len(rain_event_list),len(resample_coords)))
    #rain_event variables
    array_list = []
    for param in param_list:
        array_list.append(np.empty(shape=(len(rain_event_list),len(resample_coords))))

    #cycle through rain events
    for event_num in np.arange(0,len(rain_event_list)):
        #find times of rain onset and max sst response
        onset = rain_event_list[event_num].attrs['Rain Onset']
        t_max_δsst = (sst_event_list[event_num].attrs['Time of max δsst'] - onset).astype('timedelta64[m]').astype('float64')

        #normalized time = (current time - rain onset time)[minutes]/(minutes to max δSST)
        rain_normtime = (rain_event_list[event_num].index - onset).values.astype('timedelta64[s]').astype('float64')/60/t_max_δsst
        sst_normtime  = ((sst_event_list[event_num].index - onset).values.astype('timedelta64[s]').astype('float64')/60/t_max_δsst)

        #replace index with normalized time variable in the new list
        sst_normed_list[event_num]['index'] = sst_normtime
        rain_normed_list[event_num]['index'] = rain_normtime

        #interpolate data onto resampling coordinates
        sst_interp = sst_normed_list[event_num].interp(index=resample_coords)
        rain_interp = rain_normed_list[event_num].interp(index=resample_coords)

        #insert data into row of master array
        normed_δsst[event_num,:] = sst_interp.δsst.values
        for idx in np.arange(0,len(param_list)):
            param = param_list[idx]
            array_list[idx][event_num,:] = rain_interp[f'δ{param}'].values

    #create and fill new dataset containing normalized events
    composite_event = xr.Dataset(coords={'index':resample_coords,
                                         'event':np.arange(1,len(rain_event_list)+1)})
    #sst
    composite_event = composite_event.assign({'δsst':(['event','index'],normed_δsst)})
    #rain event variables
    for idx in np.arange(0,len(param_list)):
        param = param_list[idx]
        composite_event = composite_event.assign({f'δ{param}':(['event','index'],array_list[idx])})

    #count the number of data points at each timestep to size plots
    sizes = np.empty_like(resample_coords)
    sizes_DC = np.empty_like(resample_coords)
    for idx in np.arange(0,len(resample_coords)):
        #count how many non-nan entries there are at this timestep
        sizes[idx] = composite_event.δsst.sel(index=resample_coords[idx]).count().values.item()
        sizes_DC[idx] = composite_event.δHsC.sel(index=resample_coords[idx]).count().values.item()
        
    #takes means and stds across all events
    means = composite_event.mean(dim='event',skipna=True)
    stds = composite_event.std(dim='event',skipna=True)
    
    means_DC = composite_event.where(~np.isnan(composite_event.δHsC)).mean(dim='event',skipna=True)
    stds_DC = composite_event.where(~np.isnan(composite_event.δHsC)).mean(dim='event',skipna=True)

    if plotflag == 1:
        #plot
        fig,axx= plt.subplots(nrows=4,ncols=2,figsize=(16,12),facecolor='w')
        ticks = np.arange(start,stop+1,1)
        ticklabels = ticks.astype('str').tolist()
        ticklabels[-start] = 'Onset'
        ticklabels[-start + 1] = r'$ \delta SST_{max} $'
        plt.suptitle(title, y=0.98, fontsize=18)

        temp_ylims = [-1,0.2]
        panel_x = 0.995
        panel_y = 0.97
        panel_font = 14
        title_font = 16

        #----------LEFT SIDE----------------
        #Subplot 0,0: Rain Rate
        axx[0,0].hlines(y=0,xmin=start-1,xmax=stop+1,color='k',linewidth = 0.5,zorder=0)
        axx[0,0].grid(alpha=0.2)
        axx[0,0].errorbar(x=resample_coords, y=means.δPrecip, yerr=stds.δPrecip, fmt='-k', ecolor='k', linewidth=0.5, capsize=10*spacing, zorder=100)
        sc = axx[0,0].scatter(x=resample_coords, y=means.δPrecip, s=sizes, c='C0',zorder=200)
        axx[0,0].set_xlim([resample_coords[0],resample_coords[-1]])
        axx[0,0].set_ylim([-12,31])
        axx[0,0].set_xticks(ticks)
        axx[0,0].set_xticklabels(ticklabels, fontsize=panel_font)
        #axx[0,0].set_xlabel(r'Normalized Time $ (t-t_{onset})/t_{δSST_{max}}$', fontsize=panel_font)
        axx[0,0].set_ylabel('Rain Rate $(mm/hr)$', fontsize=panel_font)
        axx[0,0].set_title('Rain Rate',fontsize=title_font)
        axx[0,0].legend(*sc.legend_elements("sizes", num=6),title='# of Events',loc='upper left', framealpha=0.9).set_zorder(1000)
        axx[0,0].text(panel_x, panel_y, '(a)', fontsize=panel_font, horizontalalignment='right', verticalalignment='top', transform=axx[0,0].transAxes, zorder=300)
        axx[0,0].yaxis.set_tick_params(labelsize=panel_font)
        

        
        #Subplot 1,0: δSST
        axx[1,0].hlines(y=0,xmin=start-1,xmax=stop+1,color='k',linewidth = 0.5,zorder=0)
        axx[1,0].grid(alpha=0.2)
        axx[1,0].errorbar(x=resample_coords, y=means.δsst, yerr=stds.δsst, fmt='-k', ecolor='k', linewidth=0.5, capsize=10*spacing, zorder=100)
        h1 = axx[1,0].scatter(x=resample_coords, y=means.δsst, s=sizes, c='red',zorder=200)
        h2 = axx[1,0].scatter(x=resample_coords, y=means.δSST, s=sizes*0.8, marker='v', facecolors='none', edgecolors='darkred' ,zorder=150) #plots COARE skin-corrected SST
        h3 = axx[1,0].scatter(x=resample_coords, y=means.δt_int_vF_SS, s=sizes, marker='o', facecolors='none', edgecolors='darkmagenta' ,zorder=150) #plots Bellenger model (Original)
        h4 = axx[1,0].scatter(x=resample_coords, y=means.δt_int_SS, s=sizes+2, marker='*', facecolors='none', edgecolors='navy' ,zorder=200) #plots Bellenger model (Modified)

        axx[1,0].set_xlim([resample_coords[0],resample_coords[-1]])
        axx[1,0].set_xticks(ticks)
        axx[1,0].set_xticklabels(ticklabels, fontsize=panel_font)
        #axx[1,0].set_xlabel(r'Normalized Time $ (t-t_{onset})/t_{δSST_{max}}$', fontsize=panel_font)
        axx[1,0].set_ylabel('δSST ($^\circ C$)', fontsize=panel_font)
        axx[1,0].set_ylim(temp_ylims)
        axx[1,0].set_title('Skin SST',fontsize=title_font)
        axx[1,0].legend([h1,h2,h3,h4],['Observed by KT-15 (10μm)', 'COARE3.5 Model', 'Original Bellenger Model (SS)', 'Modified Bellenger Model (SS)'])
        #axx[1,0].legend([h1,h2],['Observed by KT-15 (10μm)', 'COARE'])
        axx[1,0].text(panel_x, panel_y, '(b)', fontsize=panel_font, horizontalalignment='right', verticalalignment='top', transform=axx[1,0].transAxes, zorder=300)
        axx[1,0].yaxis.set_tick_params(labelsize=panel_font)

        #Subplot 2,0: δT_bulk (sea snake)
        axx[2,0].hlines(y=0,xmin=start-1,xmax=stop+1,color='k',linewidth = 0.5,zorder=0)
        axx[2,0].grid(alpha=0.2)
        axx[2,0].errorbar(x=resample_coords, y=means.δTsea, yerr=stds.δTsea, fmt='-k', ecolor='k', linewidth=0.5, capsize=10*spacing, zorder=100)
        h3 = axx[2,0].scatter(x=resample_coords, y=means.δTsea, s=sizes, c='darkred',zorder=200)
        h4 = axx[2,0].scatter(x=resample_coords, y=means.δTseaTSG, s=sizes, c='k',marker='1',zorder=100)
        h5 = axx[2,0].scatter(x=resample_coords, y=means.δt_int_vF_SS, s=sizes, marker='o', facecolors='none', edgecolors='darkmagenta' ,zorder=200) #plots Bellenger model SST (Original)
        h6 = axx[2,0].scatter(x=resample_coords, y=means.δt_int_SS, s=sizes+5, marker='*', facecolors='none', edgecolors='navy' ,zorder=200) #plots Bellenger model SST
        axx[2,0].legend([h3,h4,h5,h6],['Sea Snake (5cm)', 'Ship TSG (5m)', 'Original Bellenger Model (SS)', 'Modified Bellenger Model (SS)'])

        axx[2,0].set_xlim([resample_coords[0],resample_coords[-1]])
        axx[2,0].set_xticks(ticks)
        axx[2,0].set_xticklabels(ticklabels, fontsize=panel_font)
        #axx[2,0].set_xlabel(r'Normalized Time $ (t-t_{onset})/t_{δSST_{max}}$', fontsize=panel_font)
        axx[2,0].set_ylabel('$\delta T_{bulk} (^\circ C)$', fontsize=panel_font)
        #axx[2,0].set_ylim(temp_ylims)
        axx[2,0].set_title('Bulk Sea Temperature',fontsize=title_font)
        axx[2,0].text(panel_x, panel_y, '(c)', fontsize=panel_font, horizontalalignment='right', verticalalignment='top', transform=axx[2,0].transAxes, zorder=300)
        axx[2,0].yaxis.set_tick_params(labelsize=panel_font)
        
         #Subplot 3,0: Sensible Fluxes
        axx[3,0].hlines(y=0,xmin=start-1,xmax=stop+1,color='k',linewidth = 0.5,zorder=0)
        axx[3,0].grid(alpha=0.2)
        axx[3,0].errorbar(x=resample_coords, y=means_DC.δHsC, yerr=stds.δHsC, fmt='-g', ecolor='g', linewidth=0.5, capsize=10*spacing, zorder=300)
        h3 = axx[3,0].scatter(x=resample_coords, y=means_DC.δHsC, s=sizes_DC, c='forestgreen',zorder=400)
        axx[3,0].errorbar(x=resample_coords, y=means_DC.δshf, yerr=stds.δshf, fmt='-m', ecolor='m', linewidth=0.5, capsize=10*spacing, zorder=100)
        h4 = axx[3,0].scatter(x=resample_coords, y=means_DC.δshf, s=sizes_DC, c='purple',zorder=200)
        axx[3,0].errorbar(x=resample_coords, y=means_DC.δrhf, yerr=stds.δrhf, fmt='-k', ecolor='k', linewidth=0.5, capsize=10*spacing, zorder=100)
        h5 = axx[3,0].scatter(x=resample_coords, y=means_DC.δrhf, s=sizes_DC, c='darkblue',zorder=200)
        axx[3,0].set_xlim([resample_coords[0],resample_coords[-1]])
        axx[3,0].set_xticks(ticks)
        axx[3,0].set_xticklabels(ticklabels, fontsize=panel_font)
        axx[3,0].legend([h3,h4,h5],['DC','COARE','Gosnell Rain'])
        #axx[3,0].legend([h4],['COARE'],loc='upper left')
        axx[3,0].set_xlabel(r'Normalized Time $ (t-t_{onset})/t_{δSST_{max}}$', fontsize=panel_font)
        axx[3,0].set_ylabel('$\delta$ Flux $(W/m^2)$', fontsize=panel_font)
        axx[3,0].set_title('Sensible Fluxes',fontsize=title_font)
        axx[3,0].text(panel_x, panel_y, '(d)', fontsize=panel_font, horizontalalignment='right', verticalalignment='top', transform=axx[3,0].transAxes, zorder=300)
        axx[3,0].yaxis.set_tick_params(labelsize=panel_font)
        
        
        #--------------RIGHT SIDE---------------
        
        #Subplot 0,1: Air Temp
        axx[0,1].hlines(y=0,xmin=start-1,xmax=stop+1,color='k',linewidth = 0.5,zorder=0)
        axx[0,1].grid(alpha=0.2)
        axx[0,1].errorbar(x=resample_coords, y=means.δT02, yerr=stds.δT02, fmt='-k', ecolor='k', linewidth=0.5, capsize=10*spacing, zorder=100)
        axx[0,1].scatter(x=resample_coords, y=means.δT02, s=sizes, c='tomato',zorder=200)
        axx[0,1].set_xlim([resample_coords[0],resample_coords[-1]])
        axx[0,1].set_ylim([-2.3,0.3])
        axx[0,1].set_xticks(ticks)
        axx[0,1].set_xticklabels(ticklabels, fontsize=panel_font)
        #axx[0,1].set_xlabel(r'Normalized Time $ (t-t_{onset})/t_{δSST_{max}}$', fontsize=panel_font)
        axx[0,1].set_ylabel('$\delta T_{air} (^\circ C)$', fontsize=panel_font)
        axx[0,1].set_title('Air Temperature',fontsize=title_font)
        axx[0,1].text(panel_x, panel_y, '(e)', fontsize=panel_font, horizontalalignment='right', verticalalignment='top', transform=axx[0,1].transAxes, zorder=300)
        
        #Subplot 1,1: Wind Speed
        axx[1,1].hlines(y=0,xmin=start-1,xmax=stop+1,color='k',linewidth = 0.5,zorder=0)
        axx[1,1].grid(alpha=0.2)
        axx[1,1].errorbar(x=resample_coords, y=means.δwspd, yerr=stds.δwspd, fmt='-k', ecolor='k', linewidth=0.5, capsize=10*spacing, zorder=100)
        axx[1,1].scatter(x=resample_coords, y=means.δwspd, s=sizes, c='darkslategray',zorder=200)
        axx[1,1].set_xlim([resample_coords[0],resample_coords[-1]])
        axx[1,1].set_xticks(ticks)
        axx[1,1].set_xticklabels(ticklabels, fontsize=panel_font)
        #axx[1,1].set_xlabel(r'Normalized Time $ (t-t_{onset})/t_{δSST_{max}}$', fontsize=panel_font)
        axx[1,1].set_ylabel('$\delta U (m/s)$', fontsize=panel_font)
        axx[1,1].set_title('Wind Speed',fontsize=title_font)
        axx[1,1].text(panel_x, panel_y, '(f)', fontsize=panel_font, horizontalalignment='right', verticalalignment='top', transform=axx[1,1].transAxes, zorder=300)
        axx[1,1].yaxis.set_tick_params(labelsize=panel_font)
        
        #Subplot 2,1: δQ
        axx[2,1].hlines(y=0,xmin=start-1,xmax=stop+1,color='k',linewidth = 0.5,zorder=0)
        axx[2,1].grid(alpha=0.2)
        axx[2,1].errorbar(x=resample_coords, y=means.δQ02, yerr=stds.δQ02, fmt='-k', ecolor='k', linewidth=0.5, capsize=10*spacing, zorder=100)
        axx[2,1].scatter(x=resample_coords, y=means.δQ02, s=sizes, c='darkslateblue',zorder=200)
        axx[2,1].set_xlim([resample_coords[0],resample_coords[-1]])
        axx[2,1].set_xticks(ticks)
        axx[2,1].set_xticklabels(ticklabels, fontsize=panel_font)
        #axx[2,1].set_xlabel(r'Normalized Time $ (t-t_{onset})/t_{δSST_{max}}$', fontsize=panel_font)
        axx[2,1].set_ylabel('$\delta Q_{2m} (g/kg)$', fontsize=panel_font)
        axx[2,1].set_title('Specific Humidity',fontsize=title_font)
        axx[2,1].text(panel_x, panel_y, '(g)', fontsize=panel_font, horizontalalignment='right', verticalalignment='top', transform=axx[2,1].transAxes, zorder=300)
        axx[2,1].yaxis.set_tick_params(labelsize=panel_font)

        
        #Subplot 3,1: Latent fluxes
        #axx[3,1].plot(resample_coords, np.zeros(len(resample_coords)),linewidth = 0.5, zorder=0)
        #axx[3,1].errorbar(x=resample_coords, y=means_DC.δHlC, yerr=stds_DC.δHlC, fmt='-g', ecolor='g', linewidth=0.5, capsize=10*spacing, zorder=100)
        #h1 = axx[3,1].scatter(x=resample_coords, y=means_DC.δHlC, s=sizes_DC, c='forestgreen',zorder=200)
        #axx[3,1].errorbar(x=resample_coords, y=means_DC.δlhf, yerr=stds.δlhf, fmt='-m', ecolor='m', linewidth=0.5, capsize=10*spacing, zorder=100)
        #h2 = axx[3,1].scatter(x=resample_coords, y=means_DC.δlhf, s=sizes_DC, c='purple',zorder=200)
        #axx[3,1].set_xlim([resample_coords[0],resample_coords[-1]])
        #axx[2,0].set_ylim([-80,300])
        #axx[3,1].set_xticks(ticks)
        #axx[3,1].set_xticklabels(ticklabels)
        #axx[3,1].legend([h1,h2],['DC','COARE'])
        #axx[3,1].set_xlabel(r'Normalized Time $ (t-t_{onset})/t_{δSST_{max}}$')
        #axx[3,1].set_ylabel('$\delta$ Flux $(W/m^2)$')
        #axx[3,1].set_title('Latent Fluxes')
        
        
        #Subplot 3,1: Bulk Salinity
        axx[3,1].hlines(y=0,xmin=start-1,xmax=stop+1,color='k',linewidth = 0.5,zorder=0)
        axx[3,1].grid(alpha=0.2)
        axx[3,1].errorbar(x=resample_coords, y=means.δSalTSG, yerr=stds.δSalTSG, fmt='-k', ecolor='k', linewidth=0.5, capsize=10*spacing, zorder=10)
        h1 = axx[3,1].scatter(x=resample_coords, y=means.δSalTSG, s=sizes, c='blue',zorder=20)
        axx[3,1].errorbar(x=resample_coords, y=means.δs_int, yerr=stds.δs_int, fmt='-k', ecolor='k', linewidth=0.5, capsize=10*spacing, zorder=10)
        h2 = axx[3,1].scatter(x=resample_coords, y=means.δs_int, s=sizes, marker='*', facecolor='None', c='k')
        axx[3,1].set_xlim([resample_coords[0],resample_coords[-1]])
        #axx[3,1].set_ylim([-0.18,0.28])
        axx[3,1].set_xticks(ticks)
        axx[3,1].set_xticklabels(ticklabels, fontsize=panel_font)
        axx[3,1].set_xlabel(r'Normalized Time $ (t-t_{onset})/t_{δSST_{max}}$', fontsize=panel_font)
        axx[3,1].set_ylabel('$\delta S_{bulk} (psu)$', fontsize=panel_font)
        axx[3,1].set_title('Salinity',fontsize=title_font)
        axx[3,1].legend([h1,h2],['Ship TSG (5m)','Bellenger Model (0m)'])
        axx[3,1].text(panel_x, panel_y, '(h)', fontsize=panel_font, horizontalalignment='right', verticalalignment='top', transform=axx[3,1].transAxes, zorder=300)
        axx[3,1].yaxis.set_tick_params(labelsize=panel_font)


        #Subplot 2,0: Sensible Heat Flux Due to Rain
        #axx[2,0].plot(resample_coords, np.zeros(len(resample_coords)),linewidth = 0.5, zorder=0)
        #axx[2,0].errorbar(x=resample_coords, y=means.δrhf, yerr=stds.δrhf, fmt='-k', ecolor='k', linewidth=0.5, capsize=10*spacing, zorder=100)
        #axx[2,0].scatter(x=resample_coords, y=means.δrhf, s=sizes, c='forestgreen',zorder=200)
        #axx[2,0].set_xlim([resample_coords[0],resample_coords[-1]])
        #axx[2,0].set_ylim([-80,300])
        #axx[2,0].set_xticks(ticks)
        #axx[2,0].set_xticklabels(ticklabels)
        #axx[2,0].set_xlabel(r'Normalized Time $ (t-t_{onset})/t_{δSST_{max}}$')
        #axx[2,0].set_ylabel('Rain Heat Flux $(W/m^2)$')
        #axx[2,0].set_title('Sensible Heat Flux Due to Rain (COARE)')
        

        plt.tight_layout()

    return composite_event




##################################################################################################

def extract_composite_event_revised(rain_event_list, sst_event_list, param_list, start, stop, spacing, plotflag, title):
    '''
    Revised extract_composite_event plotting the sea snake-initialized model results and changing some visualization.
    This function takes all of the rain events and creates a single composite event by normalizing the timescale.
    Timescale normalization is calculated as (current time - rain onset time)[minutes]/(minutes to max δSST).
    The function then plots the composite event for 4 variables (air temp, sst, bulk temp, bulk salinity).

    Inputs:
        rain_event_list  - list of xarray datasets for each rain event
        sst_event_list   - list of xarray datasets containing sst information for each rain event
        param_list       - list of data variable names for which there is a δ{param} to be found in rain_event_list (should be same as the param_list input to calculate_deltas)
        start            - desired start time in normalized coordinates (-2 recommended)
        stop             - desired end time in normalized coordinates (8-12 recommended)
        spacing          - step size ("bin size") for normalized time coordinates

    Outputs:
        composite_event  - xarray dataset with two dimensions - event number, and normalized timescale

    Also puts out a plot of the composite events in air temp, sst, bulk temp, and bulk salinity.
    '''

    #initialize array of resampling coordinates
    resample_coords = np.linspace(start, stop, num=int((stop-start)/spacing),endpoint=False)

    #create new lists of datasets to do time-normalization on
    sst_normed_list = copy.deepcopy(sst_event_list)
    rain_normed_list = copy.deepcopy(rain_event_list)

    #initialize master arrays (shaped by event # vs. normalized timescale)
    #sst
    normed_δsst = np.empty(shape=(len(rain_event_list),len(resample_coords)))
    #rain_event variables
    array_list = []
    for param in param_list:
        array_list.append(np.empty(shape=(len(rain_event_list),len(resample_coords))))

    #cycle through rain events
    for event_num in np.arange(0,len(rain_event_list)):
        #find times of rain onset and max sst response
        onset = rain_event_list[event_num].attrs['Rain Onset']
        t_max_δsst = (sst_event_list[event_num].attrs['Time of max δsst'] - onset).astype('timedelta64[m]').astype('float64')

        #normalized time = (current time - rain onset time)[minutes]/(minutes to max δSST)
        rain_normtime = (rain_event_list[event_num].index - onset).values.astype('timedelta64[s]').astype('float64')/60/t_max_δsst
        sst_normtime  = ((sst_event_list[event_num].index - onset).values.astype('timedelta64[s]').astype('float64')/60/t_max_δsst)

        #replace index with normalized time variable in the new list
        sst_normed_list[event_num]['index'] = sst_normtime
        rain_normed_list[event_num]['index'] = rain_normtime

        #interpolate data onto resampling coordinates
        sst_interp = sst_normed_list[event_num].interp(index=resample_coords)
        rain_interp = rain_normed_list[event_num].interp(index=resample_coords)

        #insert data into row of master array
        normed_δsst[event_num,:] = sst_interp.δsst.values
        for idx in np.arange(0,len(param_list)):
            param = param_list[idx]
            array_list[idx][event_num,:] = rain_interp[f'δ{param}'].values

    #create and fill new dataset containing normalized events
    composite_event = xr.Dataset(coords={'index':resample_coords,
                                         'event':np.arange(1,len(rain_event_list)+1)})
    #sst
    composite_event = composite_event.assign({'δsst':(['event','index'],normed_δsst)})
    #rain event variables
    for idx in np.arange(0,len(param_list)):
        param = param_list[idx]
        composite_event = composite_event.assign({f'δ{param}':(['event','index'],array_list[idx])})

    #count the number of data points at each timestep to size plots
    sizes = np.empty_like(resample_coords)
    sizes_DC = np.empty_like(resample_coords)
    for idx in np.arange(0,len(resample_coords)):
        #count how many non-nan entries there are at this timestep
        sizes[idx] = composite_event.δsst.sel(index=resample_coords[idx]).count().values.item()
        sizes_DC[idx] = composite_event.δHsC.sel(index=resample_coords[idx]).count().values.item()
    
    #takes means and stds across all events
    means = composite_event.mean(dim='event',skipna=True)
    stds = composite_event.std(dim='event',skipna=True)
    q10 = composite_event.quantile(0.1,dim='event')
    q90 = composite_event.quantile(0.9,dim='event')
    
    means_DC = composite_event.where(~np.isnan(composite_event.δHsC)).mean(dim='event',skipna=True)
    stds_DC = composite_event.where(~np.isnan(composite_event.δHsC)).std(dim='event',skipna=True)

    if plotflag == 1:
        #plot
        fig,axx= plt.subplots(nrows=4,ncols=2,figsize=(16,12),facecolor='w')
        ticks = np.arange(start,stop+1,1)
        ticklabels = ticks.astype('str').tolist()
        ticklabels[-start] = 'Onset'
        ticklabels[-start + 1] = r'$ \delta SST_{max} $'
        plt.suptitle(title, y=0.98, fontsize=18)

        temp_ylims = [-0.9,0.15]
        panel_x = 0.995
        panel_y = 0.97
        panel_font = 14
        title_font = 16

        #----------LEFT SIDE----------------
        #Subplot 0,0: Rain Rate
        axx[0,0].hlines(y=0,xmin=start-1,xmax=stop+1,color='k',linewidth = 0.5,zorder=10)
        axx[0,0].grid(alpha=0.8)
        axx[0,0].fill_between(x=resample_coords, y1=q10.δPrecip, y2=q90.δPrecip, color='lightgrey',zorder=1)
        axx[0,0].plot(resample_coords, means.δPrecip, color='C0',zorder=100)
        sc = axx[0,0].scatter(x=resample_coords, y=means.δPrecip, s=sizes, c='C0',zorder=200)
        
        axx[0,0].set_xlim([resample_coords[0],resample_coords[-1]])
        axx[0,0].set_ylim([-1,23])
        axx[0,0].set_xticks(ticks)
        axx[0,0].set_xticklabels(ticklabels, fontsize=panel_font)
        axx[0,0].set_ylabel('Rain Rate $(mm/hr)$', fontsize=panel_font)
        axx[0,0].set_title('Rain Rate',fontsize=title_font)
        axx[0,0].legend(*sc.legend_elements("sizes", num=5),title='# of Events',loc='upper left', framealpha=0.9).set_zorder(1000)
        axx[0,0].text(panel_x, panel_y, '(a)', fontsize=panel_font, horizontalalignment='right', verticalalignment='top', transform=axx[0,0].transAxes, zorder=300)
        axx[0,0].yaxis.set_tick_params(labelsize=panel_font)
        
        #Subplot 1,0: δSST
        axx[1,0].hlines(y=0,xmin=start-1,xmax=stop+1,color='k',linewidth = 0.5,zorder=10)
        axx[1,0].grid(alpha=0.8)
        axx[1,0].fill_between(x=resample_coords, y1=q10.δsst, y2=q90.δsst, color='lightgrey',zorder=1)
        axx[1,0].plot(resample_coords, means.δsst, color='red',zorder=100)
        h1 = axx[1,0].scatter(x=resample_coords, y=means.δsst, s=sizes, c='red',zorder=200)
        h2 = axx[1,0].scatter(x=resample_coords, y=means.δSST, s=sizes*0.8, marker='v', facecolors='none', edgecolors='darkred' ,zorder=150) #plots COARE skin-corrected SST
        h3 = axx[1,0].scatter(x=resample_coords, y=means.δt_int_vF_SS, s=sizes, marker='o', facecolors='none', edgecolors='darkmagenta' ,zorder=150) #plots Bellenger model (Original)
        h4 = axx[1,0].scatter(x=resample_coords, y=means.δt_int_SS, s=sizes+2, marker='*', facecolors='none', edgecolors='navy' ,zorder=200) #plots Bellenger model (Modified)

        axx[1,0].set_xlim([resample_coords[0],resample_coords[-1]])
        axx[1,0].set_xticks(ticks)
        axx[1,0].set_xticklabels(ticklabels, fontsize=panel_font)
        axx[1,0].set_ylabel('δSST ($^\circ C$)', fontsize=panel_font)
        axx[1,0].set_ylim(temp_ylims)
        axx[1,0].set_title('Skin SST',fontsize=title_font)
        axx[1,0].legend([h1,h2,h3,h4],['Observed by KT-15 (10μm)', 'COARE3.5 Model', 'Original Bellenger Model', 'Modified Bellenger Model'])
        axx[1,0].text(panel_x, panel_y, '(b)', fontsize=panel_font, horizontalalignment='right', verticalalignment='top', transform=axx[1,0].transAxes, zorder=300)
        axx[1,0].yaxis.set_tick_params(labelsize=panel_font)

        #Subplot 2,0: δT_bulk (sea snake)
        axx[2,0].hlines(y=0,xmin=start-1,xmax=stop+1,color='k',linewidth = 0.5,zorder=10)
        axx[2,0].grid(alpha=0.8)
        axx[2,0].fill_between(x=resample_coords, y1=q10.δTsea, y2=q90.δTsea, color='lightgrey', zorder=1)
        axx[2,0].plot(resample_coords, means.δTsea, color='darkred',zorder=100)
        h3 = axx[2,0].scatter(x=resample_coords, y=means.δTsea, s=sizes, c='darkred',zorder=200)
        h4 = axx[2,0].scatter(x=resample_coords, y=means.δTseaTSG, s=sizes, c='k',marker='^',zorder=100)
        axx[2,0].plot(resample_coords, means.δTseaTSG, color='k',linestyle='--', zorder=100)
        axx[2,0].legend([h3,h4],['Sea Snake (5cm)', 'Ship TSG (5m)'])

        axx[2,0].set_xlim([resample_coords[0],resample_coords[-1]])
        axx[2,0].set_xticks(ticks)
        axx[2,0].set_xticklabels(ticklabels, fontsize=panel_font)
        axx[2,0].set_ylabel('$\delta T_{bulk} (^\circ C)$', fontsize=panel_font)
        axx[2,0].set_title('Bulk Sea Temperature',fontsize=title_font)
        axx[2,0].text(panel_x, panel_y, '(c)', fontsize=panel_font, horizontalalignment='right', verticalalignment='top', transform=axx[2,0].transAxes, zorder=300)
        axx[2,0].yaxis.set_tick_params(labelsize=panel_font)
        
         #Subplot 3,0: Sensible Fluxes
        axx[3,0].hlines(y=0,xmin=start-1,xmax=stop+1,color='k',linewidth = 0.5,zorder=10)
        axx[3,0].grid(alpha=0.8)
        axx[3,0].fill_between(x=resample_coords, y1=q10.δrhf, y2=q90.δrhf, color='lightgrey',  zorder=1)
        h5 = axx[3,0].scatter(x=resample_coords, y=means.δrhf, s=sizes, c='darkblue',zorder=200)
        axx[3,0].plot(resample_coords, means.δrhf, color='darkblue',zorder=200)
        
        axx[3,0].fill_between(x=resample_coords, y1=q10.δshf, y2=q90.δshf, color='orchid',  alpha=0.3, zorder=8)
        h4 = axx[3,0].scatter(x=resample_coords, y=means.δshf, s=sizes, c='purple',zorder=200)
        axx[3,0].plot(resample_coords, means.δshf, color='purple',zorder=200)
        axx[3,0].set_xlim([resample_coords[0],resample_coords[-1]])
        axx[3,0].set_xticks(ticks)
        axx[3,0].set_xticklabels(ticklabels, fontsize=panel_font)
        axx[3,0].legend([h5,h4],['Gosnell Rain','COARE Turbulent'], loc='upper left')
        axx[3,0].set_xlabel(r'Normalized Time $ (t-t_{onset})/t_{δSST_{max}}$', fontsize=panel_font)
        axx[3,0].set_ylabel('$\delta$ Flux $(W/m^2)$', fontsize=panel_font)
        axx[3,0].set_title('Sensible Heat Flux',fontsize=title_font)
        axx[3,0].text(panel_x, panel_y, '(d)', fontsize=panel_font, horizontalalignment='right', verticalalignment='top', transform=axx[3,0].transAxes, zorder=300)
        axx[3,0].yaxis.set_tick_params(labelsize=panel_font)
        
        
        #--------------RIGHT SIDE---------------
        
        #Subplot 0,1: Air Temp
        axx[0,1].hlines(y=0,xmin=start-1,xmax=stop+1,color='k',linewidth = 0.5,zorder=10)
        axx[0,1].grid(alpha=0.8)
        axx[0,1].fill_between(x=resample_coords, y1=q10.δT02, y2=q90.δT02, color='lightgrey', zorder=0)
        axx[0,1].scatter(x=resample_coords, y=means.δT02, s=sizes, c='tomato',zorder=200)
        axx[0,1].plot(resample_coords, means.δT02, color='tomato',zorder=200)
        axx[0,1].set_xlim([resample_coords[0],resample_coords[-1]])
        axx[0,1].set_ylim([-2.5,0.5])
        axx[0,1].set_xticks(ticks)
        axx[0,1].set_xticklabels(ticklabels, fontsize=panel_font)
        axx[0,1].set_ylabel('$\delta T_{air} (^\circ C)$', fontsize=panel_font)
        axx[0,1].set_title('Air Temperature',fontsize=title_font)
        axx[0,1].text(panel_x, panel_y, '(e)', fontsize=panel_font, horizontalalignment='right', verticalalignment='top', transform=axx[0,1].transAxes, zorder=300)
        
        #Subplot 1,1: Wind Speed
        axx[1,1].hlines(y=0,xmin=start-1,xmax=stop+1,color='k',linewidth = 0.5,zorder=10)
        axx[1,1].grid(alpha=0.8)
        axx[1,1].fill_between(x=resample_coords, y1=q10.δwspd, y2=q90.δwspd, color='lightgrey', zorder=0)
        axx[1,1].scatter(x=resample_coords, y=means.δwspd, s=sizes, c='darkslategray',zorder=200)
        axx[1,1].plot(resample_coords, means.δwspd, color='darkslategray',zorder=200)
        axx[1,1].set_xlim([resample_coords[0],resample_coords[-1]])
        axx[1,1].set_xticks(ticks)
        axx[1,1].set_xticklabels(ticklabels, fontsize=panel_font)
        #axx[1,1].set_xlabel(r'Normalized Time $ (t-t_{onset})/t_{δSST_{max}}$', fontsize=panel_font)
        axx[1,1].set_ylabel('$\delta U (m/s)$', fontsize=panel_font)
        axx[1,1].set_title('Wind Speed',fontsize=title_font)
        axx[1,1].text(panel_x, panel_y, '(f)', fontsize=panel_font, horizontalalignment='right', verticalalignment='top', transform=axx[1,1].transAxes, zorder=300)
        axx[1,1].yaxis.set_tick_params(labelsize=panel_font)
        
        #Subplot 2,1: δQ
        axx[2,1].hlines(y=0,xmin=start-1,xmax=stop+1,color='k',linewidth = 0.5,zorder=10)
        axx[2,1].grid(alpha=0.8)
        axx[2,1].fill_between(x=resample_coords, y1=q10.δQ02, y2=q90.δQ02, color='lightgrey', zorder=0)
        axx[2,1].scatter(x=resample_coords, y=means.δQ02, s=sizes, c='darkslateblue',zorder=200)
        axx[2,1].plot(resample_coords, means.δQ02, color='darkslateblue',zorder=200)
        axx[2,1].set_xlim([resample_coords[0],resample_coords[-1]])
        axx[2,1].set_xticks(ticks)
        axx[2,1].set_xticklabels(ticklabels, fontsize=panel_font)
        #axx[2,1].set_xlabel(r'Normalized Time $ (t-t_{onset})/t_{δSST_{max}}$', fontsize=panel_font)
        axx[2,1].set_ylabel('$\delta Q_{2m} (g/kg)$', fontsize=panel_font)
        axx[2,1].set_title('Specific Humidity',fontsize=title_font)
        axx[2,1].text(panel_x, panel_y, '(g)', fontsize=panel_font, horizontalalignment='right', verticalalignment='top', transform=axx[2,1].transAxes, zorder=300)
        axx[2,1].yaxis.set_tick_params(labelsize=panel_font)
        
        #Subplot 3,1: Bulk Salinity
        axx[3,1].hlines(y=0,xmin=start-1,xmax=stop+1,color='k',linewidth = 0.5,zorder=10)
        axx[3,1].grid(alpha=0.8)
        axx[3,1].fill_between(x=resample_coords, y1=q10.δSalTSG, y2=q90.δSalTSG, color='lightgrey', zorder=0)
        h1 = axx[3,1].scatter(x=resample_coords, y=means.δSalTSG, s=sizes, c='blue',zorder=20)
        axx[3,1].plot(resample_coords, means.δSalTSG, color='blue',zorder=20)
        axx[3,1].set_xlim([resample_coords[0],resample_coords[-1]])
        #axx[3,1].set_ylim([-0.18,0.28])
        axx[3,1].set_xticks(ticks)
        axx[3,1].set_xticklabels(ticklabels, fontsize=panel_font)
        axx[3,1].set_xlabel(r'Normalized Time $ (t-t_{onset})/t_{δSST_{max}}$', fontsize=panel_font)
        axx[3,1].set_ylabel('$\delta S_{bulk} (psu)$', fontsize=panel_font)
        axx[3,1].set_title('Bulk Salinity (Ship TSG)',fontsize=title_font)
        axx[3,1].text(panel_x, panel_y, '(h)', fontsize=panel_font, horizontalalignment='right', verticalalignment='top', transform=axx[3,1].transAxes, zorder=300)
        axx[3,1].yaxis.set_tick_params(labelsize=panel_font)


        plt.tight_layout()

    return composite_event




##################################################################################################

def extract_composite_event_10panel(rain_event_list, sst_event_list, param_list, start, stop, spacing, plotflag, title):
    '''
    Revised extract_composite_event plotting the sea snake-initialized model results and changing some visualization.
    This function takes all of the rain events and creates a single composite event by normalizing the timescale.
    Timescale normalization is calculated as (current time - rain onset time)[minutes]/(minutes to max δSST).
    The function then plots the composite event for 4 variables (air temp, sst, bulk temp, bulk salinity).

    Inputs:
        rain_event_list  - list of xarray datasets for each rain event
        sst_event_list   - list of xarray datasets containing sst information for each rain event
        param_list       - list of data variable names for which there is a δ{param} to be found in rain_event_list (should be same as the param_list input to calculate_deltas)
        start            - desired start time in normalized coordinates (-2 recommended)
        stop             - desired end time in normalized coordinates (8-12 recommended)
        spacing          - step size ("bin size") for normalized time coordinates

    Outputs:
        composite_event  - xarray dataset with two dimensions - event number, and normalized timescale

    Also puts out a plot of the composite events in air temp, sst, bulk temp, and bulk salinity.
    '''

    #initialize array of resampling coordinates
    resample_coords = np.linspace(start, stop, num=int((stop-start)/spacing),endpoint=False)

    #create new lists of datasets to do time-normalization on
    sst_normed_list = copy.deepcopy(sst_event_list)
    rain_normed_list = copy.deepcopy(rain_event_list)

    #initialize master arrays (shaped by event # vs. normalized timescale)
    #sst
    normed_δsst = np.empty(shape=(len(rain_event_list),len(resample_coords)))
    #rain_event variables
    array_list = []
    for param in param_list:
        array_list.append(np.empty(shape=(len(rain_event_list),len(resample_coords))))

    #cycle through rain events
    for event_num in np.arange(0,len(rain_event_list)):
        #find times of rain onset and max sst response
        onset = rain_event_list[event_num].attrs['Rain Onset']
        t_max_δsst = (sst_event_list[event_num].attrs['Time of max δsst'] - onset).astype('timedelta64[m]').astype('float64')

        #normalized time = (current time - rain onset time)[minutes]/(minutes to max δSST)
        rain_normtime = (rain_event_list[event_num].index - onset).values.astype('timedelta64[s]').astype('float64')/60/t_max_δsst
        sst_normtime  = ((sst_event_list[event_num].index - onset).values.astype('timedelta64[s]').astype('float64')/60/t_max_δsst)

        #replace index with normalized time variable in the new list
        sst_normed_list[event_num]['index'] = sst_normtime
        rain_normed_list[event_num]['index'] = rain_normtime

        #interpolate data onto resampling coordinates
        sst_interp = sst_normed_list[event_num].interp(index=resample_coords)
        rain_interp = rain_normed_list[event_num].interp(index=resample_coords)

        #insert data into row of master array
        normed_δsst[event_num,:] = sst_interp.δsst.values
        for idx in np.arange(0,len(param_list)):
            param = param_list[idx]
            array_list[idx][event_num,:] = rain_interp[f'δ{param}'].values

    #create and fill new dataset containing normalized events
    composite_event = xr.Dataset(coords={'index':resample_coords,
                                         'event':np.arange(1,len(rain_event_list)+1)})
    #sst
    composite_event = composite_event.assign({'δsst':(['event','index'],normed_δsst)})
    #rain event variables
    for idx in np.arange(0,len(param_list)):
        param = param_list[idx]
        composite_event = composite_event.assign({f'δ{param}':(['event','index'],array_list[idx])})

    #count the number of data points at each timestep to size plots
    sizes = np.empty_like(resample_coords)
    sizes_DC = np.empty_like(resample_coords)
    for idx in np.arange(0,len(resample_coords)):
        #count how many non-nan entries there are at this timestep
        sizes[idx] = composite_event.δsst.sel(index=resample_coords[idx]).count().values.item()
        sizes_DC[idx] = composite_event.δHsC.sel(index=resample_coords[idx]).count().values.item()
    
    #takes means and stds across all events
    means = composite_event.mean(dim='event',skipna=True)
    stds = composite_event.std(dim='event',skipna=True)
    q10 = composite_event.quantile(0.1,dim='event')
    q90 = composite_event.quantile(0.9,dim='event')
    
    means_DC = composite_event.where(~np.isnan(composite_event.δHsC)).mean(dim='event',skipna=True)
    stds_DC = composite_event.where(~np.isnan(composite_event.δHsC)).std(dim='event',skipna=True)

    if plotflag == 1:
        #plot
        fig,axx= plt.subplots(nrows=5,ncols=2,figsize=(17,17),facecolor='w')
        ticks = np.arange(start,stop+1,1)
        ticklabels = ticks.astype('str').tolist()
        ticklabels[-start] = 'Onset'
        ticklabels[-start + 1] = r'$ \delta SST_{max} $'
        plt.suptitle(title, fontsize=24, y=0.92)

        temp_ylims = [-0.9,0.15]
        panel_x = 0.995
        panel_y = 0.97
        panel_font = 16
        title_font = 18

        #----------LEFT SIDE----------------
        #Subplot 0,0: Rain Rate
        axx[0,0].hlines(y=0,xmin=start-1,xmax=stop+1,color='k',linewidth = 0.5,zorder=10)
        axx[0,0].grid(alpha=0.8)
        axx[0,0].fill_between(x=resample_coords, y1=q10.δPrecip, y2=q90.δPrecip, color='lightgrey',zorder=1)
        axx[0,0].plot(resample_coords, means.δPrecip, color='C0',zorder=100)
        sc = axx[0,0].scatter(x=resample_coords, y=means.δPrecip, s=sizes, c='C0',zorder=200)
        
        axx[0,0].set_xlim([resample_coords[0],resample_coords[-1]])
        axx[0,0].set_ylim([-1,23])
        axx[0,0].set_xticks(ticks)
        axx[0,0].set_xticklabels(ticklabels, fontsize=panel_font)
        axx[0,0].set_ylabel('Rain Rate $(mm/hr)$', fontsize=panel_font)
        axx[0,0].set_title('Rain Rate',fontsize=title_font)
        axx[0,0].legend(*sc.legend_elements("sizes", num=5),title='# of Events',loc='upper left', framealpha=0.9).set_zorder(1000)
        axx[0,0].text(panel_x, panel_y, '(a)', fontsize=panel_font, horizontalalignment='right', verticalalignment='top', transform=axx[0,0].transAxes, zorder=300)
        axx[0,0].yaxis.set_tick_params(labelsize=panel_font)
        
        #Subplot 1,0: δSST
        axx[1,0].hlines(y=0,xmin=start-1,xmax=stop+1,color='k',linewidth = 0.5,zorder=10)
        axx[1,0].grid(alpha=0.8)
        axx[1,0].fill_between(x=resample_coords, y1=q10.δsst, y2=q90.δsst, color='lightgrey',zorder=1)
        axx[1,0].plot(resample_coords, means.δsst, color='red',zorder=100)
        h1 = axx[1,0].scatter(x=resample_coords, y=means.δsst, s=sizes, c='red',zorder=200)
        h2 = axx[1,0].scatter(x=resample_coords, y=means.δSST, s=sizes*0.8, marker='v', facecolors='none', edgecolors='darkred' ,zorder=150) #plots COARE skin-corrected SST
        h3 = axx[1,0].scatter(x=resample_coords, y=means.δt_int_vF_SS, s=sizes, marker='o', facecolors='none', edgecolors='darkmagenta' ,zorder=150) #plots Bellenger model (Original)
        h4 = axx[1,0].scatter(x=resample_coords, y=means.δt_int_SS, s=sizes+2, marker='*', facecolors='none', edgecolors='navy' ,zorder=200) #plots Bellenger model (Modified)

        axx[1,0].set_xlim([resample_coords[0],resample_coords[-1]])
        axx[1,0].set_xticks(ticks)
        axx[1,0].set_xticklabels(ticklabels, fontsize=panel_font)
        axx[1,0].set_ylabel('$δT_{skin}$ ($^\circ C$)', fontsize=panel_font)
        axx[1,0].set_ylim(temp_ylims)
        axx[1,0].set_title('Skin Temperature',fontsize=title_font)
        axx[1,0].legend([h1,h2,h3,h4],['Observed by KT-15 (10μm)', 'COARE3.5 Model', 'Original Bellenger Model', 'Modified Bellenger Model'])
        axx[1,0].text(panel_x, panel_y, '(b)', fontsize=panel_font, horizontalalignment='right', verticalalignment='top', transform=axx[1,0].transAxes, zorder=300)
        axx[1,0].yaxis.set_tick_params(labelsize=panel_font)

        #Subplot 2,0: δT_bulk (sea snake)
        axx[2,0].hlines(y=0,xmin=start-1,xmax=stop+1,color='k',linewidth = 0.5,zorder=10)
        axx[2,0].grid(alpha=0.8)
        axx[2,0].fill_between(x=resample_coords, y1=q10.δTsea, y2=q90.δTsea, color='lightgrey', zorder=1)
        axx[2,0].plot(resample_coords, means.δTsea, color='darkred',zorder=100)
        h3 = axx[2,0].scatter(x=resample_coords, y=means.δTsea, s=sizes, c='darkred',zorder=200)
        h4 = axx[2,0].scatter(x=resample_coords, y=means.δTseaTSG, s=sizes, c='k',marker='^',zorder=100)
        axx[2,0].plot(resample_coords, means.δTseaTSG, color='k',linestyle='--', zorder=100)
        axx[2,0].legend([h3,h4],['Sea Snake (5cm)', 'Ship TSG (5m)'])

        axx[2,0].set_xlim([resample_coords[0],resample_coords[-1]])
        axx[2,0].set_xticks(ticks)
        axx[2,0].set_xticklabels(ticklabels, fontsize=panel_font)
        axx[2,0].set_ylabel('$\delta T_{bulk}$ $(^\circ C)$', fontsize=panel_font)
        axx[2,0].set_title('Bulk Sea Temperature',fontsize=title_font)
        axx[2,0].text(panel_x, panel_y, '(c)', fontsize=panel_font, horizontalalignment='right', verticalalignment='top', transform=axx[2,0].transAxes, zorder=300)
        axx[2,0].yaxis.set_tick_params(labelsize=panel_font)
        

        #Subplot 3,0: Sensible Fluxes
        axx[3,0].hlines(y=0,xmin=start-1,xmax=stop+1,color='k',linewidth = 0.5,zorder=10)
        axx[3,0].grid(alpha=0.8)
        axx[3,0].fill_between(x=resample_coords, y1=q10.δrhf, y2=q90.δrhf, color='lightgrey',  zorder=1)
        h5 = axx[3,0].scatter(x=resample_coords, y=means.δrhf, s=sizes, c='darkblue',zorder=200)
        axx[3,0].plot(resample_coords, means.δrhf, color='darkblue',zorder=200)
        
        axx[3,0].fill_between(x=resample_coords, y1=q10.δshf, y2=q90.δshf, color='orchid',  alpha=0.3, zorder=8)
        h4 = axx[3,0].scatter(x=resample_coords, y=means.δshf, s=sizes, facecolors='none', edgecolors='purple',zorder=200)
        axx[3,0].plot(resample_coords, means.δshf, color='purple',zorder=200)
        axx[3,0].set_xlim([resample_coords[0],resample_coords[-1]])
        axx[3,0].set_xticks(ticks)
        axx[3,0].set_xticklabels(ticklabels, fontsize=panel_font)
        axx[3,0].legend([h5,h4],['Gosnell Rain','COARE Turbulent'], loc='upper left')
        axx[3,0].set_ylabel('$\delta$ Flux $(W/m^2)$', fontsize=panel_font)
        axx[3,0].set_title('Sensible Heat Flux',fontsize=title_font)
        axx[3,0].text(panel_x, panel_y, '(d)', fontsize=panel_font, horizontalalignment='right', verticalalignment='top', transform=axx[3,0].transAxes, zorder=300)
        axx[3,0].yaxis.set_tick_params(labelsize=panel_font)
        
        
        #Subplot 4,0: Tskin - Tair
        axx[4,0].hlines(y=0,xmin=start-1,xmax=stop+1,color='k',linewidth = 0.5,zorder=10)
        axx[4,0].grid(alpha=0.8)
        axx[4,0].fill_between(x=resample_coords, y1=q10.δΔT_skin_air, y2=q90.δΔT_skin_air, color='lightgrey', zorder=0)
        axx[4,0].scatter(x=resample_coords, y=means.δΔT_skin_air, s=sizes, c='C3',zorder=200)
        axx[4,0].plot(resample_coords, means.δΔT_skin_air, color='C3',zorder=200)
        axx[4,0].set_xlim([resample_coords[0],resample_coords[-1]])
        axx[4,0].set_xticks(ticks)
        axx[4,0].set_xticklabels(ticklabels, fontsize=panel_font)
        #axx[4,0].set_xlabel(r'Normalized Time $ (t-t_{onset})/t_{δSST_{max}}$', fontsize=panel_font)
        axx[4,0].set_ylabel('$\delta$ $\Delta T$ $(^\circ C)$', fontsize=panel_font)
        axx[4,0].set_title('$T_{skin} - T_{air}$',fontsize=title_font)
        axx[4,0].text(panel_x, panel_y, '(e)', fontsize=panel_font, horizontalalignment='right', verticalalignment='top', transform=axx[4,0].transAxes, zorder=300)
        axx[4,0].yaxis.set_tick_params(labelsize=panel_font)
        axx[4,0].set_xlabel(r'Normalized Time $ (t-t_{onset})/t_{δSST_{max}}$', fontsize=panel_font)


        
        #--------------RIGHT SIDE---------------
        #Subplot 0,1: Wind Speed
        axx[0,1].hlines(y=0,xmin=start-1,xmax=stop+1,color='k',linewidth = 0.5,zorder=10)
        axx[0,1].grid(alpha=0.8)
        axx[0,1].fill_between(x=resample_coords, y1=q10.δwspd, y2=q90.δwspd, color='lightgrey', zorder=0)
        axx[0,1].scatter(x=resample_coords, y=means.δwspd, s=sizes, c='darkslategray',zorder=200)
        axx[0,1].plot(resample_coords, means.δwspd, color='darkslategray',zorder=200)
        axx[0,1].set_xlim([resample_coords[0],resample_coords[-1]])
        axx[0,1].set_xticks(ticks)
        axx[0,1].set_xticklabels(ticklabels, fontsize=panel_font)
        #axx[0,1].set_xlabel(r'Normalized Time $ (t-t_{onset})/t_{δSST_{max}}$', fontsize=panel_font)
        axx[0,1].set_ylabel('$\delta U$ $(m/s)$', fontsize=panel_font)
        axx[0,1].set_title('Wind Speed',fontsize=title_font)
        axx[0,1].text(panel_x, panel_y, '(f)', fontsize=panel_font, horizontalalignment='right', verticalalignment='top', transform=axx[0,1].transAxes, zorder=300)
        axx[0,1].yaxis.set_tick_params(labelsize=panel_font)
        
        #Subplot 1,1: Air Temp
        axx[1,1].hlines(y=0,xmin=start-1,xmax=stop+1,color='k',linewidth = 0.5,zorder=10)
        axx[1,1].grid(alpha=0.8)
        axx[1,1].fill_between(x=resample_coords, y1=q10.δT02, y2=q90.δT02, color='lightgrey', zorder=0)
        axx[1,1].scatter(x=resample_coords, y=means.δT02, s=sizes, c='tomato',zorder=200)
        axx[1,1].plot(resample_coords, means.δT02, color='tomato',zorder=200)
        axx[1,1].set_xlim([resample_coords[0],resample_coords[-1]])
        axx[1,1].set_ylim([-2.5,0.5])
        axx[1,1].set_xticks(ticks)
        axx[1,1].set_xticklabels(ticklabels, fontsize=panel_font)
        axx[1,1].set_ylabel('$\delta T_{air}$ $(^\circ C)$', fontsize=panel_font)
        axx[1,1].set_title('Air Temperature',fontsize=title_font)
        axx[1,1].yaxis.set_tick_params(labelsize=panel_font)
        axx[1,1].text(panel_x, panel_y, '(g)', fontsize=panel_font, horizontalalignment='right', verticalalignment='top', transform=axx[1,1].transAxes, zorder=300)
        
        #Subplot 2,1: Bulk Salinity
        axx[2,1].hlines(y=0,xmin=start-1,xmax=stop+1,color='k',linewidth = 0.5,zorder=10)
        axx[2,1].grid(alpha=0.8)
        axx[2,1].fill_between(x=resample_coords, y1=q10.δSalTSG, y2=q90.δSalTSG, color='lightgrey', zorder=0)
        h1 = axx[2,1].scatter(x=resample_coords, y=means.δSalTSG, s=sizes, c='blue',zorder=20)
        axx[2,1].plot(resample_coords, means.δSalTSG, color='blue',zorder=20)
        axx[2,1].set_xlim([resample_coords[0],resample_coords[-1]])
        #axx[2,1].set_ylim([-0.18,0.28])
        axx[2,1].set_xticks(ticks)
        axx[2,1].set_xticklabels(ticklabels, fontsize=panel_font)
        axx[2,1].set_ylabel('$\delta S_{5m}$ $(psu)$', fontsize=panel_font)
        axx[2,1].set_title('Bulk Salinity (Ship TSG)',fontsize=title_font)
        axx[2,1].text(panel_x, panel_y, '(h)', fontsize=panel_font, horizontalalignment='right', verticalalignment='top', transform=axx[2,1].transAxes, zorder=300)
        axx[2,1].yaxis.set_tick_params(labelsize=panel_font)
        
        #Subplot 3,1: Specific Humidity
        axx[3,1].hlines(y=0,xmin=start-1,xmax=stop+1,color='k',linewidth = 0.5,zorder=10)
        axx[3,1].grid(alpha=0.8)
        axx[3,1].fill_between(x=resample_coords, y1=q10.δQ02, y2=q90.δQ02, color='lightgrey', zorder=0)
        axx[3,1].scatter(x=resample_coords, y=means.δQ02, s=sizes, c='darkslateblue',zorder=200)
        axx[3,1].plot(resample_coords, means.δQ02, color='darkslateblue',zorder=200)
        axx[3,1].set_xlim([resample_coords[0],resample_coords[-1]])
        axx[3,1].set_xticks(ticks)
        axx[3,1].set_xticklabels(ticklabels, fontsize=panel_font)
        #axx[3,1].set_xlabel(r'Normalized Time $ (t-t_{onset})/t_{δSST_{max}}$', fontsize=panel_font)
        axx[3,1].set_ylabel('$\delta q_{air}$ $(g/kg)$', fontsize=panel_font)
        axx[3,1].set_title('Specific Humidity',fontsize=title_font)
        axx[3,1].text(panel_x, panel_y, '(i)', fontsize=panel_font, horizontalalignment='right', verticalalignment='top', transform=axx[3,1].transAxes, zorder=300)
        axx[3,1].yaxis.set_tick_params(labelsize=panel_font)
        
        #Subplot 4,1: Qsat(Tskin) - Qair
        axx[4,1].hlines(y=0,xmin=start-1,xmax=stop+1,color='k',linewidth = 0.5,zorder=10)
        axx[4,1].grid(alpha=0.8)
        axx[4,1].fill_between(x=resample_coords, y1=q10.δΔq_skin_air, y2=q90.δΔq_skin_air, color='lightgrey', zorder=0)
        axx[4,1].scatter(x=resample_coords, y=means.δΔq_skin_air, s=sizes, c='k',zorder=200)
        axx[4,1].plot(resample_coords, means.δΔq_skin_air, color='k',zorder=200)
        axx[4,1].set_xlim([resample_coords[0],resample_coords[-1]])
        axx[4,1].set_xticks(ticks)
        axx[4,1].set_xticklabels(ticklabels, fontsize=panel_font)
        axx[4,1].set_xlabel(r'Normalized Time $ (t-t_{onset})/t_{δSST_{max}}$', fontsize=panel_font)
        axx[4,1].set_ylabel('$\delta$ $\Delta q$ $(g/kg)$', fontsize=panel_font)
        axx[4,1].set_title('$q_{sat}(T_{skin}) - q_{air}$',fontsize=title_font)
        axx[4,1].text(panel_x, panel_y, '(j)', fontsize=panel_font, horizontalalignment='right', verticalalignment='top', transform=axx[4,1].transAxes, zorder=300)
        axx[4,1].yaxis.set_tick_params(labelsize=panel_font)
        

        #plt.tight_layout()
        
        plt.subplots_adjust(hspace=0.4)

    return composite_event




##################################################################################################

def extract_composite_event_OLD(rain_event_list, sst_event_list, param_list, start, stop, spacing, plotflag, title):
    '''
    OLD COPY JUST IN CASE
    This function takes all of the rain events and creates a single composite event by normalizing the timescale.
    Timescale normalization is calculated as (current time - rain onset time)[minutes]/(minutes to max δSST).
    The function then plots the composite event for 4 variables (air temp, sst, bulk temp, bulk salinity).

    Inputs:
        rain_event_list  - list of xarray datasets for each rain event
        sst_event_list   - list of xarray datasets containing sst information for each rain event
        param_list       - list of data variable names for which there is a δ{param} to be found in rain_event_list (should be same as the param_list input to calculate_deltas)
        start            - desired start time in normalized coordinates (-2 recommended)
        stop             - desired end time in normalized coordinates (8-12 recommended)
        spacing          - step size ("bin size") for normalized time coordinates

    Outputs:
        composite_event  - xarray dataset with two dimensions - event number, and normalized timescale

    Also puts out a plot of the composite events in air temp, sst, bulk temp, and bulk salinity.
    '''

    #initialize array of resampling coordinates
    resample_coords = np.linspace(start, stop, num=int((stop-start)/spacing),endpoint=False)

    #create new lists of datasets to do time-normalization on
    sst_normed_list = copy.deepcopy(sst_event_list)
    rain_normed_list = copy.deepcopy(rain_event_list)

    #initialize master arrays (shaped by event # vs. normalized timescale)
    #sst
    normed_δsst = np.empty(shape=(len(rain_event_list),len(resample_coords)))
    #rain_event variables
    array_list = []
    for param in param_list:
        array_list.append(np.empty(shape=(len(rain_event_list),len(resample_coords))))

    #cycle through rain events
    for event_num in np.arange(0,len(rain_event_list)):
        #find times of rain onset and max sst response
        onset = rain_event_list[event_num].attrs['Rain Onset']
        t_max_δsst = (sst_event_list[event_num].attrs['Time of max δsst'] - onset).astype('timedelta64[m]').astype('float64')

        #normalized time = (current time - rain onset time)[minutes]/(minutes to max δSST)
        rain_normtime = (rain_event_list[event_num].index - onset).values.astype('timedelta64[s]').astype('float64')/60/t_max_δsst
        sst_normtime  = ((sst_event_list[event_num].index - onset).values.astype('timedelta64[s]').astype('float64')/60/t_max_δsst)

        #replace index with normalized time variable in the new list
        sst_normed_list[event_num]['index'] = sst_normtime
        rain_normed_list[event_num]['index'] = rain_normtime

        #interpolate data onto resampling coordinates
        sst_interp = sst_normed_list[event_num].interp(index=resample_coords)
        rain_interp = rain_normed_list[event_num].interp(index=resample_coords)

        #insert data into row of master array
        normed_δsst[event_num,:] = sst_interp.δsst.values
        for idx in np.arange(0,len(param_list)):
            param = param_list[idx]
            array_list[idx][event_num,:] = rain_interp[f'δ{param}'].values

    #create and fill new dataset containing normalized events
    composite_event = xr.Dataset(coords={'index':resample_coords,
                                         'event':np.arange(1,len(rain_event_list)+1)})
    #sst
    composite_event = composite_event.assign({'δsst':(['event','index'],normed_δsst)})
    #rain event variables
    for idx in np.arange(0,len(param_list)):
        param = param_list[idx]
        composite_event = composite_event.assign({f'δ{param}':(['event','index'],array_list[idx])})

    #count the number of data points at each timestep to size plots
    sizes = np.empty_like(resample_coords)
    sizes_DC = np.empty_like(resample_coords)
    for idx in np.arange(0,len(resample_coords)):
        #count how many non-nan entries there are at this timestep
        sizes[idx] = composite_event.δsst.sel(index=resample_coords[idx]).count().values.item()
        sizes_DC[idx] = composite_event.δHsC.sel(index=resample_coords[idx]).count().values.item()
        
    #takes means and stds across all events
    means = composite_event.mean(dim='event',skipna=True)
    stds = composite_event.std(dim='event',skipna=True)
    
    means_DC = composite_event.where(~np.isnan(composite_event.δHsC)).mean(dim='event',skipna=True)
    stds_DC = composite_event.where(~np.isnan(composite_event.δHsC)).mean(dim='event',skipna=True)

    if plotflag == 1:
        #plot
        fig,axx= plt.subplots(nrows=4,ncols=2,figsize=(16,12),facecolor='w')
        ticks = np.arange(start,stop+1,1)
        ticklabels = ticks.astype('str').tolist()
        ticklabels[-start] = 'Onset'
        ticklabels[-start + 1] = r'$ \delta SST_{max} $'
        plt.suptitle(title, y=0.98, fontsize=18)

        temp_ylims = [-1,0.2]
        panel_x = 0.995
        panel_y = 0.97
        panel_font = 14
        title_font = 16

        #----------LEFT SIDE----------------
        #Subplot 0,0: Rain Rate
        axx[0,0].hlines(y=0,xmin=start-1,xmax=stop+1,color='k',linewidth = 0.5,zorder=0)
        axx[0,0].grid(alpha=0.2)
        axx[0,0].errorbar(x=resample_coords, y=means.δPrecip, yerr=stds.δPrecip, fmt='-k', ecolor='k', linewidth=0.5, capsize=10*spacing, zorder=100)
        sc = axx[0,0].scatter(x=resample_coords, y=means.δPrecip, s=sizes, c='C0',zorder=200)
        axx[0,0].set_xlim([resample_coords[0],resample_coords[-1]])
        axx[0,0].set_ylim([-12,31])
        axx[0,0].set_xticks(ticks)
        axx[0,0].set_xticklabels(ticklabels, fontsize=panel_font)
        #axx[0,0].set_xlabel(r'Normalized Time $ (t-t_{onset})/t_{δSST_{max}}$', fontsize=panel_font)
        axx[0,0].set_ylabel('Rain Rate $(mm/hr)$', fontsize=panel_font)
        axx[0,0].set_title('Rain Rate',fontsize=title_font)
        axx[0,0].legend(*sc.legend_elements("sizes", num=6),title='# of Events',loc='upper left', framealpha=0.9).set_zorder(1000)
        axx[0,0].text(panel_x, panel_y, '(a)', fontsize=panel_font, horizontalalignment='right', verticalalignment='top', transform=axx[0,0].transAxes, zorder=300)
        axx[0,0].yaxis.set_tick_params(labelsize=panel_font)
        

        
        #Subplot 1,0: δSST
        axx[1,0].hlines(y=0,xmin=start-1,xmax=stop+1,color='k',linewidth = 0.5,zorder=0)
        axx[1,0].grid(alpha=0.2)
        axx[1,0].errorbar(x=resample_coords, y=means.δsst, yerr=stds.δsst, fmt='-k', ecolor='k', linewidth=0.5, capsize=10*spacing, zorder=100)
        h1 = axx[1,0].scatter(x=resample_coords, y=means.δsst, s=sizes, c='red',zorder=200)
        h2 = axx[1,0].scatter(x=resample_coords, y=means.δSST, s=sizes*0.8, marker='v', facecolors='none', edgecolors='darkred' ,zorder=150) #plots COARE skin-corrected SST
        h3 = axx[1,0].scatter(x=resample_coords, y=means.δSSTrain, s=sizes, facecolors='none', edgecolors='darkmagenta' ,zorder=150) #plots COARE skin-corrected SST + rain effect
        h4 = axx[1,0].scatter(x=resample_coords, y=means.δt_int, s=sizes, marker='*', facecolors='none', edgecolors='navy' ,zorder=200) #plots Bellenger model SST

        axx[1,0].set_xlim([resample_coords[0],resample_coords[-1]])
        axx[1,0].set_xticks(ticks)
        axx[1,0].set_xticklabels(ticklabels, fontsize=panel_font)
        #axx[1,0].set_xlabel(r'Normalized Time $ (t-t_{onset})/t_{δSST_{max}}$', fontsize=panel_font)
        axx[1,0].set_ylabel('δSST ($^\circ C$)', fontsize=panel_font)
        axx[1,0].set_ylim(temp_ylims)
        axx[1,0].set_title('Skin SST',fontsize=title_font)
        axx[1,0].legend([h1,h2,h3,h4],['Observed by KT-15 (10μm)', 'COARE (no rain)','COARE (with rain)', 'Bellenger Model'])
        #axx[1,0].legend([h1,h2],['Observed by KT-15 (10μm)', 'COARE'])
        axx[1,0].text(panel_x, panel_y, '(b)', fontsize=panel_font, horizontalalignment='right', verticalalignment='top', transform=axx[1,0].transAxes, zorder=300)
        axx[1,0].yaxis.set_tick_params(labelsize=panel_font)

        #Subplot 2,0: δT_bulk (sea snake)
        axx[2,0].hlines(y=0,xmin=start-1,xmax=stop+1,color='k',linewidth = 0.5,zorder=0)
        axx[2,0].grid(alpha=0.2)
        axx[2,0].errorbar(x=resample_coords, y=means.δTsea, yerr=stds.δTsea, fmt='-k', ecolor='k', linewidth=0.5, capsize=10*spacing, zorder=100)
        h3 = axx[2,0].scatter(x=resample_coords, y=means.δTsea, s=sizes, c='darkred',zorder=200)
        h4 = axx[2,0].scatter(x=resample_coords, y=means.δTseaTSG, s=sizes, c='k',marker='1',zorder=100)
        h5 = axx[2,0].scatter(x=resample_coords, y=means.δt_int, s=sizes, marker='*', facecolors='none', edgecolors='navy' ,zorder=200) #plots Bellenger model SST
        axx[2,0].legend([h3,h4,h5],['Sea Snake (5cm)', 'Ship TSG (5m)', 'Bellenger Model'])

        axx[2,0].set_xlim([resample_coords[0],resample_coords[-1]])
        axx[2,0].set_xticks(ticks)
        axx[2,0].set_xticklabels(ticklabels, fontsize=panel_font)
        #axx[2,0].set_xlabel(r'Normalized Time $ (t-t_{onset})/t_{δSST_{max}}$', fontsize=panel_font)
        axx[2,0].set_ylabel('$\delta T_{bulk} (^\circ C)$', fontsize=panel_font)
        #axx[2,0].set_ylim(temp_ylims)
        axx[2,0].set_title('Bulk Sea Temperature',fontsize=title_font)
        axx[2,0].text(panel_x, panel_y, '(c)', fontsize=panel_font, horizontalalignment='right', verticalalignment='top', transform=axx[2,0].transAxes, zorder=300)
        axx[2,0].yaxis.set_tick_params(labelsize=panel_font)
        
         #Subplot 3,0: Sensible Fluxes
        axx[3,0].hlines(y=0,xmin=start-1,xmax=stop+1,color='k',linewidth = 0.5,zorder=0)
        axx[3,0].grid(alpha=0.2)
        axx[3,0].errorbar(x=resample_coords, y=means_DC.δHsC, yerr=stds.δHsC, fmt='-g', ecolor='g', linewidth=0.5, capsize=10*spacing, zorder=300)
        h3 = axx[3,0].scatter(x=resample_coords, y=means_DC.δHsC, s=sizes_DC, c='forestgreen',zorder=400)
        axx[3,0].errorbar(x=resample_coords, y=means_DC.δshf, yerr=stds.δshf, fmt='-m', ecolor='m', linewidth=0.5, capsize=10*spacing, zorder=100)
        h4 = axx[3,0].scatter(x=resample_coords, y=means_DC.δshf, s=sizes_DC, c='purple',zorder=200)
        axx[3,0].errorbar(x=resample_coords, y=means_DC.δrhf, yerr=stds.δrhf, fmt='-k', ecolor='k', linewidth=0.5, capsize=10*spacing, zorder=100)
        h5 = axx[3,0].scatter(x=resample_coords, y=means_DC.δrhf, s=sizes_DC, c='darkblue',zorder=200)
        axx[3,0].set_xlim([resample_coords[0],resample_coords[-1]])
        axx[3,0].set_xticks(ticks)
        axx[3,0].set_xticklabels(ticklabels, fontsize=panel_font)
        axx[3,0].legend([h3,h4,h5],['DC','COARE','Gosnell Rain'])
        #axx[3,0].legend([h4],['COARE'],loc='upper left')
        axx[3,0].set_xlabel(r'Normalized Time $ (t-t_{onset})/t_{δSST_{max}}$', fontsize=panel_font)
        axx[3,0].set_ylabel('$\delta$ Flux $(W/m^2)$', fontsize=panel_font)
        axx[3,0].set_title('Sensible Fluxes',fontsize=title_font)
        axx[3,0].text(panel_x, panel_y, '(d)', fontsize=panel_font, horizontalalignment='right', verticalalignment='top', transform=axx[3,0].transAxes, zorder=300)
        axx[3,0].yaxis.set_tick_params(labelsize=panel_font)
        
        
        #--------------RIGHT SIDE---------------
        
        #Subplot 0,1: Air Temp
        axx[0,1].hlines(y=0,xmin=start-1,xmax=stop+1,color='k',linewidth = 0.5,zorder=0)
        axx[0,1].grid(alpha=0.2)
        axx[0,1].errorbar(x=resample_coords, y=means.δT02, yerr=stds.δT02, fmt='-k', ecolor='k', linewidth=0.5, capsize=10*spacing, zorder=100)
        axx[0,1].scatter(x=resample_coords, y=means.δT02, s=sizes, c='tomato',zorder=200)
        axx[0,1].set_xlim([resample_coords[0],resample_coords[-1]])
        axx[0,1].set_ylim([-2.3,0.3])
        axx[0,1].set_xticks(ticks)
        axx[0,1].set_xticklabels(ticklabels, fontsize=panel_font)
        #axx[0,1].set_xlabel(r'Normalized Time $ (t-t_{onset})/t_{δSST_{max}}$', fontsize=panel_font)
        axx[0,1].set_ylabel('$\delta T_{air} (^\circ C)$', fontsize=panel_font)
        axx[0,1].set_title('Air Temperature',fontsize=title_font)
        axx[0,1].text(panel_x, panel_y, '(e)', fontsize=panel_font, horizontalalignment='right', verticalalignment='top', transform=axx[0,1].transAxes, zorder=300)
        
        #Subplot 1,1: Wind Speed
        axx[1,1].hlines(y=0,xmin=start-1,xmax=stop+1,color='k',linewidth = 0.5,zorder=0)
        axx[1,1].grid(alpha=0.2)
        axx[1,1].errorbar(x=resample_coords, y=means.δwspd, yerr=stds.δwspd, fmt='-k', ecolor='k', linewidth=0.5, capsize=10*spacing, zorder=100)
        axx[1,1].scatter(x=resample_coords, y=means.δwspd, s=sizes, c='darkslategray',zorder=200)
        axx[1,1].set_xlim([resample_coords[0],resample_coords[-1]])
        axx[1,1].set_xticks(ticks)
        axx[1,1].set_xticklabels(ticklabels, fontsize=panel_font)
        #axx[1,1].set_xlabel(r'Normalized Time $ (t-t_{onset})/t_{δSST_{max}}$', fontsize=panel_font)
        axx[1,1].set_ylabel('$\delta U (m/s)$', fontsize=panel_font)
        axx[1,1].set_title('Wind Speed',fontsize=title_font)
        axx[1,1].text(panel_x, panel_y, '(f)', fontsize=panel_font, horizontalalignment='right', verticalalignment='top', transform=axx[1,1].transAxes, zorder=300)
        axx[1,1].yaxis.set_tick_params(labelsize=panel_font)
        
        #Subplot 2,1: δQ
        axx[2,1].hlines(y=0,xmin=start-1,xmax=stop+1,color='k',linewidth = 0.5,zorder=0)
        axx[2,1].grid(alpha=0.2)
        axx[2,1].errorbar(x=resample_coords, y=means.δQ02, yerr=stds.δQ02, fmt='-k', ecolor='k', linewidth=0.5, capsize=10*spacing, zorder=100)
        axx[2,1].scatter(x=resample_coords, y=means.δQ02, s=sizes, c='darkslateblue',zorder=200)
        axx[2,1].set_xlim([resample_coords[0],resample_coords[-1]])
        axx[2,1].set_xticks(ticks)
        axx[2,1].set_xticklabels(ticklabels, fontsize=panel_font)
        #axx[2,1].set_xlabel(r'Normalized Time $ (t-t_{onset})/t_{δSST_{max}}$', fontsize=panel_font)
        axx[2,1].set_ylabel('$\delta Q_{2m} (g/kg)$', fontsize=panel_font)
        axx[2,1].set_title('Specific Humidity',fontsize=title_font)
        axx[2,1].text(panel_x, panel_y, '(g)', fontsize=panel_font, horizontalalignment='right', verticalalignment='top', transform=axx[2,1].transAxes, zorder=300)
        axx[2,1].yaxis.set_tick_params(labelsize=panel_font)

        
        #Subplot 3,1: Latent fluxes
        #axx[3,1].plot(resample_coords, np.zeros(len(resample_coords)),linewidth = 0.5, zorder=0)
        #axx[3,1].errorbar(x=resample_coords, y=means_DC.δHlC, yerr=stds_DC.δHlC, fmt='-g', ecolor='g', linewidth=0.5, capsize=10*spacing, zorder=100)
        #h1 = axx[3,1].scatter(x=resample_coords, y=means_DC.δHlC, s=sizes_DC, c='forestgreen',zorder=200)
        #axx[3,1].errorbar(x=resample_coords, y=means_DC.δlhf, yerr=stds.δlhf, fmt='-m', ecolor='m', linewidth=0.5, capsize=10*spacing, zorder=100)
        #h2 = axx[3,1].scatter(x=resample_coords, y=means_DC.δlhf, s=sizes_DC, c='purple',zorder=200)
        #axx[3,1].set_xlim([resample_coords[0],resample_coords[-1]])
        #axx[2,0].set_ylim([-80,300])
        #axx[3,1].set_xticks(ticks)
        #axx[3,1].set_xticklabels(ticklabels)
        #axx[3,1].legend([h1,h2],['DC','COARE'])
        #axx[3,1].set_xlabel(r'Normalized Time $ (t-t_{onset})/t_{δSST_{max}}$')
        #axx[3,1].set_ylabel('$\delta$ Flux $(W/m^2)$')
        #axx[3,1].set_title('Latent Fluxes')
        
        
        #Subplot 3,1: Bulk Salinity
        axx[3,1].hlines(y=0,xmin=start-1,xmax=stop+1,color='k',linewidth = 0.5,zorder=0)
        axx[3,1].grid(alpha=0.2)
        axx[3,1].errorbar(x=resample_coords, y=means.δSalTSG, yerr=stds.δSalTSG, fmt='-k', ecolor='k', linewidth=0.5, capsize=10*spacing, zorder=100)
        axx[3,1].scatter(x=resample_coords, y=means.δSalTSG, s=sizes, c='blue',zorder=200)
        axx[3,1].set_xlim([resample_coords[0],resample_coords[-1]])
        #axx[3,1].set_ylim([-0.18,0.28])
        axx[3,1].set_xticks(ticks)
        axx[3,1].set_xticklabels(ticklabels, fontsize=panel_font)
        axx[3,1].set_xlabel(r'Normalized Time $ (t-t_{onset})/t_{δSST_{max}}$', fontsize=panel_font)
        axx[3,1].set_ylabel('$\delta S_{bulk} (psu)$', fontsize=panel_font)
        axx[3,1].set_title('Bulk Salinity (Ship TSG)',fontsize=title_font)
        axx[3,1].text(panel_x, panel_y, '(h)', fontsize=panel_font, horizontalalignment='right', verticalalignment='top', transform=axx[3,1].transAxes, zorder=300)
        axx[3,1].yaxis.set_tick_params(labelsize=panel_font)


        #Subplot 2,0: Sensible Heat Flux Due to Rain
        #axx[2,0].plot(resample_coords, np.zeros(len(resample_coords)),linewidth = 0.5, zorder=0)
        #axx[2,0].errorbar(x=resample_coords, y=means.δrhf, yerr=stds.δrhf, fmt='-k', ecolor='k', linewidth=0.5, capsize=10*spacing, zorder=100)
        #axx[2,0].scatter(x=resample_coords, y=means.δrhf, s=sizes, c='forestgreen',zorder=200)
        #axx[2,0].set_xlim([resample_coords[0],resample_coords[-1]])
        #axx[2,0].set_ylim([-80,300])
        #axx[2,0].set_xticks(ticks)
        #axx[2,0].set_xticklabels(ticklabels)
        #axx[2,0].set_xlabel(r'Normalized Time $ (t-t_{onset})/t_{δSST_{max}}$')
        #axx[2,0].set_ylabel('Rain Heat Flux $(W/m^2)$')
        #axx[2,0].set_title('Sensible Heat Flux Due to Rain (COARE)')
        

        plt.tight_layout()

    return composite_event
