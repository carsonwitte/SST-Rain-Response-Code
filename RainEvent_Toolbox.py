"""
RainEvent_Toolbox.py
"""

#import the python packages which these functions use
import pandas as pd
import xarray as xr
import numpy as np
import datetime as dt
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
            rain_start = timestep
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
                            peak_rate = rainevent.P.max().values.item()
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

    #plot all detected rain events for human review
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

        #slice out sst from this time and make new dataset
        sst_event = xr.Dataset()
        sst_event['sst'] = sst.sel(index=slice(first,last))

        #calculate mean sst of previous 'pre_onset_averaging' minutes
        pre_onset = start - dt.timedelta(minutes = pre_onset_averaging)
        sst_event.attrs['pre-onset mean'] = sst_event.sst.sel(index=slice(pre_onset,start)).mean().item()
        sst_event.attrs['pre-onset std'] = sst_event.sst.sel(index=slice(pre_onset,start)).std().item()

        #calculate max sst deviation from pre-onset mean
        sst_event['δsst'] = sst_event.sst - sst_event.attrs['pre-onset mean']
        try:
            #(if there is no data in this slice of sst, these lines will throw an error. But I want it to just pass a dataset full of nans instead.)
            sst_event.attrs['Max δsst'] = sst_event.δsst.min().values #assumes that a rain event leads to a reduction in SST
            sst_event.attrs['Time of max δsst'] = sst_event.δsst.where(sst_event.δsst==sst_event.attrs['Max δsst'], drop=True).index.values[0]
        except IndexError:
            sst_event.attrs['Max δsst'] = np.NaN
            sst_event.attrs['Time of max δsst'] = np.NaN

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
        #title
        axx[0,event_num].set_title(f'Rain Event # {event_num+1}')

        #-----BOTTOM PLOT: Winds & Heat Fluxes----
        rain_event_list[event_num].lhf.plot.line('-o',color='purple',ax=axx[1,event_num],markersize=3,fillstyle=None)
        rain_event_list[event_num].shf.plot.line('-o',color='firebrick',ax=axx[1,event_num],markersize=3,fillstyle=None)
        axx[1,event_num].set_ylabel('Latent Heat Flux [W/m^2]', color='purple')
        ax1 = axx[1,event_num].twinx()
        rain_event_list[event_num].U10.plot.line('-o',color='slategray',ax=ax1,markersize=3,fillstyle=None)
        ax1.set_ylabel('Wind Speed [m/s]',color='slategray')

    plt.tight_layout()

    return sst_event_list
