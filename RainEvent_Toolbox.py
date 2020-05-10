"""
RainEvent_Toolbox.py
"""

#import the python packages which these functions use
import pandas as pd
import xarray as xr
import numpy as np
import datetime as dt
from matplotlib import pyplot as plt

def find_rain_events(dataset, min_duration, min_separation, front_pad, end_pad):
    '''
    This function splits an xarray dataset into individual datasets for each precipitation event found in the timeseries, based on the input parameters.

    Inputs:
        dataset        - xarray dataset containing a precipitation field labeled 'P'
        min_duration   - minimum consecutive timesteps of rainfall required to be considered an event start
        min_separation - minimum consecutive timesteps of no rainfall required to be considered an event end
        front_pad      - how many minutes prior to rain onset to include in the rain event
        end_pad        - how many minutes after rain end to include in the rain event

    Outputs:
        rain_event_list - a list of xarray datasets, each of which contains one rain event and has the same format as 'dataset'
    '''

    rain_rate = dataset.P
    rain_event_list = []
    rain_event_counter = 0
    timestep = 0
    while timestep < len(rain_rate):
        #if the current rain rate isn't zero...AND...it's going to keep raining for 'min_duration':
        if (rain_rate.isel(index=timestep).values != 0 and
            np.min(rain_rate.isel(index=slice(timestep,timestep+min_duration)).values) > 0):

            #then we've got ourselves the start of a rain event!
            rain_start = timestep
            rain_event_counter = rain_event_counter + 1
            print(f'Rain Event {rain_event_counter} found at {rain_rate.isel(index=rain_start).index.values}')
            raining = True
            #we've already checked that the rain event lasts at least 'min_duration', so start searching for the end there
            timestep = timestep + min_duration - 1
            #enter a new loop that finds the end index of this rain event
            while raining:
                #if the current rain rate is zero...AND...it's going to be zero from now until 'min_separation' from now:
                if (rain_rate.isel(index=timestep).values == 0 and
                    np.mean(rain_rate.isel(index=slice(timestep,timestep+min_separation)).values) == 0):
                        #then we've located the end of the rain event
                        rain_end = timestep
                        raining = False
                        #create new dataarray containing solely this rain event, plus padding on either side
                        rainevent = dataset.isel(index=slice(rain_start-front_pad,rain_end+end_pad))
                        rainevent.attrs['Rain Event #'] = rain_event_counter
                        rainevent.attrs['Rain Onset'] = rain_rate.isel(index=rain_start).index.values
                        rainevent.attrs['Rain End'] = rain_rate.isel(index=rain_end).index.values
                        #add this dataframe to the master list of rain events
                        rain_event_list.append(rainevent)
                #as long as it keeps raining, just increment the timestep and start over
                else:
                    timestep = timestep + 1
            #[spits us out at the timestep after the event ends]
        #if we haven't found a rain event, just increment the timestep and start over
        else:
            timestep = timestep + 1

    return rain_event_list

#Max SST response
#Time of max SST response
#Cumulative rain
#Length of rain event
#Time of maximum rain rate
