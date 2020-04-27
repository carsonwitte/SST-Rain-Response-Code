"""
KT15_Toolbox.py

This python file contains functions for importing, calibrating, and sky-correcting Heitronics KT-15 infrared radiometers:
   1. KT15_importraw: import KT15 data from LabView output (Level 0) and output clean Level 1 product
   2. KT15_labcalibration: import l1 products of blackbody calibration runs and output linear calibration coefficients (slope and y-intercept)
   3. KT15_calibrate_skycorrect82: import l1 products of data runs, apply lab calibration and sky correction. Sky correction is specific to the KT15.82 8um-14um band radiometers; additional functions will be needed for other sky corrections.

To import these functions into a jupyter notebook or other python IDE, run the command:
from KT15_Toolbox import KT15_importraw, KT15_labcalibration, KT15_calibrate_skycorrect82
You will then be able to call each function as if they were defined within your notebook.
"""


def KT15_importraw(data_folder, output_path, sea_serial, sky_serial, experiment):
    # --------------------------------------------------------------------------------------------
    # This function parses an arbitrary number of raw data files that were acquired with the
    # custom-written LabView file 'KT15 UP & DOWN.vi' (stored in 'Software' folder of this
    # repository for reference), gets rid of any bad lines, and saves a single complete
    # L1 data product (.netcdf, with all relevant metadata)
    #
    # Inputs:
    #     data_folder  - a complete filepath to the folder containing the raw data. Relative paths
    #                    and paths to local files will work and may be necessary for speed, but a
    #                    path to a universally accessible server location is best.
    #     output_path  - a filepath to the desired output location of the L1 data product.
    #     sea_serial   - the serial number of the down-looking KT15
    #     sky_serial   - the serial number of the up-looking KT15
    #     experiment   - the name of the experiment (e.g., Falkor19)
    #
    # Outputs: this function does not have any python outputs. Instead it saves the imported data
    # as a netcdf in output_path. The metadata will include:
    #     - the serial numbers of the instruments
    #     - the name of the experiment
    # --------------------------------------------------------------------------------------------

    #import the python packages which this function uses
    import glob
    import pandas as pd
    import xarray as xr
    import numpy as np

    #create list of all text files in data_folder
    files = glob.glob(data_folder + '\*.txt')
    #initialize empty list to add dataframes to as we loop through files
    list = []
    #loop through files and import
    for filepath in files:
        df = pd.read_csv(filepath,                                                                                     #filename to read in
                         delimiter='\s+', skiprows=1, header=None,                                                     #treat whitespace as the delimeter, ignore the header line
                         usecols=[0,1,2,3,4,5], names=['Date','Time','SeaRef','SeaTemp','SkyRef','SkyTemp'],           #use the first 6 columns, and name them as specified
                         parse_dates={'DateTime':[0,1]}, index_col=0,                                                  #parse the first two columns as a single DateTime, and make it the index column
                         na_values=['AMB','TIMEOUT','ERROR'],                                                          #list of other things the parser might encounter in these files, that should be treated like NaNs
                         dtype={'SeaRef':np.float64, 'SeaTemp':np.float64, 'SkyRef':np.float64, 'SkyTemp':np.float64}, #explicitly specify that data columns must be 64-bit floating point numbers
                         error_bad_lines=False, warn_bad_lines=True)                                                   #if there is a bad line in the data file, drop it from the file and show a warning, but continue parsing
        df.dropna(axis='index',how='any',inplace=True)                                                                 #drop any rows that have a NaN value in them
        list.append(df)                                                                                                #append this newly imported dataframe to a list of dataframes
    #concatenate all the dataframes in this list into one dataframe along the vertical axis
    kt = pd.concat(list, axis=0)

    #convert to xarray dataset and add metadata attributes so that the data file describes itself
    kt_xr = xr.Dataset.from_dataframe(kt)
    kt_xr.attrs['experiment'] = experiment
    kt_xr.attrs['sea_serial'] = sea_serial
    kt_xr.attrs['sky_serial'] = sky_serial
    #save as netcdf in output_path
    kt_xr.to_netcdf(output_path+f'/{experiment}_KT15_{sea_serial}_{sky_serial}.cdf')

#def KT15_labcalibration(l1data_path):
    # --------------------------------------------------------------------------------------------
    # This function takes KT15 data from blackbody calibrations done in the lab, and calculates
    # a linear fit to the measured values given the blackbody reference temperature.
    # --------------------------------------------------------------------------------------------



def KT15_calibrate_skycorrect82(l1data_path, sea_slope, sea_yint, sky_slope, sky_yint):
    # --------------------------------------------------------------------------------------------
    # This function takes data from a pair of up- and down-looking KT15.82 8um-14um infrared
    # radiometers, applies a linear calibration based on mean coefficients provided from blackbody
    # calibrations in the lab, and calculates a final value for the true SST corrected for
    # reflections from the sky, following SOURCE.
    #
    # Inputs:
    #    l1data_path - path to Level 1 Data Product generated by KT15_import function
    #    sea_slope   - slope of blackbody calibration for down-looking instrument
    #    sea_yint    - y-intercept of blackbody calibration for down-looking instrument
    #    sky_slope   - slope of blackbody calibration for up-looking instrument
    #    sky_yint    - y-intercept of blackbody calibration for up-looking instrument
    #
    # Outputs:
    #    SST         - final radiometric Sea Surface Temperature for the skin layer
    #
    # Dependencies:
    #    numpy
    #    pandas
    #    matplotlib?
    # --------------------------------------------------------------------------------------------

    # open

    #apply blackbody calibration
    SkyTemp_Cal = (kt.SkyTemp - sky_yint)/sky_slope
    SeaTemp_Cal = (kt.SeaTemp - sea_yint)/sea_slope

    #apply sky correction - what is the source for these coefficients?
    K0u = -1.9151
    K1u = 0.041767
    K2u = -0.00033291
    K3u = 1.0715e-06
    K4u = -8.651e-10
    K5u = K6u = K7u = K8u = 0

    K0d = -82.135
    K1d = 2.293
    K2d = -0.020287
    K3d = 6.7301e-05
    K4d = -5.2395e-08
    K5d = K6d = K7d = K8d = 0

    y0 = 136.33
    A = 20.927
    power = 0.40045

    RadTempSKY = kt.SkyTemp_Cal + 273.16
    RadTempSEA = kt.SeaTemp_Cal + 273.16
    ReflectedSkyRad = K0u + K1u*RadTempSKY + K2u*RadTempSKY**2 + K3u*RadTempSKY**3 + K4u*RadTempSKY**4 + K5u*RadTempSKY**5 + K6u*RadTempSKY**6 + K7u*RadTempSKY**7 + K8u*RadTempSKY**8
    TotalRad = K0d + K1d*RadTempSEA + K2d*RadTempSEA**2 + K3d*RadTempSEA**3 + K4d*RadTempSEA**4 + K5d*RadTempSEA**5 + K6d*RadTempSEA**6 + K7d*RadTempSEA**7 + K8d*RadTempSEA**8
    SSTRad = TotalRad - ReflectedSkyRad
    kt['SST'] = y0 + A*(SSTRad)**power - 273.16
