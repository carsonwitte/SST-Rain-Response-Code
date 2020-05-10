"""
KT15_Toolbox.py

This python file contains functions for importing, calibrating, and sky-correcting Heitronics KT-15 infrared radiometers:
   1. KT15_reader: read a single Level 0 data file and return a pandas dataframe
   2. KT15_importraw: import KT15 data from LabView output (Level 0) and output clean Level 1 product
   3. KT15_labcalibration: import L0 products of blackbody calibration runs and save calibration data + linear calibration coefficients (slope and y-intercept)
   4. KT15_calibrate_skycorrect82: import l1 products of data runs, apply lab calibration and sky correction. Sky correction is specific to the KT15.82 8um-14um band radiometers; additional functions will be needed for other sky corrections.

To import these functions into a jupyter notebook or other python IDE, run the command:
from KT15_Toolbox import KT15_importraw, KT15_labcalibration, KT15_calibrate_skycorrect82, etc.
You will then be able to call each function as if they were defined within your notebook.
"""

#import the python packages which these functions use
import glob
import pandas as pd
import xarray as xr
import numpy as np
import datetime as dt
from matplotlib import pyplot as plt

def KT15_reader(filepath):
    '''
    This function just reads in a single Level-0 text file that was acquired with the custom-written LabView file 'KT15 UP & DOWN.vi' (stored in 'Software' folder of this repository for reference), gets rid of any bad lines, and returns a pandas dataframe

    Input: filepath - string pointing to the data file
    Output: df      - pandas dataframe containing cleaned timeseries
    '''

    # define a function to handle the importing of bad lines
    #if the data is not a valid float, return a NaN
    def valid_float(y):
        try:
            return float(y)
        except ValueError:
            return np.nan
        except TypeError:
            return np.nan

    df = pd.read_csv(filepath, engine='python',                                                                              #filename to read in, python engine is slower but better error handling (default is C)
                     delimiter='\s+', skiprows=1, header=None,                                                               #treat whitespace as the delimeter, ignore the header line
                     usecols=[0,1,2,3,4,5], names=['Date','Time','SeaRef','SeaTemp','SkyRef','SkyTemp'],                     #use the first 6 columns, and name them as specified
                     parse_dates={'DateTime':[0,1]}, index_col=0,                                                            #parse the first two columns as a single DateTime, and make it the index column
                     na_values=['AMB','TIMEOUT','ERROR','BAD','C'],                                                          #list of other things the parser might encounter in these files, that should be treated like NaNs
                     converters={'SeaRef':valid_float, 'SeaTemp':valid_float, 'SkyRef':valid_float, 'SkyTemp':valid_float},  #explicitly specify that data columns must be floats, and if not a nan is put in (should make the above line obsolete)
                     error_bad_lines=False, warn_bad_lines=True)                                                             #if there is a bad line in the data file, drop it from the file and show a warning, but continue parsing
    df.dropna(axis='index',how='any',inplace=True)                                                                           #drop any rows that have a NaN value in them
    return df

######################################################################################################

def KT15_importraw(data_folder, output_path, sea_serial, sky_serial, experiment):
    '''
    This function parses an arbitrary number of raw KT15 data files, and saves a single complete L1 data product (.netcdf, with all relevant metadata).

    Inputs:
         data_folder  - a complete filepath to the folder containing the raw data.
         output_path  - a filepath to the desired output location of the L1 data product.
         sea_serial   - the serial number of the down-looking KT15
         sky_serial   - the serial number of the up-looking KT15
         experiment   - the name of the experiment (e.g., Falkor19)

    Outputs: this function does not have any python outputs. Instead it saves the imported data as a netcdf in output_path. The metadata will include:
         - the serial numbers of the instruments
         - the name of the experiment
    '''

    #create list of all text files in data_folder
    files = glob.glob(data_folder + '\*.txt')
    #initialize empty list to add dataframes to as we loop through files
    list = []
    #loop through files and import
    for filepath in files:
        df = KT15_reader(filepath)
        list.append(df)
    #concatenate all the dataframes in this list into one dataframe along the vertical axis
    kt = pd.concat(list, axis=0, sort=True)

    #convert to xarray dataset and add metadata attributes so that the data file describes itself
    kt_xr = xr.Dataset.from_dataframe(kt)
    kt_xr.attrs['experiment'] = experiment
    kt_xr.attrs['sea_serial'] = sea_serial
    kt_xr.attrs['sky_serial'] = sky_serial
    #save as netcdf in output_path
    kt_xr.to_netcdf(output_path+f'/{experiment}_KT15_{sea_serial}_{sky_serial}.cdf')

######################################################################################################

def KT15_labcalibration(data_folder, output_path, caltemps):
    '''
    This function takes KT15 data from a blackbody calibration done in the lab, and calculates a linear fit to the measured values given the blackbody reference temperature. To apply the calibration calculated by this function, apply this formula:
        calibrated_data = (raw_data - yint)/slope

    Inputs:
        data_folder - path to the folder containing exactly one blackbody calibration
        output_path - a filepath to the desired output location of the calibration data
        caltemps    - list of temperatures used in the calibration (usually 10-40 by 5)

    This function also saves each instrument calibration as a netcdf in output_path.
    The metadata for each will include:
         - the serial number of the instrument
         - the date of the blackbody calibration
         - the slope of the fitted linear calibration
         - the y-intercept of the fitted linear calibration
    '''

    #read list of data files
    files = glob.glob(data_folder + '\*.txt')

    #the filename should always be formatted as 'YYYY_LabCalib_INST_SeaSerial_SkySerial_BBTemp_Yearday_HHMMSS.txt'
    #we want to extract the year and yearday from it once for the date, plus the sky and sea serials
    parsed_filename = files[0].replace('\\','_').split('_')
    sea_serial = parsed_filename[-5]
    sky_serial = parsed_filename[-4]
    yday = int(parsed_filename[-2])
    year = int(parsed_filename[-8])
    date = dt.datetime(year, 1, 1) + dt.timedelta(yday - 1)
    datestr = date.strftime('%Y-%m-%d')

    #initialize pandas dataframe with calibration temperatures as the index, and mean & std as the columns
    seaCal = pd.DataFrame(data=None,index=caltemps,columns=['Mean','Std'])
    skyCal = pd.DataFrame(data=None,index=caltemps,columns=['Mean','Std'])

    for filepath in files:
        #parse the temperature field in the filename
        bb_temp = int(filepath.split('_')[-3])

        #read in the file
        kt = KT15_reader(filepath)

        #populate calibration dataframe for this temperature level
        seaCal.Mean.loc[bb_temp] = kt.SeaTemp.mean()
        seaCal.Std.loc[bb_temp] = kt.SeaTemp.std()
        skyCal.Mean.loc[bb_temp] = kt.SkyTemp.mean()
        skyCal.Std.loc[bb_temp] = kt.SkyTemp.std()

    #fit a line to the bb temps vs. the mean observed temps
    [sea_slope, sea_yint] = np.polyfit(caltemps,seaCal.Mean,deg=1)
    [sky_slope, sky_yint] = np.polyfit(caltemps,skyCal.Mean,deg=1)

    #convert to xarray, add metadata, and save as netcdf in output_path
    sea_xr = xr.Dataset.from_dataframe(seaCal)
    sea_xr.attrs['serial_number'] = sea_serial
    sea_xr.attrs['calibration_date'] = datestr
    sea_xr.attrs['slope'] = sea_slope
    sea_xr.attrs['yint'] = sea_yint
    sea_xr.to_netcdf(output_path+f'/KT15_LabCalibration_{sea_serial}_{datestr}.cdf')

    sky_xr = xr.Dataset.from_dataframe(skyCal)
    sky_xr.attrs['serial_number'] = sky_serial
    sky_xr.attrs['calibration_date'] = datestr
    sky_xr.attrs['slope'] = sky_slope
    sky_xr.attrs['yint'] = sky_yint
    sky_xr.to_netcdf(output_path+f'/KT15_LabCalibration_{sky_serial}_{datestr}.cdf')

    #plot data and calibrated line for each one, with 2*Stdev error bars
    fig,axx = plt.subplots(nrows=1,ncols=2,figsize=(10,4),facecolor='w')

    axx[0].errorbar(x = caltemps, y = seaCal.Mean, yerr = 2*seaCal.Std, fmt = '.', markersize=6, markeredgewidth=2)
    axx[0].plot(caltemps, (seaCal.Mean-sea_yint)/sea_slope)
    axx[0].set_title(f'KT-15 {sea_serial} Calibration {datestr}')

    axx[1].errorbar(x = caltemps, y = skyCal.Mean, yerr = 2*skyCal.Std, fmt = '.', markersize=6, markeredgewidth=2)
    axx[1].plot(caltemps, (skyCal.Mean-sky_yint)/sky_slope)
    axx[1].set_title(f'KT-15 {sky_serial} Calibration {datestr}')

    for axis in axx:
        axis.grid()
        axis.set_ylabel('Measured Temperature ($^\circ$C)')
        axis.set_xlabel('Blackbody Temperature ($^\circ$C)')
        axis.legend(['Calibration','Data (2$\sigma$ error bars)'])

    #save the calibration figure in the output_path for reference
    plt.savefig(output_path+f'/KT15_Calibration_{datestr}_{sea_serial}_{sky_serial}.png')

######################################################################################################

def KT15_calibrate_skycorrect82(l1data_path, output_path, cal_path):
    '''
    This function takes data from a pair of up- and down-looking KT15.82 8um-14um infrared radiometers, applies a linear calibration based on mean coefficients taken from blackbody calibrations in the lab, and calculates a final value for the 'true' SST corrected for reflections from the sky, following SOURCE.

    Inputs:
        l1data_path - path to Level 1 Data Product generated by KT15_import function
        output_path - path to save new product with final calibrated SST
        cal_path    - path to folder containing all the blackbody calibrations (generated by
                      KT15_labcalibration) that you'd like to apply (for both instruments)

    Outputs: The imported L1 netcdf is saved in the output_path with a new variable added:
        SST         - final radiometric Sea Surface Temperature for the skin layer
    '''

    # open L1 data
    kt = xr.open_dataset(l1data_path)

    #extract all calibration files for each instrument from cal_path, and take the average values
    sky_files = glob.glob(cal_path+f'/*{kt.sky_serial}*.cdf')
    sky_slopes = []
    sky_yints = []
    for cal_file in sky_files:
        cal = xr.open_dataset(cal_file)
        sky_slopes.append(cal.slope)
        sky_yints.append(cal.yint)
    sky_slope = np.mean(sky_slopes)
    sky_yint = np.mean(sky_yints)

    sea_files = glob.glob(cal_path+f'/*{kt.sea_serial}*.cdf')
    sea_slopes = []
    sea_yints = []
    for cal_file in sea_files:
        cal = xr.open_dataset(cal_file)
        sea_slopes.append(cal.slope)
        sea_yints.append(cal.yint)
    sea_slope = np.mean(sea_slopes)
    sea_yint = np.mean(sea_yints)

    #apply blackbody calibration
    kt['SkyTemp_Cal'] = (kt.SkyTemp - sky_yint)/sky_slope
    kt['SeaTemp_Cal'] = (kt.SeaTemp - sea_yint)/sea_slope
    kt.attrs['sea_slope'] = sea_slope
    kt.attrs['sea_yint'] = sea_yint
    kt.attrs['sky_slope'] = sky_slope
    kt.attrs['sky_yint'] = sky_yint

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

    #plot the results and save the figure
    fig,axx = plt.subplots(facecolor='w')
    kt.SST.plot(ax=axx)
    plt.title(f'KT-15 Calibrated & Sky-Corrected Data: {kt.experiment}')
    plt.savefig(output_path+f'/{kt.experiment}_KT15_CalibratedSST.png')

    #save the final dataframe in the output path
    kt.to_netcdf(output_path+f'/{kt.experiment}_KT15_CalibratedSST.cdf')
