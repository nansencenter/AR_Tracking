'''
extract_ar.py
Loads the tracking information produced by ar_tracking.py.
Reduce to ARs that have minimum persistence in time and identify minimum domain.
Extract variable over given domain for ARs that meet criteria.

Written by Stephen Outten August 2019
'''

import numpy as np
import dill
import matplotlib.pyplot as plt
import matplotlib
from netCDF4 import Dataset


# **************************************************
FolderNameERA = '/Data/ERA5/Testing/NorthHemis/'
FolderNameAR = '/Data/ERA5/Testing/NorthHemis/'
FolderNameOut = '/Data/ERA5/Testing/NorthHemis/AR_Extract/'
FileNameAR = 'ERA5_AR_BrandsLine_Reduced.nc'
persist = 18   # Number of hours an AR must exist for to be included in data extraction
separate = 24  # Number of hours an AR must be separated from another AR to not be considered part of the same event
var_list = 'MSLP', 'Precip'    # List of variables to extract excluding IVT which is handled separately 
var_names = 'msl', 'tp'        # List of names of variables as they appear in the oroginal ERA5 NetCDF files
# **************************************************


def load_ar_data(fn):
    '''
    AR_time, ar_len, ar_lon_len, ar_pos, ar_ivt = load_ar_data(fn)
        Loads the variables from the data fiels given by the path fn.
        Variables are time in hours since 01-01-1900 that ARs occur,
        length and longitudinal extent of ARs, lat/lon position of ARs,
        and IVT along length of ARs.
    '''
    dataset = Dataset(fn)
    ar_time = np.array(dataset.variables['time'][:])
    ar_len = np.array(dataset.variables['AR_len'][:])
    ar_lon_len = np.array(dataset.variables['AR_lon_len'][:])
    ar_ivt = np.array(dataset.variables['AR_ivt'][:])
    ar_pos = np.array(dataset.variables['AR_position'][:])
    ar_pos[ar_pos==1e20] = np.nan
    dataset.close() 
    return ar_time, ar_len, ar_lon_len, ar_pos, ar_ivt




# Main Program Start Here
if __name__ == '__main__':
    # Define month lens for years and leap years
    year_times = np.append([692496],np.cumsum(np.tile([8760,8784,8760,8760],10))+692496)
    years = np.arange(1979,2019)

    # Load AR data file
    fn = FolderNameAR + FileNameAR
    AR_time, AR_len, AR_lon_len, AR_pos, AR_ivt = load_ar_data(fn)

    '''
    Checked min and max lat and lon using the following:
      np.nanmax([np.nanmax(AR_pos[ar_persist_ind[i,0]:ar_persist_ind[i,1],:,0]) for i in range(len(ar_persist_ind))])
      np.nanmax([np.nanmax(AR_pos[ar_persist_ind[i,0]:ar_persist_ind[i,1],:,1]) for i in range(len(ar_persist_ind))])
    ARs reach the boundaries on all sides except the north, where the boundary is 80 N and the most northern AR 
    only reaches 77.5 N     
    Based on this finding, extract the data from the entire domain.
    '''

    # Loop over variables to extrtact and years
    # Load data for current year and extract data for AR times
    for var in range(len(var_list)):
        fn = FolderNameERA + var_list[var] + '/' + var_list[var] + '_1979.nc'
        dataset = Dataset(fn)
        lat = np.array(dataset.variables['latitude'][:])
        lon = np.array(dataset.variables['longitude'][:])
        var_units = dataset.variables[var_names[var]].units
        var_longname = dataset.variables[var_names[var]].long_name
        dataset.close()

        for yr in range(len(years)):
            fn = FolderNameERA + var_list[var] + '/' + var_list[var] + '_' + str(years[yr]) + '.nc'
            dataset = Dataset(fn)
            var_data = np.array(dataset.variables[var_names[var]][:])
            var_time = np.array(dataset.variables['time'][:])        
            dataset.close()
            if years[yr]==1979:
                ar_data = var_data[np.isin(var_time,AR_time),:,:]
            else:
                ar_data = np.vstack((ar_data, var_data[np.isin(var_time,AR_time),:,:]))
            del var_data, var_time

        # Output each variable to a separate NetCDF file
        ##########################################################################
        ###         Write data out to NetCDF                                   ###
        ##########################################################################

        # Determine size of dimensions
        ntime, nlat, nlon = ar_data.shape

        # open a new netCDF file for writing.
        fnout = FolderNameOut + var_list[var] + '_ARs_1979_2018.nc'
        ncfile = Dataset(fnout, 'w')

        # create the lat and lon dimensions.
        time = ncfile.createDimension('time', None)
        latitude = ncfile.createDimension('latitude', nlat)
        longitude = ncfile.createDimension('longitude', nlon)

        # Define the coordinate variables.
        time_out = ncfile.createVariable('time', np.float32, ('time',))
        latitude_out = ncfile.createVariable('latitude', np.float32, ('latitude'))
        longitude_out = ncfile.createVariable('longitude', np.float32, ('longitude'))

        # Assign units attributes to coordinate variable data.
        time_out.units = 'hours since 1900-01-01 00:00:00.0'
        latitude_out.units = 'degrees_north'
        longitude_out.units = 'degrees_east'
        time_out.long_name = 'time'
        time_out.calendar = 'gregorian'
        latitude_out.long_name = 'latitude'
        longitude_out.long_name = 'longitude'

        # write data to coordinate vars.
        time_out[:] = AR_time
        latitude_out[:] = lat
        longitude_out[:] = lon

        # create main variables
        ar_data_out = ncfile.createVariable(var_names[var], np.float32, ('time','latitude','longitude'))

        # set the variable attributes.
        ar_data_out.units = var_units
        ar_data_out.long_name = var_longname

        # write data to variables along record (unlimited) dimension.
        for i in range(ntime):
           ar_data_out[i,:,::] = ar_data[i,:,:]

        # close the file.
        ncfile.close()


    for yr in range(len(years)):
        print(years[yr])
        fn = FolderNameERA + 'IVT/IVT_' + str(years[yr]) + '.dpkl'
        datafile = dill.load(open(fn,'rb'))
        ivt = datafile['ivt']
        ivt_dir = datafile['ivt_dir']
        del datafile
        ivt_time = np.arange(year_times[yr],year_times[yr+1])
        if years[yr]==1979:
            ar_ivt = ivt[np.isin(ivt_time,AR_time),:,:]
            ar_ivt_dir = ivt_dir[np.isin(ivt_time,AR_time),:,:]
        else:
            ar_ivt = np.vstack((ar_ivt, ivt[np.isin(ivt_time,AR_time),:,:]))
            ar_ivt_dir = np.vstack((ar_ivt_dir, ivt_dir[np.isin(ivt_time,AR_time),:,:]))
        del ivt, ivt_dir, ivt_time    


    ##########################################################################
    ###         Write IVT and IVT_dir out to NetCDF                        ###
    ##########################################################################

    # Determine size of dimensions
    ntime, nlat, nlon = ar_ivt.shape

    # open a new netCDF file for writing.
    fnout = FolderNameOut + 'IVT_ARs_1979_2018.nc'
    ncfile = Dataset(fnout, 'w')

    # create the lat and lon dimensions.
    time = ncfile.createDimension('time', None)
    latitude = ncfile.createDimension('latitude', nlat)
    longitude = ncfile.createDimension('longitude', nlon)

    # Define the coordinate variables.
    time_out = ncfile.createVariable('time', np.float32, ('time',))
    latitude_out = ncfile.createVariable('latitude', np.float32, ('latitude'))
    longitude_out = ncfile.createVariable('longitude', np.float32, ('longitude'))

    # Assign units attributes to coordinate variable data.
    time_out.units = 'hours since 1900-01-01 00:00:00.0'
    latitude_out.units = 'degrees_north'
    longitude_out.units = 'degrees_east'
    time_out.long_name = 'time'
    time_out.calendar = 'gregorian'
    latitude_out.long_name = 'latitude'
    longitude_out.long_name = 'longitude'

    # write data to coordinate vars.
    time_out[:] = AR_time
    latitude_out[:] = lat
    longitude_out[:] = lon

    # create main variables
    ar_ivt_out = ncfile.createVariable('IVT', np.float32, ('time','latitude','longitude'))
    ar_ivt_dir_out = ncfile.createVariable('IVT_dir', np.float32, ('time','latitude','longitude'))

    # set the variable attributes.
    ar_ivt_out.units = "kg m**-1 s**-1"
    ar_ivt_out.long_name = "Vertical integral of water vapour flux"
    ar_ivt_dir_out.units = "degrees"
    ar_ivt_dir_out.long_name = "Direction of vertical integral of water vapour flux"

    # write data to variables along record (unlimited) dimension.
    for i in range(ntime):
       ar_ivt_out[i,:,::] = ar_ivt[i,:,:]
       ar_ivt_dir_out[i,:,::] = ar_ivt_dir[i,:,:]

    # close the file.
    ncfile.close()

