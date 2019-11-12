'''
reduce_ar_persist.py
Loads the tracking information produced by ar_tracking.py.
Reduce to ARs that have minimum persistence in time and identify minimum domain.
Stores the information in reduced info file.

Written by Stephen Outten September 2019
'''

import numpy as np
import dill
import matplotlib.pyplot as plt
import matplotlib
from netCDF4 import Dataset


# **************************************************
FolderNameERA = '/Data/ERA5/Testing/NorthHemis/'
FolderNameAR = '/Home/stepoutt/Python/MUDYFEET/Data/'
FileNameAR = 'ERA5_AR_BrandsLine.nc'
FileNameAR_Out = 'ERA5_AR_BrandsLine_Reduced.nc'
persist = 18   # Number of hours an AR must exist for to be included in data extraction
separate = 24  # Number of hours an AR must be separated from another AR to not be considered part of the same event
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


def reduce_AR_persist(AR_time, persist, separate):
    '''
    AR_persist_ind = reduce_AR_persist(AR_time, persist, separate)
        Scans through the times provided and isolates those AR which persist for
        at least the minimium time (given by the persist variables), and which are
        separated by at least 'seperate' hours from other events. 
        If something is less that 'persist' hours long, but it only has a break 
        of a 'separate' hours before another AR appears, they are considered to be part 
        of the same AR.
        Returns a set of indicies for which ARs meet the persist requirement.
    '''
    diff = AR_time[1:] - AR_time[:-1]
    diff[diff<separate] = 1
    diff = np.append(np.append([0], diff), [0])
    changes = np.arange(len(diff))[diff>1]
    start = changes[diff[changes+1] == 1]
    end = changes[diff[changes-1] == 1]
    if len(end)<len(start): end = np.append(end,len(AR_time))
    ar_persist_ind = np.vstack((start,end)).T
    ar_persist_ind = ar_persist_ind[ar_persist_ind[:,1] - ar_persist_ind[:,0] >= persist]
    return ar_persist_ind


# Main Program Start Here
if __name__ == '__main__':
    # Define month lens for years and leap years
    year_times = np.append([692496],np.cumsum(np.tile([8760,8784,8760,8760],10))+692496)
    years = np.arange(1979,2019)

    # Load AR data file
    fn = FolderNameAR + FileNameAR
    AR_time, AR_len, AR_lon_len, AR_pos, AR_ivt = load_ar_data(fn)

    # Identify start and end of ARs based on hours of persistence and separation
    # And reduce AR variables accordingly
    AR_persist_ind = reduce_AR_persist(AR_time, persist, separate)
    AR_time = np.concatenate([AR_time[AR_persist_ind[i,0]:AR_persist_ind[i,1]] for i in range(len(AR_persist_ind))], axis=0)
    AR_len = np.concatenate([AR_len[AR_persist_ind[i,0]:AR_persist_ind[i,1]] for i in range(len(AR_persist_ind))], axis=0)    
    AR_lon_len = np.concatenate([AR_lon_len[AR_persist_ind[i,0]:AR_persist_ind[i,1]] for i in range(len(AR_persist_ind))], axis=0)
    AR_ivt = np.concatenate([AR_ivt[AR_persist_ind[i,0]:AR_persist_ind[i,1]] for i in range(len(AR_persist_ind))], axis=0)    
    AR_pos = np.concatenate([AR_pos[AR_persist_ind[i,0]:AR_persist_ind[i,1]] for i in range(len(AR_persist_ind))], axis=0)


    ##########################################################################
    ###         Write AR info out to NetCDF                                ###
    ##########################################################################

    # Determine size of dimensions
    ntime, nlong, ndim = AR_pos.shape

    # open a new netCDF file for writing.
    fnout = FolderNameAR + FileNameAR_Out
    ncfile = Dataset(fnout, 'w')

    # create the lat and lon dimensions.
    time = ncfile.createDimension('time', None)
    max_length = ncfile.createDimension('max_length', nlong)
    index = ncfile.createDimension('index', ndim)

    # Define the coordinate variables.
    time_out = ncfile.createVariable('time', np.float32, ('time',))
    max_len_out = ncfile.createVariable('max_length', np.float32, ('max_length'))
    index_out = ncfile.createVariable('index', np.float32, ('index'))

    # Assign units attributes to coordinate variable data.
    time_out.units = 'hours since 1900-01-01 00:00:00.0'
    max_len_out.units = ''
    index_out.units = ''
    time_out.long_name = 'Time of AR occurrence'
    max_len_out.long_name = 'Number of gridpoints of longest AR'
    index_out.long_name = 'Non-dimensional array index'

    # write data to coordinate vars.
    time_out[:] = AR_time
    max_len_out[:] = np.arange(nlong) + 1
    index_out[:] = np.arange(2) + 1

    # create main variables
    ar_len_out = ncfile.createVariable('AR_len', np.float32, ('time'))
    ar_lon_len_out = ncfile.createVariable('AR_lon_len',np.float32, ('time'))
    ar_ivt_out = ncfile.createVariable('AR_ivt',np.float32, ('time','max_length'))
    ar_pos_out = ncfile.createVariable('AR_position',np.float32, ('time','max_length','index'))

    # set the units attribute.
    ar_len_out.units = 'km'
    ar_lon_len_out.units = 'km'
    ar_ivt_out.units = 'kg m**-1 s**-1'
    ar_pos_out.units = 'lat, lon'
    ar_len_out.long_name = 'Length of AR'
    ar_lon_len_out.long_name = 'Longitudinal extent of AR at mean latitude'
    ar_ivt_out.long_name = 'Vertical integral of water vapour flux along AR'
    ar_pos_out.long_name = 'Latitude and longitude of gridpoints along AR'

    # write data to variables along record (unlimited) dimension.
    ar_len_out[:] = AR_len
    ar_lon_len_out[:] = AR_lon_len
    for i in range(ntime):
       ar_ivt_out[i,:] = AR_ivt[i,:]
       ar_pos_out[i,:,::] = AR_pos[i,:,:]

    # close the file.
    ncfile.close()

    


            
 


