'''
ar_tracking.py
Atmospheric river tracking algorithm based on Brands et al. 2017.
The details of this specific implementation are given the in the Nansen Centre 
technical report #399: DEVELOPMENT OF AN ATMOSPHERIC RIVER DETECTION ALGORITHM.

The boundary line which is checked for exceedance of the triggering thrshold was
created by comapring the line in Brands paper for Southern Noray with ERA-Interim data 
viewed in ncview. 
The boundary line consists of two lines that are:
    vertical:       lon 4.5,  lat 57 to lat 63
    horizontal:     lon 4.5 to lon 9   lat 57

Written by Stephen Outten August 2019
'''

import numpy as np
import dill
import matplotlib.pyplot as plt
import matplotlib
from latlon_disttools import dist_latlon
from datetime import date
from netCDF4 import Dataset


import time as time_pkg   # for testing only
begin_time = time_pkg.time()

# **************************************************
FolderName = '/Data/ERA5/Testing/NorthHemis/IVT/'
FileName = 'IVT_percentiles.dpkl'
FolderNameOut = '/Data/ERA5/Testing/NorthHemis/'
FileNameOut = 'ERA5_AR_BrandsLine.nc'
vmin = 57     # Lat of the vertical boundary line 
vmax = 63
hmin = 4.5    # Lon of the horizontal boundary line
hmax = 9
Pd = 90    #  Pd and Pt are percentiles that ened tob e calculated for each month and gridpoint ??
Pt = 75
minlen = 3000
# **************************************************



def load_percentiles(Pd, Pt):
    '''
    lon, lat, ivt_pd, ivt_pt = load_percentiles(Pd, Pt)
        Loads lon and lat, along with ivt percentiles at levels of Pd and Pt.
        ivt_pd and ivt_pt are maps of the monthly ivt at each grid point for 
        the percentile Pd and Pt respectively. 
    '''
    fn = FolderName + FileName
    datafile = dill.load(open(fn,'rb'))
    lon = datafile['lon']
    lat = datafile['lat']
    ivt_percentiles = datafile['ivt_percentiles']
    ivt_pd = ivt_percentiles[:,Pd,:,:].astype('float32')  # many functions and plotting don't work with float16
    ivt_pt = ivt_percentiles[:,Pt,:,:].astype('float32') 
    del datafile, ivt_percentiles
    return lon, lat, ivt_pd, ivt_pt


def define_line(vmin, vmax, hmin, hmax):
    '''
    define_line(vmin, vmax, hmin, hmax)
        Creates an aray of pairs of indicies for the points on the line.
        Returns line indicies and line lat and lon.
    '''
    # Find vertical lines points
    a = np.arange(vmin, vmax+0.25, 0.25)
    ai = np.where(lat == a[..., np.newaxis])[1]
    b = np.array([hmin] * len(a))
    bi = np.where(lon == b[..., np.newaxis])[1]
    ln_coord = np.column_stack((a,b))
    ln_ind = np.column_stack((ai,bi))
    # Find horizontal line points
    b = np.arange(hmin+0.25, hmax+0.25, 0.25)
    bi = np.where(lon == b[..., np.newaxis])[1]
    a = np.array([vmin] * len(b))
    ai = np.where(lat == a[..., np.newaxis])[1]
    ln_coord = np.vstack((ln_coord, np.column_stack((a,b))))
    ln_ind = np.vstack((ln_ind, np.column_stack((ai,bi))))
    return ln_ind, ln_coord


def load_data(FileName):
    '''
    ivt, ivt_dir = load_data(FileName)
        Loads the ivt and ivt_dir from the dpkl file for the given FileName.
        Returns ivt and ivt_dir from the file.
    '''
    fn = FolderName + FileName
    datafile = dill.load(open(fn,'rb'))
    ivt = datafile['ivt']
    ivt_dir = datafile['ivt_dir']
    del datafile
    return ivt, ivt_dir


def convert_dir(ivt_dir):
    '''
    convert_dir(ivt_dir)
        Converts all of the ivt_dir into eight directions clockwise startingt from 1 as N
        i.e. 0=N, 1=NE, 2=E, 3=SE, 4=S, 5=SW, 6=W, 7=NW
        The method for updating changes the ivt_dir without the need to return it.
    '''
    ivt_dir[(ivt_dir>337.5) | (ivt_dir<=22.5)] = 1
    ivt_dir[(ivt_dir>22.5) & (ivt_dir<=67.5)] = 2
    ivt_dir[(ivt_dir>67.5) & (ivt_dir<=112.5)] = 3
    ivt_dir[(ivt_dir>112.5) & (ivt_dir<=157.5)] = 4
    ivt_dir[(ivt_dir>157.5) & (ivt_dir<=202.5)] = 5
    ivt_dir[(ivt_dir>202.5) & (ivt_dir<=247.5)] = 6
    ivt_dir[(ivt_dir>247.5) & (ivt_dir<=292.5)] = 7
    ivt_dir[(ivt_dir>292.5) & (ivt_dir<=337.5)] = 8
    

def track_ar_step(ind, ivt_grid, ivt_dir_point):
    '''
    new_ivt, new_ind = track_ar_step(ind, ivt_grid, ivt_dir_point)
        Takes in the current location, ivt in 3x3 grid around point, and ivt_dir in grid point. 
        Returns the next grid point point along the ar. 
    '''
    ind_mod = [-1,-1], [-1,0], [-1,1], [0,1], [1,1], [1,0], [1,-1], [0,-1], [-1,-1], [-1,0]
    ind_pos = np.array(ind_mod[ivt_dir_point-1:ivt_dir_point+2])
    ivt_pos = [ivt_grid[1+ind_pos[0][0],1+ind_pos[0][1]], ivt_grid[1+ind_pos[1][0],1+ind_pos[1][1]], ivt_grid[1+ind_pos[2][0],1+ind_pos[2][1]]]
    new_ind = ind + ind_pos[np.argmax(ivt_pos)]
    new_ivt = np.max(ivt_pos)
    return new_ivt, new_ind




# Main Program Start Here
if __name__ == '__main__':
    # Define month lens for years and leap years
    year = np.array([0, 744, 1416, 2160, 2880, 3624, 4344, 5088, 5832, 6552, 7296, 8016, 8760])
    leap_year = np.array([0, 744, 1440, 2184, 2904, 3648, 4368, 5112, 5856, 6576, 7320, 8040, 8784])

    # Load percentile and calculate thresholds for line 
    lon, lat, ivt_pd, ivt_pt = load_percentiles(Pd, Pt)
    ln_ind, ln_coord = define_line(vmin, vmax, hmin, hmax)    
    ln_pd = ivt_pd[:, ln_ind[:,0], ln_ind[:,1]]

    AR_pos = []
    AR_ivt = []
    AR_len = []
    AR_lon_len = []
    AR_time = []
    # START FILE LOOP HERE
    for yy in range(1979,2019):
        start_time = (date(yy,1,1) - date(1900,1,1)).days * 24
        IVT_FileName = 'IVT_' + str(yy) + '.dpkl'
        print("Working on " + IVT_FileName)
        ivt, ivt_dir = load_data(IVT_FileName)   # change from test file to years
    #    ivt, ivt_dir = load_data('IVT_Voss.dpkl')   # change from test file to years
        convert_dir(ivt_dir)
        ln_ivt = ivt[:, ln_ind[:,0], ln_ind[:,1]]   # Extract ivt along line 
    
        # Setup thresholds along line accounting for month
        if ivt.shape[0]==8760:
            yr = year.copy()
        elif ivt.shape[0]==8784:
            yr = leap_year.copy()
    #    yr = np.array([0,24,48,81,95])     # temporary setting for test file  DELETE
        ln_pd_tiled = [np.vstack(np.tile(ln_pd[i,:], (yr[i+1]-yr[i],1)) for i in range(len(yr)-1))][0]

        # Subtract thresholds from ivt along the line, create array of inidices for time and index 
        # along line of possible AR start (i.e. max along line of IVT exceeding Pd)
        ln_ivt_ar = ln_ivt - ln_pd_tiled
        ln_ivt_ar[ln_ivt_ar<0] = 0
    #    ln_ivt_ind = np.vstack((np.arange(len(ivt)), np.argmax(ln_ivt_ar,1))).T    # This tracks from lagest exceedance
    #    ln_ivt_ind = ln_ivt_ind[ln_ivt_ar[ln_ivt_ind[:,0], ln_ivt_ind[:,1]]>0]
        ln_ivt_ind = np.vstack((np.arange(len(ivt)), np.argmax(ln_ivt,1))).T    # This tracks from largest IVT
        ln_ivt_ind = ln_ivt_ind[ln_ivt_ar[ln_ivt_ind[:,0], np.argmax(ln_ivt_ar,1)]>0]

        # Loop over identified possible ARs and confirm AR based on Pt and length
        for ii, lnii in enumerate(ln_ivt_ind):
            ar_pos = ln_coord[lnii[1]]
            ar_len = 0
            ind = ln_ind[lnii[1]]
            ar_ivt = [ivt[lnii[0], ind[0], ind[1]]]
            ar_time = start_time + lnii[0]
            found_ar = True
            
            # Loop as long as AR still passing Pt
            while found_ar is True:
                ivt_grid = ivt[lnii[0], ind[0]-1:ind[0]+2, ind[1]-1:ind[1]+2]     # pull out 3x3 grid
                ivt_dir_point = ivt_dir[lnii[0], ind[0], ind[1]].astype('int')       # pull out dir at grid point
                if ivt_grid.shape == (3,3):                    # This accounts for hittig edge of domain
                    new_ivt, new_ind = track_ar_step(ind, ivt_grid, ivt_dir_point)
                else:
                    found_ar = False
                    break

                if new_ivt > ivt_pt[np.where(yr<=lnii[0])[0][-1], new_ind[0], new_ind[1]]:
                    ar_pos = np.vstack((ar_pos, [lat[new_ind[0]], lon[new_ind[1]]]))
                    ar_len += dist_latlon(ar_pos[-2][0], ar_pos[-2][1], ar_pos[-1][0], ar_pos[-1][1])
                    ar_ivt = np.append(ar_ivt, new_ivt)
                    ivt[lnii[0], ind[0], ind[1]] = 0        # This ensures it cannnot return to previous locations
                    ind = new_ind.copy()
                else:
                    found_ar = False
                    break
        
            # Check if length of AR is greater than minimum length and store if it is
            # Use length of AR since this is always longer than longitudinal extent of AR
            # Hence more ARs stored, can be reduced based on logitudinal extent easily in post-process
            if ar_len>=minlen:
                AR_time.append(ar_time)
                AR_pos.append(ar_pos)
                AR_ivt.append(ar_ivt)
                AR_len.append(ar_len)
                AR_lon_len.append(dist_latlon(np.mean(ar_pos[:,0]), np.min(ar_pos[:,1]) , np.mean(ar_pos[:,0]),  np.max(ar_pos[:,1])))
        
         

    # Change variable length lists into arrays padded with missing value 1e20
    longest = np.max([len(AR_pos[i]) for i in range(len(AR_pos))])
    count_ar = len(AR_len)
    AR_pos_store = 1e20 * np.ones((count_ar, longest, 2)).astype('float32')
    AR_ivt_store = 1e20 * np.ones((count_ar, longest)).astype('float32')
    for ii in range(count_ar):
        AR_pos_store[ii, :len(AR_pos[ii]), :] = AR_pos[ii]
        AR_ivt_store[ii, :len(AR_pos[ii])] = AR_ivt[ii]
    AR_len = np.array(AR_len)
    AR_lon_len = np.array(AR_lon_len)
    AR_time = np.array(AR_time) 



    ##########################################################################
    ###         Write data out to NetCDF                                   ###
    ##########################################################################

    # Determine size of dimensions
    ntime, nlong, ndim = AR_pos_store.shape

    # open a new netCDF file for writing.
    fnout = FolderNameOut + FileNameOut
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
    max_len_out[:] = np.arange(longest) + 1
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
       ar_ivt_out[i,:] = AR_ivt_store[i,:]
       ar_pos_out[i,:,::] = AR_pos_store[i,:,:]

    # close the file.
    ncfile.close()





print('Final run time was {0:.6} hours'.format([time_pkg.time()-begin_time][0]/3600))

