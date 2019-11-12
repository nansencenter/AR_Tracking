'''
calc_percentiles.py
Calculates the percentiles for IVT over the reference period of
1979-2009, as per Brands et al. 2017.
These are stored for use by ar_tracking.py
This code stores the percentiles for each location and month and percentile.
i.e. Dimension of output are perceintile[month (12), percentile (100), lat, lon]
This assumes that IVTU and IVTV have been converted using calc_ivt.py.

WARNING: This code takes time to run and requires a LOT of memory.

Written by Stephen Outten August 2019
'''

import numpy as np
import dill


# **************************************************
FolderName = '/Data/ERA5/Testing/NorthHemis/IVT/'
# **************************************************


# Define month lens for years and leap years
year = [0, 744, 1416, 2160, 2880, 3624, 4344, 5088, 5832, 6552, 7296, 8016, 8760]
leap_year = [0, 744, 1440, 2184, 2904, 3648, 4368, 5112, 5856, 6576, 7320, 8040, 8784]

# Load inital file and setup monthly array and get lon and lat
print('Loading initial year...')
fn = FolderName + 'IVT_1979.dpkl'   
datafile = dill.load(open(fn,'rb'))
lon = datafile['lon']
lat = datafile['lat']
ivt = datafile['ivt'].astype('float16')     # as float32 take 160Gb memory and percentile selection adds small variations anyway
del datafile
yr = year.copy()
ivt_jan = ivt[yr[0]:yr[1], :, :]
ivt_feb = ivt[yr[1]:yr[2], :, :]
ivt_mar = ivt[yr[2]:yr[3], :, :]
ivt_apr = ivt[yr[3]:yr[4], :, :]
ivt_may = ivt[yr[4]:yr[5], :, :]
ivt_jun = ivt[yr[5]:yr[6], :, :]
ivt_jul = ivt[yr[6]:yr[7], :, :]
ivt_aug = ivt[yr[7]:yr[8], :, :]
ivt_sep = ivt[yr[8]:yr[9], :, :]
ivt_oct = ivt[yr[9]:yr[10], :, :]
ivt_nov = ivt[yr[10]:yr[11], :, :]
ivt_dec = ivt[yr[11]:yr[12], :, :]
del ivt


# Loop over other files, break into months, and append accordingly
for i in np.arange(1980,2010):
    print('Loading...')
    FileName = 'IVT_' + str(i) + '.dpkl'
    fn = FolderName + FileName
    ivt = dill.load(open(fn,'rb'))['ivt'].astype('float16')
    if ivt.shape[0]==8760:
        print(str(i) + ' is a year')
        yr = year.copy()
    elif ivt.shape[0]==8784:
        print(str(i) + ' is a leap year')
        yr = leap_year.copy()
    else:
        print('Oh crap, oh crap, this is very wrong!!')
    ivt_jan = np.vstack((ivt_jan, ivt[yr[0]:yr[1], :, :]))
    ivt_feb = np.vstack((ivt_feb, ivt[yr[1]:yr[2], :, :]))
    ivt_mar = np.vstack((ivt_mar, ivt[yr[2]:yr[3], :, :]))
    ivt_apr = np.vstack((ivt_apr, ivt[yr[3]:yr[4], :, :]))
    ivt_may = np.vstack((ivt_may, ivt[yr[4]:yr[5], :, :]))
    ivt_jun = np.vstack((ivt_jun, ivt[yr[5]:yr[6], :, :]))
    ivt_jul = np.vstack((ivt_jul, ivt[yr[6]:yr[7], :, :]))
    ivt_aug = np.vstack((ivt_aug, ivt[yr[7]:yr[8], :, :]))
    ivt_sep = np.vstack((ivt_sep, ivt[yr[8]:yr[9], :, :]))
    ivt_oct = np.vstack((ivt_oct, ivt[yr[9]:yr[10], :, :]))
    ivt_nov = np.vstack((ivt_nov, ivt[yr[10]:yr[11], :, :]))
    ivt_dec = np.vstack((ivt_dec, ivt[yr[11]:yr[12], :, :]))
    del ivt


# Create variable holding percentiles
a,b,c = ivt_jan.shape
ivt_percentiles = np.nan * np.ones((12,99,b,c),dtype='float16')
print('Sorting...')
ivt_jan.sort(0)   # sort along time axis
ivt_feb.sort(0)   
ivt_mar.sort(0)   
ivt_apr.sort(0)   
ivt_may.sort(0) 
ivt_jun.sort(0) 
ivt_jul.sort(0) 
ivt_aug.sort(0) 
ivt_sep.sort(0)  
ivt_oct.sort(0)
ivt_nov.sort(0) 
ivt_dec.sort(0)  

# Calculate the percentile values of IVT for each month
print('Percentiling...')
ivt_percentiles[0, :, :, :] = ivt_jan[np.ceil(np.arange(1,100)/100*ivt_jan.shape[0]).astype('int'), :, :]
ivt_percentiles[1, :, :, :] = ivt_feb[np.ceil(np.arange(1,100)/100*ivt_feb.shape[0]).astype('int'), :, :]
ivt_percentiles[2, :, :, :] = ivt_mar[np.ceil(np.arange(1,100)/100*ivt_mar.shape[0]).astype('int'), :, :]
ivt_percentiles[3, :, :, :] = ivt_apr[np.ceil(np.arange(1,100)/100*ivt_apr.shape[0]).astype('int'), :, :]
ivt_percentiles[4, :, :, :] = ivt_may[np.ceil(np.arange(1,100)/100*ivt_may.shape[0]).astype('int'), :, :]
ivt_percentiles[5, :, :, :] = ivt_jun[np.ceil(np.arange(1,100)/100*ivt_jun.shape[0]).astype('int'), :, :]
ivt_percentiles[6, :, :, :] = ivt_jul[np.ceil(np.arange(1,100)/100*ivt_jul.shape[0]).astype('int'), :, :]
ivt_percentiles[7, :, :, :] = ivt_aug[np.ceil(np.arange(1,100)/100*ivt_aug.shape[0]).astype('int'), :, :]
ivt_percentiles[8, :, :, :] = ivt_sep[np.ceil(np.arange(1,100)/100*ivt_sep.shape[0]).astype('int'), :, :]
ivt_percentiles[9, :, :, :] = ivt_oct[np.ceil(np.arange(1,100)/100*ivt_oct.shape[0]).astype('int'), :, :]
ivt_percentiles[10, :, :, :] = ivt_nov[np.ceil(np.arange(1,100)/100*ivt_nov.shape[0]).astype('int'), :, :]
ivt_percentiles[11, :, :, :] = ivt_dec[np.ceil(np.arange(1,100)/100*ivt_dec.shape[0]).astype('int'), :, :]

# Write percentile array to dpkl
print('Writing out...')
fnout = FolderName + 'IVT_percentiles.dpkl'
with open(fnout, 'wb') as f:
    dill.dump({'lat':lat, 'lon':lon, 'ivt_percentiles':ivt_percentiles},f)







