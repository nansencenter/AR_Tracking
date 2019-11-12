'''
identify_voss.py
Loads ar_info produced by reduce_AR_persist.py
Loop over Voss ARs.
Identify ARs with same starting position, and calculate mean and sum of distance 
between them and Voss AR.
Store index (or time), mean, and sum.
AR_table array contains 
    Time of Voss AR
    Time of extracted AR
    Length of Voss AR
    Length of extracted AR
    IVT at end of Voss AR
    IVT at end of extracted AR
    Mean of distance

Tabulate results sorted by Mean of distance
Store as dpkl or csv

Written by Stephen Outten September 2019
'''

import numpy as np
import dill
import matplotlib.pyplot as plt
import matplotlib
from extract_ar import load_ar_data
from latlon_disttools import dist_latlon
from netCDF4 import Dataset
import tabulate
from datetime import datetime, timedelta


# **************************************************
FolderNameAR = '/Data/ERA5/Testing/NorthHemis/'
FileNameAR = 'ERA5_AR_BrandsLine_Reduced.nc'
VossStart = 1006464    # Hours since 01.01.1900 until 26.10.2014:00
VossEnd = 1006559      # Hours since 01.01.1900 until 29.10.2014:23
# **************************************************


# Main Program Start Here
if __name__ == '__main__':
    # Define month lens for years and leap years
    year_times = np.append([692496],np.cumsum(np.tile([8760,8784,8760,8760],10))+692496)
    years = np.arange(1979,2019)

    # Load AR data file
    fn = FolderNameAR + FileNameAR
    AR_time, AR_len, AR_lon_len, AR_pos, AR_ivt = load_ar_data(fn)

    voss_ind = (AR_time>=VossStart) & (AR_time<=VossEnd)
    voss_ind = np.arange(len(voss_ind))[voss_ind]
    voss_times = AR_time[voss_ind]

    # Loop over times including Voss AR
    AR_table = np.array([0,0,0,0,0,0,0,0,0])    # setup array with dummy row
    for i in range(len(voss_times)):
        v_pos = AR_pos[AR_time==voss_times[i],:,:][0]
        e_ind = (AR_pos[:,0,0] == v_pos[0][0]) & (AR_pos[:,0,1] == v_pos[0][1])
        e_ind[voss_ind] = False
        e_ind = np.arange(len(e_ind))[e_ind]
        for j in range(len(e_ind)):
            e_pos = AR_pos[e_ind[j],:,:]
            compare_len = np.min((np.arange(len(e_pos))[np.isnan(e_pos[:,0])][0], np.arange(len(v_pos))[np.isnan(v_pos[:,0])][0]))
            mean_dist = np.mean([dist_latlon(v_pos[posi][0],v_pos[posi][1],e_pos[posi][0],e_pos[posi][1]) for posi in range(compare_len)])
            AR_table = np.vstack((AR_table, [voss_times[i], AR_time[e_ind[j]], AR_len[voss_ind[i]], AR_len[e_ind[j]], AR_len[e_ind[j]]-AR_len[voss_ind[i]], AR_ivt[voss_ind[i],0], AR_ivt[e_ind[j],0], np.abs(AR_ivt[e_ind[j],0]-AR_ivt[voss_ind[i],0]), mean_dist]))
    AR_table = AR_table[1:,:]     # remove dummy row


    # Create small subset of the table only with ARs similar to Voss case
    # IVT difference at last point is <= 30
    # Length difference is less than <= 500 km
    # Mean distance of AR from Voss AR  <= 500 km
    AR_table_small = AR_table[AR_table[:,7]<=8]
    AR_table_small = AR_table_small[np.abs(AR_table_small[:,4])<=500]
    AR_table_small = AR_table_small[AR_table_small[:,8]<=400]
    AR_table_small.shape

    # Tabulate and print out results
    headers = ('Voss_Time', 'AR time', 'Voss length', 'AR length', 'Length diff', 'Voss IVT', 'AR IVT', 'IVT diff', 'Mean Distance')
    AR_table_filtered = []
    for i in range(len(AR_table_small)):
        v_time_stamp = (datetime(1900,1,1,0,0,0) + timedelta(hours=int(AR_table_small[i,0]))).strftime('%Y-%m-%d:%H')
        e_time_stamp = (datetime(1900,1,1,0,0,0) + timedelta(hours=int(AR_table_small[i,1]))).strftime('%Y-%m-%d:%H')
        AR_table_filtered.extend([[v_time_stamp, e_time_stamp, AR_table_small[i,2], AR_table_small[i,3], AR_table_small[i,4], AR_table_small[i,5], AR_table_small[i,6], AR_table_small[i,7], AR_table_small[i,8]]])
    # AR_table_sort.sort(key=lambda x: x[6])    # sorted by mean distace
    AR_table_filtered.sort(key=lambda x: x[7])    # sorted by smallest difference in IVT
    print(tabulate.tabulate(AR_table_filtered, headers))
    

    # Store full resulting data in dill pickle
    fnout = FolderNameAR + 'ARs_like_Voss.dpkl'
    with open(fnout, 'wb') as f:
        dill.dump({'AR_table':AR_table, 'AR_table_small':AR_table_small, 'AR_table_filtered':AR_table_filtered, 'headers':headers},f)


    
    


