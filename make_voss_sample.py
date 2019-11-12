'''
make_voss_sample.py
Loads dpkl data for entire 2014, extracts 4 days of 26th to 29th October
around Voss event, outputs them as sample for testing.

Written by Stephen Outten August 2019
'''

import numpy as np
import dill
import matplotlib.pyplot as plt
import matplotlib

FolderName = '/Data/ERA5/Testing/NorthHemis/IVT/'
FileName = 'IVT_2014.dpkl'

# Load lon and lat from any file
fn = FolderName + FileName
datafile = dill.load(open(fn,'rb'))
lon = datafile['lon']
lat = datafile['lat']
ivt = datafile['ivt']
ivt_dir = datafile['ivt_dir']
del datafile

ivt_voss = ivt[7153:7248,:,:].copy()
ivt_dir_voss = ivt_dir[7153:7248,:,:].copy()

# Store lon, lat, ivt, and ivt_dir in dill pickle
fnout = FolderName + 'IVT_Voss.dpkl'
with open(fnout, 'wb') as f:
    dill.dump({'lat':lat, 'lon':lon, 'ivt':ivt_voss, 'ivt_dir':ivt_dir_voss},f)

