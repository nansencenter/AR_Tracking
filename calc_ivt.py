'''
calc_ivt.py
Loads the ERA5 IVTU and IVTV and converts them to magnitude and direction.
Assumes IVTU and IVTV data are in subdirectories called ITVU and IVTV repsectively, and 
the files are grouped by years with filenames of the format IVTU_year.nc (e.g. IVTU_1979).
The code also assumes there are files from 1979 to 2018 inclusive. 
The output is stored in a subdirectory called IVT and is in the format of a Dill Pickle.
The output consists of magnitude of IVT and direction of IVTThe output consists of magnitude of IVT and direction of IVT.

Written by Stephen Outten August 2019
'''

from netCDF4 import Dataset
import numpy as np
import dill

# Directory where all subdirectories are located.
FolderName = '/Data/ERA5/Testing/NorthHemis/'

# Load lon and lat from any file
fn = FolderName + 'IVTU/' + 'IVTU_1979.nc'
dataset = Dataset(fn)
lon = np.array(dataset.variables['longitude'][:])
lat = np.array(dataset.variables['latitude'][:])
dataset.close()


# Loop over files and calculate accurate magnitudes and directions
# Storing the output as dill pickles
for year in np.arange(1979,2019):
    print(year)

    # Load the IVTU for a compelte file
    fn = FolderName + 'IVTU/IVTU_' + str(year) + '.nc'
    dataset = Dataset(fn)
    ivtu = np.array(dataset.variables['p71.162'][:])
    dataset.close()

    # Load the IVTV for a compelte file
    fn = FolderName + 'IVTV/IVTV_' + str(year) + '.nc'
    dataset = Dataset(fn)
    ivtv = np.array(dataset.variables['p72.162'][:])
    dataset.close()

    # Calcualte magnitude and direction
    ivt = np.sqrt(ivtu**2 + ivtv**2)
    ivt_dir = 180 + (180/np.pi) * np.arctan2(ivtu, ivtv)
    ivt = ivt.astype('float32')
    ivt_dir = ivt_dir.astype('float32')
    del ivtu, ivtv

    # Store lon, lat, ivt, and ivt_dir in dill pickle
    fnout = FolderName + 'IVT/IVT_' + str(year) + '.dpkl'
    with open(fnout, 'wb') as f:
        dill.dump({'lat':lat, 'lon':lon, 'ivt':ivt, 'ivt_dir':ivt_dir},f)

    del ivt, ivt_dir







