# AR_Tracking
Atmospheric tracking algorithm for ERA5 based on Brands 2017 and Lavers 2012.  
This code was developed as part of the MUDYFEET project, run by the Bjerknes Centre for Climate Research (BCCR).
Details of the algorithm and examples of the results are given the the NERSC technical report number 399 provided in this repository.

This code identifies atmospheric rivers (AR) based on the integrated vapour transport (IVT) magnitude and direction.  
It was designed to work with ERA5 reanalysis to identify ARs in the North Atlantic approaching southern Norway.  
The code requires the zonal and merridional components of IVT to be stored in separate directories (e.g. IVTU/ and IVTV/), with one year of hourly data per file in NetCDF format.   
These files do not need to be global but must cover some domain that includes some of the North Atlantic and Southern Norway.  
Originally, this was used with a domain covering 20N to 80N and -80 to 40E.  
The IVTU and IVTV data can be downloaded on any chosen domain from the ECMWF web site.  
The code will combine the two components to make a single file containing IVT magnitude and direction, which it stores in an initially empty directory called IVT/.  

There are seven codes as follows:  
**calc_ivt.py**  
Converts IVTU and IVTV into files containing IVT magnitude and direction. These will be in Dill Pickle format.

**calc_percentiles.py**  
Calculates the monthly percentiles for every gripoint. This is used by the tracking algorithm for the identification thresholds. This only needs to be run once but it requires a lot of memory.

**ar_tracking.py**  
Main tracking algorithm. Uses triggering and tracking thresholds to identify the track through an AR at every time interval.
Output is NetCDF containing time of occurence of AR, length of AR, longitudinal extent of AR, IVT along AR track, and lat/lon positions of each point along track. 

**reduce_ar_persist.py**  
Reduces the identified ARs based on their persistence in time. Different literature uses different thresholds for how long an AR must persist for it to be considerd an AR. Same output as for ar_tracking.py but with reduced numbers of AR.

**extract_ar.py**  
While the reduce_AR_persist.py provides information on when and where ARs occur, this code extracts selected variables from ERA5 for those times over a specified domain, e.g. MSLP, precip, etc. These variables must be stored in a directory in advance, i.e. downlaoded from the ECMWF website. 
This code produces NetCDFs containing a subset of the ERA5 data for the selected variables limited to the specified domain and only at times when ARs are occuring. 

**identify_voss.py**  
Loads the identified ARs from reduce_AR_persist.py, isolates the times an AR was related to the Voss flooding event (26th to 29th October 2014), and identifies all other ARs that are similar in length, magnitude, and position. 
The identified ARs are tabulated. displayed, and stored.

**make_voss_sample.py**  
This code extracts the IVT magnitude and direction for the period of the Voss flooding event (26th to 29th October 2014).
Used for plotting the Voss flooding event AR.


  
