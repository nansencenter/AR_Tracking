'''
latlon_disttools.py
Set of tools for calculating distances and finding nearest grid points
based on given lat lon grid and alt lon coordiantes.
Tools include:
    dist_latlon = distance between two points
    dist_map = grid of distances between all points and specified lat/lon
    dist_map_plot = plot of distance between all points and specified lat/lon
    nearest_point = indices of grid point nearest to specified lat/lon, optionally
                    including whether the grid point is over land or sea

Written by Stephen Outten 2018
'''

import numpy as np

def dist_latlon(lat1,lon1,lat2,lon2):
    '''
    Calculates the distance between two points given their latitudes and londitudes.
    Usage:  distance = dist_latlon(lat1,lon1,lat2,lon2)
        lat1 = latitude of point 1
        lon1 = londitude of point 1
        lat2 = latitude of point 2
        lon2 = londitude of point 2

        distance = distance between the two points in kilometers
    '''
    R = 6373.0
    lat1 = np.radians(lat1)
    lat2 = np.radians(lat2)
    lon1 = np.radians(lon1)
    lon2 = np.radians(lon2)

    dlon = lon2 - lon1
    dlat = lat2 - lat1

    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))

    distance = R * c

    return distance


def dist_map(lat,lon,latgrid,longrid):
    '''
    Calculates the distance between the point specified by lat/lon, and
    every grid point in the grids latgrid and longrid.
    Usage: distances = dist_map(lat,lon,latgrid,longrid)
        lat = latitude of specified point 
        lon = longitude of specified point
        latgrid = array of latitudes
        longrid = array of longitudes

        distances = array of distances in km for each point in latgrid/longrid to point lat/lon
    '''
    if (lat<np.nanmin(latgrid)) or (lat>np.nanmax(latgrid)) or (lon<np.nanmin(longrid)) or (lon>np.nanmax(longrid)):
        print("Grid point {0:.1f}N, {0:.1f}E is outside of the grid of latitudes and longitudes provided.".format(lat,lon))
    a,b = latgrid.shape
    distances = np.nan * np.ones([a,b])
    for i in range(a):
        for j in range(b):
            distances[i,j] = dist_latlon(lat,lon,latgrid[i,j],longrid[i,j])
    return distances


def dist_map_plot(lat,lon,latgrid,longrid,latmin=None,latmax=None,lonmin=None,lonmax=None):
    '''
    Plots a map of distances between the point specified by lat/lon, and
    every grid point in the grids latgrid and longrid.
    Limits of the map may be given or are taken as limits of latgrid and lon grid.
    Uses Basemap if available.
    Usage: dist_map_plot(lat,lon,latgrid,longrid)
        lat = latitude of specified point
        lon = longitude of specified point
        latgrid = array of latitudes
        longrid = array of longitudes
        latmin = lowest latitude for plot
        latmax = highest latitude of plot
        lonmin = lowest longitude of plot
        lonmax = highest longitude of plot
    '''
    if latmin==None: latmin = np.nanmin(latgrid)
    if latmax==None: latmax = np.nanmax(latgrid)
    if lonmin==None: lonmin = np.nanmin(longrid)
    if lonmax==None: lonmax = np.nanmax(longrid)
    bmap = False
    try:    
        from mpl_toolkits.basemap import Basemap
        bmap = True
    except:
        print("Basemap not available, unable to plot orography.")
    try:
        import matplotlib.pyplot as plt
    except:
        print("Matplotlib not available. Unable to plot.")
        return

    distances = dist_map(lat, lon, latgrid, longrid)
    pltmaxdist =  int(np.nanmax(distances)/100)*100
    f1, axs = plt.subplots(1, 1)
    if bmap: 
        m = Basemap(projection='mill', llcrnrlat=latmin,urcrnrlat=latmax,llcrnrlon=lonmin,urcrnrlon=lonmax,resolution='l',ax=axs,fix_aspect=False)
        m.drawcoastlines()
        x, y = m(longrid, latgrid)
        cs = m.contourf(x,y,distances,30,vmin=0,vmax=pltmaxdist)
        m.drawparallels(np.arange(-90,90,10),labels=[True,False,False,False], linewidth=1)
        m.drawmeridians(np.arange(-180,180,10),labels=[False,False,False,True], linewidth=1)
        f1.colorbar(cs)
        plt.show()
        return
    else:
        cs = axs.contourf(longrid,latgrid,distances,30,vmin=0,vmax=pltmaxdist)
        axs.contour(longrid,latgrid,latgrid,np.arange(-90,90,10),colors='k')
        axs.contour(longrid,latgrid,longrid,np.arange(-180,180,10),colors='k')
        f1.colorbar(cs)
        plt.show()
        return


def nearest_point(lat,lon,latgrid,longrid,lsm=None,surface=None):
    '''
    Returns the indicies in the array of latgrid/longrid for the grid point closest 
    to the location specified by lat/lon.
    Usage: indx, indy = nearest_point(lat,lon,latgrid,longrid)
        lat = latitude of specified point
        lon = longitude of specified point
        latgrid = array of latitudes
        longrid = array of longitudes
        lsm = array of land sea mask (0s and 1s) or orography
        surface = "land" or "sea"   
        indx, indy = x and y index for for nearest grid point in gridlat/gridlon
    '''
    if (lat<np.nanmin(latgrid)) or (lat>np.nanmax(latgrid)) or (lon<np.nanmin(longrid)) or (lon>np.nanmax(longrid)):
        print("Grid point {0:.1f}N, {0:.1f}E is outside of the grid of latitudes and longitudes provided.".format(lat,lon))
        return None,None
    if str(surface).lower() not in ["land","sea","none"]:
        print('Surface must be "land" or "sea')
        return None, None

    a,b = latgrid.shape
    distances = np.nan * np.ones([a,b])
    for i in range(a):
        for j in range(b):
            distances[i,j] = dist_latlon(lat,lon,latgrid[i,j],longrid[i,j])

    lsm = np.array(lsm)    # avoids modifying elements of mutable object even if not returned
    if (str(lsm).lower() != 'none') and (str(surface).lower() != 'none'):
        lsm[lsm<0] = 0   # if sea have negative orography
        lsm[lsm>0] = 1     # if land is orography height
        lsm[np.isnan(lsm)] = 1   # if land is give as NaNs
        if str(surface).lower() == 'land':
            lsm[lsm==0] = np.nan
        else:
            lsm[lsm==1] = np.nan
            lsm[lsm==0] = 1
        distances *= lsm

    indx = np.where(distances==np.nanmin(distances))[0][0]
    indy = np.where(distances==np.nanmin(distances))[1][0]

    return indx, indy


