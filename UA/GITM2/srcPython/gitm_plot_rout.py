#!/usr/bin/env python
#-----------------------------------------------------------------------------
# gitm_plot_rout
#
# Author: Angeline G. Burrell, UMichigan, Jan 2013
#
# Comments: Common routine used to make GITM plots.
#
# Includes: add_colorbar                  - add a colorbar to a contour plot
#           center_polar_cap              - center radial coordinates for a
#                                           polar plot
#           find_data_limits              - find the upper and lower limits
#                                           for a list of GITM data arrays
#----------------------------------------------------------------------------

'''
Plot data from a 2D GITM file (or 3D GITM file at a single altitude) for
different geographic configurations
'''

# Import modules

def add_colorbar(contour_handle, zmin, zmax, zinc, orient, scale, name, units):
    '''
    Add a colorbar

    Input: contour_handle = handle to contour plot
           zmin           = minimum z value
           zmax           = maximum z value
           zinc           = z tick incriment (recommend 6)
           orient         = orientation of the colorbar (horizontal or vertical)
           scale          = linear or exponential?
           name           = z variable name
           units          = z variable units
    '''
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FormatStrFormatter

    w  = np.linspace(zmin, zmax, zinc, endpoint=True)
    cb = plt.colorbar(contour_handle, ticks=w, pad=.15, orientation=orient,
                      fraction=.06)
    if(scale is "exponetial"):
        cb.formatter=FormatStrFormatter('%7.2E')
    cb.set_label(r'%s ($%s$)' % (name, units))
    cb.update_ticks()
    return cb

def center_polar_cap(rcenter, redge, r):
    '''
    Adjust the radial axis in a polar plot so that it is centered about
    the northern or southern pole
    '''

    if(rcenter > redge):
        return rcenter - r
    else:
        return r

def find_data_limits(gDataList, xkey, lat_index=-1, lon_index=-1, alt_index=-2,
                     inc=6, *args, **kwargs):
    '''
    Establish the appropriate axis limits for a list of GitmBin files at
    a particular latitude/longitude index

    Input: gDataList = A list of GitmBin data structures
           xkey      = key for the desired values
           lat_index = latitude index (default -1 for no index)
           lon_index = longitude index (default -1 for no index)
           alt_index = altitude index (default -2 for no index, -1 for 2D)
           inc       = number of tick incriments (default is 6)
    '''
    import math

    hold_min = []
    hold_max = []

    for gData in gDataList:
        if(lat_index < 0 and lon_index < 0):
            if(alt_index > -2):
                flat = gData[xkey][:,:,alt_index].reshape(-1)
            else:
                flat = gData[xkey][:,:,:].reshape(-1)
        elif(lat_index < 0):
            if(alt_index > -2):
                flat = gData[xkey][lon_index,:,alt_index].reshape(-1)
            else:
                flat = gData[xkey][lon_index,:,:].reshape(-1)
        elif(lon_index < 0):
            if(alt_index > -2):
                flat = gData[xkey][:,lat_index,alt_index].reshape(-1)
            else:
                flat = gData[xkey][:,lat_index,:].reshape(-1)
        else:
            if(alt_index > -1):
                flat = gData[xkey][lon_index,lat_index,alt_index].reshape(-1)
            else:
                flat = gData[xkey][lon_index,lat_index,:].reshape(-1)

        hold_min.append(min(flat))
        hold_max.append(max(flat))

    xmin = min(hold_min)
    xmax = max(hold_max)
    xran = round((xmax-xmin)/inc)

    if(xran != 0.0):
        xmin = math.floor(float("%.14f" % (xmin / xran))) * xran
        xmax = math.ceil(float("%.14f" % (xmax / xran))) * xran

    # Consider physical limits for Latitude and Longitude keys.

    if(xkey == "dLat"):
        if(xmin < -90.0):
            xmin = -90.0
        if(xmax > 90.0):
            xmax = 90.0
    elif(xkey == "Latitude"):
        if(xmin < -np.pi / 2.0):
            xmin = -np.pi / 2.0
        if(xmax > np.pi / 2.0):
            xmax = np.pi / 2.0
    elif(xkey == "dLon" or xkey == "Longitude"):
        if(xmin < 0.0):
            xmin = 0.0
        if(xkey == "dLon" and xmax > 360.0):
            xmax = 360
        elif(xkey == "Longitude" and xmax > np.pi):
            xmax = np.pi

    return xmin, xmax

def localtime_to_glon(ut_datetime, localtime):
    '''
    Routine to compute the longitude where the local time is at a specified
    value for a specified universal time

    ut_datetime = Universal time as a datetime object
    localtime   = Local time in hours
    '''

    uth = ut_datetime.hour+(ut_datetime.minute/60.0)+(ut_datetime.second/3600.0)
    lon = (localtime - uth) * 15.0 # 15 = 360 degrees / 24 hours

    return lon
#End

