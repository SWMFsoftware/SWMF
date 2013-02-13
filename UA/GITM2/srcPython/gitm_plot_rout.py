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
#           find_zlimits                  - find the upper and lower limits
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

def find_zlimits(gDataList, zkey, aindex=-1, zinc=6, *args, **kwargs):
    '''
    Establish the appropriate z-axis limits for a list of GitmBin files

    Input: gDataList = A list of GitmBin data structures
           zkey      = key for the desired z value
           aindex    = altitude index (default -1 for 2D measurement,
                       use -2 for no index)
           zinc      = number of tick incriments (default is 6)
    '''
    import math

    hold_min = []
    hold_max = []

    for gData in gDataList:
        if(aindex > -2):
            flat = gData[zkey][:,:,aindex].reshape(-1)
        else:
            flat = gData[zkey][:,:,:].reshape(-1)

        hold_min.append(min(flat))
        hold_max.append(max(flat))

    zmin = min(hold_min)
    zmax = max(hold_max)
    zran = round((zmax-zmin)/zinc)

    if(zran != 0.0):
        zmin = math.floor(zmin / zran) * zran
        zmax = math.ceil(zmax / zran) * zran

    return zmin, zmax

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

