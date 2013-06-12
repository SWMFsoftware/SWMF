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
#           localtime_to_glon             - convert local time to longitude
#           find_lat_lon_index            - find the appropriate index for
#                                           a specified location
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

def find_lon_lat_index(gData, glon, glat, units="degrees"):
    '''
    Routine to locate the appropriate longitude and latitude indexes for a
    given location.  The location may be specified in degrees (default) or
    radians.
    '''

    import string

    # Set the keys to look for the location in the appropriate units
    if string.lower(units) == "degrees":
        latkey = "dLat"
        lonkey = "dLon"
    else:
        latkey = "Latitude"
        lonkey = "Longitude"

    # First identify the appropriate Longitude.  All longitude values are
    # the same for any latitude and altitude index.


    for (lonindex,clon) in enumerate(gData[lonkey][:,0,0]):
        if clon >= glon:
            if (clon - glon) > (glon - gData[lonkey][lonindex-1,0,0]):
                lonindex = lonindex - 1
            break
        
    # Next identify the appropriate Latitude at the specified longitude

    for (latindex,clat) in enumerate(gData[latkey][lonindex,:,0]):
        if clat >= glat:
            if (clat - glat) > (glat - gData[latkey][lonindex,latindex-1,0]):
                latindex = latindex - 1
            break

    return(lonindex, latindex)

def retrieve_key_from_web_name(name):
    '''
    A routine to retrieve a GITM key corresponding to a descriptive name.
    '''
    key_dict = {"Altitude":"Altitude", "Argon Mixing Ratio":"Ar Mixing Ratio",
                "[Ar]":"Ar","Methane Mixing Ratio":"CH4 Mixing Ratio",
                "Conduction":"Conduction", "EUV Heating":"EuvHeating",
                "[H]":"H", "[H$^+$]":"H!U+!N", "[He]":"He",
                "H$_2$ Mixing Ratio":"H2 Mixing Ratio", "[He$^+$]":"He!U+!N",
                "Hydrogen Cyanide Mixing Ratio":"HCN Mixing Ratio",
                "Heating Efficiency":"Heating Efficiency",
                "Heat Balance Total":"Heat Balance Total",
                "Latitude (rad)":"Latitude", "Longitude (rad)":"Longitude",
                "[N$_2$]":"N!D2!N", "[N$_2$$^+$]":"N!D2!U+!N",
                "[N$^+$]":"N!U+!N", "[N($^2$D)]":"N(!U2!ND)",
                "[N($^2$P)]":"N(!U2!NP)", "[N($^4$S)]":"N(!U4!NS)",
                "N$_2$ Mixing Ratio":"N2 Mixing Ratio", "[NO]":"NO",
                "[NO$^+$]":"NO!U+!N", "[O($^4$SP)$^+$]":"O_4SP_!U+!N",
                "[O($^1$D)]":"O(!U1!ND)", "[O$_2$$^+$]":"O!D2!U+!N",
                "[O($^2$D)]":"O(!U2!ND)!", "[O($^2$D)]":"O(!U2!ND)!U+!N",
                "[O($^2$P)$^+$]":"O(!U2!NP)!U+!N", "[O$_2$]":"O!D2!N",
                "[O($^2$P)]":"O(!U2!NP)!U+!N", "[O($^3$P)]":"O(!U3!NP)",
                "Radiative Cooling":"RadCooling", "Neutral Density":"Rho",
                "T$_n$":"Temperature", "v$_{East}$":"V!Di!N (east)",
                "v$_{North}$":"V!Di!N (north)", "v$_{Up}$":"V!Di!N (up)",
                "u$_{East}$":"V!Dn!N (east)","u$_{North}$":"V!Dn!N (north)",
                "u$_{Up}$":"V!Dn!N (up)", "[e-]":"e-",
                "u$_{Up, N_2}$":"V!Dn!N (up,N!D2!N              )",
                "u$_{Up, N(^4S)}$":"V!Dn!N (up,N(!U4!NS)           )",
                "u$_{Up, NO}$":"V!Dn!N (up,NO                  )",
                "u$_{Up, O_2}$":"V!Dn!N (up,O!D2!N              )",
                "u$_{Up, O(^3P)}$":"V!Dn!N (up,O(!U3!NP)           )",
                "Electron Average Energy":"Electron_Average_Energy",
                "T$_e$":"eTemperature", "T$_i$":"iTemperature",
                "Solar Zenith Angle":"Solar Zenith Angle", "[CO$_2$]":"CO!D2!N",
                "Vertical TEC":"Vertical TEC", "DivJu FL":"DivJu FL",
                "DivJuAlt":"DivJuAlt", "Field Line Length":"FL Length",
                "Electron Energy Flux":"Electron_Energy_Flux",
                "$\sigma_P$":"Pedersen FL Conductance", "Potential":"Potential",
                "$\Sigma_P$":"Pedersen Conductance", "Region 2 Current":"Je2",
                "$\sigma_H$":"Hall FL Conductance", "Region 1 Current":"Je1",
                "Hall Conductance":"$\Sigma_H$", "Ed1":"Ed1", "Ed2":"Ed2",
                "Vertical Electric Field":"E.F. Vertical",
                "Eastward Electric Field":"E.F. East", "dLat":"Latitude (deg)",
                "Northward Electric Field":"E.F. North",
                "Electric Field Magnitude":"E.F. Magnitude",
                "Magnetic Latitude":"Magnetic Latitude",
                "Magnetic Longitude":"Magnetic Longitude",
                "dLon":"Longitude (deg)", "LT":"Solar Local Time"}

    if key_dict.has_key(name):
        return key_dict[name]
    else:
        print "ERROR: unknown data type [", name, "], known names are: "
        print key_dict.keys()

#End

