#!/usr/bin/env python
#-----------------------------------------------------------------------------
# gitm_plot_rout
#
# Author: Angeline G. Burrell, UMichigan, Jan 2013
#
# Comments: Common routine used to make GITM plots.
#
# Includes: choose_contour_map         - choose a color map depending on certain
#                                        specified plot characteristics
#           add_colorbar               - add a colorbar to a contour plot
#           find_order_of_magnitude    - find the order of magnitude
#           center_polar_cap           - center radial coordinates for a
#                                        polar plot
#           find_data_limits           - find the upper and lower limits
#                                        for a list of GITM data arrays
#           find_data_limits_irange    - find the upper and lower limits
#                                        for a list of GITM data arrays and
#                                        a range of lon/lat/alt indexes
#           localtime_to_glon          - convert local time to longitude
#           find_lon_lat_index         - find the appropriate index for
#                                        a specified location
#           retrieve_key_from_web_name - a routine to retrieve the data key
#                                        from a website-friendly data name
#           find_alt_index             - A routine to find the appropriate index
#                                        for a specified altitude
#           match_cindi_key            - a routine to retrieve the CINDI data
#                                        key from a GITM key
#----------------------------------------------------------------------------

'''
Plot data from a 2D GITM file (or 3D GITM file at a single altitude) for
different geographic configurations
'''

# Import modules

def choose_contour_map(color, center):
    '''
    Choose a standard contour color map based on whether the output image
    will be black and white or color, centered around zero or not.

    Input: color  = True for color, False for black and white
           center = True for centered about zero, False if not
    '''

    if color:
        if center:
            return("seismic_r")
        else:
            return("Spectral_r")
    else:
        if center:
            return("binary")
        else:
            return("Greys")

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
    import math
    from matplotlib.ticker import FormatStrFormatter, FuncFormatter

    w  = np.linspace(zmin, zmax, zinc, endpoint=True)
    cb = plt.colorbar(contour_handle, ticks=w, pad=.15, orientation=orient,
                      fraction=.07)

    if scale.find("exponential") >= 0:
        omag = find_order_of_magnitude(zmax)

        def scaled_ticks(x, pos):
            '''
            Define ticks so that they are scaled by the order of magnitude.
            The two arguements are the value (x) and the tick position (pos)
            and are required for this function to be used by FuncFormatter
            '''
            tckstr = "{:.1f}".format(x / math.pow(10.0, omag))
            return tckstr

        # Use the previously defined function to scale the ticks
        cb.formatter = FuncFormatter(scaled_ticks)
        # Set the label
        cb.set_label(r'{:s} (${:s} \times 10^{{{:.0f}}}$)'.format(name, units,
                                                                  omag))
    else:
        zscale = max(abs(zmin), abs(zmax))
        if zscale > 1.0e3 or zscale < 1.0e-3:
            cb.formatter=FormatStrFormatter('{:.2g}')
        elif zscale < 1.0e1:
            cb.formatter=FormatStrFormatter('{:.2f}')
        else:
            cb.formatter=FormatStrFormatter('{:.0f}')
        # Set the label
        cb.set_label(r'{:s} (${:s}$)'.format(name, units))

    # Update the ticks to reflect formatting
    cb.update_ticks()
    return cb

def find_order_of_magnitude(value):
    '''
    Find the order of magnitude of a number.  Returns the exponent.
    Ex: -4000.0 = -4 x 10^3 will return 3
    '''
    import math

    return math.floor(math.log10(abs(value)))

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

def find_data_limits_irange(gDataList, xkey, min_ilat=-1, max_ilat = -1,
                            min_ilon=-1, max_ilon=-1, min_ialt=-2, max_ialt=-2,
                            inc=6, *args, **kwargs):
    '''
    Establish the appropriate axis limits for a list of GitmBin files at
    a specified range of indexes.  If you only want one index, the maximum
    index should have a value on integer higher than the minimum index, which
    contains the desired index integer.

    Input: gDataList = A list of GitmBin data structures
           xkey      = key for the desired values
           min_ilat  = minimum latitude index (default -1 for no index)
           max_ilat  = maximum latitude index (default -1 for no index)
           min_ilon  = minimum longitude index (default -1 for no index)
           max_ilon  = maximum llongitude index (default -1 for no index)
           min_ialt  = minimum altitude index (default -2 for no index, -1
                       for 2D)
           max_ialt  = maximum laltitude index (default -2 for no index, -1
                       for 2D)
           inc       = number of tick incriments (default is 6)
    '''
    import math

    hold_min = []
    hold_max = []

    for gData in gDataList:
        if(min_ilat < 0 and min_ilon < 0):
            if(min_ialt > -2):
                flat = gData[xkey][:,:,min_ialt:max_ialt].reshape(-1)
            else:
                flat = gData[xkey][:,:,:].reshape(-1)
        elif(min_ilat < 0):
            if(min_ialt > -2):
                flat = gData[xkey][min_ilon:max_ilon,:,
                                   min_ialt:max_ialt].reshape(-1)
            else:
                flat = gData[xkey][min_ilon:max_ilon,:,:].reshape(-1)
        elif(min_ilon < 0):
            if(min_ialt > -2):
                flat = gData[xkey][:,min_ilat:max_ilat,
                                     min_ialt:max_ialt].reshape(-1)
            else:
                flat = gData[xkey][:,min_ilat:max_ilat,:].reshape(-1)
        else:
            if(min_ialt > -1):
                flat = gData[xkey][min_ilon:max_ilon,min_ilat:max_ilat,
                                   min_ialt:max_ialt].reshape(-1)
            else:
                flat = gData[xkey][min_ilon:max_ilon,min_ilat:max_ilat,
                                   :].reshape(-1)

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

def find_data_limits_ivalues(gDataList, xkey, lat_indices, lon_indices,
                             alt_indices, lat_range=False, lon_range=False, 
                             alt_range=False, inc=6, rvals=True, *args,
                             **kwargs):
    '''
    Establish the appropriate axis limits for a list of GitmBin files at
    a specified list of indexes.  If you want to include the entire range of
    latitudes, longitudes, or altitudes instead of a particular index, set
    the range flag to True.  If this is set to True, the first value in the
    xxx_indices list will be used as the lower index and the second value as
    the upper index.  If you are using indexes for more than one coordinate,
    keep in mind that the indices will be paired unless they do not have the
    same length.

    Example:

    xkey = "e-"
    lat_indices = [1, 3, 5] with lat_range = False
    lon_indices = [0, 2] with lon_range = False
    alt_indices = [0, gdata.attrs['nAlt']] with alt_range = True

    Finds the maximum and minimum electron density over all altitudes at
    (lon index, lat index) locations: (0,1), (2,3), and (2,5)

    Input: gDataList   = A list of GitmBin data structures
           xkey        = key for the desired values
           lat_indices = list of latitude indices (or lower index, upper index)
           lon_indices = list of longitude indices (or lower index, upper index)
           alt_indices = list of altitude indices (or lower index, upper index)
           lat_range   = use range for lat instead of indexes (default is False)
           lat_range   = use range for lon instead of indexes (default is False)
           lat_range   = use range for alt instead of indexes (default is False)
           inc         = number of tick incriments (default is 6)
           rvals       = Round the min and max to a sensible number based
                         on the desired number of incriments?  Will never
                         decrease the maximum or increase the minimum.  Will
                         also assess the max/min for latitude and longitudes
                         to ensure they fall within physically sensible limits.
                         (default is True)

    Output: xmin = miminum data value
            xmax = maximum data value
    '''
    import math
    import sys
    import copy

    # Initialize the variables

    hold_min = []
    hold_max = []
    ilat = -1
    ilon = -1
    ialt = -1
    latlist = copy.deepcopy(lat_indices)
    lonlist = copy.deepcopy(lon_indices)
    altlist = copy.deepcopy(alt_indices)

    # Set the maximum and minimum values for the range coordinates

    if lat_range:
        ilatmax = latlist.pop()
        ilatmin = latlist.pop()

    if lon_range:
        ilonmax = lonlist.pop()
        ilonmin = lonlist.pop()

    if alt_range:
        ialtmax = altlist.pop()
        ialtmin = altlist.pop()

    # Cycle over the index coordinates
    while len(altlist) > 0 or len(lonlist) > 0 or len(latlist) > 0:
        # Initialize the index coordiates
        if len(latlist) > 0:
            ilat = latlist.pop()

        if len(lonlist) > 0:
            ilon = lonlist.pop()

        if len(altlist) > 0:
            ialt = altlist.pop()

        if((not lat_range and ilat == -1) or (not lon_range and ilon == -1)
           or (not alt_range and ialt == -1)):
            print "INPUT ERROR in find_data_limit_ivalues"
            sys.exit(1)

        # Cycle over the GITM data structures
        for gData in gDataList:
            # Make the appropriate data selection
            if lat_range:
                if lon_range:
                    if alt_range:
                        # This is silly to do, but sometimes people are silly
                        flat = gData[xkey][ilonmin:ilonmax,ilatmin:ilatmax,
                                           ialtmin:ialtmax].reshape(-1)
                    else:
                        flat = gData[xkey][ilonmin:ilonmax,ilatmin:ilatmax,
                                           ialt].reshape(-1)
                else:
                    if alt_range:
                        flat = gData[xkey][ilon,ilatmin:ilatmax,
                                           ialtmin:ialtmax].reshape(-1)
                    else:
                        flat = gData[xkey][ilon,ilatmin:ilatmax,
                                           ialt].reshape(-1)
            else:
                if lon_range:
                    if alt_range:
                        flat = gData[xkey][ilonmin:ilonmax,ilat,
                                           ialtmin:ialtmax].reshape(-1)
                    else:
                        flat = gData[xkey][ilonmin:ilonmax,ilat,
                                           ialt].reshape(-1)
                else:
                    if alt_range:
                        flat = gData[xkey][ilon,ilat,
                                           ialtmin:ialtmax].reshape(-1)
                    else:
                        # This is silly to do, but sometimes people are silly
                        flat = gData[xkey][ilon,ilat,ialt].reshape(-1)
        
            # Maintain the maximum and mimimum values for this data selection
            hold_min.append(min(flat))
            hold_max.append(max(flat))

    # Find the data max and min from the individual selection max and mins
    xmin = min(hold_min)
    xmax = max(hold_max)
    xran = round((xmax-xmin)/inc)

    if rvals:
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

def glon_to_localtime(ut_datetime, glon, units="degrees"):
    '''
    Routine to compute the local time where the longitude is at a specified
    value for a specified universal time

    ut_datetime = Universal time as a datetime object
    glon        = Longitude
    units       = units of longitude (default degrees)
    '''

    scale = 1.0 / 15.0 # (hours / degree)
    if units.find("rad") >= 0:
        scale *= (180.0 / np.pi) # Convert to hours / radians
    uth = ut_datetime.hour+(ut_datetime.minute/60.0)+(ut_datetime.second/3600.0)
    lt = glon * scale + uth

    # Ensure that the local time falls between 00:00 and 23:59
    if lt < 0.0:
        lt += 24.0

    if lt >= 24.0:
        lt -= 24.0

    return lt

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

    Input:
    gData = GitmBin data structure
    glon  = longitude
    glat  = latitude
    units = units of longitude and latitude (default = degrees)
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
                "[H]":"H", "[H+]":"H!U+!N", "[He]":"He",
                "H2 Mixing Ratio":"H2 Mixing Ratio", "[He+]":"He!U+!N",
                "Hydrogen Cyanide Mixing Ratio":"HCN Mixing Ratio",
                "Heating Efficiency":"Heating Efficiency",
                "Heat Balance Total":"Heat Balance Total",
                "Latitude (rad)":"Latitude", "Longitude (rad)":"Longitude",
                "[N2]":"N!D2!N", "[N2+]":"N!D2!U+!N",
                "[N+]":"N!U+!N", "[N(2D)]":"N(!U2!ND)",
                "[N(2P)]":"N(!U2!NP)", "[N(4S)]":"N(!U4!NS)",
                "N2 Mixing Ratio":"N2 Mixing Ratio", "[NO]":"NO",
                "[NO+]":"NO!U+!N", "[O(4SP)+]":"O_4SP_!U+!N",
                "[O(1D)]":"O(!U1!ND)", "[O2+]":"O!D2!U+!N",
                "[O(2D)]":"O(!U2!ND)!", "[O(2D)]":"O(!U2!ND)!U+!N",
                "[O(2P)+]":"O(!U2!NP)!U+!N", "[O2]":"O!D2!N",
                "[O(2P)]":"O(!U2!NP)!U+!N", "[O(3P)]":"O(!U3!NP)",
                "Radiative Cooling":"RadCooling", "Neutral Density":"Rho",
                "Tn":"Temperature", "v(East)":"V!Di!N (east)",
                "v(North)":"V!Di!N (north)", "v(Up)":"V!Di!N (up)",
                "u(East)":"V!Dn!N (east)","u(North)":"V!Dn!N (north)",
                "u(Up)":"V!Dn!N (up)", "[e-]":"e-",
                "u(Up, N_2)":"V!Dn!N (up,N!D2!N              )",
                "u(Up, N(4S))":"V!Dn!N (up,N(!U4!NS)           )",
                "u(Up, NO)":"V!Dn!N (up,NO                  )",
                "u(Up, O_2)":"V!Dn!N (up,O!D2!N              )",
                "u(Up, O(3P))":"V!Dn!N (up,O(!U3!NP)           )",
                "Electron Average Energy":"Electron_Average_Energy",
                "Te":"eTemperature", "Ti":"iTemperature",
                "Solar Zenith Angle":"Solar Zenith Angle", "[CO2]":"CO!D2!N",
                "Vertical TEC":"Vertical TEC", "DivJu FL":"DivJu FL",
                "DivJuAlt":"DivJuAlt", "Field Line Length":"FL Length",
                "Electron Energy Flux":"Electron_Energy_Flux",
                "sigmaP":"Pedersen FL Conductance", "Potential":"Potential",
                "SigmaP":"Pedersen Conductance", "Region 2 Current":"Je2",
                "sigmaH":"Hall FL Conductance", "Region 1 Current":"Je1",
                "Hall Conductance":"SigmaH", "Ed1":"Ed1", "Ed2":"Ed2",
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

def find_alt_index(gData, ilon, ilat, alt, units="km"):
    '''
    Routine to locate the appropriate altitude index in a given array.
    The altitude may be specified in km (default) or m.
    '''

    import string
    import numpy as np

    if string.lower(units) == "km":
        alt *= 1000.0

    # The GITM arrays are sorted, so we can use searchsorted to find
    # a close index

    ialt = np.searchsorted(gData['Altitude'][ilon,ilat,:], alt)

    # If the search index is zero, it is at or below the minimum altitude
    # and doesn't need to be changed.  If the distance between the desired
    # altitude is closer to the search index than the the previous index,
    # it does not need to be changed either

    if(ialt >= gData.attrs['nAlt'] or
       (ialt> 0 and abs(alt - gData['Altitude'][ilon,ilat,ialt]) >
        abs(alt + gData['Altitude'][ilon,ilat,ialt-1]))):
        # If this location is above the maximum height, return maximum index 
        # by reducing the count by one, or if this location is closer to the
        # previous altitude than the one at this index, reduce the count by one
        ialt -= 1
        
    return(ialt)
#End find_alt_index

def match_cindi_key(in_key, out_type="CINDI"):
    '''
    A routine to retrieve a CINDI/GITM data structure key from a GITM/CINDI key.

    in_key   = Input key (GITM or CINDI) that you want to pair.  May specify ALL
               to retrieve a list of all possible keys.
    out_type = Should the output key(s) be the CINDI (default) or GITM keys?
    '''
    key_dict = {"Altitude":"Alt.", "H!U+!N":"FrH", "He!U+!N":"FrHe",
                "NO!U+!N":"FrNO", "O_4SP_!U+!N":"FracO",
                "O(!U2!NP)!U+!N":"FracO", "e-":"Ni(cm^-3)", "iTemperature":"Ti",
                "dLat":"GLAT", "Magnetic Latitude":"MLAT", "dLon":"GLONG",
                "LT":"SLT", "V!Di!N (zon)":"Vzonal", "V!Di!N (par)":"Vpara",
                "V!Di!N (mer)":"Vmerid"}

    if in_key == "ALL":
        # Return a list of all keys of the desired type
        if out_type == "CINDI":
            return(key_dict.values())
        else:
            return(key_dict.keys())
    else:
        if out_type == "CINDI" and key_dict.has_key(in_key):
            return key_dict[in_key]
        elif out_type == "GITM":
            out_list = [k for k, v in key_dict.iteritems() if v == 'Alt.']
            if len(out_list) == 1:
                return(out_list[0])
            else:
                print "WARNING: unknown CINDI data type [",in_key,"]"
                return None
        elif out_type == "CINDI":
            print "WARNING: unknown GITM data type [",in_key,"], known names are:"
            print key_dict.keys()
            return None
        else:
            print "WARNING: unknown data source [", out_type, "]"
            return None
# END match_cindi_key
