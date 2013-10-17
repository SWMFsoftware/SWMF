#!/usr/bin/env python
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission
#  For more information, see http://csem.engin.umich.edu/tools/swmf
#-----------------------------------------------------------------------------
# $Id$
# gitm_loc_rout
#
# Author: Angeline G. Burrell, UMichigan, Jan 2013
#
# Comments: Common routine used to find GITM data at a specific location(s)
#           and match with observations.
#
# Includes: find_nearest_location - A routine to find the nearest neighbor in
#                                   a 1, 2, or 3D coordinate system
#           find_nearest_value - A routine to find the nearest neighbor to
#                                a specified value
#           find_nearest_datetime - A routine to find the nearest neighbor
#                                   between datetime objects
#           match_running_average - Provide running averages at specified times
#                                   and locations
#           match_running_median - Provide running medians at specified times
#                                  and locations
#           gitm_inst_loc - A routine to align GITM and instrument data.
#           gitm_net_loc - A routine to align GITM and data from a large network
#                          of instruments
#----------------------------------------------------------------------------

'''
Routines to support GITM-Data location mapping and interpolation
'''

# Import modules
import numpy as np
import math
import string
import datetime as dt
from datetime import timedelta

module_name = "gitm_loc_rout"

def find_nearest_location(x_list, y_list, z_list, location, coord="cart"):
    '''
    Routine to find the nearest neighbor to a specified 1, 2, or 3D location.
    ----------- x, y, and z should all be in appropriate units --------------
    CARTESIAN or 1D: Any units may be used, as long as all values have the
                     same units
    CYLINDRICAL or POLAR: x=r, y=theta, z=z; x (and z) should have the same
                          units, y must be in radians
    SPHERICAL: x=r, y=theta, z=phi; x may be in any units, y and z must be in
               radians
    
    Input:
    x_list = list of x values
    y_list = list of y values (for 1D, this should be empty)
    z_list = list of z values (for 1D or 2D, this should be empty)
    location = list of desired coordinates to match (x, y, z), (x, y), or (x)
    coord = type of coordinate system.  Note if 1D is specified, this flag
            is not used.  If 2D is specified, both cylindrical and spherical
            coordinate systems will default to polar.
            ('cart'=Cartesian (default), 'cyl'=Cylindrical, 'sph'=Spherical)

    Output:
    delta = distance between the desired and nearest neighbor location
    index = index of the nearest paired x, y and/or z values in the input lists
    '''
    rout_name = string.join([module_name, "find_nearest_location"], sep=" ")

    # Test the length of the input arrays
    if len(x_list) < 1:
        print rout_name, "ERROR: x-coordinate list is empty"
        return

    if len(x_list) != len(y_list) and len(y_list) > 0:
        print rout_name, "ERROR: x and y coordinate lists are mismatched"
        return

    if len(x_list) != len(z_list) and len(z_list) > 0:
        print rout_name, "ERROR: x and z coordinate lists are mismatched"
        return

    if len(location) == 3 and len(z_list) <= 0:
        print rout_name, "ERROR: z coordinate list is not specified"
        return

    if len(location) > 1 and len(y_list) <= 0:
        print rout_name, "ERROR: y coordinate list is not specified"
        return

    # Calculate the distance between each data point and the desired location
    # using the pythagorian theorem.
    r_list = list()
    for i,x in enumerate(x_list):
        # Use simple differences for 1D or cartesian systems
        if coord.find("cart") >= 0 or len(location) == 1:
            c2 = (x-location[0])**2
            if len(location) > 1:
                c2 += (y_list[i]-location[1])**2
            if len(location) > 2:
                c2 += (z_list[i]-location[2])**2
        #  Cylindrical or 2D Non-Cartesian coordinates
        elif coord.find("cyl") >= 0 or len(location) == 2:
            c2 = ((x*math.cos(y_list[i])-location[0]*math.cos(location[1]))**2
                  +(x*math.sin(y_list[i])-location[0]*math.sin(location[1]))**2)
            if len(location) > 2:
                c2 += (z_list[i]-location[2])**2
        elif coord.find("sph") >= 0 and len(location) == 3:
            c2 = ((x*math.cos(y_list[i])*math.sin(z_list[i]) -
                   location[0]*math.cos(location[1])*math.sin(location[2]))**2
                  +(x*math.sin(y_list[i])*math.sin(z_list[i]) -
                    location[0]*math.sin(location[1])*math.sin(location[2]))**2
                  +(x*math.cos(z_list[i])-location[0]*math.cos(location[2]))**2)
        else:
            print rout_name, "ERROR: didn't catch problem with input"

        r_list.append(math.sqrt(c2))

    # Find the index of the nearest location by finding the minimum distance
    # between the desired location and the list of coordinates
    delta = np.nanmin(r_list)
    index = r_list.index(delta)

    return delta, index

# End find_nearest_location

def find_nearest_value(data_list, value):
    '''
    Routine to find the nearest neighbor to a specified value

    Input:
    data_list = list to search
    value     = value to match

    Output:
    delta = difference between the nearest value in the data list and the
            desired value (eg data_list[index] - value)
    index = index of the nearest value in the data list
    '''

    # Find the index of the nearest observation time
    index = np.searchsorted(data_list, value)

    # If the value is greater than the maximum in itime, reduce the index
    # by one.  Otherwise test to see whether this index or the one
    # below it is closest
    if index >= len(data_list):
        index = len(data_list) - 1
    elif index > 0:
        if abs(data_list[index] - value) > abs(data_list[index-1] - value):
            index -= 1

    # Find the distance between the GITM and Observational times
    delta = abs(data_list[index] - value)

    return delta, index
# End find_nearest_value

def find_nearest_datetime(dt_list, value):
    '''
    Routine to find the nearest neighbor to a specified value

    Input:
    dt_list = list of datetime objects to search
    value   = datetime value to match

    Output:
    delta = difference between the nearest datetime in the data list and the
            desired datetime (eg dt_list[index] - value)
    index = index of the nearest datetime in the data list
    '''
    # Find the index of the nearest observation time
    index = np.searchsorted(dt_list, value)
    delta = np.nan

    # If the value is greater than the maximum in itime, reduce the index
    # by one.  Otherwise test to see whether this index or the one
    # below it is closest
    if index >= len(dt_list):
        index = len(dt_list) - 1

        # Find the distance between the GITM and observational times
        delta_obj = dt_list[index] - value
        delta = delta_obj.total_seconds()

    elif index >= 0:
        # Ensure that this is the right index and find the distance between
        # the GITM and Observational times
        upper_diff = dt_list[index] - value
        lower_diff = dt_list[index-1] - value

        if upper_diff.total_seconds > lower_diff.total_seconds:
            index -= 1
            delta = lower_diff.total_seconds()
        else:
            delta = upper_diff.total_seconds()

    return delta, index
# End find_nearest_datetime

def match_running_average(ydata_list, xdata, locdata, xmatch, locmatch, xbin,
                          locbin, xwrap=0, locwrap=[0,0,0], *args, **kwargs):
    '''
    Compute a running average given a directed selection criteria.  The x
    axis may be a datetime object, interger, or float.  The y axis data
    must be some sort of number.  The location data must also be numeric.

    ydata_list = list of numpy arrays denoting the dependent variable(s)
    xdata      = numpy array denoting the independent variable
    locdata    = numpy array of coordinates,
                 eg np.array([u0, v0, w0], ..., [un, vn, wn]) 
    xmatch     = numpy array containing the independent variable at the desired
                 output frequency
    locmatch   = numpy array containing the coordinates at the xmatch points
    xbin       = running average box size
    locbin     = bin size of coordinates [ubin, vbin, wbin]
    xwrap      = Wrap offset for x (eg 360 if x = Longitude), default = 0
    locwrap    = Wrap offset for coordinates, default = [0,0,0]

    Output:
    yout_list = list of numpy arrays holding the average
    ystd_list = list of numpy arrays holding the standard deviation
    nout_list = list of numpy arrays holding the number of points per average
    '''
    rname = string.join([module_name, "match_running_average"], sep=" ")

    # Test inputs
    if(len(ydata_list) <= 0 or len(xdata) != len(locdata) or
       len(locdata[0]) != len(locwrap) or len(locmatch) != len(xmatch) or
       len(locdata[0]) != len(locmatch[0]) or xbin == 0.0 or
       len(locdata[0]) != len(locbin)):
        print rname, "WARNING: input arrays don't match or xbin is zero"
        return

    # Initialize the y average, standard deviation, and count
    yout_list = [np.ndarray(shape=xmatch.shape, dtype=float) * 0.0
                 for i in ydata_list]
    ystd_list = [np.ndarray(shape=xmatch.shape, dtype=float) * 0.0
                 for i in ydata_list]
    nout_list = [np.ndarray(shape=xmatch.shape, dtype=int) * 0
                 for i in ydata_list]

    # Initialize the local x and y data inputs
    if xwrap == 0 and sum(locwrap) == 0:
        # If we don't need to buffer the data, copy the x and y data directly
        yraw_list = ydata_list
        xraw = xdata
        locraw = locdata

    else:
        # For x variables such as longitude and local time, buffer data
        # above and below the endpoints to ensure continuity

        yraw_list = list()
        for i,ydat in enumerate(ydata_list):
            yraw_list.append(np.append(ydat, ydat))
            yraw_list[i] = np.append(yraw_list[i], ydat)

        if type(xdata[0]) is dt.datetime:
            xraw = np.append(xdata, xdata)
            xraw = np.append(xraw, xdata)
        else:
            xraw = np.append(xdata - xwrap, xdata)
            xraw = np.append(xraw, xdata + xwrap)

        locraw = np.append(locdata - locwrap, locdata)
        locraw = np.append(locraw, locdata + locwrap)
        locraw = locraw.reshape((-1,len(locwrap)))

    # Cycle through the match array.  If the data's x value falls within the
    # allowable range, increase the data sums.  For this purpose, yout_list
    # will hold the Sum(y), ystd_list will hold the Sum(y**2), and nout_list
    # will hold the number of data points for each x value and data type.
    # If the y data is np.nan, it will not be included and the number of data
    # points for that type will not be incrimented.
    for i,x in enumerate(xraw):
        for m, xm in enumerate(xmatch):
            # Test to see that this data point is within range of this output
            # x location.  Also test that the observation location and output
            # position location are sufficiently similar
            good = True

            for j,lm in enumerate(locmatch[m]):
                if(locraw[i][j] >= lm + locbin[j] * 0.5 or
                   locraw[i][j] < lm - locbin[j] * 0.5):
                    good = False
                    break

            # Datetime differences must be treated differently from numbers
            if type(x) is dt.datetime:
                dt_del = x - xm
                xdel = dt_del.total_seconds()
            else:
                xdel = x - xm

            if(good and xdel < xbin * 0.5 and xdel >= -xbin * 0.5):
                for l, ydat in enumerate(yraw_list):
                    if not np.isnan(ydat[i]):
                        yout_list[l][m] += ydat[i]
                        ystd_list[l][m] += (ydat[i] * ydat[i])
                        nout_list[l][m] += 1
    del xraw, yraw_list, locraw

    # Compute the average and standard deviation
    for l,nout in enumerate(nout_list):
        for i,n in enumerate(nout):
            if(n > 1):
                ystd_list[l][i] = math.sqrt((ystd_list[l][i]-(yout_list[l][i]
                                                              * yout_list[l][i]
                                                              / n)) / (n-1))
                yout_list[l][i] /= n
            else:
                ystd_list[l][i] = 0.0

                if n == 0:
                    yout_list[l][i] = np.nan

    return yout_list, ystd_list, nout_list

# End match_running_average

def match_running_median(ydata_list, xdata, locdata, xmatch, locmatch, xbin,
                         locbin, xwrap=0, locwrap=[0,0,0], *args, **kwargs):
    '''
    Compute a running median given a directed selection criteria.  The x
    axis may be a datetime object, interger, or float.  The y axis data
    must be some sort of number.  The location data must also be numeric.

    ydata_list = list of numpy arrays denoting the dependent variable(s)
    xdata      = numpy array denoting the independent variable
    locdata    = numpy array of coordinates,
                 eg np.array([u0, v0, w0], ..., [un, vn, wn]) 
    xmatch     = numpy array containing the independent variable at the desired
                 output frequency
    locmatch   = numpy array containing the coordinates at the xmatch points
    xbin       = running average box size
    locbin     = bin size of coordinates [ubin, vbin, wbin]
    xwrap      = Wrap offset for x (eg 360 if x = Longitude), default = 0
    locwrap    = Wrap offset for coordinates, default = [0,0,0]

    Output:
    yout_list  = list of numpy arrays holding the median
    yup_list   = list of numpy arrays holding the difference between the upper
                 quartile and the median
    ydown_list = list of numpy arrays holding the difference between the median
                 and the lower quartile
    nout_list  = list of numpy arrays holding the number of points per average
    '''
    rname = string.join([module_name, "match_running_median"], sep=" ")

    # Test inputs
    if(len(ydata_list) <= 0 or len(xdata) != len(locdata) or
       len(locdata[0]) != len(locwrap) or len(locmatch) != len(xmatch) or
       len(locdata[0]) != len(locmatch[0]) or xbin == 0.0 or
       len(locdata[0]) != len(locbin)):
        print rname, "WARNING: input arrays don't match or xbin is zero"
        return

    # Initialize the data arrays
    yhold_list = [[list() for j in xmatch] for i in ydata_list]
    yout_list = [np.ndarray(shape=xmatch.shape, dtype=float) * np.nan
                 for i in ydata_list]
    yup_list = [np.ndarray(shape=xmatch.shape, dtype=float) * 0.0
                for i in ydata_list]
    ydown_list = [np.ndarray(shape=xmatch.shape, dtype=float) * 0.0
                  for i in ydata_list]
    nout_list = [np.ndarray(shape=xmatch.shape, dtype=int) * 0
                 for i in ydata_list]

    # Initialize the local x and y data inputs
    if xwrap == 0 and sum(locwrap) == 0:
        # If we don't need to buffer the data, copy the x and y data directly
        yraw_list = ydata_list
        xraw = xdata
        locraw = locdata

    else:
        # For x variables such as longitude and local time, buffer data
        # above and below the endpoints to ensure continuity
        yraw_list = list()
        for i,ydat in enumerate(ydata_list):
            yraw_list.append(np.append(ydat, ydat))
            yraw_list[i] = np.append(yraw_list[i], ydat)

        if type(xdata[0]) is dt.datetime:
            xraw = np.append(xdata, xdata)
            xraw = np.append(xraw, xdata)
        else:
            xraw = np.append(xdata - xwrap, xdata)
            xraw = np.append(xraw, xdata + xwrap)

        locraw = np.append(locdata - locwrap, locdata)
        locraw = np.append(locraw, locdata + locwrap)
        locraw = locraw.reshape((-1,len(locwrap)))

    # Cycle through the match array.  If the data's x value falls within the
    # allowable range, append the data to a temporary array. If the y data is
    # np.nan, it will not be included.
    for i,x in enumerate(xraw):
        for m,xm in enumerate(xmatch):
            good = True
            for j,lm in enumerate(locmatch[m]):
                if(locraw[i][j] >= lm + locbin[j] * 0.5 or
                   locraw[i][j] < lm - locbin[j] * 0.5):
                    good = False
                    break

            # Datetime differences must be treated differently from numbers
            if type(x) is dt.datetime:
                dt_del = x - xm
                xdel = dt_del.total_seconds()
            else:
                xdel = x - xm

            if(good and xdel < xbin * 0.5 and xdel >= -xbin * 0.5):
                for l, ydat in enumerate(yraw_list):
                    if not np.isnan(ydat[i]):
                        yhold_list[l][m].append(ydat[i])
                        nout_list[l][m] += 1
    del xraw, yraw_list, locraw

    # Compute the first, second,and third quartiles, taking the difference for
    # the upper and lower boundaries
    for l,nout in enumerate(nout_list):
        for i,n in enumerate(nout):
            if(n > 1):
                # Find the median for each x output bin
                yout_list[l][i] = np.median(yhold_list[l][i])

                # Find the upper and lower quartiles (medians of the data above
                # and below the median)

                upper = [y for y in yhold_list[l][i] if y > yout_list[l][i]]
                lower = [y for y in yhold_list[l][i] if y < yout_list[l][i]]

                yup_list[l][i] = np.median(upper) - yout_list[l][i]
                ydown_list[l][i] = yout_list[l][i] - np.median(lower)

                del upper
                del lower
            elif n > 0:
                yout_list[l][i] = yhold_list[l][i][0]

    return yout_list, yup_list, ydown_list, nout_list

# End match_running_average

def gitm_inst_loc(obs_date, obs_lat, obs_lon, obs_alt, dat_keys, obs_type,
                  gitmname_list, gitm_type="1D", mag_file=None,
                  lat_unit="degrees", lon_unit="degrees", alt_unit="km",
                  *args, **kwargs):
    '''
    Extract GITM data for (an) instrument(s) at specified location and times.
    The desired location(s) must be specified in time, latitude, longitude, and
    altitude.  If these are points along a satellite track, the GITM data will
    be output at the GITM file time and the interpolated location for this time
    along the satellite track.  If these are fixed points (eg for a ground 
    receiver), interpolation in time is not necessary.  A GitmTime data
    structure will be returned.  Since this structure requires the same number
    of latitudes, longitudes, and altitudes be available at each time, if the
    number of locations vary, times at a location with no observations will
    be filled with np.nan.

    Input:
    obs_date     = ordered numpy array of datetime objects
    obs_lat      = corresponding numpy array of latitudes 
    obs_lon      = corresponding numpy array of longitudes
    obs_alt      = corresponding numpy array of altitudes
    dat_keys     = list of keys to extract from GITM, empty list returns all
    obs_type     = observation type (satellite or ground)
    gitmbin_list = list of Gitm binary files, ordered by time
    gitm_type    = GITM binary type (1D, 2D, or 3D)
    mag_file     = 3DMAG or 3DION file (default is None)
    lat_unit     = Units of latitudes in track_lat (default degrees)
    lon_unit     = Units of longitude in track_lon (default degrees)
    alt_unit     = Units of altitude in track_alt (default km)

    Output:
    A GitmTime object, uses the PbData class
    '''
    # Local Imports
    from copy import deepcopy as dc
    from scipy import interpolate
    from spacepy.datamodel import dmarray
    import gitm
    import gitm_time
    import gitm_plot_rout as gpr

    rout_name = "gitm_inst_loc"

    # Initialize the list that will hold the desired data in GitmBin structures
    gitmbin_list = list()

    # Define the scale values for each observation location type
    alt_scale = 1.0
    if alt_unit.find("km") >= 0:
        alt_scale = 1000.0

    dlon_scale = 1.0
    rlon_scale = 180.0 / np.pi
    if lon_unit.find("rad") >= 0:
        dlon_scale = np.pi / 180.0
        rlon_scale = 1.0

    dlat_scale = 1.0
    rlat_scale = 180.0 / np.pi
    if lat_unit.find("rad") >= 0:
        dlat_scale = np.pi / 180.0
        rlat_scale = 1.0       

    # Prepare the observational data
    if obs_type.find("satellite") >= 0:
        # AGB Note: It's better to compute an average or median using a
        #           window of data at the satellite's position than it is
        #           to interpolate the satellite data so that a measurement
        #           will be available at the model time if the model and
        #           measurement are to be compared directly.  However, this
        #           will give us model output along the satellite orbit, so
        #           any windowing of the data can be done after the model
        #           estimate along the satellite track is obtained.
        #           
        # Interpolate the satellite's position in time (ms resolution)
        obs_delt = [timedelta.total_seconds(d - obs_date[0]) for d in obs_date]
        itime = list(np.arange(0.0, obs_delt[-1], 0.01))
        nlocs = [1 for d in itime] # Since there is one location for each time
        # Interpolate Altitude
        tck = interpolate.splrep(obs_delt, obs_alt * alt_scale, s=0)
        ialt = interpolate.splev(itime, tck, der=0)
        # Interpolate Latitude
        tck = interpolate.splrep(obs_delt, obs_lat * dlat_scale, s=0)
        ilat = interpolate.splev(itime, tck, der=0)
        # Interpolate Longitude
        tck = interpolate.splrep(obs_delt, obs_lon * dlon_scale, s=0)
        ilon = interpolate.splev(itime, tck, der=0)
    elif obs_type.find("ground") >= 0:
        # Group the observation locations by time.  Since we will seldom
        # require a truely vertical altitude profile (radars typically scan
        # at an angle), altitude profiles will repeat the lat/lon
        itime = list()
        nlocs = list()
        ialt = list()
        ilat = list()
        ilon = list()
        last_time = None
        j = -1

        for i,ot in enumerate(obs_date):
            if not last_time or last_time != ot:
                # This is a new time so initialize the lists
                j += 1
                itime.append(timedelta.total_seconds(ot - obs_date[0]))
                nlocs.append(0)
                ialt.append(list())
                ilat.append(list())
                ilon.append(list())

            # For new times and old times, save the locations
            nlocs[j] += 1
            ilat[j].append(obs_lat[i] * dlat_scale)
            ilon[j].append(obs_lon[i] * dlon_scale)

            if type(obs_alt) is np.ndarray:
                ialt[j].append(obs_alt[i] * alt_scale)
            
            # Update the test condition
            last_time = ot
    else:
        print rout_name, "ERROR: unknown observation data type", obs_type
        return

    # Test to ensure the input data was parsed correctly
    max_locs = max(nlocs)
    if max_locs < 1:
        print rout_name, "ERROR: there are no observation locations"
        return

    # Read the list of Gitm Binary files.  These may or may not be
    # output at the instrument location(s).
    for i, gitm_file in enumerate(gitmname_list):
        split_file = string.split(gitm_file)
        vals = dict()
        nvals = 0

        # Read in the data for this GITM binary file
        if mag_file:
            gdata = gitm.GitmBin(split_file[0], magfile=mag_file,
                                 varlist=dat_keys)
        else:
            gdata = gitm.GitmBin(split_file[0], varlist=dat_keys)

        if gdata.has_key("e-"):
            gdata.calc_2dion()

        # Find the desired location index for this time
        gdelt = timedelta.total_seconds(gdata['time'] - obs_date[0])
        tmax = timedelta.total_seconds(gdata['time'] - obs_date[-1])
        malt = gdata.attrs['nAlt']-1

        # If the GitmBin files extend beyond the observation time, stop
        # cycling through them
        if tmax > 0.0:
            print rout_name, "ADVISEMENT: GitmBin files after", split_file[0], "extend beyond the observation timeframe [", obs_date[0], " to ", obs_date[-1], "]"
            break

        if obs_type.find("satellite") >= 0:
            iobs = int(gdelt * 100)

            if iobs >= len(itime):
                iobs = -1
        else:
            # Find the index of the nearest observation time and difference
            # in time between them
            obs_dist, iobs = find_nearest_value(itime, gdelt)

            # Only use matches with less than 5 min seperation in time
            if obs_dist >= 300.0 or iobs >= len(itime):
                iobs = -1

        # Cycle through the all the locations for this time to interpolate
        # values
        if iobs >= 0 and nlocs[iobs] > 0:
            gkeys = gdata.keys()
            vkeys = list()

            # Extract and assign the values that don't require interpolation
            gkeys.pop(gkeys.index('Latitude'))
            gkeys.pop(gkeys.index('Longitude'))
            gkeys.pop(gkeys.index('dLat'))
            gkeys.pop(gkeys.index('dLon'))
            gkeys.pop(gkeys.index('Altitude'))
            gkeys.pop(gkeys.index('time'))
            gkeys.pop(gkeys.index('LT'))

            # Assign the 2D variables in the file to a seperate processing list
            if gitm_type.find("2") >= 0:
                vkeys = list(gkeys)
                gkeys = list()
            else:
                if gdata.has_key('VTEC'):
                    vkeys.append(gkeys.pop(gkeys.index('VTEC')))
                if gdata.has_key('NmF2'):
                    vkeys.append(gkeys.pop(gkeys.index('NmF2')))
                if gdata.has_key('hmF2'):
                    vkeys.append(gkeys.pop(gkeys.index('hmF2')))

            # Cycle through the observation locations at this time
            for iloc in range(nlocs[iobs]):
                # Get the location information for this point
                this_date = obs_date[0] + timedelta(0, itime[iobs])

                if obs_type.find("satellite") >= 0:
                    this_lat = ilat[iobs]
                    this_lon = ilon[iobs]
                    this_alt = ialt[iobs]
                else:
                    this_lat = ilat[iobs][iloc]
                    this_lon = ilon[iobs][iloc]
                    this_alt = 0.0

                    if type(obs_alt) is np.ndarray:
                        this_alt = ialt[iobs][iloc]

                good_interp = False
                # Perform the appropriate interpolation based on file type
                if gitm_type.find("3") >= 0 or gitm_type.find("2") >= 0:
                    # For a particular time, we need to interpolate the
                    # location in latitude, longitude, and altitude.  First,
                    # find the closest indices
                    (clon, clat) = gpr.find_lon_lat_index(gdata, this_lon,
                                                          this_lat, lat_unit)

                    if gdata["dLon"][clon,clat,0] > this_lon:
                        blon = clon - 1
                        if blon < 0:
                            blon += gdata.attrs['nLon']
                    else:
                        blon = clon
                        clon += 1
                        if clon >= gdata.attrs['nLon']:
                            clon -= gdata.attrs['nLon']

                    if gdata["dLat"][clon,clat,0] > this_lat:
                        blat = clat - 1
                        if blat < 0:
                            blat += gdata.attrs['nLat']
                    else:
                        blat = clat
                        clat += 1
                        if clat >= gdata.attrs['nLat']:
                            clat -= gdata.attrs['nLat']
                    
                    # Only set the altitude for 3D GITM files and instruments
                    if(gitm_type.find("3") >= 0
                       and type(obs_alt) is np.ndarray):
                        calt = gpr.find_alt_index(gdata,clon,clat,this_alt,"m")
                        
                        if gdata['Altitude'][clon,clat,calt] > this_alt:
                            balt = calt - 1
                            if balt < 0:
                                balt = 0
                                calt = 1
                        else:
                            balt = calt
                            calt += 1

                        if calt < malt:
                            good_interp = True
                            # Set up the grid for 3D interpolation
                            alon = gdata["dLon"][blon:clon+1,blat:clat+1,
                                                  balt:calt+1].flatten()
                            alat = gdata["dLat"][blon:clon+1,blat:clat+1,
                                                  balt:calt+1].flatten()
                            aalt = gdata['Altitude'][blon:clon+1, blat:clat+1,
                                                     balt:calt+1].flatten()

                            points = np.ndarray(shape=(len(aalt), 3), dtype=float, buffer=np.array([[l, alat[j], aalt[j]] for j,l in enumerate(alon)]))
                            xi = np.ndarray(shape=(1, 3), dtype=float,
                                            buffer=np.array([this_lon, this_lat,
                                                             this_alt]))

                            # Cycle through the keys to perform the 3D
                            # interpolation
                            for g in gkeys:
                                values = np.ndarray(shape=len(aalt), dtype=float, buffer=np.array(gdata[g][blon:clon+1,blat:clat+1, balt:calt+1].flatten()))
                                v = interpolate.griddata(points, values, xi,
                                                         method="linear")
                                if vals.has_key(g):
                                    vals[g].append(v)
                                else:
                                    vals[g] = list([v])

                    # Initialize the 2D interpolation for 2D and 3D files
                    alon = gdata["dLon"][blon:clon+1,blat:clat+1,0].flatten()
                    alat = gdata["dLat"][blon:clon+1,blat:clat+1,0].flatten()
                    points = np.ndarray(shape=(len(alat), 2), dtype=float, buffer=np.array([[l, alat[j]] for j,l in enumerate(alon)]))
                    xi = np.ndarray(shape=(1, 2), dtype=float,
                                    buffer=np.array([this_lon, this_lat]))

                    # Cycle through the keys to perform 2D interpolation
                    if len(vkeys) > 0:
                        good_interp = True

                    for g in vkeys:
                        values = np.ndarray(shape=len(alat), dtype=float, buffer=np.array(gdata[g][blon:clon+1,blat:clat+1,0].flatten()))
                        v = interpolate.griddata(points, values, xi,
                                                 method="linear")
                        
                        if vals.has_key(g):
                            vals[g].append(v)
                        else:
                            vals[g] = list([v])
                    
                elif gitm_type.find("1") >= 0:
                    # Interpolate the GITM output to the desired altitude
                    galts=np.arange(math.floor(gdata['Altitude'][0,0,0]),
                                    math.ceil(gdata['Altitude'][0,0,malt]), 1.0)
                    gindex = int(round(this_alt) - galts[0])

                    if gindex > 0 and gindex < len(galts):
                        # Assign the values that don't need interpolation
                        for g in vkeys:
                            if vals.has_key(g):
                                vals[g].append(gdata[g][0,0,0])
                            else:
                                vals[g] = list([gdata[g][0,0,0]])

                        # Assign the values that do require interpolation
                        if len(gkeys) > 0:
                            good_interp = True

                        for g in gkeys:
                            gck = interpolate.splrep(gdata['Altitude'][0,0,:],
                                                     gdata[g][0,0,:], s=0)
                            gy = interpolate.splev(galts, gck, der=0)

                            if not g in vals.keys():
                                vals[g] = list()
                            vals[g].append(gy[gindex])
                else:
                    print rout_name, "ERROR: unknown file type [",gitm_type,"]"

                if good_interp:
                    # Save the location data if interpolation was successful
                    nvals += 1

                    if not vals.has_key('Longitude'):
                        vals['Latitude'] = list([this_lat / rlat_scale])
                        vals['Longitude'] = list([this_lon / rlon_scale])
                        vals['dLon'] = list([this_lon / dlon_scale])
                        vals['dLat'] = list([this_lat / dlat_scale])
                        vals['Altitude'] = list([this_alt])
                        vals['LT'] = list([gpr.glon_to_localtime(this_date,
                                                                 this_lon,
                                                                 lon_unit)])
                    else:
                        vals['Latitude'].append(this_lat / rlat_scale)
                        vals['Longitude'].append(this_lon / rlon_scale)
                        vals['dLat'].append(this_lat / dlat_scale)
                        vals['dLon'].append(this_lon / dlon_scale)
                        vals['Altitude'].append(this_alt)
                        vals['LT'].append(gpr.glon_to_localtime(this_date,
                                                                this_lon,
                                                                lon_unit))
            
            # After interpolating at all locations for this file, save the
            # interpolated data to the list of GitmBin structures
            if nvals > 0:
                if nvals == 1 and max_locs == 1:
                    # Save this as a 1D GITM file.  Start by deleting the
                    # unneeded dimensions from the current GitmBin object and
                    # assigning the desired values
                    for k in gdata.keys():
                        if vals.has_key(k):
                            if gdata.attrs['nAlt'] > 1:
                                gdata[k] = np.delete(gdata[k], np.arange(gdata.attrs['nAlt']-1), 2)
                            if gdata.attrs['nLat'] > 1:
                                gdata[k] = np.delete(gdata[k], np.arange(gdata.attrs['nLat']-1), 1)
                            if gdata.attrs['nLon'] > 1:
                                gdata[k] = np.delete(gdata[k], np.arange(gdata.attrs['nLon']-1), 0)

                            gdata[k][0,0,0] = vals[k][0]
                        elif k.find('time') < 0:
                            del gdata[k]
                            gdata.attrs['nVars'] -= 1;
                elif max_locs > 1:
                    # Save this as a 2D GITM file.  Start by reshaping the
                    # dimensions from the current GitmBin object and 
                    # assigning the desired values
                    for k in gdata.keys():
                        if vals.has_key(k):
                            gdata[k] = dmarray(np.empty((max_locs,max_locs,1,)))
                            gdata[k][:] = np.nan

                            # Test to see if this key is lat, lon, or alt/data
                            for iloc in range(nlocs[iobs]):
                                if k.find('Lat') >= 0:
                                    gdata[k][iloc,:,0] = vals[k][iloc]
                                elif k.find('Lon') >= 0:
                                    gdata[k][:,iloc,0] = vals[k][iloc]
                                else:
                                    gdata[k][iloc,iloc,0] = vals[k][iloc]
                        elif k.find('time') < 0:
                            del gdata[k]
                            gdata.attrs['nVars'] -= 1;
                else:
                    print rout_name, "ADVISEMENT: no values interpolated at", this_lat, this_lon, this_alt, gdata['time']
                
                # Finish assigning the metadata
                gdata.attrs['nAlt'] = 1
                gdata.attrs['nLat'] = max_locs
                gdata.attrs['nLon'] = max_locs

                # Save the output
                gitmbin_list.append(dc(gdata))
        else:
            print rout_name, "ADVISEMENT: file [%s] outside of observation time range" % split_file[0]

        del gdata
                
    # Save the interpolated data in a GitmTime object
    if len(gitmbin_list) > 0:
        return(gitm_time.GitmTime(gitmbin_list))

    print rout_name, "ERROR: no data to output"
    return
#END gitm_inst_loc


def gitm_net_loc(obs_date, obs_lat, obs_lon, obs_alt, obs_dat_list,
                 obs_dat_keys, obs_dat_name, obs_dat_scale, obs_dat_units,
                 dat_keys, gitmname_list, gitm_type="3D",
                 mag_file=None, lat_unit="degrees", lon_unit="degrees",
                 alt_unit="km", *args, **kwargs):
    '''
    Extract GITM data for a network of instruments at specified locations and
    times. The desired locations must be specified in time, latitude, and
    longitude (altitude is also acceptable).  A GitmTime data structure will be
    returned.  Since this structure requires the same number of latitudes,
    longitudes, and altitudes be available at each time, if the number of
    locations vary, times at a location with no observations will be filled
    with np.nan.

    Input:
    obs_date      = ordered numpy array of datetime objects
    obs_lat       = corresponding numpy array of latitudes 
    obs_lon       = corresponding numpy array of longitudes
    obs_alt       = corresponding numpy array of altitudes (or an empty list)
    obs_dat_list  = list of numpy arrays containing observational data to append
                    to the output data structure
    obs_dat_keys  = list of key roots corresponding to the observational data
    obs_dat_name  = list of descriptive names corresponding to the
                    observational data
    obs_dat_scale = list of plotting scales corresponding to the observational
                     data
    obs_dat_units = list of units corresponding to the observational data
    dat_keys      = list of keys to extract from GITM, empty list returns all
    gitmbin_list  = list of Gitm binary files, ordered by time
    gitm_type     = GITM binary type (2D or 3D)
    mag_file      = 3DMAG or 3DION file (default is None)
    lat_unit      = Units of latitudes in track_lat (default degrees)
    lon_unit      = Units of longitude in track_lon (default degrees)
    alt_unit      = Units of altitude in track_alt (default km)

    Output:
    A GitmTime object, uses the PbData class
    '''
    # Local Imports
    from copy import deepcopy as dc
    from scipy import interpolate
    from spacepy.datamodel import dmarray
    import gitm
    import gitm_time
    import gitm_plot_rout as gpr

    rout_name = "gitm_net_loc"

    # Initialize the list that will hold the desired data in GitmBin structures
    gitmbin_list = list()

    # Define the scale values for each observation location type
    alt_scale = 1.0
    if alt_unit.find("km") >= 0:
        alt_scale = 1000.0

    dlon_scale = 1.0
    rlon_scale = 180.0 / np.pi
    if lon_unit.find("rad") >= 0:
        dlon_scale = np.pi / 180.0
        rlon_scale = 1.0

    dlat_scale = 1.0
    rlat_scale = 180.0 / np.pi
    if lat_unit.find("rad") >= 0:
        dlat_scale = np.pi / 180.0
        rlat_scale = 1.0       

    # Group the observation locations by time.  Since we will seldom
    # require a truely vertical altitude profile (radars typically scan
    # at an angle), altitude profiles will repeat the lat/lon
    itime = list()
    nlocs = list()
    ilocs = list()
    idata = list()
    last_time = None
    j = -1
    ndim = 2
    ndat = len(obs_dat_list)

    if type(obs_alt) is np.ndarray:
        ndim += 1

    # Test the validity of the data arrays that will be appended
    if ndat > len(obs_dat_keys):
        print rout_name, "WARNING: forgot to include keys for", ndat-len(obs_dat_keys), "observational data types"
    elif ndat < len(obs_dat_keys):
        print rout_name, "WARNING: forgot to include data for", len(obs_dat_keys)-ndat, "observational data keys"
        ndat = len(obs_dat_keys)

    # Sort the observational data by date
    for i,ot in enumerate(obs_date):
        if not last_time or last_time != ot:
            # If this is not the first time, save the old data
            if j >= 0:
                if ndim == 2:
                    xi = np.ndarray(shape=(nlocs[j], ndim), dtype=float, buffer=np.array([[l, ilat[k]] for k,l in enumerate(ilon)]))
                else:
                    xi = np.ndarray(shape=(nlocs[j], ndim), dtype=float, buffer=np.array([[l, ilat[k], ialt[k]] for k,l in enumerate(ilon)]))
                ilocs.append(dc(xi))

                if ndat > 0:
                    idata.append(np.ndarray(shape=(nlocs[j], ndat), dtype=float,
                                            buffer=np.array(idat)))

            # This is a new time so initialize the lists
            j += 1
            itime.append(timedelta.total_seconds(ot - obs_date[0]))
            nlocs.append(0)
            ialt = list()
            ilat = list()
            ilon = list()
            idat = list()

        # For new times and old times, save the locations
        nlocs[j] += 1
        ilat.append(obs_lat[i] * dlat_scale)
        ilon.append(obs_lon[i] * dlon_scale)

        if ndim > 2:
            ialt.append(obs_alt[i] * alt_scale)

        if ndat > 0:
            idat.append([d[i] for d in obs_dat_list])
            
        # Update the test condition
        last_time = ot

    # Save the last set of locations
    if ndim == 2:
        xi = np.ndarray(shape=(nlocs[j], ndim), dtype=float,
                        buffer=np.array([[l, ilat[i]]
                                         for i,l in enumerate(ilon)]))
    else:
        xi = np.ndarray(shape=(nlocs[j], ndim), dtype=float,
                        buffer=np.array([[l, ilat[i], ialt[i]]
                                         for i,l in enumerate(ilon)]))
    ilocs.append(dc(xi))

    if ndat > 0:
        idata.append(np.ndarray(shape=(nlocs[j], ndat), dtype=float,
                                buffer=np.array(idat)))

    # Delete the temporary arrays
    del ilat, ilon, ialt, idat, xi

    # Test to ensure the input data was parsed correctly
    max_locs = max(nlocs)
    max_saved = 0
    if max_locs < 1:
        print rout_name, "ERROR: there are no observation locations"
        return

    # Read the list of Gitm Binary files (default 2D or 3D output types).
    for i, gitm_file in enumerate(gitmname_list):
        split_file = string.split(gitm_file)
        vals = dict()

        # Read in the data for this GITM binary file
        if mag_file:
            gdata = gitm.GitmBin(split_file[0], magfile=mag_file,
                                 varlist=dat_keys)
        else:
            gdata = gitm.GitmBin(split_file[0], varlist=dat_keys)

        if gdata.has_key("e-"):
            gdata.calc_2dion()

        # Find the desired location index for this time
        gdelt = timedelta.total_seconds(gdata['time'] - obs_date[0])
        tmax = timedelta.total_seconds(gdata['time'] - obs_date[-1])
        malt = gdata.attrs['nAlt']-1

        # If the GitmBin files extend beyond the observation time, stop
        # cycling through them
        if tmax > 0.0:
            print rout_name, "ADVISEMENT: GitmBin files after", split_file[0], "extend beyond the observation timeframe [", obs_date[0], " to ", obs_date[-1], "]"
            break

        # Find the index of the nearest observation time and difference
        # in time between them
        obs_dist, iobs = find_nearest_value(itime, gdelt)

        # Only use matches with less than 5 min seperation in time.  Cycle
        # through the all the locations for this time to interpolate values
        if(obs_dist < 300.0 and iobs >= 0 and iobs < len(itime)
           and nlocs[iobs] > 0):
            gkeys = gdata.keys()
            vkeys = list()
            lkeys = dict()

            # Add the observational data to the value dictionary
            if ndat > 0:
                for j,k in enumerate(obs_dat_keys):
                    ok = "obs_{:s}".format(k)
                    vals[ok] = [d[j] for d in idata[iobs]]

            # Extract and assign the values that don't require interpolation
            lkeys[gkeys.pop(gkeys.index('Latitude'))] = 1.0 / rlat_scale
            lkeys[gkeys.pop(gkeys.index('Longitude'))] = 1.0 / rlon_scale
            lkeys[gkeys.pop(gkeys.index('dLat'))] = 1.0 / dlat_scale
            lkeys[gkeys.pop(gkeys.index('dLon'))] = 1.0 / dlon_scale
            lkeys[gkeys.pop(gkeys.index('Altitude'))] = 1.0
            lkeys[gkeys.pop(gkeys.index('LT'))] = lon_unit
            gkeys.pop(gkeys.index('time'))

            # Assign the 2D variables in the file to a seperate processing list
            if gitm_type.find("2") >= 0:
                vkeys = list(gkeys)
                gkeys = list()
            else:
                if gdata.has_key('VTEC'):
                    vkeys.append(gkeys.pop(gkeys.index('VTEC')))
                if gdata.has_key('NmF2'):
                    vkeys.append(gkeys.pop(gkeys.index('NmF2')))
                if gdata.has_key('hmF2'):
                    vkeys.append(gkeys.pop(gkeys.index('hmF2')))

            # Get the time information
            this_date = obs_date[0] + timedelta(0, itime[iobs])
            good_interp = False

            # Set up the grid for 3D interpolation
            if ndim == 3:
                good_interp = True
                alon = gdata["dLon"].flatten()
                alat = gdata["dLat"].flatten()
                aalt = gdata['Altitude'].flatten()
                points = np.ndarray(shape=(len(alat), ndim), dtype=float, buffer=np.array([[l, alat[j], aalt[j]] for j,l in enumerate(alon)]))

                # Cycle through the keys to perform the 3D interpolation
                for g in gkeys:
                    values = np.ndarray(shape=len(alat), dtype=float,
                                        buffer=np.array(gdata[g].flatten()))
                    v = interpolate.griddata(points, values, ilocs[iobs],
                                             method="linear")
                    vals[g] = v

            # Cycle through the keys to perform 2D interpolation
            if len(vkeys) > 0:
                good_interp = True
                alon = gdata["dLon"][:,:,0].flatten()
                alat = gdata["dLat"][:,:,0].flatten()
                points = np.ndarray(shape=(len(alat), ndim), dtype=float, buffer=np.array([[l, alat[j]] for j,l in enumerate(alon)]))

                for g in vkeys:
                    values = np.ndarray(shape=len(alat), dtype=float, buffer=np.array(gdata[g][:,:,0].flatten()))
                    v = interpolate.griddata(points, values, ilocs[iobs],
                                             method="linear")
                        
                    vals[g] = v

        # Save the interpolated data to the list of GitmBin structures
        if good_interp:
            # Initialize the observational data to append
            if ndat > 0:
                for j,k in enumerate(obs_dat_keys):
                    ok = "obs_{:s}".format(k)
                    gdata[ok] = dmarray(np.empty(shape=gdata['dLon'].shape,
                                                 dtype=float) * np.nan,
                                        attrs={"name":obs_dat_name[j],
                                               "scale":obs_dat_scale[j],
                                               "units":obs_dat_units[j]})
                    gdata.attrs['nVars'] += 1

            if max_locs == 1:
                # Save this as a 1D GITM file.  Start by deleting the
                # unneeded dimensions from the current GitmBin object and
                # assigning the desired values
                for k in gdata.keys():
                    if vals.has_key(k) or lkeys.has_key(k):
                        if gdata.attrs['nAlt'] > 1:
                            gdata[k] = np.delete(gdata[k], np.arange(gdata.attrs['nAlt']-1), 2)
                        if gdata.attrs['nLat'] > 1:
                            gdata[k] = np.delete(gdata[k], np.arange(gdata.attrs['nLat']-1), 1)
                        if gdata.attrs['nLon'] > 1:
                            gdata[k] = np.delete(gdata[k], np.arange(gdata.attrs['nLon']-1), 0)
                        if vals.has_key(k):
                            gdata[k][0,0,0] = vals[k][0]
                            if len(vals(k)) > max_saved:
                                max_saved = len(vals(k))
                        else:
                            if k.find("Lat") >= 0:
                                gdata[k][0,0,0] = ilocs[iobs][0][1] * lkeys[k]
                            elif k.find("Lon") >= 0:
                                gdata[k][0,0,0] = ilocs[iobs][0][0] * lkeys[k]
                            elif k.find("Altitude") >= 0:
                                if ndim == 2:
                                    gdata[k][0,0,0] = 0.0
                                else:
                                    gdata[k][0,0,0] = ilocs[iobs][0][2]
                            else:
                                gdata[k][0,0,0] = gpr.glon_to_localtime(this_date, ilocs[iobs][0][0], lon_unit)      
                    elif k.find('time') < 0:
                        # This is not an interpolated value, location, or time
                        del gdata[k]
                        gdata.attrs['nVars'] -= 1
            else:
                # Save this as a 2D or 3D GITM file.  Start by reshaping the
                # dimensions from the current GitmBin object and 
                # assigning the desired values

                for k in gdata.keys():
                    if vals.has_key(k) or lkeys.has_key(k):
                        if ndim == 2:
                            gdata[k] = dmarray(np.empty((max_locs,max_locs,1,)),
                                               attrs=gdata[k].attrs)
                        else:
                            gdata[k] = dmarray(np.empty((max_locs, max_locs,
                                                         max_locs,)),
                                               attrs=gdata[k].attrs)
                        gdata[k][:] = np.nan

                        # Assign the appropriate location or data value
                        for j,loc in enumerate(ilocs[iobs]):
                            jalt = 0
                            if ndim == 3:
                                jalt = 0
                                alt = loc[3]

                            if vals.has_key(k):
                                gdata[k][j,j,jalt] = vals[k][j]
                                if len(vals[k]) > max_saved:
                                    max_saved = len(vals[k])
                            else:
                                if k.find('Lat') >= 0:
                                    gdata[k][j,:,jalt] = loc[1] * lkeys[k]
                                elif k.find('Lon') >= 0:
                                    gdata[k][:,j,jalt] = loc[0] * lkeys[k]
                                elif k.find("Altitude") >= 0:
                                    if ndim == 2:
                                        gdata[k][:,:,:] = 0.0
                                        continue
                                    else:
                                        gdata[k][:,:,jalt] = loc[2]
                                else:
                                    gdata[k][j,j,jalt] = gpr.glon_to_localtime(this_date, loc[0], lon_unit)      
                    elif k.find('time') < 0:
                        del gdata[k]
                        gdata.attrs['nVars'] -= 1

            # Finish assigning the metadata
            gdata.attrs['nLat'] = max_locs
            gdata.attrs['nLon'] = max_locs
            gdata.attrs['nAlt'] = max_locs
            if ndim == 2:
                gdata.attrs['nAlt'] = 1

            # Save the output
            gitmbin_list.append(dc(gdata))
        else:
            print rout_name, "ADVISEMENT: file [%s] outside of observation time range" % split_file[0]

        del gdata
                
    # Save the interpolated data in a GitmTime object
    if len(gitmbin_list) > 0:
        # Remove unnecessary elements if the maximum number of dimensions is
        # larger than the saved number of dimensions
        if max_saved < max_locs:
            for gdata in gitmbin_list:
                # Re-assign the metadata
                gdata.attrs['nLat'] = max_saved
                gdata.attrs['nLon'] = max_saved
                if ndim == 3:
                    gdata.attrs['nAlt'] = max_saved

                for k in gdata.keys():
                    gdata[k] = np.delete(gdata[k], max_locs - 1 -
                                         np.arange(max_locs - max_saved), 0)
                    gdata[k] = np.delete(gdata[k], max_locs - 1 -
                                         np.arange(max_locs - max_saved), 1)
                    if ndim == 3:
                        gdata[k] = np.delete(gdata[k], max_locs - 1 -
                                             np.arange(max_locs - max_saved), 0)

        return(gitm_time.GitmTime(gitmbin_list))

    print rout_name, "ERROR: no data to output"
    return
#END gitm_net_loc
