#!/usr/bin/env python
#-----------------------------------------------------------------------------
# gitm_alt_plots
#
# Author: Angeline G. Burrell, UMichigan, Feb 2013
#
# Comments: Routine to make altitude plots of ionospheric and thermospheric
#           characteristics at a specified location from GITM.
#
# Includes: plot_single_alt_image - plots a single linear or location slice as
#                                   a function of altitude
#           plot_mult_alt_images  - plot multiple locations of linear or 3D
#                                   altitude slices
#           plot_alt_slices       - plot a single 3D altitude contour with
#                                   several linear slices
#           -----------------------------------------------------------------
#           plot_linear_alt       - plot the linear altitude dependence of a
#                                   quantity
#           plot_3D_alt           - plot the altitude dependence of a quantity
#                                   as the function of another spatiotemporal
#                                   coordinate
#----------------------------------------------------------------------------

'''
Plot data from a 3D GITM file for different spatiotemporal coordinates
'''

# Import modules
import sys
import string 
import math
import numpy as np
from spacepy.pybats import gitm
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.cm import get_cmap
from matplotlib.ticker import ScalarFormatter,FormatStrFormatter,MultipleLocator
import gitm_plot_rout as gpr

def plot_single_alt_image(plot_type, zkey, gData, lat_index=-1, lon_index=-1,
                          title=None, figname=None, xkey="dLat", color="b",
                          marker="o", line=":", *args, **kwargs):
    '''
    Creates a rectangular or polar map projection plot for a specified latitude
    range.
    Input: plot_type = key to determine plot type (linear, contour)
           zkey      = key for z variable (ie 'Vertical TEC')
           gData     = gitm bin structure
           lat_index = index of constant latitude (default -1, none)
           lon_index = index of constant longitude (default -1, none)
           title     = plot title
           figname   = file name to save figure as (default is none)
           xkey      = for contour plots specify an x key (default dLat)
                       (options dLat/dLon/Latitude/Longitude)
           color     = linear color (default blue)
           marker    = linear marker type (default circles)
           line      = linear line type (default dotted)
    '''

    # Initialize the variable limits
    zmin, zmax = gpr.find_data_limits([gData], zkey, lat_index, lon_index, -2,6)
    amin, amax = gpr.find_data_limits([gData], "Altitude", lat_index, lon_index,
                                      -2, 30)
    amin = math.ceil(amin / 10000.0) * 10.0
    amax = math.floor(amax / 10000.0) * 10.0

    # Initialize the new figure

    gf = True
    f  = plt.figure()
    ax = f.add_subplot(111)

    if(string.lower(plot_type)=="linear" and lon_index >= 0 and lat_index >= 0):
        con = plot_linear_alt(ax, zkey, gData, zmin, zmax, amin, amax, 6, 6,
                              lon_index, lat_index, title, "t", True, True,
                              color=color, marker=marker, line=line)
    elif(string.lower(plot_type)=="contour"):
        xmin, xmax = gpr.find_data_limits([gData],xkey,lat_index,lon_index,-2,6)

        con = plot_3D_alt(ax, zkey, xkey, gData, zmin, zmax, amin, amax, xmin,
                          xmax, lon_index, lat_index, 6, 10, 6, True, "r",
                          title, "t", True, True)
    else:
        print "ERROR: unknown input type [", plot_type, "]\n"
        gf = False

    if gf:
        # Draw to screen.
        if plt.isinteractive():
            plt.draw() #In interactive mode, you just "draw".
        else:
            # W/o interactive mode, "show" stops the user from typing more 
            # at the terminal until plots are drawn.
            plt.show()

        # Save output file

        if figname is not None:
            plt.savefig(figname)

        return con

# End plot_single_alt_image

def plot_mult_alt_images(plot_type, zkey, gData, lat_index, lon_index,
                         title=None, figname=None, xkey="dLat", color="b",
                         marker="o", line=":", *args, **kwargs):
    '''
    Creates a linear or contour altitude map for a specified altitude range.
    A list of latitude and longitude indexes should be specified.  They may
    be of equal length, or (for a constant value) a length of list one may be
    used.
    Input: plot_type = key to determine plot type (rectangular, polar)
           zkey      = key for z variable (ie 'Vertical TEC')
           gData     = gitm bin structure
           lat_index = list of latitude indices
           lon_index = list of longitude indices
           title     = plot title
           figname   = file name to save figure as (default is none)
           xkey      = x coordinate for contour plots (default dLat)
           color     = line color for linear plots (default blue)
           marker    = marker type for linear plots (default circle)
           line      = line type for linear plots (default dotted line)
    '''

    # Process the index lists
    lat_len = len(lat_index)
    lon_len = len(lon_index)

    # Initialize the x,y,z variable limits
    hold_xmin = []
    hold_xmax = []
    hold_amin = []
    hold_amax = []
    hold_zmin = []
    hold_zmax = []

    for ilat in lat_index:
        for ilon in lon_index:
            tmin, tmax = gpr.find_data_limits([gData], zkey, ilat, ilon, -2, 6)
            hold_zmin.append(tmin)
            hold_zmax.append(tmax)

            tmin, tmax = gpr.find_data_limits([gData], "Altitude", ilat, ilon,
                                              -2, 30)
            hold_amin.append(math.ceil(tmin / 10000.0) * 10.0)
            hold_amax.append(math.floor(tmax / 10000.0) * 10.0)

            if(string.lower(plot_type)=="contour"):
                tmin, tmax = gpr.find_data_limits([gData], xkey, ilat, ilon,
                                                  -2, 6)
                hold_xmin.append(tmin)
                hold_xmax.append(tmax)

    amin = min(hold_amin)
    amax = max(hold_amax)
    zmin = min(hold_zmin)
    zmax = max(hold_zmax)

    if(string.lower(plot_type)=="contour"):
        xmin = min(hold_xmin)
        xmax = max(hold_xmax)

    # Initialize the new figure

    pnum = max([lat_len, lon_len])

    if(pnum < 1):
        print "plot_mult_alt_images ERROR: no altitude regions specified"
        sys.exit(0)

    if(pnum == 1):
        print "plot_mult_alt_images WARNING: only one region, better to use plot_single_alt_image"

    f  = plt.figure()
    tl = " "

    if title:
        f.suptitle(title, size="medium")

    # Adjust the figure height to accomadate the number of subplots

    if(pnum > 2):
        fheight = f.get_figheight()
        f.set_figheight(fheight * 0.5 * pnum)

    for snum in reversed(range(0, pnum)):
        cl   = False
        xl   = False
        yl   = False
        fnum = (pnum * 100) + 11 + snum
        ax   = f.add_subplot(fnum)

        if(pnum == snum + 1):
            xl = True

        if(math.floor(pnum * 0.5) == snum):
            yl = True

        if(snum == 0 and string.lower(plot_type)=="contour"):
            cl = True

        ilon = lon_index[0]
        ilat = lat_index[0]

        if(lon_len > 1):
            ilon = lon_index[snum]

        if(lat_len > 1):
            ilat = lat_index[snum]

        if(string.lower(plot_type)=="linear"):

            con = plot_linear_alt(ax, zkey, gData, zmin, zmax, amin, amax, 6, 5,
                                  ilon, ilat, tl, "r", xl, yl, color, marker,
                                  line)

        elif(string.lower(plot_type)=="contour"):
            con = plot_3D_alt(ax, zkey, xkey, gData, zmin, zmax, amin, amax,
                              xmin, xmax, ilon, ilat, 6, 5, 6, False, "t", tl,
                              "r", xl, yl)

            if(cl == False):
                cpr = list(con.ax.get_position().bounds)

            else:
                # Add and adjust colorbar

                cbar = gpr.add_colorbar(con, zmin, zmax, 6, "horizontal",
                                        gData[zkey].attrs['scale'],
                                        gData[zkey].attrs['name'],
                                        gData[zkey].attrs['units'])

                bp  = list(cbar.ax.get_position().bounds)
                cp  = list(con.ax.get_position().bounds)

                cp[1] = bp[1]
                cp[3] = cpr[3]
                bp[1] = cp[1] + cp[3] + 0.05

                cbar.ax.set_position(bp)
                con.ax.set_position(cp)
        else:
            print "ERROR: unknown input type [", plot_type, "]\n"
            sys.exit(0)

    # Draw to screen.
    if plt.isinteractive():
        plt.draw() #In interactive mode, you just "draw".
    else:
        # W/o interactive mode, "show" stops the user from typing more 
        # at the terminal until plots are drawn.
        plt.show()

    # Save output file

    if figname is not None:
        plt.savefig(figname)

# End plot_mult_alt_images


def plot_alt_slices(zkey, gData, lat_index, lon_index, title=None, figname=None,
                    degrees=True, color="b", marker="o", line=":", *args,
                    **kwargs):
    '''
    Creates a contour altitude map with several linear slices as a function of
    altitude for a specified GITM variable.  A list of latitude and longitude
    indexes should be specified.  One list should consist of a single value,
    the other will be the x variable in the contour plot.  The degree flag
    determines whether x will be ploted in radians or degrees.

    Input: zkey      = key for z variable (ie 'e-')
           gData     = gitm bin structure
           lat_index = list of latitude indices
           lon_index = list of longitude indices
           title     = plot title
           figname   = file name to save figure as (default is none)
           degrees   = plot x label in radians (False) or degrees (default True)
           color     = line color for linear plots (default blue)
           marker    = marker type for linear plots (default circle)
           line      = line type for linear plots (default dotted line)
    '''

    # Process the index lists
    lat_len = len(lat_index)
    lon_len = len(lon_index)
    pnum    = max([lat_len, lon_len])

    if(pnum < 1):
        print "plot_mult_alt_slices ERROR: no altitude slices specified"
        sys.exit(0)

    elif(lat_len > 1 and lon_len > 1):
        print "plot_alt_slices ERROR: one geographic variable must be constant"
        sys.exit(0)

    else:
        ilat = -1
        ilon = -1

        if lat_len == 1:
            ilat = lat_index[0]
            xkey = "Longitude"

            if title:
                title = "%s at %5.2f$^\circ$ N" % (title,
                                                   gData['dLat'][1,ilat,1])
            else:
                title = " %5.2f$^\circ$ N" % (gData['dLat'][1,ilat,1])

            if degrees:
                xkey = "dLon"
        else:
            ilon = lon_index[0]
            xkey = "Latitude"

            if title:
                title = "%s at %5.2f$^\circ$ E" % (title,
                                                   gData['dLon'][ilon,1,1])
            else:
                title = "%5.2f$^\circ$ E" % (gData['dLon'][ilon,1,1])

            if degrees:
                xkey = "dLat"

    # Initialize the x,y,z variable limits
    hold_xmin = []
    hold_xmax = []
    hold_amin = []
    hold_amax = []
    hold_zmin = []
    hold_zmax = []

    for jlat in lat_index:
        if lat_len == 1:
            jlon = -1

        for jlon in lon_index:
            if lon_len == 1:
                jlat = -1

            tmin, tmax = gpr.find_data_limits([gData], zkey, jlat, jlon, -2, 6)
            hold_zmin.append(tmin)
            hold_zmax.append(tmax)

            tmin, tmax = gpr.find_data_limits([gData], "Altitude", jlat, jlon,
                                              -2, 30)
            hold_amin.append(math.ceil(tmin / 10000.0) * 10.0)
            hold_amax.append(math.floor(tmax / 10000.0) * 10.0)

            tmin, tmax = gpr.find_data_limits([gData], xkey, jlat, jlon, -2, 6)
            hold_xmin.append(tmin)
            hold_xmax.append(tmax)

    amin = min(hold_amin)
    amax = max(hold_amax)
    zmin = min(hold_zmin)
    zmax = max(hold_zmax)
    xmin = min(hold_xmin)
    xmax = max(hold_xmax)

    # Initialize the new figure

    f  = plt.figure()
    tl = " "

    if title:
        f.suptitle(title, size="medium")

    # Adjust the figure size to accomadate the number of subplots

    if(pnum > 2):
        fwidth = f.get_figwidth()
        f.set_figwidth(fwidth * 0.5 * pnum)

    # Display the 3D contour plot on top

    gs = gridspec.GridSpec(2,pnum)
    gs.update(hspace=.4)

    ax  = plt.subplot(gs[0,:])
    con = plot_3D_alt(ax, zkey, xkey, gData, zmin, zmax, amin, amax, xmin, xmax,
                      ilon, ilat, 6, 5, 6, True, "t", None, "r", True, True)
  
    # Determine how many x tics to include in each linear slice
    if(pnum < 4):
        zinc = 6
    else:
        zinc = 3

    # Add the linear slices

    for snum in reversed(range(0, pnum)):
        xl   = False
        yl   = False
        ax   = plt.subplot(gs[1,snum])

        if(snum == 0):
            yl = True

        if(math.floor(pnum * 0.5) == snum):
            xl = True

        ilon = lon_index[0]
        ilat = lat_index[0]

        if(lon_len > 1):
            ilon = lon_index[snum]

        if(lat_len > 1):
            ilat = lat_index[snum]

        con = plot_linear_alt(ax, zkey, gData, zmin, zmax, amin, amax, zinc, 5,
                              ilon, ilat, tl, "t", xl, yl, color, marker, line)

    # Draw to screen.
    if plt.isinteractive():
        plt.draw() #In interactive mode, you just "draw".
    else:
        # W/o interactive mode, "show" stops the user from typing more 
        # at the terminal until plots are drawn.
        plt.show()

    # Save output file

    if figname is not None:
        plt.savefig(figname)

# plot_alt_slices

def plot_linear_alt(ax, zkey, gData, zmin, zmax, amin, amax, zinc=6, ainc=6,
                    lon_index=1, lat_index=1, title=None, tloc="t", xl=True,
                    yl=True, color="b", marker="+", line="-", *args, **kwargs):
    '''
    Creates a rectangular map projection plot for a specified latitude range.
    Input: ax        = axis handle
           zkey      = key for z variable (ie 'Vertical TEC')
           gData     = gitm bin structure
           zmin      = minimum value for z variable
           zmax      = maximum value for z variable
           amin      = minimum altitude (in km)
           amax      = maximum altitude (in km)
           zinc      = number of z variable tick incriments (default 6)
           ainc      = number of alt variable tick incriments (default 6)
           lon_index = longitude index (default 1)
           lat_index = latitude index (default 1)
           title     = plot title (default is None)
           tloc      = Specify the title location (t=top, r=right, l=left,
                       b=bottom, default is top)
           xl        = Include x (z variable) label (default is True)
           yl        = Include y (altitude) label (default is True)
           color     = line/marker color (default b [blue])
           marker    = marker type (default +)
           line      = line type (default - [solid line])
    '''

    if(amin == amax or zmin == zmax):
        print "plot_linear_alt ERROR: no x-y range"
        sys.exit(0)

    # Set x and y range values
    arange = amax - amin
    awidth = arange / ainc
    zrange = zmax - zmin
    zwidth = zrange / zinc

    # Plot the values
    con  = ax.plot(gData[zkey][lon_index,lat_index,:],
                   gData['Altitude'][lon_index,lat_index,:] / 1000.0,
                   color=color, marker=marker, linestyle=line)

    # Configure axis
    ytics  = MultipleLocator(awidth)
    ax.yaxis.set_major_locator(ytics)
    if yl:
        ax.set_ylabel('Altitude ($km$)')
    plt.ylim(amin, amax)

    if gData[zkey].attrs.has_key('scale'):
        if gData[zkey].attrs['scale'] == "exponential":
            ax.set_xscale('log')
        else:
            xtics = MultipleLocator(zwidth)
            ax.xaxis.set_major_locator(xtics)
    else:
        xtics = MultipleLocator(zwidth)
        ax.xaxis.set_major_locator(xtics)

    if xl:
        if gData[zkey].attrs.has_key('name') and gData[zkey].attrs.has_key('units'):
            ax.set_xlabel(r'%s ($%s$)' % (gData[zkey].attrs['name'], 
                                        gData[zkey].attrs['units']))
        else:
            ax.set_xlabel(r'%s' % (zkey))
    plt.xlim(zmin, zmax)

    # Set the title
    if title:
        rot  = 'horizontal'
        yloc = 1.05
        xloc = .5

        if(tloc == "l" or tloc == "r"):
            xloc = -.1
            yloc = .5
            rot  = 'vertical'

            if(tloc == "r"):
                xloc = 1.1

        if(tloc == "b"):
            yloc = -.1
            
        if(title == " "):
            ax.set_title(r'%5.2f$^\circ$ E, %4.2f$^\circ$ N'
                         % (gData['dLon'][lon_index,lat_index,1],
                            gData['dLat'][lon_index,lat_index,1]),
                         size='medium', rotation=rot, y=yloc, x=xloc)
        else:
            ax.set_title(r'%s at %5.2f$^\circ$ E, %4.2f$^\circ$ N'
                         % (title, gData['dLon'][lon_index,lat_index,1],
                            gData['dLat'][lon_index,lat_index,1]),
                         size='medium', rotation=rot, y=yloc, x=xloc)

    return con

#End plot_linear_alt

def plot_3D_alt(ax, zkey, xkey, gData, zmin, zmax, amin, amax, xmin, xmax,
                lon_index=1, lat_index=-1, zinc=6, ainc=6, xinc=6, cb=True,
                cloc="r", title=None, tloc="t", xl=True, yl=True, *args,
                **kwargs):
    '''
    Creates a single polar projection, with the latitude center and range
    determined by the input.
    Input: ax         = axis handle
           zkey       = key for z variable (ie 'e-')
           xkey       = key for the x variable (dLat/Lon, Latitude, Longitude)
           gData      = gitm bin structure
           zmin       = minimum value for z variable
           zmax       = maximum value for z variable
           amin       = minimum value for altitude (in km)
           amax       = maximum value for altitude (in km)
           xmin       = minimum value for x variable
           xmax       = maximum value for x variable
           lon_index  = longitude index (default 1, -1 for Longitude type xkey)
           lat_index  = latitude index (default -1 for latitutde type xkey)
           zinc       = number of tick incriments for z variable (default 6)
           ainc       = number of tick incriments for altitude (default 6)
           xinc       = number of tick incriments for x variable (default 6)
           cb         = Add a colorbar (default is True)
           cloc       = Colorbar location (t=top, r=right, l=left, b=bottom, 
                        default is right)
           title      = plot title (default is none)
           tloc       = title location (t=top, r=right, l=left, b=bottom,
                        default is top)
           xl         = Include x (z variable) label (default is True)
           yl         = Include y (altitude) label (default is True)
    '''

    # Test latitude and longitude index input
    lat_flag = False
    lon_flag = False

    if(lat_index == -1):
        lat_flag = True

    if(lon_index == -1):
        lon_flag = True

    if (lon_flag and lat_flag) or (not lon_flag and not lat_flag):
        print "plot_3d_alt ERROR: must hold latitude or longitude constant"
        sys.exit(0)

    # Assign the X, Z, and Altitude data structures
    if lat_flag:
        x = gData[xkey][lon_index,:,:]
        a = gData['Altitude'][lon_index,:,:] / 1000.0
        z = gData[zkey][lon_index,:,:]
    else:
        x = gData[xkey][:,lat_index,:]
        a = gData['Altitude'][:,lat_index,:] / 1000.0
        z = gData[zkey][:,lat_index,:]
        
    # Set range values
    arange = xmax - xmin
    xwidth = arange / xinc
    arange = amax - amin
    awidth = arange / ainc
    arange = zmax - zmin
    zwidth = arange / zinc

    # Set the contour
    v   = np.linspace(zmin, zmax, zinc*10, endpoint=True)
    con = ax.contourf(x, a, z, v, cmap=get_cmap('Spectral_r'), vmin=zmin,
                      vmax=zmax)

    # Configure axis
    ytics  = MultipleLocator(awidth)
    ax.yaxis.set_major_locator(ytics)
    if yl:
        ax.set_ylabel('Altitude ($km$)')
    plt.ylim(amin, amax)

    xtics = MultipleLocator(xwidth)
    ax.xaxis.set_major_locator(xtics)
    if xl:
        if gData[xkey].attrs.has_key('name') and gData[xkey].attrs.has_key('units'):
            ax.set_xlabel(r'%s ($%s$)' % (gData[xkey].attrs['name'], 
                                          gData[xkey].attrs['units']))
        else:
            ax.set_xlabel(r'%s' % (xkey))
    plt.xlim(xmin, xmax)
           
    # Set the title
    if title:
        rot  = 'horizontal'
        yloc = 1.05
        xloc = 0.5

        if tloc == "b":
            yloc = -.1
        elif tloc != "t":
            rot  = 'vertical'
            yloc = 0.5
            xloc = -.2

            if tloc == "r":
                xloc = 1.1

        if lat_flag:
            if title == " ":
                title = ax.set_title(r'%5.2f$^\circ$ E' %
                                     (gData['dLon'][lon_index,1,1]), y=yloc,
                                     size='medium', x=xloc, rotation=rot)
            else:
                title = ax.set_title(r'%s at %5.2f$^\circ$ E' %
                                     (title, gData['dLon'][lon_index,1,1]),
                                     size='medium', y=yloc, x=xloc,
                                     rotation=rot)
        else:
            if title == " ":
                title = ax.set_title(r'%5.2f$^\circ$ N' %
                                     (gData['dLat'][1,lat_index,1]), y=yloc,
                                     size='medium', x=xloc, rotation=rot)
            else:
                title = ax.set_title(r'%s at %5.2f$^\circ$ N' %
                                     (title, gData['dLat'][1,lat_index,1]),
                                     size='medium', y=yloc, x=xloc,
                                     rotation=rot)
 
    # Add a colorbar
    if cb:
        orient = 'vertical'

        if(cloc == 't' or cloc == 'b'):
            orient = 'horizontal'

        cbar = gpr.add_colorbar(con, zmin, zmax, zinc, orient,
                                gData[zkey].attrs['scale'],
                                gData[zkey].attrs['name'],
                                gData[zkey].attrs['units'])

        if(cloc == 'l' or cloc == 't'):
            bp = list(cbar.ax.get_position().bounds)
            cp = list(con.ax.get_position().bounds)

            if(cloc == 't'):
                cp[1] = bp[1]
                bp[1] = cp[1] + cp[3] + 0.075
            else:
                bp[0] = 0.125
                cp[0] = bp[0] + 0.1 + bp[2]

            con.ax.set_position(cp)
            cbar.ax.set_position(bp)

    return con

#End plot_3D_alt

#End gitm_alt_plots

