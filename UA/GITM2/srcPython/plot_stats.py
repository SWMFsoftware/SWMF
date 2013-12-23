#!/usr/bin/env python
#-----------------------------------------------------------------------------
# plot_stats
#
# Author: Angeline G. Burrell, UMichigan, Dec 2013
#
# Comments: Script routine to create a plot showing different statistical
#           properties
#
# Includes: add_stat_box - adds a box with desired statistics to a plot
#           plot_hist_w_stats - routine to plot a histogram and statistics
#----------------------------------------------------------------------------

import string
import numpy as np
import operator as op
from scipy import stats
from matplotlib.ticker import FormatStrFormatter, MultipleLocator
import matplotlib.pyplot as plt
#import pylab as P <-- uncomment this when distribution lines are added

def add_stat_box(ax, x_data, x_units, moments=4, quartiles=3, modes=True,
                 loc="ir",  *args, **kwargs):
    '''
    Adds the specified statistics in a text box to a plot

    Input: ax        = plot handle
           x_data    = list or numpy array containing x-axis data
           x_units   = Units for x data
           moments   = Include moments (mean, standard deviation, skew,
                       kurtosis) of the distribution.  The number corresponds
                       to the maximum moment that will be included.  3, for
                       example, will include the mean, standard deviation, and
                       skew. (default=4)
           quartiles = Include quartiles (1st, 2nd=median, 3rd). 1=median only,
                       2=1st and 3rd only, 3=all (default=3)
           modes     = Include mode(s) (default=True)
           loc       = Text location requires two specifiers. The first letter
                       may be i=inside plot, o=outside plot.  The second letter
                       may be r=right, l=left, b=bottom, t=top. (default="ir")
    '''
    stat_text = ""

    # Add moments
    if moments >= 1:
        stat_text = "{:s}$\mu$ = {:g} ${:s}$\n".format(stat_text, round(np.mean(x_data), 2), x_units)
    if moments >= 2:
        stat_text = "{:s}$\sigma$ = {:g} ${:s}$\n".format(stat_text, round(np.std(x_data), 2), x_units)
    if moments >= 3:
        stat_text="{:s}$\gamma$ = {:g}\n".format(stat_text,round(stats.skew(x_data), 2))
    if moments >= 4:
        stat_text="{:s}$\kappa$ = {:g}\n".format(stat_text,
                                                 round(stats.kurtosis(x_data), 2))
    # Add quartiles
    if quartiles == 1:
        stat_text="{:s}$q_2$ = {:g} ${:s}$\n".format(stat_text, round(np.median(x_data), 2), x_units)
    elif quartiles == 2:
        stat_text = "{:s}$q_1$ = {:g} ${:s}$\n$q_3$ = {:g} ${:s}$\n".format(stat_text, *reduce(op.add, [[round(val, 2), x_units] for val in stats.mstats.mquantiles(x_data, prob=[.25, .75])]))
    elif quartiles == 3:
        stat_text = "{:s}$q_1$ = {:g} ${:s}$\n$q_2$ = {:g} ${:s}$\n$q_3$ = {:g} ${:s}$\n".format(stat_text, *reduce(op.add, [[round(val, 2), x_units] for val in stats.mstats.mquantiles(x_data, prob=[.25, .5, .75])]))

    # Add mode(s)
    if modes:
        m, n = stats.mode(x_data)
        n = int(n)
        if n > 0:
            i = 1
            stat_text = "{:s}$n_{{modes}}$ = {:d}\nmodes = {:g} ${:s}$".format(stat_text, n, round(m[0], 2), x_units)
            while i < n:
                stat_text = "{:s}, {:g} ${:s}$".format(stat_text,
                                                       round(m[i], 2), x_units)
                i += 1
            stat_text = "{:s}\n".format(stat_text)

    # Output the text in the desired location (initialize with loc="it")
    va = "top"
    ha = "center"
    xloc = 0.5
    yloc = 0.95

    if loc[0] == "i":
        if loc[1] == "b":
            va = "bottom"
            yloc = 0.05
        elif loc[1] == "r":
            ha = "right"
            xloc = 0.95
        elif loc[1] == "l":
            ha = "left"
            xloc = 0.05
    else:
        if loc[1] == "t":
            va = "bottom"
            yloc = 1.05
        elif loc[1] == "b":
            yloc = -0.05
        elif loc[1] == "r":
            va = "center"
            ha = "right"
            xloc = 1.05
            yloc = 0.05
        elif loc[1] == "l":
            va = "center"
            ha = "left"
            xloc = -0.05
            yloc = 0.05

    ax.text(xloc, yloc, stat_text, horizontalalignment=ha, verticalalignment=va,
            transform=ax.transAxes)
    return(stat_text)

def plot_hist_w_stats(ax, x_data, x_name, x_units, res=50, norm=1.0,
                      weights=None, cum=False, htype="bar", align="mid",
                      rwidth=None, n_label="Fractional Counts", xmin=None,
                      xmax=None, nmin=None, nmax=None, xinc=6, ninc=6,
                      color="b", title=None, tloc="t", xl=True, xt=True,
                      yl=True, yt=True, moments=4, quartiles=3, modes=False,
                      sloc="il", *args, **kwargs):
    '''
    Creates a single histogram plot with a text box showing the desired stats.

    Input: ax        = axis handle
           x_data    = list or numpy array containing x-axis data
           x_name    = Label name of x-axis
           x_units   = Units for x-axis
           res       = Histogran bin resolution. A sequence can be used to
                       provide unequally spaced bins. (default=50)
           norm      = Histogram normalization (eg False, 100), (default=1.0)
           weights   = list or numpy array of the same dimention as x_data to
                       weight each x point contribution. (default=None)
           cum       = Cumulative or not? (default=False)
           htype     = Histogram type: bar, barstacked, step, stepfilled
                       (default="bar")
           align     = Bar alignment: left, mid, right (default="mid")
           rwidth    = Specifies the bar width for htype=bar,barstacked
                       (default=None)  
           n_label   = y-axis label
           xmin      = minimum value for x axis (default=None)
           xmax      = maximum value for x axis (default=None)
           nmin      = minimum value for y axis (default=None)
           nmax      = maximum value for y axis (default=None)
           xinc      = number of tick incriments for x variable (default 6)
           ninc      = number of tick incriments for y variable (default 6)
           color     = Color value (default="b", blue)
           title     = plot title (default is none)
           tloc      = title location: t=top, r=right, l=left, b=bottom
                       (default="t")
           xl        = Include x label (default is True)
           xt        = Include x ticks (default is True)
           yl        = Include y label (default is True)
           yt        = Include y ticks (default is True)
           moments   = Include moments (mean, standard deviation, skew,
                       kurtosis) of the distribution.  The number corresponds
                       to the maximum moment that will be included.  3, for
                       example, will include the mean, standard deviation, and
                       skew. (default=4)
           quartiles = Include quartiles (1st, 2nd=median, 3rd). 1=median only,
                       2=1st and 3rd only, 3=all (default=3)
           modes     = Include mode(s) (default=False)
           sloc      = Text location requires two specifiers. The first letter
                       may be i=inside plot, o=outside plot.  The second letter
                       may be r=right, l=left, b=bottom, t=top. (default="il")
    '''
    # Add histogram to specified subplot
    n,bins,patches=ax.hist(x_data, res, normed=norm, weights=weights,
                           cumulative=cum, histtype=htype, align=align,
                           rwidth=rwidth, color=color)

    # Set the x, and y ranges if desired
    if(xmin is None):
        xmin = np.nanmin(x_data)
    if(xmax is None):
        xmax = np.nanmax(x_data)

    if(nmin is None):
        nmin = np.nanmin(n)
    if(nmax is None):
        nmax = np.nanmax(n)

    # Configure axis
    if yt:
        width = (nmax - nmin) / ninc
        ytics = MultipleLocator(width)
        ax.yaxis.set_major_locator(ytics)
        if norm == 1.0:
            ax.yaxis.set_major_formatter(FormatStrFormatter("%1.2f"))
        elif norm == 100.0:
            ax.yaxis.set_major_formatter(FormatStrFormatter("%3.1f"))
    else:
        ax.yaxis.set_major_formatter(FormatStrFormatter(""))

    if yl:
        ax.set_ylabel(n_label)
    plt.ylim(nmin, nmax)

    if xt:
        width = (xmax - xmin) / xinc
        xtics = MultipleLocator(width)
        ax.xaxis.set_major_locator(xtics)
    else:
        ax.xaxis.set_major_formatter(FormatStrFormatter(""))

    if xl:
        ax.set_xlabel(r'%s ($%s$)' % (x_name, x_units))
    plt.xlim(xmin, xmax)
           
    # Set the title
    if title:
        rot  = 'horizontal'
        yloc = 1.02
        xloc = 0.5

        if tloc == "b":
            yloc = -.1
        elif tloc != "t":
            rot  = 'vertical'
            yloc = 0.5
            xloc = -.2

            if tloc == "r":
                xloc = 1.1

        title = ax.set_title(title, x=xloc, y=yloc, size='medium', rotation=rot)
 
    # Add the statistics, if desired
    st = add_stat_box(ax, x_data, x_units, moments=moments, quartiles=quartiles,
                      modes=modes, loc=sloc)
    return(st)


def plot_lat_lt_stats(lat_data, lat_name, lat_units, sza_data, sza_units,
                      sza_east, x_data, x_name, x_units, datmin=50, 
                      latlim=[55.0, 15.0, -15.0, -55.0], szalim=[108.0, 90.0],
                      res=50, norm=False, weights=None, cum=False, htype="bar",
                      align="mid", rwidth=None, n_label="Counts", xmin=None,
                      xmax=None, nmin=None, nmax=None, xinc=6, ninc=6,
                      color="b", title=None, tloc="t", xl=True, xt=True,
                      yl=True, yt=True, moments=4, quartiles=3, modes=False,
                      sloc="il", figroot=None, draw=True, *args, **kwargs):
    '''
    Creates a single histogram plot with a text box showing the desired stats.

    Input: lat_data  = Numpy array containing lat data.  To only plot
                       global or local time regions, make this empty
           lat_name  = name of latitude data (eg Magnetic Latitude)
           lat_units = Latitude units (degrees or radians)
           sza_data  = Numpy array containing solar zenith angle.  To only plot
                       global or latitude regions, make this empty
           sza_units = local time units (radians or degrees)
           sza_east  = Numpy array of Boolian values indicating whether or not
                       the observation is east of the subsolar point (True=East)
           x_data    = Numpy array containing x-axis data
           x_name    = Label name of x-axis
           x_units   = Units for x-axis
           datmin    = Minimum number of observations needed to compute stats
                       (default=50)
           latlim    = Latitude limit in latitude units [north polar southern
                       border, equatorial northern limit, equatorial southern
                       limit, south polar northern border]
                       (default=[55.0,15.0,-15.0,-55.0])
           szalim    = Solar Zenith angle limit in SZA units ]twilight
                       nightside, twilight dayside] (default=[108.0, 90.0])
           res       = Histogran bin resolution. A sequence can be used to
                       provide unequally spaced bins. (default=50)
           norm      = Histogram normalization (eg False, 100), (default=False)
           weights   = list or numpy array of the same dimention as x_data to
                       weight each x point contribution. (default=None)
           cum       = Cumulative or not? (default=False)
           htype     = Histogram type: bar, barstacked, step, stepfilled
                       (default="bar")
           align     = Bar alignment: left, mid, right (default="mid")
           rwidth    = Specifies the bar width for htype=bar,barstacked
                       (default=None)  
           n_label   = y-axis label (default="Counts")
           xmin      = minimum value for x axis (default=None)
           xmax      = maximum value for x axis (default=None)
           nmin      = minimum value for y axis (default=None)
           nmax      = maximum value for y axis (default=None)
           xinc      = number of tick incriments for x variable (default 6)
           ninc      = number of tick incriments for y variable (default 6)
           color     = Color value (default="b", blue)
           title     = plot title (default=None)
           tloc      = title location: t=top, r=right, l=left, b=bottom
                       (default="t")
           xl        = Include x label (default is True)
           xt        = Include x ticks (default is True)
           yl        = Include y label (default is True)
           yt        = Include y ticks (default is True)
           moments   = Include moments (mean, standard deviation, skew,
                       kurtosis) of the distribution.  The number corresponds
                       to the maximum moment that will be included.  3, for
                       example, will include the mean, standard deviation, and
                       skew. (default=4)
           quartiles = Include quartiles (1st, 2nd=median, 3rd). 1=median only,
                       2=1st and 3rd only, 3=all (default=3)
           modes     = Include mode(s) (default=False)
           sloc      = Text location requires two specifiers. The first letter
                       may be i=inside plot, o=outside plot.  The second letter
                       may be r=right, l=left, b=bottom, t=top. (default="il")
           figroot   = File root name for output figures (eg "Hist.01" for
                       Hist.01.all.png, Hist.01.eq.png, etc) (default=None)
           draw      = draw figures to screen (default=True)
    '''
    rout_name = "plot_lat_lt_stats"
    f = list()

    # Start by computing statistics for all data
    if x_data.shape[0] < datmin:
        print rout_name, "WARNING: insufficient data to compute statistics"
        return f

    hist_lab = "{:s} (N={:d})".format(n_label, len(x_data))
    htitle = "All Locations"
    if title is not None:
        htitle = "{:s}\n{:s}".format(htitle, title)

    f.append(plt.figure())
    ax = f[-1].add_subplot(1,1,1)
    plot_hist_w_stats(ax, x_data, x_name, x_units, res=res, norm=norm,
                      weights=weights, cum=cum, htype=htype, align=align,
                      rwidth=rwidth, n_label=hist_lab, xmin=xmin, xmax=xmax,
                      nmin=nmin, nmax=nmax, xinc=xinc, ninc=ninc, color=color,
                      title=htitle, tloc=tloc, xl=xl, xt=xt, yl=yl, yt=yt,
                      moments=moments, quartiles=quartiles, modes=modes,
                      sloc=sloc)
    if draw and plt.isinteractive():
        plt.draw() #In interactive mode, you just "draw".

    if figroot:
        figname = figroot+'.All.png'
        plt.savefig(figname)

    # Compute statistics for all LT by latitude region:
    # NPolar (>55), NMid (55-15), Eq (+/- 15), SMid (-15- -55), SPolar (<-55)
    if lat_data.shape == x_data.shape:
        # Divide data by region and define labels
        lat_dawn_dusk = list()
        lat_index = list()
        lat_label = list()
        lat_names = list()
        lat_root = list()
        # North Polar
        lat_index.append([i for i,l in enumerate(lat_data) if l > latlim[0]])
        if len(lat_index[-1]) >= datmin:
            lat_names.append("North Polar")
            lat_root.append("NPole")
            lat_dawn_dusk.append(False)
            if lat_units.find("rad") >= 0:
                lat_label.append("{:s} (>{:.0f})\n{:s} (N={:d})".format(lat_name, latlim[0], n_label, len(lat_index[-1])))
            else:
                lat_label.append("{:s} (>{:.0f}$^\circ$)\n{:s} (N={:d})".format(lat_name, latlim[0], n_label, len(lat_index[-1])))
        else:
            lat_index.pop()
        # North Mid-Latitude
        lat_index.append([i for i,l in enumerate(lat_data)
                          if l <= latlim[0] and l >= latlim[1]])
        if len(lat_index[-1]) >= datmin:
            lat_names.append("North Mid-Latitudes")
            lat_root.append("NMid")
            lat_dawn_dusk.append(True)
            if lat_units.find("rad") >= 0:
                lat_label.append("{:s} ({:.0f} $-$ {:.0f})\n{:s} (N={:d})".format(lat_name, latlim[0], latlim[1], n_label, len(lat_index[-1])))
            else:
                lat_label.append("{:s} ({:.0f}$^\circ$ $-$ {:.0f}$^\circ$)\n{:s} (N={:d})".format(lat_name, latlim[0], latlim[1], n_label, len(lat_index[-1])))
        else:
            lat_index.pop()
        # Equatorial
        lat_index.append([i for i,l in enumerate(lat_data)
                          if l < latlim[1] and l > latlim[2]])
        if len(lat_index[-1]) >= datmin:
            lat_names.append("Equatorial")
            lat_root.append("Eq")
            lat_dawn_dusk.append(True)
            if lat_units.find("rad") >= 0:
                lat_label.append("{:s} ($\pm${:.0f})\n{:s} (N={:d})".format(lat_name, latlim[1], n_label, len(lat_index[-1])))
            else:
                lat_label.append("{:s} ($\pm${:.0f}$^\circ$)\n{:s} (N={:d})".format(lat_name, latlim[1], n_label, len(lat_index[-1])))
        else:
            lat_index.pop()
        # South Mid-Latitude
        lat_index.append([i for i,l in enumerate(lat_data)
                          if l <= latlim[2] and l >= latlim[3]])
        if len(lat_index[-1]) >= datmin:
            lat_names.append("South Mid-Latitudes")
            lat_root.append("SMid")
            lat_dawn_dusk.append(True)
            if lat_units.find("rad") >= 0:
                lat_label.append("{:s} ({:.0f} $-$ {:.0f})\n{:s} (N={:d})".format(lat_name, latlim[2], latlim[3], n_label, len(lat_index[-1])))
            else:
                lat_label.append("{:s} ({:.0f}$^\circ$ $-$ {:.0f}$^\circ$)\n{:s} (N={:d})".format(lat_name, latlim[2], latlim[3], n_label, len(lat_index[-1])))
        else:
            lat_index.pop()
        # South Polar
        lat_index.append([i for i,l in enumerate(lat_data) if l < latlim[3]])
        if len(lat_index[-1]) >= datmin:
            lat_names.append("South Polar")
            lat_root.append("SPole")
            lat_dawn_dusk.append(False)
            if lat_units.find("rad") >= 0:
                lat_label.append("{:s} (<{:.0f})\n{:s} (N={:d})".format(lat_name, latlim[3], n_label, len(lat_index[-1])))
            else:
                lat_label.append("{:s} (<{:.0f}$^\circ$)\n{:s} (N={:d})".format(lat_name, latlim[3], n_label, len(lat_index[-1])))
        else:
            lat_index.pop()

        nsub = len(lat_index)
        if nsub == 0:
            print rout_name, "ADVISEMENT: not enough data for lat divisions"
        else:
            # Initialize figure
            f.append(plt.figure())
            if nsub > 1:
                fheight = f[-1].get_figheight()
                f[-1].set_figheight(fheight * nsub)

            htitle = "{:s} Regions".format(lat_name)
            if title is not None:
                htitle = "{:s}\n{:s}".format(htitle, title)
            f[-1].suptitle(htitle)

            # Plot each subfigure
            for i, dat_index in enumerate(lat_index):
                ax = f[-1].add_subplot(nsub,1,i+1)

                try: nweights = weights[dat_index]
                except: nweights = None

                plot_hist_w_stats(ax, x_data[dat_index], x_name, x_units,
                                  res=res, norm=norm, weights=nweights, cum=cum,
                                  htype=htype, align=align, rwidth=rwidth,
                                  n_label=lat_label[i], xmin=xmin, xmax=xmax,
                                  nmin=nmin, nmax=nmax, xinc=xinc, ninc=ninc,
                                  color=color, tloc=tloc, xl=xl, xt=xt, yl=yl,
                                  yt=yt, moments=moments, quartiles=quartiles,
                                  modes=modes, sloc=sloc)

            # Adjust subplot locations
            plt.subplots_adjust(left=.15)
            if draw and plt.isinteractive():
                plt.draw() #In interactive mode, you just "draw".

            if figroot:
                figname = figroot+'.Lat.png'
                plt.savefig(figname)
        
            # Compute statistics for LT regions in each latitude region:
            # Dawn/Twilight (108-90), Day (<90), Dusk/Twilight (108-90), and
            # Night (> 108)
            if sza_data.shape==x_data.shape and sza_east.shape==x_data.shape:
                # Cycle through each of the latitude regions with data
                for i, dat_index in enumerate(lat_index):
                    # Divide into the desired LT regions using the solar
                    # zenith angle
                    lt_index = list()
                    lt_label = list()
                    # Day
                    lt_index.append([j for j in dat_index
                                     if sza_data[j] < szalim[1]])
                    if len(lt_index[-1]) >= datmin:
                        if sza_units.find("rad") >= 0:
                            lt_label.append("Day (SZA<{:.0f})\n{:s} (N={:d})".format(szalim[1], n_label, len(lt_index[-1])))
                        else:
                            lt_label.append("Day (SZA<{:.0f}$^\circ$)\n{:s} (N={:d})".format(szalim[1], n_label, len(lt_index[-1])))
                    else:
                        lt_index.pop()
                    # Night
                    lt_index.append([j for j in dat_index
                                     if sza_data[j]>szalim[0]])
                    if len(lt_index[-1]) >= datmin:
                        if sza_units.find("rad") >= 0:
                            lt_label.append("Night (SZA>{:.0f})\n{:s} (N={:d})".format(szalim[0], n_label, len(lt_index[-1])))
                        else:
                            lt_label.append("Night (SZA>{:.0f}$^\circ$)\n{:s} (N={:d})".format(szalim[0], n_label, len(lt_index[-1])))
                    else:
                        lt_index.pop()

                    if lat_dawn_dusk[i]:
                        # Dawn
                        lt_index.append([j for j in dat_index
                                         if(sza_data[j] >= szalim[1]
                                            and sza_data[j] <= szalim[0]
                                            and sza_east[j])])
                        if len(lt_index[-1]) >= datmin:
                            if sza_units.find("rad") >= 0:
                                lt_label.append("Dawn (SZA {:.0f} $-$ {:.0f})\n{:s} (N={:d})".format(szalim[1], szalim[0], n_label, len(lt_index[-1])))
                            else:
                                lt_label.append("Dawn (SZA {:.0f}$^\circ$ $-$ {:.0f}$^\circ$)\n{:s} (N={:d})".format(szalim[1], szalim[0], n_label, len(lt_index[-1])))
                        else:
                            lt_index.pop()
                        # Dusk
                        lt_index.append([j for j in dat_index
                                         if(sza_data[j] >= szalim[1]
                                            and sza_data[j] <= szalim[0]
                                            and not sza_east[j])])
                        if len(lt_index[-1]) >= datmin:
                            if sza_units.find("rad") >= 0:
                                lt_label.append("Dusk (SZA {:.0f} $-$ {:.0f})\n{:s} (N={:d})".format(szalim[1], szalim[0], n_label, len(lt_index[-1])))
                            else:
                                lt_label.append("Dusk (SZA {:.0f}$\circ$ $-$ {:.0f}$^\circ$)\n{:s} (N={:d})".format(szalim[1], szalim[0], n_label, len(lt_index[-1])))
                        else:
                            lt_index.pop()
                    else:
                        # Twilight
                        lt_index.append([j for j in dat_index
                                         if(sza_data[j] >= szalim[1]
                                            and sza_data[j] <= szalim[0])])
                        if len(lt_index[-1]) >= datmin:
                            if sza_units.find("rad") >= 0:
                                lt_label.append("Twilight (SZA {:.0f} $-$ {:.0f})\n{:s} (N={:d})".format(szalim[1], szalim[0], n_label, len(lt_index[-1])))
                            else:
                                lt_label.append("Twilight (SZA {:.0f}$^\circ$ $-$ {:.0f}$^\circ$)\n{:s} (N={:d})".format(szalim[1], szalim[0], n_label, len(lt_index[-1])))
                        else:
                            lt_index.pop()

                    nsub = len(lt_index)

                    if nsub == 0:
                        print(rout_name, "ADVISEMENT: not enough data for",
                              lat_names[i], "LT divisions")
                    else:
                        # Initialize Figure
                        f.append(plt.figure())
                        if nsub > 1:
                            fheight = f[-1].get_figheight()
                            f[-1].set_figheight(fheight * nsub)

                        htitle = "{:s} Local Time Regions".format(lat_names[i])
                        if title is not None:
                            htitle = "{:s}\n{:s}".format(htitle, title)
                        f[-1].suptitle(htitle)

                        for j, dat2_index in enumerate(lt_index):
                            ax = f[-1].add_subplot(nsub,1,j+1)

                            try: nweights = weights[dat2_index]
                            except: nweights = None

                            plot_hist_w_stats(ax, x_data[dat2_index], x_name,
                                              x_units, res=res, norm=norm,
                                              weights=nweights, cum=cum,
                                              htype=htype, align=align,
                                              rwidth=rwidth,
                                              n_label=lt_label[j], xmin=xmin,
                                              xmax=xmax, nmin=nmin, nmax=nmax,
                                              xinc=xinc, ninc=ninc, color=color,
                                              tloc=tloc, xl=xl, xt=xt, yl=yl,
                                              yt=yt, moments=moments,
                                              quartiles=quartiles, modes=modes,
                                              sloc=sloc)

                        # Adjust subplot locations
                        plt.subplots_adjust(left=.15)
                        if draw and plt.isinteractive():
                            plt.draw() #In interactive mode, you just "draw".

                        if figroot:
                            figname = '{:s}.{:s}LT.png'.format(figroot,
                                                               lat_root[i])
                            plt.savefig(figname)
    
    # Compute statistics for all latitudes by LT region:
    # Dawn/Twilight (108-90), Day (<90), Dusk/Twilight (108-90), Night (> 108)
    if sza_data.shape == x_data.shape and sza_east.shape == sza_data.shape:
        # Divide into the desired LT regions using the solar zenith angle
        lt_index = list()
        lt_label = list()
        # Day
        lt_index.append([i for i,l in enumerate(sza_data) if l < szalim[1]])
        if len(lt_index[-1]) >= datmin:
            if sza_units.find("rad") >= 0:
                lt_label.append("Day (SZA <{:.0f})\n{:s} (N={:d})".format(szalim[1], n_label, len(lt_index[-1])))
            else:
                lt_label.append("Day (SZA <{:.0f}$^\circ$)\n{:s} (N={:d})".format(szalim[1], n_label, len(lt_index[-1])))
        else:
            lt_index.pop()
        # Night
        lt_index.append([i for i,l in enumerate(sza_data) if l > szalim[0]])
        if len(lt_index[-1]) >= datmin:
            if sza_units.find("rad") >= 0:
                lt_label.append("Night (SZA >{:.0f})\n{:s} (N={:d})".format(szalim[0], n_label, len(lt_index[-1])))
            else:
                lt_label.append("Night (SZA >{:.0f}$^\circ$)\n{:s} (N={:d})".format(szalim[0], n_label, len(lt_index[-1])))
        else:
            lt_index.pop()
        # Dawn
        lt_index.append([i for i,l in enumerate(sza_data)
                         if(l >= szalim[1] and l <= szalim[0] and sza_east[i])])
        if len(lt_index[-1]) >= datmin:
            if sza_units.find("rad") >= 0:
                lt_label.append("Dawn (SZA {:.0f} $-$ {:.0f})\n{:s} (N={:d})".format(szalim[1], szalim[0], n_label, len(lt_index[-1])))
            else:
                lt_label.append("Dawn (SZA {:.0f}$^\circ$ $-$ {:.0f}$^\circ$)\n{:s} (N={:d})".format(szalim[1], szalim[0], n_label, len(lt_index[-1])))
        else:
            lt_index.pop()
        # Dusk
        lt_index.append([i for i,l in enumerate(sza_data)
                         if(l >= szalim[1] and l <= szalim[0]
                            and not sza_east[i])])
        if len(lt_index[-1]) >= datmin:
            if sza_units.find("rad") >= 0:
                lt_label.append("Dusk (SZA {:.0f} $-$ {:.0f})\n{:s} (N={:d})".format(szalim[1], szalim[0], n_label, len(lt_index[-1])))
            else:
                lt_label.append("Dusk (SZA {:.0f}$^\circ$ $-$ {:.0f}$^\circ$)\n{:s} (N={:d})".format(szalim[1], szalim[0], n_label, len(lt_index[-1])))
        else:
            lt_index.pop()

        # Divide data by region
        nsub = len(lt_index)
        if nsub == 0:
            print rout_name, "ADVISEMENT: not enough data for lt divisions"
        else:
            # Initialize figure
            f.append(plt.figure())
            if nsub > 1:
                fheight = f[-1].get_figheight()
                f[-1].set_figheight(fheight * nsub)

            htitle = "Local Time Regions"
            if title is not None:
                htitle = "{:s}\n{:s}".format(htitle, title)
            f[-1].suptitle(htitle)

            # Plot each subfigure
            for i, dat_index in enumerate(lt_index):
                ax = f[-1].add_subplot(nsub,1,i+1)

                try: nweights = weights[dat_index]
                except: nweights = None

                plot_hist_w_stats(ax, x_data[dat_index], x_name, x_units,
                                  res=res, norm=norm, weights=nweights, cum=cum,
                                  htype=htype, align=align, rwidth=rwidth,
                                  n_label=lt_label[i], xmin=xmin, xmax=xmax,
                                  nmin=nmin, nmax=nmax, xinc=xinc, ninc=ninc,
                                  color=color, tloc=tloc, xl=xl, xt=xt, yl=yl,
                                  yt=yt, moments=moments, quartiles=quartiles,
                                  modes=modes, sloc=sloc)

            # Adjust subplot locations
            plt.subplots_adjust(left=.15)
            if draw and plt.isinteractive():
                plt.draw() #In interactive mode, you just "draw".

            if figroot:
                figname = figroot+'.LT.png'
                plt.savefig(figname)

    if draw and not plt.isinteractive():
        # W/o interactive mode, "show" stops the user from typing more 
        # at the terminal until plots are drawn.
        plt.show()

    return f




