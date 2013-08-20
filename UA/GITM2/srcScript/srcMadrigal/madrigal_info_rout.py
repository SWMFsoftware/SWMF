#!/usr/bin/env python
#---------------------------------------------------------------------------
# madrigal_info_rout
#
# Author: Angeline G Burrell, UMichigan, July 2013
#
# Comments: Routines to get information from Madrigal databases
#---------------------------------------------------------------------------

import os
import string
import sys
import madrigalWeb.madrigalWeb as mWeb

module_name = "madrigal_info_rout"

#---------------------------------------------------------------------------
# Define the available madrigal sites

update = "July 12, 2013"
msites = {'haystack':'http://madrigal.haystack.mit.edu',
          'CEDAR':'http://cedar.openmadrigal.org',
          'ARO':'http://madrigal.naic.edu',
          'EISCAT':'http://www.eiscat.se/madrigal',
          'SRI':'http://isr.sri.com/madrigal',
          'Cornell':'http://landau.geo.cornell.edu/madrigal',
          'JRO':'http://jro1.igp.gob.pe/madrigal',
          'Beijing':'http://madrigal.iggcas.ac.cn/madrigal',
          'SGO':'http://madrigal.sgo.fi'}

#-----------------------------------------------------------------------------
# get_possible_instruments: retrieve the instrument names and codes from
#                           a specified Madrigal site or all sites

def get_possible_instruments(startyear, startmonth, startday, starthour,
                             startmin, startsec, endyear, endmonth, endday,
                             endhour, endmin, endsec, msite=None, *args,
                             **kwargs):
    '''
    Obtain the names and codes for instruments at a specified Madrigal site for
    a specified range of dates. If no site is specified, instrument names and
    codes from all sites will be retrieved.  The site is specified by a name
    code:

    haystack - Millstone Hill
    CEDAR - CEDAR Archive
    ARO - Arecibo
    EISCAT - EISCAT
    SRI - SRI International
    Cornell - Cornell University
    JRO - Jicamarca Radio Observatory
    Beijing - Institute of Geodesy and Geophysics, Chinese Academy of Sciences
    SGO - University of Oulu

    Input:
    startyear  = 4 digit year of starting date
    startmonth = 1-12 integer month
    startday   = integer day of month
    starthour  = 0-23 integer hour of day
    startmin   = 0-60 minutes of hour
    startsec   = 0-60 seconds of minute
    endyear    = 4 digit year of ending date
    endmonth   = 1-12 integer month
    endday     = integer day of month
    endhour    = 0-23 integer hour of day
    endmin     = 0-60 minutes of hour
    endsec     = 0-60 seconds of minute
    msite      = Madrigal site code or None (default)
    '''

    print module_name, "ADVISEMENT: Madrigal sites last updated on", update 

    mkeys = list()
    icode = dict()
    mcode = dict()
    mname = dict()
    mcat = dict()

    #------------------------------------------------
    # Determine which Madrigal database(s) to search

    if msite:
        if not msites.has_key(msite):
            print module_name, "ERROR: unknown Madrigal site [", msite, "]"
            return(mcode, mname, mcat)

        mkeys.append(msite)
        mcode[msite] = list()
        mname[msite] = list()
        mcat[msite] = list()
    else:
        for mkey in msites:
            mkeys.append(mkey)
            mcode[mkey] = list()
            mname[mkey] = list()
            mcat[mkey] = list()

    #-----------------------------------------------------------------------
    # Go through each site and determine which instruments have experiments
    # available in the desired time range
    
    for m in msites:
        mdata = mWeb.MadrigalData(msites[m])

        #-----------------------------------------------------------------
        # Get the list of instruments, saving the names, codes, and site
        #  keys.  These names are the same for all the sites, but the 
        # experiments available at each site are different.

        if len(icode.keys()) == 0:
            minst = mdata.getAllInstruments()

            # Resort instrument data by code
            for inst in minst:
                icode[inst.code] = inst

        #----------------------------------------------------------
        # Determine which experiments are available at each site

        mexp = mdata.getExperiments(0, startyear, startmonth, startday,
                                    starthour, startmin, startsec, endyear,
                                    endmonth, endday, endhour, endmin, endsec)
        new_exp = dict()

        

        #---------------------------------------------------------------
        # Save the code, name, and instrument category for each unique
        # experiment stored locally at this Madrigal database

        for exp in mexp:
            if not new_exp.has_key(exp.instcode):
                mcode[m].append(exp.instcode)
                mname[m].append(exp.instname)
                mcat[m].append(icode[exp.instcode].category)
                new_exp[exp.instcode] = 1

    #----------------------------------------------------------------------
    # Return the instrument code, instrument name, and instrument category
    # in dictionaries grouped by Madrigal site
    return(mcode, mname, mcat,msites)

# End get_instrument_codes

