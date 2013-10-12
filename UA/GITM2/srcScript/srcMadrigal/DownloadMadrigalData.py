#!/usr/bin/env python
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf
#---------------------------------------------------------------------------
# DownloadMadrigalData
#
# Author: Angeline G Burrell, UMichigan, Dec 2012
#         Nick Perlongo, UMichigan, Aug 2013
#
# Comments: A routine to download files from the Madrigal database
#---------------------------------------------------------------------------

import string
import sys
import MadrigalData_Rout as mdr
import madrigal_info_rout as mir 
import madrigalWeb.madrigalWeb as mWeb


module_name = "DownloadMadrigalData"

#---------------------------------------------------------------------------
# Get the desired starting and ending times

syear  = int(raw_input('Start Year: '))
smonth = int(raw_input('Start Month: '))
sday   = int(raw_input('Start Day: '))
shour  = int(raw_input('Start Hour: '))
smin   = int(raw_input('Start Minute: '))
ssec   = int(raw_input('Start Second: '))

eyear  = raw_input('End Year [%d]: ' % syear)
eyear  = int(eyear) if eyear else syear
emonth = raw_input('End Month [%d]: ' % (smonth + 1))
emonth = int(emonth) if emonth else smonth + 1
eday   = raw_input('End Day [%d]: ' % (sday + 1))
eday   = int(eday) if eday else sday + 1
ehour  = raw_input('End Hour [%d]: ' % shour)
ehour  = int(ehour) if ehour else shour
emin   = raw_input('End Minute [%d]: ' % smin)
emin   = int(emin) if emin else smin
esec   = raw_input('End Second [%d]: ' % ssec)
esec   = int(esec) if esec else ssec

#-----------------------------------------------------------------------------
# Get the file destination and Madrigal access input

dest_dir    = raw_input('Local Download Directory: ')
user_name   = raw_input('Madrigal User Name: ')
user_email  = raw_input('Madrigal User Email: ')
user_affil  = raw_input('Madrigal User Affiliation: [University of Michigan]')
user_affil  = user_affil if user_affil else "University of Michigan"
file_format = raw_input('Desired file format: [simple]')
file_format = file_format if file_format else "simple"

# Specify the TEC variables and download the files

mcode = dict()
mname = dict()
mcat = dict()

# Determine which instruments are available for selected paramaters

print '\n .............Retreiving available instruments................ \n'

mcode, mname, mcat, msites = mir.get_possible_instruments(syear,smonth,sday,shour,smin,ssec,
                             eyear,emonth,eday,ehour,emin,esec)

#Tell them which sites are available and ask for them to choose one
print '\nThe sites available are: \n'
count=0
for key in mname.iterkeys():
    print '%d - %s' % (count,key)
    count +=1

mSiteIn = raw_input('\n>>Choose a site[%s]: ' % '0')
if not mSiteIn:
    mSite = mname.keys()[0]
else:
    mSite =mname.keys()[int(mSiteIn)]

mSiteUrl = msites[mSite] #need the URL, not the name

#Now they need to choose the experiment from that site

print '\nThe experiments available are: \n'
for i in range(len(mname[mSite])):
    print '%d - %s' % (i, mname[mSite][i])

eIn = raw_input('\n>>Choose an experiment[%s]: ' % '0')
if not eIn:
    eName = mname[mSite][0]
else:
    eName = mname[mSite][int(eIn)]

# Download the desired experiment files from any available Madrigal database
# (local is set to false, 0)

print '\n ...........Downloading your data............ \n'

fnum  = mdr.downloadMadrigalFiles(mSiteUrl, eName, syear, smonth, sday, shour,
                                  smin, ssec, eyear, emonth, eday, ehour, emin,
                                  esec, 0, dest_dir, user_name, user_email,
                                  user_affil, file_format)

print "%s downloaded [%d] files\n" % (module_name, fnum)
