#!/usr/bin/env python
#---------------------------------------------------------------------------
# DownloadMadrigalTEC
#
# Author: Angeline G Burrell, UMichigan, Dec 2012
#
# Comments: A routine to download TEC files from the Madrigal database
#---------------------------------------------------------------------------

import string
import sys
import MadrigalData_Rout as mdr

module_name = "DownloadMadrigalTEC"

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

mSite = 'http://madrigal.haystack.mit.edu'
eName = 'World-wide GPS Receiver Network'

# Download the desired experiment files from any available Madrigal database
# (local is set to false, 0)

fnum  = mdr.downloadMadrigalFiles(mSite, eName, syear, smonth, sday, shour,
                                  smin, ssec, eyear, emonth, eday, ehour, emin,
                                  esec, 0, dest_dir, user_name, user_email,
                                  user_affil, file_format)

print "%s downloaded [%d] files\n" % (module_name, fnum)
