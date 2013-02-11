#!/usr/bin/env python
#---------------------------------------------------------------------------
# DownloadMadrigalTEC
#
# Author: Angeline G Burrell, UMichigan, Dec 2012
#
# Comments: Routines to support file downloads from the Madrigal database
#---------------------------------------------------------------------------

import madrigalWeb.madrigalWeb as mWeb
import string
import os

module_name = "MadrigalData_Rout"

def getInstrumentCode(mData, dataName, *args, **kwargs):
    '''
    Obtain the code for a given instrument name.  An open connection to the
    desired Madridal database and the instrument name must be provided.
    '''

    #------------------------------
    # Get the list of instruments

    mInst = mData.getAllInstruments()
    mCode = list()
    mName = list()

    #----------------------------------------------------------------------
    # Extract the instrument name and code from the instrument information
    # for all instruments with the desired keyword in the name
    
    for inst in mInst:
        if string.find(inst.name, dataName) >= 0:
            mCode.append(inst.code)
            mName.append(inst.name)

    return(mCode, mName)

# End getInstrumentCodes

def getFileList(mData, mCode, syear, smonth, sday, shour, smin, ssec, eyear,
                emonth, eday, ehour, emin, esec, local, *args, **kwargs):
    '''
    Get the experiment files for a particular instrument (specified by the code)
    and date range.  To search the local Madrigal database, let local = 1,
    otherwise set local = 0.
    '''

    ##
    # Get a list of experiment IDs for the specified instrument (indicated
    # by the instrument code, mCode) for the specified date range.

    mId = mData.getExperiments(mCode, syear, smonth, sday, shour, smin, ssec,
                               eyear, emonth, eday, ehour, emin, esec, local)

    mFiles = list()

    ##
    # Retireve the experiment file information for experiments with a
    # valid ID (ID is not set to -1)

    for ident in mId:
        if ident.id != -1:
            mFiles += mData.getExperimentFiles(ident.id, False)

    return(mFiles)

# End getFileList

def downloadMadrigalFiles(mSite, eName, syear, smonth, sday, shour, smin, ssec,
                          eyear, emonth, eday, ehour, emin, esec, local,
                          dest_dir, user_name, user_email, user_affil,
                          file_format, *args, **kwargs):
    '''
    Download Madrigal files from a specified Madrigal database from a 
    specified experiment for a specified date range.  The user information
    must also be specified so that Madrigal can know who is using what parts
    of the database.  Acceptable file formats include "simple (ascii)", hdf5,
    and several binary formats not recommended for general use.
    '''

    ##
    # Connect to the MIT Haystack Madrigal Database

    mData = mWeb.MadrigalData(mSite)

    ##
    # Get the file code(s) for the desired experiment(s)

    mCode, mName = getInstrumentCode(mData, eName)
    mFiles       = list()

    ##
    # Get the file information for the desired code(s) for the specified
    # range of dates.

    for codes in mCode:
        mFiles += getFileList(mData, codes, syear, smonth, sday, shour, smin,
                              ssec, eyear, emonth, eday, ehour, emin, esec,
                              local)

    if len(mFiles) <= 0:
        ##
        # Send a warning if there are no files for this instrument in the
        # desired date range
        print "%s WARNING: no [%s] files available between " % (module_name,
                                                                hName)
        print "[%d/%d/%d %d:%d:%d] and [%d/%d/%d %d:%d:%d]\n" % (smonth, sday,
                                                                 syear, shour,
                                                                 smin, ssec,
                                                                 emonth, eday,
                                                                 eyear, ehour,
                                                                 emin, esec)
        return(0)
    else:
        ##
        # Download the TEC files
        dnum = 0
        
        for fileDat in mFiles:
            # Extract the filename from the remote file information
            filelist = string.split(fileDat.name, "/")
            filename = filelist.pop()
            # Construct the destination filename
            destfile = "%s/%s" % (dest_dir, filename)
            if not os.path.exists(destfile):
                # If the destination file does not already exist, download
                # the file from Madrigal
                mData.downloadFile(fileDat.name, destfile, user_name,
                                   user_email, user_affil, file_format)

                # Test file size to see if there was an error

                statinfo = os.stat(destfile)

                if statinfo.st_size < 1000:
                    # When there is an error downloading a file, a file with
                    # the destination filename is created containing an error
                    # message.  This file is typically 910 kB so test for this
                    # size and (remove)/warn the user
                    print "%s WARNING: deleting small file [%s]\n" % (module_name, destfile)
                    #os.remove(destfile)
                else:
                    # Count the successfully downloaded files
                    dnum += 1

        return(dnum)
    
# End downloadMadrigalFiles


