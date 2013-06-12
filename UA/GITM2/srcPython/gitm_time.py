#!/usr/bin/env python
#-----------------------------------------------------------------------------
# $Id$
#
# GITM_Time.py, Angeline Burrell (AGB), UMich, June 2013
#
# Comments: Defines a class to hold data from multiple GITM binary output files,
#           allowing UT dependence to be explored
#
# Contains: class GitmTime - The class for the GITM binary, which will contain
#                            data from multiple GITM output binaries
#           def load_multiple_gitm_bin - A routine to load multiple GITM binary
#                                        files, returning a list of the
#                                        GitmBin data structures
#------------------------------------------------------------------------------

'''
PyBats submodel for handling input/output for the Global Ionosphere-Thermosphere
Model (GITM), providing data structures with universal time dependence.
'''

# Global Imports

from spacepy.pybats import PbData
from spacepy.datamodel import dmarray
from spacepy.pybats import gitm
import string
import copy

# Temporary
import gitm

class GitmTime(PbData):
    '''
    Object containing GITM data from multiple GTIM binaries, providing a 
    data structure with UT time dependence.
    '''

    def __init__(self, filelist, *args, **kwargs):
        if len(args) > 0:
            ionfile = args[0]
            newarg = args[1:len(args)]
            args = copy.deepcopy(newarg)
            del newarg
        else:
            ionfile = None

        super(GitmTime, self).__init__(*args, **kwargs) # Init as PbData
        self.attrs['nFiles'] = 0
        self.attrs['list'] = filelist

        if ionfile:
            self.attrs['ionfile'] = ionfile
        
        self._appendgitm()
        
    def __repr__(self):
        return 'File with list of GITM binaries: %s' % (self.attrs['file'])

    def _appendgitm(self):
        '''
        Append GITM binary file into a single data structure
        '''

        # Import local packages
        from datetime import datetime
        import numpy as np

        # Load the gitm binaries
        GitmList = load_multiple_gitm_bin(self.attrs['list'])
        self.attrs['nFiles'] = len(GitmList)

        #Initialize and fill the new data structure
        for i,gData in enumerate(GitmList):
            # If there is an ion file specified, compute the magnetic velocities
            if self.attrs.has_key('ionfile'):
                gData.attrs['ionfile'] = self.attrs['ionfile']
                gData.calc_magvel()

            # Determine which keys are new and which already exist
            newkeys = gData.keys()
            oldkeys = list()

            if i > 0:
                oldkeys = self.keys()

                for old in oldkeys:
                    try:
                        n = newkeys.index(old)
                        newkeys.pop(n)
                    except ValueError:
                        print "ADVISEMENT: file [",i+1,"] is missing [",old,"]"
                        

            # Initialize the first instance of each key and append the new data
            # if it is not
            for k in newkeys:
                if type(gData[k]) is dmarray:
                    scale = 0.0
                    for j in reversed(range(i+1)):
                        if j == 0:
                            scale = 1.0
                        if j == i:
                            data = np.array([gData[k]*scale])
                        else:
                            data = np.append(data, [gData[k]*scale], 0)
                elif i != 0:
                    print "WARNING: key [",k,"] not temporally aligned"
                else:
                    data = np.array([gData[k]])

                if type(gData[k]) is gitm.dmarray:
                    self[k] = gitm.dmarray(data, gData[k].attrs)
                elif type(gData[k]) is datetime:
                    self[k] = gitm.dmarray(data, attrs={"name":"Universal Time",
                                                        "units":"date",
                                                        "scale":"date"})
                else:
                    self[k] = gitm.dmarray(data)
            
            for k in oldkeys:
                if gData.has_key(k):
                    self[k] = gitm.dmarray(np.append(self[k], [gData[k]], 0),
                                           self.attrs)
                else:
                    self[k] = gitm.dmarray(np.append(self[k], [self[k][0]*0.0],
                                                     0, self.attrs))

#End Class

def load_multiple_gitm_bin(filelist, *args, **kwargs):
    '''
    Loads a list of GITM binary files into their own GitmBin data structures.
    A list of the data structures is returned.
    '''
    # import local packages
    from os import path

    func_name = "load_multiple_gitm_bin"
    fsize = path.getsize(filelist)
    outlist = list()

    if(fsize > 2.0e9):
        print func_name, "ERROR: File list size [", (fsize * 1e-9), "GB > 2 GB]"
        return outlist
    elif(fsize == 0.0):
        print func_name, "ERROR: empty file list [", filelist, "]"
        return outlist

    # Read in the list of files
    f = open(filelist, "r")
    namelist = f.readlines()

    for name in namelist:
        name = string.strip(name)
        outlist.append(gitm.GitmBin(name))

    return(outlist)
# End load_multiple_gitm_bin
