#!/usr/bin/env python
#---------------------------------------------------------------------------
# write_files
#
# Author: Angeline G Burrell, UMichigan, Feb 2014
#
# Comments: Routines to support file processing in python
#
# Contains: writeASCII_file
#---------------------------------------------------------------------------

import string
import numpy as np
from os import path

module_name = "write_files"

#---------------------------------------------------------------------------
# writeASCII_file: A routine to create an ascii file from a string or list
#                  of strings.  Will overwrite a file with the same name if
#                  such a file exists.

def writeASCII_file(filename, datalines, *args, **kwargs):
    '''
    A routine to create an ascii file from a string or list of strings.  Will
    overwrite a file with the same name if such a file exists

    Input:
    filename  = output file name
    datalines = string or list of strings to be written
    '''

    func_name = string.join([module_name, "writeASCII_file"], " ")

    #-----------------------------------------------------------------------
    # Open and test the file to ensure it can be written
    try:
        f = open(filename, 'w')

        # Print data after determinging the appropriate data type
        try:
            dlen = len(datalines)

            if type(datalines) is str:
                f.write(datalines)
            else:
                for line in datalines:
                    line = line+"\n"
                    f.write(line)
        except:
            f.write(datalines)

        f.close()
    except:
        print func_name, "ERROR: unable to open [", filename, "]"

# End writeASCII_file
