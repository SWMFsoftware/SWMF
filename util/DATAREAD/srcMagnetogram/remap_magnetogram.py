#!/usr/bin/env python

# this magnetogram remapping can either be run as a script from 
# the unix command line
# or imported into Python
# accompanying files (must be in same directory):
#    remap_magnetogram.py.README.txt 
#         by  Richard A. Frazin July-Sept 2104

from astropy.io import fits
from scipy import interpolate
from scipy import integrate
import numpy as np
import sys
import time
import argparse

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description="""
remap_magnetogram.py pre-processes the FITS format magnetograms 
into ASCII files that can by read by FDIPS.exe, BATSRUS.exe and SWMF.exe
and IDL macros. The script can read the following types of magnetograms:

   Hathaway synchronic
   ADAPT synchronic
   GONG synoptic
   GONG hourly updated
   MDI synoptic

In the case of the GONG and MDI maps, these routines only convert the
existing .fits files into the standard ASCII format with no remapping.
In case of the ADAPT and the Hathaway map, the code allows arbitrary
output resolution, specified by [nlat] and [nlon] with defaults being
180 and 360.  The code opens the .fits file and automatically
recognizes the type of map it is.  In all cases the output is a
sin(latitude) grid.  

The script uses the scipy and astropy packages that can be installed, 
for example, with MacPorts.
""")
parser.add_argument('inputfile', help='Input FITS file')
parser.add_argument('outputfile', help='Output magnetogram file')
parser.add_argument('nlat', nargs='?', type=int, default=180, help='Number of latitudes')
parser.add_argument('nlon', nargs='?', type=int, default=360, help='Number of longitudes')

args = parser.parse_args()

def remap(inputfile, outputfile, nlat = 180, nlong = 360):
    """
    This takes an ADAPT or Hathaway magnetogram on a regular spherical
    inputfile includes the path.  outputfile requires the desired suffix.
    nlat and nlong are desired number of sin(lat) point and longitude points
       in the output grid.
    Return values are new and old maps. 
    This uses an interpolation/integration scheme to conserve magnetic flux.
    Note that ADAPT files may have multiple maps.  Only the 1st is utilized.
    My MATLAB code remap_mag.m uses the same algorithm but runs much faster
      and conserves the flux 1 part in 10^5 (on Hathaway magnetograms),
      whereas here it's not very close.  Is this due to differences in the
      integrator (or maybe even the interpolator)?  
    """
    hdulist = fits.open(inputfile)
    d = hdulist[0].data
    try:
        magtype = hdulist[0].header['MODEL']
    except KeyError, err:
        try:
            magtype = hdulist[0].header['SFT_TYP']
            if magtype == 'Baseline / Assimilation' :
                magtype = 'Hathaway'
        except KeyError, er:
            magytpe = 'Unknown'

    if magtype == 'Unknown' :
        print 'Error. Unknown magnetogram type.'
        return -1
            
    nlo = hdulist[0].header['NAXIS1'] # number of longitude points
    nla = hdulist[0].header['NAXIS2'] #           latitude
    hdulist[0].header['NAXIS1'] = nlong # new number of longitude points
    hdulist[0].header['NAXIS2'] = nlat  #               latitude
    hdulist[0].header.set('REMAP','lat -> sin(lat)')  #make a new keyword for the header
    
    if magtype == 'Hathaway':
         print 'I think this is a Hathaway synchronic map. native size: ',str(nla),' X ',str(nlo)
         hdulist[0].header.set('GRID','sin(latitude)')

    if magtype == 'ADAPT':
        nim = hdulist[0].header['NAXIS3'] # number of images
        imdex = 5  #which of the 12 maps do you want?
        print 'I think this is an ADAPT synchronic map.  This file contains ', str(nim), ' images.  Using only number ' + str(imdex)
        print 'native size: ',str(nla),' X ',str(nlo)
        gridtype = hdulist[0].header['GRID'] # = 1. for a uniform spherical grid
        if gridtype != 1.0:
            print('rempap_magnetogram.py: the keyword \'GRID\' is not 1.   Either this is not an ADAPT map or this data is not on a uniform spherical mesh.')
            return(-1)
        gridtype = 'sin(latitude)'
        hdulist[0].header['GRID'] = gridtype #change the header value
        if nim > 1:  #just keep one of them for now
            d = d[imdex,:,:]


    #if nlo != nlong, first make a hybrid map that is (nla X nlong), otherwise, do nothing
    if nlo == nlong:
        hybrid = d
    else: #find the common factors of nlo and nlong for the rebin and add integration aglorithm
        ## pf = PrimeFactors(nlo)
        ## nlo_fac = np.asarray(pf.compute_prime_factors())
        ## pf = PrimeFactors(nlong)
        ## nlong_fac = np.asarray(pf.compute_prime_factors())
        ## pf = np.array(1) #array of common factors
        ## for k in np.arange(nlong_fac.size):
        ##     crap = np.where(nlong_fac[k] == nlo_fac); crap = crap[0]
        ##     if crap.size != 0: #they don't have this factor in common
        ##         pf = np.append(pf,nlong_fac[k]) # put it in the common factor array
        ##         nlong_fac[k] = 1 # replace it with 1 in the other arrays
        ##         nlo_fac[np.min(crap)] = 1 #

        ## nlong_fac = np.prod(nlong_fac) #nlong = nlong_fac*pf and nlo = nlo_fac*pf
        ## nlo_fac   = np.prod(nlo_fac)   #  thus, the least common mulitple is 
        ## pf        = np.prod(pf)        #      pf*nlo_fac*nlong_fac
        #print 'common factor: ' + str(pf) + ' upsample factor: ' + str(nlong_fac) + ' rebin factor: ' + str(nlo_fac)

        #this way is simpler ... duh!
        crap = np.arange(nlo)
        for pf in crap[::-1]:
            if ( (np.mod(nlo,pf) == 0) and (np.mod(nlong,pf) == 0)): #greatest common factor
                nlo_fac   = nlo/pf
                nlong_fac = nlong/pf
                break

        hybrid = np.zeros([nla,nlong]) #intermediate hybrid map
        for k in np.arange(nla):
            w = np.kron(d[k,:],np.ones(nlong_fac)) #this array has length pf*nlo_fac*nlong_fac
            for l in np.arange(nlong):  #take the average over nlo_fac bins of w
                hybrid[k,l] = np.sum(w[l*nlo_fac:(l+1)*nlo_fac])/nlo_fac

    newmap = np.zeros([nlat,nlong])
    
    pi = 3.141592653589793
    newlat = np.zeros(nlat)
    oldlat =  np.linspace(-pi/2 + pi/2/nla,pi/2 - pi/2/nla,nla) #old latitude grid
    oldlat = np.hstack((-pi/2-1.e-9,oldlat,pi/2+1.e-9)) #for the interpolator
    bin_boundary = np.linspace(-1.,1.,nlat+1)

    #the magnetic field value assigned is the flux divided by the area.    
    for k in np.arange(nlong):
        u = np.hstack((hybrid[0,k],hybrid[:,k],hybrid[nla-1,k]))
        crap = interpolate.interp1d(oldlat,u,kind='linear')
        fcn = lambda x : crap(x)*np.cos(x)  #x is latitude in radians
        for l in np.arange(nlat):
            newlat[l] = 0.5*(np.arcsin(bin_boundary[l]) + np.arcsin(bin_boundary[l+1]))*180./pi
            result = integrate.quad(fcn,np.arcsin(bin_boundary[l]),np.arcsin(bin_boundary[l+1]),epsabs=1.e-3,epsrel=1.e-3)/(bin_boundary[l+1] - bin_boundary[l])
            newmap[l,k] = result[0]

    #ascii output file
    fid = open(outputfile,'w')


    if magtype == 'ADAPT' :
        try:
            CRnumber = str(hdulist[0].header['CRROTEDG'])
        except KeyError, er:
            CRnumber = '0'
        try:
            mapdate = hdulist[0].header['MAPTIME']
        except KeyError, er:
            mapdate = '0000-00-00T00:00:00'

    if magtype == 'Hathaway' :
        CRnumber = '0'
        try :
            mapdate = hdulist[0].header['MAP_DATE']
        except KeyError, er:
            mapdate = '0000-00-00T00:00:00'
    
    #new output format
    #the first line is arbitary
    line0 = 'magnetogram type = '+magtype+', grid_type = sin(lat), CR'+CRnumber+ ', MapDate = '+mapdate+', units: [Deg], [G], created at: '+time.ctime()+'\n' 
    fid.write(line0)
    line0 = '       0      0.00000       2       1       1 \n'
    fid.write(line0)
    fid.write('      '+str(nlong)+'     '+str(nlat)+'\n')
    fid.write(str(0.5*360./nlong) + ' \n') #longitude shift 
    fid.write('Longitude Latitude Br LongitudeShift \n')
    
    for k in np.arange(nlat):
         for l in np.arange(nlong):
             line0 = str(l*360./nlong) + ' ' + str(newlat[k]) + ' ' + str(newmap[k,l]) + ' \n'
             fid.write(line0)

    #old output format
    ## fid.write('#CR\n')
    ## try:
    ##     fid.write(str(hdulist[0].header['CRROTEDG']) + '\n')
    ## except KeyError, er:
    ##     fid.write('-1')
    ## fid.write('#nMax\n')
    ## fid.write(str(-1) + '\n')
    ## fid.write('#ARRAYSIZE\n')
    ## fid.write(str(nlong) + '\n')
    ## fid.write(str(nlat) + '\n')
    ## fid.write('#START\n')
    ## for k in np.arange(nlat):
    ##     for l in np.arange(nlong):
    ##         fid.write(str(newmap[k,l]) + '\n')

    hdulist.close()
    fid.close()
    return(newmap,d)


def GONGFITStoASCII(inputfile,outputfile):
    """
    The GONG and MDI fits files are on a sin(latitude) grid.  All that needs
    to be done is to rewrite the data in ascii format.
    """

    g = fits.open(inputfile)
    nlong = g[0].header['NAXIS1'] # number of longitude points
    nlat  = g[0].header['NAXIS2'] #           latitude
    try:
        telescope = g[0].header['TELESCOP']
    except KeyErr, er:
        telescope = 'other'
                
    try:
        inst = g[0].header['INSTRUME']
    except KeyError,er:
        inst = 'other'

    try:
        ctyp = g[0].header['CTYPE2']
    except KeyError, er:
        ctyp = 'other'

    if telescope != 'NSO-GONG' :    
        if ( (ctyp != 'Sine Latitude') or (inst != 'MDI') ):
            print 'not a GONG magetogram, not getting the Keywords expected from an MDI magnetogram'
            return(-1)
  
    if ( (telescope != 'NSO-GONG') and (telescope != 'SOHO') ):
        print 'Not a GONG or MDI magnetogram!'
        return(-1);

    newlat = np.arcsin(np.linspace(-1.+1./nlat , 1.-1./nlat ,nlat))*180/3.141592653589793

    fid = open(outputfile,'w')
    
    try :
        CR = str(g[0].header['CAR_ROT'])
    except KeyError,er:
        CR = '0'
    
    try :
        mapdate = g[0].header['DATE']
    except KeyError, er:
        mapdate = '0000-00-00T00:00:00'

    try :
        long0 = g[0].header['LONG0']
    except KeyError, er:
        long0 = 0

    #new output format
    #the first line is arbitary
    line0 = 'magnetogram type = '+telescope+', grid_type = sin(lat), CR'+CR+ ', MapDate = '+mapdate+', units: [Deg], [G], created at: '+time.ctime()+'\n' 
    fid.write(line0)
    line0 = '       0      0.00000       2       1       1 \n'
    fid.write(line0)
    fid.write('      '+str(nlong)+'     '+str(nlat)+'\n')
    fid.write(str(long0+0.5)+' \n') #longitude shift 
    fid.write('Longitude Latitude Br LongitudeShift \n')
    
    for k in np.arange(nlat):
         for l in np.arange(nlong):
             line0 = str(l*360./nlong) + ' ' + str(newlat[k]) + ' ' + str(g[0].data[k,l]) + ' \n'
             fid.write(line0)
       

    #old output format
    ## fid.write('#CR\n')
    ## try:
    ##     fid.write(str(g[0].header['CAR_ROT']) + '\n')
    ## except KeyError, er:
    ##     fid.write('-1')
    ## fid.write('#nMax\n')
    ## fid.write(str(-1) + '\n')
    ## fid.write('#ARRAYSIZE\n')
    ## fid.write(str(nlong) + '\n')
    ## fid.write(str(nlat) + '\n')
    ## fid.write('#START\n')

    ## for k in np.arange(nlat):
    ##     for l in np.arange(nlong):
    ##         fid.write(str(g[0].data[k,l]) + '\n')

    g.close()
    fid.close()
    return(0)
    

    
if __name__ == '__main__':

    if len(sys.argv) < 3:
        print 'format: python remap_magnetogram.remap inputfile ouputfile [nlat] [nlong]'
        exit()

    inputfile = sys.argv[1] # sys.argv[0] is the 'remap_magnetogram.py' command!
    outputfile = sys.argv[2]
    print 'inputfile is: ', inputfile
    print 'outputfile will be: ', outputfile

    #check to see if it's a GONG or MDI magnetogram
    g = fits.open(inputfile)
    try:
        telescope = g[0].header['TELESCOP']
    except KeyError, er:
        telescope = 'other'

    try:
        inst = g[0].header['INSTRUME']
    except KeyError,er:
        inst = 'other'

    try:
        ctyp = g[0].header['CTYPE2']
    except KeyError, er:
        ctyp = 'other'

    g.close()
        
    if ( (ctyp == 'Sine Latitude') and (inst == 'MDI') ):
        telescope = 'NSO-GONG'
        
    if telescope == 'NSO-GONG':
        GONGFITStoASCII(inputfile,outputfile)
    else:
        if len(sys.argv) > 3:
            if len(sys.argv) != 5:
                print 'format: python remap_magnetogram.py inputfile ouputfile [nlat] [nlong]'
                exit()
            nlat = int(sys.argv[3])
            nlong = int(sys.argv[4])
            print 'the output grid have ', str(nlat), 'latitude points and ', str(nlong), ' longitude points.'
        if len(sys.argv) > 3:
            remap(inputfile,outputfile,nlat,nlong)
        else:
            remap(inputfile, outputfile)
