#!/usr/bin/env python

#this magnetogram reading can either be run as a script from the unix command 
#line or imported into Python
#
# To install astropy:
#      Go directly to https://pypi.python.org/pypi/astropy
#      Or, go to www.astrophy.org => Install from Source => PyPi
#      Download the tarball of astropy****.tar
#      Untar the tarball.
#      Go to directory astropy****/
#      python setup.py install
#

from astropy.io import fits
import numpy as np
import argparse


def readf(NameFile,BMax):
    """
    NameFile - FITS file containing original magnetogram (include path)
    """

    TypeMag = 'unknown'
    g = fits.open(NameFile)
    try:
        NameTelescope = g[0].header['TELESCOP'] #works for MDI, GONG
    except KeyError, er:
        NameTelescope = 'unknown'
        
    nLong = g[0].header['NAXIS1'] # number of longitude points
    nLat = g[0].header['NAXIS2']  #           latitude      
    Br_II = g[0].data 

    if NameTelescope.find('NSO-GONG') > -1 :
        TypeMag = 'GONG Synoptic'
        try:
            Long0 = g[0].header['LONG0']
            if float(Long0) > 0.:
                TypeMag = 'GONG Hourly'
        except KeyError, er:
            Long0 = - 1

    CRNumber = str(g[0].header['CAR_ROT'])

    if TypeMag.find('GONG') > -1:
        MapDate = g[0].header['DATE'] #works for GONG
    BUnit = g[0].header['BUNIT']  #works on GONG, MDI
                 
    g.close()
    if  ( (TypeMag == 'unknown') ):
        print "I don't recognize the type of this magnetogram."
        return(-1)
                
    print "This is a",TypeMag,"magnetogram on a",str(nLong),"X",str(nLat),"  Phi X sin(Theta) grid."
    if (Long0 > 0):  #The date is not informative for an integral synoptic map
        print "Magnetogram Date="+MapDate

    #Uniform in sinLat and longitude grid
    LatSin_I = np.arcsin(np.linspace(-1. + 1./2/nLat,1. - 1./2/nLat,nLat))
    Long_I = 2.*np.pi*np.linspace(0.,1. - 1./nLong, nLong)
    BrTransp_II = np.zeros([nLong,nLat])
 
    FileId = open('CRLong.dat','w')
    FileId.write(CRNumber + '\n')  
    FileId.write(str(Long0) + '\n')
    FileId.close()

    FileId = open('fitsfile.dat','w')

    FileId.write('#nMax\n')
    FileId.write(str(180) + '\n')
    FileId.write('#ARRAYSIZE\n')
    FileId.write(str(nLong) + '\n')
    FileId.write(str(nLat) + '\n')
    FileId.write('#START\n')
    for k in np.arange(nLat):
        for l in np.arange(nLong):
            FileId.write(str(Br_II[k,l]) + '\n')
            BrTransp_II[l,k] = Br_II[k,l]
    FileId.close()


    FileId = open('fitsfile_idl.out','w')
    
    Name0 = 'magnetogram type = '+TypeMag+', MapDate = '+MapDate+', Br['+BUnit+']'+'\n' 
    FileId.write(Name0)
    Name0 = '       0      0.00000       2       1       1 \n'
    FileId.write(Name0)
    FileId.write('      '+str(nLong)+'     '+str(nLat)+'\n')
    FileId.write(str(Long0)+'\n') #longitude shift (important for GONG Hourly)
    FileId.write('Longitude Latitude Br LongitudeShift \n')
    
    for k in np.arange(nLat):
         for l in np.arange(nLong):
             Name0 = str(Long_I[l]) + '   ' + str(LatSin_I[k]) + '   ' + str(max([-BMax,min([BMax,Br_II[k,l]])])) + ' \n'
             FileId.write(Name0)
    FileId.close()
   
    return(BrTransp_II,Long0,nLat,nLong,LatSin_I,Long_I)
        
if __name__ == '__main__':

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description="""
    pre-processes the FITS format magnetograms 
    into ASCII files. The script can read the following types of magnetograms:

       GONG Synoptic
       GONG Hourly updated
       (to be extended)
    The code opens the .fits file and automatically recognizes the type of
    map it is.
    """)
    parser.add_argument('NameFile', help='Input FITS file name including path')
    parser.add_argument('-BMax',type=float, default=1900., help='Max BrFieldAmplitude')

    args = parser.parse_args()

    readf(args.NameFile,args.BMax)

    

        

    
