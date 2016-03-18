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
import os

cPi = np.pi
def readf(NameFile,TypeOut,nSmooth,BMax):
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

    CRNumber = g[0].header['CAR_ROT']

    if TypeMag.find('GONG') > -1:
        MapDate = g[0].header['DATE'] #works for GONG
    BUnit = g[0].header['BUNIT']  #works on GONG, MDI
                 
    g.close()
    if  ( (TypeMag == 'unknown') ):
        print "I don't recognize the type of this magnetogram."
        return(-1)
                
    print "This is a",TypeMag,"magnetogram on a",str(nLong),"X",str(nLat
              ),"  Phi X sin(Theta) grid."
    if (Long0 > 0):  #The date is not informative for an integral synoptic map
        print "Magnetogram Date="+MapDate

    #Uniform in sinLat and longitude grid
    LatSin_I = np.arcsin(np.linspace(-1. + 1./2/nLat,1. - 1./2/nLat,nLat))
    Long_I = 2.*cPi*np.linspace(0.5/nLong,1. - 0.5/nLong, nLong)
    LongEarth = -1
    if (Long0>0):
        if(TypeOut=='old'):
            FileId = open('CR_Long.txt','w')
            FileId.write(CRNumber + '\n')  
            FileId.write(str(Long0) + '\n')
            FileId.close()
        LongEarth = Long0 + 60
        
    elif(os.path.isfile('Date.in')):
        FileId=open('Date.in','r')
        CRCheck = int(FileId.readline())
        if (CRCheck!=CRNumber):
            print "The Date.in is within CR"+str(CRCheck
                )+"  The synoptic map is for CR"+str(CRNumber)
            return(-1)
        LongEarth = int(FileId.readline())
        FileId.close()
    # Conservative smoothing. Boundary condition: 
    # Periodic in Longitude, reflect in Latitude.
    if (nSmooth>2):
        nSmooth2 = nSmooth//2
        Coef    = 1./(nSmooth*nSmooth)
        BrOrig_II = np.zeros([nLat,nLong+2*nSmooth2])
        for iLat in np.arange(nLat):
            BrOrig_II[iLat,:] = np.hstack((
                Br_II[iLat,nLong-nSmooth2:nLong],
                Br_II[iLat,:],Br_II[iLat,0:nSmooth2]))
        Br_II=np.zeros([nLat,nLong])
        for iLat in np.arange(nLat):
            for iLong in np.arange(nLong):
                for iSubLat in np.arange(nSmooth):
                    iLatExt  = iLat  + iSubLat  - nSmooth2
                    iLatExt  = max([-iLatExt-1,min(
                                [iLatExt, 2*nLat-1-iLatExt])])
                    Br_II[iLat,iLong] += np.sum(
                        BrOrig_II[iLatExt,iLong:iLong+nSmooth])
                Br_II[iLat,iLong]  *= Coef
    if (TypeOut=='old'):
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
        FileId.close()
    elif (TypeOut=='new'):
        FileId = open('fitsfile_idl.out','w')
    
        FileId.write('sin(lat) grid, '+TypeMag+', MapDate = '+MapDate+', Br['+BUnit+']'+'\n')
        FileId.write('       0      0.00000       2       2       1 \n')
        FileId.write('      '+str(nLong)+'     '+str(nLat)+'\n')
        FileId.write(str(Long0)+'     '+str(LongEarth)+'\n') 
        FileId.write('Longitude Latitude Br Long0 LongEarth \n')
    
        for k in np.arange(nLat):
            for l in np.arange(nLong):
                FileId.write(str(Long_I[l]*(180./cPi)) + '   ' + str(
                        LatSin_I[k]*(180./cPi)) + '   ' + str(
                        max([-BMax,min([BMax,Br_II[k,l]])])) + ' \n')
        FileId.close()
   
    return(Br_II,Long0,LongEarth,nLat,nLong,LatSin_I,Long_I)
        

#Remap from iniform in sin(Latitude) magnetogram to that uniform in latitude.
#        - by  Richard A. Frazin July 2014 - February 2015
#        - by  Igor Sokolov, 2016/March:get rid of scipy dependency. 
#          Speed up is by a factor of about 100
def remap(Br_II,Long0,LongEarth,nLat,nLong,LatSin_I,Long_I,BMax):
    #Transpose Br_II Matrix
    BrTransp_II = np.zeros([nLong,nLat])
    for k in np.arange(nLat):
        for l in np.arange(nLong):
            BrTransp_II[l,k] = Br_II[k,l]
    #
    # Centers of the bins for uniform latitude grid:
    #
    LatUniform_I  = np.linspace(-cPi/2 + cPi/2/nLat,cPi/2 - cPi/2/nLat,nLat)
    #
    # Bin boundaries for uniform latitude grid:
    #
    BinBound_I = np.linspace(-cPi/2,cPi/2,nLat+1) #boundaries of latitude grid
    #
    # We will linearly interpolate Br*cos(Latitude) given at sin(theta) grid
    # and integrate it over BinBound_I[l];BinBound_I[l+1] bin. The nodes in
    # which the magnetic field is used for the lth bin have indexes 
    # lMin_I[l] till lMax_I[l]. Construct lMin_I and lMax_I
    #
    lMin_I = np.arange(nLat)
    lMax_I = np.arange(nLat)
    lMin_I[0] = 0
    lMax_I[0] = lMin_I[0]
    while (LatSin_I[lMax_I[0]]<BinBound_I[1]):
        lMax_I[0] = lMax_I[0] + 1
    l=0
    while (l<nLat-1):
        l = l+1
        lMin_I[l]=lMin_I[l-1]
        if (lMin_I[l]<nLat-1):
            while (BinBound_I[l]>LatSin_I[ lMin_I[l]+1]):
                lMin_I[l]=lMin_I[l]+1
                if (lMin_I[l]==nLat-1):
                    break
        lMax_I[l]=lMin_I[l]
        if (lMin_I[l]<nLat-1):
            while (BinBound_I[l+1]>LatSin_I[lMax_I[l]]):
                lMax_I[l]=lMax_I[l]+1
                if (lMax_I[l]==nLat-1):
                    break
    #
    # Now, construct interpolation weights
    #
    CosLat_I      = np.cos(LatSin_I)
    SinBinBound_I = np.sin(BinBound_I)
    Weight_II = np.zeros([nLat,nLat])
    for l in np.arange(nLat):
        if (lMax_I[l]==0 and lMin_I[l]==0): #BB_I[l+1]<LS_I[0]
            Weight_II[l,0]=CosLat_I[0]*(BinBound_I[l+1]-BinBound_I[l])*(
                (BinBound_I[l+1]+BinBound_I[l])/2-BinBound_I[0])/(
                LatSin_I[0]-BinBound_I[0])
        elif (lMax_I[l]==nLat-1 and lMin_I[l]==nLat-1): #BB_I[l]>LS_I[nLat-1]
            Weight_II[l,0]=CosLat_I[nLat-1]*(
                BinBound_I[l+1]-BinBound_I[l])*(
                BinBound_I[nLat] - (BinBound_I[l+1]+BinBound_I[l])/2)/(
                BinBound_I[nLat] - LatSin_I[nLat-1])
        elif (LatSin_I[0]>BinBound_I[l]):#BB_I[l]<LS_U[0]<BB_I[l+1]<LS_I[1]
            Weight_II[l,0]=( (LatSin_I[0]-BinBound_I[l])*(
                    (LatSin_I[0] + BinBound_I[l])/2 - BinBound_I[0])/(
                    LatSin_I[0] - BinBound_I[0]) + (
                    BinBound_I[l+1] - LatSin_I[0])*(
                    LatSin_I[1]-(BinBound_I[l+1] + LatSin_I[0])/2)/(
                    LatSin_I[1]-LatSin_I[0]) )*CosLat_I[0]
            Weight_II[l,1] = (BinBound_I[l+1] - LatSin_I[0])**2/(
                2*(LatSin_I[1]-LatSin_I[0]) )*CosLat_I[1]
        elif (LatSin_I[nLat-1]<BinBound_I[l+1]):
            #LS_I[nLat-2] <BB_I[l]<LS_U[nLat-1]<BB_I[l+1]
            Weight_II[l,0] = (LatSin_I[nLat-1] - BinBound_I[l])**2/(
                2*(LatSin_I[nLat-1]-LatSin_I[nLat-2]) )*CosLat_I[nLat-2]
            Weight_II[l,1]=( (BinBound_I[l+1] - LatSin_I[nLat-1])*(
                    BinBound_I[nLat] - (LatSin_I[nLat-1] + BinBound_I[l+1])/2)/(
                    BinBound_I[nLat] - LatSin_I[nLat-1]) + (
                    LatSin_I[nLat-1] - BinBound_I[l])*(
                   (BinBound_I[l] + LatSin_I[nLat-1])/2 -LatSin_I[nLat-2])/(
                    LatSin_I[nLat-1]-LatSin_I[nLat-2]) )*CosLat_I[nLat-1]
        elif (lMax_I[l]==lMin_I[l] + 1):
             #LS_I[lMin] <BB_I[l]<BB_I[l+1]<LS_I[lMax]
            Weight_II[l,0] = CosLat_I[lMin_I[l]]*(
                BinBound_I[l+1] - BinBound_I[l])*(
                LatSin_I[lMax_I[l]]-(BinBound_I[l+1] + BinBound_I[l])/2)/(
                    LatSin_I[lMax_I[l]] - LatSin_I[lMin_I[l]]) 
            Weight_II[l,1] = CosLat_I[lMax_I[l]]*(
                BinBound_I[l+1] - BinBound_I[l])*(
                (BinBound_I[l+1] + BinBound_I[l])/2 - LatSin_I[lMin_I[l]])/(
                LatSin_I[lMax_I[l]] - LatSin_I[lMin_I[l]])
        else:
            #LS_I[lMin]<BB_I[l]<LS_I[lMin+1]..LS_I[lMax-1]<BB_I[l+1]<LS_I[lMax]
            Weight_II[l,0] = CosLat_I[lMin_I[l]]*(
                LatSin_I[lMin_I[l]+1] - BinBound_I[l])**2/(2*
                (LatSin_I[lMin_I[l]+1] - LatSin_I[lMin_I[l]]))
            Weight_II[l,1] = CosLat_I[lMin_I[l]+1]*(
                LatSin_I[lMin_I[l]+1] - BinBound_I[l])*((
                LatSin_I[lMin_I[l]+1] + BinBound_I[l])/2 - LatSin_I[lMin_I[l]])/(
                LatSin_I[lMin_I[l]+1] - LatSin_I[lMin_I[l]])
            Weight_II[l,lMax_I[l]-lMin_I[l]] = CosLat_I[lMax_I[l]]*(
                BinBound_I[l+1] - LatSin_I[lMax_I[l]-1])**2/(2*
                (LatSin_I[lMax_I[l]] - LatSin_I[lMax_I[l]-1]))
            Weight_II[l,lMax_I[l]-lMin_I[l]-1] = Weight_II[
                l,lMax_I[l]-lMin_I[l]-1] + CosLat_I[lMax_I[l]-1]*(
                BinBound_I[l+1] - LatSin_I[lMax_I[l]-1])*(
                LatSin_I[lMax_I[l]] - ( 
                BinBound_I[l+1] + LatSin_I[lMax_I[l]-1])/2)/(
                LatSin_I[lMax_I[l]] - LatSin_I[lMax_I[l]-1])
            if (lMax_I[l] - lMin_I[l]>2):
                for l1 in np.arange(lMax_I[l] - lMin_I[l]-2):
                    Weight_II[l, 1+l1] = Weight_II[l, 1+l1] + CosLat_I[
                        1+l1+lMin_I[l]]*(LatSin_I[lMin_I[l]+2+l1] - 
                                        LatSin_I[lMin_I[l]+1+l1])/2
                    Weight_II[l, 2+l1] = Weight_II[l, 2+l1] + CosLat_I[
                        2+l1+lMin_I[l]]*(LatSin_I[lMin_I[l]+2+l1] - 
                                        LatSin_I[lMin_I[l]+1+l1])/2
        Weight_II[l,0:lMax_I[l]-lMin_I[l]+1]=Weight_II[
            l,0:lMax_I[l]-lMin_I[l]+1]/(
            SinBinBound_I[l+1] - SinBinBound_I[l])
    
    BrUniform_II  = np.zeros([nLat,nLong])
    for iLong in np.arange(nLong):
        for iLat in np.arange(nLat):
            BrUniform_II[iLat,iLong]=np.sum(
                Weight_II[iLat,0:lMax_I[iLat]-lMin_I[iLat]+1]*BrTransp_II[
                iLong,lMin_I[iLat]:lMax_I[iLat]+1])
  
    FileId = open('uniform_idl.out','w')
    
    FileId.write('Uniform, non-smoothed magnetogram Br[Gauss]'+'\n')
    FileId.write('       0      0.00000       2       2       1 \n')
    FileId.write('      '+str(nLong)+'     '+str(nLat)+'\n')
    FileId.write(str(Long0) +'     '+str(LongEarth)+' \n') 
    FileId.write('Longitude Latitude Br Long0 LongEarth \n')
    
    for k in np.arange(nLat):
         for l in np.arange(nLong):
             FileId.write(str((180./cPi)*Long_I[l]) + '   ' + str(
                     (180./cPi)*LatUniform_I[k]) + '   ' + str(
                 max([-BMax,min([BMax,BrUniform_II[k,l]])])) + ' \n')
    
    FileId.close() 
    return(BrUniform_II,Long0,LongEarth,nLat,nLong,LatSin_I,Long_I)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description="""
       remap_magnetogram.py pre-processes the FITS format magnetograms 
       into ASCII. The script can read the following types of magnetograms:

       GONG Synoptic
       GONG Hourly updated

       The code opens the .fits file and automatically recognizes the type of
       map it is. Within Python, the remapping is done with the remap function contained
       in this file.
    
       The script uses the astropy packages that can be installed, 
       for example, with MacPorts.
       """)
    parser.add_argument('NameFile', help='Input FITS file name including path')
    parser.add_argument('-TypeOut',choices=['old','new','none','remap'],default='new',help=
          """
          Output file: fitsfile.dat+CR_Long.txt (old), BATSRUS standard of .out file (new),
          or remapped from uniform in sin(latitude) grid to uniform in latitude one (remap) 
          """)
    parser.add_argument(
        '-nSmooth',type=int, default=1, help='Smoothing parameter')  
    parser.add_argument(
            '-BMax',type=float, default=1900.,help='Max BrFieldAmplitude')  
    args = parser.parse_args()

 
    cc = readf(args.NameFile,args.TypeOut,args.nSmooth,args.BMax)
    Br_II      = cc[0]
    Long0      = cc[1]
    LongEarth    = cc[2]
    nLat       = cc[3]
    nLong      = cc[4]
    LatSin_I   = cc[5]
    Long_I     = cc[6]
    if(args.TypeOut=='remap'):
            remap(Br_II,Long0,LongEarth,nLat,nLong,LatSin_I,Long_I,args.BMax)
