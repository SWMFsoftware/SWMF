#!/usr/bin/env python

#Remap from iniform in sin(Latitude) magnetogram to that uniform in latitude.
#        - by  Richard A. Frazin July 2014 - February 2015
#        - by  Igor Sokolov, 2016/March:get rid of scipy dependency. 
#          Speed up is by a factor of about 100

import numpy as np
import argparse



def remap(BrTransp_II,nLat,nLong,LatSin_I,Long_I,BMax):
    cPi = np.pi
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
    l = 0
    lMin_I[l] = 0
    lMax_I[l] = lMin_I[l]
    while (LatSin_I[lMax_I[l]]<BinBound_I[l+1]):
        lMax_I[l] = lMax_I[l] + 1

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
    
    Name0 = 'Uniform, non-smoothed magnetogram Br[Gauss]'+'\n' 
    FileId.write(Name0)
    Name0 = '       0      0.00000       2       1       1 \n'
    FileId.write(Name0)
    FileId.write('      '+str(nLong)+'     '+str(nLat)+'\n')
    FileId.write(str(Long0 + 0.5*360./nLong) + ' \n') #longitude shift (important for GONG Hourly)
    FileId.write('Longitude Latitude Br LongitudeShift \n')
    
    for k in np.arange(nLat):
         for l in np.arange(nLong):
             FileId.write(str((180./cPi)*Long_I[l]) + '   ' + str(
                     (180./cPi)*LatUniform_I[k]) + '   ' + str(
                 max([-BMax,min([BMax,BrUniform_II[k,l]])])) + ' \n')
    
    FileId.close() 
    return(BrUniform_II)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description="""
    remap_magnetogram.py pre-processes the FITS format magnetograms 
    into ASCII files that can by read by FDIPS.exe, BATSRUS.exe and SWMF.exe
    and IDL macros. The script can read the following types of magnetograms:

       GONG Synoptic
       GONG Hourly updated

    The code opens the .fits file and automatically recognizes the type of
    map it is, which determines whether it is on a sin(latitude).  

    ./remap_magnetogram.py test.fits test.out
 
    Within Python, the remapping is done with the remap function contained
    in this file.
    
    The script uses the scipy and astropy packages that can be installed, 
    for example, with MacPorts.
    """)
    import read_fits as rf
    parser.add_argument('NameFile', help='Input FITS file name including path')
    parser.add_argument('-BMax',type=float, default=1900., help='Max BrFieldAmplitude')    
    args = parser.parse_args()

 
    cc = rf.readf(args.NameFile,args.BMax)
    BrTransp_II= cc[0]
    Long0      = cc[1]
    nLat       = cc[2]
    nLong      = cc[3]
    LatSin_I   = cc[4]
    Long_I     = cc[5]
 
   
    remap(BrTransp_II,nLat,nLong,LatSin_I,Long_I,args.BMax)
