#!/usr/bin/env python
import numpy as np
BMax = 1900.0

cPi = np.pi
Rad2Deg = 180/cPi

def round_my(x):
   i = int(round(x))
   return(i)

def Alg(nLong,nLat,nParam, Param_I,Long_I,Lat_I,Br_C,
        CMESpeed,UseNoARSize,GLRadius,SizeFactor, 
        GLRadiusRange_I, UsePIL, ARMag, UseCMEGrid):
   Long0     = Param_I[0]
   LongEarth = Param_I[1] 
   xPositive = Param_I[2]
   yPositive = Param_I[3]
   xNegative = Param_I[4]
   yNegative = Param_I[5]
   
   print "You chose  x/y for Positive and negative spots:" 
   print "{0:5.1f} {1:5.1f} {2:5.1f} {3:5.1f}".format(
      xPositive, yPositive,xNegative, yNegative)
   iBoxSize = round_my((16.*nLong)/360)
   
   #Calculate the weighted centers forpositive and negative spots.
   xPositiveWeight=0.
   yPositiveWeight=0.
   TotalPositiveFlux=0.
   for j in np.arange(yPositive-iBoxSize//2,yPositive+iBoxSize//2+1):
      for i in np.arange(xPositive-iBoxSize//2,xPositive+iBoxSize//2+1):
         if (Br_C[j,i] > 0.0):
            xPositiveWeight=xPositiveWeight+Br_C[j,i]*i
            yPositiveWeight=yPositiveWeight+Br_C[j,i]*j
            TotalPositiveFlux=TotalPositiveFlux+Br_C[j,i]
   xPositive=xPositiveWeight/TotalPositiveFlux
   yPositive=yPositiveWeight/TotalPositiveFlux

   xNegativeWeight=0.
   yNegativeWeight=0.
   TotalNegativeFlux=0.
   for j in np.arange(yNegative-iBoxSize//2,yNegative+iBoxSize//2+1):
      for i in np.arange(xNegative-iBoxSize//2,xNegative+iBoxSize//2+1):
         if (Br_C[j,i] < 0.0):
            xNegativeWeight=xNegativeWeight+Br_C[j,i]*i
            yNegativeWeight=yNegativeWeight+Br_C[j,i]*j
            TotalNegativeFlux=TotalNegativeFlux+Br_C[j,i]
            
   xNegative=xNegativeWeight/TotalNegativeFlux
   yNegative=yNegativeWeight/TotalNegativeFlux
   print "Weighted  x/y for Positive and negative spots:" 
   print "{0:4.1f} {1:4.1f} {2:4.1f} {3:4.1f}".format(
      xPositive, yPositive,xNegative, yNegative)
  
   #Find center of the active region as the point on the line
   #connecting the positive and negetive center in which the MF is minumal
   #(i.e as intersection of this with PIL, 
   #herewith PIL=Polarity Inversion Line 
   nProfile = max([round_my(abs(xPositive-xNegative)), 
                   round_my(abs(yPositive-yNegative))]) +1
   Misc = BMax +1.0 
   for i in np.arange(nProfile):
      xProfile = xPositive+(
         xNegative - xPositive)*i/(nProfile-1.)
      yProfile = yPositive+(
         yNegative - yPositive)*i/(nProfile-1.)
      AbsBr = abs(Br_C[round_my(yProfile),round_my(xProfile)])
      if (AbsBr < Misc):
         Misc = AbsBr
         XyARCenter_D = [xProfile,yProfile]
   print "Center for Active region:" 
   print "{0:4.1f} {1:4.1f}".format(XyARCenter_D[0],
      XyARCenter_D[1])
   iXARCenter=round_my(XyARCenter_D[0])
   iYARCenter=round_my(XyARCenter_D[1])
   #Save the physical coordinates to use as the GL flux rope center 
   GL_Latitude  = Lat_I[ iYARCenter]*Rad2Deg
   GL_Longitude = Long_I[iXARCenter]*Rad2Deg
   if Long0>0:
      GL_Longitude +=Long0
      if GL_Longitude>=360:
         lLongitude-=360
   print "GL_Longitude: {0:4.1f} GL_Latitude:{1:4.1f}".format(
      GL_Longitude, GL_Latitude)


   #Distances from the AR center and spot centers:
   DisCenter_C = np.zeros([nLat,nLong])
   PCenter_C   = np.zeros([nLat,nLong])
   NCenter_C   = np.zeros([nLat,nLong])
   for j in np.arange(nLat):
      for i in np.arange(nLong):
         DisCenter_C[j,i] = np.sqrt(
            (i-XyARCenter_D[0])**2+(j-XyARCenter_D[1])**2)
         PCenter_C[j,i]  = np.sqrt((i-xPositive)**2+(j-yPositive)**2)
         NCenter_C[j,i]  = np.sqrt((i-xNegative)**2+(j-yNegative)**2)
   #Calculate Active Region Size for determining the GL size.
   #Below, the distances are mostly compared with the distance between
   #the centers of positive and negative spots.
   NPCentersDist = np.sqrt((xPositive - xNegative)**2 + 
                           (yPositive - yNegative)**2)
    
   ARThreshold=12.*nLong/360*NPCentersDist/13.6
   PSizeMap_C = np.where(Br_C >  20., Br_C, 0.)
   PSizeMap_C = np.where(PCenter_C <= ARThreshold, PSizeMap_C, 0.)   
   NSizeMap_C = np.where(Br_C < -20., Br_C, 0.)
   NSizeMap_C = np.where(NCenter_C <= ARThreshold, NSizeMap_C, 0.)
  
 
   #Calculate the gradient of the Br field
   DDx_C = np.zeros([nLat,nLong])
   DDy_C = np.zeros([nLat,nLong]) 
   DDx_C[:,1:nLong-1] = (Br_C[:,2:nLong] - Br_C[:,0:nLong-2])/2.
   DDx_C[:,0        ] =  Br_C[:,1      ] - Br_C[:,0        ]
   DDx_C[:,nLong-1  ] =  Br_C[:,nLong-1] - Br_C[:,nLong-2  ]

   DDy_C[1:nLat-1, :] = (Br_C[2:nLat, :] - Br_C[0:nLat-2, :])/2.
   DDy_C[0       , :] =  Br_C[1     , :] - Br_C[0       , : ]
   DDy_C[nLat-1  , :] =  Br_C[nLat-1, :] - Br_C[nLat-2  , :]
   GradBr_C=np.sqrt(DDx_C**2 + DDy_C**2)
   #Calculate a Bit Map (1/0) for GradB multiplied by MF
   GradBrMap_C = np.where(GradBr_C > 0.5, Br_C, 0.)

   # Cell size is used to divide the magnetogram to sub regions 
   #in order to determine the PIL. 
   nCell = 1

   #Setup the threshold for selecting cells near the PIL.
   BThreshold = 2.0

   #Calculate the PILBitMap_C (1/0) for determining the PIL.
   M = nLong//nCell
   N = nLat//nCell
   PILBitMap_C = np.zeros([nLat,nLong])
   for j in np.arange(N-1):
      for i in np.arange(M-1): 
         if(np.amin(Br_C[j*nCell:(j+1)*nCell+1,i*nCell:(i+1)*nCell+1]) 
            < -BThreshold  and  
            np.amax(Br_C[j*nCell:(j+1)*nCell+1,i*nCell:(i+1)*nCell+1])
            >  BThreshold):
            #PIL intersects this quadrant. Among these 4 points 
            #we do not include those in which the gradient is too small: 
            PILBitMap_C[j*nCell:(j+1)*nCell+1,i*nCell:(i+1)*nCell+1
                        ] = GradBrMap_C[
               j*nCell:(j+1)*nCell+1,i*nCell:(i+1)*nCell+1]
             

    #Distance cut-off for determining the PIL.
    # 8.*(nLon/360.)*DisWeight/13.6        for calculating PIL_Length
    # 3 for calculating PIL orientation)
    # min([6.*(nLon/360.),8.*(nLon/360.)*DisWeight/13.6]) -                 
    # for calculating bt_pil 6.*(nLon/360.) and for vizualization
   DisThreshold=(8.0*nLong)/360*NPCentersDist/13.6
   DisMax      =(6.0*nLong)/360
   DisThresholdS=3
      
   if (UsePIL):
      #Calculate the orientation of the flux rope according to PIL 
      #(make it perpendicular to PIL).
       PILMap_C = np.where(DisCenter_C<=DisThresholdS, PILBitMap_C, 0.)
       
       iYPILPoints_I,iXPILPoints_I =np.where(PILMap_C > 0)
       MM = len(iYPILPoints_I)
   
       #    PIL_xx=PIL_x[sort(PIL_x)]
       #    PIL_yy=PIL_y[sort(PIL_x)]
       #    PIL_fit=ladfit(PIL_xx,PIL_yy,/double)  
              
       if PIL_fit[1] ==0: 
          if Br_C[XyARCenter_D[1]-2,XyARCenter_D[1]] < 0:

             r1=[0,-1] 
          else:
             r1=[0,1]
       else:
          aa_PIL=-1./PIL_fit[1]
          if abs(aa_PIL)<= 1.: 
             bb_PIL=ar_center[1]-aa_PIL*ar_center[0]
             xx=[ar_center[0]-2,ar_center[0]-3,ar_center[0]-4]
             yy=aa_PIL*xx+bb_PIL
             ave_field=(br_field(xx[0],yy[0])+br_field(xx[1],yy[1])+br_field(xx[2],yy[2]))/3.
             if ave_field < 0:
                r1=[-1.,-aa_PIL] 
             else:
                r1=[1.,aa_PIL]
          else:
             bb_PIL=ar_center[1]-aa_PIL*ar_center[0]
             yy=[ar_center[1]-2,AR_center[1]-3,AR_Center[1]-4]
             xx=floor((yy-bb_PIL)/aa_PIL)
             ave_field=(br_field(xx[0],yy[0])+br_field(xx[1],yy[1])+br_field(xx[2],yy[2]))/3.
             if ave_field < 0:
                r1=[-signum(aa_PIL),-signum(aa_PIL)*aa_PIL] 
             else:
                r1=[signum(aa_PIL),signum(aa_PIL)*aa_PIL]    
       if array_equal(PIL_xx,PIL_xx[0]):
          if br_field(ar_center[0]-2,ar_center[1]) < 0: 
             r1=[-1,0] 
          else:
             r1=[1,0]
          
       
   else:
      #Calculate the GL flux rope orientation from the two weighted points.
      r1=[xNegative-xPositive,yNegative-yPositive]
   r1=r1/np.sqrt(r1[0]**2+r1[1]**2) 
   r2=[1.0,0.0]
   GL_Orientation=np.arccos(r1[0]*r2[0]+r1[1]*r2[1])*Rad2Deg
   if r1[1] < 0:
      GL_Orientation=360-GL_Orientation

   if UseNoARSize:
      if (GLRadius<=0.):
       #At this moment, the PIL length is represented by degree and does not 
       #take into account the effect of different latitude. It will be 
       #improved later.
         PILMap_C = np.where(DisCenter_C<=DisThreshold, PILBitMap_C, 0.)
         PIL_Length=np.count_nonzero(PILMap_C)/2.*360./nLong
         GLRadius=PIL_Length/SizeFactor
   else:
      
       #Use Active Region size to specify the GL flux rope Radius.
      ARSize=np.count_nonzero(PSizeMap_C) +np.count_nonzero(NSizeMap_C)      
      GLRadius=0.8/280.*ARSize
      GLRadius = max([min([GLRadius,GLRadiusRange_I[1]]),GLRadiusRange_I[0]])

   PILMap_C = np.where(DisCenter_C<=min([DisMax,DisThreshold]), 
                       PILBitMap_C, 0.)
   #Calculate the poloidal flux needed for the observed CME velocity.
   #These relationships are based on the GONG magnetogram with nsmooth = 5
   if ARMag == 1:
      RegionSize_ARMag=round_my((4.0*nLong)/360)
      br_ar=np.mean(
         abs(Br_C[XyARCenter_D[1]-RegionSize_ARMag//2:
                  XyARCenter_D[1]+RegionSize_ARMag//2+1,
                  XyARCenter_D[0]-RegionSize_ARMag//2:
                  XyARCenter_D[0]+RegionSize_ARMag//2+1]))
      GL_poloidal=(CMESpeed*br_ar**0.43989278-3043.9307)/565.05018
   else:
      #Calculate the AR magnetic field strength in order to determine the
      #right poloidal flux needed.
      iYPILPoints_I, iXPILPoints_I =np.where(PILMap_C != 0.)
      NN = len(iYPILPoints_I)
      bt_pil=0.
      for i in np.arange(NN):
         bt_pil+=bt_field[iYPILPoints_I[i],iXPILPoints_I[i]]

      bt_pil=bt_pil/NN
      GL_poloidal=(CMESpeed*bt_pil**0.58148678-2814.1030)/507.60065

  
   #Print WARNING information is GL_Bstrength is negative                                               
   if GL_poloidal <= 0 :
      print '*********************************************'
      print 'WARNING: CALCULATION FAILED!USE WITH CAUTION!'
      print 'Either the active region is too weak or the'    
      print 'CME speed is too small!'
      print 'GL Poloidal Flux is set to 0!'
      print '*********************************************'
      GL_poloidal = 0.0
       #Relationship between the PIL length and the GL flux rope Radius.   
       #This factor is now based on the 2011 March 7 CME. More tests  
       #are needed in order to get a more precise value.  
   
   
      #Relationship between the GL Poloidal flux and GL Bstrength.
      #Flux rope helicity is determined by the hemisphere, northern
      #hemisphere - negative helicity.
   
   GL_Bstrength=GL_poloidal/(21.457435*GLRadius**4)
   if GL_Latitude > 0: 
      GL_Bstrength = -GL_Bstrength
      

     #Recommended GL flux rope parameters
   Distance = 1.8
   Stretch  = 0.6
   
   print '========================================'
   print 'The Recommended GL FLux Rope Parameters'
   print '========================================'
   print '#CME'
   print '%6.2f                Latitude: '%(GL_Latitude)
   print '               Longitude: %6.2f'%(GL_Longitude)
   print '             Orientation: %6.2f'%(GL_Orientation)
   print '                  Radius: %6.2f'%(GLRadius)
   print '               Bstrength: %6.2f'%(GL_Bstrength)
   print '         Stretch (FIXED): %6.2f'%(Stretch)
   print '        Distance (FIXED): %6.2f'%(Distance)
   print '             Height [Rs]: %6.2f'%(
      GLRadius + Distance - Stretch - 1.0)
   print '      Angular size [deg]: %6.2f'%(
      2*GLRadius/Distance*Rad2Deg)
   print ' Poloidal flux [1E21 Mx]: %6.2f'%(GL_poloidal)
   print '-----------------------------------------'

   if UseCMEGrid:
      #Calculate the CME grid refinement parameters based on the flux rope
      #location and size.                                                

      CMEbox_Start=[1.1,GL_Longitude-40.*GLRadius,GL_Latitude-20.*GLRadius]
      CMEbox_End=[20.0,GL_Longitude+40.*GLRadius,GL_Latitude+20.*GLRadius]
      
      print '=========================================='
      print 'The Recommended Grid Refinement Parameters'
      print '=========================================='
      print '              R_Start: %6.2f'% (CMEbox_start[0])
      print '                R_End: %6.2f'% (CMEbox_end[0])
      print '      Longitude_Start: %6.2f'% ( CMEbox_start[1])
      print '        Longitude_End: %6.2f'% ( CMEbox_end[1])
      print '       Latitude_Start: %6.2f'% ( CMEbox_start[2])
      print '         Latitude_End: %6.2f'% ( CMEbox_end[2])
      print '-----------------------------------------'
   
   if(Br_C[iYARCenter,iXARCenter]>0):
      iYPIL_I,iXPIL_I=np.where(PILMap_C>0)
   else:
      iYPIL_I,iXPIL_I=np.where(PILMap_C<0)
   nPIL=len(iXPIL_I)
  
   nParam = 8 + 2*nPIL
   Param_I = np.zeros(nParam)
   Param_I[0:8] = [Long0,LongEarth,xPositive,yPositive,xNegative,yNegative,
         XyARCenter_D[0],XyARCenter_D[1]]
   Param_I[8:8+nPIL]=iXPIL_I
   Param_I[8+nPIL:8+2*nPIL]=iYPIL_I
   FileId = open('uniform.out','w')
    
   FileId.write('Uniform, non-smoothed magnetogram Br[Gauss]'+'\n')
   FileId.write(
      '       0      0.00000       2      %2d       3 \n'% nParam)
   FileId.write('      '+str(nLong)+'     '+str(nLat)+'\n')
   FileId.write(
      ' {0:5.1f} {1:5.1f} {2:5.1f} {3:5.1f} {4:5.1f} {5:5.1f} {6:5.1f} {7:5.1f}'.format(
         Long0,LongEarth,xPositive,yPositive,xNegative,yNegative,
         XyARCenter_D[0],XyARCenter_D[1]) )
   for j in iXPIL_I:
      FileId.write(' %4.0f'%float(j))
   for j in iYPIL_I:
      FileId.write(' %4.0f'%float(j))
   FileId.write('\n')
   FileId.write(
      'Longitude Latitude Br PMap NMap Long0 LongEarth xP yP xN yN xC yC xPIL1:{0:2d} yPIL1:{0:2d}\n'.format(nPIL))
    
   for k in np.arange(nLat):
      for l in np.arange(nLong):
         FileId.write("{0:6.1f} {1:6.1f} {2:14.6e} {3:14.6e} {4:14.6e}\n".format(
               (180./cPi)*Long_I[l],(180./cPi)*Lat_I[k],
               max([-BMax,min([BMax,Br_C[k,l]])]),
                      PSizeMap_C[k,l],NSizeMap_C[k,l]))
    
   FileId.close()


   return(nLong,nLat,nParam, Param_I,Long_I,Lat_I,Br_C,PSizeMap_C,NSizeMap_C)
