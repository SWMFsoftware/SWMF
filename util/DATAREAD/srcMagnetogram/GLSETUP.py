#!/usr/bin/env python
####################### MASTER SCRIPT FOR EEGGL or SWMF_GLSETUP#######
####################### CAN BE IMPLEMENTED IN ANY SCRIPT LANGUAGE#####
import subprocess
import os
import read_fits as rf
import numpy as np
import argparse
import GLSETUPAlg as GL
BMax = 1900.0

if __name__ == '__main__':
   parser = argparse.ArgumentParser(
      formatter_class=argparse.RawTextHelpFormatter)
   parser.add_argument('NameFile', help='Input FITS file name including path')
   parser.add_argument('-Out',choices=[
         'old','new','remap','none'],action='append',dest='TypeOut',
                       default=['new'],help=
                       """
          BATSRUS standard of .out file (new),
          or remapped  to uniform in latitude one (remap) 
          """)
   parser.add_argument('-nSmooth',type=int, default=5, help=
                       """
   If nSmooth is ODD integer larger than 1, apply boxcar smoothing on the 
             magnetic field. This can help finding the PIL)
                       """)
   parser.add_argument('-CMESpeed',type=float, default=-1.0, help=
      """
        CME speed in km/s, recovered from observations
      """)
   parser.add_argument('--CMEGrid',action='store_true',
                       help='Output parameters of the refined CME grid')
   parser.add_argument('--UsePIL',action='store_true',
                       help='Use PIL orientation to find GL_Orientation')
   parser.add_argument('--ARSize_OFF',action='store_true',help=
                       """
   ==============This group of paramaters serves to find GLRadius
   There are three ways to set Radius of the GL flux rope, GLRadius:
                       """)
   parser.add_argument('-GLRadius',type=float, default=-1., help=
      """
   1. Explicitly: 
 
      --ARSize_OFF  -GLRadius=GLRadius
      """)
   parser.add_argument('-SizeFactor',type=float, default=26.25, help=
                       """
   2. In terms of Polarity Inversion Line length:

      --ARSize_OFF[ -SizeFactor=SizeFactor]

      In this case the formula, GLRadius=PILLength/SizeFactor is
      applied.  Default SizeFactor is 26.25.
                       """)
   parser.add_argument('--GLRadiusRange',nargs='+',type=float, 
                       default=[0.2,2.0], help=
                       """
   3.  In terms ofthe Active Region size (do not use --ARSize_OFF): 
       [-GLRadiusRange=GLRadiusRange]

       In this case GL is limited, using
       GLRadiusRange = 2-elements array to specify the range for GL
           Radius. Default is [0.2,2.0].
                       """)
   parser.add_argument('-ARMag',type=int, default=1., help=
                       """
   if 1, the AR magnetic parameter is calculated
   based on the average Br around the weighted center; if 2, the
   AR magnetic parameter is calculated based on the average total field
   along the PIL. The corresponding empirical relationships will be
   used accordingly.
   """)
   #PlotRadius = Set up the layer of the magnetogram for 3D input
   #      Cannot be used 2D file, requires /UseBATS. Default is 1.0.
       

   args = parser.parse_args()
   ##################OPTIONAL INPUT PARAMETERS######
   CMESpeed   = args.CMESpeed
   UseNoARSize= args.ARSize_OFF
   GLRadius   = args.GLRadius
   SizeFactor = args.SizeFactor
   GLRadiusRange_I = args.GLRadiusRange
   UsePIL     = args.UsePIL
   ARMag      = args.ARMag
   UseCMEGrid = args.CMEGrid
   TypeOut = args.TypeOut
   ##################END OF PARSER#####################
   #################SERVER SIDE, PYTHON################
   #################PROCESS MAGNETOGRAM###
   ##READ AND SMOOTH, IF DESIRED########################
   cc = rf.readf(args.NameFile,TypeOut,args.nSmooth,BMax)
   nLong        = cc[0]
   nLat         = cc[1]
   nParam       = cc[2]
   Param_I      = cc[3]
   Long0        = Param_I[0]
   LongEarth    = Param_I[1]
   Long_I       = cc[4]  
   Lat_I        = cc[5] 
   Br_C         = cc[6]  
   ###########REMAP, IF DESIRED#######################
   if(any(Type=='remap' for Type in TypeOut)):
      dd= rf.remap(nLong,nLat,nParam,Param_I,Long_I,Lat_I,Br_C,BMax)
      Lat_I      = dd[5]
      Br_C       = dd[6]
   #Info to the idl session is passed via the uniform.out file####
   ############END OF PYTHON FIRST SESSION##########
   ###IDL SESSION IN THE SWMF_GLSETUP/BROWSER SESSION IN EEGGL##
   if CMESpeed<= 0.0: 
      CMESpeed = float(raw_input(
            'Please Input the Observed CME Speed (km/s): '))

   print 'Select the CME Source Region (POSITIVE) with the left button'
   print 'Then select negative region with the right button'
   FileId=open('runidl','w')
   FileId.write(';\n;\n')
   FileId.write("      GLSETUP1,file='uniform.out',/UseBATS,CMESpeed=%5.1f  "%
                (CMESpeed))
   FileId.close()
   ########SHOW MAGNETOGRAM##########################
   ls = subprocess.Popen(["idl", "runidl"],stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT)
   #################PROCESSING STDOUT################
   stdout,stderr=ls.communicate()
   b=stdout[stdout.index('===')+4:len(stdout)]
   a=b.split()
   ###### TAKE TWO COORDINATES FROM TWO CLICKS#######
   xPositive = float(a[0])
   yPositive = float(a[1])
   xNegative = float(a[2])
   yNegative = float(a[3])
   ###The CME speed may be set within the idl/browser session:
   CMESpeed  = float(a[4])
   ## A THING TO DO: IN THE FIRST BLOWSER SESSION OF
   ## THE EEGGL SMESpeed will be provided############
   ##########SHAPE INPUTS FOR THE SECOND SERVER-SIDE SESSION####


   nParam  = 6
   Param_I = np.zeros(nParam)
     
   Param_I[0] = Long0
   Param_I[1] = LongEarth
   Param_I[2] = xPositive
   Param_I[3] = yPositive
   Param_I[4] = xNegative
   Param_I[5] = yNegative
   
   ##SECOND SERVER-SIDE SESSION (PYTHON)#######################
   CC=GL.Alg(nLong,nLat,nParam,Param_I,Long_I,Lat_I,Br_C,
                 CMESpeed,UseNoARSize,GLRadius,SizeFactor, 
                 GLRadiusRange_I, UsePIL, ARMag, UseCMEGrid )
   ##SHAPE INPUT PARAMETERS FOR THE CONCLUDING SESSION#########
   nLong =    CC[0] 
   nLat    =  CC[1] 
   nParam  =  CC[2]
   Param_I =  CC[3]
   Long0     =  Param_I[0]
   LongEarth =  Param_I[1]
   xPositive =  Param_I[2]
   yPositive =  Param_I[3]
   xnegative =  Param_I[5]
   yNegative =  Param_I[6]
   XyARCenter_D = Param_I[6:8]
   nPIL      = (nParam - 8)//2
   xPIL_I   =  Param_I[8:8+nPIL]
   yPIL_I   =  Param_I[8+nPIL:nParam]
   Long_I    = CC[4]
   Lat_I     = CC[5]
   Br_C       =  CC[6]
   PSizeMap_C =  CC[7]
   NSizeMap   =  CC[8]
   FileId=open('runidl','w')
   FileId.write(';\n;\n')
   FileId.write("GLSETUP2, file='uniform.out',/UseBATS")
   FileId.close()
   ###FINAL SESSION: SHOW MAGNETOGRAM AND BIPOLAR STRUCTURE OF AR
   subprocess.call(['idl','run2_GLSETUP'])#'runidl'])
   ###IF THE MASTER SCRIPT IS IN PYTHON, AND A CHILD PROCESS IS IN IDL
   #(1) THE TIME OF IDL SESSION SHOULD BE LIMITED (30 seconds or so) 
   #(2) WINDOWS SHOULD BE CLOSED 
   #(3) FINAL EXIT COMMAND MUST BE PRESENT IN THE IDL SCRIPT######
   #print 'GLSETUP Session is closed. Bye!!!'
  

