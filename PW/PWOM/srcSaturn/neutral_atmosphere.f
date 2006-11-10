C                                                                      C
C**********************************************************************C
C                                                                      C
C     NEW MODATM SUBROUTINE USING barostatic distribution for saturn   C
C                                                                      C
C**********************************************************************C
C     

C the MAIN program was for debuging the subroutine

C      Program MAIN
C      include 'common_variables.f'
C      REAL array(400,4),j
C
C      do i=1,100
C         j=int(i)*50.+1620.
C         call MODATM(j,a,b,c,d,e)
C         array(i,1)=a
C         array(i,2)=b
C         array(i,3)=c
C         array(i,4)=d
C      enddo
C      do i=1,100
C         write(20,*) (array(i,k),k=1,4)
C      enddo
C      end   

                                                                 
      SUBROUTINE MODATM (ALT,XNH2,XNH,XNH2O,XNCH4,TEMP)

      use ModCommonVariables
      REAL ALT1,ALT,XNH2,XNH,XNH2O,XNCH4,TEMP,Scaleheight_H2,
     ; Scaleheight_H, Scaleheight_H2O,Scaleheight_CH4,logDensity_H2,
     ; logDensity_H, logDensity_H2O, logDensity_CH4

      DATA logDensity_H2, logDensity_H, logDensity_H2O, logDensity_CH4
     ;     /9.822053667,7.661700141,3.846886683,8.2/

c      DATA logDensity_H2, logDensity_H, logDensity_H2O, logDensity_CH4
c     ;     /9.822053667,7.661700141,4.8,8.2/

c casew
c      DATA logDensity_H2, logDensity_H, logDensity_H2O, logDensity_CH4
c     ;     /10.822053667,7.661700141,5.8,8.2/

c      Temp=1000.0
c casew1
c      DATA logDensity_H2, logDensity_H, logDensity_H2O, logDensity_CH4
c     ;     /10.822053667,7.661700141,5.8,8.2/

c      Temp=800.0


C**********************************************************************
C     Case A: H2O, High CH4, Temp=800
C**********************************************************************
c      DATA Scaleheight_H2,Scaleheight_H, Scaleheight_H2O,Scaleheight_CH4
c     ; /389.7876287,797.4199227,388.7735132,44.9387829/
c      DATA logDensity_H2, logDensity_H, logDensity_H2O, logDensity_CH4
c     ; /9.909664694,7.697395955,3.884960995,4.801655892/
c      TEMP=800.

C**********************************************************************
C     Case B: H2O, High CH4, Temp=1000
C**********************************************************************
c      DATA Scaleheight_H2,Scaleheight_H, Scaleheight_H2O,Scaleheight_CH4
c     ; /475.6005082,952.5521232,487.5254497,57.47448837/
c      DATA logDensity_H2, logDensity_H, logDensity_H2O, logDensity_CH4
c     ; /9.928961423,7.731980199,3.856513957,5.423258825/
      TEMP=1000.

C**********************************************************************
C     Case C: H2O, High CH4, Temp=1500
C**********************************************************************
c      DATA Scaleheight_H2,Scaleheight_H, Scaleheight_H2O,Scaleheight_CH4
c     ; /739.97053,1544.640013,719.6795081,85.00641982/
c      DATA logDensity_H2, logDensity_H, logDensity_H2O, logDensity_CH4
c     ; /9.822053667,7.661700141,3.846886683,6.549372827/
c      TEMP=1500.

C**********************************************************************
C     Case d Temp=600
C**********************************************************************
C      Temp=600.

C**********************************************************************
C     Case E T=420
C**********************************************************************
c      Temp=420.0

C**********************************************************************
C     Case W: 100*H2O, High CH4, Temp=1500
C**********************************************************************

c scaleheight_H2O=156 before  
c      DATA Scaleheight_H2,Scaleheight_H, Scaleheight_H2O,Scaleheight_CH4
c     ; /739.97053,1544.640013,156.00,85.00641982/
c      DATA logDensity_H2, logDensity_H, logDensity_H2O, logDensity_CH4
c     ; /11.1,7.661700141,5.846886683,2.549372827/
c      TEMP=1500.


C**********************************************************************
C     Case W(a&b): H2O, Low CH4, Temp=800
C**********************************************************************

c      DATA Scaleheight_H2,Scaleheight_H, Scaleheight_H2O,Scaleheight_CH4
c     ; /475.6005082,952.5521232,156.0,57.47448837/
c      DATA logDensity_H2, logDensity_H, logDensity_H2O, logDensity_CH4
c     ; /9.928961423,7.731980199,4.856513957,3.8/
c      TEMP=1000.


      Scaleheight_H2 =1.380658e-26*Temp/(3.3452462e-27*8.52)
      Scaleheight_H  =1.380658e-26*Temp/(1.6726231e-27*8.52)
      Scaleheight_H2O=Scaleheight_H2
      Scaleheight_CH4=1.380658e-26*Temp/(2.65686432e-26*8.6)
c      Scaleheight_H2O=1.380658e-26*Temp/(2.99e-26*8.52)
c      Scaleheight_CH4=1.380658e-26*420.0/(2.65686432e-26*8.6)

C CaseW
c      Scaleheight_H2O = 144.0


      ALT1 = ALT/1.E+05
      
      
      XNH2  = (10.**logDensity_H2) *exp(-(ALT1-1400.)/Scaleheight_H2)
      XNH   = (10.**logDensity_H)  *exp(-(ALT1-1400.)/Scaleheight_H)
      XNH2O = (10.**logDensity_H2O)*exp(-(ALT1-1400.)/Scaleheight_H2O)
      XNCH4 = (10.**logDensity_CH4)*exp(-(ALT1-1000.)/Scaleheight_CH4)

      END
