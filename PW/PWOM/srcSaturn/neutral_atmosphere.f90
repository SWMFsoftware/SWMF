!                                                                      !
!**********************************************************************!
!                                                                      !
!     NEW MODATM SUBROUTINE USING barostatic distribution for saturn   !
!                                                                      !
!**********************************************************************!
!     

! the MAIN program was for debuging the subroutine

!      Program MAIN
!      include 'common_variables.f'
!      REAL array(400,4),j
!
!      do i=1,100
!         j=int(i)*50.+1620.
!         call MODATM(j,a,b,c,d,e)
!         array(i,1)=a
!         array(i,2)=b
!         array(i,3)=c
!         array(i,4)=d
!      enddo
!      do i=1,100
!         write(20,*) (array(i,k),k=1,4)
!      enddo
!      end   

                                                                 
SUBROUTINE MODATM (ALT,XNH2,XNH,XNH2O,XNCH4,TEMP)
  
  use ModCommonVariables
  REAL ALT1,ALT,XNH2,XNH,XNH2O,XNCH4,TEMP,Scaleheight_H2,&
       Scaleheight_H, Scaleheight_H2O,Scaleheight_CH4,logDensity_H2,&
       logDensity_H, logDensity_H2O, logDensity_CH4

! H2O Flux of 10^7
!      DATA logDensity_H2, logDensity_H, logDensity_H2O, logDensity_CH4
!     ;     /9.822053667,7.661700141,3.846886683,8.2/

! H2O Flux of 10^6
  DATA logDensity_H2, logDensity_H, logDensity_H2O, logDensity_CH4&
       /9.822053667,7.661700141,2.846886683,8.2/

!      DATA logDensity_H2, logDensity_H, logDensity_H2O, logDensity_CH4
!     ;     /9.822053667,7.661700141,3.8,6.2/

! casew
!      DATA logDensity_H2, logDensity_H, logDensity_H2O, logDensity_CH4
!     ;     /10.822053667,7.661700141,5.8,8.2/

!      Temp=1000.0
! casew1
!      DATA logDensity_H2, logDensity_H, logDensity_H2O, logDensity_CH4
!     ;     /10.822053667,7.661700141,5.8,8.2/

  Temp=800.0
  
  
  !**********************************************************************
  !     Case A: H2O, High CH4, Temp=800
  !**********************************************************************
  !      DATA Scaleheight_H2,Scaleheight_H, Scaleheight_H2O,Scaleheight_CH4
  !     ; /389.7876287,797.4199227,388.7735132,44.9387829/
  !      DATA logDensity_H2, logDensity_H, logDensity_H2O, logDensity_CH4
  !     ; /9.909664694,7.697395955,3.884960995,4.801655892/
  !      TEMP=800.
  
  !**********************************************************************
  !     Case B: H2O, High CH4, Temp=1000
  !**********************************************************************
  !      DATA Scaleheight_H2,Scaleheight_H, Scaleheight_H2O,Scaleheight_CH4
  !     ; /475.6005082,952.5521232,487.5254497,57.47448837/
  !      DATA logDensity_H2, logDensity_H, logDensity_H2O, logDensity_CH4
  !     ; /9.928961423,7.731980199,3.856513957,5.423258825/
  TEMP=1000.
  
  !**********************************************************************
  !     Case C: H2O, High CH4, Temp=1500
  !**********************************************************************
  !      DATA Scaleheight_H2,Scaleheight_H, Scaleheight_H2O,Scaleheight_CH4
  !     ; /739.97053,1544.640013,719.6795081,85.00641982/
  !      DATA logDensity_H2, logDensity_H, logDensity_H2O, logDensity_CH4
  !     ; /9.822053667,7.661700141,3.846886683,6.549372827/
  !       TEMP=1500.
  
  !**********************************************************************
  !     Case d Temp=600
  !**********************************************************************
  !      Temp=600.
  !
  !**********************************************************************
  !     Case E T=420
  !**********************************************************************
  !      Temp=420.0
  !
  !**********************************************************************
  !     Case W: 100*H2O, High CH4, Temp=1500
  !**********************************************************************
  !
  ! scaleheight_H2O=156 before  
  !      DATA Scaleheight_H2,Scaleheight_H, Scaleheight_H2O,Scaleheight_CH4
  !     ; /739.97053,1544.640013,156.00,85.00641982/
  !      DATA logDensity_H2, logDensity_H, logDensity_H2O, logDensity_CH4
  !     ; /11.1,7.661700141,5.846886683,2.549372827/
  !      TEMP=1500.
  !
  !
  !**********************************************************************
  !     Case W(a&b): H2O, Low CH4, Temp=800
  !**********************************************************************
  !
  !      DATA Scaleheight_H2,Scaleheight_H, Scaleheight_H2O,Scaleheight_CH4
  !     ; /475.6005082,952.5521232,156.0,57.47448837/
  !      DATA logDensity_H2, logDensity_H, logDensity_H2O, logDensity_CH4
  !     ; /9.928961423,7.731980199,4.856513957,3.8/
  !      TEMP=1000.
  
  
  Scaleheight_H2 =1.380658e-26*Temp/(3.3452462e-27*8.52)
  Scaleheight_H  =1.380658e-26*Temp/(1.6726231e-27*8.52)
  Scaleheight_H2O=Scaleheight_H2
  Scaleheight_CH4=1.380658e-26*Temp/(2.65686432e-26*8.6)
  !      Scaleheight_H2O=1.380658e-26*Temp/(2.99e-26*8.52)
  !      Scaleheight_CH4=1.380658e-26*420.0/(2.65686432e-26*8.6)
  
  ! CaseW
  !      Scaleheight_H2O = 144.0
  
  
  ALT1 = ALT/1.E+05
  
  
  XNH2  = (10.**logDensity_H2) *exp(-(ALT1-1400.)/Scaleheight_H2)
  XNH   = (10.**logDensity_H)  *exp(-(ALT1-1400.)/Scaleheight_H)
  XNH2O = (10.**logDensity_H2O)*exp(-(ALT1-1400.)/Scaleheight_H2O)
  XNCH4 = (10.**logDensity_CH4)*exp(-(ALT1-1000.)/Scaleheight_CH4)
  
END SUBROUTINE MODATM
