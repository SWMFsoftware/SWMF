!                                                                      C
!**********************************************************************C
!                                                                      C
!     NEW MODATM SUBROUTINE USING MSIS MODEL                           C
!                                                                      C
!**********************************************************************C
!                                                                      C
SUBROUTINE MODATM (ALT,XNO2,XNN2,XNO,XNH,XNHE,TEMP)
  use ModCommonVariables
  use ModLatLon,   ONLY: convert_lat_lon
  use EUA_ModMsis90, ONLY: GTD6
  
  
  REAL D(8),T(2),MT(10),ALTL(8)
  
  
  Logical :: IsUseGITM=.false.
  !      Logical :: IsUseGITM=.true.
  
  !
  ! CONVERT INPUT ALT (CM) TO KM VALUES FOR MSIS.
  !
  ALT1 = ALT/1.E05
  !      CALL GGM_PLANET(IART,GLONG,GLAT,SmLon,SmLat)
  CALL convert_lat_lon(Time,SmLat,SmLon,GLAT,GLONG,GLAT2,GLONG2,GmLat,GmLon)
  MASS=48.
  IYR=IYD/1000.
  IDAY=IYD-1000.*IYR
  UT=24.*SEC/86400.
  
  CALL GTD6(IYD,SEC,ALT1,GLAT,GLONG,STL,F107A,F107,AP,MASS,D,T)
  !      write(77,*) IYD,SEC,ALT1,GLAT,GLONG,STL,F107A,F107,AP,MASS
  XNN2 = D(3)
  XNO2 = D(4)
  XNHE = D(1)
  XNO  = D(2)
  XNH  = D(7)
  TEMP = T(2)
  
  If (IsUseGITM .and. ALT1 .le. 500.0) then
     
     call Get_Neutrals(glat,glong,Alt1,Temp,XNO,XNO2,XNN2,XXX,XXX,XXX)
     
  endif
END SUBROUTINE MODATM
