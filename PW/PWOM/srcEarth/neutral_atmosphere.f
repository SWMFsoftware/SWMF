C                                                                      C
C**********************************************************************C
C                                                                      C
C     NEW MODATM SUBROUTINE USING MSIS MODEL                           C
C                                                                      C
C**********************************************************************C
C                                                                      C
      SUBROUTINE MODATM (ALT,XNO2,XNN2,XNO,XNH,XNHE,TEMP)
      use ModCommonVariables
      use ModLatLon,   ONLY: convert_lat_lon
      use EUA_ModMsis90, ONLY: GTD6


      REAL D(8),T(2),MT(10),ALTL(8)

 
      LOGICAL REPT
      Logical :: IsUseGITM=.false.
!      Logical :: IsUseGITM=.true.
      DATA REPT/.FALSE./
C
C CONVERT INPUT ALT (CM) TO KM VALUES FOR MSIS.
C
      ALT1 = ALT/1.E+05
!      CALL GGM_PLANET(IART,GLONG,GLAT,GMLONG,GMLAT)
      CALL convert_lat_lon(Time,GMLAT,GMLONG,GLAT,GLONG)
      MASS=48.
      IYR=IYD/1000.
      IDAY=IYD-1000.*IYR
      UT=24.*SEC/86400.
      IF (REPT) GO TO 20
C      WRITE (iUnitOutput,1) IYR,IDAY,UT,STL,GMLAT,GMLONG,GLAT,GLONG,F107A,F107,
C     #AP
C1     FORMAT(/' YEAR,DAY,UT,LOCTIM,GMLAT,GMLONG,GLAT,GLONG'/
C     #I3,I4,2F5.1,4F6.1//' F10.7AV,F10.7,AP(DAILY),AP(UT),AP(UT-3)'
C     #',AP(UT-6),AP(UT-9),AP(AVER(UT-12,UT-33)),AP(AVER(UT-36,UT-59))'/
C     #2F5.0,7F5.1/)
      REPT=.TRUE.
c20    CALL GTS5(ALT1,MASS,D,T)


20    CALL GTD6(IYD,SEC,ALT1,GLAT,GLONG,STL,F107A,F107,AP,MASS,D,T)
c      write(77,*) IYD,SEC,ALT1,GLAT,GLONG,STL,F107A,F107,AP,MASS
        XNN2 = D(3)
        XNO2 = D(4)
        XNHE = D(1)
        XNO  = D(2)
        XNH  = D(7)
        TEMP = T(2)

        If (IsUseGITM .and. ALT1 .le. 500.0) then
           
          ! call Get_Neutrals(glat,glong,Alt1,XXX,XXX,XXX,XNN2,XXX,XXX,XXX)
          ! write(*,*) (XNN2-D(3))/D(3)
          call Get_Neutrals(glat,glong,Alt1,Temp,XNO,XNO2,XNN2,XXX,XXX,XXX)
           
        endif
      END
