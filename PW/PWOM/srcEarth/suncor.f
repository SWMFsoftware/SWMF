C
C
C
C
C Subroutine SUNCOR returns the declination SDEC and right ascension SRASN
C of the sun in GEI coordinates, radians, for a given date IDATE in yyddd
C format and universal time UTG in seconds.  Greenwich Sidereal Time GST
C in radians is also returned.  Reference:  C.T. Russell, Geophysical
C Coordinate Transforms
C
      SUBROUTINE SUNCOR (IDATE, UTG, SDEC, SRASN, GST)
      DATA PI/3.1415926536/
C
      FDAY=UTG/86400.
      IYR=IDATE/1000
      IDAY=IDATE-IYR*1000
      DJ=365*IYR+(IYR-1)/4+IDAY+FDAY-0.5
      T=DJ/36525.
      VL=AMOD(279.696678+.9856473354*DJ,360.)
      GST=AMOD(279.696678+.9856473354*DJ+360.*FDAY+180.,360.) * PI/180.
      G=AMOD(358.475845+.985600267*DJ,360.) * PI/180.
      SLONG=VL+(1.91946-.004789*T)*SIN(G)+.020094*SIN(2.*G)
      OBLIQ=(23.45229-0.0130125*T) *PI/180.
      SLP=(SLONG-.005686) * PI/180.
      SIND=SIN(OBLIQ)*SIN(SLP)
      COSD=SQRT(1.-SIND**2)
      SDEC=ATAN(SIND/COSD)
      SRASN=3.14159-ATAN2(1./TAN(OBLIQ)*SIND/COSD,-COS(SLP)/COSD)
      RETURN
      END
