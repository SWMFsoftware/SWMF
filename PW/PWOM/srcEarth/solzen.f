C
C
C
C
      SUBROUTINE SOLZEN (IDATE, UTG, GLAT, GLONG, SZA)
C
C Returns Solar Zenith Angle SZA in degrees for specified date in form yyddd,
C universal time in seconds, geographic latitude and longitude in degrees.
C
C S.C. Solomon, 9/88
C
      DATA PI/3.1415926536/
C
      RLAT = GLAT * PI/180.
      RLONG = GLONG * PI/180.
      CALL SUNCOR (IDATE, UTG, SDEC, SRASN, GST)
      RH = SRASN - (GST+RLONG)
      COSSZA = SIN(SDEC)*SIN(RLAT) + COS(SDEC)*COS(RLAT)*COS(RH)
CALEX sza=solar zenith angle      
      SZA = ACOS(COSSZA) * 180./PI
      RETURN
      END
