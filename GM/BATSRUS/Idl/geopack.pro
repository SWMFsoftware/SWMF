;  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
;  For more information, see http://csem.engin.umich.edu/tools/swmf
;==============================================
;;; Function jday and procedures gse_gsm, gsm_gse, gsm_gse_direction 
;;; were added by G. Toth and A. Ridley
;==============================================
function jday, year, mon, day
  if mon lt 2 then return, day
  dayofmon = [31,28,31,30,31,30,31,31,30,31,30,31]
  if year mod 4 eq 0 and (year mod 100 ne 0 or year mod 400 eq 0) $
    then dayofmon(1) = dayofmon(1) + 1
  return, total(dayofmon(0:mon-2)) + day
end
;==============================================
pro gse_gsm, field, time
; convert field array from GSE to GSM
gsm_gse_direction, field, time, -1
end
;==============================================
pro gsm_gse, field, time
; convert field array from GSM to GSE
gsm_gse_direction, field, time, +1
end
;==============================================
pro gsm_gse_direction,field,time,direction

if direction ne 1 and direction ne -1 then begin
   print,'direction should be 1 or -1'
   help,direction
   retall
endif
  
siz = size(field)
if siz(0) ne 2 then begin
    print,'field should be 2D array'
    help,field
    retall
endif

if siz(1) ne 3 then begin
    print,'field should be of size 3 in the first dimension'
    help,field
    retall
endif

n = siz(2)

siz = size(time)
if siz(0) ne 1 then begin
    print,'time should be a 1D array'
    help,time
    retall
endif

if siz(1) ne n then begin
    print,'field and time arrays should have the same size'
    help,field,time
    retall
endif

if direction eq 1 then begin
    for i = 0, n-1 do begin
        GSMGSE,field(0,i),field(1,i),field(2,i),x,y,z,+1,epoch=time(i)
        field(0,i) = x
        field(1,i) = y
        field(2,i) = z
    endfor
endif else begin
    for i = 0, n-1 do begin
        GSMGSE,x,y,z,field(0,i),field(1,i),field(2,i),-1,epoch=time(i)
        field(0,i) = x
        field(1,i) = y
        field(2,i) = z
    endfor
endelse

end
;====================================================================

; This file was updated on July 3rd, 2001 to accurately produce coordinates
; in GSM and GSE for data dates January, 2001 to 2005.  Prior to this release
; the coordinates starting on January 1, 2001 were INCORRECT.  If you were
; using a previous version of the software, please replace it immediately.
; Sorry for any inconvenience this oversight might have caused.
; T. Kovalick, Raytheon ITSS
; Further update up to 2015 by G. Toth and X. Meng.
;c--------------------------------------------------------------------
;c
      pro SUN,IYR,IDAY,IHOUR,MIN,ISEC,GST,SLONG,SRASN,SDEC,epoch=epoch
      if n_elements(epoch) ne 0 then begin
         cdf_epoch,epoch,iyr,month,idom,ihour,min,isec,/break
         iday = jday(iyr,month,idom)
;         ical,iyr,iday,month,idom,/idoy
      endif
;TJK added the next three lines for y2k - based on Scott Boardsen's changes.
      iyr = long(iyr)
      if iyr lt 50 then iyr = 2000 + iyr
      if iyr ge 50 and iyr lt 1900 then iyr = 1900 + iyr
;C
;C  CALCULATES FOUR QUANTITIES NECESSARY FOR COORDINATE TRANSFORMATIONS
;C  WHICH DEPEND ON SUN POSITION (AND, HENCE, ON UNIVERSAL TIME AND SEASON)
;C
;C-------  INPUT PARAMETERS:
;C  IYR,IDAY,IHOUR,MIN,ISEC -  YEAR, DAY, AND UNIVERSAL TIME IN HOURS, MINUTES,
;C    AND SECONDS  (IDAY=1 CORRESPONDS TO JANUARY 1).
;C
;C-------  OUTPUT PARAMETERS:
;C  GST - GREENWICH MEAN SIDEREAL TIME, SLONG - LONGITUDE ALONG ECLIPTIC
;C  SRASN - RIGHT ASCENSION,  SDEC - DECLINATION  OF THE SUN (RADIANS)
;C  THIS SUBROUTINE HAS BEEN COMPILED FROM: RUSSELL C.T., COSM.ELECTRO-
;C  DYN., 1971, V.2,PP.184-196.
;C
;C
;C                   AUTHOR: Gilbert D. Mead
;C
;C
;	IMPLICIT NONE

;	REAL GST,SLONG,SRASN,SDEC,RAD,T,VL,G,OBLIQ,SOB,SLP,SIND,
;    1       COSD,SC

;	INTEGER IYR,IDAY,IHOUR,MIN,ISEC

;      DOUBLE PRECISION DJ,FDAY
      RAD=57.295779513d0
      IF((IYR LT 1901) OR (IYR GT 2099))THEN RETURN
      FDAY=(IHOUR*3600D0+MIN*60D0+ISEC)/86400.D0
      DJ=365*(IYR-1900)+FIX((IYR-1901)/4)+IDAY-0.5D0+FDAY
      T=DJ/36525D0
      VL=(279.696678D0+0.9856473354D0*DJ) MOD 360.D0
      GST=( (279.690983D0+.9856473354D0*DJ+360.*FDAY+180.) MOD 360.D0)/RAD
      G=( (358.475845D0+0.985600267D0*DJ) MOD 360.D0)/RAD
      SLONG=(VL+(1.91946D0-0.004789D0*T)*SIN(G)+0.020094D0*SIN(2.*G))/RAD
      IF(SLONG GT 6.2831853D0)THEN SLONG=SLONG-6.2831853D0
      IF (SLONG LT 0D0)THEN SLONG=SLONG+6.2831853D0
      OBLIQ=(23.45229D0-0.0130125D0*T)/RAD
      SOB=SIN(OBLIQ)
      SLP=SLONG-9.924D-5
;C
;C   THE LAST CONSTANT IS A CORRECTION FOR THE ANGULAR ABERRATION  DUE TO
;C   THE ORBITAL MOTION OF THE EARTH
;C
      SIND=SOB*SIN(SLP)
      COSD=SQRT(1.-SIND^2)
      SC=SIND/COSD
      SDEC=ATAN(SC)
      SRASN=3.141592654D0-atan(COS(OBLIQ)/SOB*SC,-COS(SLP)/COSD)
      RETURN
      END

;C
;C----------------------------------------------------------------------
;C
      pro RECALC,IYR,IDAY,IHOUR,MIN,ISEC,epoch=epoch
      if n_elements(epoch) ne 0 then begin
         cdf_epoch,epoch,iyr,month,idom,ihour,min,isec,/break
         iday = jday(iyr,month,idom)
;         ical,iyr,iday,month,idom,/idoy
      endif
;TJK changed to below      if iyr lt 1900 then iyr=iyr+1900

      if iyr lt 50 then iyr = 2000 + iyr
      if iyr ge 50 and iyr lt 1900 then iyr = 1900 + iyr
      ;help,iyr,iday,ihour,min,isec
;C
;C
;C	THIS IS A MODIFIED VERSION OF THE SUBROUTINE RECOMP WRITTEN BY
;C	N.A. TSYGANENKO.  SINCE I WANT TO USE IT IN PLACE OF SUBROUTINE
;C	RECALC I HAVE RENAMED THIS ROUTINE RECALC AND ELIMINATED THE
;C	ORIGINAL RECALC FROM THIS VERSION OF THE GEOPACK.F PACKET. THIS
;C	WAY ALL ORIGINAL CALLS TO RECALC WILL CONTINUE TO WORK WITHOUT 
;C	HAVING TO CHANGE THEM TO CALLS TO RECOMP.
;C
;C
;C  AN ALTERNATIVE VERSION OF THE SUBROUTINE RECALC FROM THE GEOPACK PACKAGE
;C  BASED ON A DIFFERENT APPROACH TO DERIVATION OF ROTATION MATRIX ELEMENTS
;C
;C THIS SUBROUTINE WORKS BY 20% FASTER THAN RECALC AND IS EASIER TO UNDERSTAND
;C
;C   ################################################
;C   #  WRITTEN BY  N.A. TSYGANENKO ON DEC.1, 1991  #
;C   ################################################
;C
;C
;C	Modified by:	Mauricio Peredo
;C			Hughes STX at NASA/GSFC Code 695
;C			September 1992
;C
;C
;C	Modified to accept dates up to year 2000 and updated IGRF
;C	coeeficients for 1985.
;       
;       RCJ 02/2001 modified to accept dates up to year 2005 and
;	updated IGRF coefficients for 2000.
;C
;C   OTHER SUBROUTINES CALLED BY THIS ONE: SUN
;C
;C     IYR = YEAR NUMBER (FOUR DIGITS)
;C     IDAY = DAY OF YEAR (DAY 1 = JAN 1)
;C     IHOUR = HOUR OF DAY (00 TO 23)
;C     MIN = MINUTE OF HOUR (00 TO 59)
;C     ISEC = SECONDS OF DAY(00 TO 59)
;C
;	IMPLICIT NONE

;	REAL ST0,CT0,SL0,CL0,CTCL,STCL,CTSL,STSL,SFI,CFI,SPS,CPS,
;     1       SHI,CHI,HI,PSI,XMUT,A11,A21,A31,A12,A22,A32,A13,A23,
;     2       A33,DS3,F2,F1,G10,G11,H11,DT,SQ,SQQ,SQR,S1,S2,
;     3       S3,CGST,SGST,DIP1,DIP2,DIP3,Y1,Y2,Y3,Y,Z1,Z2,Z3,DJ,
;     4       T,OBLIQ,DZ1,DZ2,DZ3,DY1,DY2,DY3,EXMAGX,EXMAGY,EXMAGZ,
;     5       EYMAGX,EYMAGY,GST,SLONG,SRASN,SDEC,BA(8)

;	INTEGER IYR,IDAY,IHOUR,MIN,ISEC,K,IY,IDE,IYE,IPR
        
COMMON GEI_GSE,S1,S2,S3,DY1,DY2,DY3,DZ1,DZ2,DZ3 
COMMON C1, ST0,CT0,SL0,CL0,CTCL,STCL,CTSL,STSL,SFI,CFI,SPS,CPS, $
SHI,CHI,HI,PSI,XMUT,A11,A21,A31,A12,A22,A32,A13,A23,A33,DS3, $
K,IY,BA

common mag_dipole,g10,g11,h11
;C
;      DATA IYE,IDE,IPR/3*0/

COMMON SAVE,IYE,IDE,IPR

IF N_ELEMENTS(IYE) EQ 0 THEN IYE = 0l
IF N_ELEMENTS(IDE) EQ 0 THEN IDE = 0l
IF N_ELEMENTS(IPR) EQ 0 THEN IPR = 0l

FORMAT10 ='(//1X,"****RECOMP WARNS:  YEAR IS OUT OF INTERVAL 1965-2015: IYR=",I4,/,6X,"CALCULATIONS WILL BE DONE FOR IYR=",I4,/)'

;C
;C  IYE AND IDE ARE THE CURRENT VALUES OF YEAR AND DAY NUMBER
;C
      IY=IYR
      IDE=IDAY
      IF ((IYR EQ IYE) AND (IDAY EQ IDE))THEN GOTO, JUMP5
      IF(IY LT 1965)THEN IY=1965
      IF(IY GT 2015)THEN IY=2015
;C
;C  WE ARE RESTRICTED BY THE INTERVAL 1965-2005,
;C  FOR WHICH THE IGRF COEFFICIENTS ARE KNOWN; IF IYR IS OUTSIDE THIS INTERVAL
;C   THE SUBROUTINE GIVES A WARNING (BUT DOES NOT REPEAT IT AT THE NEXT CALLS)
;C
      IF(IY NE IYR) AND (IPR EQ 0)THEN PRINT,FORMAT=FORMAT10,IYR,IY
      IF(IY NE IYR)THEN IPR=1
      IYE=IY
;C
;C  LINEAR INTERPOLATION OF THE GEODIPOLE MOMENT COMPONENTS BETWEEN THE
;C  VALUES FOR THE NEAREST EPOCHS:
;C
	IF (IY LT 1970) THEN BEGIN			;!1965-1970
	   F2=(IY + IDAY/365D0-1965.)/5.
	   F1=1.D0-F2
	   G10=30334.*F1+30220.*F2
	   G11=-2119.*F1-2068.*F2
	   H11=5776.*F1+5737.*F2
	ENDIF ELSE IF (IY LT 1975) THEN BEGIN		;!1970-1975
	   F2=(IY+IDAY/365D0-1970.)/5.
	   F1=1.D0-F2
	   G10=30220.*F1+30100.*F2
	   G11=-2068.*F1-2013.*F2
	   H11=5737.*F1+5675.*F2
	ENDIF ELSE IF (IY LT 1980) THEN BEGIN		;!1975-1980
	   F2=(IY+IDAY/365D0-1975)/5.
	   F1=1.D0-F2
	   G10=30100.*F1+29992.*F2
	   G11=-2013.*F1-1956.*F2
	   H11=5675.*F1+5604.*F2
 	ENDIF ELSE IF (IY LT 1985) THEN BEGIN		;!1980-1985
	   F2=(IY+IDAY/365D0-1980.)/5.
	   F1=1.D0-F2
	   G10=29992.*F1+29873.*F2
	   G11=-1956.*F1-1905.*F2
	   H11=5604.*F1+5500.*F2
	ENDIF ELSE IF (IY LT 1990) THEN BEGIN		;!1985-1990
	   F2=(IY+IDAY/365D0-1985.)/5.
	   F1=1.D0-F2
	   G10=29873.*F1+29775.*F2
	   G11=-1905.*F1-1848.*F2
	   H11=5500.*F1+5406.*F2
	ENDIF ELSE IF (IY LT 1995) THEN BEGIN		;!1990-1995
;TJK added use of FLOAT	   F2=(IY+IDAY/365D0-1990.)/5.
           F2=(FLOAT(IY)+FLOAT(IDAY)/365.-1990.)/5.
	   F1=1.D0-F2
	   G10=29775.*F1+29682.*F2
	   G11=-1848.*F1-1789.*F2
	   H11=5406.*F1+5318.*F2
	ENDIF ELSE IF (IY LT 2000) THEN BEGIN		;!1995-2000
;TJK added use of FLOAT	   F2=(IY+IDAY/365D0-1995.)/5.
           F2=(FLOAT(IY)+FLOAT(IDAY)/365.-1995.)/5.
	   F1=1.D0-F2
	   G10=29682.*F1+29615.*F2
	   G11=-1789.*F1-1728.*F2
	   H11=5318.*F1+5186.*F2
        ENDIF ELSE IF (IY LT 2005) THEN BEGIN           ;!2000-2005
           ; Obtained from SWMF/share/Library/src/CON_geopack.f90
           F2=(FLOAT(IY)+FLOAT(IDAY)/365.-2000.)/5.
           F1=1.D0-F2
           G10=29615.*F1+29554.63*F2
           G11=-1728.*F1-1669.05*F2
           H11= 5186.*F1+5077.99*F2
        ENDIF ELSE IF (IY LT 2010) THEN BEGIN           ;!2005-2010
           ; Obtained from SWMF/share/Library/src/CON_geopack.f90
           F2=(FLOAT(IY)+FLOAT(IDAY)/365.-2005.)/5.
           F1=1.D0-F2
           G10=29554.63*F1+29496.5*F2
           G11=-1669.05*F1-1585.9*F2
           H11= 5077.99*F1+4945.1*F2
;
;   LINEAR EXTRAPOLATION BEYOND 2010 BY USING SECULAR VELOCITY COEFFICIENTS:
;      Use coefficient of F2 above and DT*(Coeff_F2-Coeff_F1)/5

	ENDIF ELSE BEGIN				;!2010-2015
           DT  = FLOAT(IY) + FLOAT(IDAY)/365.-2010.
	   G10 =  29496.5 - 11.6*DT
	   G11 =  -1585.9 + 16.6*DT
	   H11 =   4945.1 - 26.6*DT
	ENDELSE
;C
;C  NOW CALCULATE THE COMPONENTS OF THE UNIT VECTOR EzMAG IN GEO COORD.SYSTEM:
;C   SIN(TETA0)*COS(LAMBDA0), SIN(TETA0)*SIN(LAMBDA0), AND COS(TETA0)
;C         ST0 * CL0                ST0 * SL0                CT0
;C
      SQ=G11^2+H11^2
      SQQ=SQRT(SQ)
      SQR=SQRT(G10^2+SQ)
      SL0=-H11/SQQ
      CL0=-G11/SQQ
      ST0=SQQ/SQR
      CT0=G10/SQR
      STCL=ST0*CL0
      STSL=ST0*SL0
      CTSL=CT0*SL0
      CTCL=CT0*CL0
;C
;C      THE CALCULATIONS ARE TERMINATED IF ONLY GEO-MAG TRANSFORMATION
;C       IS TO BE DONE  (IHOUR>24 IS THE AGREED CONDITION FOR THIS CASE):
;C
JUMP5:   IF (IHOUR GT 24)THEN RETURN
;C
       SUN,IY,IDAY,IHOUR,MIN,ISEC,GST,SLONG,SRASN,SDEC
;C
;C  S1,S2, AND S3 ARE THE COMPONENTS OF THE UNIT VECTOR EXGSM=EXGSE IN THE
;C   SYSTEM GEI POINTING FROM THE EARTH'S CENTER TO THE SUN:
;C
      S1=COS(SRASN)*COS(SDEC)
      S2=SIN(SRASN)*COS(SDEC)
      S3=SIN(SDEC)
      CGST=COS(GST)
      SGST=SIN(GST)
;C
;C  DIP1, DIP2, AND DIP3 ARE THE COMPONENTS OF THE UNIT VECTOR EZSM=EZMAG
;C   IN THE SYSTEM GEI:
;C
      DIP1=STCL*CGST-STSL*SGST
      DIP2=STCL*SGST+STSL*CGST
      DIP3=CT0
;C
;C  NOW CALCULATE THE COMPONENTS OF THE UNIT VECTOR EYGSM IN THE SYSTEM GEI
;C   BY TAKING THE VECTOR PRODUCT D x S AND NORMALIZING IT TO UNIT LENGTH:
;C
      Y1=DIP2*S3-DIP3*S2
      Y2=DIP3*S1-DIP1*S3
      Y3=DIP1*S2-DIP2*S1
      Y=SQRT(Y1*Y1+Y2*Y2+Y3*Y3)
      Y1=Y1/Y
      Y2=Y2/Y
      Y3=Y3/Y
;C
;C   THEN IN THE GEI SYSTEM THE UNIT VECTOR Z = EZGSM = EXGSM x EYGSM = S x Y
;C    HAS THE COMPONENTS:
;C
      Z1=S2*Y3-S3*Y2
      Z2=S3*Y1-S1*Y3
      Z3=S1*Y2-S2*Y1
;C
;C    THE VECTOR EZGSE (HERE DZ) IN GEI HAS THE COMPONENTS (0,-SIN(DELTA),
;C     COS(DELTA)) = (0.,-0.397823,0.917462); HERE DELTA = 23.44214 DEG FOR
;C   THE EPOCH 1978 (SEE THE BOOK BY GUREVICH OR OTHER ASTRONOMICAL HANDBOOKS).
;C    HERE THE MOST ACCURATE TIME-DEPENDENT FORMULA IS USED:
;C
      DJ=365*(IY-1900)+FIX((IY-1901)/4) +IDAY -0.5 +ISEC/86400D0
      T=DJ/36525D0
      OBLIQ=(23.45229-0.0130125*T)/57.2957795
      DZ1=0.
      DZ2=-SIN(OBLIQ)
      DZ3=COS(OBLIQ)
;C
;C  THEN THE UNIT VECTOR EYGSE IN GEI SYSTEM IS THE VECTOR PRODUCT DZ x S :
;C
      DY1=DZ2*S3-DZ3*S2
      DY2=DZ3*S1-DZ1*S3
      DY3=DZ1*S2-DZ2*S1
;C
;C   THE ELEMENTS OF THE MATRIX GSE TO GSM ARE THE SCALAR PRODUCTS:
;C   CHI=EM22=(EYGSM,EYGSE), SHI=EM23=(EYGSM,EZGSE), EM32=(EZGSM,EYGSE)=-EM23,
;C     AND EM33=(EZGSM,EZGSE)=EM22
;C
      CHI=Y1*DY1+Y2*DY2+Y3*DY3
      SHI=Y1*DZ1+Y2*DZ2+Y3*DZ3
      HI=ASIN(SHI)
;C
;C    TILT ANGLE: PSI=ARCSIN(DIP,EXGSM)
;C
      SPS=DIP1*S1+DIP2*S2+DIP3*S3
      CPS=SQRT(1.-SPS^2)
      PSI=ASIN(SPS)
;C
;C    THE ELEMENTS OF THE MATRIX MAG TO SM ARE THE SCALAR PRODUCTS:
;C CFI=GM22=(EYSM,EYMAG), SFI=GM23=(EYSM,EXMAG); THEY CAN BE DERIVED AS FOLLOWS:
;C
;C IN GEO THE VECTORS EXMAG AND EYMAG HAVE THE COMPONENTS (CT0*CL0,CT0*SL0,-ST0)
;C  AND (-SL0,CL0,0), RESPECTIVELY.    HENCE, IN GEI THE COMPONENTS ARE:
;C  EXMAG:    CT0*CL0*COS(GST)-CT0*SL0*SIN(GST)
;C            CT0*CL0*SIN(GST)+CT0*SL0*COS(GST)
;C            -ST0
;C  EYMAG:    -SL0*COS(GST)-CL0*SIN(GST)
;C            -SL0*SIN(GST)+CL0*COS(GST)
;C             0
;C  THE COMPONENTS OF EYSM IN GEI WERE FOUND ABOVE AS Y1, Y2, AND Y3;
;C  NOW WE ONLY HAVE TO COMBINE THE QUANTITIES INTO SCALAR PRODUCTS:
;C
      EXMAGX=CT0*(CL0*CGST-SL0*SGST)
      EXMAGY=CT0*(CL0*SGST+SL0*CGST)
      EXMAGZ=-ST0
      EYMAGX=-(SL0*CGST+CL0*SGST)
      EYMAGY=-(SL0*SGST-CL0*CGST)
      CFI=Y1*EYMAGX+Y2*EYMAGY
      SFI=Y1*EXMAGX+Y2*EXMAGY+Y3*EXMAGZ
;C
      XMUT=(atan(SFI,CFI)+3.1415926536D0)*3.8197186342D0
;C
;C  THE ELEMENTS OF THE MATRIX GEO TO GSM ARE THE SCALAR PRODUCTS:
;C
;C   A11=(EXGEO,EXGSM), A12=(EYGEO,EXGSM), A13=(EZGEO,EXGSM),
;C   A21=(EXGEO,EYGSM), A22=(EYGEO,EYGSM), A23=(EZGEO,EYGSM),
;C   A31=(EXGEO,EZGSM), A32=(EYGEO,EZGSM), A33=(EZGEO,EZGSM),
;C
;C   ALL THE UNIT VECTORS IN BRACKETS ARE ALREADY DEFINED IN GEI:
;C
;C  EXGEO=(CGST,SGST,0), EYGEO=(-SGST,CGST,0), EZGEO=(0,0,1)
;C  EXGSM=(S1,S2,S3),  EYGSM=(Y1,Y2,Y3),   EZGSM=(Z1,Z2,Z3)
;C                                                           AND  THEREFORE:
;C
      A11=S1*CGST+S2*SGST
      A12=-S1*SGST+S2*CGST
      A13=S3
      A21=Y1*CGST+Y2*SGST
      A22=-Y1*SGST+Y2*CGST
      A23=Y3
      A31=Z1*CGST+Z2*SGST
      A32=-Z1*SGST+Z2*CGST
      A33=Z3
;C
      RETURN
      END

;C--------------------------------------------------------------------------
;C
      pro GEOMAG,XGEO,YGEO,ZGEO,XMAG,YMAG,ZMAG,J,IYR,skip=skip
;C
;C CONVERTS GEOGRAPHIC (GEO) TO DIPOLE (MAG) COORDINATES OR VICA VERSA.
;C IYR IS YEAR NUMBER (FOUR DIGITS).
;C
;C                           J>0                J<0
;C-----INPUT:  J,XGEO,YGEO,ZGEO,IYR   J,XMAG,YMAG,ZMAG,IYR
;C-----OUTPUT:    XMAG,YMAG,ZMAG        XGEO,YGEO,ZGEO
;C
;C
;C                   AUTHOR: NIKOLAI A. TSYGANENKO
;C                           INSTITUTE OF PHYSICS
;C                           ST.-PETERSBURG STATE UNIVERSITY
;C                           STARY PETERGOF 198904
;C                           ST.-PETERSBURG
;C                           RUSSIA
;C
;	IMPLICIT NONE

;	REAL XGEO,YGEO,ZGEO,XMAG,YMAG,ZMAG,ST0,CT0,SL0,CL0,CTCL,
;     1       STCL,CTSL,STSL,AB(19),BB(8)

;	INTEGER J,IYR,K,IY,II
COMMON C1, ST0,CT0,SL0,CL0,CTCL,STCL,CTSL,STSL,SFI,CFI,SPS,CPS, $
SHI,CHI,HI,PSI,XMUT,A11,A21,A31,A12,A22,A32,A13,A23,A33,DS3, $
K,IY,BA
;COMMON C1, ST0,CT0,SL0,CL0,CTCL,STCL,CTSL,STSL,AB,K,IY,BB
COMMON SAVE2,II
;if keyword_set(skip) eq 0 then begin
;      IF N_ELEMENTS(II) EQ 0 THEN II=1
;      IF(IYR EQ II)THEN GOTO, JUMP1
;      II=IYR
;      RECALC,II,0,25,0,0
;endif

JUMP1: 
      IF(J LT 0)THEN GOTO, JUMP2
      XMAG=XGEO*CTCL+YGEO*CTSL-ZGEO*ST0
      YMAG=YGEO*CL0-XGEO*SL0
      ZMAG=XGEO*STCL+YGEO*STSL+ZGEO*CT0
      RETURN
JUMP2:   XGEO=XMAG*CTCL-YMAG*SL0+ZMAG*STCL
      YGEO=XMAG*CTSL+YMAG*CL0+ZMAG*STSL
      ZGEO=ZMAG*CT0-XMAG*ST0
      RETURN
      END

;C-------------------------------------------------------------------------
;C
       pro GSMGSE,XGSM,YGSM,ZGSM,XGSE,YGSE,ZGSE,J,epoch=epoch
       common gsmgsesave,epoch0
       if n_elements(epoch) ne 0 then begin
          if n_elements(epoch0) eq 0 then epoch0 = 0d0
          if epoch ne epoch0 then recalc,epoch=epoch
          epoch0 = epoch
       endif
;C
;C CONVERTS SOLAR MAGNETOSPHERIC (GSM) TO SOLAR ECLIPTICAL (GSE) COORDS
;C   OR VICA VERSA.
;C                    J>0                J<0
;C-----INPUT: J,XGSM,YGSM,ZGSM    J,XGSE,YGSE,ZGSE
;C----OUTPUT:   XGSE,YGSE,ZGSE      XGSM,YGSM,ZGSM
;C  ATTENTION:  SUBROUTINE  RECALC  MUST BE CALLED BEFORE GSMGSE IN TWO CASES:
;C     /A/  BEFORE THE FIRST CALL OF GSMGSE
;C     /B/  IF THE CURRENT VALUES OF IYEAR,IDAY,IHOUR,MIN,ISEC ARE DIFFERENT
;C          FROM THOSE IN THE PRECEDING  CALL  OF GSMGSE
;C
;C
;C                   AUTHOR: NIKOLAI A. TSYGANENKO
;C                           INSTITUTE OF PHYSICS
;C                           ST.-PETERSBURG STATE UNIVERSITY
;C                           STARY PETERGOF 198904
;C                           ST.-PETERSBURG
;C                           RUSSIA
;C
;	IMPLICIT NONE

;	REAL XGSM,YGSM,ZGSM,XGSE,YGSE,ZGSE,SHI,CHI,
;     1       A(12),AB(13),BA(8)
;	INTEGER J,K,IY

COMMON C1, ST0,CT0,SL0,CL0,CTCL,STCL,CTSL,STSL,SFI,CFI,SPS,CPS, $
SHI,CHI,HI,PSI,XMUT,A11,A21,A31,A12,A22,A32,A13,A23,A33,DS3, $
K,IY,BA

;COMMON C1, A,SHI,CHI,AB,K,IY,BA

      IF(J LT 0)THEN GOTO, JUMP1
      XGSE=XGSM
      YGSE=YGSM*CHI-ZGSM*SHI
      ZGSE=YGSM*SHI+ZGSM*CHI
      RETURN
JUMP1:     XGSM=XGSE
      YGSM=YGSE*CHI+ZGSE*SHI
      ZGSM=ZGSE*CHI-YGSE*SHI
      RETURN
      END

;C---------------------------------------------------------------------------
;C
      PRO GEOGSM,XGEO,YGEO,ZGEO,XGSM,YGSM,ZGSM,J
;C
;C CONVERTS GEOGRAPHIC TO SOLAR MAGNETOSPHERIC COORDINATES OR VICA VERSA.
;C
;C                   J>0                   J<0
;C----- INPUT:  J,XGEO,YGEO,ZGEO    J,XGSM,YGSM,ZGSM
;C---- OUTPUT:    XGSM,YGSM,ZGSM      XGEO,YGEO,ZGEO
;C  ATTENTION:  SUBROUTINE  RECALC  MUST BE CALLED BEFORE GEOGSM IN TWO CASES:
;C     /A/  BEFORE THE FIRST USE OF GEOGSM
;C     /B/  IF THE CURRENT VALUES OF IYEAR,IDAY,IHOUR,MIN,ISEC  ARE DIFFERENT
;C            FROM THOSE IN THE PREVIOUS CALL OF THIS SUBROUTINE
;C
;C
;C                   AUTHOR: NIKOLAI A. TSYGANENKO
;C                           INSTITUTE OF PHYSICS
;C                           ST.-PETERSBURG STATE UNIVERSITY
;C                           STARY PETERGOF 198904
;C                           ST.-PETERSBURG
;C                           RUSSIA
;C
;	IMPLICIT NONE

;	REAL XGEO,YGEO,ZGEO,XGSM,YGSM,ZGSM,A11,A21,A31,A12,
;     1       A22,A32,A13,A23,A33,D,AA(17),B(8)
;	INTEGER J,K,IY
COMMON C1, ST0,CT0,SL0,CL0,CTCL,STCL,CTSL,STSL,SFI,CFI,SPS,CPS, $
SHI,CHI,HI,PSI,XMUT,A11,A21,A31,A12,A22,A32,A13,A23,A33,DS3, $
K,IY,BA


;      COMMON C1, AA,A11,A21,A31,A12,A22,A32,A13,A23,A33,D,K,IY,B
      IF (J LT 0)THEN GOTO, JUMP1
      XGSM=A11*XGEO+A12*YGEO+A13*ZGEO
      YGSM=A21*XGEO+A22*YGEO+A23*ZGEO
      ZGSM=A31*XGEO+A32*YGEO+A33*ZGEO
      RETURN
JUMP1:   XGEO=A11*XGSM+A21*YGSM+A31*ZGSM
      YGEO=A12*XGSM+A22*YGSM+A32*ZGSM
      ZGEO=A13*XGSM+A23*YGSM+A33*ZGSM
      RETURN
      END

;C--------------------------------------------------------------------------
;C
pro SMGSM,XSM,YSM,ZSM,XGSM,YGSM,ZGSM,J
;C
;C CONVERTS SOLAR MAGNETIC (SM) TO SOLAR MAGNETOSPHERIC (GSM) COORDINATES
;C   OR VICA VERSA.
;C                  J>0                 J<0
;C-----INPUT: J,XSM,YSM,ZSM        J,XGSM,YGSM,ZGSM
;C----OUTPUT:  XGSM,YGSM,ZGSM       XSM,YSM,ZSM
;C
;C  ATTENTION:  SUBROUTINE RECALC MUST BE CALLED BEFORE SMGSM IN TWO CASES:
;C     /A/  BEFORE THE FIRST USE OF SMGSM
;C     /B/  IF THE CURRENT VALUES OF IYEAR,IDAY,IHOUR,MIN,ISEC ARE DIFFERENT
;C          FROM THOSE IN THE PRECEDING CALL OF SMGSM
;C
;C
;C                AUTHOR: NIKOLAI A. TSYGANENKO
;C                        INSTITUTE OF PHYSICS
;C                        ST.-PETERSBURG  UNIVERSITY
;C                        STARY PETERGOF 198904
;C                        ST.-PETERSBURG
;C                        U.S.S.R.
;C
;	IMPLICIT NONE

;	REAL XSM,YSM,ZSM,XGSM,YGSM,ZGSM,SPS,CPS,A(10),B(15),AB(8)
;	INTEGER J,K,IY
;J=0l&k=0l&iy=0l
COMMON C1, ST0,CT0,SL0,CL0,CTCL,STCL,CTSL,STSL,SFI,CFI,SPS,CPS, $
SHI,CHI,HI,PSI,XMUT,A11,A21,A31,A12,A22,A32,A13,A23,A33,DS3, $
K,IY,BA
;      COMMON C1,A,SPS,CPS,B,K,IY,AB
      IF (J LT 0)THEN GOTO, jump1
      XGSM=XSM*CPS+ZSM*SPS
      YGSM=YSM
      ZGSM=ZSM*CPS-XSM*SPS
      RETURN
 jump1:   XSM=XGSM*CPS-ZGSM*SPS
      YSM=YGSM
      ZSM=XGSM*SPS+ZGSM*CPS
      RETURN
      END

pro geigse,xgei,ygei,zgei,xgse,ygse,zgse,j,epoch
;c added to geopack by SAB 9/9/94 to transform between
;c gei and gse coordinate system
;c j>0  gei->gse, j<0 gse > gei
;c xgse,ygse,zgse axis unit vectors in gei system
common gei_gse,x1,x2,x3,y1,y2,y3,z1,z2,z3
;c	print,' in geigse:',x1,x2,x3,y1,y2,y3,z1,z2,z3
if( j gt 0 ) then begin
	xgse = x1*xgei + x2*ygei + x3*zgei
	ygse = y1*xgei + y2*ygei + y3*zgei
	zgse = z1*xgei + z2*ygei + z3*zgei
endif else begin
	xgei = x1*xgse + y1*ygse + z1*zgse
	ygei = x2*xgse + y2*ygse + z2*zgse
	zgei = x3*xgse + y3*ygse + z3*zgse
endelse
end

pro geigeo,xgei,ygei,zgei,xgeo,ygeo,zgeo,j,epoch=epoch
SUN,IYR,IDAY,IHOUR,MIN,ISEC,GST,SLONG,SRASN,SDEC,epoch=epoch  
;help,gst
a = cos(gst)
b = sin(gst)
if( j gt 0 ) then begin
	xgeo =  a*xgei + b*ygei
	ygeo = -b*xgei + a*ygei
	zgeo =  zgei
endif else begin
	xgei = a*xgeo - b*ygeo
	ygei = b*xgeo + a*ygeo
	zgei = zgeo
endelse
end

pro scgei,xgei,ygei,zgei,xy_sunang,xsc=xsc,ysc=ysc,zsc=zsc,colat=colat
;compute direction lepedea is pointing in GEI
; or a vector direction in s/c coordinates to GEI
; colat is the angle between the vector and the s/c z-axis in degrees
common gei_gse,s1,s2,s3,dy1,dy2,dy3,dz1,dz2,dz3   
vsun = [s1,s2,s3] ;sun position in GEI ascension and declination
; declination and ascen of s/c spin axis in GEI
dec = 7.8*!pi/180 ;orientaiton of Hawkeye spin axis in GEI
ascen = 299.4*!pi/180
; compute unit vector along spin axis in GEI denoted by zsc
zsc = [cos(dec)*cos(ascen),cos(dec)*sin(ascen),sin(dec)]
dot = total(zsc*vsun)
; compute despun s/c unit vectors in GEI (xsc,ysc,zsc)
; sun lies in plane defined by xsc,zsc the sun vector projection
; xsc is positive toward the sun
; compute unit vector ysc which lies in plane defined by crossp(zsc,xsc)
xsc = vsun - dot*zsc &xsc = xsc/sqrt(total(xsc^2))
ysc = crossp(zsc,xsc)
; note xsc,ysc,zsc are orthogonal unit vectors in GEI coordinates
; that define the despun s/c coordinate system
cosa=cos(xy_sunang/!radeg)
sina=sin(xy_sunang/!radeg)
;u= xsc*cos(xy_sunang/!radeg)+ysc*sin(xy_sunang/!radeg)
;xgei = u(0)
;ygei = u(1)
;zgei = u(2)
if n_elements(colat) ne 0 then begin
	sint = sin(colat/!radeg)
	cost = cos(colat/!radeg)
	xgei = (xsc(0)*cosa + ysc(0)*sina)*sint + zsc(0)*cost
	ygei = (xsc(1)*cosa + ysc(1)*sina)*sint + zsc(1)*cost
	zgei = (xsc(2)*cosa + ysc(2)*sina)*sint + zsc(2)*cost
endif else begin
	xgei = xsc(0)*cosa + ysc(0)*sina
	ygei = xsc(1)*cosa + ysc(1)*sina
	zgei = xsc(2)*cosa + ysc(2)*sina
endelse
end

; Added S. Boardsen's terminator code   4/98  RTB
;
pro con_basis,v,x,y,z
;construct orthogonal basis from vector v
;such that x lies in the v direction
x = float(v)
x = x/norm(x)
du = max( abs(x),i)
y = fltarr(3) & z= y
y((i+1) mod 3) = x(i)
y( i) = -x( (i+1) mod 3 )
y = y/norm(y)
z = crossp(x,y)
z = z/norm(z)
;print,'x*y:',total(x*y)
;print,'z*y:',total(z*y)
;print,'z*x:',total(z*x)
end

pro cart_polar,x,y,z,r,t,p,i,radians = rad, degrees = deg
if keyword_set(deg) then f =!radeg else f = 1.
if n_elements(i) eq 0 then i=1
if (i gt 0) then begin
        r=sqrt(x^2+y^2+z^2)
        t=acos(z/r)*f
        p=atan(y,x)*f
endif  else begin
        x = r*sin(t/f)*cos(p/f)
        y = r*sin(t/f)*sin(p/f)
        z = r*cos(t/f)
endelse
end

pro termcorr,lon,lat,lat0,lon0,corr
colat0 = (90.-lat0)/!radeg
cost0 = cos(colat0)
sint0 = sin(colat0)
cosp = cos( (lon-lon0)/!radeg)
sinp = sin( (lon-lon0)/!radeg)
colat = (90.-lat)/!radeg
cost = cos(colat)
sint = sin(colat)
at = cost*sint0*cosp - sint*cost0
ap = -sint*sint0*sinp
it = where( abs(at) gt abs(ap) )
ip = where( abs(at) le abs(ap) )
if it(0) ne -1 then lat(it) = lat(it) -  corr/at(it)
if ip(0) ne -1 then lon(ip) = lon(ip) +  corr/ap(ip)
end


function terminator,lat0,lon0,check=check,corr=corr, rad = rad
; if check set then confirm orthogonality of terimator to (lat0,lon0) 
; given perpendicular to terminator by (lat0,lon0) in degrees
; compute the terminator
n=100
if n_elements(rad) ne 0 then theta =  acos(1./rad) else theta = !pi/2. 
z0 = cos(theta)
rho = sin(theta) 
;help,z0,rho
a = indgen(n)*2*!pi/n
cosa = rho*cos(a)
sina = rho*sin(a)
colat = (90.-lat0)/!radeg
v = [sin(colat)*cos(lon0/!radeg),sin(colat)*sin(lon0/!radeg),cos(colat)]
con_basis,v,zb,xb,yb
x = z0*zb(0) + cosa*xb(0) + sina*yb(0) 
y = z0*zb(1) + cosa*xb(1) + sina*yb(1) 
z = z0*zb(2) + cosa*xb(2) + sina*yb(2) 
;help,x,y,z
cart_polar,x,y,z,r,theta,lon,/degree
lat = 90. - theta
if keyword_set(corr) then termcorr,lon,lat,lat0,lon0,corr  
if keyword_set(check) then begin
  for i=0,n-1 do begin
    print, total(v*[x(i),y(i),z(i)])
  endfor
endif
return, {lat:lat,lon:lon,lat0:lat0,lon0:lon0,v:[[x],[y],[z]],u:v}
end

function terminator1,lat0,lon0,check=check,corr=corr, rad = rad
; if check set then confirm orthogonality of terimator to (lat0,lon0)
; given perpendicular to terminator by (lat0,lon0) in degrees
; compute the terminator
n=200
if n_elements(rad) ne 0 then begin
   theta =  acos(1/rad)
   z0 = cos(theta)
   rho = sin(theta)
endif else begin
   z0 = 0.
   rho = 1.
endelse
;help,z0
a = indgen(n)*2*!pi/n
xcir = rho*cos(a)
ycir = rho*sin(a)
zcir = fltarr(n) + z0
colat = (90.-lat0)/!radeg
cost = cos(colat)
sint = sin(colat)
cosp = cos(lon0/!radeg)
sinp = sin(lon0/!radeg)
x = xcir*cost*cosp - ycir*sinp + zcir*sint*cosp
y = xcir*cost*sinp + ycir*cosp + zcir*sint*sinp
z = zcir*cost - xcir*sint
lon = atan(y,x)*!radeg
lat = asin( z/sqrt(x^2+y^2+z^2) )*!radeg
u = [sint*cosp,sint*sinp,cost]
if keyword_set(corr) then termcorr,lon,lat,lat0,lon0,corr
if keyword_set(check) then begin
  for i=0,n-1 do begin
    print, total(u*[x(i),y(i),z(i)])
  endfor
endif

return, {lat:lat,lon:lon,lat0:lat0,lon0:lon0,v:[[x],[y],[z]],u:u}
end
; End terminator code
