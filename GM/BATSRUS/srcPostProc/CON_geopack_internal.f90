!^CFG COPYRIGHT UM
!BOP -------------------------------------------------------------------
!
!MODULE: CON_geopack - geopack code + heliocentric systems
!
!DESCRIPTION:
!
!Contains some subroutines of the geopack code (by N.V.Tsyganenko), 
!rewritten as the .f90 module procedures. 
!Added procedures: JulianDay(A.Ridley)and a  computation for 
!the coordinate transformation like HGI=>other systems
!In the comments and for newly introduced matrices we follow a paper by 
!M.Franz and D.Harper, Planetary and Space Scieance, V.50, 217ff(2002),
!Also see the corrected version in 
!http://www.space-plasma.qmul.ac.uk/heliocoords/systems2art.pdf
!Below we cite it as F&P, a numeration of pages follows the Internet version
!INTERFACE:
Module CON_geopack_internal

  use ModNumConst
 ! use CON_world,ONLY:CON_stop
  implicit none

  save

  private ! except
  

  !DESCRIPTION:
  !The use of CON_geopack is very simple and follows the original
  !design of the GEOPACK package. 
  !Call CON_recalc(iYear,iMonth,iDay,iHour,iMin,iSec).
  !After this a lot of global and publically accessible variables
  !and transformation matrices are constructed, related to this
  !instant of time.
  ! 
  !While using the transformation matrix, we denote HgiGse so that
  !${\bf r}_{HGI}=HgiGse_DD\cdot${\bf r}_{GSE}
  !In F&H they use T(GSE,HGI) for such the matrix
  !EOP

  !================================================================
  !Slowly varying parameters and matrices, are calculated once 

  real,dimension(3,3)::HaeGei_DD !Transformation matrix T(GEI,HAE)

  !Parameters of the Solar Equator:
  real,parameter :: cLongAscNodeSolEquator = 75.77*cDegToRad
  ! Inclination of the solar equator on the ecliptic of date
  real,parameter :: cInclinationSolEquator = 7.25*cDegToRad
  !First and second Euler angles for transformation of 
  !HAE to HGI system of coordinates.

  real,dimension(3,3) :: HgiHae_DD 
  !Transformation matrix T(HAE,HGI)=&
  !          E(cLongAscNodeSolEquator,cInclinationSolEquator)
  real,dimension(3)::  AxisMagGeo_D
  !==================================================================  
 
  !  SunE(arth)M(oon)B(arycenter) - The distance from the Sun to
  !                                 the Earth-and-Moon barycentre
  real::SunEMBDistance
  real,dimension(3)::Sun2EarthDirHae_D,Sun2EarthDirHgi_D
  !SunEMDDistance*Sun2EarthDirHae_D[1 AU]- the Sun-Earth vector in HAE
  !SunEMDDistance*Sun2EarthDirHgi_D[1 AU]- the Sun-Earth vector in HGI


  real::HPMLongitude  
  !H(eoigraphic)P(rime)M(eridian)Longitude[rad] in HGI 

  !We may want to introduce the coordinate system which is rotated by
  !a constant angle AdvanceHgr with respect to the "true" HGR system
  !The "advanced" HGR coordinate system coinsides with that advanced 
  !in time for AdvanceHgr[rad]/(2\pi)*25.38 days
  !The feature at the solar surface, having a true heliographic 
  !equal to AdvanceHgr[deg] is in xOz plane of the anvanced Hgr
  !For all the other features, Longitude_adv=Longitude-AdvanceHgr 
  logical::DoAdvanceHgr=.false.
  real::AdvanceHgr
  
  real::STPLongitude,STPColatitude 
  !S(ub)T(errestrial)P(oint) angular coordinates [rad] in HGI

  real::CarringtonLongitude 
  !=STPLongitude-HPMLongitude - Carrington longitude 
  !By definition this is a Heliographic Longitude 
  !(in rotating H(elio)G(raphic) coordinate system)
  !of the subterrestrial point

  real,dimension(3,3)::GeiGse_DD,HgiGse_DD,GeiGsm_DD,GsmGse_DD

  real,dimension(3,3)::HgrHgi_DD 
  !Transformation matrix from HGI to Heiographic coordinates
  !Definition of Heliographic coordinates (HGC) see p.7 in F&H

  integer,parameter::x_=1,y_=2,z_=3

  !PUBLIC MEMBER FUNCTIONS
  public::HaeGei_DD !Transformation matrix T(GEI,HAE) 
  public::HgiHae_DD !Transformation matrix T(HAE,HGI)
  public::AxisMagGeo_D !Vector of the Earth magnetic axis, in GEO
  public::SunEMBDistance,Sun2EarthDirHae_D,Sun2EarthDirHgi_D
  public::HPMLongitude  
  !H(eoigraphic)P(rime)M(eridian)Longitude[rad] in HGI 
  !With the use of the advanced Hgr, the value of HPMLongitude
  !is increased by a constant value of AdvanceHgr
  public::STPLongitude,STPColatitude 
  !S(ub)T(errestrial)P(oint) angular coordinates [rad] in HGI
  public::CarringtonLongitude 
  !=HPMLongitude - Carrington longitude 
  !By definition this is a Heliographic Longitude 
  !(in rotating H(elio)G(raphic) coordinate system)
  !of the subterrestrial point
  
  !Other transformational matrices
  public::GeiGse_DD,HgiGse_DD,GeiGsm_DD,GsmGse_DD,HgrHgi_DD
  
  public::rot_matrix !Rotational matrix
  public::show_matrix!prints out a matrix
  public::eulerian_matrix !Eulerian rotational matrix
  public::JulianDay !Julian day of a year
  public::CON_sun !redundat, better do not use it outside
  public::CON_recalc !(Re)calculates all the public variables
  public::CON_test_geopack !Do test
  public::set_advanced_hgr !Introduced the advanced heliographic coordinates
  public::is_advanced_hgr  !Return .true. is hgr is advanced
  public::get_advanced_hgr_longitude !Returns the longitude of
                                     !in the advanced HGR 
contains

  !BOP ========================================================================
  !IROUTINE: rotation_matrix - rotational matrix
  !INTERFACE:
  function rot_matrix(iDir,Phi)
    !INPUT ARGUMENTS:
    integer,intent(in)::iDir  !axis of rotation (1,2,3 for x-,y-,z-axis)
    real,intent(in)::Phi !Angle of rotation.
    !RETURN VALUE:
    real             :: rot_matrix(3,3) 
    !DESCRIPTION:
    ! This is the fundamental routine that calculates the rotation
    ! matrix. If the rotation of coordinate system A by the angle
    ! Phi around one of the coordinate axis of A (generally, around
    ! a unity vector {\bf n} given in the system A) gives the coordinate 
    ! system B, then the rotation matrix is defined as the 
    ! transformation matrix from  A to B:
    ! $$ 
    ! {\bf r}_B={\bf rotmatrix}({\bf n},\Phi)\cdot{\bf r}_A
    ! $$
    ! For arbitrarily directed {\bf n} the rotmatrix is as follows
    ! $$
    ! {\bf rotmatrix}({\bf n},Phi)={\bf n}{\bf n}+
    ! cos(\Phi)({\bf I}-{\bf n}{\bf n})+
    ! sin(\Phi)e_{ijk}n_k
    ! $$
    !where (\bf I} is a unity tensor and $e_{ijk}$ is a fully asymmetric 
    !object in 3D space. Usually, the rotational matrix described is used only
    !with {\bf n} being directed along one of the coordinate axis of A.
    !In this case the rotational matrices are denoted as 
    !${\bf R}_x, ${\bf R}_y, ${\bf R}_z$ correspondingly
    !
    !EOP
    !-------------------------------------------------------------------------
    !BOC
    real,dimension(3,3),parameter::&
         I_DD=reshape((/&
         !i = 1,   2,   3,     /j
            cOne, cZero, cZero,&!1
           cZero,  cOne, cZero,&!2
           cZero, cZero,  cOne &!3
           /),(/3,3/))
    real,dimension(3,3,3),parameter::&
         E_DDD=reshape((/&
         !i = 1,   2,   3,     /j  \
           cZero, cZero, cZero,&!1  | k=1
           cZero, cZero, -cOne,&!2  |
           cZero,  cOne, cZero,&!3 / 
                     !i = 1,   2,   3,     /j  \
                       cZero, cZero,  cOne,&!1  | k=2
                       cZero, cZero, cZero,&!2  |
                       -cOne, cZero, cZero,&!3 /
                                !i = 1,   2,   3,     /j  \
                                  cZero, -cOne, cZero,&!1  | k=3
                                   cOne, cZero, cZero,&!2  |
                                  cZero, cZero, cZero &!3 /
           /),(/3,3,3/))
    !E_DDD(1,2,3)=E_DDD(2,3,1)=E(3,1,2)=1
    !E_DDD(2,1,3)=E_DDD(3,2,1)=E(1,3,2)=-1
    !Otherwise E_DDD=0
    rot_matrix=I_DD
    if(Phi/=cZero)then
       rot_matrix=cos(Phi)*rot_matrix+sin(Phi)*E_DDD(:,:,iDir)
       rot_matrix(iDir,iDir)=cOne
    end if
  end function rot_matrix
  !BOP ========================================================================
  !IROUTINE: show_matrix - writes a matrix
  !INTERFACE:
  subroutine show_matrix(M_DD)
    real,dimension(3,3),intent(in)::M_DD
    integer::i,j
    write(*,'(3(3F14.10,/))')((M_DD(i,j),j=1,3),i=1,3)
  end subroutine show_matrix
  !BOP ========================================================================
  !IROUTINE: eulerian_matrix - the rotational matrix for 3 Eulerian angles
  !INTERFACE:
  function eulerian_matrix(Omega,Theta,Psi)
    !INPUT ARGUMENTS:
    real, intent(in) :: Omega,Theta,Psi
    optional::Psi
    !RETURN VALUE:
    real             :: eulerian_matrix(3,3)
    !DESCRIPTION:
    !
    ! This is the fundamental routine that calculates the rotation
    ! matrix for three Eulerian angles (see, e.g. the Appendix A.1
    ! in  www.space-plasma.qmul.ac.uk/heliocoords/systems2art.pdf)
    ! The matrix is obtained from 3 rotations:
    ! $$
    !      E(\Omega,\theta,\psi) = R_z(\psi)\cdot R_x(\theta)\cdot R_z(Omega)
    ! $$
    !EOP
    !-------------------------------------------------------------------------
    !BOC
    if(Theta==cZero)then
       eulerian_matrix=rot_matrix(z_,Omega)
    elseif(Omega==cZero)then
       eulerian_matrix=rot_matrix(x_,Theta)
    else
       eulerian_matrix=matmul(rot_matrix(x_,Theta),rot_matrix(z_,Omega))
    end if
    if(present(Psi))&
         eulerian_matrix=matmul(rot_matrix(z_,Psi),eulerian_matrix)
    !EOC
  end function eulerian_matrix

  !----------------------------------------------------------------------
  integer function JulianDay(iYear,iMonth,iDay)
    !Coded by A.Ridley
    !Comment: valid for 1900<iYear<2100
    integer,intent(in)::iYear,iMonth,iDay
    integer, dimension(1:12),parameter :: nDayInMonth_I = (/ &
         31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
    IF(iYear.LT.1901.OR.iYear.GT.2099)then
       write(*,*)'CON_geopack ERROR: No ephemers data for the year of ',iYear
       call CON_stop('CON_geopack ERROR')
    end IF
    JulianDay=iDay
    if(iMonth>1)JulianDay=JulianDay+sum(nDayInMonth_I(1:iMonth-1))
    if(iMonth>2.and.mod(iYear,4)==0)JulianDay=JulianDay+1
  end function JulianDay
  !----------------------------------------------------------------------
  subroutine CON_sun(iYear,jDay,iHour,iMin,iSec,&
       GSTime,SunLongitude,Obliq)
    !
    !  CALCULATES FOUR QUANTITIES NECESSARY FOR COORDINATE 
    !  TRANSFORMATIONS
    !  WHICH DEPEND ON SUN POSITION (AND, HENCE, ON UNIVERSAL TIME 
    !  AND SEASON)
    !  From geopack.f by N.V.Tsyganenko
    integer,intent(in)::iYear,jDay,iHour,iMin,iSec
    !-------  INPUT PARAMETERS:
    !  iYear,jDay,iHour,iMin,iSec -  YEAR, DAY OF A YEAR, AND UNIVERSAL TIME 
    !  IN HOURS, MINUTES,
    !    AND SECONDS  (IDAY=1 CORRESPONDS TO JANUARY 1).
    !
    real::Century

    real,intent(out):: GSTime,SunLongitude,Obliq   
    !Angle between the rotation axis of the Earth
    !The rotation around X-axis by the angle Obliq transforms
    !GEI system to HAE system

    !-------  OUTPUT PARAMETERS:
    !  GSTime - GREENWICH MEAN SIDEREAL TIME
    !  SunLongitude - The Sun Longitude
    !  Obliq        - The inclination of the Earth equatorial plane
    !  ORIGINAL VERSION OF THIS SUBROUTINE HAS BEEN COMPILED FROM:
    !  RUSSELL, C.T., COSMIC ELECTRODYNAMICS, 1971, V.2, PP.184-196.
    !     LAST MODIFICATION:  JAN 5, 2001 (NO ESSENTIAL CHANGES, BUT
    !     SOME REDUNDANT STATEMENTS TAKEN OUT FROM THE PREVIOUS VERSION)
    !
    !     ORIGINAL VERSION WRITTEN BY:    Gilbert D. Mead

    real*8:: DJ,FDAY
    real::VL,G !Miscellaneous
    real,parameter::cDegToRadHere=cOne/57.295779513
    logical::DoInit=.true.
    !----------------------------------------------------------------------
    IF(iYear.LT.1901.OR.iYear.GT.2099)then
       write(*,*)'CON_geopack ERROR: No ephemers data for the year of ',iYear
       call CON_stop('CON_geopack ERROR')
    end IF
    FDAY=dble(IHOUR*3600+iMIN*60+ISEC)/86400.E0       !Fraction of a day

    DJ=365*(IYear-1900)+(IYear-1901)/4+jDAY-0.5E0+FDAY
    !This is a day of epoch Jan,1,1900 

   
    Century=DJ/36525.0
    !To calculate the slowly varying parameters   

    Obliq=(23.45229-0.0130125*Century)*cDegToRadHere
    
       !Transformation matrix for Heleocentric Aries Ecliptic system
       !see F&H, p 6, the formula for T(GEI,HAE)
       HaeGei_DD=eulerian_matrix(cZero,Obliq)

       !Transformation matrix to HGI system
       !See formula for T(HAE,HCI) at p.7 in F&H, where HCI stands for HGI
       HgiHae_DD=eulerian_matrix(&
            cLongAscNodeSolEquator,cInclinationSolEquator)


    VL=MOD(279.696678+0.9856473354*DJ,360.D0)
    !Comparing this with Table 4 at p.13 in F&H, we see that this is
    !a "mean" longitude $\lambda_{mean} of the Earth+MoonBaricenter in 
    !the Heleocentric
    !Aries Ecliptic system, with 180 degree offset. Were HAE the geocentric
    !system, moving with the Earth, VL would be a "mean" longitude of the Sun.

    GSTime=&
         MOD(279.690983+0.9856473354*DJ+360.0*FDAY+180.,360.D0)&
         *cDegToRadHere
    !The G(reenwich mean) S(idereal) Time is the angle between the meridian of 
    !Greenwich and equinox

    G=MOD(358.475845+0.985600267*DJ,360.D0)*cDegToRadHere
    !g in B&H is referred to as "mean anomaly" and equals the
    !difference between the "mean longitude" and the longitude of perihelii.
    !The following formula is for the distance from the Sun to the
    !Earth+Moon barycentre. See Eq.(36) in F&H. Not that in perihelii (g=0)
    !SunEMBDistance has a minimum
    !Added by I.Sokolov&I.Roussev, 08.17.03
    SunEMBDistance=1.000140-0.016710*cos(G)-0.000140*cos(G+G)
  
    SunLongitude=(VL+(1.91946-0.004789*Century)*SIN(G)+&
         0.020094*SIN(2.*G))*cDegToRadHere
    !$\lambda$ in Eq.(36) in F&H is used instead of SunLongitude, again there
    !is difference betwen $\lambda$ and SunLongitude is a half of a full 
    !angle. At the equinox time the Earth is in -X direction from the Sun
    !in HAE system (which is the heliocentric system with X-axis along the 
    !equinox line. So $\lambda$ of F&H is 180 degrees, while SunLongitude=0
    !Note, that in perihelii (g=0) the value of d SunLongitude/ d Time 
    !has a maximum, as it should be  

  
    IF(SunLongitude.GT.6.2831853)SunLongitude=SunLongitude-6.2831853
    IF(SunLongitude.LT.0.)SunLongitude=SunLongitude+6.2831853

   
    Sun2EarthDirHae_D(1)=-cos(SunLongitude)
    Sun2EarthDirHae_D(2)=-sin(SunLongitude)
    Sun2EarthDirHae_D(3)=cZero

    Sun2EarthDirHgi_D=matmul(HgiHae_DD,Sun2EarthDirHae_D)

    !H(eliographic)P(rime)M(eridian)Longitude
    !This is a longitude of the main geliographic meridian in the HGI
    !coordinate system: see Eq.(16) in F&H
    HPMLongitude=cTwoPi*real(mod((&
         DJ + &  !Day of epoch
         2415020.0000000D0 - &!Julian day for Jan,1,1900,NoonUT (epoch)
         2398220.0000000D0)  &!Julian day For Jan,1,1854,NoonUT 
         /25.38000000000D0,&  !Sidereal rotation period for the Sun
         1.0000000000000D0))  
    
    STPLongitude=atan2(Sun2EarthDirHgi_D(2),Sun2EarthDirHgi_D(1))
    if(STPLongitude<cZero)STPLongitude=STPLongitude+cTwoPi
   
    STPColatitude=acos(Sun2EarthDirHgi_D(3))
    
    CarringtonLongitude= STPLongitude-HPMLongitude
    if(CarringtonLongitude<cZero)&
         CarringtonLongitude=CarringtonLongitude+cTwoPi

    !We can introduce, if desired, the Advanced Heliographic System
    !
    if(DoAdvanceHgr) HPMLongitude=HPMLongitude+AdvanceHgr
  end subroutine CON_sun
  !---------------------------------------------------------------------------
  subroutine set_advanced_hgr(AdvanceHgrDeg)
    real,intent(in)::AdvanceHgrDeg
    AdvanceHgr=AdvanceHgrDeg*cDegToRad
    DoAdvanceHgr=.true.
  end subroutine set_advanced_hgr
  !---------------------------------------------------------------------------
  logical function is_advanced_hgr()
    is_advanced_hgr=DoAdvanceHgr
  end function is_advanced_hgr
  !---------------------------------------------------------------------------
  subroutine get_advanced_hgr_longitude(HGRLongitute,AdvHGRLongitute)
    real,intent(in)::HGRLongitute
    real,intent(out)::AdvHGRLongitute
    AdvHGRLongitute=HGRLongitute-AdvanceHgr*cRadToDeg
  end subroutine get_advanced_hgr_longitude
  !---------------------------------------------------------------------------
  subroutine CON_mag_axis(iYear,jDay)
    !This is a part of the RECALC subroutine form geopack.f by Tsyganenko
    integer,intent(in)::iYear,jDay
    !-----INPUT PARAMETERS:
    !
    !     IYear   -  YEAR NUMBER (FOUR DIGITS)
    !     jDAY  -  DAY OF YEAR (DAY 1 = JAN 1)
    integer::IYE=0,IDE=0,IPR=0
    integer::iY
    real::F1,F2,H11,G10,G11,SQR,DT
    IF (IYear.EQ.IYE.AND.jDAY.EQ.IDE) return
    !
    !   IYE AND IDE ARE THE CURRENT VALUES OF YEAR AND DAY NUMBER
    !
    IY=IYear
    IDE=jDAY
    IF(IY.LT.1965) IY=1965
    IF(IY.GT.2005) IY=2005
    !
    !  WE ARE RESTRICTED BY THE INTERVAL 1965-2005,
    !  FOR WHICH THE IGRF COEFFICIENTS ARE KNOWN; IF IYR IS OUTSIDE THIS INTERVAL
    !  THE SUBROUTINE PRINTS A WARNING (BUT DOES NOT REPEAT IT AT NEXT INVOCATIONS)
    !
    IF(IY.NE.IYear.AND.IPR.EQ.0)&
         WRITE (*,*) 'No Igrf Coefficients are availble for the year ',&
         IYear,'We use date for year ',IY
    IF(IY.NE.IYear) IPR=1
    IYE=IY   
    !
    !  LINEAR INTERPOLATION OF THE GEODIPOLE MOMENT COMPONENTS BETWEEN THE
    !  VALUES FOR THE NEAREST EPOCHS:
    !
    IF (IY.LT.1970) THEN                            !1965-1970
       F2=(FLOAT(IY)+FLOAT(jDAY)/365.-1965.)/5.
       F1=1.E0-F2
       G10=30334.*F1+30220.*F2 ! HERE G10 HAS OPPOSITE SIGN TO THAT IN IGRF TABLES
       G11=-2119.*F1-2068.*F2
       H11=5776.*F1+5737.*F2
    ELSEIF (IY.LT.1975) THEN                        !1970-1975
       F2=(FLOAT(IY)+FLOAT(jDAY)/365.-1970.)/5.
       F1=1.E0-F2
       G10=30220.*F1+30100.*F2 ! HERE G10 HAS OPPOSITE SIGN TO THAT IN IGRF TABLES
       G11=-2068.*F1-2013.*F2
       H11=5737.*F1+5675.*F2
    ELSEIF (IY.LT.1980) THEN                        !1975-1980
       F2=(dble(IY)+dble(jDAY)/365.-1975.)/5.
       F1=1.E0-F2
       G10=30100.*F1+29992.*F2 ! HERE G10 HAS OPPOSITE SIGN TO THAT IN IGRF TABLES
       G11=-2013.*F1-1956.*F2
       H11=5675.*F1+5604.*F2
    ELSEIF (IY.LT.1985) THEN                        !1980-1985
       F2=(FLOAT(IY)+FLOAT(jDAY)/365.-1980.)/5.
       F1=1.E0-F2
       G10=29992.*F1+29873.*F2 ! HERE G10 HAS OPPOSITE SIGN TO THAT IN IGRF TABLES
       G11=-1956.*F1-1905.*F2
       H11=5604.*F1+5500.*F2
    ELSEIF (IY.LT.1990) THEN			!1985-1990
       F2=(FLOAT(IY)+FLOAT(jDAY)/365.-1985.)/5.
       F1=1.E0-F2
       G10=29873.*F1+29775.*F2 ! HERE G10 HAS OPPOSITE SIGN TO THAT IN IGRF TABLES
       G11=-1905.*F1-1848.*F2
       H11=5500.*F1+5406.*F2
    ELSEIF (IY.LT.1995) THEN			!1990-1995
       F2=(FLOAT(IY)+FLOAT(jDAY)/365.-1990.)/5.
       F1=1.E0-F2
       G10=29775.*F1+29682.*F2 ! HERE G10 HAS OPPOSITE SIGN TO THAT IN IGRF TABLES
       G11=-1848.*F1-1789.*F2
       H11=5406.*F1+5318.*F2
    ELSEIF (IY.LT.2000) THEN                        !1995-2000
       F2=(FLOAT(IY)+FLOAT(jDAY)/365.-1990.)/5.
       F1=1.E0-F2
       G10=29682.*F1+29615.*F2 ! HERE G10 HAS OPPOSITE SIGN TO THAT IN IGRF TABLES
       G11=-1789.*F1-1728.*F2
       H11=5318.*F1+5186.*F2
    ELSE                                            !2000-2005
       !
       !   LINEAR EXTRAPOLATION BEYOND 2000 BY USING SECULAR VELOCITY COEFFICIENTS:
       !
       DT=FLOAT(IY)+FLOAT(jDAY)/366.-2000.
       G10=29615.-14.6*DT      ! HERE G10 HAS OPPOSITE SIGN TO THAT IN IGRF TABLES
       G11=-1728.+10.7*DT
       H11=5186.-22.5*DT
    ENDIF
    !
    !  NOW CALCULATE THE COMPONENTS OF THE UNIT VECTOR EzMAG IN GEO COORD.SYSTEM:
    !   SIN(TETA0)*COS(LAMBDA0), SIN(TETA0)*SIN(LAMBDA0), AND COS(TETA0)
    !         ST0 * CL0                ST0 * SL0                CT0
    !
    SQR=G11**2+H11**2+G10**2
    AxisMagGeo_D(3)=G10/SQR
    AxisMagGeo_D(1)=-G11/SQR
    AxisMagGeo_D(2)=-H11/SQR
  end subroutine CON_mag_axis
  !------------------------------------------------------------------------
  subroutine CON_recalc(iYear,iMonth,iDay,iHour,iMin,iSec)
    !Updates matrices for the coordinate transformations
    !Computations for GeiGse_DD and GeiGsm_DD are from the subroutine
    !RECALC of geopack.f by N.V.Tsyganenko
    !Computations for GeiHgi_DD - i.Roussev and I.Sokolov,
    !igorsok@umich.edu, phone (734)647-4705

    ! 3/9/2005: G.Toth - corrected HgiGse_DD calculation,
    !                    which was 180 degrees off. 
    !                    NOTE: the GeiHgi_DD is only defined in the test.

    integer,intent(in)::iYear,iMonth,iDay,iHour,iMin,iSec
    integer::jDay
    real::AxisMagGei_D(3),GSTime,SunLongitude,Obliq
    real,dimension(3,3)::GseHae_DD !Transformation matrix T(HAE,GSE)
    real,dimension(3,3)::GeoGei_DD !Transformation matrix T(GEI,GEO)

    !-------------------------------------------------------------------
    jDay=JulianDay(iYear,iMonth,iDay)
    call CON_mag_axis(iYear,jDay)
    call CON_sun(iYear,jDay,iHour,iMin,iSec,GSTime,SunLongitude,Obliq)
    
    !GseHae= T(HAE,GSE)=E(\lambda+180,0,0) see F&H, p.9 
  
    GseHae_DD=eulerian_matrix(SunLongitude-9.924E-5,cZero)
    
    !Trasnformation matrix GeiGse is calculated as 
    !transpose(T(HAE,GSE)\cdot HaeGei_DD)
  
    GeiGse_DD=transpose(matmul(GseHae_DD,HaeGei_DD))

    
    !GSTime is the Eulerian angle for the transformation
    !matrix T(GEI,GEO)
    GeoGei_DD=eulerian_matrix(GSTime,cZero)
    
    !Transformation matrix HgiGse_DD is caculated as a
    !product of the above defined matrices
    HgiGse_DD = matmul(HgiHae_DD,transpose(GseHae_DD))

    !Transformation matrix to Heliographic coordinates
    !By definition (see  F&H, p.7), T(HAE,HGC)=E(\Omega,i,w_0)
    !By definition of the Eulerian Matrix 
    !E(\Omega,i,w_0)=R_z(w_0)\cdot E(\Omega,i,0)=R_z(w_0)\cdot T(HAE,HGI)
    !(above T(HAE,HGI) was mentioned to be E(\Omega,i,0)
    !w_0 is calculated above as H(eleiographic)P(rime)M(eridian)Longitude
    !Hence
    HgrHgi_DD=rot_matrix(z_,HPMLongitude)

   
    !   THE COMPONENTS OF THE UNIT VECTOR EXGSM=EXGSE IN THE
    !   SYSTEM GEI POINTING FROM THE EARTH'S CENTER TO THE SUN:
    GeiGsm_DD(:,x_)=GeiGse_DD(:,x_)


    !   THE COMPONENTS OF THE UNIT VECTOR EZSM=EZMAG
    !   IN THE SYSTEM GEI:
    AxisMagGei_D=matmul(transpose(GeoGei_DD),AxisMagGeo_D)

    !
    !  NOW CALCULATE THE COMPONENTS OF THE UNIT VECTOR EYGSM 
    !  IN THE SYSTEM GEI BY TAKING THE VECTOR PRODUCT
    !   D x S AND NORMALIZING IT TO UNIT LENGTH:

    GeiGsm_DD(1,y_)=AxisMagGei_D(2)*GeiGsm_DD(3,x_)-&
         AxisMagGei_D(3)*GeiGsm_DD(2,x_)
    GeiGsm_DD(2,y_)=AxisMagGei_D(3)*GeiGsm_DD(1,x_)-&
         AxisMagGei_D(1)*GeiGsm_DD(3,x_)
    GeiGsm_DD(3,y_)=AxisMagGei_D(1)*GeiGsm_DD(2,x_)-&
         AxisMagGei_D(2)*GeiGsm_DD(1,x_)
    GeiGsm_DD(:,y_)=GeiGsm_DD(:,y_)/&
         sqrt(dot_product(GeiGsm_DD(:,y_),GeiGsm_DD(:,y_)))

    !
    !   THEN IN THE GEI SYSTEM THE UNIT VECTOR 
    !   Z = EZGSM = EXGSM x EYGSM = S x Y
    !   HAS THE COMPONENTS:
    GeiGsm_DD(1,z_)=GeiGsm_DD(2,x_)*GeiGsm_DD(3,y_)-&
         GeiGsm_DD(3,x_)*GeiGsm_DD(2,y_)
    GeiGsm_DD(2,z_)=GeiGsm_DD(3,x_)*GeiGsm_DD(1,y_)-&
         GeiGsm_DD(1,x_)*GeiGsm_DD(3,y_)
    GeiGsm_DD(3,z_)=GeiGsm_DD(1,x_)*GeiGsm_DD(2,y_)-&
         GeiGsm_DD(2,x_)*GeiGsm_DD(1,y_)

    GsmGse_DD=matmul(transpose(GeiGsm_DD),GeiGse_DD)

  end subroutine CON_recalc
  !----------------------------------------------------------------
  subroutine CON_test_geopack
    ! Coded by I Sokolov, and I Roussev, 08.16.2003
    ! The test compares the mean position of the pole of the solar
    ! eqautor (see SUN,2001, p.C3) with the unity vector n_z of
    ! the Heliographic inertial coordinates with respect to GEI  
    ! coordinate system
    real::GeiHgi_DD(3,3)
    real,parameter::RightAssention=286.13*cDegToRad,&
         Declination=63.87*cDegToRad,&
         cQuarterPi=0.25*cPi
    integer::iYear=2000,iMonth,iDay,iHour,iMin=0,iSec=0,iDir
    do iDir=1,3
       write(*,*)'Show a rotational matrix for Phi=Pi/4 and iDir=',iDir
       call show_matrix(rot_matrix(iDir,cQuarterPi))
       if(iDir/=2)CYCLE
       write(*,*)'Check a formula R_y(Phi)=E(90,Phi,-90):'
       write(*,*)'E(90,Pi/4,-90)'
       call show_matrix(eulerian_matrix(cPi*cHalf,cQuarterPi,-cPi*cHalf))
    end do
  
    !For perihelion
    iMonth=1;iDay=3;iHour=5
    call CON_recalc(iYear,iMonth,iDay,iHour,iMin,iSec)
    write(*,'(a,f14.10,a)')'SunEMBDistance=',SunEMBDistance,&
         ', should be 0.98329'
    GeiHgi_DD=matmul(GeiGse_DD,transpose(HgiGse_DD))
    write(*,'(a,3es16.8)')&
         'Solar rotation axis vector calculated as GeiHgi_DD(:,3)',&
         GeiHgi_DD(:,3)
    write(*,'(a,3es16.8)')&
         'The vector calculated in terms of RightAss=286.13,Declin=63.87',&
         cos(RightAssention)*cos(Declination),&
         sin(RightAssention)*cos(Declination),&
         sin(Declination)
    !For aphelion
    iMonth=7;iDay=4;iHour=0
    call CON_recalc(iYear,iMonth,iDay,iHour,iMin,iSec)
    write(*,'(a,f14.10,a)')'SunEMBDistance=',SunEMBDistance,&
         ', should be 1.01671'
    GeiHgi_DD=matmul(GeiGse_DD,transpose(HgiGse_DD))
    write(*,'(a,3es16.8)')&
         'Solar rotation axis vector calculated as GeiHgi_DD(:,3)',&
         GeiHgi_DD(:,3)
    write(*,'(a,3es16.8)')&
         'The vector calculated in terms of RightAss=286.13,Declin=63.87',&
         cos(RightAssention)*cos(Declination),&
         sin(RightAssention)*cos(Declination),&
         sin(Declination)
    iYear=2003;iMonth=10;iDay=23;iHour=13
    write(*,*)&
         'For Oct,23,2003,13.00 UT - the Carrington rotation 2009 starts:'
    call CON_recalc(iYear,iMonth,iDay,iHour,iMin,iSec)
    write(*,'(a,F14.10,a)')'HMMLongitude=',HPMLongitude*cRadToDeg,' Deg'
    write(*,'(a,2(F14.10,a))')&
         'STPLongitude=',STPLongitude*cRadToDeg, &
         ' Deg,  STPLatitude=',90.0-STPColatitude*cRadToDeg,' Deg'
    write(*,'(a,F14.10,a)')'CarringtonLongitude=',&
         CarringtonLongitude*cRadToDeg,&
         ' Deg - should be zero as the Carrington rotation starts'
    call CON_recalc(iYear,iMonth,iDay,iHour+1,iMin,iSec)
    write(*,*)'Show matrix HgcHgi_DD - should be R_z(HMMLongitude)'
    call show_matrix(HgrHgi_DD)
    write(*,*)'One hour later:'
    write(*,'(a,F14.10,a)')'HMMLongitude=',HPMLongitude*cRadToDeg,' Deg'
    write(*,'(a,2(F14.10,a))')&
         'STPLongitude=',STPLongitude*cRadToDeg,&
         ' Deg,  STPLatitude=',90.0-STPColatitude*cRadToDeg,' Deg'
    write(*,'(a,F14.10,a)')'CarringtonLongitude=',&
         CarringtonLongitude*cRadToDeg,&
         ' Deg'
  end subroutine CON_test_geopack

end Module CON_geopack_internal
!=======================================================================
subroutine CON_stop(Char)
  character(LEN=*),intent(in)::Char
  write(*,*)Char
  stop
end subroutine CON_stop
