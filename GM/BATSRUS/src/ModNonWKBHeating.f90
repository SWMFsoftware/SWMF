!HISTORY:
!
!\
!Prototype: Steve Cranmer's code 
!http://www.cfa.harvard.edu/~scranmer/Data/NonWKB/heatcvb.f
!The order of lines is changed in a couple of places.
!/
!
!\
!Implementation for the SWMF
!01/10/10 Initial revision, Igor Sokolov
!/ 

module ModNonWKBHeating
  !
  !      implicit real*8 (a-h,o-z)
  !
  !C  Define necessary constants.
  !
  !      parameter (pi=3.14159265359,fourpi=4.*pi,sq3pi=3.06998016655)
  !      parameter (xRsun=6.96d+10)
  !
  use ModConst, ONLY: rSun,cPi
  implicit none
  SAVE
  PRIVATE !Except

  !PUBLIC MEMBER:
  public:: Cranmer_heating_function


  real,parameter:: rSunCGS = rSun * 1.0e2 !Solar radius in cm
  real,parameter:: rSunCGS2= rSunCGS**2

  real,parameter:: cPi4 = 4.0 * cPi, cSqrt3Pi = 3.06998016655

  !C  Wave action conservation constant.  Earlier models found values that
  !C  ranged between:
  !C     parameter (FoBconst=2.0d+04)
  !C     parameter (FoBconst=9.0d+04)
  !C  but this paper uses an intermediate value:
  !      parameter (FoBconst=5.0d+04)

  real,parameter:: cEnergyFluxPerB = 5.0e4 !erg/cm2/s/Gs

  !C  Turbulent correlation length normalization.  Earlier models found
  !C  values that ranged from a minimum of (Cranmer et al. 2007):
  !C     parameter (ellconst=2.876d+08)
  !C  to a maximum of (Cranmer and van Ballegooijen 2005):
  !C     parameter (ellconst=1.155d+09)
  !C  and this paper uses an intermediate value:
  !      parameter (ellconst=6.000d+08)
  !
  real,parameter:: cCorrelationLength = 6.0e+8   !cm/sqrt{Gs}


  !C  Do we use the dimensionless efficiency factor?  (1=yes, 0=no)
  !      parameter (iuseeff=1)
  !
  logical:: UseEfficiency = .true.

  !C  Define frequency spectrum and weighting factors from the Alfven
  !C  wave spectrum shown in Figure 3 of Cranmer et al. (2007).
  !
  !      parameter (nfreq=17)
  !      real*8 omega(nfreq),weight(nfreq)
  !
  !      data omega /
  !     >  5.2977534d-06,1.0516068d-05,2.1127114d-05,4.1765219d-05,
  !     >  8.3292362d-05,1.6558856d-04,3.3423858d-04,6.6711177d-04,
  !     >  1.3456485d-03,2.6870358d-03,5.5235168d-03,1.1217003d-02,
  !     >  2.3297412d-02,4.0542588d-02,8.0531172d-02,1.3316122d-01,
  !     >  2.3947310d-01/
  !
  !      data weight /
  !     >  1.4449414d-05,2.6764852d-05,5.4788208d-05,9.4094896d-05,
  !     >  1.5998050d-04,2.2285727d-04,3.6475333d-04,5.7441602d-04,
  !     >  1.0809003d-03,1.9369317d-03,4.8314086d-03,2.2801884d-02,
  !     >  1.5773807d-01,3.5820958d-01,3.9913118d-01,5.2719221d-02,
  !     >  3.1574298d-05/
  !
  integer,parameter:: nFrequency=17

  real,parameter:: Omega_I(nFrequency)=(/&
       5.2977534e-06,1.0516068e-05,2.1127114e-05,4.1765219e-05,&
       8.3292362e-05,1.6558856e-04,3.3423858e-04,6.6711177e-04,&
       1.3456485e-03,2.6870358e-03,5.5235168e-03,1.1217003e-02,&
       2.3297412e-02,4.0542588e-02,8.0531172e-02,1.3316122e-01,&
       2.3947310e-01/)

  real,parameter:: Weight_I(nFrequency)=(/&
       1.4449414e-05,2.6764852e-05,5.4788208e-05,9.4094896e-05,&
       1.5998050e-04,2.2285727e-04,3.6475333e-04,5.7441602e-04,&
       1.0809003e-03,1.9369317e-03,4.8314086e-03,2.2801884e-02,&
       1.5773807e-01,3.5820958e-01,3.9913118e-01,5.2719221e-02,&
       3.1574298e-05/)


contains
  !C---------------------------------------------------------------------------
  !C
  !C  Stand-alone Fortran 77 subroutine that estimates the local volumetric
  !C  coronal heating rate, using equations from Cranmer and van Ballegooijen
  !C  (2005, ApJ Supp, 156, 265), Cranmer et al. (2007, ApJ Supp, 171, 520),
  !C  and Cranmer (2010, ApJ, in press, see also
  !   An Efficient Approximation of the Coronal Heating Rate for Use in Global 
  !   Sun-Heliosphere Simulations, arXiv:0912.5333).  All calculations are done in
  !C  cgs/Gaussian units.
  !C

  !      subroutine HEATCVB (r,uwind,rho,Bmag,rcrit,Vcrit,vperp,
  !     >  Qheat,gamma)
  !
  !C  REQUIRED INPUTS:
  !C
  !C    r:      heliocentric distance, in cm.
  !C    uwind:  local solar wind speed, in cm/s.
  !C    rho:    local mass density, in g/cm^3.
  !C    Bmag:   local magnitude of magnetic field strength, in Gauss.
  !C
  !C  OPTIONAL INPUTS:
  !C  (If any of these input values are less than or equal to zero, then
  !C  the values will be re-computed using estimates from the paper.)
  !C
  !C    rcrit:  heliocentric distance of Alfven critical point, in cm.
  !C    Vcrit:  velocity at which wind speed and Alfven speed are equal
  !C            (i.e., the Alfven critical velocity), in cm/s.
  !C    vperp:  transverse Alfven wave velocity amplitude, in cm/s.
  !C
  !C  OUTPUTS:
  !C
  !C    Qheat:  local volumetric coronal heating rate, in erg/s/cm^3.
  !C    gamma:  local damping rate of waves, in 1/s.
  !C
  !C---------------------------------------------------------------------------

  subroutine Cranmer_heating_function(&
       RSI, UMagSI, RhoSI, BMagSI,&            !INPUT
       QHeatSI,                   &            !OUTPUT
       WaveEnergyDensitySI,IsFullReflection,&  !OPTIONAL INPUT
       GammaSI)                                !OPTIONAL OUTPUT

    real, intent(in) :: &
         RSI,   &         ! heliocentric distance, in m.
         UMagSI,&         ! local solar wind speed, in  m/s.
         RhoSI, &         ! local mass density, in kg/m^3.
         BMagSI           ! local magnitude of magnetic field strength, in Tesla.

    real, intent(out):: &
         QHeatSI          ! volumetric coronal heating rate, in J/s/m^3.

    real,optional,intent(in)   :: WaveEnergyDensitySI
    logical,intent(in),optional:: IsFullReflection

    real, intent(out), optional:: &            
         GammaSI            ! local damping rate of waves, in 1/s
    !-------------------------------
    real:: RCGS, UMagCGS, RhoCGS, BMagCGS, QHeatCGS   !Converted variables

    !Elzasser variables
    real:: ZPlus, ZMinus

    !L_\perp=cCorrelationLength/sqrt(BMagCGS)
    real:: CorrelationLength

    !Parameters of the efficiency factor
    real:: tEddy, tReference

    !The turbulence amplitude:
    real:: vPerp0

    !Reflection parameters

    real::ReflectionCoef, ReflectionCoef2
    real::ReflectionCoefInf, ReflectionCoef0
    logical:: IsSetReflectionCoef

    !UA:
    real::WaveEnergyDensityCGS

    real:: VAlfven, UPlusVA, UMinusVA

    !Misc
    real:: rCrit,vCrit,vNuCrit,RhoCrit
    real:: Alpha,Epsilon,e2winf
    real:: Omega0

    integer::iFrequency

    !Nu, as given in Eq.(17)
    real::vNu

    character(LEN=*),parameter::NameSub = 'Cranmer_heating_function'
    !-------------------------------

    !Inputs conversion:
    RCGS      = RSI    * 1.0e+2
    UMagCGS   = UMagSI * 1.0e+2
    RhoCGS    = RhoSI  * 1.0e-3
    BMagCGS   = BMagSI * 1.0e+4

    !--------------------------------checks-------------------------
    !C  Check that required inputs are valid.
    !
    !      if ((r.lt.(0.99*xRsun)).or.(r.gt.(4.d+04*xRsun))) then
    !        write(0,*) ' '
    !        write(0,*) ' radius out of bounds! '
    !        write(0,*) ' '
    !        stop
    !      end if
    !
    if( RCGS < 0.990 * rSunCGS .or. RCGS > 4.0e+4 * rSunCGS)&
         call stop_mpi(NameSub//&
         ': the input heliocentric distance is out of bounds')

    !      if (uwind.gt.5.0d+08) then
    !        write(0,*) ' '
    !        write(0,*) ' wind speed out of bounds! '
    !        write(0,*) ' '
    !        stop
    !      end if

    if( UMagCGS > 5.0e+8)&
         call stop_mpi(NameSub//&
         ': wind speed is out of bounds')


    !
    !      if ((rho.lt.1.0d-29).or.(rho.gt.1.0d-04)) then
    !        write(0,*) ' '
    !        write(0,*) ' density out of bounds! '
    !        write(0,*) ' '
    !        stop
    !      end if

    if(RhoCGS < 1.0e-29 .or. RhoCGS > 1.0e-04)&
         call stop_mpi(NameSub//&
         ': density is out of bounds')
    !
    !      if ((Bmag.lt.1.0d-11).or.(Bmag.gt.1.0d+04)) then
    !        write(0,*) ' '
    !        write(0,*) ' magnetic field strength out of bounds! '
    !        write(0,*) ' '
    !        stop
    !      end if

    if(BMagCGS < 1.0e-11 .or. BMagCGS > 1.0e+4)&
         call stop_mpi(NameSub//&
         ': magnetic field strength is out of bounds')



    !C  Compute some basic plasma properties at this location.
    !
    !      Valf    = Bmag / sqrt(fourpi*rho)
    !      upVA    = uwind + Valf
    !      umVA    = uwind - Valf

    !      xnu     = Valf * r/(r+xRsun)/(r-xRsun)
    !
    VAlfven = BMagCGS / sqrt(cPi4 * RhoCGS)

    UPlusVA  = UMagCGS + VAlfven
    UMinusVA = UMagCGS - VAlfven

    vNu = VAlfven * RCGS/( RCGS**2 - rSunCGS2)


    !
    !C  If any of the optional inputs aren't given, estimate them. 
    !   if (vperp0.le.0.0) then
    !     UA     = FoBconst*Valf*Bmag/upVA/upVA
    !     vperp0 = sqrt(UA/rho)
    !   end if

    if(.not.present(WaveEnergyDensitySI))then
       WaveEnergyDensityCGS = cEnergyFluxPerB * VAlfven * BMagCGS/&
            UPlusVA**2
    else
       WaveEnergyDensityCGS = WaveEnergyDensitySI * 10.0
    end if
    vPerp0 = sqrt(WaveEnergyDensityCGS / RhoCGS)

    !---------------reflection coefficient--------------------------
    !
    !C  Two things will signal "perfect" wave reflection: (1) a wind speed
    !C  that is zero or negative (which implies a closed-field region), and
    !C  (2) chromospheric plasma conditions (i.e., a high enough density),
    !C  which implies that we're below the sharp transition-region boundary).
    !
    !      isetRcavg = 0
    !      if ((uwind.le.0.0).or.(rho.ge.2.0d-13)) then
    !        uwind     = abs(uwind)
    !        Rcavg     = 1.0
    !        isetRcavg = 1
    !      end if
    !


    IsSetReflectionCoef = UMagCGS <=0.0 .or. RhoCGS >=2.0e-13
    if(present(IsFullReflection))IsSetReflectionCoef=IsSetReflectionCoef&
         .or.IsFullReflection
    if(IsSetReflectionCoef) then 
       ReflectionCoef = 1.0
    else

       !C  If we haven't already set the reflection coefficient (for the
       !C  special cases listed above), then start the process of computing it.
       !C  Variables with "zero" in their names are computed in the limit of
       !C  zero frequency.  Variables with "inf" in their names are computed
       !C  in the high-frequency limit.  For arbitrary frequencies between
       !C  these two limiting cases, the coefficients are "bridged together"
       !C  below with the "epsilon" exponent.
       !
       !      rcrit0 = rcrit
       !      Vcrit0 = Vcrit
       !      vperp0 = vperp
       !      rhocrit = rho * ((uwind/Valf)**2)
       !      if (Vcrit0.le.0.0) then
       !        alpha   = 1. / (1. + 0.3*((rho/rhocrit)**0.25))
       !        Vcrit0  = (abs(uwind)**alpha) * (Valf**(1.-alpha))
       !      end if
       !
       !      if (rcrit0.le.0.0) rcrit0 = r*Valf/Vcrit0
       !      if (rcrit0.le.(1.001*xRsun)) rcrit0 = 1.001*xRsun
       !

       !Eq.(21):
       RhoCrit = RhoCGS * (UMagCGS / VAlfven)**2

       !Eq.(23):
       Alpha = 1.0/(1.0 + 0.30 * (RhoCGS / RhoCrit)**0.250)

       !Eq.(22):
       vCrit = Abs(UMagCGS)**Alpha * VAlfven**(1.0-Alpha)

       !Eq.(24)
       rCrit = max(1.001*rSunCGS, RCGS * VAlfven/vCrit)

       !C  Compute other necessary wave properties at this location.
       !
       !      xnucrit = Vcrit0 * rcrit0/(rcrit0+xRsun)/(rcrit0-xRsun)
       !      omega0  = 0.5 * (xnu + xnucrit)

       vNuCrit = vCrit * rCrit/( rCrit**2 - rSunCGS2 )
       Omega0  = 0.5 * (vNu + vNuCrit)


       !      psifac  = upVA/umVA
       !
       !       !
       !      if (isetRcavg.eq.0) then
       !
       !        Wzero   = 0.5*log(Valf/Vcrit0)
       !        e2wzero = exp(2.*Wzero)
       !        psizero = (1. - e2wzero)/(1. + e2wzero)
       !        Rzero   = psifac * psizero
       !
       if(abs(UMinusVA)< 0.0010 * VAlfven)then
          ReflectionCoef0 = Alpha

          ReflectionCoef2 = 0.0
          do iFrequency = 1,nFrequency
             ReflectionCoefInf = Alpha * RCGS**2/(RCGS**2 - rSunCGS2) * &
                  vNu/sqrt(Omega_I(iFrequency)**2 + vNu**2)

             Epsilon = 1.0/(1.0 + (Omega0/Omega_I(iFrequency))**2)
             ReflectionCoef = ReflectionCoefInf**Epsilon * ReflectionCoef0**(1.0 - Epsilon)

             ReflectionCoef2 = ReflectionCoef2 + Weight_I(iFrequency) *&
                  ReflectionCoef**2
          end do
       else
          !Apply eq.(13)
          ReflectionCoef0 = (UPlusVA/UMinusVA) * &
               (vCrit - VAlfven)/(vCrit + VAlfven)


          !C  Loop over the frequency spectrum in order to compute the reflection
          !C  coefficients in each frequency bin.  These are then convolved
          !C  together to form a spectrum-weighted reflection coefficient.
          !
          !        Rc2avg  = 0.0
          !
          !        do 1000 i=1,nfreq
          !         Winf = .25*(xnu-xnucrit)/sqrt(omega(i)*omega(i)+xnu*xnu)
          !         e2winf  = exp(2.*Winf)
          !         psiinf  = (1. - e2winf)/(1. + e2winf)
          !         Rinf    = psifac * psiinf
          !         epsilon = 1. / (1. + (omega0/omega(i))**2)
          !         Rcoef   = (Rinf**epsilon) * (Rzero**(1.-epsilon))
          !         Rc2avg  = Rc2avg + (Rcoef*Rcoef*weight(i))
          !1000    continue
          !
          !        Rcavg = sqrt(Rc2avg)
          !
          !      end if
          !

          ReflectionCoef2 = 0.0

          do iFrequency = 1, nFrequency
             e2winf = exp(0.50 * (vNu - vNuCrit)/sqrt(Omega_I(iFrequency)**2 + vNu**2) )
             ReflectionCoefInf = (1.0 - e2winf)/(1.0 + e2winf)*UPlusVA/UMinusVA

             !Eq.(20):
             Epsilon = 1.0/(1.0 + (Omega0/Omega_I(iFrequency))**2)

             !Eq.(25):
             ReflectionCoef = ReflectionCoefInf**Epsilon * ReflectionCoef0**(1.0 - Epsilon)

             ReflectionCoef2 = ReflectionCoef2 + Weight_I(iFrequency) *&
                  ReflectionCoef**2
          end do
       end if

       ReflectionCoef = sqrt(ReflectionCoef2)
    end if
    !--------------------------heating function---------------------------
    !C  Put together the pieces into the final expression for the turbulent
    !C  dissipation rate.
    !
    !      Zminus  = 2.*vperp0 / sqrt(1. + (Rcavg*Rcavg))
    !      Zplus   = Rcavg * Zminus
    !      Zcubed  = Zminus * Zplus * (Zminus + Zplus)
    !      ellperp = ellconst / sqrt(Bmag)
    !
    !      effic   = 1.
    !      if (iuseeff.eq.1) then
    !        tref  = 1./xnu
    !        teddy = ellperp*sq3pi/(1. + uwind/Valf)/vperp0
    !        effic = 1. / (1. + (teddy/tref))
    !      end if
    !
    !      Qheat   = 0.25 * rho * effic * (Zcubed/ellperp)
    !      gamma   = Qheat / (2.*rho*vperp0*vperp0)
    !
    !      return
    !      end


    ZMinus = 2.0 * VPerp0/sqrt(1.0 + ReflectionCoef**2)
    ZPlus  = ReflectionCoef * ZMinus

    CorrelationLength = cCorrelationLength/sqrt(BMagCGS)

    QHeatCGS = 0.250 * RhoCGS * &
         ( (ZPlus**2 * ZMinus + ZMinus**2 * ZPlus)/ CorrelationLength)

    if(UseEfficiency)then
       tEddy      = CorrelationLength * cSqrt3Pi/( VPerp0 * (1.0 + UMagCGS/VAlfven) )
       tReference = 1.0/vNu

       QHeatCGS = QHeatCGS /(1.0 +  tEddy / tReference )
    end if

    if(present(GammaSI))GammaSI = QHeatCGS/(2.0 * RhoCGS * VPerp0**2)

    QHeatSI = 0.10 * QHeatCGS  !0.1= 1e-7 J/(1e-2 m)**3

  end subroutine Cranmer_heating_function

end module ModNonWKBHeating
