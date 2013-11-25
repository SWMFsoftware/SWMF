!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!This code is a copyright protected software (c) 2002- University of Michigan
!==============================================================================
!These parameters may be used to parameterize
!the coronal heating in terms of waves
module ModAlfvenWaveHeating
  implicit none

  !If this locical is set to .true. the heating function is parameterized
  !in terms of the absorption of the Alfven Waves
  logical:: UseAlfvenWaveHeating = .false.

  !-------------------Original Steve Cranmer's code-----------------------!
  !C  Wave action conservation constant.  Earlier models found values that!
  !C  ranged between:                                                     !
  !C     parameter (FoBconst=2.0d+04)                                     !
  !C     parameter (FoBconst=9.0d+04)                                     !
  !C  but this paper uses an intermediate value:                          !
  !      parameter (FoBconst=5.0d+04)                                     !
  !-----------------------------------------------------------------------!

  !------------------Comment from Igor------------------------------------!
  !The consideration of the wave action conservation (about the wave      !
  !action see for example, Sokolov et al 2009) explains only the          !
  !conservation of this invariant along the flux tube. The choice of the  !
  !ratio, (Flux Density Along The Flux tube) / (Magnetic field intensity),!
  !to be also constant over the solar surface is formulated in Suzuki,2006!
  !(see also the papers cited there) as the contition for the average of  !
  ! < \delta B_\perp \delta u_\perp >. Another reincarnation of the same  !

  real :: cEnergyFluxPerBCgs = 5.0e4 !erg/cm2/s/Gs

  !If the wave energy flux, directed along the magnetic field and equal,  !
  !in neglecting the plasma speed, to VAlfven * (Wave Energy Density      !
  !is assumed to be proportional to the magnetic field intensity, then    !
  ! Wave Energy Density [CGS] = cEnergyFluxPerBCgs * B[CGS]/VAlfvenCgs =  !
  != cEnergyFluxPerBCgs * sqrt(4\pi * RhoCgs). The "adiabatic law" comes  !
  !from here: Ew=const *\sqrt(Rho).
  real :: cAdiabaticLaw = 0.0
contains
  !===================================
  subroutine set_adiabatic_law_4_waves
    use ModPhysics, ONLY: UnitRho_, UnitEnergyDens_, Si2No_V, No2Si_V
    use ModConst  , ONLY: cTwoPi

    real, parameter:: Cgs2Si4EnergyDens = 1.0e-7 & !1 erg to J
         / 1.0e-6   !1 cm3 to m3
    real, parameter:: Si2Cgs4Dens       = 1.0e+3 & !1 kg in g
         / 1.0e+6   !1 m3 in cm3
    !---------------------------------------

    cAdiabaticLaw = & ! (Wave Energy Demsity [NoDim]/sqrt(Rho[NoDim)=    1
         ( sqrt(2.0 * cTwoPi * No2Si_V(UnitRho_) * & ! Rho[NoDim]        2
         Si2Cgs4Dens)* &                      !                   3
         cEnergyFluxPerBCgs ) &      !Wave Energy Density [CGS]          4
         *Cgs2Si4EnergyDens   &      !Wave Energy Density [SI]           5
         *Si2No_V(UnitEnergyDens_)   !Wave Energy Density [NoDim]        6
    !    /sqrt(Rho[NoDim]), this square root cancels with that in line 2
  end subroutine set_adiabatic_law_4_waves
  !===================================  
  subroutine adiabatic_law_4_wave_state(State_V, Xyz_D, B0_D)

    use ModVarIndexes, ONLY: nVar, Bx_, Bz_, Ew_, Rho_
    use ModMain, ONLY: UseB0
    use ModWaves

    !Input and output parameters:   

    !WaveFirst_:WaveLast_ components of this vector are to be filled in:
    real, intent(inout) :: State_V(nVar)  

    !If UseAlfvenWaves, the Plus or Minus waves are intialized, depending on
    !the sign of {\bf B}\cdot{\bf r}, therefore, we need the following 
    !parameters:
    real, intent(in):: Xyz_D(3), B0_D(3)


    real:: BTotal_D(3), EWaveTotal
    !--------------------------------------------------------!

    EWaveTotal = cAdiabaticLaw * sqrt(State_V(Rho_))


    BTotal_D = State_V(Bx_:Bz_)
    if(UseB0) BTotal_D = BTotal_D + B0_D

    !Figure out the sign of {\bf B}\cdot{\bf r}
    if( sum( BTotal_D*Xyz_D ) > 0) then

       State_V(WaveFirst_:WaveLast_) = EWaveTotal * SpectrumPlus_W

    else

       State_V(WaveFirst_:WaveLast_) = EWaveTotal * SpectrumMinus_W

    end if
    if( UseWavePressureLtd )&
         State_V(Ew_) = sum(State_V(WaveFirst_:WaveLast_))
  end subroutine adiabatic_law_4_wave_state

end module ModAlfvenWaveHeating


!======================!Cranmer heating function!====================!
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
!======================!NOW used with WAVES !==========================!
!While using with ALFVEN WAVES only the REFLECTION COEFFICIENT is taken!
!from the Cranmer heting function. The boundary condition for waves    !
!should be taken in accordance with the the alfven heat flux.          !
!\                                                                     !
! Igor Sokolov, 1/26/10                                                !
!/                                                                     !
!======================================================================! 

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
  use ModAlfvenWaveHeating
  use ModPhysics,ONLY: UseStar,RadiusStar
  implicit none
  SAVE
  PRIVATE !Except

  !PUBLIC MEMBER:
  logical,public::UseCranmerHeating = .false.
  public:: cranmer_heating_function

  !The Cranmers heating function.

  !The usage:

  !Use ModNonWKBHeating

  !....
  !call cranmer_heating_function(&
  !      RSI,   &       !Heliocentric distance [m]
  !      RhoSI, &       !Density, [kg/m^3]
  !      UMagSI,&       !The magnitude of the solar wind speed, [m/s]
  !      BMagSI,&       !The magnitude of the magnetic field, [T]
  !      QHeatSI,&      !The calculated heating function,[J/m3/s]
  !      WaveEnergyDensitySI,&   !Optional input, is used to calculate VPerp
  !      IsFullReflection,   &   !Optional input, sets the reflection coef to 1
  !      GammaSI)       !Optional output, nonlinear damping rate in s^{-1}

  real,parameter:: cPi4 = 4.0 * cPi, cSqrt3Pi = 3.06998016655


  !C  Turbulent correlation length normalization.  Earlier models found
  !C  values that ranged from a minimum of (Cranmer et al. 2007):
  !C     parameter (ellconst=2.876d+08)
  !C  to a maximum of (Cranmer and van Ballegooijen 2005):
  !C     parameter (ellconst=1.155d+09)
  !C  and this paper uses an intermediate value:
  !      parameter (ellconst=6.000d+08)
  !
  !Typo in comments is fixed (igorsok Dec.12,2012: dimensions of cCorreletionlength
  !should be cm*sqrt{Gs}, not "divided by"
  real,parameter:: cCorrelationLength = 6.0e+8   !cm*sqrt{Gs}


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
  !\
  ! Steve Cranmer's comments, explaining the meaning of parameters to be used 
  !/ 
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

  subroutine cranmer_heating_function(&
       RSI, RhoSI, UMagSI, BMagSI,&            !INPUT
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
    real:: RCGS, UMagCGS, RhoCGS, BMagCGS, QHeatCGS, rSunCGS, rSunCGS2   !Converted variables

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

    character(LEN=*),parameter::NameSub = 'cranmer_heating_function'
    !-------------------------------

    if(UseStar)then
       rSunCGS = rSun*RadiusStar * 1.0e2 !Stellar radius in cm
    else
       rSunCGS = rSun * 1.0e2 !Solar radius in cm
    endif
    rSunCGS2= rSunCGS**2

    !Inputs conversion:
    RCGS      = max(RSI    * 1.0e+2, 1.001 * rSunCGS) 
    !(max is to limit the denominator in the expression for rCrit)

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
       WaveEnergyDensityCGS = cEnergyFluxPerBCgs * VAlfven * BMagCGS/&
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

  end subroutine cranmer_heating_function

end module ModNonWKBHeating
!==============================================================================
module ModUnsignedFluxModel

  ! based on Abett (2007)

  implicit none
  save

  private ! except

  ! Public methods
  public :: get_coronal_heat_factor
  public :: get_coronal_heating

  !Bill Abbet model, if .true.
  logical, public :: UseUnsignedFluxModel = .false.

  ! Normalized value of Heating constant
  real, public :: HeatFactor = 0.0

  ! Cgs value of total power input from coronal heating
  real, public :: TotalCoronalHeatingCgs = 1.0e+28

  ! Exponential Scale height to truncate heating function
  real, public :: DecayLength = 1.0
  real, public :: DtUpdateFlux = -1.0
  real, public :: UnsignedFluxHeight = -99999

contains

  !============================================================================

  subroutine get_coronal_heat_factor

    use ModAdvance,     ONLY: State_VGB, Bz_
    use ModGeometry,    ONLY: true_BLK, true_cell
    use ModMagnetogram, ONLY: nTheta, nPhi, dSinTheta, dPhi, &
         get_magnetogram_field
    use ModMain,        ONLY: nI, nJ, nK, nBlock, Unused_B, Time_Simulation,z_
    use ModMpi,         ONLY: MPI_REAL, MPI_SUM
    use ModNumConst,    ONLY: cHalfPi
    use ModPhysics,     ONLY: Si2No_V, No2Si_V, UnitX_, UnitT_, &
         UnitEnergyDens_, rBody
    use ModProcMH,      ONLY: nProc, iComm
    use BATL_lib,       ONLY: CellFace_DB, CellVolume_GB

    integer :: i, j, k, iBlock
    integer :: iTheta, iPhi, iError
    real :: UnsignedFluxCgs, dAreaCgs
    real :: HeatFunction, HeatFunctionVolume, HeatFunctionVolumePe
    real :: x, y, z, Theta, Phi, SinTheta, CosTheta, SinPhi, CosPhi
    real :: B0_D(3), BrSi, BrCgs, SumUnsignedBrCgs
    real :: BzCgs(1:nI,1:nJ), SumUnsignedBzCgs, UnsignedFluxCgsPe
    real    :: TotalCoronalHeating = -1.0, TimeUpdateLast = -1.0
    logical :: DoFirst = .true.

    real, parameter :: HeatExponent = 1.1488, HeatCoef = 89.4
    !--------------------------------------------------------------------------

    if(DoFirst .and. DtUpdateFlux <= 0.0)then

       ! uniform cell area on sphere
       dAreaCgs = rBody**2*dSinTheta*dPhi*No2Si_V(UnitX_)**2*1e4
       SumUnsignedBrCgs = 0.0

       do iTheta = 0, nTheta

          Theta = cHalfPi - asin((real(iTheta) + 0.5)*dSinTheta - 1.0)
          SinTheta = sin(Theta)
          CosTheta = cos(Theta)
          do iPhi = 1, nPhi
             Phi=(real(iPhi)-0.5)*dPhi
             SinPhi = sin(Phi)
             CosPhi = cos(Phi)

             x = rBody*SinTheta*CosPhi
             y = rBody*SinTheta*SinPhi
             z = rBody*CosTheta

             call get_magnetogram_field(x, y, z, B0_D)
             BrSi = (x*B0_D(1) + y*B0_D(2) + z*B0_D(3))/rBody
             BrCgs = BrSi*1e4
             SumUnsignedBrCgs = SumUnsignedBrCgs + abs(BrCgs)

          end do
       end do

       UnsignedFluxCgs = SumUnsignedBrCgs*dAreaCgs

       TotalCoronalHeatingCgs = HeatCoef*UnsignedFluxCgs**HeatExponent

       TotalCoronalHeating = TotalCoronalHeatingCgs*1e-7 &
            *Si2No_V(UnitEnergyDens_)*Si2No_V(UnitX_)**3/Si2No_V(UnitT_)

       DoFirst = .false.

    elseif( DtUpdateFlux > 0.0 .and. &
         Time_Simulation - TimeUpdateLast > DtUpdateFlux ) then
       UnsignedFluxCgs = 0.0
       do iBlock = 1, nBlock 
          if(Unused_B(iBlock)) cycle
          if(true_BLK(iBlock)) then
             dAreaCgs = CellFace_DB(z_,iBlock)*No2Si_V(UnitX_)**2*1e4

             call get_photosphere_field(iBlock, UnsignedFluxHeight, &
                  State_VGB(Bz_,1:nI,1:nJ,0:nK+1,iBlock), BzCgs)

             SumUnsignedBzCgs = sum(abs(BzCgs))
             UnsignedFluxCgs = UnsignedFluxCgs +  SumUnsignedBzCgs*dAreaCgs
          end if
       end do
       if(nProc>1)then
          UnsignedFluxCgsPe = UnsignedFluxCgs
          call MPI_allreduce(UnsignedFluxCgsPe, UnsignedFluxCgs, 1, &
               MPI_REAL, MPI_SUM, iComm, iError)
       end if
       TotalCoronalHeatingCgs = HeatCoef*UnsignedFluxCgs**HeatExponent

       TotalCoronalHeating = TotalCoronalHeatingCgs*1e-7 &
            *Si2No_V(UnitEnergyDens_)*Si2No_V(UnitX_)**3/Si2No_V(UnitT_)
       TimeUpdateLast = Time_Simulation

    end if

    HeatFunctionVolume = 0
    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE

       if(true_BLK(iBlock)) then
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             call get_heat_function(i, j, k, iBlock, HeatFunction)
             HeatFunctionVolume = HeatFunctionVolume &
                  + HeatFunction*CellVolume_GB(i,j,k,iBlock)
          end do; end do; end do
       else
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             if(true_cell(i,j,k,iBlock))then
                call get_heat_function(i, j, k, iBlock, HeatFunction)
                HeatFunctionVolume = HeatFunctionVolume &
                     + HeatFunction*CellVolume_GB(i,j,k,iBlock)
             end if
          end do; end do; end do
       end if
    end do

    if(nProc>1)then
       HeatFunctionVolumePe = HeatFunctionVolume
       call MPI_allreduce(HeatFunctionVolumePe, HeatFunctionVolume, 1, &
            MPI_REAL, MPI_SUM, iComm, iError)
    end if

    HeatFactor = TotalCoronalHeating/HeatFunctionVolume

  end subroutine get_coronal_heat_factor

  !============================================================================

  subroutine get_coronal_heating(i, j, k, iBlock, CoronalHeating)

    integer, intent(in) :: i, j, k, iBlock
    real, intent(out) :: CoronalHeating

    real :: HeatFunction
    !--------------------------------------------------------------------------

    call get_heat_function(i, j, k, iBlock, HeatFunction)

    CoronalHeating = HeatFactor*HeatFunction

  end subroutine get_coronal_heating

  !============================================================================

  subroutine get_heat_function(i, j, k, iBlock, HeatFunction)

    use ModMain, ONLY: UseB0, z_
    use ModAdvance, ONLY: State_VGB, B0_DGB, Bx_, Bz_
    use ModGeometry, ONLY: r_BLK
    use BATL_lib, ONLY: Xyz_DGB

    integer, intent(in) :: i, j, k, iBlock
    real, intent(out) :: HeatFunction

    real :: Bmagnitude, B_D(3)
    !--------------------------------------------------------------------------

    if(UseB0) then
       B_D = B0_DGB(:,i,j,k,iBlock) + State_VGB(Bx_:Bz_,i,j,k,iBlock)
    else
       B_D = State_VGB(Bx_:Bz_,i,j,k,iBlock)
    end if

    Bmagnitude = sqrt(sum(B_D**2))

    if(DtUpdateFlux <= 0.0)then
       HeatFunction = Bmagnitude*exp(-(r_BLK(i,j,k,iBlock)-1.0)/DecayLength)
    else
       if(Xyz_DGB(z_,i,j,k,iBlock)<UnsignedFluxHeight)then
          HeatFunction = 0.0
       else
          HeatFunction = Bmagnitude
       end if
    end if

  end subroutine get_heat_function

  !===========================================================================
  subroutine get_photosphere_field(iBlock, z_cut, Bz_V, BzCgs)

    use ModMain,      ONLY: nI, nJ, nK, z_
    use ModInterpolate, ONLY: find_cell
    use ModPhysics,   ONLY: No2Si_V, UnitB_
    use BATL_lib,     ONLY: CoordMin_DB, CoordMax_DB, CellSize_DB

    integer, intent(in) :: iBlock
    real, intent(in)    :: z_cut, Bz_V(1:nI, 1:nJ, 0:nK+1)
    real, intent(out)   :: BzCgs(1:nI, 1:nJ)
    real :: MinZ, MaxZ, DxLeft, iZ
    integer :: iLeft
    !--------------------------------------------------------------------------

    MinZ = CoordMin_DB(z_,iBlock)
    MaxZ = CoordMax_DB(z_,iBlock)

    BzCgs = 0.0
    if((UnsignedFluxHeight > MaxZ).or.(UnsignedFluxHeight < MinZ))&
         return

    iZ = (UnsignedFluxHeight - MinZ)/CellSize_DB(z_,iBlock) + 0.5
    call find_cell(0, nK+1, iZ, iLeft, DxLeft)

    BzCgs = ((1.0 - DxLeft)*Bz_V(1:nI, 1:nJ, iLeft) + &
         DxLeft*Bz_V(1:nI, 1:nJ, iLeft+1))*No2Si_V(UnitB_)*1e4

  end subroutine get_photosphere_field
end module ModUnsignedFluxModel
!=========================================!Master module!======================
module ModCoronalHeating

  use ModUnsignedFluxModel
  use ModNonWKBHeating
  use ModMain,      ONLY: nBLK, nI, nJ, nK
  use ModReadParam, ONLY: lStringLine
  use ModVarIndexes,ONLY: WaveFirst_, WaveLast_

  implicit none
  SAVE

  logical,public:: UseCoronalHeating = .false.
  character(len=lStringLine) :: NameModel, TypeCoronalHeating

  ! Variables and parameters for various heating models
  !\
  ! Exponential:
  !/  

  ! quantitative parameters for exponential heating model

  real :: HeatingAmplitude, HeatingAmplitudeCgs = 6.07e-7
  real :: DecayLengthExp = 0.7  ! in Solar Radii units
  logical :: UseExponentialHeating = .false.

  ! parameters for high-B transition (in Gauss)
  ! idea is to grossly approx 'Active Region' heating values
  logical :: UseArComponent = .false.
  real :: ArHeatB0 = 30.0
  real :: DeltaArHeatB0 = 5.0
  real :: ArHeatFactorCgs = 4.03E-05  ! cgs energy density = [ergs cm-3 s-1]

  ! open closed heat is if you want to have different heating
  ! between open and closed field lines. The routines for the WSA
  ! model in the ExpansionFactors module essentially do this, so will
  ! need to call them
  logical :: DoOpenClosedHeat = .false.
  real :: WsaT0 = 3.50

  !\
  ! Abbett's model
  !/


  ! Normalization constant for Abbett Model
  real :: HeatNormalization = 1.0

  ! Alfven wave dissipation
  logical :: UseAlfvenWaveDissipation = .false.
  real    :: LperpTimesSqrtBSi = 7.5e4 ! m T^(1/2)
  real    :: LperpTimesSqrtB
  real    :: Crefl = 0.04

  ! Variables for incompressible turbulence
  logical :: UseTurbulentCascade = .false.
  logical :: UseWaveReflection = .true.
  logical :: UseScaledCorrelationLength = .true.

  logical :: IsNewBlockAlfven = .true.

  ! long scale height heating (Ch = Coronal Hole)
  logical :: DoChHeat = .false.
  real :: HeatChCgs = 5.0e-7
  real :: DecayLengthCh = 0.7

  !Arrays for the calculated heat function and dissipated wave energy
  real :: CoronalHeating_C(1:nI,1:nJ,1:nK)
  real :: WaveDissipation_VC(WaveFirst_:WaveLast_,1:nI,1:nJ,1:nK)

  ! Electron heating is a fraction of the total coronal heating
  real :: QeByQtotal = 0.4
  ! The fraction of heating for ion parallel pressure: 1/3 of Qp
  real :: QparByQtotal = 0.2

  ! Dimensionless parameter for stochastic heating
  real :: StochasticHeating = 0.17

  logical,private:: DoInit = .true. 
contains
  !==========================================================================
  subroutine read_corona_heating(NameCommand)

    use ModReadParam,  ONLY: read_var
    use ModVarIndexes, ONLY: Lperp_

    character(len=*), intent(in):: NameCommand
    !----------------------------------------------------------------------
    select case(NameCommand)
    case("#CORONALHEATING")
       call read_var('TypeCoronalHeating', TypeCoronalHeating)

       !Initialize logicals
       UseCoronalHeating = .true.
       UseUnsignedFluxModel = .false.
       UseCranmerHeating    = .false.
       UseExponentialHeating= .false.
       UseAlfvenWaveDissipation = .false.
       UseTurbulentCascade = .false.
       select case(TypeCoronalHeating)
       case('F','none')
          UseCoronalHeating = .false.
       case('exponential')
          UseExponentialHeating = .true.
          call read_var('DecayLengthExp', DecayLengthExp)
          call read_var('HeatingAmplitudeCgs', HeatingAmplitudeCgs)

       case('unsignedflux','Abbett')
          UseUnsignedFluxModel = .true.
          call read_var('DecayLength', DecayLength)
          call read_var('HeatNormalization', HeatNormalization)
       case('Cranmer','NonWKB')
          UseCranmerHeating     = .true.
       case('alfvenwavedissipation')
          UseAlfvenWaveDissipation = .true.
          call read_var('LperpTimesSqrtBSi', LperpTimesSqrtBSi)
          call read_var('Crefl', Crefl)
       case('turbulentcascade')
          UseAlfvenWaveDissipation = .true.
          UseTurbulentCascade = .true.
          call read_var('UseWaveReflection', UseWaveReflection)
          call read_var('LperpTimesSqrtBSi', LperpTimesSqrtBSi)
       case default
          call stop_mpi('Read_corona_heating: unknown TypeCoronalHeating = ' &
               // TypeCoronalHeating)
       end select
    case("#ACTIVEREGIONHEATING")
       call read_var('UseArComponent', UseArComponent)
       if(UseArComponent) then
          call read_var('ArHeatFactorCgs', ArHeatFactorCgs)
          call read_var('ArHeatB0', ArHeatB0)
          call read_var('DeltaArHeatB0', DeltaArHeatB0)
       endif
    case("#LONGSCALEHEATING")
       call read_var('DoChHeat', DoChHeat)
       if(DoChHeat)then
          call read_var('HeatChCgs', HeatChCgs)
          call read_var('DecayLengthCh', DecayLengthCh)
       end if
    case("#ELECTRONHEATING")
       call read_var('QeByQtotal', QeByQtotal)
    case("#ANISOIONHEATING")
       call read_var('QparByQtotal', QparByQtotal)
    case("#HEATPARTITIONING")
       call read_var('StochasticHeating', StochasticHeating)
    case default
       call stop_mpi('Read_corona_heating: unknown command = ' &
            // NameCommand)
    end select

  end subroutine read_corona_heating
  !============================================================================
  subroutine init_coronal_heating
    use ModPhysics, ONLY: Si2No_V, UnitEnergyDens_, UnitT_, UnitB_, UnitX_
    use ModAlfvenWaveHeating, ONLY: set_adiabatic_law_4_waves

    !--------------------------------------------------------------------------
    if(.not.DoInit) RETURN
    DoInit = .false.

    if(UseExponentialHeating)then
       HeatingAmplitude =  HeatingAmplitudeCgs*0.1 &
            *Si2No_V(UnitEnergyDens_)/Si2No_V(UnitT_)
    end if

    if(UseAlfvenWaveDissipation .and. UseScaledCorrelationLength) &
         LperpTimesSqrtB = LperpTimesSqrtBSi &
         *Si2No_V(UnitX_)*sqrt(Si2No_V(UnitB_))

    call set_adiabatic_law_4_waves

  end subroutine init_coronal_heating
  !=========================================================================
  subroutine get_cell_heating(i, j, k, iBlock, CoronalHeating)

    use ModGeometry,       ONLY: r_BLK
    use ModPhysics,        ONLY: Si2No_V, No2Si_V, UnitEnergyDens_, UnitT_, &
         No2Io_V, UnitB_,UnitRho_,UnitX_,UnitU_
    use ModExpansionFactors, ONLY: UMin
    use ModMain,       ONLY: x_, z_, UseB0
    use ModVarIndexes, ONLY: Bx_, Bz_,Rho_,RhoUx_,RhoUz_
    use ModAdvance,    ONLY: State_VGB, B0_DGB
    use BATL_lib,      ONLY: Xyz_DGB

    integer, intent(in) :: i, j, k, iBlock
    real, intent(out) :: CoronalHeating

    real :: HeatCh

    ! parameters for open/closed. This uses WSA model to determine if 
    ! cell is in an 'open' or 'closed' field region.
    !
    ! ** NOTE ** WSA does field line tracing on an auxiliary grid.
    ! should really be using the computational domain, but global 
    ! feild line tracing for this purpose is not easily implemented
    real :: Ufinal
    real :: UminIfOpen = 290.0
    real :: UmaxIfOpen = 300.0
    real :: Weight
    real :: B_D(3)

    ! local variables for ArHeating (Active Region Heating)
    real :: FractionB, Bcell

    real :: RhoSI, RSI, UMagSI, BMagSI, QHeatSI
    real :: WaveDissipation_V(WaveFirst_:WaveLast_)
    !--------------------------------------------------------------------------

    if(UseB0)then
       B_D = State_VGB(Bx_:Bz_,i,j,k,iBlock) + B0_DGB(x_:z_,i,j,k,iBlock)
    else
       B_D = State_VGB(Bx_:Bz_,i,j,k,iBlock)
    end if

    if(UseAlfvenWaveDissipation)then

       if(UseTurbulentCascade)then
          call turbulent_cascade(i, j, k, iBlock, WaveDissipation_V, &
               CoronalHeating)
       else
          call calc_alfven_wave_dissipation(i, j, k, iBlock, &
               WaveDissipation_V, CoronalHeating)
       end if

    elseif(UseUnsignedFluxModel)then

       call get_coronal_heating(i, j, k, iBlock, CoronalHeating)
       CoronalHeating = CoronalHeating * HeatNormalization

    elseif(UseExponentialHeating)then

       CoronalHeating = HeatingAmplitude &
            *exp(- max(r_BLK(i,j,k,iBlock) - 1.0, 0.0) / DecayLengthExp)

    elseif(UseCranmerHeating)then

       RSI   = r_BLK(i,j,k,iBlock) * No2Si_V(UnitX_)
       RhoSI = State_VGB(Rho_,i,j,k,iBlock) * No2Si_V(UnitRho_)
       UMagSI= sqrt( sum( State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock)**2 ) )/&
            State_VGB(Rho_,i,j,k,iBlock) * No2Si_V(UnitU_)
       BmagSI= sqrt( sum( B_D**2 ) ) * No2Si_V(UnitB_)
       if(.not.DoOpenClosedHeat)then
          call cranmer_heating_function(&
               RSI=RSI,      &
               RhoSI = RhoSI,&
               UMagSI=UMagSI,&
               BMagSI=BMagSI,&
               QHeatSI=QHeatSI)
       else
          UminIfOpen = UMin*1.05
          call get_bernoulli_integral( &
               Xyz_DGB(:,i,j,k,iBlock)/R_BLK(i,j,k,iBlock), UFinal)
          call cranmer_heating_function(&
               RSI=RSI,      &
               RhoSI  = RhoSI,&
               UMagSI =UMagSI,&
               BMagSI =BMagSI,&
               QHeatSI=QHeatSI,&
               IsFullReflection = UFinal< UMinIfOpen)
       end if

       CoronalHeating = QHeatSI *  Si2No_V(UnitEnergyDens_)/Si2No_V(UnitT_)
    else
       CoronalHeating = 0.0
    end if

    if(DoOpenClosedHeat.and.(.not.UseCranmerHeating))then
       ! If field is less than 1.05 times the minimum speed, mark as closed
       ! Interpolate between 1.05 and 1.10 for smoothness
       UminIfOpen = UMin*1.05
       UmaxIfOpen = UMin*1.1
       call get_bernoulli_integral( &
            Xyz_DGB(:,i,j,k,iBlock)/R_BLK(i,j,k,iBlock), UFinal)
       if(UFinal <= UminIfOpen) then
          Weight = 0.0
       else if (UFinal >= UmaxIfOpen) then
          Weight = 1.0
       else 
          Weight = (UFinal - UminIfOpen)/(UmaxIfOpen - UminIfOpen)
       end if

       CoronalHeating = (1.0 - Weight) * CoronalHeating
    end if

    if(DoChHeat) then
       HeatCh = HeatChCgs * 0.1 * Si2No_V(UnitEnergyDens_)/Si2No_V(UnitT_)
       CoronalHeating = CoronalHeating + HeatCh &
            *exp(- max(r_BLK(i,j,k,iBlock) - 1.0, 0.0) / DecayLengthCh)
    end if

    if(UseExponentialHeating.and.UseArComponent) then

       Bcell = No2Io_V(UnitB_) * sqrt( sum( B_D**2 ) )

       FractionB = 0.5*(1.0+tanh((Bcell - ArHeatB0)/DeltaArHeatB0))
       CoronalHeating = max(CoronalHeating, & 
            FractionB * ArHeatFactorCgs * Bcell &
            * 0.1 * Si2No_V(UnitEnergyDens_)/Si2No_V(UnitT_))

    endif

  end subroutine get_cell_heating
  !============================================================================

  subroutine get_block_heating(iBlock)

    use ModGeometry,       ONLY: r_BLK
    use ModPhysics,        ONLY: Si2No_V, No2Si_V, UnitEnergyDens_, UnitT_, &
         No2Io_V, UnitB_,UnitRho_,UnitX_,UnitU_
    use ModExpansionFactors, ONLY: UMin
    use ModMain,       ONLY: x_, z_, UseB0
    use ModVarIndexes, ONLY: Bx_, Bz_,Rho_,RhoUx_,RhoUz_
    use ModAdvance,    ONLY: State_VGB, B0_DGB
    use BATL_lib,      ONLY: Xyz_DGB

    integer, intent(in) :: iBlock

    integer             :: i, j, k
    real :: HeatCh

    ! parameters for open/closed. This uses WSA model to determine if 
    ! cell is in an 'open' or 'closed' field region.
    !
    ! ** NOTE ** WSA does field line tracing on an auxiliary grid.
    ! should really be using the computational domain, but global 
    ! feild line tracing for this purpose is not easily implemented
    real :: Ufinal
    real :: UminIfOpen = 290.0
    real :: UmaxIfOpen = 300.0
    real :: Weight
    real :: B_D(3)

    ! local variables for ArHeating (Active Region Heating)
    real :: FractionB, Bcell

    real :: RhoSI, RSI, UMagSI, BMagSI, QHeatSI
    !--------------------------------------------------------------------------

    if(UseAlfvenWaveDissipation)then

       if(UseTurbulentCascade)then
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             call turbulent_cascade(i, j, k, iBlock, &
                  WaveDissipation_VC(:,i,j,k), CoronalHeating_C(i,j,k))
          end do; end do; end do
       else
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             call calc_alfven_wave_dissipation(i, j, k, iBlock, &
                  WaveDissipation_VC(:,i,j,k), CoronalHeating_C(i,j,k))
          end do; end do; end do
       end if

    elseif(UseUnsignedFluxModel)then

       do k=1,nK;do j=1,nJ; do i=1,nI

          call get_coronal_heating(i, j, k, iBlock, CoronalHeating_C(i,j,k))
          CoronalHeating_C(i,j,k) = CoronalHeating_C(i,j,k) * HeatNormalization

       end do; end do; end do

    elseif(UseExponentialHeating)then
       do k=1,nK;do j=1,nJ; do i=1,nI

          CoronalHeating_C(i,j,k) = HeatingAmplitude &
               *exp(- max(r_BLK(i,j,k,iBlock) - 1.0, 0.0) / DecayLengthExp)

       end do; end do; end do

    elseif(UseCranmerHeating)then
       UminIfOpen = UMin*1.05

       do k=1,nK;do j=1,nJ; do i=1,nI

          if(UseB0)then
             B_D = State_VGB(Bx_:Bz_,i,j,k,iBlock) + B0_DGB(x_:z_,i,j,k,iBlock)
          else
             B_D = State_VGB(Bx_:Bz_,i,j,k,iBlock)
          end if

          RSI   = r_BLK(i,j,k,iBlock) * No2Si_V(UnitX_)
          RhoSI = State_VGB(Rho_,i,j,k,iBlock) * No2Si_V(UnitRho_)
          UMagSI= sqrt( sum( State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock)**2 ) )/&
               State_VGB(Rho_,i,j,k,iBlock) * No2Si_V(UnitU_)
          BmagSI= sqrt( sum( B_D**2 ) ) * No2Si_V(UnitB_)
          if(.not.DoOpenClosedHeat)then
             call cranmer_heating_function(&
                  RSI=RSI,      &
                  RhoSI = RhoSI,&
                  UMagSI=UMagSI,&
                  BMagSI=BMagSI,&
                  QHeatSI=QHeatSI)
          else

             call get_bernoulli_integral( &
                  Xyz_DGB(:,i,j,k,iBlock)/R_BLK(i,j,k,iBlock), UFinal)
             call cranmer_heating_function(&
                  RSI=RSI,      &
                  RhoSI  = RhoSI,&
                  UMagSI =UMagSI,&
                  BMagSI =BMagSI,&
                  QHeatSI=QHeatSI,&
                  IsFullReflection = UFinal< UMinIfOpen)
          end if

          CoronalHeating_C(i,j,k) = QHeatSI *  Si2No_V(UnitEnergyDens_)/Si2No_V(UnitT_)

       end do; end do; end do
    else
       CoronalHeating_C = 0.0
    end if

    if(DoOpenClosedHeat.and.(.not.UseCranmerHeating))then
       ! If field is less than 1.05 times the minimum speed, mark as closed
       ! Interpolate between 1.05 and 1.10 for smoothness
       UminIfOpen = UMin*1.05
       UmaxIfOpen = UMin*1.1
       do k=1,nK; do j=1,nJ; do i=1,nI
          call get_bernoulli_integral( &
               Xyz_DGB(:,i,j,k,iBlock)/R_BLK(i,j,k,iBlock), UFinal)
          if(UFinal <= UminIfOpen) then
             Weight = 0.0
          else if (UFinal >= UmaxIfOpen) then
             Weight = 1.0
          else 
             Weight = (UFinal - UminIfOpen)/(UmaxIfOpen - UminIfOpen)
          end if

          CoronalHeating_C(i,j,k) = (1.0 - Weight) * CoronalHeating_C(i,j,k)
       end do; end do; end do
    end if

    if(DoChHeat) then
       HeatCh = HeatChCgs * 0.1 * Si2No_V(UnitEnergyDens_)/Si2No_V(UnitT_)
       do k=1,nK; do j=1,nJ; do i=1,nI
          CoronalHeating_C(i,j,k) = CoronalHeating_C(i,j,k) + HeatCh &
               *exp(- max(r_BLK(i,j,k,iBlock) - 1.0, 0.0) / DecayLengthCh)
       end do; end do; end do
    end if

    if(UseExponentialHeating.and.UseArComponent) then
       do k=1,nK; do j=1,nJ; do i=1,nI

          if(UseB0)then
             B_D = State_VGB(Bx_:Bz_,i,j,k,iBlock) + B0_DGB(x_:z_,i,j,k,iBlock)
          else
             B_D = State_VGB(Bx_:Bz_,i,j,k,iBlock)
          end if

          Bcell = No2Io_V(UnitB_) * sqrt( sum( B_D**2 ) )

          FractionB = 0.5*(1.0+tanh((Bcell - ArHeatB0)/DeltaArHeatB0))
          CoronalHeating_C(i,j,k) = max(CoronalHeating_C(i,j,k), & 
               FractionB * ArHeatFactorCgs * Bcell &
               * 0.1 * Si2No_V(UnitEnergyDens_)/Si2No_V(UnitT_))
       end do; end do; end do
    endif

  end subroutine get_block_heating
  !============================================================================

  subroutine calc_alfven_wave_dissipation(i, j, k, iBlock, WaveDissipation_V, &
       CoronalHeating)

    use ModAdvance, ONLY: State_VGB, B0_DGB
    use ModMain, ONLY: UseB0
    use ModVarIndexes, ONLY: Rho_, Bx_, Bz_

    integer, intent(in) :: i, j, k, iBlock
    real, intent(out)   :: WaveDissipation_V(WaveFirst_:WaveLast_), &
         CoronalHeating

    real :: EwavePlus, EwaveMinus, FullB_D(3), FullB, Coef
    !--------------------------------------------------------------------------

    if(UseB0)then
       FullB_D = B0_DGB(:,i,j,k,iBlock) + State_VGB(Bx_:Bz_,i,j,k,iBlock)
    else
       FullB_D = State_VGB(Bx_:Bz_,i,j,k,iBlock)
    end if
    FullB = sqrt(sum(FullB_D**2))

    Coef = sqrt(FullB/State_VGB(Rho_,i,j,k,iBlock))/LperpTimesSqrtB

    EwavePlus  = State_VGB(WaveFirst_,i,j,k,iBlock)
    EwaveMinus = State_VGB(WaveLast_,i,j,k,iBlock)

    WaveDissipation_V(WaveFirst_) = Coef*EwavePlus &
         *sqrt(max(EwaveMinus,Crefl**2*EwavePlus))

    WaveDissipation_V(WaveLast_) = Coef*EwaveMinus &
         *sqrt(max(EwavePlus,Crefl**2*EwaveMinus))

    CoronalHeating = sum(WaveDissipation_V)

  end subroutine calc_alfven_wave_dissipation

  !============================================================================

  subroutine turbulent_cascade(i, j, k, iBlock, WaveDissipation_V, &
       CoronalHeating)

    use BATL_size, ONLY: nDim
    use ModAdvance, ONLY: State_VGB
    use ModB0, ONLY: B0_DGB
    use ModMain, ONLY: UseB0
    use ModVarIndexes, ONLY: Rho_, Bx_, Bz_, Lperp_

    integer, intent(in) :: i, j, k, iBlock
    real, intent(out)   :: WaveDissipation_V(WaveFirst_:WaveLast_), &
         CoronalHeating

    real :: FullB_D(3), FullB, Rho, Coef
    real :: ReflectionRate, GradLogAlfven_D(nDim), CurlU_D(3), b_D(3)
    real :: EwavePlus, EwaveMinus

    character(len=*), parameter :: &
         NameSub = 'ModCoronalHeating::turbulent_cascade'
    !--------------------------------------------------------------------------

    ! Low-frequency cascade due to small-scale nonlinearities

    if(UseB0)then
       FullB_D = B0_DGB(:,i,j,k,iBlock) + State_VGB(Bx_:Bz_,i,j,k,iBlock)
    else
       FullB_D = State_VGB(Bx_:Bz_,i,j,k,iBlock)
    end if
    FullB = sqrt(sum(FullB_D**2))

    if(Lperp_ > 1 .and. .not.UseScaledCorrelationLength)then
       ! Note that Lperp is multiplied with the density
       Coef = 2.0*sqrt(State_VGB(Rho_,i,j,k,iBlock)) &
            /State_VGB(Lperp_,i,j,k,iBlock)
    else
       Coef = 2.0*sqrt(FullB/State_VGB(Rho_,i,j,k,iBlock))/LperpTimesSqrtB
    end if

    EwavePlus  = State_VGB(WaveFirst_,i,j,k,iBlock)
    EwaveMinus = State_VGB(WaveLast_,i,j,k,iBlock)

    if(UseWaveReflection)then
       WaveDissipation_V(WaveFirst_) = Coef*sqrt(EwaveMinus)*EwavePlus
       WaveDissipation_V(WaveLast_) = Coef*sqrt(EwavePlus)*EwaveMinus
    else
       b_D = FullB_D/max(1e-15, FullB)
       Rho = State_VGB(Rho_,i,j,k,iBlock)

       call get_grad_log_alfven_speed(i, j, k, iBlock, GradLogAlfven_D)

       call get_curl_u(i, j, k, iBlock, CurlU_D)

       ! Reflection rate driven by Alfven speed gradient and
       ! vorticity along the field lines
       ReflectionRate = sqrt( (sum(b_D*CurlU_D))**2 &
            + (sum(FullB_D(:nDim)*GradLogAlfven_D))**2/Rho )

       ! Local dissipation of Alfven waves.
       ! Correct at top of closed field lines, where turbulence is balanced
       WaveDissipation_V(WaveFirst_) = &
            max(ReflectionRate, Coef*sqrt(EwaveMinus))*EwavePlus
       WaveDissipation_V(WaveLast_) = &
            max(ReflectionRate, Coef*sqrt(EwavePlus))*EwaveMinus
    end if

    CoronalHeating = sum(WaveDissipation_V)

  end subroutine turbulent_cascade

  !============================================================================

  subroutine get_wave_reflection(iBlock)

    use BATL_size, ONLY: nDim, nI, nJ, nK
    use ModAdvance, ONLY: State_VGB, Source_VC
    use ModB0, ONLY: B0_DGB
    use ModChromosphere,  ONLY: DoExtendTransitionRegion, extension_factor, &
         get_tesi_c, TeSi_C
    use ModGeometry, ONLY: true_cell
    use ModMain, ONLY: UseB0
    use ModVarIndexes, ONLY: Rho_, Bx_, Bz_

    integer, intent(in) :: iBlock

    integer :: i, j, k
    real :: GradLogAlfven_D(nDim), CurlU_D(3), b_D(3)
    real :: FullB_D(3), FullB, Rho, DissipationRateMax, ReflectionRate
    real :: EwavePlus, EwaveMinus
    !--------------------------------------------------------------------------

    if(DoExtendTransitionRegion) call get_tesi_c(iBlock, TeSi_C)

    do k = 1, nK; do j = 1, nJ; do i = 1, nI
       if(.not.true_cell(i,j,k,iBlock)) CYCLE

       call get_grad_log_alfven_speed(i, j, k, iBlock, GradLogAlfven_D)

       call get_curl_u(i, j, k, iBlock, CurlU_D)

       if(UseB0)then
          FullB_D = B0_DGB(:,i,j,k,iBlock) + State_VGB(Bx_:Bz_,i,j,k,iBlock)
       else
          FullB_D = State_VGB(Bx_:Bz_,i,j,k,iBlock)
       end if
       FullB = sqrt(sum(FullB_D**2))
       b_D = FullB_D/max(1e-15, FullB)

       Rho = State_VGB(Rho_,i,j,k,iBlock)

       EwavePlus  = State_VGB(WaveFirst_,i,j,k,iBlock)
       EwaveMinus = State_VGB(WaveLast_,i,j,k,iBlock)

       DissipationRateMax = &
            2.0*sqrt(max(EwavePlus,EwaveMinus)*FullB/Rho)/LperpTimesSqrtB

       if(DoExtendTransitionRegion) DissipationRateMax = &
            DissipationRateMax/extension_factor(TeSi_C(i,j,k))

       ! Reflection rate driven by Alfven speed gradient and
       ! vorticity along the field lines
       ReflectionRate = sqrt( (sum(b_D*CurlU_D))**2 &
            + (sum(FullB_D(:nDim)*GradLogAlfven_D))**2/Rho )

       ! Clip the reflection rate from above with maximum dissipation rate
       ReflectionRate = min(ReflectionRate, DissipationRateMax)

       ! No reflection when turbulence is balanced (waves are then
       ! assumed to be uncorrelated)
       if(4.0*EwaveMinus <= EwavePlus)then
          ReflectionRate = ReflectionRate*(1.0-2.0*sqrt(EwaveMinus/EwavePlus))
       elseif(4.0*EwavePlus <= EwaveMinus)then
          ReflectionRate = ReflectionRate*(2.0*sqrt(EwavePlus/EwaveMinus)-1.0)
       else
          ReflectionRate = 0.0
       end if

       Source_VC(WaveFirst_,i,j,k) = Source_VC(WaveFirst_,i,j,k) &
            - ReflectionRate*sqrt(EwavePlus*EwaveMinus)
       Source_VC(WaveLast_,i,j,k) = Source_VC(WaveLast_,i,j,k) &
            + ReflectionRate*sqrt(EwavePlus*EwaveMinus)

    end do; end do; end do

  end subroutine get_wave_reflection

  !============================================================================

  subroutine get_grad_log_alfven_speed(i, j, k, iBlock, GradLogAlfven_D)

    use BATL_lib, ONLY: IsCartesianGrid, x_, y_, z_, &
         CellSize_DB, FaceNormal_DDFB, CellVolume_GB
    use BATL_size, ONLY: nDim, nI, j0_, nJp1_, k0_, nKp1_

    integer, intent(in) :: i, j, k, iBlock
    real, intent(out) :: GradLogAlfven_D(nDim)

    real, save :: LogAlfven_FD(0:nI+1,j0_:nJp1_,k0_:nKp1_,nDim)
    !--------------------------------------------------------------------------
    if(IsNewBlockAlfven)then
       call get_log_alfven_speed

       IsNewBlockAlfven = .false.
    end if

    if(IsCartesianGrid)then
       GradLogAlfven_D(x_) = 1.0/CellSize_DB(x_,iBlock) &
            *(LogAlfven_FD(i+1,j,k,x_) - LogAlfven_FD(i,j,k,x_))
       if(nJ > 1) GradLogAlfven_D(y_) = 1.0/CellSize_DB(y_,iBlock) &
            *(LogAlfven_FD(i,j+1,k,y_) - LogAlfven_FD(i,j,k,y_))
       if(nK > 1) GradLogAlfven_D(z_) = 1.0/CellSize_DB(z_,iBlock) &
            *(LogAlfven_FD(i,j,k+1,z_) - LogAlfven_FD(i,j,k,z_))
    else
       GradLogAlfven_D = &
            LogAlfven_FD(i+1,j,k,x_)*FaceNormal_DDFB(:,x_,i+1,j,k,iBlock) &
            - LogAlfven_FD(i,j,k,x_)*FaceNormal_DDFB(:,x_,i,j,k,iBlock)
       if(nJ > 1) GradLogAlfven_D = GradLogAlfven_D + &
            LogAlfven_FD(i,j+1,k,y_)*FaceNormal_DDFB(:,y_,i,j+1,k,iBlock) &
            - LogAlfven_FD(i,j,k,y_)*FaceNormal_DDFB(:,y_,i,j,k,iBlock)
       if(nK > 1) GradLogAlfven_D = GradLogAlfven_D + &
            LogAlfven_FD(i,j,k+1,z_)*FaceNormal_DDFB(:,z_,i,j,k+1,iBlock) &
            - LogAlfven_FD(i,j,k,z_)*FaceNormal_DDFB(:,z_,i,j,k,iBlock)

       GradLogAlfven_D = GradLogAlfven_D/CellVolume_GB(i,j,k,iBlock)
    end if

  contains

    subroutine get_log_alfven_speed

      use ModAdvance, ONLY: &
           LeftState_VX,  LeftState_VY,  LeftState_VZ,  &
           RightState_VX, RightState_VY, RightState_VZ
      use ModB0, ONLY: B0_DX, B0_DY, B0_DZ
      use ModMain, ONLY: UseB0
      use ModVarIndexes, ONLY: Rho_, Bx_, Bz_

      integer :: i, j, k
      real :: Rho, FullB_D(3)
      !------------------------------------------------------------------------
      do k = 1, nK; do j = 1, nJ; do i = 1, nI+1
         FullB_D = 0.5*(LeftState_VX(Bx_:Bz_,i,j,k) &
              + RightState_VX(Bx_:Bz_,i,j,k))
         if(UseB0) FullB_D = FullB_D + B0_DX(:,i,j,k)
         Rho = 0.5*(LeftState_VX(Rho_,i,j,k) + RightState_VX(Rho_,i,j,k))
         LogAlfven_FD(i,j,k,x_) = log(sqrt(max(sum(FullB_D**2), 1e-30)/Rho))
      end do; end do; end do

      if(nJ > 1)then
         do k = 1, nK; do j = 1, nJ+1; do i = 1, nI
            FullB_D = 0.5*(LeftState_VY(Bx_:Bz_,i,j,k) &
                 + RightState_VY(Bx_:Bz_,i,j,k))
            if(UseB0) FullB_D = FullB_D + B0_DY(:,i,j,k)
            Rho = 0.5*(LeftState_VY(Rho_,i,j,k) + RightState_VY(Rho_,i,j,k))
            LogAlfven_FD(i,j,k,y_) = log(sqrt(max(sum(FullB_D**2), 1e-30)/Rho))
         end do; end do; end do
      end if

      if(nK > 1)then
         do k = 1, nK+1; do j = 1, nJ; do i = 1, nI
            FullB_D = 0.5*(LeftState_VZ(Bx_:Bz_,i,j,k) &
                 + RightState_VZ(Bx_:Bz_,i,j,k))
            if(UseB0) FullB_D = FullB_D + B0_DZ(:,i,j,k)
            Rho = 0.5*(LeftState_VZ(Rho_,i,j,k) + RightState_VZ(Rho_,i,j,k))
            LogAlfven_FD(i,j,k,z_) = log(sqrt(max(sum(FullB_D**2), 1e-30)/Rho))
         end do; end do; end do
      end if

    end subroutine get_log_alfven_speed

  end subroutine get_grad_log_alfven_speed

  !============================================================================

  subroutine get_curl_u(i, j, k, iBlock, CurlU_D)

    use BATL_lib, ONLY: IsCartesianGrid, CellSize_DB, FaceNormal_DDFB, &
         CellVolume_GB, x_, y_, z_
    use ModAdvance, ONLY: &
         LeftState_VX,  LeftState_VY,  LeftState_VZ,  &
         RightState_VX, RightState_VY, RightState_VZ
    use ModCoordTransform, ONLY: cross_product
    use ModSize, ONLY: MaxDim
    use ModVarIndexes, ONLY: Ux_, Uy_, Uz_

    integer, intent(in) :: i, j, k, iBlock
    real, intent(out) :: CurlU_D(MaxDim)

    real :: DxInvHalf, DyInvHalf, DzInvHalf
    !--------------------------------------------------------------------------

    if(IsCartesianGrid)then
       DxInvHalf = 0.5/CellSize_DB(x_,iBlock)
       DyInvHalf = 0.5/CellSize_DB(y_,iBlock)
       DzInvHalf = 0.5/CellSize_DB(z_,iBlock)

       CurlU_D(x_) = &
            DyInvHalf*(LeftState_VY(Uz_,i,j+1,k) + RightState_VY(Uz_,i,j+1,k) &
            - LeftState_VY(Uz_,i,j,k) - RightState_VY(Uz_,i,j,k)) - &
            DzInvHalf*(LeftState_VZ(Uy_,i,j,k+1) + RightState_VZ(Uy_,i,j,k+1) &
            - LeftState_VZ(Uy_,i,j,k) - RightState_VZ(Uy_,i,j,k))

       CurlU_D(y_) = &
            DzInvHalf*(LeftState_VZ(Ux_,i,j,k+1) + RightState_VZ(Ux_,i,j,k+1) &
            - LeftState_VZ(Ux_,i,j,k) - RightState_VZ(Ux_,i,j,k)) - &
            DxInvHalf*(LeftState_VX(Uz_,i+1,j,k) + RightState_VX(Uz_,i+1,j,k) &
            - LeftState_VX(Uz_,i,j,k) - RightState_VX(Uz_,i,j,k))

       CurlU_D(z_) = &
            DxInvHalf*(LeftState_VX(Uy_,i+1,j,k) + RightState_VX(Uy_,i+1,j,k) &
            - LeftState_VX(Uy_,i,j,k) - RightState_VX(Uy_,i,j,k)) - &
            DyInvHalf*(LeftState_VY(Ux_,i,j+1,k) + RightState_VY(Ux_,i,j+1,k) &
            - LeftState_VY(Ux_,i,j,k) - RightState_VY(Ux_,i,j,k))
    else
       CurlU_D(:) = &
            + cross_product( FaceNormal_DDFB(:,1,i+1,j,k,iBlock),           &
            LeftState_VX(Ux_:Uz_,i+1,j,k) + RightState_VX(Ux_:Uz_,i+1,j,k)) &
            - cross_product( FaceNormal_DDFB(:,1,i  ,j,k,iBlock),           &
            LeftState_VX(Ux_:Uz_,i  ,j,k) + RightState_VX(Ux_:Uz_,i  ,j,k)) &
            + cross_product( FaceNormal_DDFB(:,2,i,j+1,k,iBlock),           &
            LeftState_VY(Ux_:Uz_,i,j+1,k) + RightState_VY(Ux_:Uz_,i,j+1,k)) &
            - cross_product( FaceNormal_DDFB(:,2,i,j  ,k,iBlock),           &
            LeftState_VY(Ux_:Uz_,i,j  ,k) + RightState_VY(Ux_:Uz_,i,j  ,k)) &
            + cross_product( FaceNormal_DDFB(:,3,i,j,k+1,iBlock),           &
            LeftState_VZ(Ux_:Uz_,i,j,k+1) + RightState_VZ(Ux_:Uz_,i,j,k+1)) &
            - cross_product( FaceNormal_DDFB(:,3,i,j,k  ,iBlock),           &
            LeftState_VZ(Ux_:Uz_,i,j,k  ) + RightState_VZ(Ux_:Uz_,i,j,k  ))

       CurlU_D(:) = 0.5*CurlU_D(:)/CellVolume_GB(i,j,k,iBlock)
    end if

  end subroutine get_curl_u

  !============================================================================

  subroutine apportion_coronal_heating(i, j, k, iBlock, &
       CoronalHeating, QeFraction, QparFraction) 

    ! Apportion the coronal heating to the electrons and protons based on
    ! how the Alfven waves dissipate at length scales << Lperp

    use ModMain, ONLY: UseB0
    use ModPhysics, ONLY: IonMassPerCharge
    use ModAdvance, ONLY: State_VGB, B0_DGB, UseAnisoPressure, &
         Rho_, Bx_, Bz_, Pe_, p_, Ppar_, WaveFirst_, WaveLast_

    integer, intent(in) :: i, j, k, iBlock
    real, intent(in) :: CoronalHeating
    real, intent(out) :: QeFraction, QparFraction

    real :: TeByTp, B2, BetaElectron, BetaProton, Pperp, LperpInvGyroRad
    real :: WaveLarge
    real :: DampingElectron, DampingPar, DampingPerp, DampingTotal

    character(len=*), parameter :: &
         NameSub = 'ModCoronalHeating::apportion_coronal_heating'
    !--------------------------------------------------------------------------

    if(UseTurbulentCascade)then
       ! Damping rates and wave energy partition based on Chandran et al.[2011]

       TeByTp = State_VGB(Pe_,i,j,k,iBlock) &
            /max(State_VGB(p_,i,j,k,iBlock), 1e-15)

       if(UseB0) then
          B2 = sum((B0_DGB(:,i,j,k,iBlock)+State_VGB(Bx_:Bz_,i,j,k,iBlock))**2)
       else
          B2 = sum(State_VGB(Bx_:Bz_,i,j,k,iBlock)**2)
       end if
       BetaElectron = 2.0*State_VGB(Pe_,i,j,k,iBlock)/max(B2, 1e-30)
       BetaProton = 2.0*State_VGB(p_,i,j,k,iBlock)/max(B2, 1e-30)

       ! Linear Landau damping and transit-time damping of kinetic Alfven
       ! waves contributes to electron and parallel proton heating
       DampingElectron = 0.01*sqrt(TeByTp/max(BetaProton, 1.0e-8)) &
            *(1.0 + 0.17*BetaProton**1.3)/(1.0 +(2800.0*BetaElectron)**(-1.25))
       DampingPar = 0.08*sqrt(sqrt(TeByTp))*BetaProton**0.7 &
            *exp(-1.3/max(BetaProton, 1.0e-8))

       ! Nonlinear damping/stochastic heating to perpendicular proton heating
       if(UseAnisoPressure)then
          Pperp = 0.5*(3*State_VGB(p_,i,j,k,iBlock) &
               - State_VGB(Ppar_,i,j,k,iBlock))
       else
          Pperp = State_VGB(p_,i,j,k,iBlock)
       end if

       ! The following is the fourth sqrt of Lperp by proton gyroradius
       LperpInvGyroRad = sqrt(sqrt(LperpTimesSqrtB* &
            sqrt(sqrt(B2)*State_VGB(Rho_,i,j,k,iBlock)/(2.0*Pperp)) &
            /IonMassPerCharge))

       ! Extract the Alfven wave energy density of the dominant wave
       WaveLarge = maxval(State_VGB(WaveFirst_:WaveLast_,i,j,k,iBlock))

       ! We note that if anistropic pressure is correctly included in
       ! the Alfven wave transport equation than the relationship between
       ! the wave energy density and velocity perturbation squared changes
       ! to: E_w = 0.5*rho*<\delta u^2> (1 + 1/A), where <.> is the
       ! ensemble average and A = 1 - (Ppar - Pperp)/B^2 in normalized
       ! units, resulting in a slightly different expression below.

       DampingPerp = 0.18*WaveLarge*sqrt(WaveLarge*B2/(2.0*Pperp)) &
            /IonMassPerCharge/max(CoronalHeating*LperpInvGyroRad**3,1e-30) &
            *exp(-StochasticHeating &
            *sqrt(2.0*Pperp/max(WaveLarge,1e-15))*LperpInvGyroRad)

       ! The 1+ is due to the fraction of the cascade power that succeeds
       ! to cascade to the smallest scale (<< proton gyroradius),
       ! where the dissipation is via interactions with the electrons
       DampingTotal = 1.0 + DampingElectron + DampingPar + DampingPerp

       QeFraction = (1.0 + DampingElectron)/DampingTotal
       QparFraction = DampingPar/Dampingtotal
    else
       QeFraction = QeByQtotal
       QparFraction = QparByQtotal
    end if

  end subroutine apportion_coronal_heating

end module ModCoronalHeating
