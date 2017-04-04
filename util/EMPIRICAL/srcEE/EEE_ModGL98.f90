!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!==============================================================================
module EEE_ModGL98
  use EEE_ModCommonVariables
  implicit none
  save

  private

  public :: set_parameters_GL98
  public :: get_GL98_fluxrope
  public :: adjust_GL98_fluxrope

  real :: cme_a=0.0, cme_r1=0.0, cme_r0=0.0, cme_a1=0.0, cme_alpha=0.0
  real :: pBackgroundDim =0.0, cme_rho2=0.0

  real :: ModulationRho=0.0, ModulationP=0.0
  !\
  !The derivative of current over flux function
  !(\mu_0)dI/d\psi) has the dimentions of inverse length
  !/
  real :: Alpha0
  !\
  ! Characteristic value of magnetic field of the spheromak
  ! configuration: the field in the center of configuration
  ! equals 2(1/3 -\beta_0)B_0 \approx 0.7 B_0:
  !/ 
  real :: B0
  !\
  ! Characteristic ratio of plasma pressure to B_0^2/\mu_0
  ! Exact definition in terms of the pressure derivative over
  ! the flux function: \beta_0=\mu_0/(B_0 Alpha_0^2) dp/d\psi 
  !/
  real :: Beta0
  !\
  ! Dimensionless product of R0 by Alpha0. For any given \beta_0,
  ! this product is found from the equation, 
  !
  ! j_1(\alpha_0r_0)/(\alpha_0r_0)=\beta_0, (*)
  !
  ! where j_1(x)=sin(x)/x^2 -cos(x)/x is the spherical Bessel function.
  ! With this equation satisfied, the radial and phi- components of the
  ! spheromak field as well as the flux function, current function and
  ! pressure all turn to zero at the external boundary of configuration, 
  ! at r=r_0, while theta-component of the field is non-zero and 
  ! proportional to j_2(\alpha_0r_0). The GL paper chooses a special case of
  ! (negative) \beta_0. Specifically, the product \alpha_0r_0 is chosen
  ! to be a first zero of j_2 function, and then \beta_0 is chosen to satisfy
  ! Eq (*) with this value of .     
  !/
  real :: Alpha0R0
  !\
  ! Vector characteristic of the configuration: radius vector of the
  ! configuration center and B0 multiplied by unit vector along
  ! the axis of symmetry
  !/
  real,dimension(3) :: XyzCenterConf_D, BConf_D
  !\
  !Misc
  !/
  real :: delta
  real :: pBackground,rho2scl 
contains
  !=========================================
  subroutine set_parameters_GL98(NameCommand)
    use ModReadParam, ONLY: read_var
    character(len=*), intent(in):: NameCommand

    character(len=*), parameter:: NameSub = 'set_parameters_GL98'
    !--------------------------------------------------------------------------
    select case(NameCommand)
    case("#GL98FLUXROPE")
       call read_var('UseFluxRope',     DoAddFluxRope)
       call read_var('cme_a',           cme_a)
       call read_var('cme_r1',          cme_r1)
       call read_var('cme_r0',          cme_r0)
       call read_var('cme_a1',          cme_a1)
       call read_var('cme_alpha',       cme_alpha)
       call read_var('pBackground',     pBackgroundDim)
       call read_var('cme_rho2',        cme_rho2)
       call read_var('ModulationRho',   ModulationRho)
       call read_var('ModulationP',     ModulationP)
       call read_var('LongitudeCme',   LongitudeCme)
       call read_var('LatitudeCme',    LatitudeCme)
       call read_var('OrientationCme', OrientationCme)
    case("#CME")
       call read_var('Stretch',     cme_a)
       call read_var('Distance',    cme_r1)
       call read_var('Radius',      cme_r0)
       call read_var('Bstrength',   cme_a1)
       call read_var('pBackground',     pBackgroundDim)
       call read_var('ModulationRho',   ModulationRho)
       call read_var('ModulationP',     ModulationP)
    case default
       call CON_stop(NameSub//' unknown NameCommand='//NameCommand)
    end select

  end subroutine set_parameters_GL98

  !============================================================================

  subroutine get_GL98_fluxrope(XyzIn_D,rho_GL98,p_GL98,B_D,U_D)
    !--------------------------------------------------------------------------
    ! Definition of Parameters used for the initial state
    !   cme_a    = contraction distance as in   r --> r -a
    !   cme_r1   = distance of flux rope from sun center = 1.2
    !   cme_r0   = radius of flux rope
    !   cme_a1   = constant for setting pressure in flux rope
    !
    ! OrientationCME : The counter-clockwise angle between the fluxrope
    !                  axis and the solar equator.
    ! cme_a1 : The sign of cme_a1 (magnitude of the fluxrope field strength)
    !          sets the helicity Positive (negative) value results in
    !          postive or left handed or sinistral (negative or right handed
    !          or dextral) helicity.
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////
    !=====================================================================
    !\
    ! Calculates magnetic field, pressure and density for a coronal flux
    ! rope capable of self-similar expansion and erupting in a CME.
    ! The analytical solution is taken from Gibson and Low 
    ! Astrophysical Journal, Vol 493, p. 460.
    !
    ! Written by Chip Manchester Jan 18 2001
    ! Rewritten by Chip Nov 29 2001 for flux rope injection
    !/
    !   Bug fixes
    !   March 18       dpres_1dr1 = cme_a1*(dA2dr*dr2dr1 + dA2dth*dth2dr1)
    !   change to..... dpres_1dr1 = a1scl*(dA2dr*dr2dr1 + dA2dth*dth2dr1)
    !   without above fix, same as used for runs 12, 13, 14
    !   Feb  07, 2002 Br1 is changed to Br1**2 in density calc thanks to Ilia
    !   Feb  12, 2002 expression for ga0r0 is fixed     
    !   Feb  22, 2002 derivative of B2*2/dr1 is fixed
    !   Feb  22, 2002 angles in 2nd coordinate system are redefined
    !
    !Igor 2017-03-03 Bug is found. Formulae from the Gibson-Low paper 
    !were derived in CGS system with the use of 1/4\pi factor in the 
    !magnetic pressure and magnetic force density. In the existing 
    !version of the eruptive event generator, the dimensionless formulation 
    !of the MHD equations are employed, with no 1/4\pi coefficient. To keep 
    !some backward compatibility, the coefficient 4\pi is retained in the
    !definition of a1scl, hence old input parameters will provide same
    !distribution for the magnetic field, however, the perturbations of
    !pressure and density will be a factor of 4\pi times larger.
    !
    !Igor 2017-03-04 Fix a bug  Before the 'helicity switch' was
    !implemented as follows:
    !if(cme_a1 > 0.0) then;Br2 = -Br2; Btheta2 = -Btheta2; end if
    !However, below, the magnetic pressure gradient
    !is calculated including Br2*dBr2dr1 + Btheta2*dBtheta2dr1.
    !Therefore, the sign in the magnetic field derivatives should
    !be changed too. Add two more operators prior to end if:
    !dBr2dr1     = -dBr2dr1; dBtheta2dr1 = -dBtheta2dr1
    !------------------------------------------------------------------------
    use ModNumConst,       ONLY: cPi, cDegToRad
    use ModCoordTransform, ONLY: rot_matrix_y, rot_matrix_z, cross_product
    implicit none
    !\
    ! Coordinates of the input point
    !/
    real, intent(in) :: XyzIn_D(3)
    !\
    ! OUTPUTS
    !/
    !\
    ! Magnetic field
    !/
    real, intent(out) :: B_D(3)
    !\
    ! Optional: velocity of self-similar expansion
    !/
    real, intent(out),optional :: U_D(3)
    !\
    ! Density, pressure
    !/
    real, intent(out) :: rho_GL98,p_GL98
    !\
    ! User declared local variables go here::
    !/
    !\
    ! Radial coordinate of the input point
    !/
    real :: R
    !\
    ! Coordinates of the input point after
    ! stretching transformation
    !/
    real :: RTransformed, XyzTransformed_D(3)
    !\
    ! Ccordinate vector originating at the center of
    ! the magnetic configuration and its module
    !/
    real :: Distance2ConfCenter,XyzConf_D(3)
    !\
    ! Local pressure in the magnetic configuration
    !/ 
    real :: PConf
    !\
    ! Variables needed to calculate parameters of the
    ! stretched magnetic configuration:
    !/
    ! 1. Unit vector of radial direction and the magnetic 
    ! field projection onto it:
    real :: eBoldR_D(3), Br1
    ! 2. Magnetic field squared, total pressure gradient
    ! and its projection on the radial direction 
    real ::  BSquared, DPTotalDR1, GradPTotal_D(3)
    ! 3. Radial tension of stretched field
    real :: RadialTension
    ! 4. Local acceleration of the gravitational force
    real :: gGravity
    !\ 
    ! Misc
    !/
    real :: Alpha0R2
    real, dimension(3) :: R2CrossB0_D
    real, dimension(3,3):: RotateGL98_DD
    logical, save :: DoFirst_GL=.true.
    !------------------------------------------------------------------------
    if (DoFirst_GL) then
       DoFirst_GL=.false.
       if(iProc==0)then
          write(*,*) prefix
          write(*,*) prefix, &
               '>>>>>>>>>>>>>>>>>>>                  <<<<<<<<<<<<<<<<<<<<<'
          write(*,*) prefix
          write(*,*) prefix, &
               '            Gibson and Low CME is initiated!!!'
          write(*,*) prefix
          write(*,*) prefix, &
               '>>>>>>>>>>>>>>>>>>>                  <<<<<<<<<<<<<<<<<<<<<'
          write(*,*) prefix
          write(*,*) prefix, 'cme_a  = ',cme_a, '[rSun]'
          write(*,*) prefix, 'cme_r1 = ',cme_r1,'[rSun]'
          write(*,*) prefix, 'cme_r0 = ',cme_r0,'[rSun]'
          write(*,*) prefix, 'cme_a1 = ',cme_a1,'[Gauss/rSun^2]'
          write(*,*) prefix, 'pBackground = ', pBackgroundDim, '[dyne/cm^2]'
          write(*,*) prefix, 'cme_rho2 = ',cme_rho2,'[kg/m^3]'
          write(*,*) prefix, 'ModulationRho = ',ModulationRho
          write(*,*) prefix, 'ModulationP   = ',ModulationP
          write(*,*) prefix, 'LongitudeCme = ',LongitudeCme,'[degrees]'
          write(*,*) prefix, 'LatitudeCme = ',LatitudeCme,'[degrees]'
          write(*,*) prefix, 'OrientationCme = ',OrientationCme,'[degrees]'
          write(*,*) prefix
       end if
       pBackground = pBackgroundDim*Si2No_V(UnitP_)
       rho2scl = cme_rho2*Si2No_V(UnitRho_)
       delta = 0.1
       Alpha0R0 = 5.763854
       alpha0 = Alpha0R0/cme_r0
       Beta0 = (sin(Alpha0R0) - Alpha0R0*cos(Alpha0R0))/Alpha0R0**3
       !\
       ! The constant coefficient, Beta0 = -2.8723629148938019E-02 
       ! This is Beta0 parameter for the GL particular configuration
       !/
       !/
       !\
       ! 4\pi in the formula below is incorrect. The GL98 paper is in
       ! CGS system, while in dimensionless formulae used below it may not
       ! appear. To mantain the capability with EEGGL, the multiplier 4\pi
       ! is included into a definition of A1Scaled
       !/    
       B0 = - cme_a1*Io2No_V(UnitB_)/(Beta0*Alpha0**2) &
            *4.0*cPi
       !\
       ! The costant coefficient,
       ! -4.0*cPi/(Beta0*Alpha0R0**2) = 13.1687517342067082
       !/
       !\
       ! Construct the rotational matrix RotateGL98_DD,
       !/
       RotateGL98_DD  = matmul( &
            rot_matrix_z(-OrientationCme*cDegToRad),&
            rot_matrix_y((LatitudeCme-90)*cDegToRad))
       RotateGL98_DD = matmul(RotateGL98_DD, &
            rot_matrix_z(-LongitudeCme*cDegToRad))
       !\
       ! In the rotated coordinates the coordinate vector from 
       ! the heiocenter to the center of configuration is
       ! (/0.0, 0.0, cme_r1/). Find this vector in the original
       ! coodinate frame.
       !/
       XyzCenterConf_D = matmul((/0.0, 0.0, cme_r1/), RotateGL98_DD)
       !\
       ! The same for the magnetic field of configuration
       !/
       BConf_D = matmul((/B0, 0.0, 0.0/), RotateGL98_DD)
    end if

    R = sqrt(sum(XyzIn_D**2))
    ! Unit vector of radial direction
    eBoldR_D = XyzIn_D/R
    !\
    ! Stretched CARTESIAN COORDINATES 
    !/
    RTransformed = R + cme_a
    XyzTransformed_D = eBoldR_D*RTransformed
    !COORDINATES RELATIVE TO CME FLUX ROPE CENTER
    !----------------------------------------------------------------
    XyzConf_D = XyzTransformed_D - XyzCenterConf_D 
    Distance2ConfCenter = sqrt(sum(XyzConf_D**2))
    if (Distance2ConfCenter <= delta) then
       XyzConf_D = XyzConf_D*(delta/Distance2ConfCenter)
       Distance2ConfCenter = delta
    end if
    if (Distance2ConfCenter <= cme_r0) then
       !\
       ! INSIDE THE MAGNETIC FLUX ROPE REGION
       !/
       !\
       !An argument of spherical Bessel functions
       !/
       Alpha0R2 = Alpha0*Distance2ConfCenter 
       !\
       ! Magnetic field components
       !/ 
       R2CrossB0_D = cross_product(XyzConf_D,BConf_D) 
       B_D = (2.0*BConf_D &
            + sign(Alpha0,cme_a1)* &  !Helicity
            R2CrossB0_D)*(spher_bessel_1_over_x(Alpha0R2) - Beta0) &
            + spher_bessel_2(Alpha0R2)/Distance2ConfCenter**2*&
            cross_product(XyzConf_D,R2CrossB0_D) 
       !\
       ! Compute kinetic gas pressure
       !/
       PConf     = pBackground + Beta0*(Alpha0**2)*&
            sum(R2CrossB0_D**2)*(spher_bessel_1_over_x(Alpha0R2) - Beta0)
       !If we use stretch:
       Br1 = sum(eBoldR_D*B_D) 
       BSquared = sum(B_D**2)
       GradPTotal_D =&
            Alpha0**2*(B0**2*XyzConf_D - BConf_D*sum(BConf_D*XyzConf_D))*&
            (spher_bessel_1_over_x(Alpha0R2)**2 - Beta0**2)&
            -Alpha0*(XyzConf_D/Distance2ConfCenter**3)*sum(R2CrossB0_D**2)*&
            spher_bessel_2(Alpha0R2)*spher_bessel_1(Alpha0R2) &
            -4.0*B0**2*(XyzConf_D/Distance2ConfCenter**2)*&
            spher_bessel_2(Alpha0R2)*(spher_bessel_1_over_x(Alpha0R2) - Beta0)&
            +(XyzConf_D*sum(BConf_D*XyzConf_D)**2/Distance2ConfCenter**4 - &
            BConf_D*sum(BConf_D*XyzConf_D)/Distance2ConfCenter**2&
            )*(spher_bessel_2(Alpha0R2)**2 - 4.0*spher_bessel_2(Alpha0R2)*&
            (spher_bessel_1_over_x(Alpha0R2) - Beta0)) +&
            Alpha0*(XyzConf_D/Distance2ConfCenter**3)*sum(R2CrossB0_D**2)*(&
            (spher_bessel_1(Alpha0R2) - 3.0*spher_bessel_2(Alpha0R2)/Alpha0R2)*&
            (spher_bessel_2(Alpha0R2) - 2.0*&
            (spher_bessel_1_over_x(Alpha0R2) - Beta0)) + &
            2.0*spher_bessel_2(Alpha0R2)**2/Alpha0R2 )        
       DPTotalDR1 = sum(GradPTotal_D*eBoldR_D)
       !\
       ! MAGNETIC FIELD transformed with stretching transformation
       !/
       B_D = (RTransformed/r)*B_D + cme_a*XyzTransformed_D*Br1/r**2
       !\
       ! PLASMA DENSITY with stretching transformation
       !/
       gGravity   = (abs(Gbody)/(r**2) + cme_alpha*r)
       RadialTension = (cme_a/r)*(RTransformed/r)**2*(&
            (2.0 + cme_a/r)*(BSquared/RTransformed + DPTotalDR1) &
           + 2.0*PConf/RTransformed - (3.0 + 2.0*cme_a/r)*Br1**2/r)
       rho_GL98 = RadialTension/gGravity
       !\
       ! Add background density
       !/
       rho_GL98 = rho_GL98 + rho2scl/r**3
       !\
       ! PLASMA PRESSURE with contraction transformation
       !/
       p_GL98 = PConf*(RTransformed/r)**2 &
            - (1.0/2.0)*((RTransformed/r)**2)*((RTransformed/r)**2 - 1.0)*Br1**2 
       !\
       ! Add background pressure
       !/
       p_GL98 = p_GL98 + abs(Gbody)*rho2scl/(4.0*r**4)

       rho_GL98 = rho_GL98*No2Si_V(UnitRho_)
       p_GL98 = p_GL98*No2Si_V(UnitP_)
       B_D = B_D*No2Si_V(UnitB_)

    else
       B_D = 0.0; rho_GL98 = 0.0; p_GL98 = 0.0
    endif
  contains
    real function spher_bessel_0(x)
      real, intent(in) :: x
      spher_bessel_0 = sin(x)/x
    end function spher_bessel_0
    !===================
    real function spher_bessel_1(x)
      real, intent(in) :: x
      spher_bessel_1 = (sin(x) - x*cos(x))/x**2
    end function spher_bessel_1
    !===================
    real function spher_bessel_1_over_x(x)
      real, intent(in) :: x
      if(x==0.0)then
         spher_bessel_1_over_x = 1.0/3.0
      else
         spher_bessel_1_over_x = (sin(x) - x*cos(x))/x**3
      end if
    end function spher_bessel_1_over_x
    !===================
    real function spher_bessel_2(x)
      real, intent(in) :: x
      spher_bessel_2 = &
           3.0*spher_bessel_1_over_x(x) - spher_bessel_0(x)
    end function spher_bessel_2
    !===================
  end subroutine get_GL98_fluxrope

  !===================================================

  subroutine adjust_GL98_fluxrope(Rho,p)

    real, intent(inout) :: Rho,p
    !--------------------------------------------------------------------------

    !\
    ! Add just `ModulationRho' times of the CME mass
    ! to the mass density::
    !/
    if(ModulationRho*Rho <= 0.0)then
       Rho = Rho
    else
       Rho = ModulationRho*Rho
    end if

    !\
    ! Add just `ModulationP' times of the CME pressure
    ! to the kinetic pressure::
    !/
    if(ModulationP*p <= 0.0) then
       p = p
    else
       p = ModulationP*p
    end if

  end subroutine adjust_GL98_fluxrope

end module EEE_ModGL98
