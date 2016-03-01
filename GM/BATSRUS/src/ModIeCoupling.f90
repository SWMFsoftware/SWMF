!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!============================================================================

module ModIeCoupling

  use ModProcMH

  implicit none
  save

  private ! except

  public:: init_ie_grid                ! set the ionosphere grid
  public:: clean_mod_ie_coupling       ! deallocate arrays
  public:: read_ie_velocity_param      ! parameters for velocity nudging
  public:: get_ie_potential            ! interpolate potential to x,y,z
  public:: logvar_ionosphere           ! get cross polar cap potential
  public:: calc_grad_ie_potential      ! calculate gradient of iono. potential
  public:: calc_inner_bc_velocity      ! get ExB/B^2 drift velocity at inner BC
  public:: apply_ie_velocity           ! nudge MHD velocity in a volume
  public:: map_jouleheating_to_inner_bc! map iono. Joule heating to inner BC
  public:: get_inner_bc_jouleheating   ! get mapped Joule heating at inner BC
  public:: calc_ie_mag_perturb         ! calculate dB due to IE currents

  ! Ionosphere grid description
  integer, public              :: nThetaIono = -1, nPhiIono = -1
  real,    public              :: rIonosphere
  real,    public, allocatable :: ThetaIono_I(:), PhiIono_I(:)

  ! Ionosphere potential and its gradient
  real, allocatable, public :: IonoPotential_II(:,:)
  real, allocatable, public :: dIonoPotential_DII(:,:,:)

  ! Joule heating
  real, public, allocatable :: IonoJouleHeating_II(:,:)

  ! Conductances
  real, public, allocatable :: SigmaHall_II(:,:), SigmaPedersen_II(:,:)

  ! Hall and Pedersen currents
  real, public, allocatable :: jHall_DII(:,:,:), jPedersen_DII(:,:,:)

  ! Local variables
  real:: dThetaIono, dPhiIono
  real, allocatable:: SinTheta_I(:), CosTheta_I(:), SinPhi_I(:), CosPhi_I(:)

  ! Velocity nudging
  logical :: UseIonoVelocity = .false.
  real    :: rCoupleUiono    = 3.5
  real    :: TauCoupleUiono  = 20.0

  ! Debugging parameters
  logical, parameter:: DoDebug = .false.
  integer, parameter:: iDebug = 10, jDebug = 10

contains
  !============================================================================
  subroutine init_ie_grid(Theta_I, Phi_I, rIono)

    use ModNumConst, ONLY: cPi, cTwoPi
    use ModPhysics, ONLY: Si2No_V, UnitX_

    ! Set internal variables for the ionosphere grid
    ! The ionosphere works on two hemispheres with a node based grid
    ! iSize is the number of latitude nodes from the pole to the equator.
    ! jSize is the number of longitude nodes (including a periodic overlap)
    
    real, intent(in):: Theta_I(:)  ! co-latitudes in radians
    real, intent(in):: Phi_I(:)    ! longitudes in radians
    real, intent(in):: rIono       ! Ionosphere radius in meters

    logical:: DoTest, DoTestMe
    character(len=*), parameter :: NameSub='init_ie_grid'
    !-----------------------------------------------------------------------
    if(nThetaIono > 1) RETURN

    call set_oktest(NameSub, DoTest, DoTestMe)

    nThetaIono = size(Theta_I)
    nPhiIono   = size(Phi_I)

    allocate(ThetaIono_I(nThetaIono), PhiIono_I(nPhiIono), &
         SinTheta_I(nThetaIono), CosTheta_I(nThetaIono), &    
         SinPhi_I(nPhiIono), CosPhi_I(nPhiIono))

    ThetaIono_I = Theta_I
    PhiIono_I   = Phi_I
    rIonosphere = rIono * Si2No_V(UnitX_)

    ! Save sin and cos of coordinates for sake of speed
    SinTheta_I = sin(ThetaIono_I)
    CosTheta_I = cos(ThetaIono_I)
    SinPhi_I   = sin(PhiIono_I)
    CosPhi_I   = cos(PhiIono_I)

    ! This only works for the uniform ionoosphere grid
    dThetaIono = cPi    / (nThetaIono - 1)
    dPhiIono   = cTwoPi / (nPhiIono - 1)

    if(DoTestMe)write(*,*) NameSub,': nThetaIono, nPhiIono=', &
         nThetaIono, nPhiIono

  end subroutine init_ie_grid
  !============================================================================
  subroutine clean_mod_ie_coupling

    if(allocated(IonoPotential_II)) &
         deallocate(IonoPotential_II, dIonoPotential_DII)

    if(allocated(IonoJouleHeating_II)) &
         deallocate(IonoJouleHeating_II)

    if(allocated(SigmaHall_II)) &
         deallocate(SigmaHall_II, SigmaPedersen_II)

    if(allocated(ThetaIono_I)) deallocate( &
         ThetaIono_I, PhiIono_I, SinTheta_I, CosTheta_I, SinPhi_I, CosPhi_I)

    nThetaIono = -1
    nPhiIono = -1

  end subroutine clean_mod_ie_coupling

  !===========================================================================
  subroutine get_ie_grid_index(Theta, Phi, ThetaNorm, PhiNorm)
    real, intent(in) :: Theta, Phi
    real, intent(out):: ThetaNorm, PhiNorm

    ThetaNorm = Theta / dThetaIono
    PhiNorm   = Phi   / dPhiIono
  end subroutine get_ie_grid_index

  !============================================================================
  subroutine get_ie_potential(Xyz_D, Potential)

    use ModMain,           ONLY: Time_Simulation, TypeCoordSystem
    use ModInterpolate,    ONLY: bilinear
    use CON_planet_field,  ONLY: map_planet_field
    use ModCoordTransform, ONLY: xyz_to_dir

    ! Interpolate IE potential to Xyz_D location
    ! Assume equipotential dipole field lines

    real, intent(in) :: Xyz_D(3)
    real, intent(out):: Potential

    integer:: iHemisphere
    real :: XyzIono_D(3)         ! Mapped point on the ionosphere
    real :: Theta, Phi           ! Mapped point colatitude, longitude

    logical, parameter :: DoTestMe = .false.
    character(len=*), parameter :: NameSub = 'get_ie_potential'
    !-------------------------------------------------------------------------
    ! call set_oktest(NameSub, DoTest, DoTestMe)

    if(.not.allocated(IonoPotential_II))then
       Potential = 0.0
!       Potential = sign(10.0, Xyz_D(2))
       RETURN
    end if

    ! Map down to the ionosphere at radius rIonosphere. Result is in SMG.
    call map_planet_field(Time_Simulation, Xyz_D, TypeCoordSystem//' NORM', &
         rIonosphere, XyzIono_D, iHemisphere, DoNotConvertBack=.true.)

    ! Calculate angular coordinates
    call xyz_to_dir(XyzIono_D, Theta, Phi)

    ! Interpolate the potential on the ionospheric grid
    Potential = bilinear(IonoPotential_II, 1, nThetaIono, 1, nPhiIono, &
         (/ Theta/dThetaIono+1, Phi/dPhiIono+1 /))

  end subroutine get_ie_potential

  !============================================================================

  subroutine calc_grad_ie_potential

    ! Calculate gradient of ionosphere potential on the IE grid

    integer, parameter :: Theta_=1, Phi_=2
    integer :: i, j
    !-------------------------------------------------------------------------
    if(.not.allocated(dIonoPotential_DII)) &
         allocate(dIonoPotential_DII(2,nThetaIono,nPhiIono))

    ! Calculate the gradients for the internal points with central differences
    do j = 1, nPhiIono; do i = 2, nThetaIono-1
       dIonoPotential_DII(Theta_, i, j) = &
            (IonoPotential_II(i+1, j) - IonoPotential_II(i-1, j)) &
            / (ThetaIono_I(i+1) - ThetaIono_I(i-1))
    end do; end do

    do j = 2, nPhiIono-1; do i = 1, nThetaIono
       dIonoPotential_DII(Phi_, i, j) = &
            (IonoPotential_II(i, j+1) - IonoPotential_II(i, j-1)) &
            / (PhiIono_I(j+1)-PhiIono_I(j-1))
    end do; end do

    ! Calculate the theta gradient at the poles
    ! with one sided second order approximations

    ! df/dx = (4f(x+dx)-3f(x)-f(x+2dx))/(2dx)
    dIonoPotential_DII(Theta_, 1, :) = &
         ( 4*IonoPotential_II(2,:) &
         - 3*IonoPotential_II(1,:) &
         -   IonoPotential_II(3,:) ) / (ThetaIono_I(3)-ThetaIono_I(1))

    ! df/dx = (3f(x)-4f(x-dx)+f(x-2dx))/(2dx)
    dIonoPotential_DII(Theta_, nThetaIono, :) = &
         ( 3*IonoPotential_II(nThetaIono  ,:) &
         - 4*IonoPotential_II(nThetaIono-1,:) &
         +   IonoPotential_II(nThetaIono-2,:) ) / &
         (ThetaIono_I(nThetaIono)-ThetaIono_I(nThetaIono-2))

    ! Calculate the phi gradient at the edges from the periodicity
    dIonoPotential_DII(Phi_, :, 1) = &
         (IonoPotential_II(:, 2) - IonoPotential_II(:, nPhiIono - 1)) &
         / (2*(PhiIono_I(2)-PhiIono_I(1)))

    dIonoPotential_DII(Phi_,:,nPhiIono) = dIonoPotential_DII(Phi_,:,1)

  end subroutine calc_grad_ie_potential

  !============================================================================
  real function logvar_ionosphere(NameLogvar)

    ! Calculate cross polar cap potentials for write_logvar

    use ModProcMH,  ONLY: iProc
    use ModIO,      ONLY: write_myname

    character (len=*), intent(in) :: NameLogvar
    integer :: nWarn = 0 ! warn multiple times to catch multiple log variables
    !--------------------------------------------------------------------------
    if(.not.allocated(IonoPotential_II))then
       logvar_ionosphere = 0.0
       return
    endif

    select case(NameLogvar)
    case('cpcpn','cpcp_n','cpcp_north','cpcpnorth')
       logvar_ionosphere = &
            maxval(IonoPotential_II(1:(nThetaIono+1)/2,:)) - &
            minval(IonoPotential_II(1:(nThetaIono+1)/2,:))
    case('cpcps','cpcp_s','cpcp_south','cpcpsouth')
       logvar_ionosphere = &
            maxval(IonoPotential_II((nThetaIono+1)/2:nThetaIono,:)) - &
            minval(IonoPotential_II((nThetaIono+1)/2:nThetaIono,:))
    case default
       if(nWarn < 2 .and. iProc==0)then
          call write_myname;
          write(*,'(a)')'WARNING in logvar_ionosphere: unknown NameLogvar='//&
               trim(NameLogvar)//', returning -777.77'
          nWarn = nWarn + 1
       end if
       logvar_ionosphere = -777.77
    end select

  end function logvar_ionosphere

  !============================================================================
  subroutine calc_inner_bc_velocity(tSimulation, Xyz_D, b_D, u_D)

    !USES:
    use ModMain,           ONLY: TypeCoordSystem, MaxDim
    use ModCoordTransform, ONLY: xyz_to_dir, cross_product
    use CON_planet_field,  ONLY: map_planet_field

    !INPUT ARGUMENTS:
    real, intent(in)    :: tSimulation      ! Simulation time
    real, intent(in)    :: Xyz_D(MaxDim)    ! Position vector
    real, intent(in)    :: b_D(MaxDim)      ! Magnetic field 

    !OUTPUT ARGUMENTS:
    real, intent(out)   :: u_D(MaxDim)      ! Velocity vector

    !DESCRIPTION:
    ! This subroutine calculates the velocity vector derived from
    ! the electric field of the ionosphere. The location is given by
    ! the Xyz\_D coordinates and the corresponding magnetic field is
    ! passed because it has been calculated already.
    ! The algorithm is the following: the input location is mapped down
    ! to the ionosphere where the Theta and Phi gradients of the potential
    ! are interpolated to the mapped point. This gradient is multiplied by
    ! the 2 by 3 Jacobian matrix of the mapping which converts the 
    ! Theta,Phi gradient to the X,Y,Z gradient of the potential, which is
    ! the electric field at the required location. The velocity is
    ! determined from the electric field and the magnetic field using
    ! the fact that the electric field is orthogonal to the magnetic field.
    !EOP

    real :: XyzIono_D(MaxDim)    ! Mapped point on the ionosphere
    real :: Theta, Phi           ! Mapped point colatitude, longitude
    real :: ThetaNorm, PhiNorm   ! Normalized colatitude, longitude
    real :: Dist1, Dist2         ! Distance from ionosphere grid point

    real :: dPotential_D(2)      ! Gradient of potential at the mapped position
    real :: DdirDxyz_DD(2,3)     ! Jacobian matrix between Theta, Phi and Xyz_D
    real :: eField_D(MaxDim)     ! Electric field
    real :: B2                   ! Magnetic field squared

    integer :: iTheta, iPhi, iHemisphere

    logical, parameter :: DoTestMe = .false.
    character(len=*), parameter :: NameSub = 'calc_inner_bc_velocity'
    !-------------------------------------------------------------------------
    ! call set_oktest(NameSub, DoTest, DoTestMe)

    if(DoTestMe)write(*,*)NameSub,' Xyz_D=',Xyz_D

    ! Map down to the ionosphere at radius rIonosphere. Result is in SMG.
    ! Also obtain the Jacobian matrix between Theta,Phi and Xyz_D
    call map_planet_field(tSimulation, Xyz_D, TypeCoordSystem//' NORM', &
         rIonosphere, XyzIono_D, iHemisphere, .true., DdirDxyz_DD)

    ! Calculate angular coordinates
    call xyz_to_dir(XyzIono_D, Theta, Phi)

    ! Interpolate the spherical gradients of the electric potential
    call get_ie_grid_index(Theta, Phi, ThetaNorm, PhiNorm)

    iTheta    = floor(ThetaNorm) + 1
    iPhi      = floor(PhiNorm)   + 1

    if(iTheta<1 .or. iTheta > nThetaIono .or. &
         iPhi < 1 .or. iPhi > nPhiIono)then
       write(*,*)NameSub,' PhiNorm, ThetaNorm=',PhiNorm,ThetaNorm
       write(*,*)NameSub,' Phi, Theta=',Phi,Theta
       write(*,*)NameSub,' nPhi, nTheta=',nPhiIono,nThetaIono
       write(*,*)NameSub,' iPhi, iTheta=',iPhi,iTheta
       call stop_mpi(NameSub//' index out of bounds')
    end if

    Dist1     = ThetaNorm - (iTheta - 1)
    Dist2     = PhiNorm   - (iPhi   - 1)

    dPotential_D = &
         (1 - Dist1)*( (1-Dist2) * dIonoPotential_DII(:, iTheta  , iPhi  )  &
         +             Dist2     * dIonoPotential_DII(:, iTheta  , iPhi+1)) &
         + Dist1    *( (1-Dist2) * dIonoPotential_DII(:, iTheta+1, iPhi  )  &
         +             Dist2     * dIonoPotential_DII(:, iTheta+1, iPhi+1))

    ! E = -grad(Potential) = - dPotential/d(Theta,Phi) * d(Theta,Phi)/d(x,y,z)
    eField_D = - matmul( dPotential_D, DdirDxyz_DD)

    ! Magnetic field
    B2  = sum(b_D**2)

    ! U = (E x B) / B^2
    u_D = cross_product(eField_D, b_D) / B2

    if(DoTestMe)then
       write(*,*)NameSub,' Xyz_D        =',Xyz_D
       write(*,*)NameSub,' XyzIono_D    =',XyzIono_D
       write(*,*)NameSub,' Theta, Phi   =',Theta,Phi
       write(*,*)NameSub,' iTheta, iPhi =',iTheta,iPhi
       write(*,*)NameSub,' Dist1, Dist2 =',Dist1,Dist2
       write(*,*)NameSub,' dPotential_D =',dPotential_D
       write(*,*)NameSub,' DdirDxyz_DD  =',DdirDxyz_DD
       write(*,*)NameSub,' E_D=',eField_D
       write(*,*)NameSub,' b_D=',b_D
       write(*,*)NameSub,' u_D=',u_D
    endif

    if(DoTestMe)write(*,*)NameSub,' Final u_D=',u_D

  end subroutine calc_inner_bc_velocity
  !==========================================================================
  subroutine map_jouleheating_to_inner_bc

    use ModMain,    ONLY: Time_Simulation
    use ModPhysics, ONLY: rBody
    use CON_planet_field, ONLY: get_planet_field, map_planet_field
    use ModCoordTransform, ONLY: sph_to_xyz

    integer :: i, j, iHemisphere
    real, dimension(3) :: XyzIono_D, bIono_D, B_D, Xyz_tmp
    real    :: bIono, b
    !-------------------------------------------------------------------------
    do j = 1, nPhiIono; do i = 1, nThetaIono
       call sph_to_xyz(rIonosphere, ThetaIono_I(i), PhiIono_I(j), XyzIono_D)
       call get_planet_field(Time_Simulation,XyzIono_D, 'SMG NORM', bIono_D)
       bIono = sqrt(sum(bIono_D**2))

       ! map out to GM (caution!, not like map down to the ionosphere, 
       ! there is always a corresponding position.)
       call map_planet_field(Time_Simulation, XyzIono_D, 'SMG NORM', &
            rBody, Xyz_tmp, iHemisphere)

       if (iHemisphere == 0) then 
          ! not a mapping in the dipole, but to the equator
          ! assume no outflow to GM inner boundary
          IonoJouleHeating_II(i,j) = 0
       else
          call get_planet_field(Time_Simulation, Xyz_tmp, 'SMG NORM', B_D)
          b = sqrt(sum(B_D**2))

          ! scale the jouleheating
          IonoJouleHeating_II(i,j) = IonoJouleHeating_II(i,j) * b/bIono

       endif
    end do; end do

  end subroutine map_jouleheating_to_inner_bc

  !==========================================================================
  subroutine get_inner_bc_jouleheating(XyzSm_D, JouleHeating)

    ! Get the Joule heating at the position XyzSm_D (given in SMG) 
    ! at the inner boundary.

    use ModCoordTransform, ONLY: xyz_to_dir
    use ModMain,           ONLY: Time_Simulation
    use CON_planet_field,  ONLY: map_planet_field

    !INPUT ARGUMENTS:
    real, intent(in)    :: XyzSm_D(3) ! Position vector in SMG
    real, intent(out)   :: JouleHeating    ! Houle heating

    real :: XyzIono_D(3)         ! Location mapped to ionosphere
    real :: Theta, Phi           ! Mapped point colatitude, longitude
    real :: ThetaNorm, PhiNorm   ! Normalized colatitude, longitude
    integer :: iTheta, iPhi      ! Indexes on the IE grid

    integer:: iHemisphere

    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub = 'map_inner_bc_jouleheating'
    !---------------------------------------------------------------------
    call set_oktest(NameSub, DoTest, DoTestMe)

    ! Map down to the ionosphere at radius rIonosphere
    call map_planet_field(Time_Simulation, XyzSm_D, 'SMG NORM', rIonosphere, &
         XyzIono_D, iHemisphere)

    ! Calculate angular coordinates
    call xyz_to_dir(XyzIono_D, Theta, Phi)

    ! Interpolate the spherical jouleheating
    call get_ie_grid_index(Theta, Phi, ThetaNorm, PhiNorm)
    iTheta    = floor(ThetaNorm) + 1
    iPhi      = floor(PhiNorm)   + 1
    if(iTheta<1 .or. iTheta > nThetaIono .or. &
         iPhi < 1 .or. iPhi > nPhiIono)then
       write(*,*)NameSub,' PhiNorm, ThetaNorm=',PhiNorm,ThetaNorm
       write(*,*)NameSub,' Phi, Theta=',Phi,Theta
       write(*,*)NameSub,' nPhi, nTheta=',nPhiIono,nThetaIono
       write(*,*)NameSub,' iPhi, iTheta=',iPhi,iTheta
       call stop_mpi(NameSub//' index out of bounds')
    end if

    ! Note: this is a first order accurate algorithm
    JouleHeating = IonoJouleHeating_II(iTheta, iPhi)

  end subroutine get_inner_bc_jouleheating

  !===========================================================================
  subroutine calc_ie_currents

    use CON_planet_field,  ONLY: get_planet_field
    use ModCoordTransform, ONLY: cross_product, sph_to_xyz
    use ModMain,           ONLY: Time_Simulation

    use ModPhysics, ONLY: No2Si_V, UnitX_, UnitElectric_, UnitJ_

    ! Calculate the ionospheric Hall and Pedersen currents 
    ! from the Hall and Pedersen conductivities and the electric field.

    real:: XyzIono_D(3), eTheta, ePhi, eIono_D(3), bUnit_D(3)

    integer:: i, j

    character(len=*), parameter:: NameSub = 'calc_ie_currents'
    !------------------------------------------------------------------------
    ! If the current arrays are allocated it means they are already calculated
    if(allocated(jHall_DII)) RETURN

    ! check if iono grid is defined
    if(.not.allocated(ThetaIono_I)) RETURN

    ! check if iono potential and conductances are known
    if(.not.allocated(dIonoPotential_DII)) RETURN
    if(.not.allocated(SigmaHall_II)) RETURN

    allocate( &
         jHall_DII(3,nThetaIono,nPhiIono), &
         jPedersen_DII(3,nThetaIono,nPhiIono))

    do j =1, nPhiIono; do i = 2, nThetaIono-1

       ! Calculate magnetic field direction for the Hall current
       call sph_to_xyz(rIonosphere, ThetaIono_I(i), PhiIono_I(j), XyzIono_D)
       call get_planet_field(Time_Simulation, XyzIono_D, 'SMG NORM', bUnit_D)
       bUnit_D = bUnit_D/sqrt(sum(bUnit_D**2))

       ! Calculate spherical components of the electric field from the
       ! derivatives of the potential in dIonoPotential_DII: E = -grad(Phi)
       eTheta = -dIonoPotential_DII(1,i,j)/rIonosphere
       ePhi   = -dIonoPotential_DII(2,i,j)/(rIonosphere*SinTheta_I(i))

       ! Convert to Cartesian components
       eIono_D(1) =  eTheta*CosTheta_I(i)*CosPhi_I(j) - ePhi*SinPhi_I(j)
       eIono_D(2) =  eTheta*CosTheta_I(i)*SinPhi_I(j) + ePhi*CosPhi_I(j)
       eIono_D(3) = -eTheta*SinTheta_I(i)

       ! Hall current for negative charge carriers: jH = -sigmaH*(E x B)
       jHall_DII(:,i,j) = -SigmaHall_II(i,j)*cross_product(eIono_D, bUnit_D)

       ! Perdersen current: jP = sigmaP*E
       jPedersen_DII(:,i,j) = SigmaPedersen_II(i,j)*eIono_D

       if(DoDebug.and.iProc==0.and.i==iDebug.and.j==jDebug)then
          write(*,*)NameSub,': i,j,Theta,Phi=', &
               i,j,ThetaIono_I(i),PhiIono_I(j)
          write(*,*)NameSub,': rIonosphere, dTheta, dPhi, sin(theta)=', &
               rIonosphere*No2Si_V(UnitX_), dThetaIono, dPhiIono, SinTheta_I(i)
          write(*,*)NameSub,': SigmaH,SigmaP=', &
               SigmaHall_II(i,j) &
               *No2Si_V(UnitJ_)*No2Si_V(UnitX_)/No2Si_V(UnitElectric_), &
               SigmaPedersen_II(i,j) &
               *No2Si_V(UnitJ_)*No2Si_V(UnitX_)/No2Si_V(UnitElectric_)
          write(*,*)NameSub,':  IonoPotential=', &
               IonoPotential_II(max(1,i-1):i+1,max(1,j-1):j+1) &
               *No2Si_V(UnitElectric_)*No2Si_V(UnitX_)
          write(*,*)NameSub,': dIonoPotential=', &
               dIonoPotential_DII(:,i,j) &
               *No2Si_V(UnitElectric_)/No2Si_V(UnitX_)
          write(*,*)NameSub,': eTheta, ePhi =', &
               eTheta*No2Si_V(UnitElectric_), &
               ePhi*No2Si_V(UnitElectric_)
          write(*,*)NameSub,': SinTheta,CosTheta,SinPhi,CosPhi=', &
               SinTheta_I(i), CosTheta_I(i), SinPhi_I(j), CosPhi_I(j)
          write(*,*)NameSub,': bUnit_D unit =', bUnit_D
          write(*,*)NameSub,': eIono_D      =', eIono_D*No2Si_V(UnitElectric_)
          ! These are height integrated currents in SI units of A/m
          write(*,*)NameSub,': Jhall        =', jHall_DII(:,i,j) &
               *No2Si_V(UnitJ_)*No2Si_V(UnitX_)
          write(*,*)NameSub,': JPedersen    =', jPedersen_DII(:,i,j) &
               *No2Si_V(UnitJ_)*No2Si_V(UnitX_)
       end if

    end do; end do

  end subroutine calc_ie_currents

  !===========================================================================

  subroutine calc_ie_mag_perturb(nMag, XyzSm_DI, dBHall_DI, dBPedersen_DI)

    ! For nMag points at locations XyzSm_DI (in SMG coordinates)
    ! calculate the magnetic pertubations dBHall_DI and dBPedersen_DI
    ! from the ionospheric Pedersen and Hall currents, respectively
    ! using the Biot-Savart integral. The output is given in SMG coordinates.

    use ModNumConst, ONLY: cPi
    use ModCoordTransform, ONLY: cross_product, sph_to_xyz
    use ModPhysics, ONLY: No2Si_V, UnitB_, UnitX_
    use ModMain,    ONLY: Time_Simulation

    integer,intent(in) :: nMag
    real,   intent(in) :: XyzSm_DI(3,nMag)
    real,   intent(out):: dBHall_DI(3,nMag), dBPedersen_DI(3,nMag)

    real:: XyzIono_D(3), dXyz_D(3)
    integer:: iMag, i, j, iLine
    real:: Coef0, Coef
    
    logical:: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'calc_ie_mag_perturb'
    !========================================================================
    call timing_start(NameSub)
    call set_oktest(NameSub, DoTest, DoTestMe)

    if(DoTestMe)write(*,*) NameSub,' starting with XyzSm(iMag=1)=', &
         XyzSm_DI(:,1)*No2Si_V(UnitX_)

    ! Initialize magnetic perturbations
    dBHall_DI     = 0.0
    dBPedersen_DI = 0.0

    ! Check if IE has provided the necessary information
    if(.not.allocated(SigmaHall_II)) RETURN

    ! Calculate the height integrated Hall and Pedersen current densities
    call calc_ie_currents

    ! Loop over the ionosphere grid, but skip the poles
    iLine = 0
    do j = 1, nPhiIono; do i = 2, nThetaIono-1

       iLine = iLine + 1
       ! distribute the work among the BATSRUS processors
       if(mod(iLine, nProc) /= iProc)CYCLE

       ! Get Cartesian coordinates for the ionospheric grid point
       call sph_to_xyz(rIonosphere, ThetaIono_I(i), PhiIono_I(j), XyzIono_D)

       ! 1/4pi times the area of a surface element
       Coef0 = 1/(4*cPi)*rIonosphere**2*dThetaIono*dPhiIono*SinTheta_I(i)

       do iMag = 1, nMag
          ! Distance vector between magnetometer position 
          ! and ionosphere surface element
          dXyz_D = XyzSm_DI(:,iMag) - XyzIono_D

          ! Surface element area divided by (4pi*distance cubed)
          Coef = Coef0/sqrt(sum(dXyz_D**2))**3

          ! Do Biot-Savart integral: dB = j x d/(4pi|d|^3) dA  (mu0=1)
          dBHall_DI(:,iMag)     = dBHall_DI(:,iMag) + &
               Coef*cross_product(jHall_DII(:,i,j), dXyz_D)
          dBPedersen_DI(:,iMag) = dBPedersen_DI(:,iMag) + &
               Coef*cross_product(jPedersen_DII(:,i,j), dXyz_D)

          if(DoDebug.and.iProc==0.and.i==iDebug.and.j==jDebug.and.iMag==1)then
             write(*,*)NameSub,': Time=', Time_Simulation
             write(*,*)NameSub,': XyzSm,XyzIono=', &
                  XyzSm_DI(:,iMag)*No2Si_V(UnitX_), XyzIono_D*No2Si_V(UnitX_)
             write(*,*)NameSub,': rIono**2, dTheta, dPhi, sinTheta=', &
                   rIonosphere**2*No2Si_V(UnitX_)**2, &
                   dThetaIono, dPhiIono, SinTheta_I(i)
             write(*,*)NameSub,': dArea, r3    =', &
                  rIonosphere**2*dThetaIono*dPhiIono*SinTheta_I(i)&
                  *No2Si_V(UnitX_)**2, &
                  sqrt(sum(dXyz_D**2))**3*No2Si_V(UnitX_)**3
             write(*,*)NameSub,': dBHall,  sum =', &
                  Coef*cross_product(jHall_DII(:,i,j), dXyz_D) &
                  *No2Si_V(UnitB_), dBHall_DI(:,iMag)*No2Si_V(UnitB_)
             write(*,*)NameSub,': dBPeders,sum=', &
                  Coef*cross_product(jPedersen_DII(:,i,j), dXyz_D) &
                  *No2Si_V(UnitB_), dBPedersen_DI(:,iMag)*No2Si_V(UnitB_)
          end if

       end do; end do
    end do
    call timing_stop(NameSub)

    if(DoTestMe)write(*,*) NameSub,': dBHall(iMag=1), dBPede(iMag=1)=', &
         dBHall_DI(:,1)*No2Si_V(UnitB_), &
         dBPedersen_DI(:,1)*No2Si_V(UnitB_)

  end subroutine calc_ie_mag_perturb
  !===========================================================================
  subroutine read_ie_velocity_param

    use ModReadParam, ONLY: read_var

    call read_var('UseIonoVelocity', UseIonoVelocity)
    if(UseIonoVelocity)then
       call read_var('rCoupleUiono'   , rCoupleUiono)
       call read_var('TauCoupleUiono',  TauCoupleUiono)
    end if

  end subroutine read_ie_velocity_param

  !===========================================================================
  subroutine apply_ie_velocity

    use ModMain,    ONLY: nI, nJ, nK, nBlock, time_accurate, time_simulation, &
         Dt, Unused_B, UseB0, UseRotatingBc
    use ModAdvance, ONLY: State_VGB, Rho_, RhoUx_, RhoUz_, Bx_, Bz_
    use ModGeometry,ONLY: r_BLK, Rmin_BLK
    use ModB0,      ONLY: B0_DGB
    use ModPhysics, ONLY: Si2No_V, UnitT_, rBody
    use ModEnergy,  ONLY: calc_energy_cell
    use BATL_lib,   ONLY: Xyz_DGB

    real :: Factor, RhoUdotB
    real, dimension(3) :: Xyz_D, Uiono_D, Urot_D, RhoUiono_D, RhoU_D, b_D

    integer :: i,j,k,iBlock
    character (len=*), parameter :: NameSub='apply_ie_velocity'
    logical :: DoTest, DoTestMe
    !--------------------------------------------------------------------------
    if(.not.UseIonoVelocity) RETURN

    call set_oktest(NameSub, DoTest, DoTestMe)

    if(time_accurate)then
       ! Ramp up is based on physical time: u' = u + dt/tau * (uIE - u)
       ! A typical value might be 10 sec, to get close to the IE velocity
       ! in 20 seconds

       Factor = min(1.0, Dt/(TauCoupleUiono*Si2No_V(UnitT_)))

    else
       ! Ramp up is based on number of iterations: u' = u + (uIE-u)/(1+nTau)
       ! A typical value might be 10, to get close to the E x B velocity
       ! in 20 iterations

       Factor = 1.0/(1.0 + TauCoupleUiono)

    end if

    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       if(rMin_BLK(iBlock) > rCoupleUiono) CYCLE

       do k = 1, nK; do j = 1, nJ; do i = 1, nI

          if( r_BLK(i,j,k,iBlock) > rCoupleUiono) CYCLE
          if( r_BLK(i,j,k,iBlock) < rBody) CYCLE

          Xyz_D = Xyz_DGB(:,i,j,k,iBlock)

          b_D = State_VGB(Bx_:Bz_,i,j,k,iBlock)

          if(UseB0)b_D = b_D + b0_DGB(:,i,j,k,iBlock)

          ! Calculate E x B velocity
          call calc_inner_bc_velocity(Time_Simulation, Xyz_D, b_D, uIono_D)

          ! Add rotational velocity if required
          if (UseRotatingBc) then
             call calc_corotation_velocities(Xyz_D, uRot_D)
             uIono_D = uIono_D + uRot_D
          end if

          ! Convert to momentum
          RhoUIono_D = State_VGB(Rho_,i,j,k,iBlock)*uIono_D

          ! The original momentum
          RhoU_D = State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock)

          ! Store field aligned momentum component
          RhoUdotB = sum(RhoU_D*b_D)

          ! Push momenta towards the Rho*(E x B + Omega x r) value
          RhoU_D = RhoU_D + Factor*(RhoUIono_D - RhoU_D)

          ! Restore field aligned momentum component
          RhoU_D = RhoU_D + b_D*(RhoUdotB - sum(RhoU_D*b_D))/sum(b_D**2)

          ! Store result
          State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) = RhoU_D

       end do; end do; end do

       ! Recalculate the energy
       call calc_energy_cell(iBlock)

    end do

  end subroutine apply_ie_velocity

end module ModIeCoupling
