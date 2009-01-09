!^CFG COPYRIGHT UM
!^CMP FILE IE
!==========================================================================
module ModIonoPotential

  use ModIeGrid
  implicit none
  save

  real, allocatable :: IonoPotential_II(:,:)
  real, allocatable :: dIonoPotential_DII(:,:,:)

contains

  subroutine init_mod_iono_potential(iSize, jSize)

    integer, intent(in) :: iSize, jSize
    character(len=*), parameter :: NameSub='init_mod_iono_potential'
    !-------------------------------------------------------------------------

    if(allocated(IonoPotential_II)) RETURN

    call init_mod_ie_grid(iSize, jSize)

    allocate( IonoPotential_II(nThetaIono, nPhiIono), &
         dIonoPotential_DII(2, nThetaIono, nPhiIono) )

  end subroutine init_mod_iono_potential

  !============================================================================

  subroutine calc_grad_iono_potential

    integer, parameter :: Theta_=1, Phi_=2
    integer :: i, j

    !write(*,*)'calc_grad_iono_potential: overwriting potential !!!'
    !do j = 1, nPhiIono;do i = 1, nThetaIono
    !   IonoPotential_II(i,j) = 0.1*i**2 + 0.01*j**2
    !end do; end do

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

  end subroutine calc_grad_iono_potential

end module ModIonoPotential

!=============================================================================
subroutine GM_put_from_ie(Buffer_II,iSize,jSize)

  use ModIonoPotential, ONLY: IonoPotential_II, &
       init_mod_iono_potential, calc_grad_iono_potential
  use ModPhysics,       ONLY: Si2No_V, UnitX_, UnitElectric_
  implicit none
  character(len=*), parameter :: NameSub='GM_put_from_ie'

  integer, intent(in) :: iSize,jSize
  real, intent(in) :: Buffer_II(iSize,jSize)
!  character(len=*), intent(in) :: NameVar

  logical :: DoTest, DoTestMe
  !---------------------------------------------------------------------------
  call CON_set_do_test(NameSub,DoTest,DoTestMe)
!  if(DoTest)write(*,*)NameSub,': NameVar,iSize,jSize=',NameVar,iSize,jSize

  if(.not. allocated(IonoPotential_II)) &
       call init_mod_iono_potential(iSize,jSize)

  IonoPotential_II = Buffer_II * Si2No_V(UnitElectric_)*Si2No_V(UnitX_)
  call calc_grad_iono_potential

  if(DoTest)write(*,*)NameSub,': done'

end subroutine GM_put_from_ie

!==========================================================================
!BOP
!ROUTINE: calc_inner_bc_velocity - calculate velocity at the inner boundary
!INTERFACE:
subroutine calc_inner_bc_velocity(tSimulation, Xyz_D, b_D, u_D)

  !USES:
  use ModIonoPotential
  use ModMain,           ONLY: TypeCoordSystem, nDim
  use CON_axes,          ONLY: transform_matrix
  use ModCoordTransform, ONLY: xyz_to_dir, cross_product
  use CON_planet_field,  ONLY: map_planet_field

  implicit none

  !INPUT ARGUMENTS:
  real, intent(in)    :: tSimulation    ! Simulation time
  real, intent(in)    :: Xyz_D(nDim)    ! Position vector
  real, intent(in)    :: b_D(nDim)      ! Magnetic field 

  !OUTPUT ARGUMENTS:
  real, intent(out)   :: u_D(nDim)      ! Velocity vector


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

  real :: XyzIono_D(nDim)      ! Mapped point on the ionosphere
  real :: Theta, Phi           ! Mapped point colatitude, longitude
  real :: ThetaNorm, PhiNorm   ! Normalized colatitude, longitude
  real :: Dist1, Dist2         ! Distance from ionosphere grid point

  real :: dPotential_D(2)      ! Gradient of potential at the mapped position
  real :: DdirDxyz_DD(2,3)     ! Jacobian matrix between Theta, Phi and Xyz_D
  real :: eField_D(nDim)       ! Electric field
  real :: B2                   ! Magnetic field squared

  integer :: iTheta, iPhi, iHemisphere

  character(len=*), parameter :: NameSub = 'calc_inner_bc_velocity'
  logical :: DoTest = .false., DoTestMe = .false.
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
subroutine calc_inner_bc_velocity1(tSimulation,Xyz_D,B1_D,B0_D,u_D)

  use ModIonoPotential
  use ModMain,           ONLY: TypeCoordSystem, nDim
  use CON_axes,          ONLY: transform_matrix
  use ModCoordTransform, ONLY: xyz_to_dir, cross_product
  use CON_planet_field,  ONLY: map_planet_field

  implicit none
  real, intent(in)    :: tSimulation
  real, intent(in)    :: Xyz_D(nDim)    ! Position vector
  real, intent(in)    :: B1_D(nDim)     ! Magnetic field perturbation
  real, intent(in)    :: B0_D(nDim)     ! Magnetic field of planet
  real, intent(out)   :: u_D(nDim)      ! Velocity vector (output)

  real, parameter :: Epsilon = 0.01 ! Perturbation of X, Y or Z

  real :: XyzEpsilon_D(nDim)       ! Points shifted by Epsilon
  real :: XyzIono_D(nDim)          ! Mapped point on the ionosphere
  real :: Theta, Phi               ! Mapped point colatitude, longitude
  real :: ThetaNorm, PhiNorm       ! Normalized colatitude, longitude
  real :: Dist1, Dist2             ! Distance from ionosphere grid point

  real :: Potential_DI(nDim, 2)    ! Potential at the shifted positions
  real :: eField_D(nDim)           ! Electric field
  real :: b_D(nDim)                ! Magnetic field
  real :: B2                       ! Magnetic field squared

  integer :: iDim, iSide, iTheta, iPhi, iHemisphere

  character(len=*), parameter :: NameSub = 'calc_inner_bc_velocity1'
  logical :: DoTest, DoTestMe
  real :: tSimulationLast=-1.0
  real, save :: SmgGm_DD(nDim,nDim)
  !-------------------------------------------------------------------------

  call set_oktest(NameSub, DoTest, DoTestMe)

  if(DoTestMe)write(*,*)NameSub,' Xyz_D=',Xyz_D

  ! Calculate conversion matrix between GM and SMG coordinates
  if( tSimulationLast /= tSimulation ) then
     tSimulationLast = tSimulation
     SmgGm_DD = transform_matrix(tSimulation, TypeCoordSystem, 'SMG')
  end if

  ! Map points to obtain potential
  do iSide = 1, 2
     do iDim = 1, nDim

        ! Perturb the iDim coordinate
        XyzEpsilon_D = Xyz_D
        if(iSide == 1)then
           XyzEpsilon_D(idim) = XyzEpsilon_D(idim) - Epsilon
        else
           XyzEpsilon_D(idim) = XyzEpsilon_D(idim) + Epsilon
        end if

        ! Transform into SMG coordinates
        XyzEpsilon_D = matmul(SmgGm_DD, XyzEpsilon_D)

        ! Map down to the ionosphere at radius rIonosphere
        call map_planet_field(tSimulation, XyzEpsilon_D, 'SMG NORM', &
             rIonosphere, XyzIono_D, iHemisphere)

        ! Calculate angular coordinates
        call xyz_to_dir(XyzIono_D, Theta, Phi)

        ! Interpolate potential

        ! Get normalized coordinates
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

        Potential_DI(iDim, iSide) = &
             (1 - Dist1)*( (1-Dist2) * IonoPotential_II(iTheta  , iPhi  )  &
             +             Dist2     * IonoPotential_II(iTheta  , iPhi+1)) &
             + Dist1    *( (1-Dist2) * IonoPotential_II(iTheta+1, iPhi  )  &
             +             Dist2     * IonoPotential_II(iTheta+1, iPhi+1))

        if(DoTestMe)then
           write(*,*)NameSub,' iDim, iSide  =',iDim, iSide
           write(*,*)NameSub,' XyzEpsilon_D =',XyzEpsilon_D
           write(*,*)NameSub,' XyzIono_D    =',XyzIono_D
           write(*,*)NameSub,' Theta, Phi   =',Theta,Phi
           write(*,*)NameSub,' iTheta, iPhi =',iTheta,iPhi
           write(*,*)NameSub,' Dist1, Dist2 =',Dist1,Dist2
           write(*,*)NameSub,' Potential_DI =',Potential_DI(iDim, iSide)
        end if

     end do
  end do

  b_D = B1_D + B0_D
  B2  = sum(b_D**2)
  
  ! E = -grad(Potential)
  eField_D = - (Potential_DI(:,2) - Potential_DI(:,1))/(2*Epsilon)

  ! U = (E x B) / B^2
  u_D = cross_product(eField_D, b_D) / B2

  if(DoTestMe)then
     write(*,*)NameSub,' b_D=',b_D
     write(*,*)NameSub,' E_D=',eField_D
     write(*,*)NameSub,' u_D=',u_D
  endif

  ! Subtract the radial component of the velocity
  u_D = u_D - Xyz_D * sum(Xyz_D * u_D) / sum(Xyz_D**2)

  if(DoTestMe)then
     write(*,*)NameSub,' Final u_D=',u_D
  end if

end subroutine calc_inner_bc_velocity1

!=============================================================================

real function logvar_ionosphere(NameLogvar)
  use ModProcMH,  ONLY: iProc
  use ModIO,      ONLY: write_myname
  use ModIonoPotential, ONLY: nThetaIono, IonoPotential_II
  implicit none
  character (len=*), intent(in) :: NameLogvar
  integer :: nWarn = 0 ! warn multiple times to catch multiple log variables
  !---------------------------------------------------------------------------
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
