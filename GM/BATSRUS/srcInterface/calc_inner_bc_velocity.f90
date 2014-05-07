!==========================================================================
!BOP
!ROUTINE: calc_inner_bc_velocity - calculate velocity at the inner boundary
!INTERFACE:
		  subroutine calc_inner_bc_velocity(tSimulation, Xyz_D, b_D, u_D)

  !USES:
  use ModIonoPotential
  use ModMain,           ONLY: TypeCoordSystem, MaxDim
  use ModCoordTransform, ONLY: xyz_to_dir, cross_product
  use CON_planet_field,  ONLY: map_planet_field

  implicit none

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

  character(len=*), parameter :: NameSub = 'calc_inner_bc_velocity'
  logical :: DoTestMe = .false.
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

