!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

!============================================================================

module GM_couple_ie

  use ModProcMH

  implicit none
  save

  private ! except

  public:: GM_get_for_ie
  public:: GM_put_from_ie
  public:: GM_put_mag_from_ie
  public:: print_iono_potential
  public:: clean_mod_field_aligned_current
  public:: get_ie_grid_index

  ! Joule heating
  real, public, allocatable :: IonoJouleHeating_II(:,:)

  ! Ionosphere potential
  real, allocatable, public :: IonoPotential_II(:,:)
  real, allocatable, public :: dIonoPotential_DII(:,:,:)

  ! Ionosphere grid description
  integer, public   :: nThetaIono = -1, nPhiIono = -1
  real, public      :: rIonosphere

  ! Local variables
  real, allocatable :: ThetaIono_I(:), PhiIono_I(:)
  real              :: dThetaIono, dPhiIono

  ! Field aligned current
  real, allocatable :: FieldAlignedCurrent_II(:,:), &
       bCurrentLocal_VII(:,:,:), bCurrent_VII(:,:,:)

contains

  !========================================================================== 
  subroutine print_iono_potential

    use ModIoUnit, ONLY: UnitTmp_

    integer:: i, j

    do j=1,nPhiIono
       do i=1,nThetaIono
          write(UNITTMP_,'(2i4,3G14.6)')i,j,&
               ThetaIono_I(i),PhiIono_I(j),IonoPotential_II(i,j)
       end do
    end do

  end subroutine print_iono_potential

  !==========================================================================
  subroutine init_mod_ie_grid(iSize, jSize)
    ! The ionosphere works on two hemispheres with a node based grid
    ! iSize is the number of latitude nodes from the pole to the equator.
    ! jSize is the number of longitude nodes (including a periodic overlap)

    use ModNumConst, ONLY: cPi, cTwoPi
    use ModPhysics,  ONLY: Si2No_V, UnitX_
    use CON_coupler, ONLY: Grid_C, IE_

    integer, intent(in) :: iSize, jSize
    character(len=*), parameter :: NameSub='init_mod_ie_grid'
    !-------------------------------------------------------------------------

    if(nThetaIono > 0) RETURN

    nThetaIono = Grid_C(IE_) % nCoord_D(1)
    nPhiIono   = Grid_C(IE_) % nCoord_D(2)

    if(nThetaIono /= iSize .or. nPhiIono /= jSize)then
       write(*,*)NameSub,': Grid_C(IE_)%nCoord_D(1:2)=',&
            Grid_C(IE_) % nCoord_D(1:2)
       write(*,*)NameSub,': iSize,2*iSize-1,jSize=',iSize,2*iSize-1,jSize
       call stop_mpi(NameSub//' ERROR: Inconsistent IE grid sizes')
    endif

    allocate(ThetaIono_I(nThetaIono), PhiIono_I(nPhiIono))
    ThetaIono_I = Grid_C(IE_) % Coord1_I
    PhiIono_I   = Grid_C(IE_) % Coord2_I
    rIonosphere = Grid_C(IE_) % Coord3_I(1) * Si2No_V(UnitX_)

    dThetaIono = cPi    / (nThetaIono - 1)
    dPhiIono   = cTwoPi / (nPhiIono - 1)

  end subroutine init_mod_ie_grid
  !===========================================================================
  subroutine get_ie_grid_index(Theta, Phi, ThetaNorm, PhiNorm)
    real, intent(in) :: Theta, Phi
    real, intent(out):: ThetaNorm, PhiNorm

    ThetaNorm = Theta / dThetaIono
    PhiNorm   = Phi   / dPhiIono
  end subroutine get_ie_grid_index

  !============================================================================

  subroutine init_mod_field_aligned_current(iSize,jSize)
    integer, intent(in) :: iSize, jSize

    if(allocated(FieldAlignedCurrent_II)) RETURN

    call init_mod_ie_grid(iSize, jSize)

    allocate( bCurrent_VII(0:6, nThetaIono, nPhiIono), &
         bCurrentLocal_VII(0:6, nThetaIono, nPhiIono), &
         FieldAlignedCurrent_II(nThetaIono, nPhiIono))

  end subroutine init_mod_field_aligned_current

  !============================================================================
  subroutine clean_mod_field_aligned_current

    if(allocated(FieldAlignedCurrent_II)) &
         deallocate( bCurrent_VII, bCurrentLocal_VII, FieldAlignedCurrent_II)
  end subroutine clean_mod_field_aligned_current

  !============================================================================

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

  !============================================================================

  subroutine init_mod_iono_jouleheating(iSize, jSize)

    integer, intent(in) :: iSize, jSize
    character(len=*), parameter :: NameSub='init_mod_iono_jouleheating'
    !-------------------------------------------------------------------------

    if(allocated(IonoJouleHeating_II)) RETURN

    call init_mod_ie_grid(iSize, jSize)

    allocate( IonoJouleHeating_II(nThetaIono, nPhiIono))

  end subroutine init_mod_iono_jouleheating
  !============================================================================

  !============================================================================

  subroutine GM_get_for_ie(Buffer_IIV,iSize,jSize,nVar)

    use CON_axes, ONLY: transform_matrix
    use ModMain, ONLY: Time_Simulation, TypeCoordSystem
    use ModCurrent,            ONLY: calc_field_aligned_current
    use ModRaytrace, ONLY: DoTraceIE, RayResult_VII, RayIntegral_VII, &
         InvB_, RhoInvB_, pInvB_, xEnd_, yEnd_, zEnd_, CLOSEDRAY, GmSm_DD
    use ModNumConst, ONLY: cRadToDeg
    use ModProcMH,         ONLY: iProc
    use ModPhysics, ONLY: No2Si_V, UnitX_, UnitP_, UnitRho_, UnitB_, UnitJ_
    use ModCoordTransform, ONLY: sph_to_xyz, xyz_to_sph

    character (len=*), parameter :: NameSub='GM_get_for_ie'

    integer, intent(in) :: iSize, jSize, nVar
    real, intent(out), dimension(iSize, jSize, nVar) :: Buffer_IIV
    integer :: i, j
    real :: Radius, Phi, Theta, LatBoundary
    real, allocatable, dimension(:), save :: IE_lat, IE_lon
    real, allocatable, dimension(:,:,:)   :: bSm_DII 
    real :: XyzIono_D(3), RtpIono_D(3), Lat,Lon, dLat,dLon
    logical :: DoTest, DoTestMe
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest, DoTestMe)
    if(DoTest)write(*,*)NameSub,': starting'

    if(.not.allocated(FieldAlignedCurrent_II)) &
         call init_mod_field_aligned_current(iSize, jSize)

    if(.not.allocated(IE_lat)) &
         allocate(IE_lat(iSize), IE_lon(jSize))

    if(.not. allocated(bSm_DII))&
         allocate(bSm_DII(3,iSize,jSize))

    ! initialize all elements to zero on proc 0, others should not use it
    if(iProc == 0) Buffer_IIV = 0.0

    ! Put field aligned currents into the first "slice" of the buffer
    call calc_field_aligned_current(nThetaIono, nPhiIono, rIonosphere, &
         FieldAlignedCurrent_II, bSm_DII, LatBoundary, ThetaIono_I, PhiIono_I)

    ! Take radial component of FAC and put it into the buffer sent to IE
    ! The resulting FAC_r will be positive for radially outgoing current
    ! and negative for radially inward going currents.
    if(iProc==0)then
       do j=1, nPhiIono
          Phi = PhiIono_I(j)

          do i=1, nThetaIono
             Theta = ThetaIono_I(i)
             call sph_to_xyz(rIonosphere, Theta, Phi, XyzIono_D)
             FieldAlignedCurrent_II(i,j) = FieldAlignedCurrent_II(i,j) &
                  * sum(bSm_DII(:,i,j)*XyzIono_D) &
                  / (sqrt(sum((bSm_DII(:,i,j))**2)) * rIonosphere)
          end do
       end do

       ! Save the latitude boundary information to the equator 
       Buffer_IIV(:,:,1) = FieldAlignedCurrent_II(:,:)*No2Si_V(UnitJ_)
       Buffer_IIV(nThetaIono/2:nThetaIono/2+1,1,1) = LatBoundary                          

    end if

    deallocate(bSm_DII)

    if(DoTraceIE) then
       ! Load grid and convert to lat-lon in degrees
       IE_lat = 90.0 - cRadToDeg * ThetaIono_I(1:iSize)
       IE_lon =        cRadToDeg * PhiIono_I
       Radius = (6378.+100.)/6378.
       call integrate_ray_accurate(iSize, jSize, IE_lat, IE_lon, Radius, &
            'InvB,RhoInvB,pInvB,Z0x,Z0y,Z0b')

       ! Only processor 0 has the resulting integrals,
       !   the others do not participate in the coupling
       if(iProc == 0)then
          where(RayResult_VII(InvB_,:,:) > 0.)
             RayResult_VII(RhoInvB_,:,:) = RayResult_VII(RhoInvB_,:,:) &
                  /RayResult_VII(InvB_,:,:)
             RayResult_VII(pInvB_,:,:)   = RayResult_VII(pInvB_,:,:) &
                  /RayResult_VII(InvB_,:,:)
          end where
          where(RayResult_VII(xEnd_,:,:) <= CLOSEDRAY)
             RayResult_VII(   InvB_,:,:) = -1.
             RayResult_VII(rhoInvB_,:,:) = 0.
             RayResult_VII(  pInvB_,:,:) = 0.
          end where
          Buffer_IIV(:,:,2) = RayResult_VII(   InvB_,:,:) &
               * No2Si_V(UnitX_)/No2Si_V(UnitB_)
          Buffer_IIV(:,:,3) = RayResult_VII(RhoInvB_,:,:) * No2Si_V(UnitRho_)
          Buffer_IIV(:,:,4) = RayResult_VII(  pInvB_,:,:) * No2Si_V(UnitP_)

          ! Transformation matrix from default (GM) to SM coordinates
          GmSm_DD = transform_matrix(time_simulation,TypeCoordSystem,'SMG')

          ! Loop to compute deltas
          do i=1,iSize
             do j=1,jSize
                Lat = -IE_Lat(i)
                Lon =  IE_Lon(j)
                if(RayResult_VII(InvB_,i,j)>1.e-10)then
                   XyzIono_D(1)=RayResult_VII(xEnd_,i,j)
                   XyzIono_D(2)=RayResult_VII(yEnd_,i,j)
                   XyzIono_D(3)=RayResult_VII(zEnd_,i,j)
                   XyzIono_D = matmul(GmSm_DD,XyzIono_D)
                   call xyz_to_sph(XyzIono_D, RtpIono_D)
                   Lat=90.-RtpIono_D(2)*cRadToDeg
                   Lon=    RtpIono_D(3)*cRadToDeg
                end if

                ! Shift to closer to local value, even if outside 0-360
                if( (IE_Lon(j)-Lon)<-180. ) Lon = Lon-360.
                if( (IE_Lon(j)-Lon)> 180. ) Lon = Lon+360.

                ! Compute deltas
                dLat = abs(Lat) - abs(IE_Lat(i))
                dLon =     Lon  -     IE_Lon(j)

                ! Put into exchange buffer
                Buffer_IIV(i,j,5) = dLat
                Buffer_IIV(i,j,6) = dLon
             end do
          end do
       end if
       deallocate(RayIntegral_VII, RayResult_VII)
    end if

    if(DoTest)write(*,*)NameSub,': finished'
  end subroutine GM_get_for_ie

  !============================================================================
  subroutine GM_put_from_ie(Buffer_IIV,iSize,jSize)

    use ModPhysics,       ONLY: Si2No_V, UnitX_, UnitElectric_, UnitPoynting_
    use ModProcMH

    character(len=*), parameter :: NameSub='GM_put_from_ie'

    integer, intent(in) :: iSize,jSize
    integer, parameter  :: nVar = 2
    real, intent(in) :: Buffer_IIV(iSize,jSize,nVar)
    !  character(len=*), intent(in) :: NameVar

    logical :: DoTest, DoTestMe
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)
    !  if(DoTest)write(*,*)NameSub,': NameVar,iSize,jSize=',NameVar,iSize,jSize

    if(.not. allocated(IonoPotential_II)) &
         call init_mod_iono_potential(iSize,jSize)

    IonoPotential_II = Buffer_IIV(:,:,1)*Si2No_V(UnitElectric_)*Si2No_V(UnitX_)
    call calc_grad_iono_potential

    if (.not. allocated(IonoJouleHeating_II)) &
         call init_mod_iono_jouleheating(iSize,jSize)

    ! Add the iono. jouleheating - Yiqun Sep 2008
    IonoJouleHeating_II = Buffer_IIV(:,:,2) * Si2No_V(UnitPoynting_)
    call map_jouleheating_to_inner_bc

    if(DoTest)write(*,*)NameSub,': done'

  end subroutine GM_put_from_ie

  !==========================================================================
  subroutine GM_put_mag_from_ie(Buffer_DII, iSize)
    ! Get magnetometer "measurements" from IE.  They are received via input
    ! array Buffer_DII, which has dimensions (3x2xiSize) where iSize should
    ! equal the total number of shared magnetometers.

    use ModGroundMagPerturb, ONLY: nMagnetometer, IeMagPerturb_DII
    use ModGmGeoindices,     ONLY: nIndexMag, MagPerbIE_DI

    integer, intent(in) :: iSize
    real, intent(in)    :: Buffer_DII(3,2,iSize)

    character(len=*), parameter :: NameSub='GM_put_mag_from_ie'
    logical :: DoTest, DoTestMe
    !--------------------------------------------------------------------------
    call set_oktest(NameSub, DoTest, DoTestMe)

    if(DoTestMe)write(*,*)NameSub, ' nIndexMag, nMag, iSize = ', &
         nIndexMag, nMagnetometer, iSize

    if( (nIndexMag+nMagnetometer) /= iSize)call CON_stop(NameSub// &
         ' Number of shared magnetometers does not match!')
    
    ! Place geomagnetic index data into right place (combine hall+pederson):
    if (nIndexMag>0) &
         MagPerbIE_DI(:,:) = Buffer_DII(:,1,1:nIndexMag) &
         + Buffer_DII(:,2,1:nIndexMag)

    ! Place regular mag data into right place (keep hall/pederson separate):
    if (nMagnetometer>0) &
         IeMagPerturb_DII = Buffer_DII(:,:,nIndexMag+1:)

  end subroutine GM_put_mag_from_ie

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
    do i = 1, nThetaIono; do j =1, nPhiIono
       call sph_to_xyz(rIonosphere, ThetaIono_I(i), PhiIono_I(j), XyzIono_D)
       call get_planet_field(Time_Simulation,XyzIono_D, 'SMG NORM', bIono_D)
       bIono = sqrt(sum(bIono_D**2))

       ! map out to GM (caution!, not like map down to the ionosphere, 
       ! there is always a corresponding position.)
       call map_planet_field(Time_Simulation, XyzIono_D, 'SMG NORM', &
            rBody, Xyz_tmp, iHemisphere)

       if (iHemisphere == 0) then 
          ! not a mapping in the dipole, but to the equator
          ! assume not outflow to GM inner boundary
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
  subroutine calc_inner_bc_velocity1(tSimulation,Xyz_D,B1_D,B0_D,u_D)

    use ModMain,           ONLY: TypeCoordSystem, MaxDim
    use CON_axes,          ONLY: transform_matrix
    use ModCoordTransform, ONLY: xyz_to_dir, cross_product
    use CON_planet_field,  ONLY: map_planet_field

    real, intent(in)    :: tSimulation
    real, intent(in)    :: Xyz_D(MaxDim)    ! Position vector
    real, intent(in)    :: B1_D(MaxDim)     ! Magnetic field perturbation
    real, intent(in)    :: B0_D(MaxDim)     ! Magnetic field of planet
    real, intent(out)   :: u_D(MaxDim)      ! Velocity vector (output)

    real, parameter :: Epsilon = 0.01 ! Perturbation of X, Y or Z

    real :: XyzEpsilon_D(MaxDim)     ! Points shifted by Epsilon
    real :: XyzIono_D(MaxDim)        ! Mapped point on the ionosphere
    real :: Theta, Phi               ! Mapped point colatitude, longitude
    real :: ThetaNorm, PhiNorm       ! Normalized colatitude, longitude
    real :: Dist1, Dist2             ! Distance from ionosphere grid point

    real :: Potential_DI(MaxDim, 2)  ! Potential at the shifted positions
    real :: eField_D(MaxDim)         ! Electric field
    real :: b_D(MaxDim)              ! Magnetic field
    real :: B2                       ! Magnetic field squared

    integer :: iDim, iSide, iTheta, iPhi, iHemisphere

    character(len=*), parameter :: NameSub = 'calc_inner_bc_velocity1'
    logical :: DoTest, DoTestMe
    real :: tSimulationLast=-1.0
    real, save :: SmgGm_DD(MaxDim,MaxDim)
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
       do iDim = 1, MaxDim

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
  !============================================================================

end module GM_couple_ie

!==========================================================================

subroutine map_inner_bc_jouleheating(tSimulation, Xyz_D, JouleHeating)

  use GM_couple_ie, ONLY: get_ie_grid_index, IonoJouleHeating_II
  use ModCoordTransform, ONLY: xyz_to_dir
  use ModMain,           ONLY: TypeCoordSystem, MaxDim
  use CON_axes,          ONLY: transform_matrix

  !INPUT ARGUMENTS:
  real, intent(in)    :: tSimulation    ! Simulation time
  real, intent(in)    :: Xyz_D(MaxDim)    ! Position vector
  real, intent(out)   :: JouleHeating
  real :: Theta, Phi           ! Mapped point colatitude, longitude
  real :: ThetaNorm, PhiNorm   ! Normalized colatitude, longitude
  integer :: iTheta, iPhi
  character(len=*), parameter :: NameSub = 'map_inner_bc_jouleheating'
  logical :: DoTest, DoTestMe
  real :: Xyz_D_tmp(3)
  ! -----------------------------------------------------

  call set_oktest(NameSub, DoTest, DoTestMe)

  ! Convert Xyz_D into SMG coordinates
  Xyz_D_tmp = matmul(transform_matrix(tSimulation, TypeCoordSystem, 'SMG'), Xyz_D)

  ! Calculate angular coordinates
  call xyz_to_dir(Xyz_D_tmp, Theta, Phi)

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

  JouleHeating = IonoJouleHeating_II(iTheta, iPhi)

end subroutine map_inner_bc_jouleheating

