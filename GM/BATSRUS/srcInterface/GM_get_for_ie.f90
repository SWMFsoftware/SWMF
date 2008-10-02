!^CFG COPYRIGHT UM
!^CMP FILE IE
!============================================
!                                           |
!     Magnetosphere/Ionosphere Coupling     |
!     Subroutines here are used with MPI    |
!============================================
module ModFieldAlignedCurrent

  use ModIeGrid
  
  implicit none
  save
  real, allocatable :: bCurrentLocal_VII(:,:,:), bCurrent_VII(:,:,:), &
       FieldAlignedCurrent_II(:,:)
  real :: LatBoundary

contains

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

  subroutine calc_field_aligned_current

    use ModVarIndexes,     ONLY: Bx_, Bz_, nVar
    use ModMain,           ONLY: Time_Simulation, TypeCoordSystem, nBlock
    use ModPhysics,        ONLY: rCurrents, No2Si_V, Si2No_V, UnitB_, UnitJ_
    use ModCoordTransform, ONLY: sph_to_xyz
    use CON_planet_field,  ONLY: get_planet_field, map_planet_field
    use CON_axes,          ONLY: transform_matrix
    use ModProcMH,         ONLY: iProc, iComm
    use ModMpi
    use ModNumConst, ONLY: cHalfPi

    ! Map the grid points from the ionosphere radius to rCurrents.
    ! Calculate the field aligned currents there, use the ratio of the
    ! magnetic field strength, and project to the radial direction.
    ! The result is saved into

    integer :: i, j, iHemisphere, iError
    real    :: Phi, Theta, XyzIono_D(3), Xyz_D(3)
    real    :: bIono_D(3), bIono, B0_D(3), b_D(3), b, j_D(3), Fac
    real    :: GmSmg_DD(3,3)
    real    :: State_V(Bx_-1:nVar+3)
    !-------------------------------------------------------------------------

    GmSmg_DD = transform_matrix(Time_Simulation, 'SMG', TypeCoordSystem)

    LatBoundary = 100.0
    do j = 1, nPhiIono
       Phi = PhiIono_I(j)
       do i = 1, nThetaIono

          Theta = ThetaIono_I(i)
          call sph_to_xyz(rIonosphere,Theta,Phi, XyzIono_D)

          call map_planet_field(Time_Simulation, XyzIono_D, 'SMG NORM', &
               rCurrents, Xyz_D, iHemisphere)

          if(iHemisphere == 0) then
             ! Assign weight 1, magnetic field of 1,0,0 and current 0,0,0
             bCurrentLocal_VII(:,i,j) = (/1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
             CYCLE
          end if

          LatBoundary = min( abs(Theta - cHalfPi), LatBoundary )

          ! Get the B0 field at the mapped position
          call get_planet_field(Time_Simulation, Xyz_D,'SMG NORM',B0_D)
          B0_D = B0_D*Si2No_V(UnitB_)

          ! Convert to GM coordinates
          Xyz_D = matmul(GmSmg_DD, Xyz_D)

          ! Extract currents and magnetic field for this position
          call get_point_data(0.0, Xyz_D, 1, nBlock, Bx_, nVar+3, State_V)

          bCurrentLocal_VII(0  ,i,j) = State_V(Bx_-1)        ! Weight
          bCurrentLocal_VII(1:3,i,j) = State_V(Bx_:Bz_) + &  ! B1 and B0
               State_V(Bx_-1)*B0_D
          bCurrentLocal_VII(4:6,i,j) = State_V(nVar+1:nVar+3) ! Currents

          if(.false. .and. i==6 .and. j==6)then
             write(*,*)'iHemispher=',iHemisphere
             write(*,*)'Phi,Theta=',Phi,Theta
             write(*,*)'XyzIono_D=',XyzIono_D
             write(*,*)'Xyz_D    =',Xyz_D
             write(*,*)'rCurrents=',rCurrents,sqrt(sum(Xyz_D**2))
             write(*,*)'b0_D     =',b0_D
             write(*,*)'bCurrentLocal_VII =',bCurrentLocal_VII(:,i,j)
             call stop_mpi('DEBUG')
          end if
       end do
    end do

    ! Add up contributions from all PE-s to processor 0
    call MPI_reduce(bCurrentLocal_VII, bCurrent_VII, nThetaIono*nPhiIono*7,&
         MPI_REAL, MPI_SUM, 0, iComm,iError)

    if(iProc /= 0) RETURN

    ! Calculate the field aligned current
    do j = 1, nPhiIono
       Phi = PhiIono_I(j)
       do i = 1, nThetaIono
          Theta = ThetaIono_I(i)

          ! Calculate magnetic field strength at the ionosphere grid point
          call sph_to_xyz(rIonosphere, Theta, Phi, XyzIono_D)
          call get_planet_field(Time_Simulation, XyzIono_D,'SMG NORM',bIono_D)

          ! Convert to GM units and get magnitude
          bIono_D = bIono_D*Si2No_V(UnitB_)
          bIono   = sqrt(sum(bIono_D**2))

          ! Divide MHD values by the total weight if it exceeds 1.0
          if(bCurrent_VII(0,i,j) > 1.0) bCurrent_VII(:,i,j) = &
               bCurrent_VII(:,i,j) / bCurrent_VII(0,i,j)

          ! Extract magnetic field and current
          b_D = bCurrent_VII(1:3,i,j)
          j_D = bCurrent_VII(4:6,i,j)


          ! The strength of the field
          b = sqrt(sum(b_D**2))


          ! Convert b_D into a unit vector
          b_D = b_D / b

          ! Calculate field aligned current
          Fac = sum(b_D*j_D)

          ! Multiply by the ratio of the magnetic field strengths
          Fac = bIono / b * Fac

          ! Take the radial component
          Fac = Fac * sum(b_D*XyzIono_D) / rIonosphere

          ! Store the result in SI units
          FieldAlignedCurrent_II(i,j) = Fac * No2Si_V(UnitJ_)

       end do
    end do

    ! Save the latitude boundary information to the equator
    FieldAlignedCurrent_II(nThetaIono/2:nThetaIono/2+1,1) = LatBoundary

  end subroutine calc_field_aligned_current

end module ModFieldAlignedCurrent

!==============================================================================

subroutine GM_get_for_ie(Buffer_IIV,iSize,jSize,nVar)

  use ModProcMH, ONLY: iProc
  use ModIeGrid
  use ModFieldAlignedCurrent,ONLY: FieldAlignedCurrent_II, &
       init_mod_field_aligned_current, calc_field_aligned_current
  use ModRaytrace, ONLY: DoTraceIE, RayResult_VII, RayIntegral_VII, &
       InvB_, RhoInvB_, pInvB_, xEnd_, CLOSEDRAY
  use ModNumConst, ONLY: cRadToDeg
  use ModPhysics, ONLY: No2Si_V, UnitP_, UnitRho_, UnitB_

  implicit none

  character (len=*), parameter :: NameSub='GM_get_for_ie'

  integer, intent(in) :: iSize, jSize, nVar
  real, intent(out), dimension(iSize, jSize, nVar) :: Buffer_IIV

  real :: Radius
  real, allocatable, dimension(:), save :: IE_lat, IE_lon

  logical :: DoTest, DoTestMe
  !--------------------------------------------------------------------------
  call CON_set_do_test(NameSub,DoTest, DoTestMe)
  if(DoTest)write(*,*)NameSub,': starting'

  if(.not.allocated(FieldAlignedCurrent_II)) &
       call init_mod_field_aligned_current(iSize, jSize)

  if(.not.allocated(IE_lat)) &
       allocate(IE_lat(iSize), IE_lon(jSize))

  ! initialize all elements to zero on proc 0, others should not use it
  if(iProc == 0) Buffer_IIV = 0.0

  ! Put field aligned currents into the first "slice" of the buffer
  call calc_field_aligned_current
  if(iProc == 0)Buffer_IIV(:,:,1) = FieldAlignedCurrent_II(:,:)

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
        where(RayResult_VII(InvB_,:,:)>0.)
           RayResult_VII(RhoInvB_,:,:)=RayResult_VII(RhoInvB_,:,:)/RayResult_VII(InvB_,:,:)
           RayResult_VII(pInvB_,:,:)=RayResult_VII(pInvB_,:,:)/RayResult_VII(InvB_,:,:)
        end where
        where(RayResult_VII(xEnd_,:,:) <= CLOSEDRAY)
           RayResult_VII(   InvB_,:,:) = -1.*No2Si_V(UnitB_)
           RayResult_VII(rhoInvB_,:,:) = 0.
           RayResult_VII(  pInvB_,:,:) = 0.
        end where
        Buffer_IIV(:,:,2) = RayResult_VII(   InvB_,:,:) / No2Si_V(UnitB_)
        Buffer_IIV(:,:,3) = RayResult_VII(rhoInvB_,:,:) * No2Si_V(UnitRho_)
        Buffer_IIV(:,:,4) = RayResult_VII(  pInvB_,:,:) * No2Si_V(UnitP_)
     end if
     deallocate(RayIntegral_VII, RayResult_VII)
  end if

  if(DoTest)write(*,*)NameSub,': finished'

end subroutine GM_get_for_ie
