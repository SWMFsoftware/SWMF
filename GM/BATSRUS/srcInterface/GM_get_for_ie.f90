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
  real, allocatable ::FieldAlignedCurrentLocal_II(:,:), &
       FieldAlignedCurrent_II(:,:)!, bCurrentLocal_VII(:,:,:), bCurrent_VII(:,:,:)

  real :: LatBoundary

contains

  !============================================================================

  subroutine init_mod_field_aligned_current(iSize,jSize)
    use ModProcMH,         ONLY: iProc, iComm
    integer, intent(in) :: iSize, jSize

    if(allocated(FieldAlignedCurrent_II)) RETURN

    call init_mod_ie_grid(iSize, jSize)

    allocate( FieldAlignedCurrentLocal_II(nThetaIono, nPhiIono), &
         FieldAlignedCurrent_II(nThetaIono, nPhiIono))

  end subroutine init_mod_field_aligned_current

  !============================================================================
  subroutine clean_mod_field_aligned_current

    if(allocated(FieldAlignedCurrent_II)) &
         deallocate( FieldAlignedCurrentLocal_II, FieldAlignedCurrent_II)
  end subroutine clean_mod_field_aligned_current
 
end module ModFieldAlignedCurrent

!==============================================================================

subroutine GM_get_for_ie(Buffer_IIV,iSize,jSize,nVar)

  use ModIeGrid
  use ModCurrent,            ONLY: calc_field_aligned_current
  use ModFieldAlignedCurrent,ONLY: FieldAlignedCurrent_II, &
       init_mod_field_aligned_current, FieldALignedCurrentLocal_II
  use ModRaytrace, ONLY: DoTraceIE, RayResult_VII, RayIntegral_VII, &
       InvB_, RhoInvB_, pInvB_, xEnd_, CLOSEDRAY
  use ModNumConst, ONLY: cRadToDeg
  use ModProcMH,         ONLY: iProc, iComm
  use ModPhysics, ONLY: No2Si_V, UnitP_, UnitRho_, UnitB_,rCurrents, UnitJ_
  use ModCoordTransform, ONLY: sph_to_xyz
  implicit none

  character (len=*), parameter :: NameSub='GM_get_for_ie'

  integer, intent(in) :: iSize, jSize, nVar
  real, intent(out), dimension(iSize, jSize, nVar) :: Buffer_IIV
  integer :: i, j
  real :: Radius,iError, Phi, Theta
  real, allocatable, dimension(:), save :: IE_lat, IE_lon
  real, allocatable, dimension(:,:,:)   :: B_DII 
  real :: XyzIono_D(3)
  logical :: DoTest, DoTestMe
  !--------------------------------------------------------------------------
  call CON_set_do_test(NameSub,DoTest, DoTestMe)
  if(DoTest)write(*,*)NameSub,': starting'

  if(.not.allocated(FieldAlignedCurrent_II)) &
       call init_mod_field_aligned_current(iSize, jSize)

  if(.not.allocated(IE_lat)) &
       allocate(IE_lat(iSize), IE_lon(jSize))

 if(.not. allocated(B_DII))&
      allocate(B_DII(3,iSize,jSize))

  ! initialize all elements to zero on proc 0, others should not use it
  if(iProc == 0) Buffer_IIV = 0.0

  ! Put field aligned currents into the first "slice" of the buffer
  call calc_field_aligned_current(iSize,jSize,rIonosphere,FieldAlignedCurrent_II,B_DII)

  ! Convert Fac into radial direction
  if(iProc==0)then
     do j=1, nPhiIono
         Phi = PhiIono_I(j)

        do i=1, nThetaIono
           Theta = ThetaIono_I(i)
           ! the Latboundary point
           if (j==1 .and. i== nThetaIono/2 .or. j==1 .and. i==nThetaIono/2+1)CYCLE

           call sph_to_xyz(rIonosphere, Theta, Phi, XyzIono_D)
           ! in radial direction
           FieldAlignedCurrent_II(i,j) = FieldAlignedCurrent_II(i,j) &
                * sum(B_DII(:,i,j)*XyzIono_D) / rIonosphere /sqrt(sum((B_DII(:,i,j))**2))
        end do
     end do
  end if

  if(iProc == 0)Buffer_IIV(:,:,1) = FieldAlignedCurrent_II(:,:)*No2Si_V(UnitJ_)

  deallocate(B_DII)

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
