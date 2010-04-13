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
  real, allocatable :: &
       FieldAlignedCurrent_II(:,:), bCurrentLocal_VII(:,:,:), bCurrent_VII(:,:,:)

contains

  !============================================================================

  subroutine init_mod_field_aligned_current(iSize,jSize)
    use ModProcMH,         ONLY: iProc, iComm
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
 
end module ModFieldAlignedCurrent

!==============================================================================

subroutine GM_get_for_ie(Buffer_IIV,iSize,jSize,nVar)

  use CON_axes, ONLY: transform_matrix
  use ModMain, ONLY: Time_Simulation, TypeCoordSystem
  use ModIeGrid
  use ModCurrent,            ONLY: calc_field_aligned_current
  use ModFieldAlignedCurrent,ONLY: FieldAlignedCurrent_II, &
       init_mod_field_aligned_current!, calc_field_aligned_current
  use ModRaytrace, ONLY: DoTraceIE, RayResult_VII, RayIntegral_VII, &
       InvB_, RhoInvB_, pInvB_, xEnd_, yEnd_, zEnd_, CLOSEDRAY, GmSm_DD
  use ModNumConst, ONLY: cRadToDeg, cDegToRad
  use ModProcMH,         ONLY: iProc, iComm
  use ModPhysics, ONLY: No2Si_V, UnitP_, UnitRho_, UnitB_,rCurrents, UnitJ_
  use ModCoordTransform, ONLY: sph_to_xyz, xyz_to_sph
  implicit none

  character (len=*), parameter :: NameSub='GM_get_for_ie'

  integer, intent(in) :: iSize, jSize, nVar
  real, intent(out), dimension(iSize, jSize, nVar) :: Buffer_IIV
  integer :: i, j
  real :: Radius,iError, Phi, Theta, LatBoundary
  real, allocatable, dimension(:), save :: IE_lat, IE_lon
  real, allocatable, dimension(:,:,:)   :: bSm_DII 
  real :: XyzIono_D(3), RtpIono_D(3), Xyz_D(3), Lat,Lon, dLat,dLon
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
        where(RayResult_VII(InvB_,:,:)>0.)
           RayResult_VII(RhoInvB_,:,:)=RayResult_VII(RhoInvB_,:,:)/RayResult_VII(InvB_,:,:)
           RayResult_VII(pInvB_,:,:)=RayResult_VII(pInvB_,:,:)/RayResult_VII(InvB_,:,:)
        end where
        where(RayResult_VII(xEnd_,:,:) <= CLOSEDRAY)
           RayResult_VII(   InvB_,:,:) = -1.
           RayResult_VII(rhoInvB_,:,:) = 0.
           RayResult_VII(  pInvB_,:,:) = 0.
        end where
        Buffer_IIV(:,:,2) = RayResult_VII(   InvB_,:,:) / No2Si_V(UnitB_)
        Buffer_IIV(:,:,3) = RayResult_VII(rhoInvB_,:,:) * No2Si_V(UnitRho_)
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
