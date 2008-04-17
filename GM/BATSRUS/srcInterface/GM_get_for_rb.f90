!^CFG COPYRIGHT UM
!^CMP FILE RB

!==========================================================================
subroutine GM_get_for_rb_trace(iSizeIn, jSizeIn, NameVar, nVarLine, nPointLine)

  ! Do ray tracing for RB. 
  ! Provide total number of points along rays 
  ! and the number of variables to pass to RB
  use ModGmRbCoupling, ONLY: allocate_gm_rb, RB_lat, RB_lon
  use ModRayTrace, ONLY: DoExtractUnitSi

  use CON_line_extract, ONLY: line_get
  implicit none
  integer, intent(in)           :: iSizeIn, jSizeIn
  character (len=*), intent(in) :: NameVar
  integer, intent(out)          :: nVarLine, nPointLine
  real :: Radius

  character (len=*), parameter :: NameSub='GM_get_for_rb_trace'
  !---------------------------------------------------------------------

  if(NameVar /= 'Z0x:Z0y:Z0b:I_I:S_I:R_I:B_I:IMF') &
       call CON_stop(NameSub//' invalid NameVar='//NameVar)

  ! Allocate arrays
  call allocate_gm_rb(iSizeIn, jSizeIn)

  ! The RB ionosphere radius in normalized units
  Radius = (6378.+100.)/6378.  !!! could be derived from Grid_C ?
  DoExtractUnitSi = .true.
  call integrate_ray_accurate(iSizeIn, jSizeIn, RB_lat, RB_lon, Radius, &
       NameVar)

  call line_get(nVarLine, nPointLine)

  nVarLine = 4 ! We only pass line index, length, B and radial distance to RB

end subroutine GM_get_for_rb_trace

!==========================================================================

subroutine GM_get_for_rb(Buffer_IIV, iSizeIn, jSizeIn, nVar, &
     BufferLine_VI, nVarLine, nPointLine, NameVar)

  !call stop_mpi('RAYTRACE is OFF') !^CFG UNCOMMENT IF NOT RAYTRACE
  !^CFG IF RAYTRACE BEGIN

  use ModProcMH, ONLY: iProc

  use ModMain, ONLY: Time_Simulation

  use ModGmRbCoupling, ONLY: &
       RB_lat, RB_lon, &
       write_integrated_data_tec, write_integrated_data_idl, &
       MHD_SUM_vol, MHD_Xeq, MHD_Yeq, MHD_Beq, MHD_SUM_rho, MHD_SUM_p, NoValue

  use ModRaytrace, ONLY: RayResult_VII, RayIntegral_VII, &
       InvB_, Z0x_, Z0y_, Z0b_, RhoInvB_, pInvB_, xEnd_, CLOSEDRAY

  use ModVarIndexes, ONLY: Bx_, Bz_, MassFluid_I, IonFirst_

  use ModPhysics, ONLY: No2Si_V, UnitN_, UnitU_, UnitB_, UnitP_

  use CON_line_extract, ONLY: line_get, line_clean

  implicit none

  character (len=*), parameter :: NameSub='GM_get_for_rb'

  integer, intent(in)                                :: iSizeIn, jSizeIn, nVar
  real, intent(out), dimension(iSizeIn,jSizeIn,nVar) :: Buffer_IIV

  integer, intent(in) :: nPointLine, nVarLine
  real, intent(out)   :: BufferLine_VI(nVarLine, nPointLine)
  character (len=*), intent(in):: NameVar

  integer :: nVarExtract, nPoint, iPoint
  real, allocatable :: Buffer_VI(:,:)
  real :: Rho, Ux, Uy, Uz, Bx, By, Bz, p

  logical :: DoTestTec, DoTestIdl
  logical :: DoTest, DoTestMe

  integer :: iLat,iLon,iLine
  !--------------------------------------------------------------------------

  if(NameVar /= 'Z0x:Z0y:Z0b:I_I:S_I:R_I:B_I:IMF') &
       call CON_stop(NameSub//' invalid NameVar='//NameVar)

  if(iProc /= 0)then
     ! Clean and return
     deallocate(RayIntegral_VII, RayResult_VII)
     call line_clean
     RETURN
  end if

  call CON_set_do_test(NameSub//'_tec', DoTestTec, DoTestMe)
  call CON_set_do_test(NameSub//'_idl', DoTestIdl, DoTestMe)
  call CON_set_do_test(NameSub, DoTest, DoTestMe)

  ! Copy RayResult into small arrays used in old algorithm
  MHD_SUM_vol = RayResult_VII(InvB_   ,:,:)
  MHD_Xeq     = RayResult_VII(Z0x_    ,:,:)
  MHD_Yeq     = RayResult_VII(Z0y_    ,:,:)
  MHD_Beq     = RayResult_VII(Z0b_    ,:,:)
  MHD_SUM_rho = RayResult_VII(RhoInvB_,:,:)
  MHD_SUM_p   = RayResult_VII(pInvB_  ,:,:)

  ! Put impossible values if ray was not found for a lat-lon grid cell
  where(RayResult_VII(xEnd_,:,:) <= CLOSEDRAY)
     MHD_Xeq     = NoValue
     MHD_Yeq     = NoValue
     MHD_SUM_vol = 0.0
     MHD_SUM_rho = 0.0
     MHD_SUM_p   = 0.0
     MHD_Beq     = NoValue
  end where

  ! Put the extracted data into BufferLine_VI

  call line_get(nVarExtract, nPoint)
  if(nPoint /= nPointLine)call stop_mpi(NameSub//': nPointLine error')
  if(nVarExtract < nVarLine)call stop_mpi(NameSub//': nVarLine error')
  allocate(Buffer_VI(0:nVarExtract, nPoint))
  call line_get(nVarExtract, nPoint, Buffer_VI, DoSort=.true.)
  
  do iPoint = 1, nPoint

     iLine =  Buffer_VI(0,iPoint)     ! line index
     iLat = mod(iLine-1, iSizeIn) + 1
     iLon = (iLine-1)/iSizeIn + 1

     ! exclude open field lines by setting impossible line index
     if(MHD_Xeq(iLat, iLon) == NoValue)then
        !write(*,*)'!!!iPoint,iSizeIn,iLine,iLat,iLon=',iPoint,iSizeIn,iLine,iLat,iLon
        iLine = -1
     endif

     BufferLine_VI(1,iPoint) = iLine
     BufferLine_VI(2,iPoint) = Buffer_VI(1,iPoint)                 ! Length
     BufferLine_VI(3,iPoint) = sqrt(sum(Buffer_VI(2:4,iPoint)**2)) ! |r|
     BufferLine_VI(4,iPoint) = &
          sqrt(sum(Buffer_VI(4+Bx_:4+Bz_,iPoint)**2))       ! |B|
  end do

!  write(*,*) 'iLat, iLon, RB_lat(iLat), RB_lon(iLon),MHD_Xeq(iLat,iLon),MHD_Yeq(iLat,iLon), RayResult_VII(Z0x_:Z0y_,iLat,iLon)'
!  do iLon=1,jSizeIn
!     do iLat=1,iSizeIn
!        write(*,'(2i3,6f13.5)') iLat, iLon, RB_lat(iLat), RB_lon(iLon), &
!             MHD_Xeq(iLat,iLon),MHD_Yeq(iLat,iLon), RayResult_VII(Z0x_:Z0y_,iLat,iLon)
!     enddo
!  enddo
     
  deallocate(RayIntegral_VII, RayResult_VII, Buffer_VI)
  call line_clean

  ! Output before processing
  if(DoTest .or. DoTestTec)call write_integrated_data_tec
  if(DoTest .or. DoTestIdl)call write_integrated_data_idl

  ! Put results into output buffer
  Buffer_IIV(:,:,1)  = MHD_Xeq
  Buffer_IIV(:,:,2)  = MHD_Yeq
  Buffer_IIV(:,:,3)  = MHD_Beq * No2Si_V(UnitB_)
  
!  write(*,*) 'iLat, iLon, RB_lat(iLat), RB_lon(iLon),MHD_Xeq(iLat,iLon),MHD_Yeq(iLat,iLon), Buffer_IIV(:,:,1:2)'
!  do iLon=1,jSizeIn
!     do iLat=1,iSizeIn
!        write(*,'(2i3,6f13.5)') iLat, iLon, RB_lat(iLat), RB_lon(iLon), &
!             MHD_Xeq(iLat,iLon),MHD_Yeq(iLat,iLon), Buffer_IIV(iLat,iLon,1:2)
!     enddo
!  enddo

  ! Send solar wind values in the array of the extra integral
  ! This is a temporary solution. RB should use MHD_SUM_rho and MHD_SUM_p

  call get_solar_wind_point(Time_Simulation, 0.0, 0.0, &
       Rho, Ux, Uy, Uz, Bx, By, Bz, p)

  Buffer_IIV(1,:,4) = Rho/MassFluid_I(IonFirst_)  * No2Si_V(UnitN_)
  Buffer_IIV(2,:,4) = Ux * No2Si_V(UnitU_)
  Buffer_IIV(3,:,4) = Uy * No2Si_V(UnitU_)
  Buffer_IIV(4,:,4) = Uz * No2Si_V(UnitU_)
  Buffer_IIV(5,:,4) = Bx * No2Si_V(UnitB_)
  Buffer_IIV(6,:,4) = By * No2Si_V(UnitB_)
  Buffer_IIV(7,:,4) = Bz * No2Si_V(UnitB_)
  Buffer_IIV(8,:,4) = p  * No2Si_V(UnitP_)

  !^CFG END RAYTRACE
end subroutine GM_get_for_rb
