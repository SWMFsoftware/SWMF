!^CFG COPYRIGHT UM
!^CMP FILE IM

!==========================================================================
subroutine GM_get_for_im_trace(nRadius, nLon, NameVar, nVarLine, nPointLine)

  ! Do ray tracing for IM/RAM_SCB. 
  ! Provide total number of points along rays 
  ! and the number of variables to pass to IM

  use CON_line_extract, ONLY: line_get
  use CON_coupler, ONLY: Grid_C, IM_

  implicit none
  integer, intent(in)           :: nRadius, nLon
  character (len=*), intent(in) :: NameVar
  integer, intent(out)          :: nVarLine, nPointLine

  real, allocatable, save :: RadiusIm_I(:), LongitudeIm_I(:)

  character(len=*), parameter :: NameSub='GM_get_for_im_trace'
  !---------------------------------------------------------------------

  if(NameVar /= 'l:x:y:z:rho:ux:uy:uz:bx:by:bz:p') &
       call CON_stop(NameSub//' invalid NameVar='//NameVar)

  if(.not.allocated(RadiusIm_I))then
     if(   nRadius /= Grid_C(IM_) % nCoord_D(1) .or. &
           nLon    /= Grid_C(IM_) % nCoord_D(2) )then
        write(*,*)NameSub//' grid sizes do not agree nRadius,nLon,nCells=',&
             nRadius, nLon, Grid_C(IM_) % nCoord_D(1:2)
        call CON_stop(NameSub//' ERROR')
     end if
     allocate(RadiusIm_I(nRadius), LongitudeIm_I(nLon))
     RadiusIm_I    = Grid_C(IM_) % Coord1_I
     LongitudeIm_I = Grid_C(IM_) % Coord2_I
  end if

  call trace_ray_equator(nRadius, nLon, RadiusIm_I, LongitudeIm_I, .true.)

  call line_get(nVarLine, nPointLine)

  ! We only pass index, length, 3 coordinates and 8 MHD variables
  nVarLine = 13 

end subroutine GM_get_for_im_trace

!==========================================================================

subroutine GM_get_for_im_line(BufferLine_VI, nVarLine, nPointLine, NameVar)

  !call stop_mpi('RAYTRACE is OFF') !^CFG UNCOMMENT IF NOT RAYTRACE
  !^CFG IF RAYTRACE BEGIN

  use ModProcMH,  ONLY: iProc

  use ModPhysics, ONLY: rBody
  use ModVarIndexes, ONLY: Rho_, Ux_, Uz_, Bx_, Bz_, p_

  use CON_line_extract, ONLY: line_get, line_clean

  implicit none

  character (len=*), parameter :: NameSub='GM_get_for_im_line'

  integer, intent(in) :: nPointLine, nVarLine
  real, intent(out)   :: BufferLine_VI(nVarLine, nPointLine)
  character (len=*), intent(in):: NameVar

  integer :: nVarExtract, nPoint, iPoint, iLine
  real, allocatable :: Buffer_VI(:,:)

  logical :: DoTest, DoTestMe
  !--------------------------------------------------------------------------

  if(NameVar /= 'l:x:y:z:rho:ux:uy:uz:bx:by:bz:p') &
       call CON_stop(NameSub//' invalid NameVar='//NameVar)

  if(iProc /= 0)then
     ! Clean and return
     call line_clean
     RETURN
  end if

  call CON_set_do_test(NameSub, DoTest, DoTestMe)

  ! Put the extracted data into BufferLine_VI

  call line_get(nVarExtract, nPoint)
  if(nPoint /= nPointLine)call stop_mpi(NameSub//': nPointLine error')
  if(nVarExtract+1 < nVarLine)then
     write(*,*)NameSub, ' nVarExtract, nVarLine=',nVarExtract, nVarLine
     call stop_mpi(NameSub//': nVarLine error')
  end if
  allocate(Buffer_VI(0:nVarExtract, nPoint))
  call line_get(nVarExtract, nPoint, Buffer_VI, DoSort=.true.)

  do iPoint = 1, nPoint

     iLine = Buffer_VI(0,iPoint)     ! line index

     ! exclude open field lines by setting impossible line index
     ! open field lines end somewhere far from the inner boundary at rBody
     if(maxval(abs(Buffer_VI(2:4,iPoint))) > 1.5*rBody) iLine = -1

     BufferLine_VI(1,iPoint)    = iLine
     BufferLine_VI(2:5,iPoint)  = Buffer_VI(1:4        ,iPoint) ! Length,X,Y,Z
     BufferLine_VI(6,iPoint)    = Buffer_VI(Rho_+4     ,iPoint) ! density
     BufferLine_VI(7:9,iPoint)  = Buffer_VI(Ux_+4:Uz_+4,iPoint) ! velocity
     BufferLine_VI(10:12,iPoint)= Buffer_VI(Bx_+4:Bz_+4,iPoint) ! magnetic field
     BufferLine_VI(13,iPoint)   = Buffer_VI(p_         ,iPoint) ! pressure
  end do
  
  deallocate(Buffer_VI)
  call line_clean

  !^CFG END RAYTRACE
end subroutine GM_get_for_im_line

!==========================================================================

subroutine GM_get_for_im(Buffer_IIV,iSizeIn,jSizeIn,nVar,NameVar)

  !call stop_mpi('RCM is OFF') !^CFG UNCOMMENT IF NOT RCM
  !^CFG IF RCM BEGIN
  use ModProcMH, ONLY: iProc

  use ModGmImCoupling, ONLY: &
       allocate_gm_im, &
       process_integrated_data, DoTestTec, DoTestIdl, &
       write_integrated_data_tec, write_integrated_data_idl, &
       RCM_lat, RCM_lon, &
       MHD_SUM_vol, MHD_Xeq, MHD_Yeq, MHD_Beq, MHD_SUM_rho, MHD_SUM_p, &
       NoValue

  use ModRaytrace, ONLY: RayResult_VII, RayIntegral_VII, &
       InvB_, Z0x_, Z0y_, Z0b_, RhoInvB_, pInvB_, xEnd_, CLOSEDRAY

  implicit none

  character (len=*), parameter :: NameSub='GM_get_for_im'

  integer, intent(in)                                :: iSizeIn, jSizeIn, nVar
  real, intent(out), dimension(iSizeIn,jSizeIn,nVar) :: Buffer_IIV
  character (len=*), intent(in)                      :: NameVar

  real :: Radius

  logical :: DoTest, DoTestMe
  !--------------------------------------------------------------------------
  if(NameVar /= 'vol:z0x:z0y:bmin:rho:p') &
       call CON_stop(NameSub//' invalid NameVar='//NameVar)

  call CON_set_do_test(NameSub//'_tec', DoTestTec, DoTestMe)
  call CON_set_do_test(NameSub//'_idl', DoTestIdl, DoTestMe)
  call CON_set_do_test(NameSub, DoTest, DoTestMe)

  ! Allocate arrays
  call allocate_gm_im(iSizeIn, jSizeIn)

  ! The RCM ionosphere radius in normalized units
  Radius = (6378.+100.)/6378.  !!! could be derived from Grid_C ?
  call integrate_ray_accurate(iSizeIn, jSizeIn, RCM_lat, RCM_lon, Radius, &
       'InvB,RhoInvB,pInvB,Z0x,Z0y,Z0b')

  if(iProc==0)then
     ! Copy RayResult into small arrays
     MHD_SUM_vol = RayResult_VII(InvB_   ,:,:)
     MHD_Xeq     = RayResult_VII(Z0x_    ,:,:)
     MHD_Yeq     = RayResult_VII(Z0y_    ,:,:)
     MHD_Beq     = RayResult_VII(Z0b_    ,:,:)
     MHD_SUM_rho = RayResult_VII(RhoInvB_,:,:)
     MHD_SUM_p   = RayResult_VII(pInvB_  ,:,:)

     ! Put impossible values if the ray is not closed
     where(RayResult_VII(xEnd_,:,:) <= CLOSEDRAY)
        MHD_Xeq     = NoValue
        MHD_Yeq     = NoValue
        MHD_SUM_vol = 0.0
        MHD_SUM_rho = 0.0
        MHD_SUM_p   = 0.0
        MHD_Beq     = NoValue
     end where
  end if

  deallocate(RayIntegral_VII, RayResult_VII)

  if (iProc == 0) then
     ! Output before processing
     if(DoTestTec)call write_integrated_data_tec
     if(DoTestIdl)call write_integrated_data_idl
     call process_integrated_data
     ! Output after processing
     if(DoTestTec)call write_integrated_data_tec
     if(DoTestIdl)call write_integrated_data_idl

     ! Put results into output buffer
     Buffer_IIV(:,:,InvB_)    = MHD_SUM_vol
     Buffer_IIV(:,:,Z0x_)     = MHD_Xeq
     Buffer_IIV(:,:,Z0y_)     = MHD_Yeq
     Buffer_IIV(:,:,Z0b_)     = MHD_Beq
     Buffer_IIV(:,:,RhoInvB_) = MHD_SUM_rho
     Buffer_IIV(:,:,pInvB_)   = MHD_SUM_p
  end if

end subroutine GM_get_for_im

!==========================================================================
subroutine GM_satinit_for_im(nSats)

  !This subroutine collects the number of satellite files for use in 
  !SWMF GM and IM coupling.   !!!DTW 2007

  !Module variables to use:
  use ModMain,   ONLY: DoImSatTrace
  use ModSatelliteFile, ONLY: nSatellite

  implicit none

  !Subroutine Arguments:
  integer,           intent(out) :: nSats
!--------------------------------------------------------------------------
  
  !If IM sat tracing is on, collect the number of satellites to trace.
  !If IM sat tracing is off, set nSats to zero.
  if (DoImSatTrace) then
     nSats = nSatellite
  else 
     nSats = 0
  endif

end subroutine GM_satinit_for_im

!==========================================================================
subroutine GM_get_sat_for_im(Buffer_III, Buffer_I, nSats)

  ! Subroutine to update and collect satellite locations for IM tracing
  ! !!!DTW 2007
  
  !Modules
  use ModProcMH, ONLY: iProc, iComm
  use ModSatelliteFile, ONLY: Satellite_name, Xsatellite, &
       get_satellite_ray, set_satellite_flags
  use ModMPI
  
  implicit none
  
  !Arguments
  integer, intent(in)               :: nSats
  real, intent(out)                 :: Buffer_III(3,2,nSats)
  character (len=100), intent(out)  :: Buffer_I(nSats)
  
  !Internal variables
  character (len=*), parameter :: NameSub='GM_get_sat_for_im'

  real :: sat_RayVars(5), sat_RayVarsSum(5)
  
  integer :: iSat, iError
  !--------------------------------------------------------------------------
  ! Store satellite names in Buffer_I
  if (iProc == 0) then
     Buffer_I = Satellite_name(1:nSats)
  end if

  do iSat=1, nSats
     ! Update satellite position.
     call set_satellite_flags(iSat)
     call get_satellite_ray(iSat, sat_RayVars)

     ! Reduce values from all 
     call MPI_reduce(sat_RayVars, sat_RayVarsSum, 5, MPI_REAL, MPI_SUM, &
          0, iComm, iError)
     
     ! Store results in Buffer_III
     if (iProc == 0) then 
        Buffer_III(:,1,iSat) = Xsatellite(iSat,:)
        Buffer_III(:,2,iSat) = sat_RayVarsSum(1:3)
     end if
  end do

  !^CFG END RCM      
end subroutine GM_get_sat_for_im

!==========================================================================
