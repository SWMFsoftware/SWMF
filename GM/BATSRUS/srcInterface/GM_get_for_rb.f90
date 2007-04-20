!^CFG COPYRIGHT UM
!^CMP FILE RB

!==========================================================================
subroutine GM_get_for_rb_trace(iSizeIn, jSizeIn, NameVar, nVarLine, nPointLine)

  ! Do ray tracing for RB. 
  ! Provide total number of points along rays 
  ! and the number of variables to pass to RB
  use ModGmRbCoupling, ONLY: allocate_gm_rb, RB_lat, RB_lon

  use CON_line_extract, ONLY: line_get
  implicit none
  integer, intent(in)           :: iSizeIn, jSizeIn
  character (len=*), intent(in) :: NameVar
  integer, intent(out)          :: nVarLine, nPointLine
  real :: Radius

  character (len=*), parameter :: NameSub='GM_get_for_rb_trace'
  !---------------------------------------------------------------------

  if(NameVar /= 'Z0x:Z0y:Z0b:B_I:S_I:IMF') &
       call CON_stop(NameSub//' invalid NameVar='//NameVar)

  ! Allocate arrays
  call allocate_gm_rb(iSizeIn, jSizeIn)

  ! The RB ionosphere radius in normalized units
  Radius = (6378.+100.)/6378.  !!! could be derived from Grid_C ?
  call integrate_ray_accurate(iSizeIn, jSizeIn, RB_lat, RB_lon, Radius, &
       NameVar)

  call line_get(nVarLine, nPointLine)
  nVarLine = 2 ! We only pass B and Length to RB

end subroutine GM_get_for_rb_trace

!==========================================================================

subroutine GM_get_for_rb(Buffer_IIV, iSizeIn, jSizeIn, nVar, &
     BufferLine_VI, nVarLine, nPointLine)

  !call stop_mpi('RAYTRACE is OFF') !^CFG UNCOMMENT IF NOT RAYTRACE
  !^CFG IF RAYTRACE BEGIN

  use ModProcMH, ONLY: iProc

  use ModGmRbCoupling, ONLY: &
       process_integrated_data, &
       write_integrated_data_tec, write_integrated_data_idl, &
       MHD_SUM_vol, MHD_Xeq, MHD_Yeq, MHD_Beq, MHD_SUM_rho, MHD_SUM_p, NoValue

  use ModRaytrace, ONLY: RayResult_VII, RayIntegral_VII, &
       InvB_, Z0x_, Z0y_, Z0b_, RhoInvB_, pInvB_, xEnd_, CLOSEDRAY

  use ModVarIndexes, ONLY: Bx_, Bz_

  use ModPhysics, ONLY: No2Si_V, UnitN_, UnitU_, UnitB_, UnitP_, &
       SW_n, SW_Ux, SW_Uy, SW_Uz, SW_p, SW_Bx, SW_By, SW_Bz

  use CON_line_extract, ONLY: line_get, line_clean

  implicit none

  character (len=*), parameter :: NameSub='GM_get_for_rb'

  integer, intent(in)                                :: iSizeIn, jSizeIn, nVar
  real, intent(out), dimension(iSizeIn,jSizeIn,nVar) :: Buffer_IIV

  integer, intent(in) :: nPointLine, nVarLine
  real, intent(out)   :: BufferLine_VI(nVarLine, nPointLine)

  integer :: nVarExtract, nPoint, iPoint
  real, allocatable :: Buffer_VI(:,:)

  logical :: DoTestTec, DoTestIdl
  logical :: DoTest, DoTestMe
  !--------------------------------------------------------------------------

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
     BufferLine_VI(1,iPoint) = Buffer_VI(0,iPoint)    ! Length
     BufferLine_VI(2,iPoint) = &
          sqrt(sum(Buffer_VI(3+Bx_:3+Bz_,iPoint)**2)) ! |B|
  end do
     
  deallocate(RayIntegral_VII, RayResult_VII, Buffer_VI)
  call line_clean

  ! Output before processing
  if(DoTest .or. DoTestTec)call write_integrated_data_tec
  if(DoTest .or. DoTestIdl)call write_integrated_data_idl
  call process_integrated_data
  ! Output after processing
  if(DoTest .or. DoTestTec)call write_integrated_data_tec
  if(DoTest .or. DoTestIdl)call write_integrated_data_idl

  ! Put results into output buffer
  Buffer_IIV(:,:,InvB_)    = MHD_SUM_vol
  Buffer_IIV(:,:,Z0x_)     = MHD_Xeq
  Buffer_IIV(:,:,Z0y_)     = MHD_Yeq
  Buffer_IIV(:,:,Z0b_)     = MHD_Beq
  Buffer_IIV(:,:,RhoInvB_) = MHD_SUM_rho
  Buffer_IIV(:,:,pInvB_)   = MHD_SUM_p
  Buffer_IIV(:,:,pInvB_+1) = 0.

  ! Send solar wind values in the array of the extra integral
  ! This is a temporary solution. RB should use MHD_SUM_rho and MHD_SUM_p
  Buffer_IIV(1,:,pInvB_+1) = SW_n  * No2Si_V(UnitN_)
  Buffer_IIV(2,:,pInvB_+1) = SW_Ux * No2Si_V(UnitU_)
  Buffer_IIV(3,:,pInvB_+1) = SW_Uy * No2Si_V(UnitU_)
  Buffer_IIV(4,:,pInvB_+1) = SW_Uz * No2Si_V(UnitU_)
  Buffer_IIV(5,:,pInvB_+1) = SW_Bx * No2Si_V(UnitB_)
  Buffer_IIV(6,:,pInvB_+1) = SW_By * No2Si_V(UnitB_)
  Buffer_IIV(7,:,pInvB_+1) = SW_Bz * No2Si_V(UnitB_)
  Buffer_IIV(8,:,pInvB_+1) = SW_p  * No2Si_V(UnitP_)

  !^CFG END RAYTRACE
end subroutine GM_get_for_rb
