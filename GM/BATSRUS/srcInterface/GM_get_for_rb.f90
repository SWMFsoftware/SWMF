!^CFG COPYRIGHT UM


!==========================================================================
!==========================================================================
!==========================================================================
subroutine GM_get_for_rb(Buffer_IIV,iSizeIn,jSizeIn,nVar,NameVar)

  use ModProcMH, ONLY: iProc

  use ModGmRbCoupling, ONLY: &
       allocate_gm_rb, &
       process_integrated_data, &
       write_integrated_data_tec, write_integrated_data_idl, &
       RB_lat, RB_lon, &
       MHD_SUM_vol, MHD_Xeq, MHD_Yeq, MHD_Beq, MHD_SUM_rho, MHD_SUM_p

  use ModRaytrace, ONLY: UseAccurateIntegral, RayResult_VII, RayIntegral_VII, &
       InvB_, Z0x_, Z0y_, Z0b_, RhoInvB_, pInvB_, xEnd_, CLOSEDRAY

  implicit none

  character (len=*), parameter :: NameSub='GM_get_for_rb'

  integer, intent(in)                                :: iSizeIn, jSizeIn, nVar
  real, intent(out), dimension(iSizeIn,jSizeIn,nVar) :: Buffer_IIV
  character (len=*), intent(in)                      :: NameVar

  real :: Radius

  logical :: DoTest, DoTestMe
  !--------------------------------------------------------------------------
  if(NameVar /= 'vol:z0x:z0y:bmin:rho:p') &
       call CON_stop(NameSub//' invalid NameVar='//NameVar)

  call CON_set_do_test(NameSub, DoTest, DoTestMe)

  ! Allocate arrays
  call allocate_gm_rb(iSizeIn, jSizeIn)

  ! Use accurate integration only

  ! The RB ionosphere radius in normalized units
  Radius = (6378.+100.)/6378.  !!! could be derived from Grid_C ?
  call integrate_ray_accurate(iSizeIn, jSizeIn, RB_lat, RB_lon, Radius)

  if(iProc==0)then
     ! Copy RayResult into small arrays used in old algorithm
     MHD_SUM_vol = RayResult_VII(InvB_   ,:,:)
     MHD_Xeq     = RayResult_VII(Z0x_    ,:,:)
     MHD_Yeq     = RayResult_VII(Z0y_    ,:,:)
     MHD_Beq     = RayResult_VII(Z0b_    ,:,:)
     MHD_SUM_rho = RayResult_VII(RhoInvB_,:,:)
     MHD_SUM_p   = RayResult_VII(pInvB_  ,:,:)

     ! Put impossible values if ray was not found for a lat-lon grid cell
     where(RayResult_VII(xEnd_,:,:) <= CLOSEDRAY)
        MHD_Xeq     = -99999.0
        MHD_Yeq     = -99999.0
        MHD_SUM_vol = 0.0
        MHD_SUM_rho = 0.0
        MHD_SUM_p   = 0.0
        MHD_Beq     = -99999.0
     end where

  end if

  deallocate(RayIntegral_VII, RayResult_VII)

  if (iProc == 0) then
     if(DoTest)call write_integrated_data_tec  ! TecPlot output before processing
     if(DoTest)call write_integrated_data_idl  ! IDL output before processing
     call process_integrated_data
     if(DoTest)call write_integrated_data_tec  ! TecPlot output
     if(DoTest)call write_integrated_data_idl  ! IDL     output

     ! Put results into output buffer
     Buffer_IIV(:,:,InvB_)    = MHD_SUM_vol
     Buffer_IIV(:,:,Z0x_)     = MHD_Xeq
     Buffer_IIV(:,:,Z0y_)     = MHD_Yeq
     Buffer_IIV(:,:,Z0b_)     = MHD_Beq
     Buffer_IIV(:,:,RhoInvB_) = MHD_SUM_rho
     Buffer_IIV(:,:,pInvB_)   = MHD_SUM_p
  end if

end subroutine GM_get_for_rb
