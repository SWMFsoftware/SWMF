!^CFG COPYRIGHT UM
!^CMP FILE IM

!==========================================================================
!==========================================================================
!==========================================================================
subroutine GM_get_for_im(Buffer_IIV,iSizeIn,jSizeIn,nVar,NameVar)

  !call stop_mpi('RCM is OFF') !^CFG UNCOMMENT IF NOT RCM
  !^CFG IF RCM BEGIN
  use ModProcMH, ONLY: iProc

  use ModGmImCoupling, ONLY: &
       allocate_gm_im, MHD_compute_from_raytrace, &
       process_integrated_data, DoTestTec, DoTestIdl, &
       write_integrated_data_tec, write_integrated_data_idl, &
       RCM_lat, RCM_lon, &
       MHD_SUM_vol, MHD_Xeq, MHD_Yeq, MHD_Beq, MHD_SUM_rho, MHD_SUM_p, &
       NoValue

  use ModRaytrace, ONLY: UseAccurateIntegral, RayResult_VII, RayIntegral_VII, &
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

  if(UseAccurateIntegral)then
     ! The RCM ionosphere radius in normalized units
     Radius = (6378.+100.)/6378.  !!! could be derived from Grid_C ?
     call integrate_ray_accurate(iSizeIn, jSizeIn, RCM_lat, RCM_lon, Radius)

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

  else
     ! Make sure that we have ray tracing
     call ray_trace

     !compute various values
     call MHD_compute_from_raytrace
  end if

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

  !^CFG END RCM
end subroutine GM_get_for_im
