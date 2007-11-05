
! ---------------------------------------------------

subroutine IO_set_IMF_Bz_single(BzIn)

  use ModIndices
  implicit none
  real, intent(in) :: BzIn

  nIndices_V(imf_bz_) = 1
  Indices_TV(1,imf_bz_) = BzIn

end subroutine IO_set_IMF_Bz_single

! ---------------------------------------------------

subroutine IO_set_IMF_By_single(ByIn)

  use ModIndices
  implicit none
  real, intent(in) :: ByIn

  nIndices_V(imf_by_) = 1
  Indices_TV(1,imf_by_) = ByIn

end subroutine IO_set_IMF_By_single

! ---------------------------------------------------

subroutine IO_set_IMF_Bx_single(BxIn)

  use ModIndices
  implicit none
  real, intent(in) :: BxIn

  nIndices_V(imf_bx_) = 1
  Indices_TV(1,imf_bx_) = BxIn

end subroutine IO_set_IMF_Bx_single

! ---------------------------------------------------

subroutine IO_set_SW_v_single(VIn)

  use ModIndices
  implicit none
  real, intent(in) :: VIn

  nIndices_V(sw_v_) = 1
  Indices_TV(1,sw_v_) = VIn

end subroutine IO_set_SW_v_single

! ---------------------------------------------------

subroutine IO_set_hpi_single(HpiIn)

  use ModIndices
  implicit none
  real, intent(in) :: HpiIn

  nIndices_V(Hpi_) = 1
  Indices_TV(1,Hpi_) = HpiIn

end subroutine IO_set_hpi_single

! ---------------------------------------------------

subroutine IO_set_f107_single(F107In)

  use ModIndices
  implicit none
  real, intent(in) :: F107In

  nIndices_V(f107_) = 1
  Indices_TV(1,f107_) = F107In

end subroutine IO_set_f107_single

! ---------------------------------------------------

subroutine IO_set_f107a_single(F107aIn)

  use ModIndices
  implicit none
  real, intent(in) :: F107aIn

  nIndices_V(f107a_) = 1
  Indices_TV(1,f107a_) = F107aIn

end subroutine IO_set_f107a_single

! ---------------------------------------------------

subroutine IO_set_kp_single(kpIn)

  use ModIndices
  implicit none
  real, intent(in) :: kpIn

  nIndices_V(kp_) = 1
  Indices_TV(1,kp_) = kpIn

end subroutine IO_set_kp_single


! ---------------------------------------------------

subroutine IO_set_ap_single(apIn)

  use ModIndices
  implicit none
  real, intent(in) :: apIn

  nIndices_V(ap_) = 1
  Indices_TV(1,ap_) = apIn

end subroutine IO_set_ap_single


