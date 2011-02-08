!^CFG FILE DISSFLUX
!^CFG COPYRIGHT UM

! This interface subroutine is needed to avoid circular dependency between
! ModResistivity.f90 and ModUser.f90

subroutine set_resistivity(iBlock)

  use ModResistivity
  use ModUser, ONLY: user_set_resistivity
  implicit none

  integer, intent(in) :: iBlock
  character (len=*), parameter :: NameSub = 'set_resistivity'
  !--------------------------------------------------------------------------
  select case(TypeResistivity)
  case('constant')
     Eta_GB(:,:,:,iBlock) = Eta0
  case('spitzer')
     call spitzer_resistivity(iBlock, Eta_GB(:,:,:,iBlock))
  case('anomalous')
     call anomalous_resistivity(iBlock, Eta_GB(:,:,:,iBlock))
  case('user')
     call user_set_resistivity(iBlock, Eta_GB(:,:,:,iBlock))
  case default
     call stop_mpi(NameSub//' : invalid TypeResistivity='//TypeResistivity)
  end select

  call mask_resistivity(iBlock, Eta_GB(:,:,:,iBlock))

end subroutine set_resistivity

