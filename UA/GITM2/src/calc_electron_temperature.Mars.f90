
subroutine calc_electron_temperature(iBlock)

  use ModGITM

  implicit none

  integer, intent(in) :: iBlock

  call report("Electron Density", 2)

  eTemperature(:,:,:,iBlock) = Temperature(:,:,:,iBlock) * TempUnit(:,:,:)

end subroutine calc_electron_temperature







