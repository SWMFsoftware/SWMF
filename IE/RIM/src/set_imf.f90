!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

subroutine set_imf

  use ModRIM
  use ModParamRIM
  use ModIndicesInterfaces

  implicit none

  real :: temp
  integer :: iError

  call get_IMF_Bz(CurrentTime, temp, iError)
  call IO_SetIMFBz(temp)

  if (iError /= 0) then
     write(*,*) &
          "Code Error in get_IMF_Bz called from set_imf, part of RIM"
     write(*,*) "Code : ",iError
     call stop_RIM("Stopping in set_imf, part of RIM")
  endif

  if (iDebugLevel > 1) write(*,*) "==> Bz : ",temp

  call get_IMF_By(CurrentTime, temp, iError)
  call IO_SetIMFBy(temp)

  if (iError /= 0) then
     write(*,*) &
          "Code Error in get_IMF_By called from set_imf, part of RIM"
     call stop_RIM("Stopping in set_imf, part of RIM")
  endif

  if (iDebugLevel > 1) write(*,*) "==> IMF By : ",temp

  call get_SW_V(CurrentTime, temp, iError)
  call IO_SetSWV(temp)

  if (iError /= 0) then
     write(*,*) "Code Error in get_sw_v called from set_imf, part of RIM"
     call stop_RIM("Stopping in set_imf, part of RIM")
  endif

  if (iDebugLevel > 1) write(*,*) "==> Solar Wind Velocity : ",temp

  call get_SW_N(CurrentTime, temp, iError)
  if (iError /= 0) then
     temp = 5.0
     iError = 0
  endif
  call IO_SetSWN(temp)

  if (iError /= 0) then
    write(*,*) "Code Error in get_sw_n called from set_imf, part of RIM"
    call stop_RIM("Stopping in set_imf, part of RIM")
  endif

  call get_HPI(CurrentTime, temp, iError)
  call IO_SetHPI(temp)

  if (iError /= 0) then
     write(*,*) "Code Error in get_HPI called from set_imf, part of RIM"
     call stop_RIM("Stopping in set_imf, part of RIM")
  endif

  if (iDebugLevel > 1) write(*,*) "==> Hemispheric Power Index : ",temp

end subroutine set_imf
