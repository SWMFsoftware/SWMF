
subroutine stop_RIM(str)

  implicit none

  character (len=*), intent(in) :: str

  call CON_stop("IE/RIM Error: "//str)

end subroutine stop_RIM
