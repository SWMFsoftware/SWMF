
subroutine DGCPM_read_inputs(cFile)

  use ModIoDGCPM
  implicit none

  character (len=*), intent(in) :: cFile

  character (len=iCharLen_) :: line
  logical :: IsThere
  integer :: iError, i

  nInputLines = 1

  inquire(file=cFile,EXIST=IsThere)
  if (.not.IsThere) &
       call stop_dgcpm(cFile//" cannot be found by read_inputs")

  open(iInputUnit_,file=cFile,status="old")

  iError = 0
  do while (iError == 0)

     read(iInputUnit_,'(a)',iostat=iError) line

     if (nInputLines > nInputMaxLines) &
          call stop_dgcpm("Too many lines of input in read_inputs")

     cInputText(nInputLines) = line
     nInputLines = nInputLines + 1

  enddo
           
  close(iInputUnit_)

  if (nInputLines==0) &
       call stop_dgcpm("No lines of input read by read_inputs")
  
end subroutine DGCPM_read_inputs
