!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------

subroutine report(str, iLevel)

  use ModInputs, only : iDebugLevel
  implicit none

  character (len=*), intent(in) :: str
  integer, intent(in) :: iLevel
  character (len = 11) :: cArrow
  integer :: i

  if (iDebugLevel < iLevel .or. iDebugLevel > 10) return

  do i=1,iLevel
     cArrow(i:i) = "="
  enddo
  cArrow(iLevel+1:iLevel+1) = ">"

  write(*,*) cArrow(1:iLevel+1), " ",str

end subroutine report

!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------

subroutine stop_gitm(str)

  use ModGITM
  use ModInputs, only: IsFramework
  use ModMpi
  implicit none

  character (len=*), intent(in) :: str
  integer :: ierror, erno

  if (IsFramework) then
     call CON_stop("UA/GITM Error: "//str)
  else
     write(*,*)'Stopping execution! iProc=',iProc,' with msg=',str
     call MPI_abort(iCommGITM, erno, ierror)
     stop
  endif

end subroutine stop_gitm

!!---------------------------------------------------------------------------
!!
!!---------------------------------------------------------------------------
!
!subroutine merge_str(str1, str2)
!
!  character (len=100) :: str1, str2, temp
!  integer :: i, j, k
!
!  i = 1
!  do while (iachar(str1(i:i)) /= 32 .and. &
!            iachar(str1(i:i)) /= 9  .and. &
!            i < 100) 
!     i=i+1
!  enddo
!
!  j = 1
!  do while (iachar(str2(j:j)) /= 32 .and. &
!            iachar(str2(j:j)) /= 9  .and. &
!            j < 100) 
!     j=j+1
!  enddo
!
!  temp = str1
!  do k = i,100
!     temp(k:k) = ' '
!  enddo
!
!  if (i+j-1 > 100) j = 100 - i + 1
!
!  temp(i:i+j-1) = str2(1:j)
!
!  str2 = temp
!
!end subroutine merge_str

!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------

subroutine i2s(iValue, cOut, iLength)

  integer, intent(in)            :: iValue, iLength
  character (len=*), intent(out) :: cOut
  character (len=10)             :: cFormat
  integer                        :: i

  if (iLength < 10) then
     write(cFormat,"('(I',I1,')')") iLength
  else
     write(cFormat,"('(I',I2,')')") iLength
  endif

  write(cOut, cFormat) iValue

  do i=1, iLength
    if (cOut(i:i) == ' ') cOut(i:i) = '0'
  enddo

end subroutine i2s
