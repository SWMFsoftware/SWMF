!!BOP -------------------------------------------------------------------
!
! MODULE: ModIO_f77_to_swmf - specific handling of F77 units.
!
! DESCRIPTION:
!
! This routine assigns new file name *iFile* in SWMF and is compatible with
! a program written in f77.
!
! In order to comply with SWMF requirements, you must do the following:
!
! Step 1: Create a file named *stdout.h* that contains the following three
! lines of code:
!
!       common/stdout/iStdOut,prefix
!       integer iStdOut
!       character*4 prefix
!
! Here the number of characters in the last line is assumed to be 4, so
! adjust this number to meet your needs. In the case of solar particles
! module, prefix='SP: ', i.e. it contains 4 characters. It is simple as
! that.
! NOTE: The convention is that the first two letters correspond to
! the two-symbol name of the component in SWMF (like IH,GM, etc.). 
!
! Step 2: If anywhere in your f77 program there is a write statement
! of the kind
!
!           write(*,*)'B_I(j) =',B_I(j)
!           write(6,*)'B_I(j) =',B_I(j)
!
! replace that line by
!
!           write(iStdout,*)prefix,'B_I(j) =',B_I(j)
!
! However, remember to add the following line in the header of your
! program or subroutine, right after the line *implicit none*
!
!           include 'stdout.h'
!
! Also, in the "INIT" section of your program, add a call to
!
!           call set_stdout('SP: ')
!
! in order to assign the value of variable *prefix*.
!
! Step 3: If anywhere in your f77 program there is a handling of file
! unit, before you open the unit make the call:
!
!           call get_io_unit_new(iFile)
!
! Remembder to declare the integer variable iFile before that:
!
!           integer::iFile
!
! Then, instead of opening or closing the file unit (assumed to be 11
! here) by:
!
!           open(11,file='SP_dist.dat',status='unknown',form='formatted')
!           close(11)
!
! USE the following way:
!
!           open(iFile,file='SP_dist.dat',status='unknown',form='formatted')
!           close(iFile)
!
! Then, if you are reading/writing from/to the file, instead of:
!
!           read(11,*) B_I(j)
!           write(11,*) 'B_I(j) =',B_I(j)
!
! USE the following:
!
!           read(iFile,*) B_I(j)
!           write(iFile,*) 'B_I(j) =',B_I(j)
!
! That is it you have to do. Good luck with SWMF!
! 
!INTERFACE:
      subroutine get_io_unit_new(iFileNew)
!EOP
      implicit none
  !--------------------------------------------------------------------------!
      integer iFileNew
      integer :: iUnit
      logical :: IsExisting, IsOpened
      integer :: iError

      do iUnit = 10, 10000
         inquire (
     1   unit   = iUnit,       
     2   exist  = IsExisting,  
     3   opened = IsOpened,    
     4   iostat = iError)
         if (IsExisting .and. .not. IsOpened
     1           .and. iError == 0) then
            iFileNew = iUnit
            return
         end if
      end do

      iFileNew = -1
      end subroutine get_io_unit_new
!BOP
!INTERFACE:
      subroutine set_stdout(NameComp)
!DESCRIPTION:
! This routine sets the prefix and assigns a standard/non-standard number    
! for STD output, compatible with f77. The f77 model is assumed not to be     
! parallel so that the prefix consists of two-symbol component name followed 
! by a semicolon and space. 
!EOP                                                  
      implicit none
      include 'stdout.h'
      character*4  NameComp
      iStdOut = 6
      prefix = NameComp
      DoWriteAll=.false.
      end subroutine set_stdout
!============================================================================! 
