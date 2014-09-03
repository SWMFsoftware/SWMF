module ModRandomNumber

  implicit none

  public :: random_integer ! Returns a random integer number between 1 and N
  public :: random_real    ! Returns a random real number between 0 and 1

contains
  !============================================================================

  real function random_real(iSeed)

    ! returns a pseudo_random_number between 0 and 1
    ! iSeed needs to be initialized and saved by the calling routine

    integer, intent(inout)  :: iSeed
    !---------------------------------------------------------
    iSeed = iSeed*48828125

    if(iSeed < 0) iSeed = (iSeed + 2147483647) + 1
    if(iSeed == 0) iSeed = 1

    random_real = iSeed/2147483647.0

  end function random_real

  !============================================================================
  integer function random_integer(n, iSeed1, iSeed2, iSeed3)

    integer, intent(in)   :: n                 ! maximum value to be returned
    integer, intent(inout):: iSeed1, iSeed2, iSeed3  ! integer seeds

    ! Returns a random integer between 1 and N.
    ! The maximum N is 2147483647 for a default-sized 4-byte integer.
    ! The three integer seeds should be initialized to 
    ! values in the range 1 to 30,000 before the first call and 
    ! not altered between subsequent calls (unless a sequence of random 
    ! numbers is to be repeated by reinitializing the seeds).
    !
    !  Author:
    !
    !    Robert Renka,
    !    Department of Computer Science,
    !    University of North Texas,
    !    renka@cs.unt.edu
    !
    !  Reference:  
    !
    !    B. A. Wichmann and I. D. Hill, 
    !    An Efficient and Portable Pseudo-random Number Generator,
    !    Applied Statistics, 
    !    Volume 31, Number 2, 1982, pages 188-190.

    real:: x ! Pseudo-random number in the range 0 to 3
    real:: u ! Pseudo-random number uniformly distributed in the interval (0,1)
    !--------------------------------------------------------------------------
    iSeed1 = mod(171*iSeed1, 30269)
    iSeed2 = mod(172*iSeed2, 30307)
    iSeed3 = mod(170*iSeed3, 30323)

    x = iSeed1/30269.0 + iSeed2/30307.0 + iSeed3/30323.0

    u = x - int(x)
    random_integer = n*u + 1.0

  end function random_integer

end module ModRandomNumber
