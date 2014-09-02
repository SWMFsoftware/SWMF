module ModRandomNumber

implicit none

public :: ipseudo_random_number    ! Returns a random integer number between 1 and N
public :: rpseudo_random_number    ! Returns a random real number between 0 and 1

contains
  !============================================================================

  real function rpseudo_random_number(iSeedIn, iSeedOut)
    ! rpseudo_random_number returns a pseudo_random_number between 0 and 1

    integer, intent(in)  :: iSeedIn
    integer, intent(out) :: iSeedOut

    integer :: iSeed

    iSeed = iSeedIn*48828125

    IF(iSeed < 0) iSeed=(iSeed+2147483647)+1
    if(iSeed==0) iSeed=1

    iSeedOut = iSeed

    rpseudo_random_number=FLOAT(iSeed)/2147483647

  end function rpseudo_random_number

  !=============================================================================
  integer function ipseudo_random_number( n, ix, iy, iz )
    !
    !*******************************************************************************
    !
    !! ipseudo_random_number returns a random integer between 1 and N.
    !! The maximum N is 2147483647 for a default-sized integer
    !
    !  Discussion:
    !
    !   This function returns a uniformly distributed pseudo-
    !   random integer in the range 1 to N.
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
    !
    !  Parameters:
    !
    !    Input, integer N, the maximum value to be returned.
    !
    !    Input/output, integer IX, IY, IZ = Integer seeds initialized to 
    !    values in the range 1 to 30,000 before the first call to ipseudo_random_number, and 
    !    not altered between subsequent calls (unless a sequence of random 
    !    numbers is to be repeated by reinitializing the seeds).
    !
    !    Output, integer ipseudo_random_number, a random integer in the range 1 to N.
    !
    !  Local parameters:
    !
    !    U = Pseudo-random number uniformly distributed in the interval (0,1).
    !    X = Pseudo-random number in the range 0 to 3 whose fractional part is U.
    !
    implicit none
    !
    integer:: ix, iy, iz
    integer:: n
    real:: u
    real:: x
    !
    ix = mod ( 171 * ix, 30269 )
    iy = mod ( 172 * iy, 30307 )
    iz = mod ( 170 * iz, 30323 )

    x = ( real ( ix ) / 30269.0E+00 ) &
         + ( real ( iy ) / 30307.0E+00 ) &
         + ( real ( iz ) / 30323.0E+00 )

    u = x - int ( x )
    ipseudo_random_number = real ( n ) * u + 1.0E+00

    return
  end function ipseudo_random_number

end module ModRandomNumber
