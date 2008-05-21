!^CFG COPYRIGHT UM
!
!Data is given in electron volts (eV)
!Data taken from the book "Allen' Astrophysical Quantities, fourth edition" pg. 36
!by Allen, Clabon W. Editor: Cox, Arthur N.
!Published in 2000 by Springer
!=================================

module ModIonizPotential
  implicit none
  PRIVATE !Except
 ! integer,parameter,dimension(3) :: nCell2_D=(/nI/2,nJ/2,nK/2/)
 !/////////////////////////////////////////////////////////////////////////////
 
  real, parameter, dimension(10,10) :: cPotentials10_II = reshape(  (/   &
!     1    |  2     :     3     :    4    : 5      :   6   :   7    :  8      : 9      :  10       :
  13.59844, 0.       , 0.      , 0.      , 0.     , 0.     , 0.     , 0.     , 0.     ,  0.       ,&  ! 1 - H
  24.58741, 54.41778, 0.      , 0.      , 0.     , 0.     , 0.     , 0.     , 0.     ,  0.       ,&  ! 2 - He
  5.39172 , 75.64018, 122.454, 0.      , 0.     , 0.     , 0.     , 0.     , 0.     ,  0.       ,&  ! 3 - Li
  9.32263 , 18.21116, 153.897, 217.713, 0.     , 0.     , 0.     , 0.     , 0.     ,  0.       ,&  ! 4 - Be
  8.29803 , 25.15484, 37.931 , 259.366, 340.22, 0.     , 0.     , 0.     , 0.     ,  0.       ,&  ! 5 - B
  11.26030, 24.38332, 47.888 , 64.492 , 392.08, 489.98, 0.     , 0.     , 0.     ,  0.       ,&  ! 6 - C
  14.53414, 29.6013 , 47.449 , 77.472 , 97.89 , 552.06, 667.03, 0.     , 0.     ,  0.       ,&  ! 7 - N
  13.61806, 35.11730, 54.936 , 77.413 , 113.90, 138.12, 739.29, 871.41, 0.     ,  0.       ,&  ! 8 - O
  17.42282, 34.97082, 62.708 , 87.140 , 114.24, 157.17, 185.19, 953.91, 1103.1,  0.       ,&  ! 9 - F
  21.56454, 40.96328, 63.45  , 97.12  , 126.21, 157.93, 207.28, 239.10, 1195.8, 1362.2 /),&  !10 - Ne
                                                                       (/10,10/))
!public :: get_ionization_potential10

  
contains
  subroutine get_ioniz_potential10( nZ, cPotential_I)
         integer,intent(in) :: nZ

         real,intent(out), dimension(nZ) :: cPotential_I

         cPotential_I = cPotentials10_II( 1:nZ, nZ)

  end subroutine get_ioniz_potential10  
end module ModIonizPotential
  


