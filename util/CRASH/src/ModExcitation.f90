!^CFG COPYRIGHT UM
!
!Data are provided in electron volts (eV)
!=================================



module ModIonizPotential
  implicit none
  PRIVATE !Except
 !/////////////////////////////////////////////////////////////////////////////
!Ground Term
!The data for the first 3 columns are taken from the book 
!"Allen' Astrophysical Quantities, Edn iV" p. 33
!by Allen, Clabon W. Editor: Cox, Arthur N.
!Publisher: Springer (2000) (g_0 for Y=1, Y=2, Y=3

!First 10 elements - full ionizations 

integer, parameter, dimension(10,10) :: GroundTerm_II = reshape(  (/   &
  !   I   !   II    !   III  !   IV   !   V   !  VI   !  VII  ! VIII  !  IX   !   X   !
      2 ,    9*1                                                                     ,&  ! 1 - H
      1 ,      2 ,     8*1                                                           ,&  ! 2 - He
      2 ,      1 ,       2 ,    7*1                                                  ,&  ! 3 - Li
      1 ,      2 ,       1 ,      2 ,   6*1                                          ,&  ! 4 - Be
      6 ,      1 ,       2 ,      1 ,     2 .   5*1                                  ,&  ! 5 - B
      9 ,      6 ,       1 ,      2 ,     1 ,     2 ,   4*1 ,                        ,&  ! 6 - C
      4 ,      9 ,       6 ,      1 ,     2 ,     1 ,     2 ,   3*1                  ,&  ! 7 - N
      9 ,      4 ,       9 ,      6 ,     1 ,     2 ,     1 ,     2 ,   2*1          ,&  ! 8 - O
      6 ,      9 ,       4 ,      9 ,     6 ,     1 ,     2 ,     1 ,     2 ,      1 ,&  ! 9 - F
      1 ,      6 ,       9 ,      4 ,     9 ,     6 ,     1 ,     2 ,     1 ,      2     /),&  !10 - Ne
                                                                                          (/10,10/))



contains
    
end module ModIonizPotential


  

