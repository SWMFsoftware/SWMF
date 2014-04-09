module ModBlk
  use ModPar,   ONLY: inmax, jnmax
  implicit none

  ! blktri parameters
  real, dimension(1:inmax-5) :: ar, br ,cr
  real, dimension(1:jnmax-5) :: ath, bth, cth, bthk

end module ModBlk
