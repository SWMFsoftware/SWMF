Module ModTsyInput
  ! max no. of SW-IMF, dst, Kp data pts 
  integer, parameter :: nswmax=40000, ndstmax=45000, nKpmax=250  
  
  real    :: tf,dsth(ndstmax),tdst(ndstmax),bxw(nswmax),byw(nswmax), &
             bzw(nswmax),timf(nswmax),xnswa(nswmax),vswa(nswmax),tsw(nswmax),&
             ndst,nimf,nsw,js
  integer :: iyear,iday,imod=3
  logical :: UseSmooth=.false.


end Module ModTsyInput
