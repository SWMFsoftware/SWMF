 !-------
 ! \
module M_expTab
  ! /
  !  tabulated exponential, 10 times faster that hard wired exp()
  !-------
  implicit none
  integer,save :: nb_Exp=-1	! flag of preparation
  integer,parameter :: ex_nb=15000
  real,parameter :: ex_du=1e-2,ex_mm=1e-25 &
       ,ex_min=1e-4,ex_max=ex_nb*ex_du &
       ,ex_zero=0.,ex_one=1.
  real,save :: ex_sdu,ex_c,ex_cm,ex_tab(0:ex_nb)
  real :: ex_u,ex_y
  integer :: ex_i
  !-------
contains
  !-------
  subroutine exp_tab8()
    implicit none
    if(nb_exp.eq.ex_nb) return
    nb_exp=ex_nb
    do ex_i=0,ex_nb
       ex_u=ex_i*ex_du
       ex_tab(ex_i)=exp(-ex_u)
    end do
    ex_sdu=ex_one/ex_du
    ex_c=(ex_one-ex_tab(1))*ex_sdu
    ex_cm=(ex_one/ex_tab(1)-ex_one)*ex_sdu
    ex_y=0.
    write(*,*)'.. tabulation of 1/exp(0..',ex_max,') is done.'
  end subroutine exp_tab8
 !-------
 ! \
end module M_expTab
 ! /
 !-------
 
