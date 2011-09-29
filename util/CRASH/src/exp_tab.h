!*** exp_tab.h ***
!*
!* IN : y , OUT: ey=exp(-Y) , Y=min(y,ex_max) ; 0 <= y 
!*
! ey=exp(-y)	! pour verification
! ex_y=ex_one-ey
!*		
!c	ex_y=min(ex_max,max(ex_min,y))
	ex_y=min(ex_max,max(ex_zero,u))
	ex_u=ex_y*ex_sdu
	ex_i=int(ex_u)
!c	ex_u=frac(ex_u)
	ex_u=ex_y -ex_i*ex_du
	ey=ex_tab(ex_i)*(ex_one-ex_c*ex_u)
