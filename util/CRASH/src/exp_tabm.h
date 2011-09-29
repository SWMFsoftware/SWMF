 ! *** exp_tabm.h ***
 ! *
 ! * IN : y , OUT: ey=exp(-Y), ex_y=1-ey , |Y|=min(|y|,ex_max) 
 ! *
	if(y.ge.ex_zero) then
	 ex_y=min(ex_max,y)
	 ex_u=ex_y*ex_sdu
	 ex_i=int(ex_u)
	 ex_u=ex_y -ex_i*ex_du
	 ey=ex_tab(ex_i)*(ex_one-ex_c*ex_u)
	else
	 ex_y=min(ex_max,-y)
	 ex_u=ex_y*ex_sdu
	 ex_i=int(ex_u)
	 ex_u=ex_y -ex_i*ex_du
	 ey=(ex_one+ex_cm*ex_u)/ex_tab(ex_i)
	end if
	if(ex_i.ne.0) ex_y=ex_one-ey
