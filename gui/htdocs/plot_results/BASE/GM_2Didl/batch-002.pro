;  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
;  For more information, see http://csem.engin.umich.edu/tools/swmf
filename='file.out'
!x.range = [-20, 10]
!y.range = [-15, 15]
func='rho'
plottitle='rho'
plotmode='contbargrid'
transform='n'
set_device,'image.ps',/port
.r animate
close_device

