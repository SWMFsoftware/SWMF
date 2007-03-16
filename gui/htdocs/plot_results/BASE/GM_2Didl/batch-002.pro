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

