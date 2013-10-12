;  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
;  For more information, see http://csem.engin.umich.edu/tools/swmf
filename='file.out'
func='p'
plottitle='p'
plotmode='contbar'
transform='n'
set_device,'image.ps',/port
.r animate
close_device

