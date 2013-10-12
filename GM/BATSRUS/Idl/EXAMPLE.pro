;  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
;  For more information, see http://csem.engin.umich.edu/tools/swmf
;^CFG COPYRIGHT VAC_UM
; This file plots the 10-th picture from run/IO2/example.out into the
; PostScript file example.ps
; Usage    : idl EXAMPLE.pro
; or in IDL: @EXAMPLE
filename='run/IO2/example.out'
npict=10
.r getpict
headline='Example PostScript Plot'
func='rho bx;by'
plotmode='contfill stream'
multiplot=1
plottitle='rho and B'
set_plot,'PS'
device,filename='example.ps',/color,bits=8,xsize=18,ysize=18,$
	xoffset=1,yoffset=2
loadct,3
.r plotfunc
device,/close
set_plot,'X'
