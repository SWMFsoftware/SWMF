;  Copyright (C) 2002 Regents of the University of Michigan, 
; portions used with permission 
;  For more information, see http://csem.engin.umich.edu/tools/swmf
; This file plots the 10-th picture from run/IO2/y=0*.outs into the
; PostScript file example.eps
; Usage    : idl EXAMPLE.pro
; or in IDL: @EXAMPLE

filename='run/IO2/y=0*.outs'
npict=10
read_data
func='rho bx;by'
plotmode='contfill streamover'
plottitle='rho and B'
set_device,'example.eps'
loadct,3
plot_data
close_device,/pdf
