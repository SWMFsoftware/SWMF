#^CFG FILE _FALSE_
$tree = [{'attrib' => {'name' => 'SEP with pitch-angle dependence: SP Component'},'content' => [{'attrib' => {'name' => 'Testing'},'content' => [{'attrib' => {'name' => 'VERBOSE'},'content' => [{'attrib' => {'min' => '0','default' => '0','type' => 'integer','name' => 'iVerbose'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

#VERBOSE
0                       iVerbose

For test purpose set iVerbose>0. SP code will provide a huge
output of the information to STDOUT
','type' => 't'}],'type' => 'e','name' => 'command'}],'type' => 'e','name' => 'commandgroup'},{'attrib' => {'name' => 'Output'},'content' => [{'content' => '

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! OUTPUT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
','type' => 't'},{'attrib' => {'name' => 'PLOT'},'content' => [{'attrib' => {'min' => '1','default' => '5','type' => 'integer','name' => 'nPlot'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

#PLOT
5                       nPlot

The iDtPlot parameter defines the frequency for saving the information .
to file. Default is 5

','type' => 't'}],'type' => 'e','name' => 'command'}],'type' => 'e','name' => 'commandgroup'},{'attrib' => {'name' => 'MAGNETIC FIELD LINES'},'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!    MAGNETIC FIELD LINES       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
','type' => 't'},{'attrib' => {'name' => 'LINE'},'content' => [{'attrib' => {'default' => '0','type' => 'real','name' => 'XyzLine_Dx'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '0','type' => 'real','name' => 'XyzLine_Dy'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '0','type' => 'real','name' => 'XyzLine_Dz'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '1','default' => '1.1','type' => 'real','name' => 'RBoundSC'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '1','default' => '21','type' => 'real','name' => 'RBoundIH'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

 #LINE
 0.0 		        XyzLine_D(x_)
 0.0                    XyzLine_D(y_)
 0.0                    XyzLine_D(z_)
 1.1 	                RBoundSC             
 21.0                   RBoundIH            

XyzLine_D are three coordinates of the starting point of the magnetic field line.
Coordinates are expressed in Rsun, measured in SGI frame of reference.
Default values are Zero, at these values and at any other satisfying the condition
sum(XyzLine_D**2).lt.1 the XyzLine_D is set to coinside with the Earth position at the time the
SP component is switched on.
RBoundSC is the radius of sphere, inside of which the points of the line are considered as
not belonging to the Solar Corona. The integration of the original magnetic field line is from
XyzLine_D towards the Sun. On reaching the surface of the radius of RBoundSC the integration
inside the SC component stops. Expressed in Rsun, default value is 1.1 [Rsun]
RBoundIH is the radius of sphere, inside of which the points of the line are considered as
not belonging to the Inner Heliosphere. The integration of the original magnetic field line is from
XyzLine_D towards the Sun. On reaching the surface of the radius of RBoundIH the integration
inside the IH component stops.  Expressed in Rsun, default value is 21.0 [Rsun]

','type' => 't'}],'type' => 'e','name' => 'command'}],'type' => 'e','name' => 'commandgroup'},{'attrib' => {'name' => 'CONTROL'},'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!               CONTROL         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

','type' => 't'},{'attrib' => {'name' => 'DORUN'},'content' => [{'attrib' => {'default' => 'T','type' => 'logical','name' => 'DoRun'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '
#DORUN
T                       DoRun

When DoRun is set to non-default value .false., the component skips all computations

','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'SAVEMHDATA'},'content' => [{'attrib' => {'default' => 'F','type' => 'logical','name' => 'SaveMhData'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

#SAVEMHDATA
F                       SaveMhData

Default .false.. When .true. saves all the MHD data coming from the coupler(s). Setting
DoRun=.false. (or .true.) and SaveMhData=.true. one can the perform (or repeat) complete 
run for SP only, without running the other components.

','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'DOREADMHDATA'},'content' => [{'attrib' => {'default' => 'F','type' => 'logical','name' => 'DoReadMhData'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

#DOREADMHDATA
F                       DoReadMhData

Default .false.. When .true. the SP component is expected to work as alone component or 
at least it is assumed not to be coupled with MHD model. The MHD data are read from
the flat files in this case, rather then from the couplers

','type' => 't'}],'type' => 'e','name' => 'command'}],'type' => 'e','name' => 'commandgroup'},{'attrib' => {'name' => 'SHOCK WAVE SENSOR'},'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!   SHOCK WAVE SENSOR       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
','type' => 't'},{'attrib' => {'name' => 'RTRANSIENT'},'content' => [{'attrib' => {'default' => '1','type' => 'real','name' => 'rTransient'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

#RTRANSIENT
1.0                    rTransient

Default=1.0, Expressed in Rsun. The code searches for the shock wave at the
distances from the Solar centre larger than rTransient

','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'NSMOOTH'},'content' => [{'attrib' => {'default' => '0','type' => 'integer','name' => 'nSmooth'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

#NSMOOTH
0                      nSmooth

Default=0. Usually the MHD data from the coupler are messed with some spurious occilation.
If nSmooth>0, the SP code nSmooth times repeats some smoothing precedure for the
data data obtained from the coupler.
','type' => 't'}],'type' => 'e','name' => 'command'}],'type' => 'e','name' => 'commandgroup'}],'type' => 'e','name' => 'commandList'}];