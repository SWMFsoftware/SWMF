#^CFG FILE _FALSE_
$tree = [{'type' => 'e','attrib' => {'name' => 'SEP with pitch-angle dependence: SP Component'},'content' => [{'type' => 'e','attrib' => {'name' => 'Testing'},'content' => [{'type' => 'e','attrib' => {'name' => 'VERBOSE'},'content' => [{'type' => 'e','attrib' => {'min' => '0','type' => 'integer','default' => '0','name' => 'iVerbose'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '

#VERBOSE
0                       iVerbose

For testing purposes set iVerbose>0 and the SP code will write a lot of
information to STDOUT.
'}],'name' => 'command'}],'name' => 'commandgroup'},{'type' => 'e','attrib' => {'name' => 'Output'},'content' => [{'type' => 't','content' => '

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! OUTPUT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'type' => 'e','attrib' => {'name' => 'PLOT'},'content' => [{'type' => 'e','attrib' => {'min' => '1','type' => 'integer','default' => '5','name' => 'DnPlot'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '

#PLOT
5                       DnPlot

The DnPlot parameter defines the frequency for saving the information
to file. Default value is 5.
'}],'name' => 'command'}],'name' => 'commandgroup'},{'type' => 'e','attrib' => {'name' => 'MAGNETIC FIELD LINES'},'content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!    MAGNETIC FIELD LINES       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'type' => 'e','attrib' => {'name' => 'LINE'},'content' => [{'type' => 'e','attrib' => {'type' => 'real','default' => '0','name' => 'xLine'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '0','name' => 'yLine'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '0','name' => 'zLine'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '1','type' => 'real','default' => '1.1','name' => 'rBoundSC'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '1','type' => 'real','default' => '21','name' => 'rBoundIH'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '

#LINE
0.0 		        xLine    [Rs]
0.0                     yLine    [Rs]
0.0                     zLine    [Rs]
1.1 	                rBoundSc [Rs]
21.0                    rBoundIh [Rs]

The xLine, yLine and zLine are the coordinates of the starting point of the 
magnetic field line given in units of solar radii (Rs) in the 
Heliographic Inertial (HGI) frame of reference.
If the starting point is given as all zeros, which is the 
center of the Sun (or any position inside the Sun),
then the starting point is set to coincide with the position of the Earth 
at the time the SP component is switched on.
The default values are zero for these 3 variables.

The rBoundSc variable defines the radius of a sphere, 
inside which the points of the magnetic field
line are considered as not belonging to the Solar Corona. 
The integration of the original magnetic field line is from
the starting point towards the Sun. On reaching the surface of the radius of 
rBoundSc, the integration inside the SC component stops. 
The default value is rBoundSc = 1.1 [Rs].

The rBoundIh variable defines the radius of a sphere, inside which 
the points of the magnetic field line are considered as not belonging 
to the Inner Heliosphere. 
The integration of the original magnetic field line is from
the starting point towards the Sun. On reaching the surface of the radius of 
rBoundIh, the integration inside the IH component stops.  
The default value is rBoundIh = 21.0 [Rs].
'}],'name' => 'command'}],'name' => 'commandgroup'},{'type' => 'e','attrib' => {'name' => 'CONTROL'},'content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!               CONTROL         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

'},{'type' => 'e','attrib' => {'name' => 'DORUN'},'content' => [{'type' => 'e','attrib' => {'type' => 'logical','default' => 'T','name' => 'DoRun'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '
#DORUN
T                       DoRun

When DoRun is set to the value .false., the 
component skips all computations.
The default value is true, i.e. computations are performed.
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'SAVEMHDATA'},'content' => [{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'SaveMhData'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '

#SAVEMHDATA
F                       DoSaveMhData

When DoSaveMhData is true, the SEP code saves all the MHD data coming from 
the coupler(s). 
Setting DoRun=.false. (or .true.) and SaveMhData=.true. one can 
perform (or repeat) complete run for SP only, without running the other 
components.

Default value is false. 
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'DOREADMHDATA'},'content' => [{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'DoReadMhData'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '

#DOREADMHDATA
F                       DoReadMhData

When .true. the SP component is expected to work as a single component or 
at least it is assumed not to be coupled with an MHD model. 
The MHD data are read from the flat files in this case, 
rather then received thrhough the couplers.
Default value is false. 
'}],'name' => 'command'}],'name' => 'commandgroup'},{'type' => 'e','attrib' => {'name' => 'SHOCK WAVE SENSOR'},'content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!   SHOCK WAVE SENSOR       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

'},{'type' => 'e','attrib' => {'name' => 'RTRANSIENT'},'content' => [{'type' => 'e','attrib' => {'type' => 'real','default' => '1','name' => 'rTransient'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '

#RTRANSIENT
1.0                    rTransient [Rs]

The code searches for the shock wave at a
distance larger than rTransient from the Solar center.
The rTransient is given in units of solar radii. 
Default value is 1.
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'NSMOOTH'},'content' => [{'type' => 'e','attrib' => {'type' => 'integer','default' => '0','name' => 'nSmooth'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '

#NSMOOTH
0                      nSmooth

Usually the MHD data from the coupler contains some spurious occilation.
If nSmooth > 0, the SP code repeats some smoothing procedure nSmooth times 
for the data data obtained from the coupler.
Default value is 0. 
'}],'name' => 'command'}],'name' => 'commandgroup'},{'type' => 'e','attrib' => {'expr' => '-d \'SP\' or not $_IsFirstSession'},'content' => [{'type' => 't','content' => '
	Directory SP should exist!
'}],'name' => 'rule'}],'name' => 'commandList'}];