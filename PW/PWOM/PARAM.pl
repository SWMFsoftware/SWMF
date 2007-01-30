#^CFG FILE _FALSE_
$tree = [{'type' => 'e','attrib' => {'name' => 'PWOM: PW Component'},'name' => 'commandList','content' => [{'type' => 't','content' => '

List of PW commands used in the PARAM.in file

'},{'type' => 'e','attrib' => {'type' => 'string','value' => '$_NameComp/restartOUT','name' => 'NameRestartOutDir'},'name' => 'set','content' => []},{'type' => 'e','attrib' => {'type' => 'string','value' => '$_NameComp/plots','name' => 'NamePlotDir'},'name' => 'set','content' => []},{'type' => 'e','attrib' => {'name' => 'NUMERICAL SCHEME'},'name' => 'commandgroup','content' => [{'type' => 'e','attrib' => {'name' => 'SCHEME'},'name' => 'command','content' => [{'type' => 'e','attrib' => {'type' => 'string','input' => 'select','name' => 'TypeSolver'},'name' => 'parameter','content' => [{'type' => 'e','attrib' => {'default' => 'T','name' => 'Godunov'},'name' => 'option','content' => []},{'type' => 'e','attrib' => {'name' => 'Rusanov'},'name' => 'option','content' => []}]},{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'IsImplicit'},'name' => 'parameter','content' => []},{'type' => 'e','attrib' => {'type' => 'real','min' => '0','default' => '0.05','name' => 'DtVertical'},'name' => 'parameter','content' => []},{'type' => 't','content' => '

#SCHEME
Godunov			TypeSolver
F			IsImplicit
0.05			DtVertical

TypeSolver sets the type of solver which is to be used. Currently only two
options are available: Godunov and Rusanov. IsImplicit determines whether
collision terms are handled implicitly or not. DtVertical sets the time step
for plasma propagating along the field line.

The default values are shown.
'}]},{'type' => 'e','attrib' => {'name' => 'TIMESTEP'},'name' => 'command','content' => [{'type' => 'e','attrib' => {'type' => 'real','min' => '0','default' => '50','name' => 'DtHorizontal'},'name' => 'parameter','content' => []},{'type' => 't','content' => '

#TIMESTEP
50		DtHorizontal

DtHorizontal is the timestep for horizontal motion of the field line.
Default value is shown.
'}]},{'type' => 'e','attrib' => {'name' => 'MOTION'},'name' => 'command','content' => [{'type' => 'e','attrib' => {'type' => 'logical','default' => 'T','name' => 'DoMoveLine'},'name' => 'parameter','content' => []},{'type' => 't','content' => '
#MOTION
T			DoMoveLine

This command determines which to move the field lines as determined by
the horizontal convection, or to hold them in their initial locations. The 
default value is shown.
'}]},{'type' => 'e','attrib' => {'name' => 'ROTATION'},'name' => 'command','content' => [{'type' => 'e','attrib' => {'type' => 'logical','default' => 'T','name' => 'UseCentrifugal'},'name' => 'parameter','content' => []},{'type' => 't','content' => '

#ROTATION
T		UseCentrifugalForce

This command determines whether centrifugal forces should be included. The
default is shown. 
'}]},{'type' => 'e','attrib' => {'name' => 'FAC'},'name' => 'command','content' => [{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'UseJr'},'name' => 'parameter','content' => []},{'type' => 't','content' => '

#FAC
F               UseJr

UseJr determines whether to use field aligned currents to affect the ion
outflow. The default is shown.
'}]}]},{'type' => 'e','attrib' => {'name' => 'INPUT/OUTPUT'},'name' => 'commandgroup','content' => [{'type' => 'e','attrib' => {'name' => 'FIELDLINE'},'name' => 'command','content' => [{'type' => 'e','attrib' => {'type' => 'integer','default' => '1','name' => 'nTotalLine'},'name' => 'parameter','content' => []},{'type' => 't','content' => '

#FIELDLINE
1		nTotalLine

nTotalLine sets the number of field lines included in the simulation.
The default is shown.
'}]},{'type' => 'e','attrib' => {'name' => 'RESTART'},'name' => 'command','content' => [{'type' => 'e','attrib' => {'type' => 'logical','default' => 'T','name' => 'IsRestart'},'name' => 'parameter','content' => []},{'type' => 't','content' => '

#RESTART
T                       IsRestart

If the IsRestart variable is true, then the PWOM uses a restart file.
Otherwise, the PWOM uses a cold start routine. The default is shown.
'}]},{'type' => 'e','attrib' => {'name' => 'SAVEPLOT'},'name' => 'command','content' => [{'type' => 'e','attrib' => {'type' => 'real','default' => '10.0','name' => 'DtSavePlot'},'name' => 'parameter','content' => []},{'type' => 't','content' => '

#SAVEPLOT
10.0			DtSavePlot

The frequency which plot files are written out are defined here. The 
default value is given.
'}]},{'type' => 'e','attrib' => {'name' => 'END'},'name' => 'command','content' => [{'type' => 't','content' => '

#END

The #END command signals the end of the included file or the
end of the PARAM.in file. Lines following the #END command are
ignored. It is not required to use the #END command. The end
of the included file or PARAM.in file is equivalent with an 
#END command in the last line.
'}]}]}]}];