#^CFG FILE _FALSE_
$tree = [{'content' => [{'content' => '

List of PW commands used in the PARAM.in file

','type' => 't'},{'content' => [],'name' => 'set','type' => 'e','attrib' => {'value' => '$_NameComp/restartOUT','type' => 'string','name' => 'NameRestartOutDir'}},{'content' => [],'name' => 'set','type' => 'e','attrib' => {'value' => '$_NameComp/plots','type' => 'string','name' => 'NamePlotDir'}},{'content' => [{'content' => [{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'default' => 'T','name' => 'Godunov'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Rusanov'}}],'name' => 'parameter','type' => 'e','attrib' => {'type' => 'string','name' => 'TypeSolver','input' => 'select'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'default' => 'F','type' => 'logical','name' => 'IsImplicit'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'min' => '0','default' => '0.05','type' => 'real','name' => 'DtVertical'}},{'content' => '

#SCHEME
Godunov			TypeSolver
F			IsImplicit
0.05			DtVertical

TypeSolver sets the type of solver which is to be used. Currently only two
options are available: Godunov and Rusanov. IsImplicit determines whether
collision terms are handled implicitly or not. DtVertical sets the time step
for plasma propagating along the field line.

The default values are shown.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'SCHEME'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'min' => '0','default' => '50','type' => 'real','name' => 'DtHorizontal'}},{'content' => '

#TIMESTEP
50		DtHorizontal

DtHorizontal is the timestep for horizontal motion of the field line.
Default value is shown.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'TIMESTEP'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'default' => 'T','type' => 'logical','name' => 'DoMoveLine'}},{'content' => '
#MOTION
T			DoMoveLine

This command determines which to move the field lines as determined by
the horizontal convection, or to hold them in their initial locations. The 
default value is shown.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'MOTION'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'default' => 'T','type' => 'logical','name' => 'UseCentrifugal'}},{'content' => '

#ROTATION
T		UseCentrifugalForce

This command determines whether centrifugal forces should be included. The
default is shown. 
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'ROTATION'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'default' => 'F','type' => 'logical','name' => 'UseJr'}},{'content' => '

#FAC
F               UseJr

UseJr determines whether to use field aligned currents to affect the ion
outflow. The default is shown.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'FAC'}}],'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'NUMERICAL SCHEME'}},{'content' => [{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'length' => '100','type' => 'string','name' => 'StringTest'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'min' => '0','default' => '0','type' => 'integer','name' => 'iProcTest'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'min' => '1','default' => '1','type' => 'integer','name' => 'iLineTest'}},{'content' => '
#TEST
PW_move_line		StringTest
0			iProcTest
2			iLineTest

Set test parameters. The subroutines to be tested are listed in StringTest,
the tested processor is iProcTest, and the tested field line is iLineTest.
If iLineTest is 0, all field lines produce test output.

Default is an empty StringTest, i.e. no test output is produced.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'TEST'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'default' => '1','type' => 'integer','name' => 'nTotalLine'}},{'content' => '

#FIELDLINE
1		nTotalLine

nTotalLine sets the number of field lines included in the simulation.
The default is shown.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'FIELDLINE'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'default' => 'T','type' => 'logical','name' => 'IsRestart'}},{'content' => '

#RESTART
T                       IsRestart

If the IsRestart variable is true, then the PWOM uses a restart file.
Otherwise, the PWOM uses a cold start routine. The default is shown.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'RESTART'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'default' => '10.0','type' => 'real','name' => 'DtSavePlot'}},{'content' => '

#SAVEPLOT
10.0			DtSavePlot

The frequency which plot files are written out are defined here. The 
default value is given.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'SAVEPLOT'}},{'content' => [{'content' => '

#END

The #END command signals the end of the included file or the
end of the PARAM.in file. Lines following the #END command are
ignored. It is not required to use the #END command. The end
of the included file or PARAM.in file is equivalent with an 
#END command in the last line.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'END'}}],'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'INPUT/OUTPUT'}}],'name' => 'commandList','type' => 'e','attrib' => {'name' => 'PWOM: PW Component'}}];