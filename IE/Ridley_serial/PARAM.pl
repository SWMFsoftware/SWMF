#^CFG FILE _FALSE_
$tree = [{'content' => [{'content' => [{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'default' => 'T','type' => 'logical','name' => 'UseStrict'}},{'content' => '

#STRICT
T                       UseStrict

If true then stop when parameters are incompatible. If false, try to
correct parameters and continue. Default is true, i.e. strict mode.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'STRICT'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'default' => '-1','min' => '-1','type' => 'integer','name' => 'iDebugLevel'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'default' => '0','min' => '0','max' => '$_nProc','type' => 'integer','name' => 'iDebugProc'}},{'content' => '

#DEBUG
2			iDebugLevel
3			iDebugProc

The iDebugLevel variable sets the level of debug information
for the processor selected by iDebugProc.

Default is iDebugLevel=-1 which is no debug info on any and iDebugProc=0.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'DEBUG'}}],'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'Testing'}},{'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! OUTPUT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
','type' => 't'},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'default' => 'ionosphere','length' => '100','type' => 'string','name' => 'NameIonoDir'}},{'content' => '

#IONODIR
IE/Plots		NameIonoDir

The NameIonoDir variable contains the name of the directory to store
output files. Default value is "IE/ionosphere".

','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'IONODIR'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'default' => '0','min' => '0','max' => '10','type' => 'integer','name' => 'nPlotFile'}},{'content' => [{'content' => [{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'default' => 'T','name' => 'IDL','value' => 'idl'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'TecPlot','value' => 'tec'}}],'name' => 'part','type' => 'e','attrib' => {'required' => 'T','type' => 'string','name' => 'TypePlotForm','input' => 'select'}},{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'default' => 'T','name' => 'Minimum','value' => 'min'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Maximum','value' => 'max'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Aurora','value' => 'aur'}}],'name' => 'part','type' => 'e','attrib' => {'required' => 'T','type' => 'string','name' => 'TypePlotVar','input' => 'select'}}],'name' => 'parameter','type' => 'e','attrib' => {'min' => '2','max' => '2','type' => 'strings','name' => 'StringPlot'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'min' => '-1','type' => 'integer','name' => 'DnOutput'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'min' => '-1','type' => 'real','name' => 'DtOutput'}}],'name' => 'for','type' => 'e','attrib' => {'to' => '$nPlotFile','from' => '1'}},{'content' => '

#SAVEPLOT
2			nFile
min idl			StringPlot
-1			DnOuput
10.0			DtOutput [sec]
max tec			StringPlot
-1			DnOuput
20.0			DtOutput [sec]

The StringPlot variable consists of two string parts:
the TypePlotVar string can be \'min\', \'max\' or \'aur\'
corresponding to a minimum, maximum or auroral set of
plot variables. The other part TypePlotForm can be \'idl\'
or \'tec\' corresponding to plot files for IDL or TecPlot.
The DnOuput and DtOutput variables determine the frequency
of saves in terms of time step or physical time.

The default is that no plotfiles are saved.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'alias' => 'IE_SAVEPLOT','name' => 'SAVEPLOT'}}],'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'Output'}},{'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!! Physical parameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

','type' => 't'},{'content' => [{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => '(0) constant pedersen and 0 Hall','value' => '0'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => '(1) constant pedersen and Hall','value' => '1'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => '(2) F107 flux and and constant Hall','value' => '2'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => '(3) Simple oval','value' => '3'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'default' => 'T','name' => '(5) Realistic oval','value' => '5'}}],'name' => 'parameter','type' => 'e','attrib' => {'type' => 'integer','name' => 'ConductanceModel','input' => 'select'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'default' => 'F','type' => 'logical','name' => 'UseFullCurrent'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'default' => 'F','type' => 'logical','name' => 'UseFakeRegion2'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'default' => '150','min' => '0','type' => 'real','name' => 'F107Flux'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'default' => '0.25','min' => '0.001','type' => 'real','name' => 'StarLightPedConductance'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'default' => '0.25','min' => '0','type' => 'real','name' => 'PolarCapPedConductance'}},{'content' => [{'content' => '
		F107Flux must be positive for ConductanceModel = 3 or 5
	','type' => 't'}],'name' => 'rule','type' => 'e','attrib' => {'expr' => 'not ($ConductanceModel=~/3|5/ and $F107Flux==0)'}},{'content' => '
#IONOSPHERE
5			iConductanceModel
F			UseFullCurrent
F			UseFakeRegion2
150.0			F107Flux
0.25			StarLightPedConductance
0.25			PolarCapPedConductance

The iConductanceModel variable determines which ionosphere model is used:\\\\
  0 - uses a constant Pedersen conductance which is set by 
      StarLightPedConductance\\\\
  1 - uses a constant Pedersen conductance which is set by 
      StarLightPedConductance, and a constant Hall conductance
      which is set by PolarCapPedConductance\\\\
  2 - uses a solar EUV combined with a nightside conductance, so
      it uses F107Flux and StarLightPedConductance\\\\
  3 - uses solar EUV, nightside, and crude oval, so uses
      F107Flux, StarLightPedConductance, and PolarCapPedConductance,
      since a polar cap is defined with the oval.\\\\
  4 - doesn\'t work\\\\
  5 - more realistic oval, uses same variables as 3.\\\\

The UseFullCurrent and UseFakeRegion2 logicals were used in the past
to test various algorithmic choices. They should have a false value now.
The default values are shown by the example above.

','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'IONOSPHERE'}},{'content' => [{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'default' => 'T','name' => 'north'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'south'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'cpcpmin'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'average'}}],'name' => 'parameter','type' => 'e','attrib' => {'type' => 'string','name' => 'TypeImCouple','input' => 'select'}},{'content' => '

#IM
average				TypeImCouple

The TypeImCouple parameter determines which hemisphere the IM component 
is coupled to. If the value is \'north\' or \'south\', the potential and radial
current are sent from the corresponding magnetic hemisphere. For \'cpcpmin\'
the hemisphere with the lower cross polar cap potential is selected.
For TypeImCouple=\'average\' the potential and radial current are averaged
for the north and south hemispheres.

The default value is \'north\', which is backward compatible,
and it requires no communication between the IE processors.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'IM'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'default' => 'F','type' => 'logical','name' => 'DoCoupleUaCurrent'}},{'content' => '

#UA
T		DoCoupleUaCurrent

The DoCoupleUaCurrent parameter determines if the field aligned current
calculated by the UA component should be used. Usually the currents are
dominated by the field aligned currents provided by the GM component.
The coupling with the UA currents is still experimental.
The default value is false, i.e. the UA currents are not included.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'UA'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'length' => '100','type' => 'string','name' => 'AMIEFileNorth'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'length' => '100','type' => 'string','name' => 'AMIEFileSouth'}},{'content' => '

#AMIEFILES
IE/amie.north
IE/amie.south

Set the files to read the AMIE data from.

Default is not reading AMIE files.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'AMIEFILES'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'default' => 'F','type' => 'logical','name' => 'UseSPS'}},{'content' => '

#SPS
T			UseSPS

The UseSPS parameter indicates if the serial potential solver is used.
This is the default in the SWMF and it cannot be modified in the SWMF.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'SPS'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'length' => '100','type' => 'string','name' => 'NameOfModelDir'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'length' => '100','type' => 'string','name' => 'NameOfEFieldModel'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'length' => '100','type' => 'string','name' => 'NameOfAuroralModel'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'length' => '100','type' => 'string','name' => 'NameOfSolarModel'}},{'content' => '

#BACKGROUND
dir			NameOfModelDir
weimer96		NameOfEFieldModel
ihp			NameOfAuroralModel
xyz			NameOfSolarModel

This command cannot be used in the SWMF.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_IsStandAlone','name' => 'BACKGROUND'}}],'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'Physical parameters'}},{'content' => [{'content' => '
	Output directory IE/ionosphere should exist!
','type' => 't'}],'name' => 'rule','type' => 'e','attrib' => {'expr' => '-d \'IE/ionosphere\' or not $_IsFirstSession'}}],'name' => 'commandList','type' => 'e','attrib' => {'name' => 'Ridley Ionosphere Model: IE Component'}}];