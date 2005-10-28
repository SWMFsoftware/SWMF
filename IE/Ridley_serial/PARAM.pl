#^CFG FILE _FALSE_
$tree = [{'content' => [{'content' => [{'content' => [{'content' => [],'attrib' => {'default' => 'T','type' => 'logical','name' => 'UseStrict'},'type' => 'e','name' => 'parameter'},{'content' => '

#STRICT
T                       UseStrict

If true then stop when parameters are incompatible. If false, try to
correct parameters and continue. Default is true, i.e. strict mode.
','type' => 't'}],'attrib' => {'name' => 'STRICT'},'type' => 'e','name' => 'command'},{'content' => [{'content' => [],'attrib' => {'default' => '-1','type' => 'integer','name' => 'iDebugLevel','min' => '-1'},'type' => 'e','name' => 'parameter'},{'content' => [],'attrib' => {'default' => '0','type' => 'integer','name' => 'iDebugProc','max' => '$_nProc','min' => '0'},'type' => 'e','name' => 'parameter'},{'content' => '

#DEBUG
2			iDebugLevel
3			iDebugProc

The iDebugLevel variable sets the level of debug information
for the processor selected by iDebugProc.

Default is iDebugLevel=-1 which is no debug info on any and iDebugProc=0.
','type' => 't'}],'attrib' => {'name' => 'DEBUG'},'type' => 'e','name' => 'command'}],'attrib' => {'name' => 'Testing'},'type' => 'e','name' => 'commandgroup'},{'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! OUTPUT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
','type' => 't'},{'content' => [{'content' => [],'attrib' => {'default' => 'ionosphere','length' => '100','type' => 'string','name' => 'NameIonoDir'},'type' => 'e','name' => 'parameter'},{'content' => '

#IONODIR
IE/Plots		NameIonoDir

The NameIonoDir variable contains the name of the directory to store
output files. Default value is "IE/ionosphere".

','type' => 't'}],'attrib' => {'name' => 'IONODIR'},'type' => 'e','name' => 'command'},{'content' => [{'content' => [],'attrib' => {'default' => '0','type' => 'integer','name' => 'nPlotFile','max' => '10','min' => '0'},'type' => 'e','name' => 'parameter'},{'content' => [{'content' => [{'content' => [{'content' => [],'attrib' => {'default' => 'T','name' => 'IDL','value' => 'idl'},'type' => 'e','name' => 'option'},{'content' => [],'attrib' => {'name' => 'TecPlot','value' => 'tec'},'type' => 'e','name' => 'option'}],'attrib' => {'required' => 'T','type' => 'string','input' => 'select','name' => 'TypePlotForm'},'type' => 'e','name' => 'part'},{'content' => [{'content' => [],'attrib' => {'default' => 'T','name' => 'Minimum','value' => 'min'},'type' => 'e','name' => 'option'},{'content' => [],'attrib' => {'name' => 'Maximum','value' => 'max'},'type' => 'e','name' => 'option'},{'content' => [],'attrib' => {'name' => 'UA','value' => 'uam'},'type' => 'e','name' => 'option'},{'content' => [],'attrib' => {'name' => 'Aurora','value' => 'aur'},'type' => 'e','name' => 'option'}],'attrib' => {'required' => 'T','type' => 'string','input' => 'select','name' => 'TypePlotVar'},'type' => 'e','name' => 'part'}],'attrib' => {'type' => 'strings','name' => 'StringPlot','max' => '2','min' => '2'},'type' => 'e','name' => 'parameter'},{'content' => [],'attrib' => {'type' => 'integer','name' => 'DnOutput','min' => '-1'},'type' => 'e','name' => 'parameter'},{'content' => [],'attrib' => {'type' => 'real','name' => 'DtOutput','min' => '-1'},'type' => 'e','name' => 'parameter'}],'attrib' => {'from' => '1','to' => '$nPlotFile'},'type' => 'e','name' => 'for'},{'content' => '

#SAVEPLOT
2			nFile
min idl			StringPlot
-1			DnOuput
10.0			DtOutput [sec]
max tec			StringPlot
-1			DnOuput
20.0			DtOutput [sec]

The StringPlot variable consists of two string parts:
the TypePlotVar string can be \'min\', \'max\', \'uam\' or \'aur\'
corresponding to a minimum, maximum, upper atmosphere or auroral set of
plot variables. The other part TypePlotForm can be \'idl\'
or \'tec\' corresponding to plot files for IDL or TecPlot.
The DnOuput and DtOutput variables determine the frequency
of saves in terms of time step or physical time.

The default is that no plotfiles are saved.
','type' => 't'}],'attrib' => {'alias' => 'IE_SAVEPLOT','name' => 'SAVEPLOT'},'type' => 'e','name' => 'command'}],'attrib' => {'name' => 'Output'},'type' => 'e','name' => 'commandgroup'},{'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!! Physical parameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

','type' => 't'},{'content' => [{'content' => [{'content' => [],'attrib' => {'name' => '(0) constant pedersen and 0 Hall','value' => '0'},'type' => 'e','name' => 'option'},{'content' => [],'attrib' => {'name' => '(1) constant pedersen and Hall','value' => '1'},'type' => 'e','name' => 'option'},{'content' => [],'attrib' => {'name' => '(2) F107 flux and and constant Hall','value' => '2'},'type' => 'e','name' => 'option'},{'content' => [],'attrib' => {'name' => '(3) Simple oval','value' => '3'},'type' => 'e','name' => 'option'},{'content' => [],'attrib' => {'default' => 'T','name' => '(5) Realistic oval','value' => '5'},'type' => 'e','name' => 'option'}],'attrib' => {'type' => 'integer','input' => 'select','name' => 'ConductanceModel'},'type' => 'e','name' => 'parameter'},{'content' => [],'attrib' => {'default' => 'F','type' => 'logical','name' => 'UseFullCurrent'},'type' => 'e','name' => 'parameter'},{'content' => [],'attrib' => {'default' => 'F','type' => 'logical','name' => 'UseFakeRegion2'},'type' => 'e','name' => 'parameter'},{'content' => [],'attrib' => {'default' => '150','type' => 'real','name' => 'F107Flux','min' => '0'},'type' => 'e','name' => 'parameter'},{'content' => [],'attrib' => {'default' => '0.25','type' => 'real','name' => 'StarLightPedConductance','min' => '0.001'},'type' => 'e','name' => 'parameter'},{'content' => [],'attrib' => {'default' => '0.25','type' => 'real','name' => 'PolarCapPedConductance','min' => '0'},'type' => 'e','name' => 'parameter'},{'content' => [{'content' => '
		F107Flux must be positive for ConductanceModel = 3 or 5
	','type' => 't'}],'attrib' => {'expr' => 'not ($ConductanceModel=~/3|5/ and $F107Flux==0)'},'type' => 'e','name' => 'rule'},{'content' => '
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

','type' => 't'}],'attrib' => {'name' => 'IONOSPHERE'},'type' => 'e','name' => 'command'},{'content' => [{'content' => [{'content' => [],'attrib' => {'default' => 'T','name' => 'north'},'type' => 'e','name' => 'option'},{'content' => [],'attrib' => {'name' => 'south'},'type' => 'e','name' => 'option'},{'content' => [],'attrib' => {'name' => 'cpcpmin'},'type' => 'e','name' => 'option'},{'content' => [],'attrib' => {'name' => 'average'},'type' => 'e','name' => 'option'}],'attrib' => {'type' => 'string','input' => 'select','name' => 'TypeImCouple'},'type' => 'e','name' => 'parameter'},{'content' => '

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
','type' => 't'}],'attrib' => {'name' => 'IM'},'type' => 'e','name' => 'command'},{'content' => [{'content' => [],'attrib' => {'default' => 'F','type' => 'logical','name' => 'DoCoupleUaCurrent'},'type' => 'e','name' => 'parameter'},{'content' => '

#UA
T		DoCoupleUaCurrent

The DoCoupleUaCurrent parameter determines if the field aligned current
calculated by the UA component should be used. Usually the currents are
dominated by the field aligned currents provided by the GM component.
The coupling with the UA currents is still experimental.
The default value is false, i.e. the UA currents are not included.
','type' => 't'}],'attrib' => {'name' => 'UA'},'type' => 'e','name' => 'command'},{'content' => [{'content' => [],'attrib' => {'length' => '100','type' => 'string','name' => 'AMIEFileNorth'},'type' => 'e','name' => 'parameter'},{'content' => [],'attrib' => {'length' => '100','type' => 'string','name' => 'AMIEFileSouth'},'type' => 'e','name' => 'parameter'},{'content' => '

#AMIEFILES
IE/amie.north
IE/amie.south

Set the files to read the AMIE data from.

Default is not reading AMIE files.
','type' => 't'}],'attrib' => {'name' => 'AMIEFILES'},'type' => 'e','name' => 'command'},{'content' => [{'content' => [],'attrib' => {'default' => 'F','type' => 'logical','name' => 'UseSPS'},'type' => 'e','name' => 'parameter'},{'content' => '

#SPS
T			UseSPS

The UseSPS parameter indicates if the serial potential solver is used.
This is the default in the SWMF and it cannot be modified in the SWMF.
','type' => 't'}],'attrib' => {'name' => 'SPS'},'type' => 'e','name' => 'command'},{'content' => [{'content' => [],'attrib' => {'length' => '100','type' => 'string','name' => 'NameOfModelDir'},'type' => 'e','name' => 'parameter'},{'content' => [],'attrib' => {'length' => '100','type' => 'string','name' => 'NameOfEFieldModel'},'type' => 'e','name' => 'parameter'},{'content' => [],'attrib' => {'length' => '100','type' => 'string','name' => 'NameOfAuroralModel'},'type' => 'e','name' => 'parameter'},{'content' => [],'attrib' => {'length' => '100','type' => 'string','name' => 'NameOfSolarModel'},'type' => 'e','name' => 'parameter'},{'content' => '

#BACKGROUND
dir			NameOfModelDir
weimer96		NameOfEFieldModel
ihp			NameOfAuroralModel
xyz			NameOfSolarModel

This command cannot be used in the SWMF.
','type' => 't'}],'attrib' => {'if' => '$_IsStandAlone','name' => 'BACKGROUND'},'type' => 'e','name' => 'command'}],'attrib' => {'name' => 'Physical parameters'},'type' => 'e','name' => 'commandgroup'},{'content' => [{'content' => '
	Output directory IE/ionosphere should exist!
','type' => 't'}],'attrib' => {'expr' => '-d \'IE/ionosphere\' or not $_IsFirstSession'},'type' => 'e','name' => 'rule'}],'attrib' => {'name' => 'Ridley Ionosphere Model: IE Component'},'type' => 'e','name' => 'commandList'}];