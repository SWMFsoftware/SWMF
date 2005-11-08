#^CFG FILE _FALSE_
$tree = [{'name' => 'commandList','content' => [{'name' => 'commandgroup','content' => [{'name' => 'command','content' => [{'name' => 'parameter','content' => [],'attrib' => {'name' => 'UseStrict','default' => 'T','type' => 'logical'},'type' => 'e'},{'content' => '

#STRICT
T                       UseStrict

If true then stop when parameters are incompatible. If false, try to
correct parameters and continue. Default is true, i.e. strict mode.
','type' => 't'}],'attrib' => {'name' => 'STRICT'},'type' => 'e'},{'name' => 'command','content' => [{'name' => 'parameter','content' => [],'attrib' => {'name' => 'iDebugLevel','default' => '-1','min' => '-1','type' => 'integer'},'type' => 'e'},{'name' => 'parameter','content' => [],'attrib' => {'name' => 'iDebugProc','default' => '0','min' => '0','type' => 'integer','max' => '$_nProc'},'type' => 'e'},{'content' => '

#DEBUG
2			iDebugLevel
3			iDebugProc

The iDebugLevel variable sets the level of debug information
for the processor selected by iDebugProc.

Default is iDebugLevel=-1 which is no debug info on any and iDebugProc=0.
','type' => 't'}],'attrib' => {'name' => 'DEBUG'},'type' => 'e'}],'attrib' => {'name' => 'Testing'},'type' => 'e'},{'name' => 'commandgroup','content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! OUTPUT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
','type' => 't'},{'name' => 'command','content' => [{'name' => 'parameter','content' => [],'attrib' => {'name' => 'NameIonoDir','default' => 'ionosphere','length' => '100','type' => 'string'},'type' => 'e'},{'content' => '

#IONODIR
IE/Plots		NameIonoDir

The NameIonoDir variable contains the name of the directory to store
output files. Default value is "IE/ionosphere".

','type' => 't'}],'attrib' => {'name' => 'IONODIR'},'type' => 'e'},{'name' => 'command','content' => [{'name' => 'parameter','content' => [],'attrib' => {'name' => 'nPlotFile','default' => '0','min' => '0','type' => 'integer','max' => '10'},'type' => 'e'},{'name' => 'for','content' => [{'name' => 'parameter','content' => [{'name' => 'part','content' => [{'name' => 'option','content' => [],'attrib' => {'name' => 'IDL','value' => 'idl','default' => 'T'},'type' => 'e'},{'name' => 'option','content' => [],'attrib' => {'name' => 'TecPlot','value' => 'tec'},'type' => 'e'}],'attrib' => {'name' => 'TypePlotForm','required' => 'T','input' => 'select','type' => 'string'},'type' => 'e'},{'name' => 'part','content' => [{'name' => 'option','content' => [],'attrib' => {'name' => 'Minimum','value' => 'min','default' => 'T'},'type' => 'e'},{'name' => 'option','content' => [],'attrib' => {'name' => 'Maximum','value' => 'max'},'type' => 'e'},{'name' => 'option','content' => [],'attrib' => {'name' => 'UA','value' => 'uam'},'type' => 'e'},{'name' => 'option','content' => [],'attrib' => {'name' => 'Aurora','value' => 'aur'},'type' => 'e'}],'attrib' => {'name' => 'TypePlotVar','required' => 'T','input' => 'select','type' => 'string'},'type' => 'e'}],'attrib' => {'name' => 'StringPlot','min' => '2','type' => 'strings','max' => '2'},'type' => 'e'},{'name' => 'parameter','content' => [],'attrib' => {'name' => 'DnOutput','min' => '-1','type' => 'integer'},'type' => 'e'},{'name' => 'parameter','content' => [],'attrib' => {'name' => 'DtOutput','min' => '-1','type' => 'real'},'type' => 'e'}],'attrib' => {'to' => '$nPlotFile','from' => '1'},'type' => 'e'},{'content' => '

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
','type' => 't'}],'attrib' => {'name' => 'SAVEPLOT','alias' => 'IE_SAVEPLOT'},'type' => 'e'}],'attrib' => {'name' => 'Output'},'type' => 'e'},{'name' => 'commandgroup','content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!! Physical parameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

','type' => 't'},{'name' => 'command','content' => [{'name' => 'parameter','content' => [{'name' => 'option','content' => [],'attrib' => {'name' => '(0) constant pedersen and 0 Hall','value' => '0'},'type' => 'e'},{'name' => 'option','content' => [],'attrib' => {'name' => '(1) constant pedersen and Hall','value' => '1'},'type' => 'e'},{'name' => 'option','content' => [],'attrib' => {'name' => '(2) F107 flux and and constant Hall','value' => '2'},'type' => 'e'},{'name' => 'option','content' => [],'attrib' => {'name' => '(3) Simple oval','value' => '3'},'type' => 'e'},{'name' => 'option','content' => [],'attrib' => {'name' => '(5) Realistic oval','value' => '5','default' => 'T'},'type' => 'e'}],'attrib' => {'name' => 'ConductanceModel','input' => 'select','type' => 'integer'},'type' => 'e'},{'name' => 'parameter','content' => [],'attrib' => {'name' => 'UseFullCurrent','default' => 'F','type' => 'logical'},'type' => 'e'},{'name' => 'parameter','content' => [],'attrib' => {'name' => 'UseFakeRegion2','default' => 'F','type' => 'logical'},'type' => 'e'},{'name' => 'parameter','content' => [],'attrib' => {'name' => 'F107Flux','default' => '150','min' => '0','type' => 'real'},'type' => 'e'},{'name' => 'parameter','content' => [],'attrib' => {'name' => 'StarLightPedConductance','default' => '0.25','min' => '0.001','type' => 'real'},'type' => 'e'},{'name' => 'parameter','content' => [],'attrib' => {'name' => 'PolarCapPedConductance','default' => '0.25','min' => '0','type' => 'real'},'type' => 'e'},{'name' => 'rule','content' => [{'content' => '
		F107Flux must be positive for ConductanceModel = 3 or 5
	','type' => 't'}],'attrib' => {'expr' => 'not ($ConductanceModel=~/3|5/ and $F107Flux==0)'},'type' => 'e'},{'content' => '
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

','type' => 't'}],'attrib' => {'name' => 'IONOSPHERE'},'type' => 'e'},{'name' => 'command','content' => [{'name' => 'parameter','content' => [{'name' => 'option','content' => [],'attrib' => {'name' => 'north','default' => 'T'},'type' => 'e'},{'name' => 'option','content' => [],'attrib' => {'name' => 'south'},'type' => 'e'},{'name' => 'option','content' => [],'attrib' => {'name' => 'cpcpmin'},'type' => 'e'},{'name' => 'option','content' => [],'attrib' => {'name' => 'average'},'type' => 'e'}],'attrib' => {'name' => 'TypeImCouple','input' => 'select','type' => 'string'},'type' => 'e'},{'content' => '

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
','type' => 't'}],'attrib' => {'name' => 'IM'},'type' => 'e'},{'name' => 'command','content' => [{'name' => 'parameter','content' => [],'attrib' => {'name' => 'DoCoupleUaCurrent','default' => 'F','type' => 'logical'},'type' => 'e'},{'name' => 'parameter','content' => [],'attrib' => {'name' => 'LatBoundary','default' => '45','if' => '$DoCoupleUaCurrent','min' => '0','type' => 'real','max' => '60'},'type' => 'e'},{'content' => '

#UA
T		DoCoupleUaCurrent
45.0		LatBoundary [deg] (only read if DoCoupleUaCurrent is true)

The DoCoupleUaCurrent parameter determines if the field aligned current
calculated by the UA component should be used. Usually the currents are
dominated by the field aligned currents provided by the GM component.
The coupling with the UA currents is still experimental.
If DoCoupleUaCurrent is set to true, the lower latitude boundary for
the potential solver should be given with the LatBoundary paramter.

The default value is DoCoupleUaCurrent=.false, i.e. the UA currents 
are not included.
','type' => 't'}],'attrib' => {'name' => 'UA'},'type' => 'e'},{'name' => 'command','content' => [{'name' => 'parameter','content' => [],'attrib' => {'name' => 'AMIEFileNorth','length' => '100','type' => 'string'},'type' => 'e'},{'name' => 'parameter','content' => [],'attrib' => {'name' => 'AMIEFileSouth','length' => '100','type' => 'string'},'type' => 'e'},{'content' => '

#AMIEFILES
IE/amie.north
IE/amie.south

Set the files to read the AMIE data from.

Default is not reading AMIE files.
','type' => 't'}],'attrib' => {'name' => 'AMIEFILES'},'type' => 'e'},{'name' => 'command','content' => [{'name' => 'parameter','content' => [],'attrib' => {'name' => 'UseSPS','default' => 'F','type' => 'logical'},'type' => 'e'},{'content' => '

#SPS
T			UseSPS

The UseSPS parameter indicates if the serial potential solver is used.
This is the default in the SWMF and it cannot be modified in the SWMF.
','type' => 't'}],'attrib' => {'name' => 'SPS'},'type' => 'e'},{'name' => 'command','content' => [{'name' => 'parameter','content' => [],'attrib' => {'name' => 'NameOfModelDir','length' => '100','type' => 'string'},'type' => 'e'},{'name' => 'parameter','content' => [],'attrib' => {'name' => 'NameOfEFieldModel','length' => '100','type' => 'string'},'type' => 'e'},{'name' => 'parameter','content' => [],'attrib' => {'name' => 'NameOfAuroralModel','length' => '100','type' => 'string'},'type' => 'e'},{'name' => 'parameter','content' => [],'attrib' => {'name' => 'NameOfSolarModel','length' => '100','type' => 'string'},'type' => 'e'},{'content' => '

#BACKGROUND
dir			NameOfModelDir
weimer96		NameOfEFieldModel
ihp			NameOfAuroralModel
xyz			NameOfSolarModel

This command cannot be used in the SWMF.
','type' => 't'}],'attrib' => {'name' => 'BACKGROUND','if' => '$_IsStandAlone'},'type' => 'e'}],'attrib' => {'name' => 'Physical parameters'},'type' => 'e'},{'name' => 'rule','content' => [{'content' => '
	Output directory IE/ionosphere should exist!
','type' => 't'}],'attrib' => {'expr' => '-d \'IE/ionosphere\' or not $_IsFirstSession'},'type' => 'e'}],'attrib' => {'name' => 'Ridley Ionosphere Model: IE Component'},'type' => 'e'}];