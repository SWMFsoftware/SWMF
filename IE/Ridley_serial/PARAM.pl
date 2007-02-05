#^CFG FILE _FALSE_
$tree = [{'attrib' => {'name' => 'Ridley Ionosphere Model: IE Component'},'content' => [{'attrib' => {'name' => 'Testing'},'content' => [{'attrib' => {'name' => 'STRICT'},'content' => [{'attrib' => {'default' => 'T','type' => 'logical','name' => 'UseStrict'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

#STRICT
T                       UseStrict

If true then stop when parameters are incompatible. If false, try to
correct parameters and continue. Default is true, i.e. strict mode.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'DEBUG'},'content' => [{'attrib' => {'min' => '-1','default' => '-1','type' => 'integer','name' => 'iDebugLevel'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','max' => '$_nProc','default' => '0','type' => 'integer','name' => 'iDebugProc'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

#DEBUG
2			iDebugLevel
3			iDebugProc

The iDebugLevel variable sets the level of debug information
for the processor selected by iDebugProc.

Default is iDebugLevel=-1 which is no debug info on any and iDebugProc=0.
','type' => 't'}],'type' => 'e','name' => 'command'}],'type' => 'e','name' => 'commandgroup'},{'attrib' => {'name' => 'Output'},'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! OUTPUT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
','type' => 't'},{'attrib' => {'name' => 'IONODIR'},'content' => [{'attrib' => {'length' => '100','default' => 'ionosphere','type' => 'string','name' => 'NameIonoDir'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

#IONODIR
IE/Plots		NameIonoDir

The NameIonoDir variable contains the name of the directory to store
output files. Default value is "IE/ionosphere".

','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'alias' => 'IE_SAVEPLOT','name' => 'SAVEPLOT'},'content' => [{'attrib' => {'min' => '0','max' => '10','default' => '0','type' => 'integer','name' => 'nPlotFile'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'to' => '$nPlotFile','from' => '1'},'content' => [{'attrib' => {'min' => '2','max' => '2','type' => 'strings','name' => 'StringPlot'},'content' => [{'attrib' => {'input' => 'select','required' => 'T','type' => 'string','name' => 'TypePlotForm'},'content' => [{'attrib' => {'value' => 'idl','default' => 'T','name' => 'IDL'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'tec','name' => 'TecPlot'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'part'},{'attrib' => {'input' => 'select','required' => 'T','type' => 'string','name' => 'TypePlotVar'},'content' => [{'attrib' => {'value' => 'min','default' => 'T','name' => 'Minimum'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'max','name' => 'Maximum'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'uam','name' => 'UA'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'aur','name' => 'Aurora'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'part'}],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '-1','type' => 'integer','name' => 'DnOutput'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '-1','type' => 'real','name' => 'DtOutput'},'content' => [],'type' => 'e','name' => 'parameter'}],'type' => 'e','name' => 'for'},{'content' => '

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
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'SAVELOGFILE'},'content' => [{'attrib' => {'default' => 'F','type' => 'logical','name' => 'DoSaveIELogfile'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

#SAVELOGFILE
F			DoSaveIELogfile

If true, every time that iono_solve is called, iteration, time, and soltution
information is written to a logfile.
','type' => 't'}],'type' => 'e','name' => 'command'}],'type' => 'e','name' => 'commandgroup'},{'attrib' => {'name' => 'Physical parameters'},'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!! Physical parameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

','type' => 't'},{'attrib' => {'name' => 'IONOSPHERE'},'content' => [{'attrib' => {'input' => 'select','type' => 'integer','name' => 'ConductanceModel'},'content' => [{'attrib' => {'value' => '0','name' => '(0) constant pedersen and 0 Hall'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '1','name' => '(1) constant pedersen and Hall'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '2','name' => '(2) F107 flux and and constant Hall'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '3','name' => '(3) Simple oval'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '5','default' => 'T','name' => '(5) Realistic oval'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => 'F','type' => 'logical','name' => 'UseFullCurrent'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => 'F','type' => 'logical','name' => 'UseFakeRegion2'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','default' => '150','type' => 'real','name' => 'F107Flux'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0.001','default' => '0.25','type' => 'real','name' => 'StarLightPedConductance'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','default' => '0.25','type' => 'real','name' => 'PolarCapPedConductance'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'expr' => 'not ($ConductanceModel=~/3|5/ and $F107Flux==0)'},'content' => [{'content' => '
		F107Flux must be positive for ConductanceModel = 3 or 5
	','type' => 't'}],'type' => 'e','name' => 'rule'},{'content' => '
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

','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'IM'},'content' => [{'attrib' => {'input' => 'select','type' => 'string','name' => 'TypeImCouple'},'content' => [{'attrib' => {'default' => 'T','name' => 'north'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'south'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'cpcpmin'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'average'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'parameter'},{'content' => '

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
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'UA'},'content' => [{'attrib' => {'default' => 'F','type' => 'logical','name' => 'DoCoupleUaCurrent'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','max' => '60','if' => '$DoCoupleUaCurrent','default' => '45','type' => 'real','name' => 'LatBoundary'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

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
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'AMIEFILES'},'content' => [{'attrib' => {'length' => '100','type' => 'string','name' => 'AMIEFileNorth'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'length' => '100','type' => 'string','name' => 'AMIEFileSouth'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

#AMIEFILES
IE/amie.north
IE/amie.south

Set the files to read the AMIE data from.

Default is not reading AMIE files.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'SPS'},'content' => [{'attrib' => {'default' => 'F','type' => 'logical','name' => 'UseSPS'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

#SPS
T			UseSPS

The UseSPS parameter indicates if the serial potential solver is used.
This is the default in the SWMF and it cannot be modified in the SWMF.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'if' => '$_IsStandAlone','name' => 'BACKGROUND'},'content' => [{'attrib' => {'length' => '100','type' => 'string','name' => 'NameOfModelDir'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'length' => '100','type' => 'string','name' => 'NameOfEFieldModel'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'length' => '100','type' => 'string','name' => 'NameOfAuroralModel'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'length' => '100','type' => 'string','name' => 'NameOfSolarModel'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

#BACKGROUND
dir			NameOfModelDir
weimer96		NameOfEFieldModel
ihp			NameOfAuroralModel
xyz			NameOfSolarModel

This command cannot be used in the SWMF.
','type' => 't'}],'type' => 'e','name' => 'command'}],'type' => 'e','name' => 'commandgroup'},{'attrib' => {'expr' => '-d \'IE/ionosphere\' or not $_IsFirstSession'},'content' => [{'content' => '
	Output directory IE/ionosphere should exist!
','type' => 't'}],'type' => 'e','name' => 'rule'}],'type' => 'e','name' => 'commandList'}];