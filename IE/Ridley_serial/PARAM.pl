#^CFG FILE _FALSE_
$tree = [{'attrib' => {'name' => 'Ridley Ionosphere Model: IE Component'},'content' => [{'attrib' => {'name' => 'Testing'},'content' => [{'attrib' => {'name' => 'STRICT'},'content' => [{'attrib' => {'default' => 'T','type' => 'logical','name' => 'UseStrict'},'content' => [],'name' => 'parameter','type' => 'e'},{'content' => '

#STRICT
T                       UseStrict

If true then stop when parameters are incompatible. If false, try to
correct parameters and continue. Default is true, ie. strict mode.
','type' => 't'}],'name' => 'command','type' => 'e'},{'attrib' => {'name' => 'DEBUG'},'content' => [{'attrib' => {'min' => '-1','default' => '-1','type' => 'integer','name' => 'iDebugLevel'},'content' => [],'name' => 'parameter','type' => 'e'},{'attrib' => {'min' => '0','default' => '0','max' => '$_nProc','type' => 'integer','name' => 'iDebugProc'},'content' => [],'name' => 'parameter','type' => 'e'},{'content' => '

#DEBUG
2			iDebugLevel
3			iDebugProc

The iDebugLevel variable sets the level of debug information
for the processor selected by iDebugProc.

Default is iDebugLevel=-1 which is no debug info on any and iDebugProc=0.
','type' => 't'}],'name' => 'command','type' => 'e'}],'name' => 'commandgroup','type' => 'e'},{'attrib' => {'name' => 'Output'},'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! OUTPUT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
','type' => 't'},{'attrib' => {'name' => 'IONODIR'},'content' => [{'attrib' => {'default' => 'ionosphere','length' => '100','type' => 'string','name' => 'NameIonoDir'},'content' => [],'name' => 'parameter','type' => 'e'},{'content' => '

#IONODIR
IE/Plots		NameIonoDir

The NameIonoDir variable contains the name of the directory to store
ouput files. Default value is "IE/ionosphere".

','type' => 't'}],'name' => 'command','type' => 'e'},{'attrib' => {'alias' => 'IE_SAVEPLOT','name' => 'SAVEPLOT'},'content' => [{'attrib' => {'min' => '0','default' => '0','max' => '10','type' => 'integer','name' => 'nPlotFile'},'content' => [],'name' => 'parameter','type' => 'e'},{'attrib' => {'from' => '1','to' => '$nPlotFile'},'content' => [{'attrib' => {'min' => '2','max' => '2','type' => 'strings','name' => 'StringPlot'},'content' => [{'attrib' => {'input' => 'select','required' => 'T','type' => 'string','name' => 'TypePlotForm'},'content' => [{'attrib' => {'default' => 'T','value' => 'idl','name' => 'IDL'},'content' => [],'name' => 'option','type' => 'e'},{'attrib' => {'value' => 'tec','name' => 'TecPlot'},'content' => [],'name' => 'option','type' => 'e'}],'name' => 'part','type' => 'e'},{'attrib' => {'input' => 'select','required' => 'T','type' => 'string','name' => 'TypePlotVar'},'content' => [{'attrib' => {'default' => 'T','value' => 'min','name' => 'Minimum'},'content' => [],'name' => 'option','type' => 'e'},{'attrib' => {'value' => 'max','name' => 'Maximum'},'content' => [],'name' => 'option','type' => 'e'},{'attrib' => {'value' => 'aur','name' => 'Aurora'},'content' => [],'name' => 'option','type' => 'e'}],'name' => 'part','type' => 'e'}],'name' => 'parameter','type' => 'e'},{'attrib' => {'min' => '-1','type' => 'integer','name' => 'DnOutput'},'content' => [],'name' => 'parameter','type' => 'e'},{'attrib' => {'min' => '-1','type' => 'real','name' => 'DtOutput'},'content' => [],'name' => 'parameter','type' => 'e'}],'name' => 'for','type' => 'e'},{'content' => '

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
','type' => 't'}],'name' => 'command','type' => 'e'}],'name' => 'commandgroup','type' => 'e'},{'attrib' => {'name' => 'Physical parameters'},'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!! Physical parameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

','type' => 't'},{'attrib' => {'name' => 'IONOSPHERE'},'content' => [{'attrib' => {'input' => 'select','type' => 'integer','name' => 'ConductanceModel'},'content' => [{'attrib' => {'value' => '0','name' => '(0) constant pedersen and 0 Hall'},'content' => [],'name' => 'option','type' => 'e'},{'attrib' => {'value' => '1','name' => '(1) constant pedersen and Hall'},'content' => [],'name' => 'option','type' => 'e'},{'attrib' => {'value' => '2','name' => '(2) F107 flux and and constant Hall'},'content' => [],'name' => 'option','type' => 'e'},{'attrib' => {'value' => '3','name' => '(3) Simple oval'},'content' => [],'name' => 'option','type' => 'e'},{'attrib' => {'default' => 'T','value' => '5','name' => '(5) Realistic oval'},'content' => [],'name' => 'option','type' => 'e'}],'name' => 'parameter','type' => 'e'},{'attrib' => {'default' => 'F','type' => 'logical','name' => 'UseFullCurrent'},'content' => [],'name' => 'parameter','type' => 'e'},{'attrib' => {'default' => 'F','type' => 'logical','name' => 'UseFakeRegion2'},'content' => [],'name' => 'parameter','type' => 'e'},{'attrib' => {'min' => '0','default' => '150','type' => 'real','name' => 'F107Flux'},'content' => [],'name' => 'parameter','type' => 'e'},{'attrib' => {'min' => '0.001','default' => '0.25','type' => 'real','name' => 'StarLightPedConductance'},'content' => [],'name' => 'parameter','type' => 'e'},{'attrib' => {'min' => '0','default' => '0.25','type' => 'real','name' => 'PolarCapPedConductance'},'content' => [],'name' => 'parameter','type' => 'e'},{'attrib' => {'expr' => 'not ($ConductanceModel=~/3|5/ and $F107Flux==0)'},'content' => [{'content' => '
		F107Flux must be positive for ConductanceModel = 3 or 5
	','type' => 't'}],'name' => 'rule','type' => 'e'},{'content' => '

The ConductanceModel variable determines which ionosphere model is used:
  0 - uses a constant Pedersen conductance which is set by 
      StarLightPedConductance
  1 - uses a constant Pedersen conductance which is set by 
      StarLightPedConductance, and a constant Hall conductance
      which is set by PolarCapPedConductance
  2 - uses a solar EUV combined with a nightside conductance, so
      it uses f107_flux and StarLightPedConductance
  3 - uses solar EUV, nightside, and crude oval, so uses
      f107_flux, StarLightPedConductance, and PolarCapPedConductance,
      since a polar cap is defined with the oval.
  4 - doesn\'t work
  5 - more realistic oval, uses same variables as 3.
','type' => 't'}],'name' => 'command','type' => 'e'},{'attrib' => {'name' => 'IM'},'content' => [{'attrib' => {'input' => 'select','type' => 'string','name' => 'TypeImCouple'},'content' => [{'attrib' => {'default' => 'T','name' => 'north'},'content' => [],'name' => 'option','type' => 'e'},{'attrib' => {'name' => 'south'},'content' => [],'name' => 'option','type' => 'e'},{'attrib' => {'name' => 'cpcpmin'},'content' => [],'name' => 'option','type' => 'e'},{'attrib' => {'name' => 'average'},'content' => [],'name' => 'option','type' => 'e'}],'name' => 'parameter','type' => 'e'},{'content' => '

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
','type' => 't'}],'name' => 'command','type' => 'e'},{'attrib' => {'name' => 'AMIEFILES'},'content' => [{'attrib' => {'length' => '100','type' => 'string','name' => 'AMIEFileNorth'},'content' => [],'name' => 'parameter','type' => 'e'},{'attrib' => {'length' => '100','type' => 'string','name' => 'AMIEFileSouth'},'content' => [],'name' => 'parameter','type' => 'e'},{'content' => '

#AMIEFILES
IE/amie.north
IE/amie.south

Set the files to read the AMIE data from.

Default is not reading AMIE files.
','type' => 't'}],'name' => 'command','type' => 'e'},{'attrib' => {'name' => 'SPS'},'content' => [{'attrib' => {'default' => 'F','type' => 'logical','name' => 'UseSPS'},'content' => [],'name' => 'parameter','type' => 'e'},{'content' => '

#SPS
T			UseSPS

The UseSPS parameter indicates if the serial potential solver is used.
','type' => 't'}],'name' => 'command','type' => 'e'},{'attrib' => {'name' => 'BACKGROUND'},'content' => [{'attrib' => {'length' => '100','type' => 'string','name' => 'NameOfModelDir'},'content' => [],'name' => 'parameter','type' => 'e'},{'attrib' => {'length' => '100','type' => 'string','name' => 'NameOfFieldModel'},'content' => [],'name' => 'parameter','type' => 'e'},{'attrib' => {'length' => '100','type' => 'string','name' => 'NameOfAuroralModel'},'content' => [],'name' => 'parameter','type' => 'e'},{'attrib' => {'length' => '100','type' => 'string','name' => 'NameOfSolarModel'},'content' => [],'name' => 'parameter','type' => 'e'},{'content' => '

#BACKGROUND
???

???
','type' => 't'}],'name' => 'command','type' => 'e'}],'name' => 'commandgroup','type' => 'e'},{'attrib' => {'expr' => '-d \'IE/ionosphere\' or not $_IsFirstSession'},'content' => [{'content' => '
	Output directory IE/ionosphere should exist!
','type' => 't'}],'name' => 'rule','type' => 'e'}],'name' => 'commandList','type' => 'e'}];