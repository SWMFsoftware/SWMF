#^CFG FILE _FALSE_
$tree = [{'attrib' => {'name' => 'Ridley Ionosphere Model: IE Component'},'content' => [{'attrib' => {'name' => 'Testing'},'content' => [{'attrib' => {'name' => 'STRICT'},'content' => [{'attrib' => {'type' => 'logical','default' => 'T','name' => 'UseStrict'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

#STRICT
T                       UseStrict

If true then stop when parameters are incompatible. If false, try to
correct parameters and continue. Default is true, ie. strict mode.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'DEBUG'},'content' => [{'attrib' => {'min' => '-1','type' => 'integer','default' => '-1','name' => 'iDebugLevel'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','max' => '$_nProc','type' => 'integer','default' => '0','name' => 'iDebugProc'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

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
','type' => 't'},{'attrib' => {'name' => 'IONODIR'},'content' => [{'attrib' => {'length' => '100','type' => 'string','default' => 'ionosphere','name' => 'NameIonoDir'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

#IONODIR
IE/Plots		NameIonoDir

The NameIonoDir variable contains the name of the directory to store
ouput files. Default value is "IE/ionosphere".

','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'alias' => 'IE_SAVEPLOT','name' => 'SAVEPLOT'},'content' => [{'attrib' => {'min' => '0','max' => '10','type' => 'integer','default' => '0','name' => 'nPlotFile'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'to' => '$nPlotFile','from' => '1'},'content' => [{'attrib' => {'min' => '2','max' => '2','type' => 'strings','name' => 'StringPlot'},'content' => [{'attrib' => {'required' => 'T','input' => 'select','type' => 'string','name' => 'TypePlotForm'},'content' => [{'attrib' => {'default' => 'T','value' => 'idl','name' => 'IDL'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'tec','name' => 'TecPlot'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'part'},{'attrib' => {'required' => 'T','input' => 'select','type' => 'string','name' => 'TypePlotVar'},'content' => [{'attrib' => {'default' => 'T','value' => 'min','name' => 'Minimum'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'max','name' => 'Maximum'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'aur','name' => 'Aurora'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'part'}],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '-1','type' => 'integer','name' => 'DnOutput'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '-1','type' => 'real','name' => 'DtOutput'},'content' => [],'type' => 'e','name' => 'parameter'}],'type' => 'e','name' => 'for'},{'content' => '

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
','type' => 't'}],'type' => 'e','name' => 'command'}],'type' => 'e','name' => 'commandgroup'},{'attrib' => {'name' => 'Physical parameters'},'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!! Physical parameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

','type' => 't'},{'attrib' => {'name' => 'IONOSPHERE'},'content' => [{'attrib' => {'input' => 'select','type' => 'integer','name' => 'ConductanceModel'},'content' => [{'attrib' => {'value' => '0','name' => '(0) constant pedersen and 0 Hall'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '1','name' => '(1) constant pedersen and Hall'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '2','name' => '(2) F107 flux and and constant Hall'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '3','name' => '(3) Simple oval'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'default' => 'T','value' => '5','name' => '(5) Realistic oval'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'parameter'},{'attrib' => {'type' => 'logical','default' => 'F','name' => 'UseFullCurrent'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'type' => 'logical','default' => 'F','name' => 'UseFakeRegion2'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','type' => 'real','default' => '150','name' => 'F107Flux'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0.001','type' => 'real','default' => '0.25','name' => 'StarLightPedConductance'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','type' => 'real','default' => '0.25','name' => 'PolarCapPedConductance'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'expr' => 'not ($ConductanceModel=~/3|5/ and $F107Flux==0)'},'content' => [{'content' => '
		F107Flux must be positive for ConductanceModel = 3 or 5
	','type' => 't'}],'type' => 'e','name' => 'rule'},{'content' => '

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
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'AMIEFILES'},'content' => [{'attrib' => {'length' => '100','type' => 'string','name' => 'AMIEFileNorth'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'length' => '100','type' => 'string','name' => 'AMIEFileSouth'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

#AMIEFILES
IE/amie.north
IE/amie.south

Set the files to read the AMIE data from.

Default is not reading AMIE files.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'SPS'},'content' => [{'attrib' => {'type' => 'logical','default' => 'F','name' => 'UseSPS'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

#SPS
T			UseSPS

The UseSPS parameter indicates if the serial potential solver is used.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'BACKGROUND'},'content' => [{'attrib' => {'length' => '100','type' => 'string','name' => 'NameOfModelDir'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'length' => '100','type' => 'string','name' => 'NameOfFieldModel'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'length' => '100','type' => 'string','name' => 'NameOfAuroralModel'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'length' => '100','type' => 'string','name' => 'NameOfSolarModel'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

#BACKGROUND
???

???
','type' => 't'}],'type' => 'e','name' => 'command'}],'type' => 'e','name' => 'commandgroup'}],'type' => 'e','name' => 'commandList'}];