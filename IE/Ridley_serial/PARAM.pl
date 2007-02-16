#^CFG FILE _FALSE_
$tree = [{'content' => [{'content' => [{'content' => [{'content' => [],'name' => 'parameter','attrib' => {'name' => 'UseStrict','type' => 'logical','default' => 'T'},'type' => 'e'},{'content' => '

#STRICT
T                       UseStrict

If true then stop when parameters are incompatible. If false, try to
correct parameters and continue. Default is true, i.e. strict mode.
','type' => 't'}],'name' => 'command','attrib' => {'name' => 'STRICT'},'type' => 'e'},{'content' => [{'content' => [],'name' => 'parameter','attrib' => {'name' => 'iDebugLevel','type' => 'integer','default' => '-1','min' => '-1'},'type' => 'e'},{'content' => [],'name' => 'parameter','attrib' => {'name' => 'iDebugProc','max' => '$_nProc','type' => 'integer','default' => '0','min' => '0'},'type' => 'e'},{'content' => '

#DEBUG
2			iDebugLevel
3			iDebugProc

The iDebugLevel variable sets the level of debug information
for the processor selected by iDebugProc.

Default is iDebugLevel=-1 which is no debug info on any and iDebugProc=0.
','type' => 't'}],'name' => 'command','attrib' => {'name' => 'DEBUG'},'type' => 'e'}],'name' => 'commandgroup','attrib' => {'name' => 'Testing'},'type' => 'e'},{'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! OUTPUT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
','type' => 't'},{'content' => [{'content' => [],'name' => 'parameter','attrib' => {'name' => 'NameIonoDir','length' => '100','type' => 'string','default' => 'ionosphere'},'type' => 'e'},{'content' => '

#IONODIR
IE/Plots		NameIonoDir

The NameIonoDir variable contains the name of the directory to store
output files. Default value is "IE/ionosphere".

','type' => 't'}],'name' => 'command','attrib' => {'name' => 'IONODIR'},'type' => 'e'},{'content' => [{'content' => [],'name' => 'parameter','attrib' => {'name' => 'nPlotFile','max' => '10','type' => 'integer','default' => '0','min' => '0'},'type' => 'e'},{'content' => [{'content' => [{'content' => [{'content' => [],'name' => 'option','attrib' => {'name' => 'IDL','value' => 'idl','default' => 'T'},'type' => 'e'},{'content' => [],'name' => 'option','attrib' => {'name' => 'TecPlot','value' => 'tec'},'type' => 'e'}],'name' => 'part','attrib' => {'name' => 'TypePlotForm','input' => 'select','type' => 'string','required' => 'T'},'type' => 'e'},{'content' => [{'content' => [],'name' => 'option','attrib' => {'name' => 'Minimum','value' => 'min','default' => 'T'},'type' => 'e'},{'content' => [],'name' => 'option','attrib' => {'name' => 'Maximum','value' => 'max'},'type' => 'e'},{'content' => [],'name' => 'option','attrib' => {'name' => 'UA','value' => 'uam'},'type' => 'e'},{'content' => [],'name' => 'option','attrib' => {'name' => 'Aurora','value' => 'aur'},'type' => 'e'}],'name' => 'part','attrib' => {'name' => 'TypePlotVar','input' => 'select','type' => 'string','required' => 'T'},'type' => 'e'}],'name' => 'parameter','attrib' => {'name' => 'StringPlot','max' => '2','type' => 'strings','min' => '2'},'type' => 'e'},{'content' => [],'name' => 'parameter','attrib' => {'name' => 'DnOutput','type' => 'integer','min' => '-1'},'type' => 'e'},{'content' => [],'name' => 'parameter','attrib' => {'name' => 'DtOutput','type' => 'real','min' => '-1'},'type' => 'e'}],'name' => 'for','attrib' => {'to' => '$nPlotFile','from' => '1'},'type' => 'e'},{'content' => '

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
','type' => 't'}],'name' => 'command','attrib' => {'name' => 'SAVEPLOT','alias' => 'IE_SAVEPLOT'},'type' => 'e'},{'content' => [{'content' => [],'name' => 'parameter','attrib' => {'name' => 'DoSaveIELogfile','type' => 'logical','default' => 'F'},'type' => 'e'},{'content' => '

#SAVELOGFILE
F			DoSaveIELogfile

If true, every time that iono_solve is called, iteration, time, and soltution
information is written to a logfile.
','type' => 't'}],'name' => 'command','attrib' => {'name' => 'SAVELOGFILE'},'type' => 'e'}],'name' => 'commandgroup','attrib' => {'name' => 'Output'},'type' => 'e'},{'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!! Physical parameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

','type' => 't'},{'content' => [{'content' => [{'content' => [],'name' => 'option','attrib' => {'name' => '(0) constant Pedersen and 0 Hall','value' => '0'},'type' => 'e'},{'content' => [],'name' => 'option','attrib' => {'name' => '(1) constant Pedersen and Hall','value' => '1'},'type' => 'e'},{'content' => [],'name' => 'option','attrib' => {'name' => '(2) F107 flux and and constant Hall','value' => '2'},'type' => 'e'},{'content' => [],'name' => 'option','attrib' => {'name' => '(3) Simple oval','value' => '3'},'type' => 'e'},{'content' => [],'name' => 'option','attrib' => {'name' => '(4) Restricted oval','value' => '4'},'type' => 'e'},{'content' => [],'name' => 'option','attrib' => {'name' => '(5) Realistic oval','value' => '5','default' => 'T'},'type' => 'e'}],'name' => 'parameter','attrib' => {'name' => 'iConductanceModel','input' => 'select','type' => 'integer'},'type' => 'e'},{'content' => [],'name' => 'parameter','attrib' => {'name' => 'UseFullCurrent','type' => 'logical','default' => 'F'},'type' => 'e'},{'content' => [],'name' => 'parameter','attrib' => {'name' => 'UseFakeRegion2','type' => 'logical','default' => 'F'},'type' => 'e'},{'content' => [],'name' => 'parameter','attrib' => {'name' => 'F107Flux','type' => 'real','default' => '150','min' => '0'},'type' => 'e'},{'content' => [],'name' => 'parameter','attrib' => {'name' => 'StarLightPedConductance','type' => 'real','default' => '0.25','min' => '0.001'},'type' => 'e'},{'content' => [],'name' => 'parameter','attrib' => {'name' => 'PolarCapPedConductance','type' => 'real','default' => '0.25','min' => '0'},'type' => 'e'},{'content' => [{'content' => '
		F107Flux must be positive for ConductanceModel = 3, 4 or 5
	','type' => 't'}],'name' => 'rule','attrib' => {'expr' => 'not ($iConductanceModel=~/3|4|5/ and $F107Flux==0)'},'type' => 'e'},{'content' => '
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
  4 - restricted oval, uses same variables as 3.\\\\
  5 - more realistic oval, uses same variables as 3.\\\\

Model 4 and 5 differ in the way the conductances are limited to the
fitted oval. Model 4 is more restrictive, while model 5 is somewhat
more relaxed.

The UseFullCurrent and UseFakeRegion2 logicals were used in the past
to test various algorithmic choices. They should have a false value now.
The default values are shown by the example above.

','type' => 't'}],'name' => 'command','attrib' => {'name' => 'IONOSPHERE'},'type' => 'e'},{'content' => [{'content' => [{'content' => [],'name' => 'option','attrib' => {'name' => 'north','default' => 'T'},'type' => 'e'},{'content' => [],'name' => 'option','attrib' => {'name' => 'south'},'type' => 'e'},{'content' => [],'name' => 'option','attrib' => {'name' => 'cpcpmin'},'type' => 'e'},{'content' => [],'name' => 'option','attrib' => {'name' => 'average'},'type' => 'e'}],'name' => 'parameter','attrib' => {'name' => 'TypeImCouple','input' => 'select','type' => 'string'},'type' => 'e'},{'content' => '

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
','type' => 't'}],'name' => 'command','attrib' => {'name' => 'IM'},'type' => 'e'},{'content' => [{'content' => [],'name' => 'parameter','attrib' => {'name' => 'DoCoupleUaCurrent','type' => 'logical','default' => 'F'},'type' => 'e'},{'content' => [],'name' => 'parameter','attrib' => {'name' => 'LatBoundary','max' => '60','if' => '$DoCoupleUaCurrent','type' => 'real','default' => '45','min' => '0'},'type' => 'e'},{'content' => '

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
','type' => 't'}],'name' => 'command','attrib' => {'name' => 'UA'},'type' => 'e'},{'content' => [{'content' => [],'name' => 'parameter','attrib' => {'name' => 'AMIEFileNorth','length' => '100','type' => 'string'},'type' => 'e'},{'content' => [],'name' => 'parameter','attrib' => {'name' => 'AMIEFileSouth','length' => '100','type' => 'string'},'type' => 'e'},{'content' => '

#AMIEFILES
IE/amie.north
IE/amie.south

Set the files to read the AMIE data from.

Default is not reading AMIE files.
','type' => 't'}],'name' => 'command','attrib' => {'name' => 'AMIEFILES'},'type' => 'e'},{'content' => [{'content' => [],'name' => 'parameter','attrib' => {'name' => 'UseSPS','type' => 'logical','default' => 'F'},'type' => 'e'},{'content' => '

#SPS
T			UseSPS

The UseSPS parameter indicates if the serial potential solver is used.
This is the default in the SWMF and it cannot be modified in the SWMF.
','type' => 't'}],'name' => 'command','attrib' => {'name' => 'SPS'},'type' => 'e'},{'content' => [{'content' => [],'name' => 'parameter','attrib' => {'name' => 'NameOfModelDir','length' => '100','type' => 'string'},'type' => 'e'},{'content' => [],'name' => 'parameter','attrib' => {'name' => 'NameOfEFieldModel','length' => '100','type' => 'string'},'type' => 'e'},{'content' => [],'name' => 'parameter','attrib' => {'name' => 'NameOfAuroralModel','length' => '100','type' => 'string'},'type' => 'e'},{'content' => [],'name' => 'parameter','attrib' => {'name' => 'NameOfSolarModel','length' => '100','type' => 'string'},'type' => 'e'},{'content' => '

#BACKGROUND
dir			NameOfModelDir
weimer96		NameOfEFieldModel
ihp			NameOfAuroralModel
xyz			NameOfSolarModel

This command cannot be used in the SWMF.
','type' => 't'}],'name' => 'command','attrib' => {'name' => 'BACKGROUND','if' => '$_IsStandAlone'},'type' => 'e'},{'content' => [{'content' => [],'name' => 'parameter','attrib' => {'name' => 'OvalWidthFactor','type' => 'real','default' => '1.0','min' => '0'},'type' => 'e'},{'content' => [],'name' => 'parameter','attrib' => {'name' => 'OvalStrengthFactor','type' => 'real','default' => '1.7','min' => '0'},'type' => 'e'},{'content' => '

#CONDUCTANCE
1.0			OvalWidthFactor
1.7			OvalStrengthFactor

Modifies the conductance by adjusting the oval width/strength.
Only works for conductance model 4!
','type' => 't'}],'name' => 'command','attrib' => {'name' => 'CONDUCTANCE'},'type' => 'e'}],'name' => 'commandgroup','attrib' => {'name' => 'Physical parameters'},'type' => 'e'},{'content' => [{'content' => [{'content' => [],'name' => 'parameter','attrib' => {'name' => 'UsePreconditioner','type' => 'logical','default' => 'T'},'type' => 'e'},{'content' => [],'name' => 'parameter','attrib' => {'name' => 'UseInitialGuess','type' => 'logical','default' => 'T'},'type' => 'e'},{'content' => [],'name' => 'parameter','attrib' => {'name' => 'Tolerance','type' => 'real','default' => '0.01','min' => '0'},'type' => 'e'},{'content' => [],'name' => 'parameter','attrib' => {'name' => 'MaxIteration','type' => 'integer','default' => '100','min' => '1'},'type' => 'e'},{'content' => '

#KRYLOV
T			UsePreconditioner
T			UseInitialGuess
0.01			Tolerance
100			MaxIteration

This command controls the parameters for the Krylov solver used to
solve the Poisson type equation for the electric potential.
If UsePreconditioner is true the solver uses a preconditioner.
If UseInitialGuess is true, the previous solution is used as an 
initial guess. The Tolerance parameter sets the second norm of
the final (preconditioned) residual. The MaxIteration parameter sets the
maximum number of iterations before the linear solver gives up.
In most cases the default values should work fine.

The default values are shown above.
','type' => 't'}],'name' => 'command','attrib' => {'name' => 'KRYLOV'},'type' => 'e'}],'name' => 'commandgroup','attrib' => {'name' => 'Scheme parameters'},'type' => 'e'},{'content' => [{'content' => '
	Output directory IE/ionosphere should exist!
','type' => 't'}],'name' => 'rule','attrib' => {'expr' => '-d \'IE/ionosphere\' or not $_IsFirstSession'},'type' => 'e'}],'name' => 'commandList','attrib' => {'name' => 'Ridley Ionosphere Model: IE Component'},'type' => 'e'}];