#^CFG FILE _FALSE_
$tree = [{'type' => 'e','attrib' => {'name' => 'Rice Convection Model 2: IM Component'},'name' => 'commandList','content' => [{'type' => 'e','attrib' => {'name' => 'Testing'},'name' => 'commandgroup','content' => [{'type' => 'e','attrib' => {'name' => 'STRICT'},'name' => 'command','content' => [{'type' => 'e','attrib' => {'type' => 'logical','name' => 'UseStrict','default' => 'T'},'name' => 'parameter','content' => []},{'type' => 't','content' => '

#STRICT
T                       UseStrict

If true then stop when parameters are incompatible. If false, try to
correct parameters and continue. Default is true, ie. strict mode.
'}]}]},{'type' => 'e','attrib' => {'name' => 'Output'},'name' => 'commandgroup','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! OUTPUT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'type' => 'e','attrib' => {'name' => 'ASCII'},'name' => 'command','content' => [{'type' => 'e','attrib' => {'type' => 'logical','name' => 'IsAscii','default' => 'T'},'name' => 'parameter','content' => []},{'type' => 't','content' => '

#ASCII
T			IsAscii

The input and output files for RCM can be either ascii or binary.
Default value is true for ascii.

'}]},{'type' => 'e','attrib' => {'name' => 'RCMDIR'},'name' => 'command','content' => [{'type' => 'e','attrib' => {'type' => 'string','length' => '100','name' => 'NameRcmDir','default' => 'IM'},'name' => 'parameter','content' => []},{'type' => 't','content' => '

#RCMDIR
IM/Plots		NameRcmDir

The NameRcmDir variable contains the name of the directory to store
output files. Default value is "IM".

'}]},{'type' => 'e','attrib' => {'name' => 'SAVEPLOT'},'name' => 'command','content' => [{'type' => 'e','attrib' => {'type' => 'integer','max' => '9','min' => '0','name' => 'nFilesPlot','default' => '0'},'name' => 'parameter','content' => []},{'type' => 'e','attrib' => {'from' => '1','to' => '$nFilesPlot'},'name' => 'for','content' => [{'type' => 'e','attrib' => {'type' => 'strings','max' => '3','min' => '3','name' => 'StringPlot'},'name' => 'parameter','content' => [{'type' => 'e','attrib' => {'input' => 'select','type' => 'string','required' => 'T','name' => 'plotarea'},'name' => 'part','content' => [{'type' => 'e','attrib' => {'value' => '2d/2D','name' => '2D vars.','default' => 'T'},'name' => 'option','content' => []},{'type' => 'e','attrib' => {'value' => '3d/2D','name' => '3D vars.'},'name' => 'option','content' => []}]},{'type' => 'e','attrib' => {'input' => 'select','type' => 'string','required' => 'T','name' => 'plotvar'},'name' => 'part','content' => [{'type' => 'e','attrib' => {'value' => 'min/MIN','name' => 'minimum set','default' => 'T'},'name' => 'option','content' => []},{'type' => 'e','attrib' => {'value' => 'max/MAX/rcm/RCM','name' => 'maximum set'},'name' => 'option','content' => []}]},{'type' => 'e','attrib' => {'input' => 'select','type' => 'string','required' => 'T','name' => 'plotformat'},'name' => 'part','content' => [{'type' => 'e','attrib' => {'value' => 'tec/TEC','name' => 'TECPLOT','default' => 'T'},'name' => 'option','content' => []},{'type' => 'e','attrib' => {'value' => 'idl/IDL','name' => 'IDL','if' => '$plotarea =~ /2d/i'},'name' => 'option','content' => []}]}]},{'type' => 'e','attrib' => {'type' => 'integer','min' => '-1','name' => 'DnSavePlot'},'name' => 'parameter','content' => []},{'type' => 'e','attrib' => {'type' => 'integer','min' => '-1','name' => 'DtSavePlot'},'name' => 'parameter','content' => []}]},{'type' => 't','content' => '
#SAVEPLOT
2			nFilesPlot
2d min idl		StringPlot
100			DnSavePlot
-1			DtSavePlot [sec]
3d max tec		StringPlot
-1			DnSavePlot
60			DtSavePlot [sec]

Default is nFilesPlot=0. \\\\

\\noindent
StringPlot must contain the following 3 parts (separated with space)
in arbitrary order:
\\begin{verbatim}
plotarea   = \'2d\', \'3d\'
plotvar    = \'min\', \'max\', \'rcm\'
plotformat = \'tec\', \'idl\'
\\end{verbatim}

\\noindent
The plotformat specifies whether to output files for processing in either
Tecplot or Idl.

The plotarea string defines the region for the plots, 2d or 3d.
Currently the 3D output is implemented for Tecplot format only.

The plotvar string specifies which variables to write, either a minimum set,
or a maximum set (the \'rcm\' string is preserved for backwards
compatibility only, and it is the same as the \'max\' string).
The \'min\' value corresponds to the variables 
\'rho T P rho(MHD) T(MHD) P(MHD\' in 2D and to \'eeta veff\' in 3D.

The \'max\' value corresponds to a lot of variables in 2D. For Tecplot format
they include the \'min\' variables, but for IDL format the \'max\' and
\'min\' variables are distinct, so you need to specify two IDL plot files
to get all the variables. In 3D only the Tecplot format is available,
and the \'max\' variable set contains \'bndloc Vm |B| V birk alam eeta veff w\'.

The DnSavePlot and DtSavePlot integers define the plotting frequency
in number of time steps and number of seconds, respectively.
A negative value -1 means that the frequency parameter is ignored.
Note that DtSavePlot must be a multiple of the time step iDtRcm
(see the #TIMESTEP command).
'}]}]},{'type' => 'e','attrib' => {'name' => 'TIME INTEGRATION'},'name' => 'commandgroup','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  TIME INTEGRATION PARAMETERS  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'type' => 'e','attrib' => {'name' => 'RESTART'},'name' => 'command','content' => [{'type' => 'e','attrib' => {'type' => 'logical','name' => 'DoRestart','default' => 'F'},'name' => 'parameter','content' => []},{'type' => 'e','attrib' => {'expr' => '-f \'IM/restartIN/RCMrestart.rst\' or not $DoRestart'},'name' => 'rule','content' => [{'type' => 't','content' => '
		File IM/restartIN/RCMrestart.rst should exist!
	'}]},{'type' => 't','content' => '

#RESTART
T			DoRestart

The DoRestart parameter determines if the RCM starts from the restart
files saved from a previous run, or from scratch via the standard input files.
The default is DoRestart = .false.
'}]},{'type' => 'e','attrib' => {'name' => 'TIMESTEP'},'name' => 'command','content' => [{'type' => 'e','attrib' => {'type' => 'integer','min' => '1','name' => 'iDtRcm','default' => '5'},'name' => 'parameter','content' => []},{'type' => 't','content' => '

#TIMESTEP
5			iDtRcm

The iDtRcm parameter defines the time step in seconds. The default
value is 5 seconds.
'}]}]},{'type' => 'e','attrib' => {'name' => 'PHYSICS PARAMETERS'},'name' => 'commandgroup','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! PHYSICS PARAMETERS  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'type' => 'e','attrib' => {'name' => 'COMPOSITION'},'name' => 'command','content' => [{'type' => 'e','attrib' => {'type' => 'real','name' => 'FractionH','default' => '0.7'},'name' => 'parameter','content' => []},{'type' => 'e','attrib' => {'type' => 'real','name' => 'FractionO','default' => '0.3'},'name' => 'parameter','content' => []},{'type' => 't','content' => '

#COMPOSITION
0.7      		FractionH
0.3      		FractionO

Fractional composition in RCM of H+ and O+.  The two need to add up to 1.0
'}]},{'type' => 'e','attrib' => {'name' => 'CHARGEEXCHANGE'},'name' => 'command','content' => [{'type' => 'e','attrib' => {'type' => 'logical','name' => 'UseChargeExchange','default' => 'T'},'name' => 'parameter','content' => []},{'type' => 'e','attrib' => {'expr' => '$UseChargeExchange'},'name' => 'if','content' => [{'type' => 'e','attrib' => {'type' => 'real','min' => '0','name' => 'SunspotNumber','default' => '125.'},'name' => 'parameter','content' => []},{'type' => 'e','attrib' => {'type' => 'real','min' => '0','name' => 'F107MonthlyMean','default' => '169.'},'name' => 'parameter','content' => []},{'type' => 'e','attrib' => {'type' => 'real','max' => '366','min' => '0','name' => 'DayOfYear','default' => '90.'},'name' => 'parameter','content' => []}]},{'type' => 't','content' => '

#CHARGEEXCHANGE
T			UseChargeExchange
125.			SunspotNumber
169.			F107MonthlyMean
90.			DayOfYear

Activate charge exchange in RCM and specify solar conditions.
'}]}]},{'type' => 'e','attrib' => {'expr' => '-d \'IM/restartOUT\' or not $_IsFirstSession'},'name' => 'rule','content' => [{'type' => 't','content' => '
	Directory IM/restartOUT should exist!
'}]}]}];