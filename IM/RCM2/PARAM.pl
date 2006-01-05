#^CFG FILE _FALSE_
$tree = [{'attrib' => {'name' => 'Rice Convection Model 2: IM Component'},'content' => [{'attrib' => {'name' => 'Testing'},'content' => [{'attrib' => {'name' => 'STRICT'},'content' => [{'attrib' => {'default' => 'T','type' => 'logical','name' => 'UseStrict'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

#STRICT
T                       UseStrict

If true then stop when parameters are incompatible. If false, try to
correct parameters and continue. Default is true, ie. strict mode.
','type' => 't'}],'type' => 'e','name' => 'command'}],'type' => 'e','name' => 'commandgroup'},{'attrib' => {'name' => 'Output'},'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! OUTPUT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
','type' => 't'},{'attrib' => {'name' => 'ASCII'},'content' => [{'attrib' => {'default' => 'T','type' => 'logical','name' => 'IsAscii'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

#ASCII
T			IsAscii

The input and output files for RCM can be either ascii or binary.
Default value is true for ascii.

','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'RCMDIR'},'content' => [{'attrib' => {'length' => '100','default' => 'IM','type' => 'string','name' => 'NameRcmDir'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

#RCMDIR
IM/Plots		NameRcmDir

The NameRcmDir variable contains the name of the directory to store
output files. Default value is "IM".

','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'SAVEPLOT'},'content' => [{'attrib' => {'min' => '0','max' => '9','default' => '0','type' => 'integer','name' => 'nFilesPlot'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'to' => '$nFilesPlot','from' => '1'},'content' => [{'attrib' => {'min' => '3','max' => '3','type' => 'strings','name' => 'StringPlot'},'content' => [{'attrib' => {'input' => 'select','required' => 'T','type' => 'string','name' => 'plotarea'},'content' => [{'attrib' => {'value' => '2d','name' => '2D vars.'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '2D','name' => '2D vars.'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '3d','name' => '3D vars.'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '3D','name' => '3D vars.'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'part'},{'attrib' => {'input' => 'select','required' => 'T','type' => 'string','name' => 'plotvar'},'content' => [{'attrib' => {'value' => 'min','name' => 'minimum variables'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'MIN','name' => 'minimum variables'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'max','name' => 'maximum variables'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'MAX','name' => 'maximum variables'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'rcm','name' => 'normal variables'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'RCM','name' => 'normal variables'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'part'},{'attrib' => {'input' => 'select','required' => 'T','type' => 'string','name' => 'plotformat'},'content' => [{'attrib' => {'value' => 'tec','name' => 'TECPLOT'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'idl','name' => 'IDL'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'part'}],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '-1','type' => 'integer','name' => 'DnSavePlot'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '-1','type' => 'integer','name' => 'DtSavePlot'},'content' => [],'type' => 'e','name' => 'parameter'}],'type' => 'e','name' => 'for'},{'content' => '
#SAVEPLOT
2			nFilesPlot
RCM 2d tec		StringPlot
100			DnSavePlot
-1			DtSavePlot
RCM 3d idl		StringPlot
-1			DnSavePlot
60			DtSavePlot

! Default is nFilesPlot=0. \\\\
!
! \\noindent
! StringPlot must contain the following 3 parts in arbitrary order
!\\begin{verbatim}
! plotvar plotarea plotformat
!
! plotarea   = \'2d\', \'3d\'
! plotvar    = \'min\', \'max\', \'rcm\'
! plotformat = \'tec\', \'idl\'
!\\end{verbatim}
!
!\\noindent
! The plotarea string defines the region for the plots, 2d or 3d.
! The plotvar string specifies which variables to write, either a minimum,
! maximum, or typical values.
! The plotformat specifies whether to output files for processing in either
! Tecplot or Idl.
! The plot_string is always followed by the plotting frequency
! DnSavePlot and DtSavePlot, whether the run it time accurate or not.\\\\
!

','type' => 't'}],'type' => 'e','name' => 'command'}],'type' => 'e','name' => 'commandgroup'},{'attrib' => {'name' => 'TIME INTEGRATION'},'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  TIME INTEGRATION PARAMETERS  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
','type' => 't'},{'attrib' => {'name' => 'RESTART'},'content' => [{'attrib' => {'default' => 'F','type' => 'logical','name' => 'DoRestart'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'expr' => '-f \'IM/restartIN/RCMrestart.rst\' or not $DoRestart'},'content' => [{'content' => '
		File IM/restartIN/RCMrestart.rst should exist!
	','type' => 't'}],'type' => 'e','name' => 'rule'},{'content' => '

#RESTART
T			DoRestart

The DoRestart parameter determines if the RCM starts from the restart
files saved from a previous run, or from scratch via the standard input files.
The default is DoRestart = .false.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'TIMESTEP'},'content' => [{'attrib' => {'min' => '1','default' => '5','type' => 'integer','name' => 'iDtRcm'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

#TIMESTEP
5			iDtRcm

The iDtRcm parameter defines the time step in seconds. The default
value is 5 seconds.
','type' => 't'}],'type' => 'e','name' => 'command'}],'type' => 'e','name' => 'commandgroup'},{'attrib' => {'name' => 'CONTROL PARAMETERS'},'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  CONTROL PARAMETERS  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
','type' => 't'},{'attrib' => {'name' => 'COMPOSITION'},'content' => [{'attrib' => {'default' => '0.7','type' => 'real','name' => 'FractionH'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '0.3','type' => 'real','name' => 'FractionO'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

#COMPOSITION
0.7      		FractionH
0.3      		FractionO

Fractional composition in RCM of H+ and O+.  The two need to add up to 1.0
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'CHARGEEXCHANGE'},'content' => [{'attrib' => {'default' => 'T','type' => 'logical','name' => 'L_dktime'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'expr' => '$L_dktime'},'content' => [{'attrib' => {'default' => '125.','type' => 'real','name' => 'SunspotNumber'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '169.','type' => 'real','name' => 'F107MonthlyMean'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '90.','type' => 'real','name' => 'DayOfYear'},'content' => [],'type' => 'e','name' => 'parameter'}],'type' => 'e','name' => 'if'},{'content' => '


#CHARGEEXCHANGE
T			L_dktime
125.			SunspotNumber
169.			F107MonthlyMean
90.			DayOfYear

Activate charge exchange in RCM and specify solar conditions.
','type' => 't'}],'type' => 'e','name' => 'command'}],'type' => 'e','name' => 'commandgroup'},{'attrib' => {'expr' => '-d \'IM/restartOUT\' or not $_IsFirstSession'},'content' => [{'content' => '
	Directory IM/restartOUT should exist!
','type' => 't'}],'type' => 'e','name' => 'rule'}],'type' => 'e','name' => 'commandList'}];