#^CFG FILE _FALSE_
$tree = [{'name' => 'commandList','content' => [{'name' => 'commandgroup','content' => [{'name' => 'command','content' => [{'name' => 'parameter','content' => [],'type' => 'e','attrib' => {'name' => 'UseStrict','type' => 'logical','default' => 'T'}},{'content' => '

#STRICT
T                       UseStrict

If true then stop when parameters are incompatible. If false, try to
correct parameters and continue. Default is true, ie. strict mode.
','type' => 't'}],'type' => 'e','attrib' => {'name' => 'STRICT'}}],'type' => 'e','attrib' => {'name' => 'Testing'}},{'name' => 'commandgroup','content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! OUTPUT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
','type' => 't'},{'name' => 'command','content' => [{'name' => 'parameter','content' => [],'type' => 'e','attrib' => {'name' => 'IsAscii','type' => 'logical','default' => 'T'}},{'content' => '

#ASCII
T			IsAscii

The input and output files for RCM can be either ascii or binary.
Default value is true for ascii.

','type' => 't'}],'type' => 'e','attrib' => {'name' => 'ASCII'}},{'name' => 'command','content' => [{'name' => 'parameter','content' => [],'type' => 'e','attrib' => {'length' => '100','name' => 'NameRcmDir','type' => 'string','default' => 'IM'}},{'content' => '

#RCMDIR
IM/Plots		NameRcmDir

The NameRcmDir variable contains the name of the directory to store
output files. Default value is "IM".

','type' => 't'}],'type' => 'e','attrib' => {'name' => 'RCMDIR'}},{'name' => 'command','content' => [{'name' => 'parameter','content' => [],'type' => 'e','attrib' => {'name' => 'nFilesPlot','type' => 'integer','default' => '0','min' => '0','max' => '9'}},{'name' => 'for','content' => [{'name' => 'parameter','content' => [{'name' => 'part','content' => [{'name' => 'option','content' => [],'type' => 'e','attrib' => {'name' => '2D vars.','value' => '2d'}},{'name' => 'option','content' => [],'type' => 'e','attrib' => {'name' => '2D vars.','value' => '2D'}},{'name' => 'option','content' => [],'type' => 'e','attrib' => {'name' => '3D vars.','value' => '3d'}},{'name' => 'option','content' => [],'type' => 'e','attrib' => {'name' => '3D vars.','value' => '3D'}}],'type' => 'e','attrib' => {'input' => 'select','name' => 'plotarea','type' => 'string','required' => 'T'}},{'name' => 'part','content' => [{'name' => 'option','content' => [],'type' => 'e','attrib' => {'name' => 'minimum variables','value' => 'min'}},{'name' => 'option','content' => [],'type' => 'e','attrib' => {'name' => 'minimum variables','value' => 'MIN'}},{'name' => 'option','content' => [],'type' => 'e','attrib' => {'name' => 'maximum variables','value' => 'max'}},{'name' => 'option','content' => [],'type' => 'e','attrib' => {'name' => 'maximum variables','value' => 'MAX'}},{'name' => 'option','content' => [],'type' => 'e','attrib' => {'name' => 'normal variables','value' => 'rcm'}},{'name' => 'option','content' => [],'type' => 'e','attrib' => {'name' => 'normal variables','value' => 'RCM'}}],'type' => 'e','attrib' => {'input' => 'select','name' => 'plotvar','type' => 'string','required' => 'T'}},{'name' => 'part','content' => [{'name' => 'option','content' => [],'type' => 'e','attrib' => {'name' => 'TECPLOT','value' => 'tec'}},{'name' => 'option','content' => [],'type' => 'e','attrib' => {'name' => 'IDL','value' => 'idl'}}],'type' => 'e','attrib' => {'input' => 'select','name' => 'plotformat','type' => 'string','required' => 'T'}}],'type' => 'e','attrib' => {'name' => 'StringPlot','type' => 'strings','min' => '3','max' => '3'}},{'name' => 'parameter','content' => [],'type' => 'e','attrib' => {'name' => 'DnSavePlot','type' => 'integer','min' => '-1'}},{'name' => 'parameter','content' => [],'type' => 'e','attrib' => {'name' => 'DtSavePlot','type' => 'integer','min' => '-1'}}],'type' => 'e','attrib' => {'from' => '1','to' => '$nFilesPlot'}},{'content' => '
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

','type' => 't'}],'type' => 'e','attrib' => {'name' => 'SAVEPLOT'}}],'type' => 'e','attrib' => {'name' => 'Output'}},{'name' => 'commandgroup','content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  TIME INTEGRATION PARAMETERS  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
','type' => 't'},{'name' => 'command','content' => [{'name' => 'parameter','content' => [],'type' => 'e','attrib' => {'name' => 'DoRestart','type' => 'logical','default' => 'F'}},{'name' => 'rule','content' => [{'content' => '
		File IM/restartIN/RCMrestart.rst should exist!
	','type' => 't'}],'type' => 'e','attrib' => {'expr' => '-f \'IM/restartIN/RCMrestart.rst\' or not $DoRestart'}},{'content' => '

#RESTART
T			DoRestart

The DoRestart parameter determines if the RCM starts from the restart
files saved from a previous run, or from scratch via the standard input files.
The default is DoRestart = .false.
','type' => 't'}],'type' => 'e','attrib' => {'name' => 'RESTART'}},{'name' => 'command','content' => [{'name' => 'parameter','content' => [],'type' => 'e','attrib' => {'name' => 'iDtRcm','type' => 'integer','default' => '5','min' => '1'}},{'content' => '

#TIMESTEP
5			iDtRcm

The iDtRcm parameter defines the time step in seconds. The default
value is 5 seconds.
','type' => 't'}],'type' => 'e','attrib' => {'name' => 'TIMESTEP'}}],'type' => 'e','attrib' => {'name' => 'TIME INTEGRATION'}},{'name' => 'commandgroup','content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! PHYSICS PARAMETERS  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
','type' => 't'},{'name' => 'command','content' => [{'name' => 'parameter','content' => [],'type' => 'e','attrib' => {'name' => 'FractionH','type' => 'real','default' => '0.7'}},{'name' => 'parameter','content' => [],'type' => 'e','attrib' => {'name' => 'FractionO','type' => 'real','default' => '0.3'}},{'content' => '

#COMPOSITION
0.7      		FractionH
0.3      		FractionO

Fractional composition in RCM of H+ and O+.  The two need to add up to 1.0
','type' => 't'}],'type' => 'e','attrib' => {'name' => 'COMPOSITION'}},{'name' => 'command','content' => [{'name' => 'parameter','content' => [],'type' => 'e','attrib' => {'name' => 'UseChargeExchange','type' => 'logical','default' => 'T'}},{'name' => 'if','content' => [{'name' => 'parameter','content' => [],'type' => 'e','attrib' => {'name' => 'SunspotNumber','type' => 'real','default' => '125.','min' => '0'}},{'name' => 'parameter','content' => [],'type' => 'e','attrib' => {'name' => 'F107MonthlyMean','type' => 'real','default' => '169.','min' => '0'}},{'name' => 'parameter','content' => [],'type' => 'e','attrib' => {'name' => 'DayOfYear','type' => 'real','default' => '90.','min' => '0','max' => '366'}}],'type' => 'e','attrib' => {'expr' => '$UseChargeExchange'}},{'content' => '

#CHARGEEXCHANGE
T			UseChargeExchange
125.			SunspotNumber
169.			F107MonthlyMean
90.			DayOfYear

Activate charge exchange in RCM and specify solar conditions.
','type' => 't'}],'type' => 'e','attrib' => {'name' => 'CHARGEEXCHANGE'}}],'type' => 'e','attrib' => {'name' => 'PHYSICS PARAMETERS'}},{'name' => 'rule','content' => [{'content' => '
	Directory IM/restartOUT should exist!
','type' => 't'}],'type' => 'e','attrib' => {'expr' => '-d \'IM/restartOUT\' or not $_IsFirstSession'}}],'type' => 'e','attrib' => {'name' => 'Rice Convection Model 2: IM Component'}}];