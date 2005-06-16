#^CFG FILE _FALSE_
$tree = [{'type' => 'e','content' => [{'type' => 't','content' => '

CON reads input parameters from the PARAM.in file and the files
included into PARAM.in. All commands interpreted by CON start
with the # character followed by capital letters and numbers.
The commands can have an arbitrary number of parameters, 
which are written into the lines following the command. 
Other lines are ignored, and can be used for remarks.
The general format of the parameter file is

#COMMANDNAME1
variable1
variable2

#COMMANDNAME2

#COMMANDNAME3
variable3

The #BEGIN_COMP ID and #END_COMP ID commands are exceptional in the sense 
that their parameters are written in the same line as the command itself.
This exception makes the parameter file more readable. The parameter
is the two-character component ID. There must be exactly one space
between the #BEGIN_COMP or #END_COMP string and the ID.
The lines between the #BEGIN_COMP ID and the matching 
#END_COMP ID commands are passed to the component with the
corresponding ID. For example 

#BEGIN_COMP GM

#AMR
-1

#END_COMP GM

The parameters passed to the components can be of arbitrary format.
The only restriction is that the length of the lines cannot exceed
100 characters (extra characters will be ignored).


'},{'type' => 'e','content' => [],'name' => 'set','attrib' => {'type' => 'logical','value' => 'T','name' => 'DoTimeAccurate'}},{'type' => 'e','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!! GENERAL COMMANDS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'string','length' => '100','name' => 'NameIncludeFile'}},{'type' => 't','content' => '

#INCLUDE
run/GM/restartIN/restart.H		NameIncludeFile

The NameIncludeFile parameter contains the name of the file to be included.
The file name may be followed with a trailing comment if it
is separated with at least 3 spaces or one TAB character.
The #INCLUDE command can be used anywhere in the parameter file,
even in the sections which contain the component specific parameters.
For example the information in the run/GM/restartIN/restart.H 
file or parameters specific to a component can be included.
'}],'name' => 'command','attrib' => {'name' => 'INCLUDE'}},{'type' => 'e','content' => [{'type' => 't','content' => '

#END

The #END command signals the end of the included file or the
end of the PARAM.in file. Lines following the #END command are
ignored. It is not required to use the #END command. The end
of the included file or PARAM.in file is equivalent with an 
#END command in the last line.
'}],'name' => 'command','attrib' => {'name' => 'END'}},{'type' => 'e','content' => [{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'logical','default' => 'T','name' => 'UseStrict'}},{'type' => 't','content' => '
#STRICT
T                       UseStrict

If UseStrict is true, the SWMF does not attempt to correct problems,
but it stops after the warning message. If it is set to false, SWMF
attempts to correct the problems after the warning message is issued.
It is possible to switch back and forth between strict and non-strict mode. 
The default is strict mode.
'}],'name' => 'command','attrib' => {'name' => 'STRICT'}},{'type' => 'e','content' => [{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'string','length' => '100','name' => 'StringDescription'}},{'type' => 't','content' => '

#DESCRIPTION
This is a test run for GM-IE-UA coupling.

The StringDescription string can be used to describe the simulation 
for which the parameter file is written. The #DESCRIPTION command and 
the StringDescription string are saved into the restart file, 
which helps in identifying the restart files.

The default value is "Please describe me!", which is self explanatory.
'}],'name' => 'command','attrib' => {'multiple' => 'T','name' => 'DESCRIPTION'}}],'name' => 'commandgroup','attrib' => {'name' => 'GENERAL COMMANDS'}},{'type' => 'e','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! TIME AND SESSION CONTROL  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

The execution of SWMF is done in consecutive sessions. Each session is
executed with a set of parameters read from the parameter file.
After the session is finished, CON reads and distributes the parameters for
the next session. Parameters from previous sessions are carried over,
so only the changes relative to the previous session need to be given.

'},{'type' => 'e','content' => [{'type' => 't','content' => '

#RUN

The #RUN command does not have any parameters. It signals the end
of the current session, and makes CON execute the session with
the current set of parameters. The parameters for the next session
start after the #RUN command. For the last session there is no
need to use the #RUN command, since the #END command or simply
the end of the PARAM.in file makes CON execute the last session.
'}],'name' => 'command','attrib' => {'name' => 'RUN'}},{'type' => 'e','content' => [{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'logical','default' => 'T','name' => 'DoTimeAccurate'}},{'type' => 't','content' => '

#TIMEACCURATE
F		DoTimeAccurate

If DoTimeAccurate is set to true, the SWMF is solving
a time dependent problem. If DoTimeAccurate is false, a steady-state
solution is sought for. It is possible to use steady-state mode
in the first few sessions to obtain a steady state solution,
and then to switch to time accurate mode in the following sessions.
In time accurate mode the frequency of coupling, saving restart files,
or stopping conditions are taken in simulation time, which is the
time relative to the initial time. In steady state mode the simulation
time is not advanced at all, instead the time step or iteration number
is used to control the frequencies of various actions.

The steady-state mode also allows the components to use various
acceleration techniques. For example the BATSRUS code (in the GM, IH, SC
components) can use local time stepping to accelerate convergence
to steady state. It is also envisioned that in steady state mode
the various components will be allowed to use different number
of iterations per global iteration, thus they can converge to steady state
at the same rate and an optimal global convergence can be achieved.

The default value is the time accurate mode.
'}],'name' => 'command','attrib' => {'name' => 'TIMEACCURATE'}},{'type' => 'e','content' => [{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'integer','default' => '2000','name' => 'iYear'}},{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'integer','default' => '3','max' => '12','name' => 'iMonth','min' => '1'}},{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'integer','default' => '21','max' => '31','name' => 'iDay','min' => '1'}},{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'integer','default' => '0','max' => '23','name' => 'iHour','min' => '0'}},{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'integer','default' => '0','max' => '59','name' => 'iMinute','min' => '0'}},{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'integer','default' => '0','max' => '59','name' => 'iSecond','min' => '0'}},{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'real','default' => '0','max' => '1','name' => 'FracSecond','min' => '0'}},{'type' => 't','content' => '

#STARTTIME
2000			iYear
   3			iMonth
  21			iDay
  10			iHour
  45			iMinute
   0			iSecond
   0.0			FracSecond

The #STARTTIME command sets the initial date and time for the
simulation in Greenwich Mean Time (GMT) or Universal Time (UT).
This time is stored in the restart files. 

The default values are shown above.
This is a date and time when, for the Earth, both the rotational and the magnetic axes
have approximately zero tilt towards the Sun.
'}],'name' => 'command','attrib' => {'if' => '$_IsFirstSession','name' => 'STARTTIME'}},{'type' => 'e','content' => [{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'real','default' => '0.0','name' => 'tSimulation','min' => '0'}},{'type' => 't','content' => '

#TIMESIMULATION
3600.0			tSimulation [sec]

The tSimulation variable contains the simulation time in seconds
relative to the initial time set by the #STARTTIME command.
The #TIMESIMULATION command and tSimulation are saved into the restart file, 
so that the simulation can continue from the same time when the 
restart was saved.  
The default value is tSimulation=0.
'}],'name' => 'command','attrib' => {'if' => '$_IsFirstSession','name' => 'TIMESIMULATION'}},{'type' => 'e','content' => [{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'integer','default' => '0','name' => 'nStep','min' => '0'}},{'type' => 't','content' => '

#NSTEP
100                     nStep

The nStep variable contains the number of time steps since the
beginning of the simulation (including all restarts). 
The #NSTEP command and the nStep variable  are saved into the restart file, 
so that the simulation can continue from the same time step when
the restart was saved.
The default value is nStep=0.
'}],'name' => 'command','attrib' => {'if' => '$_IsFirstSession','name' => 'NSTEP'}},{'type' => 'e','content' => [{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'integer','default' => '-1','name' => 'MaxIter','min' => '-1'}},{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'real','default' => '-1','name' => 'TimeMax','min' => '-1'}},{'type' => 't','content' => '

#STOP
100			MaxIteration
10.0                    tSimulationMax [sec]

The MaxIteration variable contains the
maximum number of iterations {\\it since the beginning of the current run}
(in case of a restart, the time steps done before the restart do not count).
If nIteration reaches this value the session is finished.
The tSimulationMax variable contains the maximum simulation time
relative to the initial time determined by the #STARTTIME command.
If tSimulation reaches this value the session is finished.

Using a negative value for either variables means that the
corresponding condition is  not checked. The default values
are MaxIteration=0 and tSimulationMax = 0.0, so the #STOP command
must be used in every session.
'}],'name' => 'command','attrib' => {'required' => 'T','name' => 'STOP'}},{'type' => 'e','content' => [{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'logical','default' => 'F','name' => 'DoCheckStop'}},{'type' => 'e','content' => [{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'integer','default' => '-1','name' => 'DnCheckStop','min' => '-1'}},{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'real','default' => '-1','name' => 'DtCheckStop','min' => '-1'}}],'name' => 'if','attrib' => {'expr' => '$DoCheckStop'}},{'type' => 't','content' => '

#CHECKSTOP
T			DoCheckStop
-1			DnCheckStop
10.0			DtCheckStop

The DoCheckStop variable controls whether CON should
check the CPU time or the existence of the SWMF.STOP file in the
run directory. If it is set to false, there is no checking.
If it is set to true, the stop conditions are checked at the
frequencies given by the DnCheckStop and DtCheckStop variables.
The DnCheckStop variable determines the frequency in terms of
the time step number nStep, while  DtCheckStop determines the
frequency in terms of the simulation time tSimulation.
Negative values for either variable mean that the corresponding condition
is not checked. For time accurate mode DtCheckStop, for
steady-state mode DnCheckStop is the relevant frequency.

The default value is DoCheckStop=.false., because the checks require
synchronization of the components. The more frequent the checks
are the less efficient the execution.  On the other hand the
less frequent the checks are, the less control the user has to stop
the code at a given time.
'}],'name' => 'command','attrib' => {'name' => 'CHECKSTOP'}},{'type' => 'e','content' => [{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'logical','default' => 'T','name' => 'DoCheckStopFile'}},{'type' => 't','content' => '

#CHECKSTOPFILE
T			DoCheckStopFile

If DoCheckStopFile is true (and DoCheckStop is set
to true in the #CHECKSTOP command) then the code checks if the
SWMF.STOP file exists in the run directory. This file is deleted at
the beginning of the run, so the user must explicitly create the file
with, for example, the "touch SWMF.STOP" UNIX command. 
If the file is found in the run directory,
the execution stops in a graceful manner.
Restart files and plot files are saved as required by the
appropriate parameters.

The default is DoCheckStopFile=.true. (but the default for DoCheckStop
is .false.).
'}],'name' => 'command','attrib' => {'name' => 'CHECKSTOPFILE'}},{'type' => 'e','content' => [{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'real','default' => '-1','name' => 'CpuTimeMax','min' => '-1'}},{'type' => 't','content' => '

#CPUTIMEMAX
3600                    CpuTimeMax [sec]

The CpuTimeMax variable contains the maximum allowed CPU time (wall clock
time) for the execution of the current run. If the CPU time reaches
this time, the execution stops in a graceful manner.
Restart files and plot files are saved as required by the
appropriate parameters.
This command is very useful when the code is submitted to a batch
queue with a limited wall clock time.

The default value is -1.0, which means that the CPU time is not checked.
To do the check the CpuTimeMax variable has to be set to a positive
value and the DoCheckStop variable also must be set to .true. 
in the #CHECKSTOP command.
'}],'name' => 'command','attrib' => {'name' => 'CPUTIMEMAX'}}],'name' => 'commandgroup','attrib' => {'name' => 'TIME AND SESSION CONTROL'}},{'type' => 'e','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!  TESTING AND TIMING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'string','length' => '100','name' => 'StringTest'}},{'type' => 't','content' => '

#TEST
CON_session::do_session		StringTest

A space separated list of subroutine and function names to be tested. 
Only subroutines containing the \'call CON_set_do_test(...)\' statement
can be tested. The first argument is the name of the subroutine,
usually defined as the string parameter \'NameSub\'.
Default is an empty string.
'}],'name' => 'command','attrib' => {'name' => 'TEST'}},{'type' => 'e','content' => [{'type' => 'e','content' => [{'type' => 'e','content' => [],'name' => 'option','attrib' => {'value' => '0','name' => 'errors and warnings only'}},{'type' => 'e','content' => [],'name' => 'option','attrib' => {'default' => 'T','value' => '1','name' => 'normal'}},{'type' => 'e','content' => [],'name' => 'option','attrib' => {'value' => '10','name' => 'calls on test processor'}},{'type' => 'e','content' => [],'name' => 'option','attrib' => {'value' => '100','name' => 'calls on all processors'}}],'name' => 'parameter','attrib' => {'type' => 'integer','input' => 'select','name' => 'iVarTest'}},{'type' => 't','content' => '

#VERBOSE
100			lVerbose

The lVerbose variable sets the verbosity of CON.\\\\
If lVerbose=0 the verbose information is minimized.\\\\
If lVerbose=1 some verbose information is printed in CON_main and CON_io
    from processor zero.\\\\
If lVerbose=10, processor zero will produce a line on the standard output
    with the name of the subroutine and the iteration number 
    for all subroutines which call CON_set_do_test.\\\\
If lVerbose=100, all processors and all subroutines which call CON_set_do_test
    will produce a line on the standard output with the 
    name of the subroutine, the iteration number and the processor number.\\\\
The default value is lVerbose=1.
'}],'name' => 'command','attrib' => {'name' => 'VERBOSE'}},{'type' => 'e','content' => [{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'logical','default' => 'T','name' => 'UseTiming'}},{'type' => 'e','content' => [{'type' => 'e','content' => [{'type' => 'e','content' => [],'name' => 'option','attrib' => {'value' => '-3','name' => 'none'}},{'type' => 'e','content' => [],'name' => 'option','attrib' => {'default' => 'T','value' => '-2','name' => 'final only'}},{'type' => 'e','content' => [],'name' => 'option','attrib' => {'value' => '-1','name' => 'end of sessions'}},{'type' => 'e','content' => [],'name' => 'optioninput','attrib' => {'default' => '100','name' => 'every X steps','min' => '1'}}],'name' => 'parameter','attrib' => {'type' => 'integer','input' => 'select','name' => 'Frequency'}},{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'integer','default' => '-1','name' => 'nDepthTiming','min' => '-1'}},{'type' => 'e','content' => [{'type' => 'e','content' => [],'name' => 'option','attrib' => {'default' => 'T','value' => 'cumu','name' => 'cumulative'}},{'type' => 'e','content' => [],'name' => 'option','attrib' => {'name' => 'list'}},{'type' => 'e','content' => [],'name' => 'option','attrib' => {'name' => 'tree'}}],'name' => 'parameter','attrib' => {'type' => 'string','input' => 'select','name' => 'TypeTimingReport'}}],'name' => 'if','attrib' => {'expr' => '$UseTiming'}},{'type' => 't','content' => '

#TIMING
T                       UseTiming      (rest of parameters read if true)
-2                      DnTiming       (-3 none, -2 final, -1 each session/AMR)
-1                      nDepthTiming   (-1 for arbitrary depth)
tree                    TypeTimingReport   (\'cumu\', \'list\', or \'tree\')

If UseTiming=.true., the execution is timed by the TIMING utility.\\\\
If UseTiming=.false., the execution is not timed.\\\\

\\noindent
The Dntiming parameter determines the frequency of timing reports:\\\\
If DnTiming is positive, a timing report is produced every DnTiming steps.\\\\
If DnTiming is -1, a timing report is shown at the end of each session.\\\\
If DnTiming is -2, a timing report is shown at the end of the whole run.\\\\
If DnTiming is -3, no timing report is shown.\\\\

\\noindent
The nDepthTiming parameters defines the depth of the timing tree.\\\\
A negative value means unlimited depth.\\\\
If nDepthTiming is 1, only the total SWMF execution is timed.\\\\

\\noindent
The TypeTimingReport parameter determines the format of the timing reports:
\\begin{verbatim}
\'cumu\' - cumulative list sorted by timings
\'list\' - list based on caller and sorted by timings
\'tree\' - tree based on calling sequence.
\\end{verbatim}
The default values are shown above.
'}],'name' => 'command','attrib' => {'name' => 'TIMING'}},{'type' => 'e','content' => [{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'integer','default' => '10','name' => 'DnProgressShort','min' => '-1'}},{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'integer','default' => '100','name' => 'DnProgressLong','min' => '-1'}},{'type' => 't','content' => '
#PROGRESS
10			DnProgressShort
100			DnProgressLong

The DnShowProgressShort and DnShowProgressLong variables determine
the frequency of showing short and long progress reports in terms
of the number of time steps nStep. The short progress report
consists of a single line which shows the number of time steps,
simulation time and CPU time. The long progress report also shows
a small timing report on the root processor.
Negative values indicate that no report is requested.

The default values are DnShowProgressShort=10 and DnShowProgressLong=100.
'}],'name' => 'command','attrib' => {'name' => 'PROGRESS'}},{'type' => 'e','content' => [{'type' => 'e','content' => [{'type' => 'e','content' => [],'name' => 'option','attrib' => {'default' => '$nByteReal==4','value' => '4','name' => 'single precision (4)'}},{'type' => 'e','content' => [],'name' => 'option','attrib' => {'default' => '$nByteReal==8','value' => '8','name' => 'double precision (8)'}}],'name' => 'parameter','attrib' => {'type' => 'integer','input' => 'select','name' => 'nByteReal'}},{'type' => 'e','content' => [{'type' => 't','content' => '
		nByteReal in file must agree with _nByteReal.
	'}],'name' => 'rule','attrib' => {'expr' => '$nByteReal==$_nByteReal'}},{'type' => 't','content' => '

#PRECISION
8			nByteReal

The nByteReal variable gives the number of bytes in the default real variable. 
Possible values are 4 and 8. This command serves to check consistency with 
binary data, such as binary restart files. 
The #PRECISION command and the nByteReal variable
are saved into the restart file. 

There is no default value. If the command is not used, the precision of 
the real numbers is not checked.
'}],'name' => 'command','attrib' => {'if' => '$_IsFirstSession','name' => 'PRECISION'}},{'type' => 'e','content' => [{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'real','default' => '1.0','name' => 'CodeVersion','min' => '0'}},{'type' => 't','content' => '

#VERSION
1.0                    Version

The Version variable contains the version number of SWMF. This command
serves to check consistency for restarted runs. The #VERSION command
and the Version variable are saved into the restart file.

There is no default value. If the command is not used, the version
number is not checked.
'}],'name' => 'command','attrib' => {'if' => '$_IsFirstSession','name' => 'VERSION'}}],'name' => 'commandgroup','attrib' => {'name' => 'TESTING AND TIMING'}},{'type' => 'e','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!! COMPONENT CONTROL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

The session model allows concurrent execution of the components.
The components are synchronized only when conditions for the actions of 
coupling, saving restart files, or stopping execution are met. 
This is possible because it is known in advance by each component when 
these actions will occur. In time accurate mode the components are 
required to make a time step which does not exceed the next synchronization
time. In steady state mode there is no limit on the time step, and 
components can be called at different frequencies.
Components which are not used in a session can be switched off.

'},{'type' => 'e','content' => [{'type' => 'e','content' => [{'type' => 'e','content' => [],'name' => 'option','attrib' => {'name' => 'SC'}},{'type' => 'e','content' => [],'name' => 'option','attrib' => {'default' => 'T','name' => 'IH'}},{'type' => 'e','content' => [],'name' => 'option','attrib' => {'name' => 'SP'}},{'type' => 'e','content' => [],'name' => 'option','attrib' => {'name' => 'GM'}},{'type' => 'e','content' => [],'name' => 'option','attrib' => {'name' => 'IM'}},{'type' => 'e','content' => [],'name' => 'option','attrib' => {'name' => 'RB'}},{'type' => 'e','content' => [],'name' => 'option','attrib' => {'name' => 'IE'}},{'type' => 'e','content' => [],'name' => 'option','attrib' => {'name' => 'UA'}}],'name' => 'parameter','attrib' => {'type' => 'string','input' => 'options','name' => 'NameComp'}},{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'logical','default' => 'T','name' => 'UseComp'}},{'type' => 'e','content' => [{'type' => 't','content' => '
		Component must be registered.
	'}],'name' => 'rule','attrib' => {'expr' => '$_Registered{$NameComp} or not $_Components'}},{'type' => 'e','content' => [],'name' => 'set','attrib' => {'type' => 'logical','value' => '$UseComp','name' => '_UsedComp{$NameComp}'}},{'type' => 't','content' => '

#COMPONENT
IE			NameComp
F			UseComp

The NameComp variable contains the two-character component ID, 
while the UseComp variable defines if the component should be 
used in the current session or not. It is possible that in the 
initial sessions some components are not used, and they are turned on later.
Unused components should not be coupled, and they should not
be given any parameters.

The default is that all the components registered in the LAYOUT.in
file are used.
'}],'name' => 'command','attrib' => {'multiple' => 'T','name' => 'COMPONENT'}},{'type' => 'e','content' => [{'type' => 'e','content' => [{'type' => 'e','content' => [],'name' => 'option','attrib' => {'default' => 'T','name' => 'SC'}},{'type' => 'e','content' => [],'name' => 'option','attrib' => {'name' => 'IH'}},{'type' => 'e','content' => [],'name' => 'option','attrib' => {'name' => 'SP'}},{'type' => 'e','content' => [],'name' => 'option','attrib' => {'name' => 'GM'}},{'type' => 'e','content' => [],'name' => 'option','attrib' => {'name' => 'IM'}},{'type' => 'e','content' => [],'name' => 'option','attrib' => {'name' => 'RB'}},{'type' => 'e','content' => [],'name' => 'option','attrib' => {'name' => 'IE'}},{'type' => 'e','content' => [],'name' => 'option','attrib' => {'name' => 'UA'}}],'name' => 'parameter','attrib' => {'type' => 'string','input' => 'options','name' => 'NameComp'}},{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'integer','default' => '1','name' => 'DnRun','min' => '1'}},{'type' => 't','content' => '

#CYCLE
IH			NameComp
10			DnRun

The DnRun variable defines the frequency of calling component NameComp
during a steady state run. In the example IH will be called for 
nStep = 10, 20, 30, ... For time accurate runs this command has no effect.
The default is DnRun = 1 for all the active components.
'}],'name' => 'command','attrib' => {'multiple' => 'T','name' => 'CYCLE'}}],'name' => 'commandgroup','attrib' => {'name' => 'COMPONENT CONTROL'}},{'type' => 'e','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!! COUPLING CONTROL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'integer','default' => '13','max' => '13','name' => 'nCouple','min' => '0'}},{'type' => 'e','content' => [{'type' => 'e','content' => [{'type' => 'e','content' => [],'name' => 'option','attrib' => {'default' => '$iCouple==1','name' => 'SC IH'}},{'type' => 'e','content' => [],'name' => 'option','attrib' => {'default' => '$iCouple==2','name' => 'IH SC'}},{'type' => 'e','content' => [],'name' => 'option','attrib' => {'default' => '$iCouple==3','name' => 'SC SP'}},{'type' => 'e','content' => [],'name' => 'option','attrib' => {'default' => '$iCouple==4','name' => 'IH SP'}},{'type' => 'e','content' => [],'name' => 'option','attrib' => {'default' => '$iCouple==5','name' => 'IH GM'}},{'type' => 'e','content' => [],'name' => 'option','attrib' => {'default' => '$iCouple==6','name' => 'GM IE'}},{'type' => 'e','content' => [],'name' => 'option','attrib' => {'default' => '$iCouple==7','name' => 'GM IM'}},{'type' => 'e','content' => [],'name' => 'option','attrib' => {'default' => '$iCouple==8','name' => 'GM RB'}},{'type' => 'e','content' => [],'name' => 'option','attrib' => {'default' => '$iCouple==9','name' => 'UA IE'}},{'type' => 'e','content' => [],'name' => 'option','attrib' => {'default' => '$iCouple==10','name' => 'IE IM'}},{'type' => 'e','content' => [],'name' => 'option','attrib' => {'default' => '$iCouple==11','name' => 'IM GM'}},{'type' => 'e','content' => [],'name' => 'option','attrib' => {'default' => '$iCouple==12','name' => 'IE UA'}},{'type' => 'e','content' => [],'name' => 'option','attrib' => {'default' => '$iCouple==13','name' => 'IE GM'}}],'name' => 'parameter','attrib' => {'type' => 'string','input' => 'select','name' => 'NameSourceTarget'}}],'name' => 'for','attrib' => {'from' => '1','name' => 'iCouple','to' => '$nCouple'}},{'type' => 't','content' => '

#COUPLEORDER
13                nCouple
SC IH             NameSourceTarget
IH SC             NameSourceTarget
SC SP             NameSourceTarget
IH SP             NameSourceTarget
IH GM             NameSourceTarget
GM IE             NameSourceTarget
GM IM             NameSourceTarget
GM RB             NameSourceTarget
UA IE             NameSourceTarget
IE IM             NameSourceTarget
IM GM             NameSourceTarget
IE UA             NameSourceTarget
IE GM             NameSourceTarget

The nCouple variable determines the maximum number of couplings
among the components. The NameSourceTarget string contains the
two-character ID for the source and target components separated by
spaces. The order of the couplings is most important for the
couplings done at the beginning of the session. Some components need
to get information from the others before they can initialize properly.

The default coupling order is shown above. It is based on the propagation 
of information from component to component.

NOTE: The order of the \'SC SP\' and \'IH SP\' couplings should not be reversed!
'}],'name' => 'command','attrib' => {'name' => 'COUPLEORDER'}},{'type' => 'e','content' => [{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'string','length' => '2','name' => 'NameSource'}},{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'string','length' => '2','name' => 'NameTarget'}},{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'integer','default' => '-1','name' => 'DnCouple','min' => '-1'}},{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'real','default' => '-1.0','name' => 'DtCouple','min' => '-1'}},{'type' => 'e','content' => [{'type' => 't','content' => '
		Source component must be registered and ON.
	'}],'name' => 'rule','attrib' => {'expr' => '$_UsedComp{$NameSource} or not $_Components'}},{'type' => 'e','content' => [{'type' => 't','content' => '
		Target component must be registered and ON.
	'}],'name' => 'rule','attrib' => {'expr' => '$_UsedComp{$NameTarget} or not $_Components'}},{'type' => 'e','content' => [{'type' => 't','content' => '
		Source and target components must be different.
	'}],'name' => 'rule','attrib' => {'expr' => '$NameSource ne $NameTarget'}},{'type' => 't','content' => '

#COUPLE1
IE			NameSource
IM			NameTarget
-1			DnCouple
10.0			DtCouple

The NameSource and NameTarget variables contain the two-character 
component ID-s for the source and target components. 
The DnCouple variable determines the frequency
of couplings in terms of the number of time steps nStep for steady state
runs, while DtCouple defines the frequency in terms of the simulation
time tSimulation in seconds for time-accurate runs. Setting both
frequencies to a negative value means that there is no coupling.

The default is no coupling between the components.
'}],'name' => 'command','attrib' => {'multiple' => 'T','name' => 'COUPLE1'}},{'type' => 'e','content' => [{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'string','length' => '2','name' => 'NameComp1'}},{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'string','length' => '2','name' => 'NameComp2'}},{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'integer','default' => '-1','name' => 'DnCouple','min' => '-1'}},{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'real','default' => '-1.0','name' => 'DtCouple','min' => '-1'}},{'type' => 'e','content' => [{'type' => 't','content' => '
		Component 1 must be registered and ON.
	'}],'name' => 'rule','attrib' => {'expr' => '$_UsedComp{$NameComp1} or not $_Components'}},{'type' => 'e','content' => [{'type' => 't','content' => '
		Component 2 must be registered and ON.
	'}],'name' => 'rule','attrib' => {'expr' => '$_UsedComp{$NameComp2} or not $_Components'}},{'type' => 'e','content' => [{'type' => 't','content' => '
		Components 1 and 2 must be different.
	'}],'name' => 'rule','attrib' => {'expr' => '$NameComp1 ne $NameComp2'}},{'type' => 't','content' => '

#COUPLE2
GM			NameComp1
IE			NameComp2
-1			DnCouple
10.0			DtCouple

The NameComp1 and NameComp2
variables contain the two-character component ID-s for the two
components which are coupled both ways.
The DnCouple variable determines the frequency
of couplings in terms of the number of time steps nStep for steady state runs, 
while DtCouple defines the frequency in terms of the simulation
time tSimulation in seconds for time-accurate runs. Setting both
frequencies to a negative value means that there is no coupling.

The default is no coupling between the components.
'}],'name' => 'command','attrib' => {'multiple' => 'T','name' => 'COUPLE2'}},{'type' => 'e','content' => [{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'string','length' => '2','name' => 'NameSource'}},{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'string','length' => '2','name' => 'NameTarget'}},{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'integer','default' => '-1','name' => 'DnCouple','min' => '-1'}},{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'real','default' => '-1.0','name' => 'DtCouple','min' => '-1'}},{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'integer','default' => '0.0','max' => '$DnCouple','name' => 'nNext12','min' => '-1'}},{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'real','default' => '0.0','max' => '$DtCouple','name' => 'tNext12','min' => '-1'}},{'type' => 'e','content' => [{'type' => 't','content' => '
		Source component must be registered and ON.
	'}],'name' => 'rule','attrib' => {'expr' => '$_UsedComp{$NameSource} or not $_Components'}},{'type' => 'e','content' => [{'type' => 't','content' => '
		Target component must be registered and ON.
	'}],'name' => 'rule','attrib' => {'expr' => '$_UsedComp{$NameTarget} or not $_Components'}},{'type' => 'e','content' => [{'type' => 't','content' => '
		Source and target components must be different.
	'}],'name' => 'rule','attrib' => {'expr' => '$NameSource ne $NameTarget'}},{'type' => 't','content' => '

#COUPLE1SHIFT
IH			NameSource
GM			NameTarget
-1			DnCouple
10.0			DtCouple
-1			nNext12
3.0			tNext12

The NameSource and NameTarget
variables contain the two-character component ID-s for the source
and target components. The DnCouple variable determines the frequency
of couplings in terms of the number of time steps nStep for steady state
runs, while DtCouple defines the frequency in terms of the simulation
time tSimulation in seconds for time-accurate runs.

For steady-state simulations the nNext12 variable determines in which time
step the first coupling occurs after the initial coupling, namely
when mod(nStep,DnCouple) equals nNext12. For time accurate
simulations the tNext12 variable determines at what simulation
time the first coupling occurs after the initial coupling, namely
when mod(tSimulation,DtCouple) equals tNext12.

The above example will couple IH to GM at simulation times 3, 13, 23, etc.

The default is no shifting.
'}],'name' => 'command','attrib' => {'multiple' => 'T','name' => 'COUPLE1SHIFT'}},{'type' => 'e','content' => [{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'string','length' => '2','name' => 'NameComp1'}},{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'string','length' => '2','name' => 'NameComp2'}},{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'integer','default' => '-1','name' => 'DnCouple','min' => '-1'}},{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'real','default' => '-1.0','name' => 'DtCouple','min' => '-1'}},{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'integer','default' => '0.0','max' => '$DnCouple','name' => 'nNext12','min' => '-1'}},{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'real','default' => '0.0','max' => '$DtCouple','name' => 'tNext12','min' => '-1'}},{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'integer','default' => '0.0','max' => '$DnCouple','name' => 'nNext21','min' => '-1'}},{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'real','default' => '0.0','max' => '$DtCouple','name' => 'tNext21','min' => '-1'}},{'type' => 'e','content' => [{'type' => 't','content' => '
		Component 1 must be registered and ON.
	'}],'name' => 'rule','attrib' => {'expr' => '$_UsedComp{$NameComp1} or not $_Components'}},{'type' => 'e','content' => [{'type' => 't','content' => '
		Component 2 must be registered and ON.
	'}],'name' => 'rule','attrib' => {'expr' => '$_UsedComp{$NameComp2} or not $_Components'}},{'type' => 'e','content' => [{'type' => 't','content' => '
		Components 1 and 2 must be different.
	'}],'name' => 'rule','attrib' => {'expr' => '$NameComp1 ne $NameComp2'}},{'type' => 't','content' => '

#COUPLE2SHIFT
GM                 NameComp1
IE                 NameComp2
-1                 DnCouple
10.0               DtCouple
-1                  nNext12
3.0                tNext12
-1                 nNext21
6.0                tNext21

The NameComp1 and NameComp2
variables contain the two-character component ID-s for the two
components which are coupled both ways.
The DnCouple variable determines the frequency
of couplings in terms of the number of time steps nStep for steady state
runs, while DtCouple defines the frequency in terms of the simulation
time tSimulation in seconds for time-accurate runs.

For steady-state simulations the nNext12 variable determines in which time
step the first coupling occurs from NameComp1 to NameComp2
after the initial coupling, namely
when mod(nStep,DnCouple) equals nNext12. For time accurate
simulations the tNext12 variable determines at what simulation
time the first coupling occurs from NameComp1 to NameComp2
after the initial coupling, namely
when mod(tSimulation,DtCouple) equatls tNext12.

The first coupling step and time for the NameComp2 to NameComp1
coupling is determined by the nNext21 and tNext21 variables
in a similar fashion. This command allows to shift the couplings
relative to each other. 

The above example will couple GM to IE at simulation times 
3, 13, 23, etc, while IE will be coupled to GM at simulation times 
6, 16, 26 etc. 
This way IE can solve the potential problem while GM advances by 3 seconds.
That can improve the parallelization and efficiency.

The default is no shifting.
'}],'name' => 'command','attrib' => {'multiple' => 'T','name' => 'COUPLE2SHIFT'}},{'type' => 'e','content' => [{'type' => 'e','content' => [{'type' => 'e','content' => [],'name' => 'option','attrib' => {'default' => 'T','name' => 'SC'}},{'type' => 'e','content' => [],'name' => 'option','attrib' => {'name' => 'IH'}},{'type' => 'e','content' => [],'name' => 'option','attrib' => {'name' => 'SP'}},{'type' => 'e','content' => [],'name' => 'option','attrib' => {'name' => 'GM'}},{'type' => 'e','content' => [],'name' => 'option','attrib' => {'name' => 'IM'}},{'type' => 'e','content' => [],'name' => 'option','attrib' => {'name' => 'RB'}},{'type' => 'e','content' => [],'name' => 'option','attrib' => {'name' => 'IE'}},{'type' => 'e','content' => [],'name' => 'option','attrib' => {'name' => 'UA'}}],'name' => 'parameter','attrib' => {'type' => 'string','input' => 'options','name' => 'NameComp'}},{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'logical','default' => 'T','name' => 'DoCoupleOnTime'}},{'type' => 't','content' => '

#COUPLETIME
GM			NameComp
F			DoCoupleOnTime

The NameComp variable contains the two-character component ID, 
while the DoCoupleOnTime parameter defines if the time step of the
component should be limited such that it does not exceed the next 
coupling time. If DoCoupleOnTime is true, the component will limit the
time step for couplings, if it is false, the time step is only limited
by the final time, the time of saving restarts and the time of checking
stop conditions. 

The default is that all components limit their time steps to match
the coupling time.
'}],'name' => 'command','attrib' => {'multiple' => 'T','name' => 'COUPLETIME'}}],'name' => 'commandgroup','attrib' => {'name' => 'COUPLING CONTROL'}},{'type' => 'e','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!! RESTART CONTROL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CON needs to coordinate the saving of restart information for itself
and all the components. It is important that all components save
the necessary information at the same simulation time.

'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'logical','default' => 'T','name' => 'DoSaveRestart'}},{'type' => 'e','content' => [{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'integer','default' => '-1','name' => 'DnSaveRestart','min' => '-1'}},{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'real','default' => '-1','name' => 'DtSaveRestart','min' => '-1'}}],'name' => 'if','attrib' => {'expr' => '$SaveRestart'}},{'type' => 't','content' => '

#SAVERESTART
T			DoSaveRestart (Rest of parameters read if true)
-1			DnSaveRestart
-1.			DtSaveRestart

The DoSaveRestart variable determines whether restart information 
should be saved or not. The rest of the parameters are read only 
if DoSaveRestart is true.
For steady state runs the frequency of saving restart information is given by
the DnSaveRestart variable in terms of the number of time steps nStep,
while for time accurate run, the DtSaveRestart variable determines
the frequency in terms of the simulation time tSimulation in seconds.
Negative frequencies mean that no restart file is saved during the run,
but they do not exclude the saving of the restart information at
the very end. 

The default is given above, which means that restart information
is saved at the end of the run.
'}],'name' => 'command','attrib' => {'name' => 'SAVERESTART'}},{'type' => 'e','content' => [{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'string','default' => 'RESTART.out','length' => '100','name' => 'NameRestartFile'}},{'type' => 't','content' => '

#RESTARTFILE
RESTART_test.in

The NameRestartFile variable contains the file name
for the restart file for CON itself. This file contains information
such as initial time, simulation time and time step, version number,
description of the simulation, name of the planet, etc. The file
is in the same format as PARAM.in, so it can be simply included
with the #INCLUDE command. To avoid unpleasant surprises,
do not include the file which is being written.

The default value for NameRestartFile is "RESTART.out".
'}],'name' => 'command','attrib' => {'name' => 'RESTARTFILE'}}],'name' => 'commandgroup','attrib' => {'name' => 'RESTART CONTROL'}},{'type' => 'e','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!! OUTPUT CONTROL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'logical','default' => 'F','name' => 'DoEcho'}},{'type' => 't','content' => '

#ECHO
T			DoEcho

If the DoEcho variable is true, the input parameters are echoed back.
The echoing either goes to the standard output or into the log files
depending on the UseStdout variable, which can be set in the
#STDOUT command.

The default value for DoEcho is .false., but it is a good idea to
set it to true at the beginning of the PARAM.in file.
'}],'name' => 'command','attrib' => {'name' => 'ECHO'}},{'type' => 'e','content' => [{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'logical','default' => 'T','name' => 'UseStdout'}},{'type' => 't','content' => '

#STDOUT
F			UseStdout

If the UseStdout variable is true, the echoed input parameters and other
verbose information produced by the components are written to the
standard output. To distinguish between the output produced by the
various components, a prefix string is written at the beginning
of the lines. Usually the prefix string contains the component ID, the
processor number if the line is written by a processor which
is not the root processor of the component, and a colon 
(for example "GM0001:").
Even with the prefix, it may be difficult to collect the output
from the various components and processors. The order of the output
depends on how the MPI library buffers that, which is platform
dependent.

If the UseStdout variable is false, the echoed input parameters and other
verbose information produced by the components are written
into separate files in the STDOUT directory (the name of
this directory can be changed by the #STDOUTDIR command).
The files are named similarly to the prefix string:
the component ID is followed by the global processor number
and a ".log" extension is added. For example the root processor
of GM may write into "STDOUT/GM0014.log" if GM\'s root processor
has global rank 14.

Note that warnings and error messages should always be written to
the standard output, and they should always have a prefix which
identifies the component issuing the warning or error.
CON itself always writes to the standard output and it does not
use a string prefix.

The default value for UseStdout is true.
'}],'name' => 'command','attrib' => {'name' => 'STDOUT'}},{'type' => 'e','content' => [{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'string','length' => '100','name' => 'NameStdoutDir'}},{'type' => 't','content' => '

#STDOUTDIR
STDOUT/Test		NameStdoutDir

The NameStdoutDir variable contains the name of the directory where 
the log files with the redirected standard output of the components 
are written if UseStdout is set to .false. in the #STDOUT command.
The directory must exist before the run is started.

The default value of NameStdoutDir is "STDOUT".
'}],'name' => 'command','attrib' => {'name' => 'STDOUTDIR'}}],'name' => 'commandgroup','attrib' => {'name' => 'OUTPUT CONTROL'}},{'type' => 'e','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! SOLAR COORDINATE COMMANDS !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

We allow an offset for the HGR and HGI/HGC systems so that the
interesting features are aligned with the primary axis.
One common option is to have the planet in the -X,Z plane.
Another option would be to move an active region into an appropriate plane.

'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'real','default' => '0','max' => '360','name' => 'dLongitudeHgr','min' => '-360'}},{'type' => 't','content' => '
#ROTATEHGR
145.6			dLongitudeHgr [deg]

Rotate the HGR system by dLongitudeHgr degrees around the Z axis.
A negative value is interpreted as an offset angle which moves the
planet into the -X, Z plane (so roughly towards the -X axis).
Default value is 0, i.e. the true HGR system is used.
'}],'name' => 'command','attrib' => {'if' => '$_IsFirstSession','name' => 'ROTATEHGR'}},{'type' => 'e','content' => [{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'real','default' => '0','max' => '360','name' => 'dLongitudeHgi','min' => '-360'}},{'type' => 't','content' => '
#ROTATEHGI
-1.0			dLongitudeHgi [deg]

Rotate the HGI and the related rotating HGC systems by 
dLongitudeHgr degrees around the Z axis.
A negative value is interpreted as an offset angle which moves the
planet into the -X, Z plane (so roughly towards the -X axis).
Default value is 0, i.e. the true HGI system is used.
'}],'name' => 'command','attrib' => {'if' => '$_IsFirstSession','name' => 'ROTATEHGI'}}],'name' => 'commandgroup','attrib' => {'name' => 'SOLAR COORDINATE COMMANDS'}},{'type' => 'e','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! PLANET COMMANDS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Several components share information about the planet for which the
simulation is done. It is important that the various components
use compatible information about the planet. It is also useful
for the couplers that they can globally access this information,
such as radius and orientation of the planet, or its magnetic
field. The SWMF is designed to work for an arbitrary planet.
It also allows to change some parameters of the planet relative
to the real values.

By default the SWMF works with Earth and its real parameters.
Another planet can be selected with the #PLANET command.
The real planet parameters can be modified and simplified
with the other planet commands listed in this subsection.
These modifier commands cannot proceed the #PLANET command!

'},{'type' => 'e','content' => [{'type' => 'e','content' => [{'type' => 'e','content' => [],'name' => 'option','attrib' => {'default' => 'T','value' => 'EARTH/Earth/earth','name' => 'Earth'}},{'type' => 'e','content' => [],'name' => 'option','attrib' => {'value' => 'SATURN/Saturn/saturn','name' => 'Saturn'}},{'type' => 'e','content' => [],'name' => 'option','attrib' => {'name' => 'New'}}],'name' => 'parameter','attrib' => {'type' => 'string','input' => 'select','name' => 'NamePlanet'}},{'type' => 'e','content' => [{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'real','name' => 'RadiusPlanet','min' => '0'}},{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'real','name' => 'MassPlanet','min' => '0'}},{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'real','name' => 'OmegaPlanet','min' => '0'}},{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'real','name' => 'TiltRotation','min' => '0'}},{'type' => 'e','content' => [{'type' => 'e','content' => [],'name' => 'option','attrib' => {'name' => 'NONE'}},{'type' => 'e','content' => [],'name' => 'option','attrib' => {'default' => 'T','name' => 'DIPOLE'}}],'name' => 'parameter','attrib' => {'type' => 'string','input' => 'select','name' => 'TypeBField'}}],'name' => 'if','attrib' => {'expr' => '$NamePlanet eq \'New\''}},{'type' => 'e','content' => [{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'real','max' => '180','name' => 'MagAxisThetaGeo','min' => '0'}},{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'real','max' => '360','name' => 'MagAxisPhiGeo','min' => '0'}},{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'real','name' => 'DipoleStrength'}}],'name' => 'if','attrib' => {'expr' => '$TyepBField eq \'DIPOLE\''}},{'type' => 'e','content' => [{'type' => 't','content' => '
		PLANET should precede $PlanetCommand
	'}],'name' => 'rule','attrib' => {'expr' => 'not $PlanetCommand'}},{'type' => 't','content' => '

#PLANET
New			NamePlanet (rest of parameters read for unknown planet)
6300000.0		RadiusPlanet [m]
5.976E+24		MassPlanet   [kg]
0.000000199		OmegaPlanet  [radian/s]
23.5			TiltRotation [degree]
DIPOLE			TypeBField
11.0			MagAxisThetaGeo [degree]
289.1			MagAxisPhiGeo   [degree]
-31100.0E-9		DipoleStrength  [T]

The NamePlanet parameter contains the name of the planet
with arbitrary capitalization. In case the name of the planet
is not recognized, the following variables are read:
RadiusPlanet is the radius of the planet,
MassPlanet is the mass of the planet, 
OmegaPlanet is the angular speed relative to an inertial frame,
TiltRotation is the tilt of the rotation axis relative to ecliptic North,
TypeBField, which can be "NONE" or "DIPOLE". 
TypeBField="NONE" means that the planet does not have magnetic field. 
It TypeBField is set to "DIPOLE" than the following variables are read:
MagAxisThetaGeo and MagAxisPhiGeo are the colatitude and longitude
of the north magnetic pole in corotating planetocentric coordinates.
Finally DipoleStrength is the equatorial strength of the magnetic dipole
field. The units are indicated in the above example, which shows the
Earth values approximately.

The default value is NamePlanet="Earth", which is currently
the only recognized planet.
'}],'name' => 'command','attrib' => {'if' => '$_IsFirstSession','name' => 'PLANET'}},{'type' => 'e','content' => [{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'logical','default' => 'T','name' => 'IsRotAxisPrimary'}},{'type' => 'e','content' => [{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'real','max' => '180','name' => 'RotAxisTheta','min' => '0'}},{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'real','max' => '360','name' => 'RotAxisPhi','min' => '0'}}],'name' => 'if','attrib' => {'expr' => '$IsRotAxisPrimary'}},{'type' => 'e','content' => [],'name' => 'set','attrib' => {'type' => 'string','value' => 'ROTATIONAXIS','name' => 'PlanetCommand'}},{'type' => 't','content' => '

#ROTATIONAXIS
T			IsRotAxisPrimary (rest of parameters read if true)
23.5			RotAxisTheta
198.3			RotAxisPhi

If the IsRotAxisPrimary variable is false, the rotational axis
is aligned with the magnetic axis. If it is true, the other two variables
are read, which give the position of the rotational axis at the
initial time in the GSE coordinate system. Both angles are read in degrees
and stored internally in radians.

The default is to use the true rotational axis determined by the
date and time given by #STARTTIME.
'}],'name' => 'command','attrib' => {'if' => '$_IsFirstSession','name' => 'ROTATIONAXIS'}},{'type' => 'e','content' => [{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'logical','default' => 'T','name' => 'UseRotation'}},{'type' => 'e','content' => [{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'real','name' => 'RotationPeriod'}}],'name' => 'if','attrib' => {'expr' => '$UseRotation'}},{'type' => 'e','content' => [],'name' => 'set','attrib' => {'type' => 'string','value' => 'ROTATION','name' => 'PlanetCommand'}},{'type' => 't','content' => '

#ROTATION
T			UseRotation
24.06575		RotationPeriod [hour] (read if UseRotation is true)

If UseRotation is false, the planet is assumed to stand still, 
and the OmegaPlanet variable is set to zero. 
If UseRotation is true, the RotationPeriod variable is read in hours, 
and it is converted to the angular speed OmegaPlanet given in radians/second.
Note that OmegaPlanet is relative to an inertial coordinate system,
so the RotationPeriod is not 24 hours for the Earth, but the
length of the astronomical day.

The default is to use rotation with the real rotation period of the planet.
'}],'name' => 'command','attrib' => {'if' => '$_IsFirstSession','name' => 'ROTATION'}},{'type' => 'e','content' => [{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'logical','default' => 'T','name' => 'IsMagAxisPrimary'}},{'type' => 'e','content' => [{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'real','max' => '180','name' => 'MagAxisTheta','min' => '0'}},{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'real','max' => '360','name' => 'MagAxisPhi','min' => '0'}}],'name' => 'if','attrib' => {'expr' => '$IsMagAxisPrimary'}},{'type' => 'e','content' => [],'name' => 'set','attrib' => {'type' => 'string','value' => 'MAGNETICAXIS','name' => 'PlanetCommand'}},{'type' => 't','content' => '

#MAGNETICAXIS
T			IsMagAxisPrimary (rest of parameters read if true)
34.5			MagAxisTheta [degree]
0.0			MagAxisPhi   [degree]

If the IsMagAxisPrimary variable is false, the magnetic axis
is aligned with the rotational axis. If it is true, the other two variables
are read, which give the position of the magnetic axis at the
initial time in the GSE coordinate system. Both angles are read in degrees
and stored internally in radians.

The default is to use the true magnetic axis determined by the
date and time given by #STARTTIME.
'}],'name' => 'command','attrib' => {'if' => '$_IsFirstSession','name' => 'MAGNETICAXIS'}},{'type' => 'e','content' => [{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'real','name' => 'DipoleStrength'}},{'type' => 'e','content' => [],'name' => 'set','attrib' => {'type' => 'string','value' => 'DIPOLE','name' => 'PlanetCommand'}},{'type' => 't','content' => '

#DIPOLE
-3.11e-4		DipoleStrength [Tesla]

The DipoleStrength variable contains the
magnetic equatorial strength of the dipole magnetic field in Tesla.

The default value is the real dipole strength for the planet.
For the Earth the default is taken to be -31100 nT.
The sign is taken to be negative so that the magnetic axis can
point northward as usual.
'}],'name' => 'command','attrib' => {'if' => '$_IsFirstSession','name' => 'DIPOLE'}},{'type' => 'e','content' => [{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'real','default' => '0.0001','name' => 'DtUpdateB0','min' => '-1'}},{'type' => 't','content' => '

The DtUpdateB0 variable determines how often the position of
the magnetic axis is recalculated. A negative value indicates that
the motion of the magnetic axis during the course of the simulation
is neglected. This is an optimization parameter, since recalculating
the values which depend on the orientation of the magnetic
field can be costly. Since the magnetic field moves relatively
slowly as the planet rotates around, it may not be necessary
to continuously update the magnetic field orientation.

The default value is 0.0001, which means that the magnetic axis
is continuously followed.
'}],'name' => 'command','attrib' => {'name' => 'UPDATEB0'}},{'type' => 'e','content' => [{'type' => 'e','content' => [],'name' => 'set','attrib' => {'type' => 'string','value' => 'IDEALAXES','name' => 'PlanetCommand'}},{'type' => 't','content' => '

#IDEALAXES

The #IDEALAXES command has no parameters. It sets both the rotational
and magnetic axes parallel with the ecliptic North direction. In fact
it is identical with

#ROTATIONAXIS
T               IsRotAxisPrimary
0.0             RotAxisTheta
0.0             RotAxisPhi

#MAGNETICAXIS
F               IsMagAxisPrimary

but much shorter.
'}],'name' => 'command','attrib' => {'if' => '$_IsFirstSession','name' => 'IDEALAXES'}}],'name' => 'commandgroup','attrib' => {'name' => 'PLANET COMMANDS'}},{'type' => 'e','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!! STUB COMPONENTS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
If SWMF is compiled with the interface in srcCON/Stubs,
the stub components recognize only one command #TIMESTEP.

'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'real','name' => 'DtRun','min' => '0'}},{'type' => 'e','content' => [],'name' => 'parameter','attrib' => {'type' => 'real','name' => 'DtCpu','min' => '0'}},{'type' => 't','content' => '

#TIMESTEP
0.01       DtRun (the typical time step of the component)
0.12       DtCpu (the CPU time needed for 1 time step)

The DtRun variable defines the typical time step of
the component in terms of simulation time. The DtCpu variable
determines the CPU time needed to execute one time step for the component.
Both variables are given in seconds.

Of course it is not necessary to put in the actual CPU times.
One can take the same fraction for all components to accelerate
the run.
'}],'name' => 'command','attrib' => {'name' => 'TIMESTEP'}}],'name' => 'commandgroup','attrib' => {'name' => 'STUB COMPONENTS'}}],'name' => 'commandList','attrib' => {'name' => 'Control Module'}}];