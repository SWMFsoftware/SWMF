#^CFG FILE _FALSE_
$tree = [{'attrib' => {'name' => 'BATSRUS: GM, SC and IH Components'},'content' => [{'content' => '

List of MH (GM, IH and SC) commands used in the PARAM.in file


','type' => 't'},{'attrib' => {'value' => '$_GridSize[0]','type' => 'integer','name' => 'nI'},'content' => [],'type' => 'e','name' => 'set'},{'attrib' => {'value' => '$_GridSize[1]','type' => 'integer','name' => 'nJ'},'content' => [],'type' => 'e','name' => 'set'},{'attrib' => {'value' => '$_GridSize[2]','type' => 'integer','name' => 'nK'},'content' => [],'type' => 'e','name' => 'set'},{'attrib' => {'value' => '$_GridSize[3]','type' => 'integer','name' => 'MaxBlock'},'content' => [],'type' => 'e','name' => 'set'},{'attrib' => {'value' => '$_GridSize[4]','type' => 'integer','name' => 'MaxImplBlock'},'content' => [],'type' => 'e','name' => 'set'},{'attrib' => {'value' => '$_nProc and $MaxBlock and $_nProc*$MaxBlock','type' => 'integer','name' => 'MaxBlockALL'},'content' => [],'type' => 'e','name' => 'set'},{'attrib' => {'value' => '$_NameComp/restartOUT','type' => 'string','name' => 'NameRestartOutDir'},'content' => [],'type' => 'e','name' => 'set'},{'attrib' => {'value' => '$_NameComp/IO2','type' => 'string','name' => 'NamePlotDir'},'content' => [],'type' => 'e','name' => 'set'},{'attrib' => {'name' => 'STAND ALONE MODE'},'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! STAND ALONE PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

','type' => 't'},{'attrib' => {'name' => 'COMPONENT'},'content' => [{'attrib' => {'input' => 'select','type' => 'string','name' => 'NameComp'},'content' => [{'attrib' => {'if' => '$_NameComp eq \'SC\'','default' => 'T','name' => 'SC'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'if' => '$_NameComp eq \'IH\'','default' => 'T','name' => 'IH'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'if' => '$_NameComp eq \'GM\'','default' => 'T','name' => 'GM'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'parameter'},{'content' => '

#COMPONENT
GM			NameComp

This command is only used in the stand alone mode.

The NameComp variable contains the two-character component ID
for the component which BATSRUS is representing.
If NameComp does not agree with the value of the NameThisComp
variable, BATSRUS stops with an error message.
This command is saved into the restart header file for consistency check.

There is no default value: if the command is not given, the component ID is not checked.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'if' => '$_IsStandAlone','name' => 'DESCRIPTION'},'content' => [{'attrib' => {'length' => '100','type' => 'string','name' => 'StringDescription'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

#DESCRIPTION
This is a test run for Jupiter with no rotation.

This command is only used in the stand alone mode.

The StringDescription string can be used to describe the simulation
for which the parameter file is written. The #DESCRIPTION command and
the StringDescription string are saved into the restart file,
which helps in identifying the restart files.

The default value is ``Please describe me!", which is self explanatory.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'if' => '$_IsStandAlone','name' => 'ECHO'},'content' => [{'attrib' => {'default' => 'F','type' => 'logical','name' => 'DoEcho'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

#ECHO
T                       DoEcho

This command is only used in the stand alone mode.

If the DoEcho variable is true, the input parameters are echoed back.
The default value for DoEcho is .false., but it is a good idea to
set it to true at the beginning of the PARAM.in file.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'if' => '$_IsStandAlone','name' => 'PROGRESS'},'content' => [{'attrib' => {'min' => '-1','default' => '10','type' => 'integer','name' => 'DnProgressShort'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '-1','default' => '100','type' => 'integer','name' => 'DnProgressLong'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '
#PROGRESS
10			DnProgressShort
100			DnProgressLong

The frequency of short and long progress reports for BATSRUS in
stand alone mode. These are the defaults. Set -1-s for no progress reports.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'if' => '$_IsStandAlone','name' => 'TIMEACCURATE'},'content' => [{'attrib' => {'default' => 'T','type' => 'logical','name' => 'DoTimeAccurate'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

#TIMEACCURATE
F               DoTimeAccurate

This command is only used in stand alone mode.

If DoTimeAccurate is set to true, BATSRUS solves
a time dependent problem. If DoTimeAccurate is false, a steady-state
solution is sought for. It is possible to use steady-state mode
in the first few sessions to obtain a steady state solution,
and then to switch to time accurate mode in the following sessions.
In time accurate mode saving plot files, log files and restart files,
or stopping conditions are taken in simulation time, which is the
time relative to the initial time. In steady state mode the simulation
time is not advanced at all, instead the time step or iteration number
is used to control the frequencies of various actions.

The steady-state mode allows BATSRUS to use local time stepping
to accelerate the convergence towards steady state.

The default value depends on how the stand alone code was installed.
See the description of the NEWPARAM command.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'multiple' => 'T','if' => '$_IsStandAlone','name' => 'BEGIN_COMP'},'content' => [{'content' => '

This command is allowed in stand alone mode only for the sake of the 
test suite, which contains these commands when the framework is tested.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'multiple' => 'T','if' => '$_IsStandAlone','name' => 'END_COMP'},'content' => [{'content' => '

This command is allowed in stand alone mode only for the sake of the 
test suite, which contains these commands when the framework is tested.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'if' => '$_IsStandAlone','name' => 'RUN'},'content' => [{'content' => '

#RUN

This command is only used in stand alone mode.

The #RUN command does not have any parameters. It signals the end
of the current session, and makes BATSRUS execute the session with
the current set of parameters. The parameters for the next session
start after the #RUN command. For the last session there is no
need to use the #RUN command, since the #END command or simply
the end of the PARAM.in file makes BATSRUS execute the last session.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'END'},'content' => [{'content' => '

#END

The #END command signals the end of the included file or the
end of the PARAM.in file. Lines following the #END command are
ignored. It is not required to use the #END command. The end
of the included file or PARAM.in file is equivalent with an 
#END command in the last line.
','type' => 't'}],'type' => 'e','name' => 'command'}],'type' => 'e','name' => 'commandgroup'},{'attrib' => {'name' => 'PLANET COMMANDS'},'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! PLANET COMMANDS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

The planet commands can only be used in stand alone mode.
The commands allow to work with an arbitrary planet.
It is also possible to change some parameters of the planet relative
to the real values.

By default Earth is assumed with its real parameters.
Another planet can be selected with the #PLANET command.
The real planet parameters can be modified and simplified
with the other planet commands listed in this subsection.
These modified commands cannot precede the #PLANET command!

','type' => 't'},{'attrib' => {'if' => '$_IsFirstSession and $_IsStandAlone','name' => 'PLANET'},'content' => [{'attrib' => {'input' => 'select','type' => 'string','name' => 'NamePlanet'},'content' => [{'attrib' => {'value' => 'EARTH/Earth/earth','default' => 'T','name' => 'Earth'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'SATURN/Saturn/saturn','name' => 'Saturn'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'New'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'parameter'},{'attrib' => {'expr' => '$NamePlanet eq \'New\''},'content' => [{'attrib' => {'min' => '0','type' => 'real','name' => 'RadiusPlanet'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','type' => 'real','name' => 'MassPlanet'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','type' => 'real','name' => 'OmegaPlanet'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','type' => 'real','name' => 'TiltRotation'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'input' => 'select','type' => 'string','name' => 'TypeBField'},'content' => [{'attrib' => {'name' => 'NONE'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'default' => 'T','name' => 'DIPOLE'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'parameter'}],'type' => 'e','name' => 'if'},{'attrib' => {'expr' => '$TyepBField eq \'DIPOLE\''},'content' => [{'attrib' => {'min' => '0','max' => '180','type' => 'real','name' => 'MagAxisThetaGeo'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','max' => '360','type' => 'real','name' => 'MagAxisPhiGeo'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'type' => 'real','name' => 'DipoleStrength'},'content' => [],'type' => 'e','name' => 'parameter'}],'type' => 'e','name' => 'if'},{'attrib' => {'expr' => 'not $PlanetCommand'},'content' => [{'content' => '
		PLANET should precede $PlanetCommand
	','type' => 't'}],'type' => 'e','name' => 'rule'},{'content' => '

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
OmegaPlanet is the angular speed relative to an inertial frame, and
TiltRotation is the tilt of the rotation axis relative to ecliptic North,
TypeBField, which can be "NONE" or "DIPOLE". 
TypeBField="NONE" means that the planet does not have magnetic field. 
If TypeBField is set to "DIPOLE" then the following variables are read:
MagAxisThetaGeo and MagAxisPhiGeo are the colatitude and longitude
of the north magnetic pole in corotating planetocentric coordinates.
Finally DipoleStrength is the equatorial strength of the magnetic dipole
field. The units are indicated in the above example, which shows the
Earth values approximately.

The default value is NamePlanet="Earth", which is currently
the only recognized planet.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'if' => '$_IsStandAlone','name' => 'ROTATIONAXIS'},'content' => [{'attrib' => {'default' => 'T','type' => 'logical','name' => 'IsRotAxisPrimary'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'expr' => '$IsRotAxisPrimary'},'content' => [{'attrib' => {'min' => '0','max' => '180','type' => 'real','name' => 'RotAxisTheta'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','max' => '360','type' => 'real','name' => 'RotAxisPhi'},'content' => [],'type' => 'e','name' => 'parameter'}],'type' => 'e','name' => 'if'},{'attrib' => {'value' => 'ROTATIONAXIS','type' => 'string','name' => 'PlanetCommand'},'content' => [],'type' => 'e','name' => 'set'},{'content' => '

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
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'if' => '$_IsStandAlone','name' => 'ROTATION'},'content' => [{'attrib' => {'default' => 'T','type' => 'logical','name' => 'UseRotation'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'expr' => '$UseRotation'},'content' => [{'attrib' => {'type' => 'real','name' => 'RotationPeriod'},'content' => [],'type' => 'e','name' => 'parameter'}],'type' => 'e','name' => 'if'},{'attrib' => {'value' => 'MAGNETICAXIS','type' => 'string','name' => 'PlanetCommand'},'content' => [],'type' => 'e','name' => 'set'},{'content' => '

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
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'if' => '$_IsStandAlone','name' => 'MAGNETICAXIS'},'content' => [{'attrib' => {'default' => 'T','type' => 'logical','name' => 'IsMagAxisPrimary'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'expr' => '$IsMagAxisPrimary'},'content' => [{'attrib' => {'min' => '0','max' => '180','type' => 'real','name' => 'MagAxisTheta'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','max' => '360','type' => 'real','name' => 'MagAxisPhi'},'content' => [],'type' => 'e','name' => 'parameter'}],'type' => 'e','name' => 'if'},{'attrib' => {'value' => 'MAGNETICAXIS','type' => 'string','name' => 'PlanetCommand'},'content' => [],'type' => 'e','name' => 'set'},{'content' => '

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
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'if' => '$_IsStandAlone','name' => 'DIPOLE'},'content' => [{'attrib' => {'type' => 'real','name' => 'DipoleStrength'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

#DIPOLE
-3.11e-5		DipoleStrength [Tesla]

The DipoleStrength variable contains the
magnetic equatorial strength of the dipole magnetic field in Tesla.

The default value is the real dipole strength for the planet.
For the Earth the default is taken to be -31100 nT.
The sign is taken to be negative so that the magnetic axis can
point northward as usual.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'if' => '$_IsStandAlone','name' => 'UPDATEB0'},'content' => [{'attrib' => {'min' => '-1','default' => '0.0001','type' => 'real','name' => 'DtUpdateB0'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

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
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'if' => '$_IsStandAlone','name' => 'IDEALAXES'},'content' => [{'content' => '

#IDEALAXES

The #IDEALAXES command has no parameters. It sets both the rotational
and magnetic axes parallel with the ecliptic North direction. In fact
it is identical with the commands:

#ROTATIONAXIS
T               IsRotAxisPrimary
0.0             RotAxisTheta
0.0             RotAxisPhi

#MAGNETICAXIS
F               IsMagAxisPrimary

but much shorter.
','type' => 't'}],'type' => 'e','name' => 'command'}],'type' => 'e','name' => 'commandgroup'},{'attrib' => {'name' => 'USER DEFINED INPUT'},'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!  USER DEFINED INPUT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

','type' => 't'},{'attrib' => {'name' => 'USER_FLAGS'},'content' => [{'attrib' => {'default' => 'F','type' => 'logical','name' => 'UseUserInnerBcs'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => 'F','type' => 'logical','name' => 'UseUserSource'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => 'F','type' => 'logical','name' => 'UseUserPerturbation'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => 'F','type' => 'logical','name' => 'UseUserOuterBcs'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => 'F','type' => 'logical','name' => 'UseUserICs'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => 'F','type' => 'logical','name' => 'UseUserSpecifyRefinement'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => 'F','type' => 'logical','name' => 'UseUserLogFiles'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => 'F','type' => 'logical','name' => 'UseUserWritePlot'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => 'F','type' => 'logical','name' => 'UseUserAMR'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => 'F','type' => 'logical','name' => 'UseUserEchoInput'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => 'F','type' => 'logical','name' => 'UseUserB0'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => 'F','type' => 'logical','name' => 'UseUserSetPhysConst'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => 'F','type' => 'logical','name' => 'UseUserUpdateStates'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

#USER_FLAGS
F			UseUserInnerBcs
F			UseUserSource
F			UseUserPerturbation
F                       UseUserOuterBcs
F                       UseUserICs
F                       UseUserSpecifyRefinement
F                       UseUserLogFiles
F                       UseUserWritePlot
F                       UseUserAMR
F                       UseUserEchoInput
F                       UseUserB0
F                       UseUserSetPhysConst
F                       UseUserUpdateStates

This command controls the use of user defined routines in user_routines.f90.
For each flag that is set, an associated routine will be called in 
user_routines.f90.  Default is .false. for all flags.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'USERINPUTBEGIN'},'content' => [{'content' => '

This command signals the beginning of the section of the file which 
is read by the subroutine user\\_read\\_inputs in the user\\_routines.f90 file.
The section ends with the #USERINPUTEND command. There is no XML based parameter
checking in the user section.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'USERINPUTEND'},'content' => [{'content' => '

This command signals the end of the section of the file which 
is read by the subroutine user\\_read\\_inputs in the user\\_routines.f90 file.
The section begins with the #USERINPUTBEGIN command. There is no XML based parameter
checking in the user section.
','type' => 't'}],'type' => 'e','name' => 'command'}],'type' => 'e','name' => 'commandgroup'},{'attrib' => {'name' => 'TESTING AND TIMING'},'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  TESTING AND TIMING PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
','type' => 't'},{'attrib' => {'name' => 'TEST'},'content' => [{'attrib' => {'length' => '100','type' => 'string','name' => 'TestString'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '
#TEST
read_inputs

! A space separated list of subroutine names. Default is empty string.
!
! Examples:\\\\
!   read_inputs  - echo the input parameters following the #TEST line\\\\
!   project_B    - info on projection scheme\\\\   
!   implicit     - info on implicit scheme\\\\     
!   krylov       - info on the Krylov solver\\\\   
!   message_count- count messages\\\\
!   initial_refinement\\\\
!   ...
!
! Check the subroutines for call setoktest("...",oktest,oktest_me) to
! see the appropriate strings.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'TESTIJK'},'content' => [{'attrib' => {'min' => '-2','max' => '$nI+2','type' => 'integer','name' => 'iTest'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '-2','max' => '$nJ+2','type' => 'integer','name' => 'jTest'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '-2','max' => '$nK+2','type' => 'integer','name' => 'kTest'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '1','max' => '$MaxBlock','type' => 'integer','name' => 'iBlockTest'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','type' => 'integer','name' => 'iProcTest'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '
#TESTIJK
1                       iTest           (cell index for testing)
1                       jTest           (cell index for testing)
1                       kTest           (cell index for testing)
1                       iBlockTest      (block index for testing)
0                       iProcTest       (processor index for testing)

! The location of test info in terms of indices, block and processor number.
! Note that the user should set #TESTIJK or #TESTXYZ, not both.  If both
! are set, the final one in the session will set the test point.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'TESTXYZ'},'content' => [{'attrib' => {'min' => '$xMin','max' => '$xMax','type' => 'real','name' => 'xTest'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '$yMin','max' => '$yMax','type' => 'real','name' => 'yTest'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '$zMin','max' => '$zMax','type' => 'real','name' => 'zTest'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '
#TESTXYZ
1.5                     xTest           (X coordinate of cell for testing)
-10.5                   yTest           (Y coordinate of cell for testing)
-10.                    zTest           (Z coordinate of cell for testing)

! The location of test info in terms of coordinates.
! Note that the user should set #TESTIJK or #TESTXYZ, not both.  If both
! are set, the final one in the session will set the test point.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'TESTTIME'},'content' => [{'attrib' => {'min' => '-1','default' => '-1','type' => 'integer','name' => 'nIterTest'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '-1','default' => '1e30','type' => 'real','name' => 'TimeTest'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

#TESTTIME
-1                      nIterTest       (iteration number to start testing)
10.5                    TimeTest        (time to start testing in seconds)

! The time step and physical time to start testing.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'TESTVAR'},'content' => [{'attrib' => {'input' => 'select','type' => 'integer','name' => 'iVarTest'},'content' => [{'attrib' => {'value' => '1','default' => 'T','name' => 'Rho'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '2','name' => 'RhoUx'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '3','name' => 'RhoUy'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '4','name' => 'RhoUz'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '5','name' => 'Bx'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '6','name' => 'By'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '7','name' => 'Bz'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '8','name' => 'e'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '9','name' => 'p'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'parameter'},{'content' => '
#TESTVAR
1                       iVarTest

! Index of variable to be tested. Default is rho_="1", i.e. density.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'TESTDIM'},'content' => [{'attrib' => {'input' => 'select','type' => 'integer','name' => 'iVarTest'},'content' => [{'attrib' => {'value' => '0','name' => 'all'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '1','default' => 'T','name' => 'x'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '2','name' => 'y'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '3','name' => 'z'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'parameter'},{'content' => '
#TESTDIM
1                       iDimTest

! Index of dimension/direction to be tested. Default is X dimension.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'STRICT'},'content' => [{'attrib' => {'default' => 'T','type' => 'logical','name' => 'UseStrict'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '
#STRICT
T                       UseStrict

! If true then stop when parameters are incompatible. If false, try to
! correct parameters and continue. Default is true, i.e. strict mode
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'VERBOSE'},'content' => [{'attrib' => {'input' => 'select','type' => 'integer','name' => 'lVerbose'},'content' => [{'attrib' => {'value' => '-1','name' => 'errors and warnings only'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '0','name' => 'start and end of sessions'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '1','default' => 'T','name' => 'normal'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '10','name' => 'calls on test processor'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '100','name' => 'calls on all processors'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'parameter'},{'content' => '
#VERBOSE
-1                      lVerbose

! Verbosity level controls the amount of output to STDOUT. Default level is 1.
!\\\\
!   lVerbose $\\leq -1$ only warnings and error messages are shown.\\\\
!   lVerbose $\\geq  0$ start and end of sessions is shown.\\\\
!   lVerbose $\\leq  1$ a lot of extra information is given.\\\\
!   lVerbose $\\leq 10$ all calls of set_oktest are shown for the test processor.\\\\
!   lVerbose $\\leq 100$ all calls of set_oktest are shown for all processors.\\\\
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'DEBUG'},'content' => [{'attrib' => {'default' => 'F','type' => 'logical','name' => 'DoDebug'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => 'F','type' => 'logical','name' => 'DoDebugGhost'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '
#DEBUG
F                       DoDebug         (use it as if(okdebug.and.oktest)...)
F                       DoDebugGhost    (parameter for show_BLK in library.f90)

! Excessive debug output can be controlled by the global okdebug parameter
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'if' => '$_IsFirstSession','name' => 'CODEVERSION'},'content' => [{'attrib' => {'min' => '0','type' => 'real','default' => '7.50','name' => 'CodeVersion'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '
#CODEVERSION
7.50                    CodeVersion

! Checks CodeVersion. Prints a WARNING if it differs from the CodeVersion
! defined in ModMain.f90. Used in newer restart header files. 
! Should be given in PARAM.in when reading old restart files, 
! which do not have version info in the header file.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'if' => '$_IsFirstSession','name' => 'EQUATION'},'content' => [{'attrib' => {'length' => '100','default' => 'MHD','type' => 'string','name' => 'NameEquation'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '8','type' => 'integer','name' => 'nVar'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '
#EQUATION
MHD			NameEquation
8			nVar

! Define the equation name and the number of variables.
! If any of these do not agree with the values determined 
! by the code, BATSRUS stops with an error. Used in restart
! header files and can be given in PARAM.in as a check
! and as a description.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'if' => '$_IsFirstSession','name' => 'PRECISION'},'content' => [{'attrib' => {'input' => 'select','type' => 'integer','name' => 'nByteReal'},'content' => [{'attrib' => {'value' => '4','default' => '$_nByteReal==4','name' => 'single precision (4)'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '8','default' => '$_nByteReal==8','name' => 'double precision (8)'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'parameter'},{'attrib' => {'expr' => '$nByteReal==$_nByteReal'},'content' => [{'content' => '
		nByteReal in file must agree with _nByteReal.
	','type' => 't'}],'type' => 'e','name' => 'rule'},{'content' => '

#PRECISION
8                       nByteReal

! Define the number of bytes in a real number. If it does not agree
! with the value determined by the code, BATSRUS stops with an error.
! This is used in latest restart header files to check binary 
! compatibility between the restart file and the compiled code
! The command may also be given in PARAM.in to enforce a certain precision.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'if' => '$_IsFirstSession','name' => 'CHECKGRIDSIZE'},'content' => [{'attrib' => {'min' => '$nI','max' => '$nI','default' => '$nI','type' => 'integer','name' => 'nI'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '$nJ','max' => '$nJ','default' => '$nJ','type' => 'integer','name' => 'nJ'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '$nK','max' => '$nK','default' => '$nK','type' => 'integer','name' => 'nK'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '1','max' => '$MaxBlockALL','type' => 'integer','name' => 'MinBlockALL'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

#CHECKGRIDSIZE
       4                        nI
       4                        nJ
       4                        nK
     576                        MinBlockALL

! Checks block size and number of blocks. Stops with an error message,
! if nI, nJ, or nK differ from those set in ModSize. 
! Also stops if number_of_blocks exceeds nBLK*numprocs, where nBLK 
! is defined in ModSize and numprocs is the number of processors.
! This command is used in the restart headerfile to check consistency,
! and it is also useful to check if the executable is consistent with the 
! requirements of the problem described in the PARAM.in file.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'BLOCKLEVELSRELOADED'},'content' => [{'content' => '
#BLOCKLEVELSRELOADED

This command means that the restart file contains the information about
the minimum and maximum allowed refinement levels for each block.
This command is only used in the restart header file.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'if' => '$_IsStandAlone','name' => 'TIMING'},'content' => [{'attrib' => {'default' => 'T','type' => 'logical','name' => 'UseTiming'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'expr' => '$UseTiming'},'content' => [{'attrib' => {'input' => 'select','type' => 'integer','name' => 'Frequency'},'content' => [{'attrib' => {'value' => '-3','name' => 'none'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '-2','default' => 'T','name' => 'final only'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '-1','name' => 'end of sessions'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'min' => '1','default' => '100','name' => 'every X steps'},'content' => [],'type' => 'e','name' => 'optioninput'}],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '-1','default' => '-1','type' => 'integer','name' => 'nDepthTiming'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'input' => 'select','type' => 'string','name' => 'TypeTimingReport'},'content' => [{'attrib' => {'value' => 'cumu','default' => '1','name' => 'cumulative'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'list'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'tree'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'parameter'}],'type' => 'e','name' => 'if'},{'content' => '
#TIMING
T                       UseTiming      (rest of parameters read if true)
-2                      DnTiming       (-3 none, -2 final, -1 each session/AMR)
-1                      nDepthTiming   (-1 for arbitrary depth)
tree                    TypeTimingReport   (\'cumu\', \'list\', or \'tree\')

! The default values are shown.
!
! This command can only be used in stand alone mode. In the SWMF the
! #TIMING command should be issued for CON.
!
! If UseTiming=.true., the TIMING module must be on.
! If UseTiming=.false., the execution is not timed.
!
! Dntiming determines the frequency of timing reports.
! If DnTiming .ge.  1, a timing report is produced every dn_timing step.
! If DnTiming .eq. -1, a timing report is shown at the end of each session,
!                    before each AMR, and at the end of the whole run.
! If DnTiming .eq. -2, a timing report is shown at the end of the whole run.
! If DnTiming .eq. -3, no timing report is shown.
!
! nDepthTiming determines the depth of the timing tree. A negative number
! means unlimited depth. If TimingDepth is 1, only the full BATSRUS execution
! is timed.
!
! TypeTimingReport determines the format of the timing reports:
! \'cumu\' - cumulative list sorted by timings
! \'list\' - list based on caller and sorted by timings
! \'tree\' - tree based on calling sequence
','type' => 't'}],'type' => 'e','name' => 'command'}],'type' => 'e','name' => 'commandgroup'},{'attrib' => {'name' => 'INITIAL AND BOUNDARY CONDITIONS'},'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! MAIN INITIAL AND BOUNDARY CONDITION PARAMETERS  !!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
','type' => 't'},{'attrib' => {'required' => 'T','if' => '$_IsFirstSession','name' => 'PROBLEMTYPE'},'content' => [{'attrib' => {'input' => 'select','type' => 'integer','name' => 'iProblem'},'content' => [{'attrib' => {'value' => '1','name' => 'Uniform'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '2','name' => 'Shock tube'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '3','name' => 'Heliosphere'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '5','name' => 'Comet'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '6','name' => 'Rotation'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '7','name' => 'Diffusion'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '11','default' => 'T','name' => 'Earth'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '12','name' => 'Saturn'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '13','name' => 'Jupiter'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '14','name' => 'Venus'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '21','name' => 'Cylinder'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '22','name' => 'Sphere'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '25','name' => 'Arcade'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '26','name' => 'CME'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '30','name' => 'Dissipation'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'parameter'},{'attrib' => {'length' => '20','if' => '$iProblem==30','type' => 'string','name' => 'TypeDissipation'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '
#PROBLEMTYPE
30			iProblem
heat_test1		TypeProblemDiss

! select a problem type which defines defaults for a lot of parameters
!
! Problem type has to be defined as the first item after #TEST..#DEBUG items!
!\\begin{verbatim}
!       iProblem:     1=MHD Uniform Flow
!                     2=Shock tube
!                     3=Solar Wind and Inner Heliosphere
!                     5=Mass-Loaded Comet
!                     6=Rotation test
!                     7=Diffusion test
!                    11=Earth Magnetosphere
!                    12=Saturn Magnetosphere
!                    13=Jupiter Magnetosphere
!                    14=Venus Ionosphere
!                    21=Conducting Cylinder (2-D)
!                    22=Conducting Sphere   (3-D)
!                    25=Arcade
!                    26=CME
!		     30=Test Dissipative MHD
!\\end{verbatim}
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'multiple' => 'T','alias' => 'COORDINATESYSTEM','name' => 'COORDSYSTEM'},'content' => [{'attrib' => {'input' => 'select','type' => 'string','name' => 'TypeCoordSystem'},'content' => [{'attrib' => {'if' => '$_NameComp eq \'GM\'','default' => 'T','name' => 'GSM'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'if' => '$_NameComp eq \'IH\'','default' => 'T','name' => 'HGI'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'if' => '$_NameComp eq \'IH\'','name' => 'HGC'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'if' => '$_NameComp eq \'SC\'','default' => 'T','name' => 'HGR'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'if' => '$_NameComp eq \'SC\'','name' => 'HGI'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'parameter'},{'content' => '

#COORDSYSTEM
GSM			TypeCoordSystem

! TypeCoordSystem defines the coordinate system for the component.
! Currently only one coordinate system is available for GM ("GSM")
! and two for IH ("HGI" or "HGC") and two for SC ("HGR" or "HGI"). 
! In the future "GSE" should be also an option for GM.
! The coordinate systems are defined in share/Library/src/CON_axes.
!
! Default is component dependent: "GSM" for GM, "HGI" for IH, and "HGR" for SC.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'if' => '$_IsFirstSession','name' => 'RESTARTINDIR'},'content' => [{'attrib' => {'length' => '100','default' => 'GM/restartIN','type' => 'string','name' => 'NameRestartInDir'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'expr' => '-d $NameRestartInDir'},'content' => [{'content' => '
		Restart input directory $NameRestartInDir must exist!
	','type' => 't'}],'type' => 'e','name' => 'rule'},{'content' => '

#RESTARTINDIR
GM/restart_n5000	NameRestartInDir

! The NameRestartInDir variable contains the name of the directory
! where restart files are saved relative to the run directory.
! The directory should be inside the subdirectory with the name 
! of the component.
!
! Default value is "GM/restartIN".
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'if' => '$_IsFirstSession','name' => 'NEWRESTART'},'content' => [{'attrib' => {'default' => 'F','type' => 'logical','name' => 'DoRestartBFace'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => ' 

#NEWRESTART
T		DoRestartBFace 

! The RESTARTINDIR/restart.H file always contains the #NEWRESTART command.
! This command is really used only in the restart headerfile.  Generally
! it is not inserted in a PARAM.in file by the user.
!
! The #NEWRESTART command sets the following global variables:
! DoRestart=.true. (read restart files),
! DoRestartGhost=.false.  (no ghost cells are saved into restart file)
! DoRestartReals=.true.   (only real numbers are saved in blk*.rst files).

! The DoRestartBFace parameter tells if the face centered magnetic field
! is saved into the restart files. These values are used by the 
! Constrained Transport scheme.

','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'required' => 'T','if' => '$_IsFirstSession','name' => 'GRID'},'content' => [{'attrib' => {'min' => '1','default' => '2','type' => 'integer','name' => 'nRootBlockX'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '1','default' => '1','type' => 'integer','name' => 'nRootBlockY'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '1','default' => '1','type' => 'integer','name' => 'nRootBlockZ'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '-192.0','type' => 'real','name' => 'xMin'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '$xMin','default' => '  64.0','type' => 'real','name' => 'xMax'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => ' -64.0','type' => 'real','name' => 'yMin'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '$yMin','default' => '  64.0','type' => 'real','name' => 'yMax'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => ' -64.0','type' => 'real','name' => 'zMin'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '$zMin','default' => '  64.0','type' => 'real','name' => 'zMax'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '
#GRID
2                       nRootBlockX
1                       nRootBlockY
1                       nRootBlockZ
-224.                   xMin
 32.                    xMax
-64.                    yMin
 64.                    yMax
-64.                    zMin
 64.                    zMax

! The nRootBlockX, nRootBlockY and nRootBlockZ parameters define the 
! number of blocks of the base grid, i.e. the roots of the octree. 
! By varying these parameters, one can setup a grid which is elongated
! in some direction. The xMin, ..., zMax parameters define the physical
! size of the grid.
!
! There is no default value, the grid size must always be given.
! The #GRID command should be used before the #SAVEPLOT command.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'OUTERBOUNDARY'},'content' => [{'attrib' => {'values' => 'TypeBcEast,TypeBcWest,TypeBcSouth,TypeBcNorth,TypeBcBot,TypeBcTop','name' => 'Side'},'content' => [{'attrib' => {'input' => 'select','type' => 'string','name' => '$Side'},'content' => [{'attrib' => {'name' => 'coupled'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'default' => '$Side ne \'TypeBcEast\'','name' => 'fixed/inflow'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'default' => '$Side eq \'TypeBcEast\'','name' => 'float/outflow'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'heliofloat'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'reflect'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'periodic'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'vary'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'shear'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'linetied'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'raeder'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'arcadetop'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'arcadebot'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'arcadebotcont'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'user'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'parameter'}],'type' => 'e','name' => 'foreach'},{'attrib' => {'expr' => 'not($TypeBcEast eq \'periodic\' xor $TypeBcWest eq \'periodic\')'},'content' => [{'content' => '
		East and west BCs must be both periodic or neither
	','type' => 't'}],'type' => 'e','name' => 'rule'},{'attrib' => {'expr' => 'not($TypeBcSouth eq \'periodic\' xor $TypeBcNorth eq \'periodic\')'},'content' => [{'content' => '
		South and North BCs must be both periodic or neither
	','type' => 't'}],'type' => 'e','name' => 'rule'},{'attrib' => {'expr' => 'not($TypeBcBot eq \'periodic\' xor $TypeBcTop eq \'periodic\')'},'content' => [{'content' => '
		Bottom and top BCs must be both periodic or neither
	','type' => 't'}],'type' => 'e','name' => 'rule'},{'content' => '
#OUTERBOUNDARY
outflow                 TypeBcEast
inflow                  TypeBcWest
float                   TypeBcSouth
float                   TypeBcNorth
float                   TypeBcBottom
float                   TypeBcTop

! Default depends on problem type.\\\\
! Possible values:
!\\begin{verbatim}
! coupled       - GM coupled to the IH component (at the \'west\' boundary)
! fixed/inflow  - fixed solarwind values
! fixedB1       - fixed solarwind values without correction for the dipole B0
! float/outflow - zero gradient
! heliofloat    - floating for the SC component (requires #FACEOUTERBC)
! linetied      - float P, rho, and B, reflect all components of U
! raeder        - Jimmy Raeder\'s BC
! reflect       - reflective
! periodic      - periodic
! vary          - time dependent BC (same as fixed for non time_accurate)
! shear         - sheared (intended for shock tube problem only)
! arcadetop     - intended for arcade problem only
! arcadebot     - intended for arcade problem only
! arcadebotcont - intended for arcade problem only
! user          - user defined
!\\end{verbatim}
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'INNERBOUNDARY'},'content' => [{'attrib' => {'input' => 'select','type' => 'string','name' => 'TypeBcInner'},'content' => [{'attrib' => {'name' => 'reflect'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'float'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'fixed'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'default' => 'T','name' => 'ionosphere'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'ionosphereB0/ionosphereb0','name' => 'ionosphereB0'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'ionospherefloat'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'coronatoih'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'user'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'parameter'},{'attrib' => {'expr' => '$UseBody2'},'content' => [{'attrib' => {'input' => 'select','type' => 'string','name' => 'TypeBcBody2'},'content' => [{'attrib' => {'default' => 'T','name' => 'reflect'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'float'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'fixed'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'ionosphere'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'ionosphereB0'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'ionospherefloat'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'parameter'}],'type' => 'e','name' => 'if'},{'attrib' => {'expr' => 'not($TypeBcBody2 =~ /ionosphere/)'},'content' => [{'content' => ' 
	For the second body COROTATION AND AN IONOSPHERIC BOUNDARY DO NOT WORK.
	','type' => 't'}],'type' => 'e','name' => 'rule'},{'content' => '

#INNERBOUNDARY
ionosphere              TypeBcInner
ionosphere              TypeBcBody2  !read only if UseBody2=T 


This command should appear after the #SECONDBODY command when using 
two bodies. Note: for the second body COROTATION AND AN IONOSPHERIC 
BOUNDARY DO NOT WORK.

Default value for TypeBcBody2 is \'reflect\'.


Possible values for TypeBcInner are:
\\begin{verbatim}
\'reflect\'         - reflect Vr, reflect Vphi to rotation, float Vtheta,
                    reflect Br, float Bphi, float Btheta, float rho, float P
\'float\'           - float Vr, reflect Vphi to rotation, float Vtheta,
                    float B, float rho, float P
\'fixed\'           - Vr=0, Vphi=rotation, Vtheta=0
                    B=B0 (ie B1=0), fix rho, fix P
\'ionosphere\'      - set V as if ionosphere gave V_iono=0
                    float B, fix rho, fix P
\'ionospherefloat\' - set V as if ionosphere gave V_iono=0
                    float B, float rho, float P
\'coronatoih\'      - IH component obtains inner boundary from the SC component
\'user\'            - user defined
\\end{verbatim}
For \'ionosphere\' and \'ionospherefloat\' types and a coupled GM-IE run,
the velocity at the inner boundary is determined by the ionosphere model.

Default value for TypeBcInner is \'ionosphere\' for problem types
Earth, Saturn, Jupiter, and rotation.
For all other problems with an inner boundary the default is \'unknown\',
so the inner boundary must be set.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'EXTRABOUNDARY'},'content' => [{'attrib' => {'default' => 'F','type' => 'logical','name' => 'UseExtraBoundary'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'expr' => '$UseExtraBoundary'},'content' => [{'attrib' => {'type' => 'string','name' => 'TypeBcExtra'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => 'F','type' => 'logical','name' => 'DoFixExtraboundary'},'content' => [],'type' => 'e','name' => 'parameter'}],'type' => 'e','name' => 'if'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'FACEOUTERBC'},'content' => [{'attrib' => {'min' => '0','max' => '6','default' => '0','type' => 'integer','name' => 'MaxBoundary'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'expr' => '$MaxBoundary >= 1'},'content' => [{'attrib' => {'default' => 'F','type' => 'logical','name' => 'DoFixOuterBoundary'},'content' => [],'type' => 'e','name' => 'parameter'}],'type' => 'e','name' => 'if'},{'content' => '
#FACEOUTERBC
0              MaxBoundary            
F              DoFixOuterBoundary)    !read only for MaxBoundary>=East_(=1).
! If MaxBoundary is East_(=1) or more then the outer boundaries with
! the number of boundary being between East_ and MaxBoundary
! are treated using set_BCs.f90 subroutines instead of set_outerBCs.f90 
! if DoFixOuterBoundary is .true., there is no resolution
! change along the outer boundaries with the number of
! of boundary being between East_ and MaxBoundary
','type' => 't'}],'type' => 'e','name' => 'command'}],'type' => 'e','name' => 'commandgroup'},{'attrib' => {'name' => 'INITIAL TIME'},'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! INITIAL TIME AND STEP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
','type' => 't'},{'attrib' => {'alias' => 'SETREALTIME','if' => '$_IsFirstSession','name' => 'STARTTIME'},'content' => [{'attrib' => {'default' => '2000','type' => 'integer','name' => 'iYear'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '1','max' => '12','default' => '3','type' => 'integer','name' => 'iMonth'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '1','max' => '31','default' => '21','type' => 'integer','name' => 'iDay'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','max' => '23','default' => '0','type' => 'integer','name' => 'iHour'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','max' => '59','default' => '0','type' => 'integer','name' => 'iMinute'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','max' => '59','default' => '0','type' => 'integer','name' => 'iSecond'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '
#STARTTIME
2000                    iYear
3                       iMonth
21                      iDay
10                      iHour
45                      iMinute
0                       iSecond

The #STARTTIME command sets the initial date and time for the
simulation in Greenwich Mean Time (GMT) or Universal Time (UT)
in stand alone mode. 
In the SWMF this command checks start times against the SWMF start time 
and warns if the difference exceeds 1 millisecond.
This time is stored in the BATSRUS restart header file.

The default values are shown above.
This is a date and time when both the rotational and the magnetic axes
have approximately zero tilt towards the Sun.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'if' => '$_IsFirstSession','name' => 'TIMESIMULATION'},'content' => [{'attrib' => {'min' => '0','default' => '0.0','type' => 'real','name' => 'tSimulation'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

#TIMESIMULATION
3600.0			tSimulation [sec]

The tSimulation variable contains the simulation time in seconds
relative to the initial time set by the #STARTTIME command.
The #TIMESIMULATION command and tSimulation are saved into the restart 
header file, which provides human readable information about the restart state.

In SWMF the command is ignored (SWMF has its own #TIMESIMULATION command).
In stand alone mode time\\_simulation is set, but in case of a restart,
it gets overwritten by the binary value saved into the .rst binary files. 

The default value is tSimulation=0.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'if' => '$_IsFirstSession','name' => 'NSTEP'},'content' => [{'attrib' => {'min' => '0','default' => '0','type' => 'integer','name' => 'nStep'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

#NSTEP
100			nStep

! Set nStep for the component. Typically used in the restart.H header file.
! Generally it is not inserted in a PARAM.in file by the user.
!
! The default is nStep=0 as the starting time step with no restart.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'if' => '$_IsFirstSession','name' => 'NPREVIOUS'},'content' => [{'attrib' => {'min' => '-1','default' => '-1','type' => 'integer','name' => 'nPrevious'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

#NPREVIOUS
100			nPrev
1.5			DtPrev

! This command should only occur in the restart.H header file.
! If it is present, it indicates that the restart file contains
! the state variables for the previous time step.
! nPrev is the time step number and DtPrev is the length of the previous 
! time step in seconds.
! The previous time step is needed for a second order in time restart 
! with the implicit scheme. 
!
! The default is that the command is not present and no previous time step 
! is saved into the restart files.
','type' => 't'}],'type' => 'e','name' => 'command'}],'type' => 'e','name' => 'commandgroup'},{'attrib' => {'name' => 'TIME INTEGRATION'},'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  TIME INTEGRATION PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
','type' => 't'},{'attrib' => {'name' => 'TIMESTEPPING'},'content' => [{'attrib' => {'input' => 'select','type' => 'integer','name' => 'nStage'},'content' => [{'attrib' => {'value' => '1','default' => 'T'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '2'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','max' => '1','default' => '0.8','type' => 'real','name' => 'CflExpl'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

#TIMESTEPPING
2                       nStage
0.80                    CflExpl

! Parameters for explicit time integration.
! Default is 1 stage and CflExpl=0.8
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'FIXEDTIMESTEP'},'content' => [{'attrib' => {'default' => 'F','type' => 'logical','name' => 'UseDtFixed'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','if' => '$UseDtFixed','default' => '1.0','type' => 'real','name' => 'DtFixedDim'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '
#FIXEDTIMESTEP
T                       UseDtFixed
10.                     DtFixedDim [sec] (read if UseDtFixed is true)

! Default is UseDtFixed=.false. Effective only if DoTimeAccurate is true.
! If UseDtFixed is true, the time step is fixed to DtFixedDim.
!
! This is useful for debugging explicit schemes.

! The real application is, however, for implicit and partially
! implicit/local schemes. The time step is set to DtFixedDim unless the
! update checking decides to reduce the time step for the sake of robustness.

','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'PARTSTEADY'},'content' => [{'attrib' => {'default' => 'F','type' => 'logical','name' => 'UsePartSteady'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

! If UsePartSteady is true, the partially steady state algorithm is used.
! Only blocks which are changing or next to changing blocks are evolved.
! This scheme can speed up the calculation if part 
! of the domain is in a numerical steady state. 
! In steady state runs the code stops when a full steady state is
! achieved. The conditions for checking the numerical steady state are set 
! by the #PARTSTEADYCRITERIA command.
! Default value is UsePartSteady = .false.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'alias' => 'STEADYCRITERIA','name' => 'PARTSTEADYCRITERIA'},'content' => [{'attrib' => {'min' => '1','default' => '1','type' => 'integer','name' => 'MinCheckVar'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '$MinCheckVar','default' => '8','type' => 'integer','name' => 'MaxCheckVar'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'to' => '$MaxCheckVar','from' => '$MinCheckVar'},'content' => [{'attrib' => {'min' => '0','max' => '1','default' => '0.001','type' => 'real','name' => 'RelativeEps'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','default' => '0.0001','type' => 'real','name' => 'AbsoluteEps'},'content' => [],'type' => 'e','name' => 'parameter'}],'type' => 'e','name' => 'for'},{'content' => '
#PARTSTEADYCRITERIA
5               MinCheckVar
8               MaxCheckVar
0.001           RelativeEps(bx)
0.0001          AbsoluteEps(bx)
0.001           RelativeEps(by)
0.0001          AbsoluteEps(by)
0.001           RelativeEps(bz)
0.0001          AbsoluteEps(bz)
0.001           RelativeEps(p)
0.0001          AbsoluteEps(p)

The part steady scheme only evolves blocks which are changing,
or neighbors of changing blocks. The scheme checks the neighbor blocks
every time step if their state variable has changed significantly.
This command allows the user to select the variables to be checked,
and to set the relative and absolute limits for each variable.
Only the state variables indexed from MinCheckVar to MaxCheckVar are checked.
The change in the block is significant if 
\\begin{verbatim}
max(abs(State - StateOld)) / (RelativeEps*abs(State) + AbsoluteEps)
\\end{verbatim}
exceeds 1.0 for any of the checked variables in any cells of the block.
(including body cells but excluding ghost cells).
The RelativeEps variable determines the maximum ratio of the change
and the norm of the old state. The AbsoluteEps variable is only needed
if the old state is very close to zero. It should be set to a positive
value which is much smaller than the typical significantly non-zero
value of the variable.

Default values are such that all variables are checked with
relative error 0.001 and absolute error 0.0001.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'PARTLOCAL'},'content' => [{'attrib' => {'default' => 'F','type' => 'logical','name' => 'UsePartLocal'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '
#PARTLOCAL
T               UsePartLocal

! Default is UsePartLocal=.false. If UsePartLocal is true and the
! run is time accurate, then the blocks selected as "implicit"
! by the criteria defined in #STEPPINGCRITERIA are not used to
! calculate the time step, and all cells are advanced with the
! smaller of the stable and the global time steps.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'IMPLICIT'},'content' => [{'attrib' => {'default' => 'F','type' => 'logical','name' => 'UsePointImplicit'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => 'F','type' => 'logical','name' => 'UsePartImplicit'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => 'F','type' => 'logical','name' => 'UseFullImplicit'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','if' => '$UsePartImplicit or $UseFullImplicit','default' => '100','type' => 'real','name' => 'CflImpl'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'expr' => '$UsePointImplicit + $UsePartImplicit + $UseFullImplicit <= 1'},'content' => [{'content' => '
		At most one of these logicals can be true!
	','type' => 't'}],'type' => 'e','name' => 'rule'},{'content' => '

#IMPLICIT
F               UsePointImplicit   
F               UsePartImplicit
F               UseFullImplicit
100.0           CflImpl (read if UsePartImplicit or UseFullImplicit is true)

! Default is false for all logicals. Only one of them can be set to true!
! The CFL number is used in the implicit blocks of the fully or partially
! implicit schemes. Ignored if UseDtFixed is true.
','type' => 't'}],'type' => 'e','name' => 'command'}],'type' => 'e','name' => 'commandgroup'},{'attrib' => {'name' => 'IMPLICIT PARAMETERS'},'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! PARAMETERS FOR FULL AND PART IMPLICIT TIME INTEGRATION !!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
','type' => 't'},{'attrib' => {'alias' => 'STEPPINGCRITERIA','name' => 'IMPLICITCRITERIA'},'content' => [{'attrib' => {'input' => 'select','type' => 'string','name' => 'TypeImplCrit'},'content' => [{'attrib' => {'value' => 'dt','default' => 'T','name' => 'Time step'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'r/R','name' => 'Radial distance'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'test','name' => 'Test block'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','if' => '$TypeImplCrit eq \'R\'','type' => 'real','name' => 'rImplicit'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

#IMPLICITCRITERIA
r		TypeImplCrit (dt or r or test)
10.0		rImplicit    (only read for TypeImplCrit = r)

! Both #IMPLICITCRITERIA and #STEPPINGCRITERIA are acceptable.
! Only effective if PartImplicit or PartLocal is true in a time accurate run.
! Default value is ImplCritType=\'dt\'.
!
! The options are
!\\begin{verbatim}
! if     (TypeImplCrit ==\'dt\'  ) then blocks with DtBLK .gt. DtFixed
! elseif (TypeImplCrit ==\'r\'   ) then blocks with rMinBLK .lt. rImplicit
! elseif (TypeImplCrit ==\'test\') then block iBlockTest on processor iProcTest
!\\end{verbatim}
! and are handled with local/implicit scheme. 
! Here DtBlock is the time step
! allowed by the CFL condition for a given block, while rMinBLK is the
! smallest radial distance for all the cells in the block.\\\\
!
! \\noindent
! The iBlockTest and iProcTest can be defined in the #TESTIJK command.\\\\
! DtFixed must be defined in the #FIXEDTIMESTEP command.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'alias' => 'IMPLICITSTEP','name' => 'IMPLSTEP'},'content' => [{'attrib' => {'min' => '0','max' => '1','default' => '1','type' => 'real','name' => 'ImplCoeff'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => 'T','type' => 'logical','name' => 'UseBdf2'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => 'F','type' => 'logical','name' => 'UseSourceImpl'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '
! For steady state run the default values are shown. For second order
! time accurate run the default is UseBdf2=T, since
! BDF2 is a 3 level second order stable implicit scheme.
! This can be overwritten with #IMPLSTEP after the #TIMESTEPPING command.
! For example one could use the 2-level trapezoid scheme with
! ImplCoeff=0.5 and UseBDF2=F. ImplCoeff is the coefficient for $R^{n+1}$.
! For BDF2 scheme ImplCoeff is used in the first time step only, later on it
! is overwritten by the BDF2 scheme.
! UseSourceImpl true means that the preconditioner should take point
! source terms into account. Default is false.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'alias' => 'IMPLICITSCHEME','name' => 'IMPLSCHEME'},'content' => [{'attrib' => {'input' => 'select','type' => 'integer','name' => 'nOrderImpl'},'content' => [{'attrib' => {'default' => 'T','name' => '1'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => '2'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'parameter'},{'attrib' => {'input' => 'select','type' => 'string','name' => 'TypeFluxImpl'},'content' => [{'attrib' => {'value' => 'Roe/roe/1','name' => 'Roe'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'Rusanov/rusanov/2/TVDLF','default' => 'T','name' => 'Rusanov'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'Linde/linde/3/HLLEL','name' => 'Linde'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'Sokolov/sokolov/4/AW','name' => 'Sokolov'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'parameter'},{'content' => '
#IMPLSCHEME
1               nOrderImpl
Rusanov         TypeFluxImpl

! This command defines the scheme used in the implicit part (\'left hand side\').
! The default order is first order. The default scheme is the same as the
! scheme selected for the explicit part. 
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'alias' => 'IMPLICITCHECK','name' => 'IMPLCHECK'},'content' => [{'attrib' => {'min' => '0','max' => '0.9','default' => '0.3','type' => 'real','name' => 'RejectStepLevel'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','max' => '0.9','default' => '0.5','type' => 'real','name' => 'RejectStepFactor'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','max' => '0.9','default' => '0.6','type' => 'real','name' => 'ReduceStepLevel'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','max' => '1','default' => '0.9','type' => 'real','name' => 'ReduceStepFactor'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','max' => '1','default' => '0.8','type' => 'real','name' => 'IncreaseStepLevel'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '1','max' => '2','default' => '1.05','type' => 'real','name' => 'IncreaseStepFactor'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

#IMPLCHECK
0.3		RejectStepLevel
0.5		RejectStepFactor
0.6		ReduceStepLevel
0.95		ReduceStepFactor
0.8		IncreaseStepLevel
1.05		IncreaseStepFactor

The update checking of the implicit scheme can be tuned with this command.
Update checking is done unless it is switched off (see UPDATECHECK command).
After each (partially) implicit time step, the code computes pRhoRelMin,
which is the minimum of the relative pressure and density drops over 
the whole computational domain. The algorithm is the following:

If pRhoRelMin is less than RejectStepLevel,
the step is rejected, and the time step is reduced by RejectStepFactor;
else if pRhoRelMin is less than ReduceStepLevel,
the step is accepted, but the next time step is reduced by ReduceStepFactor;
else if pRhoRelMin is greater than IncreaseStepFactor,
the step is accepted and the next time step is increased by IncreaseStepFactor,
but it is never increased above the value given in the FIXEDTIMESTEP command.

Assigning ReduceStepFactor=1.0 means that the
time step is not reduced unless the step is rejected.
Assigning IncreaseStepFactor=1.0 means that the 
time step is never increased, only reduced.

Default values are shown.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'NEWTON'},'content' => [{'attrib' => {'default' => 'F','type' => 'logical','name' => 'UseConservativeImplicit'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => 'F','type' => 'logical','name' => 'UseNewton'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'expr' => '$UseNewton'},'content' => [{'attrib' => {'default' => 'T','type' => 'logical','name' => 'UseNewMatrix'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '1','default' => '10','type' => 'integer','name' => 'MaxIterNewton'},'content' => [],'type' => 'e','name' => 'parameter'}],'type' => 'e','name' => 'if'},{'content' => '
#NEWTON
F		UseConservativeImplicit
T               UseNewton
F               UseNewMatrix  (only read if UseNewton is true)
10              MaxIterNewton (only read if UseNewton is true)

! Default is UseConservativeImplicit=F and UseNewton=F, i.e. 
! no conservative fix is used and only one "Newton" iteration is done.
! UseNewMatrix decides whether the Jacobian should be recalculated
! for every Newton iteration. MaxIterNewton is the maximum number
! of Newton iterations before giving up.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'JACOBIAN'},'content' => [{'attrib' => {'input' => 'select','type' => 'string','name' => 'TypeJacobian'},'content' => [{'attrib' => {'value' => 'prec','default' => 'T','name' => 'Preconditioned'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'free','name' => 'No preconditioning'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','max' => '1.e-5','default' => '$doublePrecision ? 1.e-12 : 1.e-6','type' => 'real','name' => 'JacobianEps'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '
#JACOBIAN
prec            TypeJacobian (prec, free)
1.E-12          JacobianEps

! The Jacobian matrix is always calculated with a matrix free approach,
! however it can be preconditioned  (\'prec\'), or not (\'free\').  The
! Default value is TypeJacobian=\'prec\'.
! JacobianEps contains the machine round off error for numerical derivatives.
! The default value is 1.E-12 for 8 byte reals and 1.E-6 for 4 byte reals.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'PRECONDITIONER'},'content' => [{'attrib' => {'input' => 'select','type' => 'string','name' => 'TypePrecondSide'},'content' => [{'attrib' => {'name' => 'left'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'default' => 'T','name' => 'symmetric'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'right'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'parameter'},{'attrib' => {'input' => 'select','type' => 'string','name' => 'TypePrecond'},'content' => [{'attrib' => {'default' => 'T','name' => 'MBILU'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','max' => '1','default' => '0.5','type' => 'real','name' => 'GustafssonPar'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '
#PRECONDITIONER
symmetric       TypePrecondSide (left, symmetric, right)
MBILU           TypePrecond (MBILU)
0.5             GustafssonPar (0. no modification, 1. full modification)

! Default parameters are shown. Right preconditioning does not affect
! the normalization of the residual. The Gustafsson parameter determines
! how much the MBILU preconditioner is modified. The default 0.5 value
! means a relaxed modification.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'KRYLOV'},'content' => [{'attrib' => {'input' => 'select','type' => 'string','name' => 'TypeKrylov'},'content' => [{'attrib' => {'default' => 'T','name' => 'gmres'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'bicgstab'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'parameter'},{'attrib' => {'input' => 'select','type' => 'string','name' => 'TypeInitKrylov'},'content' => [{'attrib' => {'value' => 'nul','default' => 'T','name' => '0'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'old','name' => 'previous'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'explicit'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'explicit','name' => 'scaled explicit'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','max' => '0.1','default' => '0.001','type' => 'real','name' => 'ErrorMaxKrylov'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '1','default' => '100','type' => 'integer','name' => 'MaxMatvecKrylov'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '
#KRYLOV
gmres           TypeKrylov  (gmres, bicgstab)
nul             TypeInitKrylov (nul, old, explicit, scaled)
0.001           ErrorMaxKrylov
100             MaxMatvecKrylov

! Default values are shown. Initial guess for the Krylov type iterative scheme
! can be 0 (\'nul\'), the previous solution (\'old\'), the explicit solution
! (\'explicit\'), or the scaled explicit solution (\'scaled\'). The iterative
! scheme stops if the required accuracy is achieved or the maximum number
! of matrix-vector multiplications is exceeded.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'KRYLOVSIZE'},'content' => [{'attrib' => {'min' => '1','default' => 'MaxMatvecKrylov','type' => 'integer','name' => 'nKrylovVector'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '
#KRYLOVSIZE
10		nKrylovVector

! The number of Krylov vectors only matters for GMRES (TypeKrylov=\'gmres\').
! If GMRES does not converge within nKrylovVector iterations, it needs
! a restart, which usually degrades its convergence rate and robustness.
! So nKrylovVector should exceed the number of iterations, but
! it should not exceed the maximum number of iterations MaxMatvecKrylov.
! On the other hand the dynamically allocated memory is also proportional 
! to nKrylovVector. The default is nKrylovVector=MaxMatvecKrylov (in #KRYLOV)
! which can be overwritten by #KRYLOVSIZE after the #KRYLOV command (if any).
','type' => 't'}],'type' => 'e','name' => 'command'}],'type' => 'e','name' => 'commandgroup'},{'attrib' => {'name' => 'STOPPING CRITERIA'},'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! STOPPING CRITERIA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

The commands in this group only work in stand alone mode.

','type' => 't'},{'attrib' => {'required' => '$_IsStandAlone','if' => '$_IsStandAlone','name' => 'STOP'},'content' => [{'attrib' => {'min' => '-1','default' => '-1','type' => 'integer','name' => 'MaxIteration'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '-1','default' => '-1','type' => 'real','name' => 'tSimulationMax'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

#STOP
100			MaxIteration
10.0			tSimulationMax [sec]

This command is only used in stand alone mode.

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
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'if' => '$_IsStandAlone','name' => 'CHECKSTOPFILE'},'content' => [{'attrib' => {'default' => 'T','type' => 'logical','name' => 'DoCheckStopFile'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

#CHECKSTOPFILE
T			DoCheckStopFile

This command is only used in stand alone mode.

If DoCheckStopFile is true then the code checks if the
BATSRUS.STOP file exists in the run directory. This file is deleted at
the beginning of the run, so the user must explicitly create the file
with e.g. the "touch BATSRUS.STOP" UNIX command.
If the file is found in the run directory,
the execution stops in a graceful manner.
Restart files and plot files are saved as required by the
appropriate parameters.

The default is DoCheckStopFile=.true.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'if' => '$_IsStandAlone','name' => 'CPUTIMEMAX'},'content' => [{'attrib' => {'min' => '-1','default' => '-1','type' => 'real','name' => 'CpuTimeMax'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

#CPUTIMEMAX
3600                    CpuTimeMax [sec]

This command is only used in stand alone mode.

The CpuTimeMax variable contains the maximum allowed CPU time (wall clock
time) for the execution of the current run. If the CPU time reaches
this time, the execution stops in a graceful manner.
Restart files and plot files are saved as required by the
appropriate parameters.
This command is very useful when the code is submitted to a batch
queue with a limited wall clock time.

The default value is -1.0, which means that the CPU time is not checked.
To do the check the CpuTimeMax variable has to be set to a positive value.
','type' => 't'}],'type' => 'e','name' => 'command'}],'type' => 'e','name' => 'commandgroup'},{'attrib' => {'name' => 'OUTPUT PARAMETERS'},'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  OUTPUT PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
','type' => 't'},{'attrib' => {'name' => 'RESTARTOUTDIR'},'content' => [{'attrib' => {'length' => '100','default' => 'GM/restartOUT','type' => 'string','name' => 'NameRestartOutDir'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'expr' => '-d $NameRestartOutDir'},'content' => [{'content' => '
		Restart output directory $NameRestartOutDir must exist
	','type' => 't'}],'type' => 'e','name' => 'rule'},{'content' => '

#RESTARTOUTDIR
GM/restart_n5000	NameRestartOutDir

! The NameRestartOutDir variable contains the name of the directory
! where restart files are saved relative to the run directory.
! The directory should be inside the subdirectory with the name 
! of the component.
!
! Default value is "GM/restartOUT".
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'SAVERESTART'},'content' => [{'attrib' => {'default' => 'T','type' => 'logical','name' => 'DoSaveRestart'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'expr' => '$DoSaveRestart'},'content' => [{'attrib' => {'min' => '-1','default' => '-1','type' => 'integer','name' => 'DnSaveRestart'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '-1','default' => '-1','type' => 'real','name' => 'DtSaveRestart'},'content' => [],'type' => 'e','name' => 'parameter'}],'type' => 'e','name' => 'if'},{'content' => '
#SAVERESTART
T			DoSaveRestart Rest of parameters read if true
100			DnSaveRestart
-1.			DtSaveRestart [seconds]

! Default is DoSaveRestart=.true. with DnSaveRestart=-1 and 
! DtSaveRestart=-1. This results in the restart file being 
! saved only at the end.  A binary restart file is produced for every 
! block and named as RESTARTOUTDIR/blkGLOBALBLKNUMBER.rst.
! In addition the grid is described by RESTARTOUTDIR/octree.rst
! and an ASCII header file is produced with timestep and time info:
! RESTARTOUTDIR/restart.H
!
! The restart files are overwritten every time a new restart is done,
! but one can change the name of the RESTARTOUTDIR with the #RESTARTOUTDIR
! command from session to session. The default directory name is \'restartOUT\'.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'PLOTDIR'},'content' => [{'attrib' => {'length' => '100','default' => 'GM/IO2','type' => 'string','name' => 'NamePlotDir'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'expr' => '-d $NamePlotDir'},'content' => [{'content' => '
		Plot directory $NamePlotDir must exist
	','type' => 't'}],'type' => 'e','name' => 'rule'},{'content' => '

The NamePlotDir variable contains the name of the directory
where plot files and logfiles are saved relative to the run directory.
The directory should be inside the subdirectory with the name
of the component.

Default value is "GM/IO2".
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'SAVELOGFILE'},'content' => [{'attrib' => {'default' => 'F','type' => 'logical','name' => 'DoSaveLogfile'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'expr' => '$DoSaveLogfile'},'content' => [{'attrib' => {'min' => '1','max' => '4','type' => 'strings','name' => 'StringLog'},'content' => [{'attrib' => {'input' => 'select','required' => 'T','type' => 'string','name' => 'TypeLogVar'},'content' => [{'attrib' => {'value' => 'MHD','default' => 'T','name' => 'MHD vars. dimensional'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'FLX','name' => 'Flux vars. dimensional'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'RAW','name' => 'Raw vars. dimensional'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'VAR','name' => 'Set vars. dimensional'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'mhd','default' => 'T','name' => 'MHD vars. scaled'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'flx','name' => 'Flux vars. scaled'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'raw','name' => 'Raw vars. scaled'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'var','name' => 'Set vars. scaled'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'part'},{'attrib' => {'multiple' => 'T','input' => 'select','required' => 'F','type' => 'string','name' => 'TypeTime'},'content' => [{'attrib' => {'exclusive' => 'T','name' => 'none'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'step'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'date'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'time'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'part'}],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '-1','default' => '1','type' => 'integer','name' => 'DnSaveLogfile'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '-1','default' => '-1','type' => 'real','name' => 'DtSaveLogfile'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'length' => '100','if' => '$TypeLogVar =~ /var/i','type' => 'string','name' => 'NameLogVars'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '1','length' => '100','max' => '10','if' => '($TypeLogVar=~/flx/i or $NameLogVars=~/flx/i)','type' => 'strings','name' => 'StringLogRadii'},'content' => [{'attrib' => {'min' => '$rBody','multiple' => 'T','type' => 'real','name' => 'LogRadii'},'content' => [],'type' => 'e','name' => 'part'}],'type' => 'e','name' => 'parameter'}],'type' => 'e','name' => 'if'},{'content' => '
#SAVELOGFILE
T                       DoSaveLogfile, rest of parameters read if true
VAR step date           StringLog
100                     DnSaveLogfile
-1.                     DtSaveLogfile [sec]
rho p rhoflx            NameLogVars (read if StrigLog is \'var\' or \'VAR\')
4.0  10.0               rLog  (radii for the flux. Read if vars include \'flx\')

! Default is DoSaveLogfile=.false.
! The logfile can contain averages or point values and other scalar
! quantities.  It is written into an ASCII file named as\\\\
!
! NAMEPLOTDIR/log_TIMESTEP.log\\\\
!
! \\noindent
! where NAMEPLOTDIR can be defined with the #PLOTDIR command (default is IO2).\\\\
! The StringLog can contain two groups of information in arbitrary order.
! The first is LogVar which is a single 3 character string that indicates
! the type of variables that are to be written.  The second group indicates
! the type of time/iteration output format to use.  This second group is
! not required and defaults to something standard for each logvar case.\\\\
! Any of the identifiers for the timetype can be included in arbitrary order.
!
!\\begin{verbatim}
! logvar   = \'mhd\', \'raw\', \'flx\' or \'var\' - unitless output
! logvar   = \'MHD\', \'RAW\', \'FLX\' or \'VAR\' - dimensional output
! timetype = \'none\', \'step\', \'time\', \'date\'
!\\end{verbatim}
!
! The logvar string is not optional and must be found on the line.
! The timetype is optional - when not specified a logical choice is made
!  by the code.
!
! The log_var string defines the variables to print in the log file
! It also controls whether or not the variables will come out in
! dimensional or non-dimensional form by the capitalization of the log_var
! string.
!\\begin{verbatim}
! ALL CAPS  - dimensional
! all lower - dimensionless
!
! \'raw\' - vars: dt rho rhoUx rhoUy rhoUz Bx By Bz E Pmin Pmax
!       - time: step time
! \'mhd\' - vars: rho rhoUx rhoUy rhoUz Bx By Bz E Pmin Pmax
!       - time: step date time
! \'flx\' - vars: rho Pmin Pmax rhoflx pvecflx e2dflx
!       - time: step date time
! \'var\' - vars: READ FROM PARAMETER FILE
!       - time: step time
!\\end{verbatim}
! log_vars is read only when the log_string contains var or VAR.  The choices
! for variables are currently:
!\\begin{verbatim}
! Average value on grid:   rho rhoUx rhoUy rhoUz Ux Uy Uz Bx By Bz P E
! Value at the test point: rhopnt rhoUxpnt rhoUypnt rhoUxpnt Uxpnt Uypnt Uzpnt
!                          Bxpnt Bypnt Bzpnt B1xpnt B1ypnt B1zpnt
!                          Epnt Ppnt Jxpnt Jypnt Jzpnt
!                          theta1pnt theta2pnt phi1pnt phi2pnt statuspnt
! Ionosphere values:       cpcpn cpcps                  
!
! Max or Min on grid:  Pmin Pmax
! Flux values:         Aflx rhoflx Bflx B2flx pvecflx e2dflx
! Other variables:     dt
!\\end{verbatim}
! timetype values mean the following:
!\\begin{verbatim}
!  none  = there will be no indication of time in the logfile (not even an
!                # of steps)
!  step  = # of time steps (n_steps)
!  date  = time is given as an array of 7 integers:  year mo dy hr mn sc msc
!  time  = time is given as a real number - elapsed time since the start of
!          the run.  Units are determined by log_var and unitUSER_t
!\\end{verbatim}
! these can be listed in any combination in the log_string line.\\\\
! R_log is read only when one of the variables used is a \'flx\' variable.  R_log
! is a list of radii at which to calculate the flux through a sphere.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'if' => '$_IsFirstSession','name' => 'SATELLITE'},'content' => [{'attrib' => {'min' => '0','default' => '0','type' => 'integer','name' => 'nSatellite'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'to' => '$nSatellite','from' => '1'},'content' => [{'attrib' => {'min' => '1','max' => '5','type' => 'strings','name' => 'StringSatellite'},'content' => [{'attrib' => {'input' => 'select','required' => 'T','type' => 'string','name' => 'TypeSatelliteVar'},'content' => [{'attrib' => {'value' => 'MHD','default' => 'T','name' => 'MHD vars. dimensional'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'FUL','name' => 'All vars. dimensional'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'VAR','name' => 'Set vars. dimensional'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'mhd','name' => 'MHD vars. scaled'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'ful','name' => 'All vars. scaled'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'var','name' => 'Set vars. scaled'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'part'},{'attrib' => {'input' => 'select','required' => 'F','type' => 'string','name' => 'TypeTrajectory'},'content' => [{'attrib' => {'default' => 'T','name' => 'file'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'eqn','name' => 'equation'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'part'},{'attrib' => {'multiple' => 'T','input' => 'select','required' => 'F','type' => 'string','name' => 'TypeTime'},'content' => [{'attrib' => {'exclusive' => 'T','name' => 'none'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'step'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'date'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'time'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'part'}],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '-1','default' => '1','type' => 'integer','name' => 'DnOutput'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '-1','default' => '-1','type' => 'real','name' => 'DtOutput'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'length' => '100','type' => 'string','name' => 'NameTrajectoryFile'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'expr' => '-f $NameTrajectoryFile'},'content' => [{'content' => '
			Trajectory file $NameTrajectoryFile must exist
		','type' => 't'}],'type' => 'e','name' => 'rule'},{'attrib' => {'length' => '100','if' => '$TypeSatelliteVar =~ /\\bvar\\b/i','type' => 'string','name' => 'NameSatelliteVars'},'content' => [],'type' => 'e','name' => 'parameter'}],'type' => 'e','name' => 'for'},{'content' => '
#SATELLITE
2                       nSatellite
MHD file                StringSatellite (variables and traj type)
100                     DnOutput
-1.                     DtOutput [sec]
satellite1.dat          NameTrajectoryFile
VAR eqn step date       StringSatellite
100                     DnOutput
-1.                     DtOutput [sec]
satellite2.dat          NameTrajectoryFile
rho p                   NameSatelliteVars ! Read if StringSatellite 
                                          ! contains \'var\' or \'VAR\'

! The numerical solution can be extracted along one or more satellite
! trajectories. The number of satellites is defined by the 
! nSatellite parameter (default is 0).
!
! For each satellite the StringSatellite parameter determines what
! is saved into the satellite file(s).
! The StringSatellite can contain the following 3 parts in arbitrary order
!\\begin{verbatim}
! satellitevar   = \'mhd\', \'ful\' or \'var\' (unitless output)
!                  \'MHD\', \'FUL\' or \'VAR\' (dimensional output)
! trajectorytype = \'file\' or \'eqn\'
! timetype       = \'none\', \'step\', \'time\', \'date\'
!\\end{verbatim}
! The \'satellitevar\' part is required, 
! the \'trajectorytype\' part is optional (defaults to \'file\'), and
! the \'timetype\' part is also optional (default depends on satellitevar)
!
! The \'satellitevar\' string defines the variables to print in the satellite
! output file.  It also controls whether or not the variables will come out in
! dimensional or non-dimensional form by the capitalization of the
! satellitevars string: ALL CAPS means dimensional, all lower means 
! dimensionless. 
!
! If \'satellitevar\' is set to \'mhd\', the variables 
! \'rho ux uy uz bx by bz p jx jy jz\' will be saved, while \'ful\' implies
! \'rho ux uy uz bx by bz b1x b1y b1z p jx jy jz\'.
!
! If satellitevar is set to \'var\' then the list of variables is read 
! from the NameSatelliteVar parameter as a space separated list. 
! The choices for variables are currently:
!\\begin{verbatim}
! rho, rho, rhouy, rhouz, ux, uy, uz,Bx, By, Bz, B1x, B1y, B1z,
! E, P, Jx, Jy, Jz, theta1, theta2, phi1, phi2, status.
!\\end{verbatim}
!
! If \'trajectorytype\' is \'file\' (default) than the trajectory of the 
! satellite is read from the file given by the NameTrajectoryFile parameter.
! If \'trajectorytype\' is \'eqn\' then the trajectory is defined by an
! equation, which is hard coded in subroutine satellite_trajectory_formula
! in satellites.f90.
!
! The \'timetype\' values mean the following:
!\\begin{verbatim}
!  none  = there will be no indication of time in the logfile 
!          (not even the number of steps),
!  step  = number of time steps (n_steps),
!  date  = time is given as an array of 7 integers:  year mo dy hr mn sc msc,
!  time  = time is given as a real number - elapsed time since the start of
!          the run.  Units are determined by satellitevar and unitUSER_t.
!\\end{verbatim}
!  More than one \'timetype\' can be listed. They can be put together in any
!  combination.\\\\
!
! \\noindent
! The DnOutput and DtOutput parameters determine the frequency of extracting
! values along the satellite trajectories. \\\\
!
! \\noindent
! The extracted satellite information is saved into the files named
!\\begin{verbatim}
! PLOTDIR/sat_TRAJECTORYNAME_nTIMESTEP.sat
!\\end{verbatim}
! where TIMESTEP is the number of time steps (e.g. 000925), 
! and TRAJECTORYNAME is the name of the trajectory file.\\\\
!
! \\noindent
! The default is nSatellite=0, i.e. no satellite data is saved.
!
! Satellite input files contain the trajectory of the satellite.  They should
! have to following format:
!\\begin{verbatim}
! #COOR
! GSM
!
! #START
!  2004  6  24   0   0  58   0  2.9  -3.1 - 3.7  
!  2004  6  24   0   1  58   0  2.8  -3.2 - 3.6  
!\\end{verbatim}
!
! The #COOR command is optional.  It indicates which coordinate system the data
! represents.  The default is GSM, but others are possible.

! The file containing the upstream conditions should include data in the 
! following order:
!\\begin{verbatim}
! yr mn dy hr min sec msec x y z
!\\end{verbatim}
! with the position variables in units of the body radii or the length scale
! normalization.
!
! The maximum number of lines of data allowed in the input file is 50,000.  
! However, this can be modified by changing the variable Max_Satellite_Npts 
! in the file GM/BATSRUS/ModIO.f90.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'SAVEPLOT'},'content' => [{'attrib' => {'min' => '0','max' => '100','default' => '0','type' => 'integer','name' => 'nPlotFile'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'to' => '$nPlotFile','from' => '1'},'content' => [{'attrib' => {'min' => '3','max' => '3','type' => 'strings','name' => 'StringPlot'},'content' => [{'attrib' => {'input' => 'select','required' => 'T','type' => 'string','name' => 'plotform'},'content' => [{'attrib' => {'value' => 'tec','name' => 'TECPLOT'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'idl','name' => 'IDL'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'part'},{'attrib' => {'input' => 'select','required' => 'T','type' => 'string','name' => 'plotarea'},'content' => [{'attrib' => {'value' => '3d/3d_','name' => '3D'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'x=0'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'y=0','default' => 'T'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'z=0'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'sph'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'los'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'lin'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'ptn','if' => '$plotform =~ /\\btec\\b/'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'dpl','if' => '$plotform =~ /\\btec\\b/'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'cut','if' => '$plotform =~ /\\bidl\\b/'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'part'},{'attrib' => {'input' => 'select','required' => 'T','type' => 'string','name' => 'plotvar'},'content' => [{'attrib' => {'value' => 'MHD','name' => 'MHD vars. dimensional'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'FUL','name' => 'All vars. dimensional'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'RAW','name' => 'Raw vars. dimensional'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'RAY','name' => 'Ray tracing vars. dim.'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'FLX','name' => 'Flux vars. dimensional'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'SOL','name' => 'Solar vars. dimensional'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'POS','if' => '$plotarea eq \'lin\'','name' => 'Position vars. dimensional'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'VAR','name' => 'Select dimensional vars.'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'mhd','name' => 'MHD vars. scaled'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'ful','name' => 'All vars. scaled'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'raw','name' => 'Raw vars. scaled'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'ray','name' => 'Ray tracing vars. scaled'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'flx','name' => 'Flux vars. scaled'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'sol','name' => 'Solar vars. scaled'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'pos','if' => '$plotarea eq \'lin\'','name' => 'Position vars. scaled'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'var','name' => 'Select scaled vars.'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'part'}],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '-1','type' => 'integer','name' => 'DnSavePlot'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '-1','type' => 'real','name' => 'DtSavePlot'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'expr' => '$plotarea =~ /\\bptn\\b/'},'content' => [{'attrib' => {'type' => 'real','name' => 'xMinCut'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '$xMinCut','type' => 'real','name' => 'xMaxCut'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'type' => 'real','name' => 'yMinCut'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '$yMinCut','type' => 'real','name' => 'yMaxCut'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'type' => 'real','name' => 'zMinCut'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '$zMinCut','type' => 'real','name' => 'zMaxCut'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'type' => 'real','name' => 'xPoint'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'type' => 'real','name' => 'yPoint'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'type' => 'real','name' => 'zPoint'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'type' => 'real','name' => 'xNormal'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'type' => 'real','name' => 'yNormal'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'type' => 'real','name' => 'zNormal'},'content' => [],'type' => 'e','name' => 'parameter'}],'type' => 'e','name' => 'if'},{'attrib' => {'expr' => '$plotarea =~ /\\bdpl\\b/'},'content' => [{'attrib' => {'type' => 'real','name' => 'xMinCut'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '$xMinCut','type' => 'real','name' => 'xMaxCut'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'type' => 'real','name' => 'yMinCut'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '$yMinCut','type' => 'real','name' => 'yMaxCut'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'type' => 'real','name' => 'zMinCut'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '$zMinCut','type' => 'real','name' => 'zMaxCut'},'content' => [],'type' => 'e','name' => 'parameter'}],'type' => 'e','name' => 'if'},{'attrib' => {'expr' => '$plotarea =~ /\\bcut\\b/'},'content' => [{'attrib' => {'type' => 'real','name' => 'xMinCut'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '$xMinCut','type' => 'real','name' => 'xMaxCut'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'type' => 'real','name' => 'yMinCut'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '$yMinCut','type' => 'real','name' => 'yMaxCut'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'type' => 'real','name' => 'zMinCut'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '$zMinCut','type' => 'real','name' => 'zMaxCut'},'content' => [],'type' => 'e','name' => 'parameter'}],'type' => 'e','name' => 'if'},{'attrib' => {'min' => '0','if' => '$plotarea =~ /\\bsph\\b/','default' => '10','type' => 'real','name' => 'Radius'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'expr' => '$plotarea =~ /\\blos\\b/'},'content' => [{'attrib' => {'default' => '0','type' => 'real','name' => 'LosVectorX'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '0.0001','type' => 'real','name' => 'LosVectorY'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '1','type' => 'real','name' => 'LosVectorZ'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','default' => '20','type' => 'real','name' => 'xSizeImage'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','default' => '20','type' => 'real','name' => 'ySizeImage'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '10','type' => 'real','name' => 'xOffset'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '10','type' => 'real','name' => 'yOffset'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '1','default' => '2.5','type' => 'real','name' => 'rOccult'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','max' => '1','default' => '0.5','type' => 'real','name' => 'MuLimbDarkening'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '2','default' => '200','type' => 'integer','name' => 'nPixX'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '2','default' => '200','type' => 'integer','name' => 'nPixY'},'content' => [],'type' => 'e','name' => 'parameter'}],'type' => 'e','name' => 'if'},{'attrib' => {'expr' => '$plotarea =~ /\\blin\\b/'},'content' => [{'attrib' => {'input' => 'select','type' => 'string','name' => 'NameLine'},'content' => [{'attrib' => {'value' => 'A','name' => 'Advected B'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'default' => 'T','name' => 'B'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'U'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'J'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => 'F','type' => 'logical','name' => 'IsSingleLine'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '1','max' => '20','default' => '1','type' => 'integer','name' => 'nLine'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'to' => '$nLine','from' => '1'},'content' => [{'attrib' => {'type' => 'real','name' => 'xStartLine'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'type' => 'real','name' => 'yStartLine'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'type' => 'real','name' => 'zStartLine'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'type' => 'logical','name' => 'IsParallel'},'content' => [],'type' => 'e','name' => 'parameter'}],'type' => 'e','name' => 'for'}],'type' => 'e','name' => 'if'},{'attrib' => {'min' => '-1.0','if' => '($plotform=~/\\bidl\\b/ and $plotarea!~/\\b(sph|los|lin)\\b/)','default' => '-1.0','type' => 'real','name' => 'DxSavePlot'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'expr' => '$plotvar =~ /\\bvar\\b/i'},'content' => [{'attrib' => {'length' => '100','type' => 'string','name' => 'NameVars'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'length' => '100','type' => 'string','name' => 'NamePars'},'content' => [],'type' => 'e','name' => 'parameter'}],'type' => 'e','name' => 'if'}],'type' => 'e','name' => 'for'},{'content' => '
#SAVEPLOT
8			nPlotfile
3d  MHD tec		StringPlot ! 3d MHD data
100			DnSavePlot
-1.			DtSavePlot
y=0 VAR idl		StringPlot ! y=0 cut
-1			DnSavePlot
100.			DtSavePlot
2.			DxSavePlot ! Read only for format \'idl\'
jx jy jz		NameVars   ! Read only for content \'var\'
g unitx unitv unitn	NamePars   ! Read only for content \'var\'
cut ray idl		StringPlot ! ray tracing plot
1			DnSavePlot
-1.			DtSavePlot
-10.			xMinCut    ! Read only for area \'cut\'
10.			xMaxCut    ! Read only for area \'cut\'
-10.			yMinCut    ! Read only for area \'cut\'
10.			yMaxCut    ! Read only for area \'cut\'
-10.			zMinCut    ! Read only for area \'cut\'
10.			zMaxCut    ! Read only for area \'cut\'
1.			DxSavePlot ! Read only for format \'idl\'
sph flx idl		StringPlot ! spherical plot
-1			DnSavePlot
100.			DtSavePlot
4.			Radius     ! of spherical cut, Read only for area \'sph\'
los sol idl             StringPlot ! line of sight plot
-1			DnSavePlot
100.			DtSavePlot
1.			xLosVector
0.			yLosVector
0.			zLosVector
30.			xSizeImage
50.			ySizeImage
10.			xOffset
20.			yOffset
5.			rOccult
0.5			MuLimbDarkening
256			nPixX
256			nPixY
lin mhd idl		StringPlot  ! field line plot
-1			DnSavePlot
10.			DtSavePlot
B			NameLine ! B - magnetic field line, U - stream line
F			IsSingleLine
2			nLine
-2.0			xStartLine
0.0			yStartLine
3.5			zStartLine
F			IsParallel
-1.0			xStartLine
1.0			yStartLine
-3.5			zStartLine
T			IsParallel
dpl MHD tec		StringPlot  ! dipole tilt cut plot, ONLY WORKS WITH tec
-1			DnSavePlot
10.			DtSavePlot
-10.			xMinCut
 10.			xMaxCut
-10.			yMinCut
 10.			yMaxCut
-10.			zMinCut
 10.			zMaxCut
ptn MHD tec		StringPlot  ! point-normal cut plot, ONLY WORKS WITH tec
-1			DnSavePlot
10.			DtSavePlot
-10.			xMinCut
 10.			xMaxCut
-10.			yMinCut
 10.			yMaxCut
-10.			zMinCut
 10.			zMaxCut
 0.			xPoint
 0.			yPoint
 0.			zPoint
  0.			xNormal
  0.			yNormal
  1.			zNormal

! Default is nPlotFile=0. \\\\
!
! \\noindent
! StringPlot must contain the following 3 parts in arbitrary order
!\\begin{verbatim}
! plotarea plotvar plotform
!
! plotarea = \'3d\' , \'x=0\', \'y=0\', \'z=0\', \'cut\', \'dpl\', \'ptn\', \'sph\', \'los\', \'lin\'
! plotvar  = \'mhd\', \'ful\',\'raw\', \'ray\', \'flx\', \'sol\', \'var\' - unitless output
! plotvar  = \'MHD\', \'FUL\',\'RAW\', \'RAY\', \'FLX\', \'SOL\', \'VAR\' - dimensional
! plotform = \'tec\', \'idl\'
!\\end{verbatim}
! NOTES: The plotvar option \'sol\' is only valid for plotarea \'los\'.\\\\
!
!\\noindent
! The plotarea string defines the 1, 2, or 3D volume of the plotting area:
!\\begin{verbatim}
! x=0	- full x=0 plane: xmin=-0.001, xmax=0.001, average for symmetry plane
! y=0	- full y=0 plane: ymin=-0.001, ymax=0.001, average for symmetry plane
! z=0	- full z=0 plane: zmin=-0.001, zmax=0.001, average for symmetry plane
! 3d	- full 3D volume
! cut	- READ PLOTRANGE FROM PARAM.in, slightly different behavior for idl / Tec
! dpl   - cut at dipole \'equator\', uses PLOTRANGE to clip plot
! ptn   - cut plane given point and normal, uses PLOTRANGE to clip plot
! sph   - spherical cut at radius R_plot, READ FROM PARAM.in
! los   - line of sight integrated plot
! lin   - one dimensional plot along a field or stream line
!\\end{verbatim}
! The plotvar string defines the plot variables and the equation parameters.
! It also controls whether or not the variables will be plotted in dimensional
! values or as non-dimensional values:
!\\begin{verbatim}
! ALL CAPS  - dimensional
! all lower - dimensionless
!
! \'mhd\' - vars: rho Ux Uy Uz E Bx By Bz P Jx Jy Jz
!         pars: g eta
! \'ful\' - vars: rho Ux Uy Uz E Bx By Bz B1x B1y B1z P Jx Jy Jz
!         pars: g eta
! \'raw\' - vars: rho rhoUx rhoUy rhoUz E Bx By Bz P b1x b1y b1z divb
!         pars: g eta
! \'ray\' - vars: bx by bz theta1 phi1 theta2 phi2 status blk
!         pars: R_ray
! \'flx\' - vars: rho rhoUr Br jr pvecr
!         pars: g eta
! \'var\' - vars: READ FROM PARAMETER FILE
!         pars: READ FROM PARAMETER FILE
! \'sol\' - vars: wl pb
!         pars: mu
!\\end{verbatim}
! The plot_string is always followed by the plotting frequency
! DnSavePlot and for time accurate runs by DtSavePlot.\\\\
!
!\\noindent
! Depending on StringPlot, further information is read from the parameter file
! in this order:
!\\begin{verbatim}
! PlotRange		if plotarea is \'cut\', \'dpl\', or \'ptn\'
! Point,Normal		if plotarea is \'ptn\'
! DxSavePlot		if plotform is \'idl\' and plotarea is not sph, ion, los
! Radius		if plotarea is \'sph\'
! NameVars		if plotform is \'var\'
! NamePars		if plotform is \'var\'
!\\end{verbatim}
! The PlotRange is described by 6 coordinates. For IDL plots, If the width in 
! one or two 
! dimensions is less than the smallest cell size within the plotarea, 
! then the plot file will be 2 or 1 dimensional. If the range is thin but
! symmetric about one of the x=0, y=0, or z=0 planes, data will be averaged
! in the postprocessing.\\\\
!
! For Tecplot (tec) file and \'cut\', Plotrange is read but only 1 dimension is used.  
! Cuts are entire x, y, or z = constant planes (2D only, 1D or 3D cuts are not
! implemented.  For x=constant, for example, the y and z ranges 
! do not matter as long at they are "wider" than the x range.  The slice will be 
! located at the average of the two x ranges.  So, for example to save a plot in
! a x=-5 constant plane 
! cut in tec.  The following would work for the plot range:
!\\begin{verbatim}
! -5.01			xMinCut
! -4.99			xMaxCut
! -10.			yMinCut
!  10.			yMaxCut
! -10.			zMinCut
!  10.			zMaxCut
!\\end{verbatim}
!
! The \'dpl\' and \'ptn\' Tecplot plots use Plotrange like the IDL plots, and will
! clip the cut plane when it exits the defined box.
!
! \\noindent
! Point is described by the coordinate of any point on the cut plane, often the origin.
! Normal is the coordinate of a vector normal to the plane.  If the normal in any given
! coordinate direction is less than 0.01, then no cuts are computed for cell edges
! parallel to that coordinate direction.  For example, the following would result in
! only computing cuts on cell edges parallel to the Z axis.
!\\begin{verbatim}
!  0.0			xNormal
!  0.0			yNormal
!  1.0			zNormal
!\\end{verbatim}
!
! \\noindent
! Possible values for DxSavePlot (for IDL files):
!\\begin{verbatim}
!  0.5	- fixed resolution (any positive value)
!  0.	- fixed resolution based on the smallest cell in the plotting area
! -1.	- unstructured grid will be produced by PostIDL.exe
!\\end{verbatim}
! Radius is the radius of the spherical cut for plotarea=\'sph\'
!
! With #SAVEPLOT it is possible to create plots which are the integral along
! the line of site of some quantity.  The variables which control this are the 
! following:
!\\begin{verbatim}
!    LosVectorX,Y,Z define the direction of the line of sight integration
!    xSizeImage, ySizeImage defines the size of the LOS image
!    xOffset, yOffset defines the offset relative to the origin (Sun)
!    rOccult defines the minimum distance of the line from the origin (Sun)
!    MuLimbDarkening is the limb darkening parameter for the \'wl\' (white light)
!                 and \'pb\' (polarization brightness) plot variables.
!\\end{verbatim}
!
!\\noindent
! The possible values for NameVars with plotarea \'los\' 
!       are listed in subroutine set_plotvar_los in write_plot_los.f90. \\\\
!
! \\noindent
! The possible values for NameVars for other plot areas
!       are listed in subroutine set_plotvar in write_plot_common.f90.\\\\
!
! \\noindent
! The possible values for NamePars 
!       are listed in subroutine set_eqpar in write_plot_common.f90\\\\
!
! A plot file is produced by each processor.  This file is ASCII in \'tec\'
! format and can be either binary or ASCII in \'idl\' format as chosen under
! the #SAVEBINARY flag.  The name of the files are
!\\begin{verbatim}
! IO2/plotarea_plotvar_plotnumber_timestep_PEnumber.extenstion
!\\end{verbatim}
! where extension is \'tec\' for the TEC and \'idl\' for the IDL file formats.
! The plotnumber goes from 1 to nplot in the order of the files in PARAM.in.
! After all processors wrote their plot files, processor 0 writes a small 
! ASCII header file named as
!\\begin{verbatim}
! IO2/plotarea_plotvar_plotnumber_timestep.headextension
!\\end{verbatim}
! where headextension is:
!\\begin{verbatim}
!           \'T\' for TEC file format
!           \'S\' for TEC and plot_area \'sph\' 
!           \'h\' for IDL file format       
!\\end{verbatim}
!
!\\noindent
! The line of sight integration produces TecPlot and IDL files directly:
!\\begin{verbatim}
! IO2/los_plotvar_plotnumber_timestep.extension
\\end{verbatim}
! where extension is \'dat\' for TecPlot and \'out\' for IDL file formats.
! The IDL output from line of sight integration is always in ASCII format.

','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'SAVEBINARY'},'content' => [{'attrib' => {'default' => 'T','type' => 'logical','name' => 'DoSaveBinary'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '
#SAVEBINARY
T			DoSaveBinary   used only for \'idl\' plot file

! Default is .true. Saves unformatted IO2/*.idl files if true. 
! This is the recommended method, because it is fast and accurate.
! The only advantage of saving IO2/*.idl in formatted text files is
! that it can be processed on another machine or with a different 
! (lower) precision. For example PostIDL.exe may be compiled with 
! single precision to make IO2/*.out files smaller, while BATSRUS.exe is 
! compiled in double precision, to make results more accurate.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'SAVEPLOTSAMR'},'content' => [{'attrib' => {'default' => 'F','type' => 'logical','name' => 'DoSavePlotsAmr'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '
#SAVEPLOTSAMR
F			DoSavePlotsAmr

! Save plots before each AMR. Default is DoSavePlotsAMR=.false.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'FLUSH'},'content' => [{'attrib' => {'default' => 'T','type' => 'logical','name' => 'DoFlush'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

#FLUSH
F			DoFlush

If the DoFlush variable is true, the output is flushed when
subroutine ModUtility::flush_unit is called. This is used in the 
log and satellite files. The flush is useful to see the output immediately, 
and to avoid truncated files when the code crashes,
but on some systems the flush may be very slow. 

The default is to flush the output, i.e. DoFlush=T.
','type' => 't'}],'type' => 'e','name' => 'command'}],'type' => 'e','name' => 'commandgroup'},{'attrib' => {'name' => 'AMR PARAMETERS'},'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  AMR PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
','type' => 't'},{'attrib' => {'if' => '$_IsFirstSession','name' => 'AMRINIT'},'content' => [{'attrib' => {'input' => 'select','type' => 'string','name' => 'InitialRefineType'},'content' => [{'attrib' => {'default' => '1','name' => 'default'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'all'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'none'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => '3Dbodyfocus'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'spherefocus'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'magnetosphere'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'points'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'coupledhelio'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'helio_init'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'helio_z=4'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'all_then_focus'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'cme'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'points'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'mag_new'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'magnetosphere'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'magneto_fine'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'magneto12'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'magnetosaturn'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'magnetojupiter'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'paleo'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'comet'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','default' => '4','type' => 'integer','name' => 'InitialRefineLevel'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '
#AMRINIT
default			TypeRefineInit
4			nRefineLevelInit

! These are the default values for the initial refinement.\\\\
! Possible values for InitialRefineType:\\\\
! Default depends on problem_type. 
!\\begin{verbatim}
! \'none\'		- Refine no blocks
! \'all\' 		- Refine all blocks
! \'3Dbodyfocus\'		- Refinement focusing on body
! \'spherefocus\'		- Refinement focusing on the origin, does not require 
!                           a body
! \'points\'      	- Refine around given points
! \'magnetosphere\'	- Refine for generic magnetosphere
! *			- any other value will use default value by ProblemType
!\\end{verbatim}
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'if' => '$_IsFirstSession','name' => 'AMRINITPHYSICS'},'content' => [{'attrib' => {'min' => '0','default' => '0','type' => 'integer','name' => 'nRefineLevelIC'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '
#AMRINITPHYSICS
3			nRefineLevelIC

! Defines number of physics (initial condition) based AMR-s AFTER the 
! geometry based initial AMR-s defined by #AMRINIT were done.
! Only useful if the initial condition has a non-trivial analytic form.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'multiple' => 'T','alias' => 'GRIDLEVEL','name' => 'GRIDRESOLUTION'},'content' => [{'attrib' => {'value' => '0','name' => 'RotateArea'},'content' => [],'type' => 'e','name' => 'set'},{'attrib' => {'min' => '0','if' => '$_command =~ /RESOLUTION/','type' => 'real','name' => 'Resolution'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','if' => '$_command =~ /LEVEL/','type' => 'integer','name' => 'nLevel'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '1','max' => '2','type' => 'strings','name' => 'StringArea'},'content' => [{'attrib' => {'input' => 'select','required' => 'T','type' => 'string','name' => 'NameArea'},'content' => [{'attrib' => {'value' => 'init/initial','default' => 'T','name' => 'initial'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'all'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'box'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'brick'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'brick0'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'sphere'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'sphere0'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'shell'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'shell0'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'cylinderx'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'cylinderx0'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'cylindery'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'cylindery0'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'cylinderz'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'cylinderz0'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'ringx'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'ringx0'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'ringy'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'ringy0'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'ringz'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'ringz0'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'part'},{'attrib' => {'input' => 'select','required' => 'F','type' => 'string','name' => 'RotateArea'},'content' => [{'attrib' => {'name' => 'rotated'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'part'}],'type' => 'e','name' => 'parameter'},{'attrib' => {'expr' => '$NameArea !~ /box|all|init|0/'},'content' => [{'attrib' => {'default' => '0','type' => 'real','name' => 'xCenter'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '0','type' => 'real','name' => 'yCenter'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '0','type' => 'real','name' => 'zCenter'},'content' => [],'type' => 'e','name' => 'parameter'}],'type' => 'e','name' => 'if'},{'attrib' => {'expr' => '$NameArea =~ /box/'},'content' => [{'attrib' => {'type' => 'real','name' => 'xMinBox'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'type' => 'real','name' => 'yMinBox'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'type' => 'real','name' => 'zMinBox'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'type' => 'real','name' => 'xMaxBox'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'type' => 'real','name' => 'yMaxBox'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'type' => 'real','name' => 'zMaxBox'},'content' => [],'type' => 'e','name' => 'parameter'}],'type' => 'e','name' => 'if'},{'attrib' => {'expr' => '$NameArea =~ /brick/i'},'content' => [{'attrib' => {'min' => '0','type' => 'real','name' => 'xSizeBrick'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','type' => 'real','name' => 'ySizeBrick'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','type' => 'real','name' => 'zSizeBrick'},'content' => [],'type' => 'e','name' => 'parameter'}],'type' => 'e','name' => 'if'},{'attrib' => {'expr' => '$NameArea =~ /sphere/i'},'content' => [{'attrib' => {'min' => '0','type' => 'real','name' => 'Radius'},'content' => [],'type' => 'e','name' => 'parameter'}],'type' => 'e','name' => 'if'},{'attrib' => {'expr' => '$NameArea =~ /shell/i'},'content' => [{'attrib' => {'min' => '0','type' => 'real','name' => 'Radius1'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','type' => 'real','name' => 'Radius2'},'content' => [],'type' => 'e','name' => 'parameter'}],'type' => 'e','name' => 'if'},{'attrib' => {'expr' => '$NameArea =~ /cylinder/i'},'content' => [{'attrib' => {'min' => '0','type' => 'real','name' => 'LengthCylinder'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','type' => 'real','name' => 'Radius'},'content' => [],'type' => 'e','name' => 'parameter'}],'type' => 'e','name' => 'if'},{'attrib' => {'expr' => '$NameArea =~ /ring/i'},'content' => [{'attrib' => {'min' => '0','type' => 'real','name' => 'HeightRing'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','type' => 'real','name' => 'Radius1'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','type' => 'real','name' => 'Radius2'},'content' => [],'type' => 'e','name' => 'parameter'}],'type' => 'e','name' => 'if'},{'attrib' => {'expr' => '$RotateArea =~ /rotated/i'},'content' => [{'attrib' => {'min' => '-360','max' => '360','default' => '0','type' => 'real','name' => 'xRotate'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '-360','max' => '360','default' => '0','type' => 'real','name' => 'yRotate'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '-360','max' => '360','default' => '0','type' => 'real','name' => 'zRotate'},'content' => [],'type' => 'e','name' => 'parameter'}],'type' => 'e','name' => 'if'},{'content' => '

#GRIDRESOLUTION
2.0			Resolution
initial			NameArea

#GRIDLEVEL
3			nLevel
all			NameArea

#GRIDLEVEL
4			nLevel
box			NameArea
-64.0			xMinBox
-16.0			yMinBox
-16.0			zMinBox
-32.0			xMaxBox
 16.0			yMaxBox
  0.0			zMaxBox

#GRIDLEVEL
4			nLevel
brick			NameArea
-48.0			xCenter
  0.0			yCenter
 -8.0			zCenter
 32.0			xSizeBrick
 32.0			ySizeBrick
 16.0			zSizeBrick

#GRIDRESOLUTION
1/8			Resolution
shell0			NameArea
3.5			Radius1
4.5			Radius2

#GRIDRESOLUTION
0.5			Resolution
sphere			NameArea
-10.0			xCenterSphere
 10.0			yCenterSphere
  0.0			zCenterSphere
 20.0			rSphere

#GRIDRESOLUTION
1/8			Resolution
cylinderx		NameArea
-30.0			xCenter
  0.0			yCenter
  0.0			zCenter
 60.0			LengthCylinder
 20.0			rCylinder

#GRIDRESOLUTION
1/8			Resolution
ringz0 rotated		NameArea
  5.0			HeightRing
 20.0			Radius1
 25.0			Radius2
 10.0			xRotate
 10.0			yRotate
  0.0			zRotate

The #GRIDRESOLUTION and #GRIDLEVEL commands allow to set the grid resolution
or refinement level, respectively, in a given area. The Resolution parameter
refers to the size of the cell in the X direction (Dx).
The nLevel parameter is an integer with level 0 meaning no refinement relative
to the root block, while level N is a refinement by 2 to the power N.

The next parameter NameArea defines the shape of the area. 
If NameArea is set to \'initial\' or \'init\', the highest initial resolution
or level is set by the command. This is similar to the #AMRINIT command.
For other values of NameArea, the command specifies the blocks to be refined.
All computational blocks that intersect the area and have a coarser
resolution than the resolution set for the area are selected for refinement.
There are the following basic shapes: 
\'all\', \'box\', \'brick\', \'sphere\', \'shell\', \'cylinderx\', \'cylindery\', 
\'cylinderz\', \'ringx\', \'ringy\' and \'ringz\'.

The area \'all\' refers to the whole computational domain, and it can be
used to set the overall minimum resolution. The area \'box\' is a box
aligned with the X, Y and Z axes, and it is given with the coordinates
of two diagonally opposite corners. The area \'brick\' has the same shape
as \'box\', but it is defined with the center of the brick and the 
size of the brick. The area \'sphere\' is a sphere around an arbitrary point,
which is defined with the center point and the radius of the sphere.
The area \'shell\' consists of the volume between two concentric spherical
surfaces, which is given with the center point and the two radii.
The area \'cylinderx\' is a cylinder with an axis parallel with the X axis,
and it is given with the center, the length of the axis and the radius,
The areas \'cylindery\' and \'cylinderz\' are cylinders parallel with the 
Y and Z axes, respectively, and are defined analogously as \'cylinderx\'.
The area \'ringx\', \'ringy\' and \'ringz\' are the volumes between 
two cylindrical surfaces parallel with the X, Y and Z axes, respectively.
The ring area is given with the center, the height and the two radii.

If the area name contains a zero, the center is taken to be at the origin.
If the word rotated is added, the area can be rotated by 3 angles around
the X, Y and Z axes in this order.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'AMRLEVELS'},'content' => [{'attrib' => {'min' => '-1','default' => '0','type' => 'integer','name' => 'MinBlockLevel'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '-1','default' => '99','type' => 'integer','name' => 'MaxBlockLevel'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => 'F','type' => 'logical','name' => 'DoFixBodyLevel'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '
#AMRLEVELS
0			MinBlockLevel
99			MaxBlockLevel
F			DoFixBodyLevel

! Set the minimum/maximum levels that can be affected by AMR.  The usage is as
! follows:
!\\begin{verbatim}
! MinBlockLevel .ge.0 Cells can be coarsened up to the listed level but not
!                       further.
! MinBlockLevel .lt.0 The current grid is ``frozen\'\' for coarsening such that
!                       blocks are not allowed to be coarsened to a size
!                       larger than their current one.
! MaxBlockLevel .ge.0 Any cell at a level greater than or equal to
!                       MaxBlockLevel is unaffected by AMR (cannot be coarsened
!                       or refined).
! MaxBlockLevel .lt.0 The current grid is ``frozen\'\' for refinement such that
!                       blocks are not allowed to be refined to a size
!                       smaller than their current one.
! DoFixBodyLevel = T  Blocks touching the body cannot be coarsened or refined.
!\\end{verbatim}
! This command has no effect when DoAutoRefine is .false. in the #AMR command.
!
! Note that the user can set either #AMRLEVELS or #AMRRESOLUTION but not
! both.  If both are set, the final one in the session will set the values
! for AMR.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'AMRRESOLUTION'},'content' => [{'attrib' => {'min' => '-1','default' => '0','type' => 'real','name' => 'DxCellMin'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '-1','default' => '99999','type' => 'real','name' => 'DxCellMax'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => 'F','type' => 'logical','name' => 'DoFixBodyLevel'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '
#AMRRESOLUTION
0.			DxCellMin
99999.			DxCellMax
F			DoFixBodyLevel

! Serves the same function as AMRLEVELS. The DxCellMin and DxCellMmax
! parameters are converted into MinBlockLevel and MaxBlockLevel 
! when they are read.
! Note that MinBlockLevel corresponds to DxCellMax and MaxBlockLevel
! corresponds to DxCellMin.  See details above.
!
! This command has no effect when DoAutoRefine is .false. in the #AMR command.
!
! Note that the user can set either #AMRLEVELS or #AMRRESOLUTION but not
! both.  If both are set, the final one in the session will set the values
! for AMR.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'AMR'},'content' => [{'attrib' => {'min' => '-1','default' => '-1','type' => 'integer','name' => 'DnRefine'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'expr' => '$DnRefine>0'},'content' => [{'attrib' => {'default' => 'F','type' => 'logical','name' => 'DoAutoRefine'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'expr' => '$DoAutoRefine'},'content' => [{'attrib' => {'min' => '0','max' => '100','default' => '20','type' => 'real','name' => 'PercentCoarsen'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','max' => '100','default' => '20','type' => 'real','name' => 'PercentRefine'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '1','default' => '99999','type' => 'integer','name' => 'MaxTotalBlocks'},'content' => [],'type' => 'e','name' => 'parameter'}],'type' => 'e','name' => 'if'}],'type' => 'e','name' => 'if'},{'content' => '
#AMR
2001			DnRefine
T			DoAutoRefine   ! read if DnRefine is positive
0.			PercentCoarsen ! read if DoAutoRefine is true
0.			PercentRefine  ! read if DoAutoRefine is true
99999			MaxTotalBlocks ! read if DoAutoRefine is true

! The DnRefine parameter determines the frequency of adaptive mesh refinements
! in terms of total steps nStep.
!
! When DoAutoRefine is false, the grid is refined by one more level
! based on the TypeRefineInit parameter given in the #AMRINIT command. 
! If the number of blocks is not sufficient for this pre-specified refinement, 
! the code stops with an error.
!
! When DoAutoRefine is true, the grid is refined or coarsened 
! based on the criteria given in the #AMRCRITERIA command.
! The number of blocks to be refined or coarsened are determined by
! the PercentRefine and PercentCoarsen parameters. These percentages
! are approximate only, because the constraints of the block adaptive
! grid may result in more or fewer blocks than prescribed.
! The total number of blocks will not exceed the smaller of the 
! MaxTotalBlocks parameter and the total number of blocks available on all 
! the PE-s (which is determined by the number of PE-s and 
! the MaxBlocks parameter in ModSize.f90).
! 
! Default for DnRefine is -1, i.e. no run time refinement.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'AMRCRITERIA'},'content' => [{'attrib' => {'input' => 'select','type' => 'integer','name' => 'nRefineCrit'},'content' => [{'attrib' => {'name' => '1'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => '2'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'default' => '1','name' => '3'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'parameter'},{'attrib' => {'to' => '$nRefineCrit','from' => '1'},'content' => [{'attrib' => {'input' => 'select','type' => 'string','name' => 'TypeRefine'},'content' => [{'attrib' => {'value' => 'gradt/gradT','name' => 'grad T'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'gradp/gradP','name' => 'grad P'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'gradlogrho','name' => 'grad log(Rho)'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'gradlogP/gradlogp','name' => 'grad log(p)'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'gradE','name' => 'grad E'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'curlV/curlv/curlU/curlu','name' => 'curl U'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'curlB/curlb','name' => 'curl B'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'divU/divu/divV/divv','name' => 'div U'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'divb/divB','name' => 'divB'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'Valfven/vAlfven/valfven','name' => 'vAlfven'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'heliobeta','name' => 'heliospheric beta'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'flux'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'heliocurrentsheet','name' => 'heliospheric current sheet'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'rcurrents/Rcurrents','name' => 'rCurrents'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'transient/Transient','name' => 'Transient'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'parameter'},{'attrib' => {'expr' => '$TypeRefine =~ /transient/i'},'content' => [{'attrib' => {'input' => 'select','type' => 'string','name' => 'TypeTransient'},'content' => [{'attrib' => {'value' => 'p_dot/P_dot','name' => 'P_dot'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 't_dot/T_dot','name' => 'T_dot'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'rho_dot/Rho_dot','default' => 'T','name' => 'Rho_dot'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'RhoU_dot/rhou_dot','name' => 'RhoU_dot'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'Rho_2nd_1/rho_2nd_1','name' => 'Rho_2nd_1'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'Rho_2nd_2/rho_2nd_2','name' => 'Rho_2nd_2'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => 'F','type' => 'logical','name' => 'UseSunEarth'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'expr' => '$UseSunEarth'},'content' => [{'attrib' => {'type' => 'real','name' => 'xEarth'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'type' => 'real','name' => 'yEarth'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'type' => 'real','name' => 'zEarth'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'type' => 'real','name' => 'InvD2Ray'},'content' => [],'type' => 'e','name' => 'parameter'}],'type' => 'e','name' => 'if'}],'type' => 'e','name' => 'if'}],'type' => 'e','name' => 'for'},{'content' => '
#AMRCRITERIA
3			nRefineCrit (number of refinement criteria: 1,2 or 3)
gradlogP		TypeRefine
divB			TypeRefine
Transient		TypeRefine
Rho_dot			TypeTransient ! Only if \'Transient\' or \'transient\'
T			UseSunEarth   ! Only if \'Transient\'
0.00E+00		xEarth        ! Only if UseSunEarth
2.56E+02 		yEarth        ! Only if UseSunEarth
0.00E+00		zEarth        ! Only if UseSunEarth
5.00E-01		InvD2Ray      ! Only if UseSunEarth

! The default values depend on problem_type.\\\\ 
! At most three criteria can be given. Possible criteria:
!\\begin{verbatim}
! \'gradT\'		- gradient of temperature
! \'gradP\'		- gradient of pressure
! \'gradlogrho\'		- gradient of log(rho)
! \'gradlogP\'		- gradient of log(P)
! \'gradE\'		- gradient of electric field magnitude
! \'curlV\',\'curlU\' 	- magnitude of curl of velocity
! \'curlB\'		- magnitude of current
! \'divU\', \'divV\'	- divergence of velocity
! \'divB\'		- div B
! \'vAlfven\',\'Valfven\'	- Alfven speed
! \'heliobeta\' 		- special function for heliosphere $R^2 B^2/rho$
! \'flux\'		- radial mass flux
! \'heliocurrentsheet\'	- refinement in the currentsheet of the heliosphere
! \'Rcurrents\'		- refinement near Rcurrents value
!\\end{verbatim}
! All the names can also be spelled with all small case letters.\\\\
!
!\\noindent
! The possible choices for TypeTransient:
!\\begin{verbatim}
! \'P_dot\' (same as \'p_dot\')
! \'T_dot\' (same as \'t_dot\')
! \'Rho_dot\' (same as \'rho_dot\')
! \'RhoU_dot\' (same as \'rhou_dot\')
! \'B_dot\' (same as \'b_dot\')
! \'Rho_2nd_1\' (same as \'rho_2nd_1\')
! \'Rho_2nd_2\' (same as \'rho_2nd_2\')
!\\end{verbatim}
!
! Also, (xEarth,yEarth,zEarth) are the coordinates of the Earth. InvD2Ray is
! a factor that defines how close to the ray Sun-Earth to refine the grid.
! Note that the AMR occurs in a cylinder around the ray.
! Example for InvD2Ray =
\\begin{verbatim}
!   1 - refine_profile = 0.3679 at distance Rsun/10 from the ray
!   2 - refine_profile = 0.0183 at distance Rsun/10 from the ray
!   3 - refine_profile = 0.0001 at distance Rsun/10 from the ray
\\end{verbatim}
','type' => 't'}],'type' => 'e','name' => 'command'}],'type' => 'e','name' => 'commandgroup'},{'attrib' => {'name' => 'SCHEME PARAMETERS'},'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  SCHEME PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
','type' => 't'},{'attrib' => {'name' => 'SCHEME'},'content' => [{'attrib' => {'input' => 'select','type' => 'integer','name' => 'nOrder'},'content' => [{'attrib' => {'default' => 'T','name' => '1'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => '2'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'parameter'},{'attrib' => {'input' => 'select','type' => 'string','name' => 'TypeFlux'},'content' => [{'attrib' => {'value' => 'Roe/roe/1','name' => 'Roe'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'Rusanov/rusanov/2/TVDLF','default' => 'T','name' => 'Rusanov'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'Linde/linde/3/HLLEL','name' => 'Linde'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'Sokolov/sokolov/4/AW','name' => 'Sokolov'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'parameter'},{'attrib' => {'expr' => '$nOrder == 2'},'content' => [{'attrib' => {'input' => 'select','type' => 'string','name' => 'TypeLimiter'},'content' => [{'attrib' => {'default' => 'T','name' => 'minmod'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'mc'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'beta'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '1','max' => '2','if' => '$TypeLimiter ne \'minmod\'','default' => '1.2','type' => 'real','name' => 'LimiterBeta'},'content' => [],'type' => 'e','name' => 'parameter'}],'type' => 'e','name' => 'if'},{'content' => '
#SCHEME
2			nOrder (1 or 2)
Rusanov			TypeFlux
mc			TypeLimiter ! Only for nOrder=2
1.2			LimiterBeta ! Only if LimiterType is NOT \'minmod\'

! Default values are shown above.\\\\
!
!\\noindent
! Possible values for TypeFlux:
!\\begin{verbatim}
! \'Rusanov\'     - Rusanov or Lax-Friedrichs flux     
! \'Linde        - Linde\'s HLLEL flux                   
! \'Sokolov\'     - Sokolov\'s Local Artificial Wind flux 
! \'Roe\'         - Roe\'s approximate Riemann flux       
!\\end{verbatim}
! Possible values for TypeLimiter:
!\\begin{verbatim}
! \'minmod\'	- minmod limiter is the most robust limiter
! \'mc\'          - monotonized central limiter with a beta parameter
! \'beta\'        - beta limiter is less robust than the mc limiter for 
!                 the same beta value
!\\end{verbatim}
! Possible values for LimiterBeta (for \'beta\' and \'mc\' limiters only)
! are between 1.0 and 2.0 : 
!\\begin{verbatim}
!  LimiterBeta = 1.0 is the same as the minmod limiter
!  LimiterBeta = 2.0 for the beta limiter is the same as the superbee limiter
!  LimiterBeta = 1.5 is a typical value for the mc limiter
!  LimiterBeta = 1.2 is the recommended value for the beta limiter
!\\end{verbatim}
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'NONCONSERVATIVE'},'content' => [{'attrib' => {'default' => 'T','type' => 'logical','name' => 'UseNonConservative'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '
#NONCONSERVATIVE
T		UseNonConservative

! For Earth the default is using non-conservative equations 
! (close to the body).
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'CONSERVATIVECRITERIA'},'content' => [{'attrib' => {'min' => '0','max' => '3','default' => '1','type' => 'integer','name' => 'nConservCrit'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'to' => '$nConservCrit','from' => '1'},'content' => [{'attrib' => {'input' => 'select','type' => 'string','name' => 'TypeConservCrit'},'content' => [{'attrib' => {'value' => 'r/R/radius/Radius','default' => 'T','name' => 'radius'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'parabola/paraboloid','name' => 'parabola'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'p/P','name' => 'p'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'gradp/GradP','name' => 'grad P'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '$rBody','if' => '$TypeConservCrit =~ /^r|radius$/i','default' => '6','type' => 'real','name' => 'rConserv'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','if' => '$TypeConservCrit =~ /^parabol/i','default' => '6','type' => 'real','name' => 'xParabolaConserv'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','if' => '$TypeConservCrit =~ /^parabol/i','default' => '36','type' => 'real','name' => 'yParabolaConserv'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','if' => '$TypeConservCrit =~ /^p$/i','default' => '0.05','type' => 'real','name' => 'pCoeffConserv'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','if' => '$TypeConservCrit =~ /gradp/i','default' => '0.1','type' => 'real','name' => 'GradPCoeffConserv'},'content' => [],'type' => 'e','name' => 'parameter'}],'type' => 'e','name' => 'for'},{'content' => '

#CONSERVATIVECRITERIA
3		nConservCrit
r		TypeConservCrit
6.		rConserv             ! read if TypeConservCrit is \'r\'
parabola        TypeConservCrit
6.		xParabolaConserv     ! read if TypeConservCrit is \'parabola\'
36.		yParabolaConserv     ! read if TypeConservCrit is \'parabola\'
p		TypeConservCrit
0.05		pCoeffConserv	     ! read if TypeConservCrit is \'p\'
GradP		TypeConservCrit
0.1		GradPCoeffConserv    ! read if TypeConservCrit is \'GradP\'

! Select the parts of the grid where the conservative vs. non-conservative
! schemes are applied. The number of criteria is arbitrary, although 
! there is no point applying the same criterion more than once.
!
! If no criteria is used, the whole domain will use conservative or
! non-conservative equations depending on UseNonConservative set in
! command #NONCONSERVATIVE.
!
! The physics based conservative criteria (\'p\' and \'GradP\')
! select cells which use the non-conservative scheme if ALL of them are true:
!\\begin{verbatim}
! \'p\'      - the pressure is smaller than fraction pCoeffConserv of the energy
! \'GradP\'  - the relative gradient of pressure is less than GradPCoeffConserv
!\\end{verbatim}
! The geometry based criteria are applied after the physics based criteria 
! (if any) and they select the non-conservative scheme if ANY of them is true:
!\\begin{verbatim}
! \'r\'        - radial distance of the cell is less than rConserv
! \'parabola\' - x less than xParabolaConserv - (y**2+z**2)/yParabolaConserv
!\\end{verbatim}
! Default values are nConservCrit = 1 with TypeConservCrit = \'r\'
! and rConserv=2*rBody, where rBody has a problem dependent default.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'UPDATECHECK'},'content' => [{'attrib' => {'default' => 'T','type' => 'logical','name' => 'UseUpdateCheck'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'expr' => '$UseUpdateCheck'},'content' => [{'attrib' => {'min' => '0','max' => '100','default' => '40','type' => 'real','name' => 'RhoMinPercent'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '100','default' => '400','type' => 'real','name' => 'RhoMaxPercent'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','max' => '100','default' => '40','type' => 'real','name' => 'pMinPercent'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '100','default' => '400','type' => 'real','name' => 'pMaxPercent'},'content' => [],'type' => 'e','name' => 'parameter'}],'type' => 'e','name' => 'if'},{'content' => '
#UPDATECHECK
T			UseUpdateCheck
40.			RhoMinPercent
400.			RhoMaxPercent
40.			pMinPercent
400.			pMaxPercent

! Default values are shown.  This will adjust the timestep so that
! density and pressure cannot change by more than the given percentages
! in a single timestep.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'PROLONGATION'},'content' => [{'attrib' => {'input' => 'select','type' => 'integer','name' => 'nOrderProlong'},'content' => [{'attrib' => {'default' => 'T','name' => '1'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => '2'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'parameter'},{'attrib' => {'input' => 'select','if' => '$nOrderProlong==2','type' => 'string','name' => 'TypeProlong'},'content' => [{'attrib' => {'value' => 'lr','default' => 'T','name' => 'left-right'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'central'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'name' => 'minmod'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'lr2','name' => 'left-right extrapolate'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'central2','name' => 'central    extrapolate'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'minmod2','name' => 'minmod     extrapolate'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'parameter'},{'content' => '
#PROLONGATION
2			nOrderProlong (1 or 2 for ghost cells)
lr			TypeProlong  ! Only for nOrderProlong=2

! Default is prolong_order=1. \\\\
! Possible values for prolong_type:\\\\
!
!\\noindent
! 1. in message_pass_dir (used if limiter_type is not \'LSG\')
!\\begin{verbatim}
!    \'lr\'          - interpolate only with left and right slopes 
!    \'central\'	   - interpolate only with central difference slope
!    \'minmod\' 	   - interpolate only with minmod limited slope
!    \'lr2\'	   -  like \'lr\' but extrapolate when necessary
!    \'central2\'	   - like \'central\' but extrapolate when necessary
!    \'minmod2\'	   - like \'minmod\' but extrapolate when necessary
!    \'lr3\'         - only experimental
!\\end{verbatim}
!
!\\noindent
! 2. in messagepass_all (used if limiter_type is \'LSG\')
\\begin{verbatim}
!    \'lr\',\'lr2\'		    - left and right slopes (all interpolation)
!    \'central\',\'central2\'   - central differences (all interpolation)
!    \'minmod\',\'minmod2\'	    - to be implemented
!\\end{verbatim}
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'alias' => 'OPTIMIZE','name' => 'MESSAGEPASS'},'content' => [{'attrib' => {'input' => 'select','type' => 'string','name' => 'TypeMessagePass'},'content' => [{'attrib' => {'value' => 'allopt','default' => 'T','name' => 'm_p_cell FACES ONLY'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'all','name' => 'm_p_cell'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'opt','name' => 'm_p_dir FACES ONLY'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'dir','name' => 'm_p_dir group by directions'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'face','name' => 'm_p_dir group by faces     '},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'min','name' => 'm_p_dir group by kind and face'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'parameter'},{'content' => '
#MESSAGEPASS
allopt			TypeMessagePass

! Default value is shown above.\\\\
! Possible values for optimize_message_pass:
!\\begin{verbatim}
! \'dir\'		- message_pass_dir: group messages direction by direction
! \'face\'	- message_pass_dir: group messages face by face
! \'min\'		- message_pass_dir: send equal, restricted and prolonged 
!				    messages face by face
!
! \'opt\'		- message_pass_dir: do not send corners, send one layer for
!				    first order, send direction by direction
!
! \'all\'		- message_pass_cell: corners, edges and faces in single message
!
! \'allopt\'      - message_pass_cell:  faces only in a single message
!\\end{verbatim}
! Constrained transport requires corners, default is \'all\'!\\\\ 
! Diffusive control requires corners, default is \'all\'!\\\\
! Projection uses message_pass_dir for efficiency!\\\\
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'TVDRESCHANGE'},'content' => [{'attrib' => {'default' => 'F','type' => 'logical','name' => 'UseTvdAtReschange'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '
#TVDRESCHANGE
T		UseTvdAtResChange

! For UseTvdAtResChange=T a second order TVD limited scheme is used 
! at the resolution changes. 
!
! Default value is false, which results in first order 
! prolongation and restriction operators at the resolution changes.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'BORIS'},'content' => [{'attrib' => {'default' => 'F','type' => 'logical','name' => 'UseBorisCorrection'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','max' => '1','if' => '$UseBorisCorrection','default' => '1','type' => 'real','name' => 'BorisClightFactor'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '
#BORIS
T			UseBorisCorrection
1.0			BorisClightFactor !Only if UseBorisCorrection is true

! Default is boris_correction=.false.
! Use semi-relativistic MHD equations with speed of light reduced by
! the BorisClightFactor. Set BorisClightFactor=1.0 for true semi-relativistic
! MHD. Gives the same steady state as normal MHD analytically, but there
! can be differences due to discretization errors. 
! You can use either Boris or BorisSimple but not both. 
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'BORISSIMPLE'},'content' => [{'attrib' => {'default' => 'F','type' => 'logical','name' => 'UseBorisSimple'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','max' => '1','if' => '$UseBorisSimple','default' => '1','type' => 'real','name' => 'BorisClightFactor'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '
#BORISSIMPLE
T			UseBorisSimple
0.05			BorisClightFactor !Only if UseBorisSimple is true

! Default is UseBorisSimple=.false. 
! Use simplified semi-relativistic MHD with speed of light reduced by the
! BorisClightFactor. This is only useful with BorisClightFactor less than 1.
! Should give the same steady state as normal MHD, but there can be a
! difference due to discretization errors.
! You can use either Boris or BorisSimple but not both. 
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'DIVB'},'content' => [{'attrib' => {'default' => 'T','type' => 'logical','name' => 'UseDivbSource'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => 'F','type' => 'logical','name' => 'UseDivbDiffusion'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => 'F','type' => 'logical','name' => 'UseProjection'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => 'F','type' => 'logical','name' => 'UseConstrainB'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'expr' => '$UseDivbSource or $UseDivbDiffusion or $UseProjection or $UseConstrainB'},'content' => [{'content' => '
		At least one of the options should be true.
	','type' => 't'}],'type' => 'e','name' => 'rule'},{'attrib' => {'expr' => 'not($UseProjection and ($UseDivbSource or $UseDivbDiffusion or $UseConstrainB))'},'content' => [{'content' => '
		If UseProjection is true, all others should be false.
	','type' => 't'}],'type' => 'e','name' => 'rule'},{'attrib' => {'expr' => 'not($UseConstrainB and ($UseDivbSource or $UseDivbDiffusion or $UseProjection))'},'content' => [{'content' => '
		If UseConstrainB is true, all others should be false.
	','type' => 't'}],'type' => 'e','name' => 'rule'},{'content' => '
	
#DIVB
T			UseDivbSource
F			UseDivbDiffusion	
F			UseProjection           
F			UseConstrainB           

! Default values are shown above.
! If UseProjection is true, all others should be false.
! If UseConstrainB is true, all others should be false.
! At least one of the options should be true.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'DIVBSOURCE'},'content' => [{'attrib' => {'default' => 'T','type' => 'logical','name' => 'UseB0Source'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '
#DIVBSOURCE
T			UseB0Source

! Add extra source terms related to the non-zero divergence and curl of B0.
! Default is true.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'DIVBDIFFUSION'},'content' => [{'attrib' => {'min' => '0','max' => '1','default' => '0.1666667','type' => 'real','name' => 'DivbDiffCoeff'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '
#DIVBDIFFUSION
0.1666667		DivbDiffCoeff

! Default value is shown above. 1.0/6.0
! If divb_diffcoeff .gt. 0.5 then cfl .lt. 0.5/DivbDiffCoeff is required!
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'PROJECTION'},'content' => [{'attrib' => {'input' => 'select','type' => 'string','name' => 'TypeProjectIter'},'content' => [{'attrib' => {'value' => 'cg','default' => 'T','name' => 'Conjugate Gradients'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'bicgstab','name' => 'BiCGSTAB'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'parameter'},{'attrib' => {'input' => 'select','type' => 'string','name' => 'TypeProjectStop'},'content' => [{'attrib' => {'value' => 'rel','default' => 'T','name' => 'Relative norm'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'max','name' => 'Maximum error'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','max' => '1','default' => '0.1','type' => 'real','name' => 'RelativeLimit'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','default' => '0.0','type' => 'real','name' => 'AbsoluteLimit'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '1','default' => '50','type' => 'integer','name' => 'MaxMatvec'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '
#PROJECTION
cg			TypeProjectIter:\'cg\' or \'bicgstab\' for iterative scheme
rel			TypeProjectStop:\'rel\' or \'max\' error for stop condition
0.1			RelativeLimit
0.0			AbsoluteLimit 
50			MaxMatvec (upper limit on matrix.vector multipl.)

! Default values are shown above.\\\\
!
!\\noindent
! For symmetric Laplacian matrix TypeProjectIter=\'cg\' (Conjugate Gradients)
! should be used, as it is faster than BiCGSTAB. In current applications
! the Laplacian matrix is always symmetric.\\\\
!
!\\noindent 
! The iterative scheme stops when the stopping condition is fulfilled:
!\\begin{verbatim}
!   TypeProjectStop = \'rel\': 
!        stop if ||div B||    < RelativeLimit*||div B0||
!   TypeProjectStop = \'max\' and RelativeLimit is positive:
!        stop if max(|div B|) < RelativeLimit*max(|div B0|)
!   TypeProjectStop = \'max\' and RelativeLimit is negative: 
!        stop if max(|div B|) < AbsoluteLimit
!\\end{verbatim}
!   where {\\tt ||.||} is the second norm, and B0 is the magnetic
!   field before projection. In words \'rel\' means that the norm of the error
!   should be decreased by a factor of RelativeLimit, while 
!   \'max\' means that the maximum error should be less than either
!   a fraction of the maximum error in div B0, or less than the constant 
!   AbsoluteLimit.
! 
!   Finally the iterations stop if the number of matrix vector
!   multiplications exceed MaxMatvec. For the CG iterative scheme
!   there is 1 matvec per iteration, while for BiCGSTAB there are 2/iteration.
!
!  In practice reducing the norm of the error by a factor of 10 to 100 in 
!  every iteration works well.
!

!
!  Projection is also used when the scheme switches to constrained transport.
!  It is probably a good idea to allow many iterations and require an
!  accurate projection, because it is only done once, and the constrained
!  transport will carry along the remaining errors in div B. An example is
!
#PROJECTION
cg			TypeProjIter
rel			TypeProjStop
0.0001			RelativeLimit
0.0			AbsoluteLimit 
500			MaxMatvec


','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'CORRECTP'},'content' => [{'attrib' => {'min' => '0','max' => '1','default' => '0.01','type' => 'real','name' => 'pRatioLow'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '$pRatioLow','max' => '1','default' => '0.1','type' => 'real','name' => 'pRatioHigh'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '
#CORRECTP
0.01			pRatioLow
0.1			pRatioHigh

Default values are shown. 

The purpose of the correctP subroutine is to remove any discrepancies between
pressure stored as the primitive variable P and the pressure calculated 
from the total energy E. Such discrepancies can be caused by the 
constrained transport scheme and by the projection scheme which modify 
the magnetic energy. The algorithm is the following:
\\begin{verbatim}
Define the rato of thermal and total energies q = eThermal/e and

If              q < pRatioLow   then E is recalculated from P
If pRatioLow  < q < pRatioHigh  then both P and E are modified depending on q
If pratioHigh < q               then P is recalculated from E
\\end{verbatim}
The second case is a linear interpolation between the first and third cases.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'RAYTRACE'},'content' => [{'attrib' => {'default' => 'T','type' => 'logical','name' => 'UseAccurateIntegral'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => 'F','type' => 'logical','name' => 'UseAccurateTrace'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0.01','max' => '60','default' => '0.1','type' => 'real','name' => 'DtExchangeRay'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '1','default' => '1','type' => 'integer','name' => 'DnRaytrace'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

#RAYTRACE
T			UseAccurateIntegral
T			UseAccurateTrace
0.1			DtExchangeRay [sec]
1			DnRaytrace

Raytracing (field-line tracing) is needed to couple the GM and IM components.
It can also be used to create plot files with open-closed field line 
information. There are two algorithms implemented for integrating rays
and for tracing rays.

If UseAccurateIntegral is true (default), the field line integrals
are calculated with the accurate algorithm, which follows the lines
all the way. If UseAccurateIntegral is false, the block-wise algorithm
is used, which actually needs the face values computed by the 
block-wise tracing algorithm (UseAccurateTrace must be false).

If UseAccurateTrace is false (default), the block-wise algorithm is used,
which interpolates at block faces. This algorithm is fast, but less 
accurate than the other algorithm. If UseAccurateTrace is true, 
the field lines are followed all the way. It is more accurate but 
potentially slower than the other algorithm.

In the accurate tracing algorithms, when the ray exits the domain that belongs 
to the PE, its information is sent to the other PE where the ray continues. 
The information is buffered for sake of efficiency and to synchronize
communication. The frequency of the information exchanges 
(in terms of CPU seconds) is given by the DtExchangeRay parameter. 
This is an optimization parameter for speed. Very small values of DtExchangeRay
result in many exchanges with few rays, while very large values result
in infrequent exchanges thus some PE-s may become idle (no more work to do).
The optimal value is problem dependent. A typically acceptable value is 
DtExchangeRay = 0.1 seconds (default).

The DnRaytrace parameter contains the minimum number of iterations between
two ray tracings. The default value 1 means that every new step requires
a new trace (since the magnetic field is changing). A larger value implies
that the field does not change significantly in that many time steps.
The ray tracing is always redone if the grid changes due to an AMR.

Default values are UseAccurateIntegral = .true., UseAccurateTrace = .false.,
DtExchangeRay = 0.1 and DnRaytrace=1.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'IM'},'content' => [{'attrib' => {'min' => '1','type' => 'real','name' => 'TauCoupleIm'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

#IM
20.0			TauCoupleIm

Same as command IMCOUPLING, except it does not read the two logicals.
See description for command IMCOUPLING.

The default value is TauCoupleIm=20.0, which corresponds to typical nudging.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'IMCOUPLING'},'content' => [{'attrib' => {'min' => '1','type' => 'real','name' => 'TauCoupleIm'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => 'T','type' => 'logical','name' => 'DoCoupleImPressure'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => 'F','type' => 'logical','name' => 'DoCoupleImDensity'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

#IMCOUPLING
20.0			TauCoupleIm
T			DoCoupleImPressure
F			DoCoupleImDensity

TauCoupleIm:

Determine how fast the GM pressure p (and possibly density rho) 
should be nudged towards the IM pressure pIm (and density dIM). 
The RCM prssure and denstiy are updated every time IM->GM coupling occurs,
but the nudging towards these values is done in every GM time step.

If not time accurate, a weighted average is taken: 

p\' = (p*TauCoupleIm + pIm)/(TauCoupleIm+1)

Therefore the larger TauCoupleIm is the slower the adjustment will be.
It takes approximately 2*TauCoupleIm time steps to get p close to pIm.

If time accurate, the nudging is based on physical time:

p\' = p + max(1.0, dt/TauCoupleIm)*(pIm - p)

where dt is the time step. It takes about 2*TauCoupleIm seconds
to get p close to pIm. If the time step dt exceeds TauCoupleIm, 
p\' = pIm is set in a single step.

The default value is TauCoupleIm=20.0, which corresponds to typical nudging.

DoCoupleImPressure:

Logical which sets whether GM pressure is driven by IM pressure, default 
is true, and it should always be true, because pressure is the dominant
variable in the IM to GM coupling.

DoCoupleImDensity:

Logical which sets whether GM density is driven by IM density, default is
false. This is a new feature in the IM-GM coupling, which is not fully
tested.
','type' => 't'}],'type' => 'e','name' => 'command'}],'type' => 'e','name' => 'commandgroup'},{'attrib' => {'name' => 'PHYSICS PARAMETERS'},'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  PHYSICS PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
','type' => 't'},{'attrib' => {'if' => '$_IsFirstSession','name' => 'GAMMA'},'content' => [{'attrib' => {'min' => '1','default' => '1.6666666667','type' => 'real','name' => 'Gamma'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '
#GAMMA
1.6666666667		Gamma

! The adiabatic index (ratio of the specific heats for fixed pressure
! and volume. The default value is 5.0/3.0, which is valid for
! monoatomic gas or plasma.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'SHOCKTUBE'},'content' => [{'attrib' => {'min' => '0','default' => '1','type' => 'real','name' => 'RhoLeft'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '0','type' => 'real','name' => 'UnLeft'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '0','type' => 'real','name' => 'Ut1Left'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '0','type' => 'real','name' => 'Ut2Left'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '0.75','type' => 'real','name' => 'BnLeft'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '1','type' => 'real','name' => 'Bt1Left'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '0','type' => 'real','name' => 'Bt2Left'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','default' => '1','type' => 'real','name' => 'pRight'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','default' => '0.125','type' => 'real','name' => 'RhoRight'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '0','type' => 'real','name' => 'UnRight'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '0','type' => 'real','name' => 'Ut1Right'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '0','type' => 'real','name' => 'Ut2Right'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '0.75','type' => 'real','name' => 'BnRight'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '-1','type' => 'real','name' => 'Bt1Right'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '0','type' => 'real','name' => 'Bt2Right'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','default' => '0.1','type' => 'real','name' => 'pRight'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'input' => 'select','type' => 'real','name' => 'ShockSlope'},'content' => [{'attrib' => {'value' => '0','default' => 'T','name' => 'no rotation'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '0.25'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '0.3333333333333'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '0.5'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '1'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '2'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '3'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '4'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'parameter'},{'content' => '
#SHOCKTUBE
1.		rho (left state)
0.		Ux (Un)
0.		Uy (Ut1)
0.		Uz (Ut2)
0.75		Bx (Bn)
1.		By (Bt1)
0.		Bz (Bt2)
1.		P
0.125		rho (right state)
0.		Ux (Un)
0.		Uy (Ut1)
0.		Uz (Ut2)
0.75		Bx (Bn)
-1.		By (Bt1)
0.		Bz (Bt2)
0.1		P
0.0		ShockSlope

! Default values are shown (Brio-Wu problem).
! The shock is rotated if ShockSlope is not 0, and the tangent of 
! the rotation angle is ShockSlope. 
! When the shock is rotated, it is best used in combination
! with sheared outer boundaries, but then only
!\\begin{verbatim}
! ShockSlope = 1., 2., 3., 4., 5.      .....
! ShockSlope = 0.5, 0.33333333, 0.25, 0.2, .....
!\\end{verbatim}
! can be used, because these angles can be accurately represented
! on the grid.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'if' => '$_IsFirstSession','name' => 'SOLARWIND'},'content' => [{'attrib' => {'min' => '-1','default' => '5','type' => 'real','name' => 'SwRhoDim'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '-1','default' => '181712.175','type' => 'real','name' => 'SwTDim'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'max' => '0','default' => '-400','type' => 'real','name' => 'SwUxDim'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '0','type' => 'real','name' => 'SwUyDim'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '0','type' => 'real','name' => 'SwUzDim'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '0','type' => 'real','name' => 'SwBxDim'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '0','type' => 'real','name' => 'SwByDim'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '5','type' => 'real','name' => 'SwBzDim'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '
#SOLARWIND
5.0			SwRhoDim [n/cc]
181712.175		SwTDim [K]
-400.0			SwUxDim [km/s]
0.0			SwUyDim [km/s]
0.0			SwUzDim [km/s]
0.0			SwBxDim [nT]
0.0			SwByDim [nT]
5.0			SwBzDim [nT]

! This command defines the solar wind parameters for the GM component.
! It also defines the normalization for all the variables therefore
! it is saved into the restart header file.
! One of the #SOLARWIND command and the #UPSTREAM_INPUT_FILE command
! (with UseUpstreamInputFile = .true.) is required by the GM component.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'UPSTREAM_INPUT_FILE'},'content' => [{'attrib' => {'default' => 'F','type' => 'logical','name' => 'UseUpstreamInputFile'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'expr' => '$UseUpstreamInputFile'},'content' => [{'attrib' => {'length' => '100','type' => 'string','name' => 'NameUpstreamFile'},'content' => [],'type' => 'e','name' => 'parameter'}],'type' => 'e','name' => 'if'},{'attrib' => {'expr' => '-f $NameUpstreamFile'},'content' => [{'content' => '
		Upstream file $NameUpstreamFile must exist
	','type' => 't'}],'type' => 'e','name' => 'rule'},{'content' => '
#UPSTREAM_INPUT_FILE
T			UseUpstreamInputFile (rest of parameters read if true)
IMF.dat                 NameUpstreamFile

! Default is UseUpstreamInputFile = .false.
!
! Read IMF data from file NameUpstreamFile if UseUpstreamInputFile is true.
! The data file contains all information required for setting the upstream
! boundary conditions. Parameter TypeBcEast should be set to \'vary\' for
! the time dependent boundary condition.
!
! If the #SOLARWIND command is not provided then the first time read from
! the upstream input file will set the normalization of all variables
! in the GM component. Consequently either the #SOLARWIND command or
! the #UPSTREAM_INPUT_FILE command with UseUpstreamInputFile=.true.
! is required by the GM component.
!
! The input files are strutured similar to the PARAM.in file.  There are 
! {\\tt #commands} that can be  inserted as well as the data.
! The file containing the upstream conditions should include data in the 
! following order:
! \\begin{verbatim}
! yr mn dy hr min sec msec bx by bz vx vy vz dens temp
! \\end{verbatim}
! The units of the variables should be:
! \\begin{verbatim}
! Magnetic field (b)     nT
! Velocity (v)           km/s
! Number Density (dens)  cm^-3
! Temperature (Temp)     K
! \\end{verbatim}
!
! The input files  can have the following optional commands at the beginning
! \\begin{verbatim}
! #COOR
! GSM          The coordinate system of the data (GSM or GSE)
!
! #PLANE       The input data represents values on a tilted plane
! 30.0         Angle to rotate in the XY plane
! 30.0         Angle to rotate in the XZ plane
!
! #POSITION    Origin for the plane rotation (see #PLANE)
! 30.0         Y location
! 30.0         Z location
!
! #ZEROBX
! T            T means Bx is ignored and set to zero
!
! #TIMEDELAY
! 3600.0       A real number in seconds by which to delay the input
! \\end{verbatim}
!
! Finally, the data should be preceded by a {\\tt #START}.  The beginning of 
! a typical upstream input file might look like:
! \\begin{verbatim}
! #COOR
! GSM
! 
! #ZEROBX
! T
!
! #START
!  2004  6  24   0   0  58   0  2.9  -3.1 - 3.7  -300.0  0.0  0.0  5.3  2.00E+04
!  2004  6  24   0   1  58   0  3.0  -3.2 - 3.6  -305.0  0.0  0.0  5.4  2.01E+04
! \\end{verbatim}
!
! The maximum number of lines of data allowed in the input file is 50,000.  
! However, this can be modified by changing the variable Max_Upstream_Npts 
! in the file GM/BATSRUS/get_solar_wind_point.f90.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'alias' => 'MAGNETOSPHERE','if' => '$_IsFirstSession','name' => 'BODY'},'content' => [{'attrib' => {'default' => 'F','type' => 'logical','name' => 'UseBody'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'expr' => '$UseBody'},'content' => [{'attrib' => {'min' => '0','default' => '3','type' => 'real','name' => 'rBody'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'expr' => '$_NameComp eq \'GM\''},'content' => [{'attrib' => {'min' => '-1','default' => '4','type' => 'real','name' => 'rCurrents'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','default' => '1','type' => 'real','name' => 'BodyRhoDim'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','default' => '10000','type' => 'real','name' => 'BodyTDim'},'content' => [],'type' => 'e','name' => 'parameter'}],'type' => 'e','name' => 'if'}],'type' => 'e','name' => 'if'},{'content' => '
#BODY
T			UseBody (rest of parameters read if true)
3.0			rBody (user units)
4.0			rCurrents
1.0			BodyRhoDim (/ccm) density for fixed BC for rho_BLK
10000.0			BodyTDim (K) temperature for fixed BC for P_BLK

! If UseBody is true, the inner boundary is a spherical surface
! with radius rBody. The rBody is defined in units of the planet/solar
! radius. It can be 1.0, in which case the simulation extends all the
! way to the surface of the central body. In many cases it is more
! economic to use an rBody larger than 1.0. 
!
! The rCurrents parameter defines where the currents are calculated for
! the GM-IE coupling. This only matters if BATSRUS is running as GM
! and it is coupled to IE.
!
! The BodyRhoDim and BodyTDim parameters define the density and temperature
! at the inner boundary. The exact effect of these parameters depends 
! on the settings in the #INNERBOUNDARY command.
! 
! The default values depend on the problem type defined 
! in the #PROBLEMTYPE command.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'if' => '$_IsFirstSession','name' => 'GRAVITY'},'content' => [{'attrib' => {'default' => 'F','type' => 'logical','name' => 'UseGravity'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'input' => 'select','if' => '$UseGravity','type' => 'integer','name' => 'iDirGravity'},'content' => [{'attrib' => {'value' => '0','default' => 'T','name' => 'central mass'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '1','name' => 'X direction'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '2','name' => 'Y direction'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => '3','name' => 'Z direction'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'parameter'},{'content' => '
#GRAVITY
T			UseGravity (rest of parameters read if true)
0			iDirGravity(0 - central, 1 - X, 2 - Y, 3 - Z direction)

! If UseGravity is false, the gravitational force of the central body
! is neglected. If UseGravity is true and iDirGravity is 0, the
! gravity points towards the origin. If iDirGravity is 1, 2 or 3,
! the gravitational force is parallel with the X, Y or Z axes, respectively.
!
! Default values depend on problem_type.

! When a second body is used the gravity direction for the second body
! is independent of the GravityDir value.  Gravity due to the second body
! is radially inward toward the second body.

','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'MASSLOADING'},'content' => [{'attrib' => {'default' => 'F','type' => 'logical','name' => 'UseMassLoading'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => 'F','type' => 'logical','name' => 'DoAccelerateMassLoading'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '
#MASSLOADING
F			UseMassLoading
F			DoAccelerateMassLoading
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'HEATFLUX'},'content' => [{'attrib' => {'default' => 'F','type' => 'logical','name' => 'UseHeatFlux'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => 'T','type' => 'logical','name' => 'UseSpitzerForm'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'expr' => 'not $UseSpitzerForm'},'content' => [{'attrib' => {'default' => '1.23E-11','type' => 'real','name' => 'Kappa0Heat'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '2.5','type' => 'real','name' => 'Kappa0Heat'},'content' => [],'type' => 'e','name' => 'parameter'}],'type' => 'e','name' => 'if'},{'content' => '
#HEATFLUX
T		UseHeatFlux
F		UseSpitzerForm
1.23E-11	Kappa0Heat [W/m/K]	! Only if not UseSpitzerForm
2.50E+00	ExponentHeat [-]	! Only if not UseSpitzerForm
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'name' => 'RESISTIVEFLUX'},'content' => [{'attrib' => {'default' => 'F','type' => 'logical','name' => 'UseResistFlux'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => 'T','type' => 'logical','name' => 'UseSpitzerForm'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'expr' => 'not $UseSpitzerForm'},'content' => [{'attrib' => {'input' => 'select','type' => 'string','name' => 'TypeResist'},'content' => [{'attrib' => {'value' => 'Localized/localized','name' => 'localized'},'content' => [],'type' => 'e','name' => 'option'},{'attrib' => {'value' => 'Constant/constant','default' => 'T','name' => 'constant'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '9.69953E+8','type' => 'real','name' => 'Eta0Resist'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'expr' => '$TypeResist =~ /localized/i'},'content' => [{'attrib' => {'default' => '150','type' => 'real','name' => 'Alpha0Resist'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '0.5','type' => 'real','name' => 'yShiftResist'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '0.05','type' => 'real','name' => 'TimeInitRise'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '1','type' => 'real','name' => 'TimeConstLev'},'content' => [],'type' => 'e','name' => 'parameter'}],'type' => 'e','name' => 'if'}],'type' => 'e','name' => 'if'},{'attrib' => {'default' => 'F','type' => 'logical','name' => 'UseAnomResist'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'expr' => '$UseAnomResist'},'content' => [{'attrib' => {'default' => '1.93991E+09','type' => 'real','name' => 'Eta0AnomResist'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '1.93991E+10','type' => 'real','name' => 'EtaAnomMaxResist'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '1','type' => 'real','name' => 'ThresholdFactorResist'},'content' => [],'type' => 'e','name' => 'parameter'}],'type' => 'e','name' => 'if'},{'content' => '
#RESISTIVEFLUX
T		UseResistFlux
F		UseSpitzerForm
Localized	TypeResist		! Only if not UseSpitzerForm
9.69953E+08	Eta0Resist [m^2/s]	! Only if not UseSpitzerForm
1.50000E+02	Alpha0Resist [-]	! Only if TypeResist==\'Localized\'
5.00000E-01	yShiftResist [-]	! Only if TypeResist==\'Localized\'
5.00000E-02	TimeInitRise [-]	! Only if TypeResist==\'Localized\'
1.00000E+00	TimeConstLev [-]	! Only if TypeResist==\'Localized\'
T		UseAnomResist
1.93991E+09	Eta0AnomResist [m^2/s]		! Only if UseAnomResist
1.93991E+10	EtaAnomMaxResist [m^2/s]	! Only if UseAnomResist
1.00000E+00	ThresholdFactorResist [-]	! Only if UseAnomResist

! Note: ResistType = `Constant\'  (the same as `constant\')
!		     \'Localized\' (the same as \'localized\')
!
! The first choice results in a uniform resistivity of value Eta0Resist.
! The second choice represents localized in space magnetic diffusion 
! with a peak value Eta0Resist. The enhanced resistivity has a Gaussian 
! shape with HWHM of 1/sqrt(Alpha0Resist), shifted along the y-axis on 
! -yShistResist*y2.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'if' => '$_IsFirstSession','name' => 'TESTDISSMHD'},'content' => [{'attrib' => {'default' => 'T','type' => 'logical','name' => 'UseDefaultUnits'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '2.635620E-02','type' => 'real','name' => 'Grav0Diss'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '1.640000E-01','type' => 'real','name' => 'Beta0Diss'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '1.500000E+06','type' => 'real','name' => 'Length0Diss'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '1.159850E+01','type' => 'real','name' => 'Time0Diss'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '5.019000E-11','type' => 'real','name' => 'Rho0Diss'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '1.000000E+05','type' => 'real','name' => 'Tem0Diss'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '6.000000E-01','type' => 'real','name' => 'Theta0Diss'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '2.500000E+01','type' => 'real','name' => 'Delta0Diss'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '7.000000E+00','type' => 'real','name' => 'EpsilonDiss'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '4.500000E+00','type' => 'real','name' => 'RhoDifDiss'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '4.000000E-01','type' => 'real','name' => 'yShiftDiss'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '5.000000E-01','type' => 'real','name' => 'ScaleHeightDiss'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '1.159850E+01','type' => 'real','name' => 'ScaleFactorDiss'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '5.000000E-01','type' => 'real','name' => 'BZ0iss'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '
#TESTDISSMHD
F                       UseDefaultUnits
2.635620E-02            Grav0Diss
1.640000E-01            Beta0Diss
1.500000E+06            Length0Diss
1.159850E+01            Time0Diss
5.019000E-11            Rho0Diss
1.000000E+05            Tem0Diss
6.000000E-01            ThetaDiss
2.500000E+01            DeltaDiss
7.000000E+00            EpsilonDiss
4.500000E+00            RhoDifDiss
4.000000E-01            yShiftDiss
5.000000E-01            scaleHeightDiss
1.000000E+00            scaleFactorDiss
0.000000E-01            BZ0Diss

! Default values are shown. Parameters for problem_dissipation
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'if' => '$_IsFirstSession','name' => 'SECONDBODY'},'content' => [{'attrib' => {'default' => 'F','type' => 'logical','name' => 'UseBody2'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'expr' => '$UseBody2'},'content' => [{'attrib' => {'min' => '0','default' => '0.1','type' => 'real','name' => 'rBody2'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '$xMin','max' => '$xMax','default' => '-40','type' => 'real','name' => 'xBody2'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '$yMin','max' => '$yMax','default' => '0','type' => 'real','name' => 'yBody2'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '$zMin','max' => '$zMax','default' => '0','type' => 'real','name' => 'zBody2'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '$rBody2','default' => '1.3*$rBody2','type' => 'real','name' => 'rCurrents2'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','default' => '5','type' => 'real','name' => 'RhoDimBody2'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','default' => '25000','type' => 'real','name' => 'tDimBody2'},'content' => [],'type' => 'e','name' => 'parameter'}],'type' => 'e','name' => 'if'},{'content' => '

#SECONDBODY
T			UseBody2   ! Rest of the parameters read if .true.
0.01			rBody2 
-40.			xBody2
0.			yBody2
0.			zBody2
0.011                   rCurrents2  !This is unused currently 
5.0			RhoDimBody2 (/ccm) density for fixed BC for rho_BLK
25000.0			TDimBody2 (K) temperature for fixed BC for P_BLK

! Default for UseBody2=.false.   -   All others have no defaults!
! This command should appear before the #INNERBOUNDARY command when using
! a second body.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'if' => '$_IsFirstSession','name' => 'DIPOLEBODY2'},'content' => [{'attrib' => {'type' => 'real','name' => 'BdpDimBody2x'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'type' => 'real','name' => 'BdpDimBody2y'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'type' => 'real','name' => 'BdpDimBody2z'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

#DIPOLEBODY2
0.0			BdpDimBody2x [nT]
0.0			BdpDimBody2y [nT]
-1000.0			BdpDimBody2z [nT]

! The BdpDimBody2x, BdpDimBody2y and BdpDimBody2z variables contain
! the 3 components of the dipole vector in the GSE frame.
! The absolute value of the dipole vector is the equatorial field strength
! in nano Tesla.
!
! Default is no dipole field.

! For now the dipole of the second body can only be aligned with the z-axis
','type' => 't'}],'type' => 'e','name' => 'command'}],'type' => 'e','name' => 'commandgroup'},{'attrib' => {'name' => 'SOLAR PROBLEM TYPES'},'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!! SOLAR PROBLEM TYPES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
','type' => 't'},{'attrib' => {'if' => '$_IsFirstSession and $_NameComp ne \'GM\'','name' => 'HELIOSPHERE'},'content' => [{'attrib' => {'min' => '0','default' => '2.85E06','type' => 'real','name' => 'BodyTDim'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','default' => '1.50E8','type' => 'real','name' => 'BodyRhoDim'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','default' => '25.0','type' => 'real','name' => 'qSun'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','default' => '1.75','type' => 'real','name' => 'tHeat'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','default' => '1.0','type' => 'real','name' => 'rHeat'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','default' => '4.5','type' => 'real','name' => 'SigmaHeat'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => 'F','type' => 'logical','name' => 'DoInitRope'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'expr' => '$DoInitRope'},'content' => [{'attrib' => {'min' => '0','default' => '0.7','type' => 'real','name' => 'CmeA'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','default' => '1.2','type' => 'real','name' => 'CmeR1'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','default' => '1.0','type' => 'real','name' => 'CmeR0'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '0.23','type' => 'real','name' => 'CmeA1'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '0.0','type' => 'real','name' => 'CmeAlpha'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','default' => '2.5E-12','type' => 'real','name' => 'CmeRho1'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','default' => '2.0E-13','type' => 'real','name' => 'CmeRho2'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','max' => '10','default' => '0.0','type' => 'real','name' => 'ModulationRho'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','max' => '10','default' => '0.0','type' => 'real','name' => 'ModulationP'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '-360','max' => '360','default' => '0.0','type' => 'real','name' => 'OrientationGL98'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '-90','max' => '90','default' => '0.0','type' => 'real','name' => 'LatitudeGL98'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '-360','max' => '360','default' => '0.0','type' => 'real','name' => 'LongitudeGL98'},'content' => [],'type' => 'e','name' => 'parameter'}],'type' => 'e','name' => 'if'},{'content' => '

#HELIOSPHERE
2.85E06			BodyTDim	[K]
1.50E8			BodyRhoDim	[N/ccm]
25.00			qSun		
1.75			tHeat
1.00			rHeat
4.50			SIGMAheat
F			InitRope
0.7     		CmeA    [scaled] contraction distance
1.2     		CmeR1   [scaled] distance of spheromac from sun center
1.0     		CmeR0   [scaled] diameter of spheromac
0.23    		CmeA1   [Gauss]  spheromac B field strength
0.0     		CmeAlpha[scaled] cme acceleration rate
2.5E-12 		CmeRho1 [kg/m^3] density of background corona before contract
2.0E-13 		CmeRho2 [kg/m^3] density of background corona after contract 
0.0                     ModulationRho
0.0                     ModulationP
0.0			OrientationGL98 [deg]
0.0			LatitudeGL98 [deg]
0.0			LongitudeGL98 [deg]

This command defines the heliosphere parameters with a CME model.
The coronal eruptive event generator is based on the
Gibson and Low (GL) analytical solution prescribing a
three-dimensional twisted magnetic flux rope in
hydrostatic equilibrium in the presence of gravity.
The GL solution is described in the Astrophysical
Journal, volume 493, page 460.
This flux rope is formed by applying a mathematical
stretching operation to axisymmetric  flux
rope.  The flux rope is of radius Cme_R0 and is
placed Cme_R1 from the origin (solar center).  The
stretching transformation draws space radially inward
toward the origin by a distance of Cme_A, which
distorts the flux rope to have a tear-drop shape.
The parameter Cme_A1 modulates the magnetic field strength
and negative values of Cme_A1 reverse the overall field
direction.  For the GL flux rope to be in equilibrium, requires
both dense plasma in the form of a filament inside the rope,
(prescribed by the GL solution) as well as plasma pressure
outside the flux rope which tends to be large than the
solar corona can provide.  To initiate an eruption (the CME)
we linearly superimpose the GL flux rope in the solar
corona within the streamer belt.  The location of the flux
rope is determined by the parameters cRotxGl98, cRotYGl98
and cRotZGl98.  The flux rope is line-tied with both ends
attached to the inner boundary.  The eruption follows from
the flux rope being out of equilibrium, owing to a reduction
in filament mass (set with ModulationRho) and from pressure
of the corona being unable to balance the magnetic pressure
of the flux rope.  The eruption takes the form of the flux
rope being bodily expelled from the corona.  Eruption energy
increases with flux rope size, field strength, stretching
deformation and the buoyancy of the flux rope.

The flux rope can be rotated to an arbitrary position.
The LatitudeGM98 and LongitudeGL98 parameters define the position
of the center of the fluxrope in the coordinate system of the Solar
Corona component. The OrientationGL98 parameter determines the 
orientation of the fluxrope relative to the East-West direction
(clock-wise).
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'if' => '$_NameComp ne \'GM\'','name' => 'HELIOUPDATEB0'},'content' => [{'attrib' => {'min' => '-1','default' => '0.0001','type' => 'real','name' => 'DtUpdateB0'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

#HELIOUPDATEB0
-1.0			DtUpdateB0 [s]

Set the frequency of updating the B0 field for the solar corona.
A negative value means that the B0 field is not updated.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'if' => '$_NameComp ne \'GM\'','name' => 'HELIODIPOLE'},'content' => [{'attrib' => {'type' => 'real','name' => 'HelioDipoleStrength'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '-90','max' => '90','default' => '0','type' => 'real','name' => 'HelioDipoleTilt'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

#HELIODIPOLE
-3.0                    HelioDipoleStrength [G]
 0.0                    HelioDipoleTilt     [deg]

! Variable HelioDipoleStrength defines the equatorial field strength in Gauss,
! while HelioDipoleTilt is the tilt relative to the ecliptic North 
! (negative sign means towards the planet) in degrees.
!
! Default value is HelioDipoleStrength = 0.0.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'if' => '$_NameComp ne \'GM\'','name' => 'HELIOTEST'},'content' => [{'attrib' => {'default' => 'F','type' => 'logical','name' => 'DoSendMHD'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

#HELIOTEST
F			DoSendMHD

! If DoSendMHD is true, IH sends the real MHD solution to GM in the coupling.
! If DoSendMHD is false then the values read from the IMF file are sent,
! so there is no real coupling. Mostly used for testing the framework.
!
! Default value is true, i.e. real coupling.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'if' => '$_IsFirstSession and $_NameComp eq \'IH\'','name' => 'HELIOBUFFERGRID'},'content' => [{'attrib' => {'min' => '1','default' => '19','type' => 'real','name' => 'rBuffMin'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '1','default' => '21','type' => 'real','name' => 'rBuffMax'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '18','default' => '45','type' => 'integer','name' => 'nThetaBuff'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '36','default' => '90','type' => 'integer','name' => 'nPhiBuff'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '
#HELIOBUFFERGRID
19.0		rBuffMin
21.0		rBuffMax
45		nThetaBuff
90		nPhiBuff

Define the radius and the grid resolution for the uniform 
spherical buffer grid which passes information 
from the SC component to the IH component. The resolution should
be similar to the grid resolution of the coarser of the SC and IH grids.
This command can only be used in the first session by the IH component. 
Default values are shown above.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'if' => '$_IsFirstSession and $_NameComp ne \'GM\'','name' => 'CME'},'content' => [{'attrib' => {'input' => 'select','type' => 'string','name' => 'TypeCme'},'content' => [{'attrib' => {'default' => 'T','name' => 'Low'},'content' => [],'type' => 'e','name' => 'option'}],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','default' => '0.7','type' => 'real','name' => 'CmeA'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','default' => '1.2','type' => 'real','name' => 'CmeR1'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','default' => '1.0','type' => 'real','name' => 'CmeR0'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '0.23','type' => 'real','name' => 'CmeA1'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '0.0','type' => 'real','name' => 'CmeAlpha'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','default' => '2.5E-12','type' => 'real','name' => 'CmeRho1'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','default' => '2.0E-13','type' => 'real','name' => 'CmeRho2'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '1.0','type' => 'real','name' => 'CmeB1Dim'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','default' => '4.0E5','type' => 'real','name' => 'CmeUErupt'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '
#CME
Low		TypeCme   model type (\'Low\')
0.7		CmeA    [scaled] contraction distance
1.2             CmeR1   [scaled] distance of spheromac from sun center
1.0             CmeR0   [scaled] diameter of spheromac
0.23		CmeA1   [Gauss]  sets spheromac B strength which can be negative
0.0		Cmealpha   [scaled] cme acceleration rate
2.5E-12		CmeRho1 [kg/m^3] density of background corona before contract
2.0E-13		CmeRho2 [kg/m^3] density of background corona after contract 
1.0             CmeB1Dim [Gauss] field strength of dipole-type B field
4.0E5           CmeUErupt  [m/s] cme velocity

The coronal eruptive event generator (TypeCme Low) is based on the
Gibson and Low (GL) analytical solution prescribing a
three-dimensional twisted magnetic flux rope in
hydrostatic equilibrium in the presence of gravity.
The GL solution is described in the Astrophysical
Journal, volume 493, page 460.
This flux rope is formed by applying a mathematical
stretching operation to axisymmetric spheromak flux
rope.  The flux rope is of radius Cme_R0 and is
placed Cme_R1 from the origin (solar center).  The
stretching transformation draws space radially inward
toward the origin by a distance of Cme_A, which
distorts the flux rope to have a tear-drop shape.
The parameter Cme_A1 modulates the magnetic field strength
and negative values of Cme_A1 reverse the overall field
direction.  For the GL flux rope to be in equilibrium, requires
both dense plasma in the form of a filament inside the rope,
(prescribed by the GL solution) as well as plasma pressure
outside the flux rope which tends to be large than the
solar corona can provide.  To initiate an eruption (the CME)
we linearly superimpose the GL flux rope in the solar
corona within the streamer belt.  The location of the flux
rope is determined by the parameters cRotxGl98, cRotYGl98
and cRotZGl98.  The flux rope is line-tied with both ends
attached to the inner boundary.  The eruption follows from
the flux rope being out of equilibrium, owing to a reduction
in filament mass (set with ModulationRho) and from pressure
of the corona being unable to balance the magnetic pressure
of the flux rope.  The eruption takes the form of the flux
rope being bodily expelled from the corona.  Eruption energy
increases with flux rope size, field strength, stretching
deformation and the buoyancy of the flux rope.
Default values are shown above for the GL flux rope CME model.
','type' => 't'}],'type' => 'e','name' => 'command'},{'attrib' => {'if' => '$_IsFirstSession and $_NameComp ne \'GM\'','name' => 'ARCADE'},'content' => [{'attrib' => {'min' => '0','default' => '1.0E6','type' => 'real','name' => 'tArcDim'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','default' => '1.0E-12','type' => 'real','name' => 'RhoArcDim'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','default' => '0.718144','type' => 'real','name' => 'bArcDim'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '1.0E6','type' => 'real','name' => 'ByArcDim'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '5.0E3','type' => 'real','name' => 'UzArcDim'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '0.5','type' => 'real','name' => 'Phi0Arc'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'default' => '1.3','type' => 'real','name' => 'MuArc'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','default' => '3','type' => 'real','name' => 'ExpArc'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','default' => '0.5','type' => 'real','name' => 'WidthArc'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '
#ARCADE
1.0E6                   tArcDim   [K]      1.0E6
1.0E-12                 RhoArcDim [kg/m^3] 1.0E-12
0.71814                 bArcDim   [Gauss]  0.718144
0.0                     ByArcDim  [Gauss]
5.0E3                   UzArcDim  [5.0E3 m/s]
0.5                     Phi0Arc
1.3                     MuArc
3                       ExpArc
0.5                     WidthArc

! Default values are shown. Parameters for problem_arcade
','type' => 't'}],'type' => 'e','name' => 'command'}],'type' => 'e','name' => 'commandgroup'},{'attrib' => {'name' => 'COMET PROBLEM TYPE'},'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!! COMET PROBLEM TYPE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
','type' => 't'},{'attrib' => {'if' => '$_IsFirstSession','name' => 'COMET'},'content' => [{'attrib' => {'min' => '0','type' => 'real','name' => 'ProdRate'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','type' => 'real','name' => 'UrNeutral'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','type' => 'real','name' => 'AverageMass'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','type' => 'real','name' => 'IonizationRate'},'content' => [],'type' => 'e','name' => 'parameter'},{'attrib' => {'min' => '0','type' => 'real','name' => 'kFriction'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '
#COMET
1.0E28		ProdRate    - Production rate (#/s)
1.0		UrNeutral   - neutral radial outflow velocity (km/s)
17.0		AverageMass - average particle mass (amu)
1.0E-6		IonizationRate (1/s)
1.7E-9		kFriction - ion-neutral friction rate coefficient (cm^3/s)

! Only used by problem_comet.  Defaults are as shown.
','type' => 't'}],'type' => 'e','name' => 'command'}],'type' => 'e','name' => 'commandgroup'},{'attrib' => {'name' => 'SCRIPT COMMANDS'},'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! SCRIPT COMMANDS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
','type' => 't'},{'attrib' => {'name' => 'INCLUDE'},'content' => [{'attrib' => {'length' => '100','default' => 'Param/','type' => 'string','name' => 'NameIncludeFile'},'content' => [],'type' => 'e','name' => 'parameter'},{'content' => '

#INCLUDE
Param/SSS_3000		NameIncludeFile

! Include a library file from Param/ or any file from anywhere else.
','type' => 't'}],'type' => 'e','name' => 'command'}],'type' => 'e','name' => 'commandgroup'},{'attrib' => {'expr' => '($SwRhoDim > 0) or $UseUpstreamInputFile or $_NameComp ne \'GM\''},'content' => [{'content' => '
	Either command #SOLARWIND or #UPSTREAM_INPUT_FILE must be used!
','type' => 't'}],'type' => 'e','name' => 'rule'},{'attrib' => {'expr' => '$MaxImplBlock>1 or not $UsePartImplicit or not $MaxImplBlock'},'content' => [{'content' => '
	Part implicit scheme requires more than 1 implicit block!
','type' => 't'}],'type' => 'e','name' => 'rule'},{'attrib' => {'expr' => '$MaxImplBlock==$MaxBlock or not $UseFullImplicit'},'content' => [{'content' => '
	Full implicit scheme should be used with equal number of 
	explicit and implicit blocks!
','type' => 't'}],'type' => 'e','name' => 'rule'},{'attrib' => {'expr' => '-d $NameRestartOutDir or not $_IsFirstSession'},'content' => [{'content' => '
	Output restart directory $NameRestartOutDir should exist!
','type' => 't'}],'type' => 'e','name' => 'rule'},{'attrib' => {'expr' => '-d $NamePlotDir or not $_IsFirstSession'},'content' => [{'content' => '
	Plot directory $NamePlotDir should exist!
','type' => 't'}],'type' => 'e','name' => 'rule'}],'type' => 'e','name' => 'commandList'}];