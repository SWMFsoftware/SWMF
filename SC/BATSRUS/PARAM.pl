#^CFG FILE _FALSE_
$tree = [{'type' => 'e','content' => [{'type' => 't','content' => '

List of MH (GM, IH and SC) commands used in the PARAM.in file




'},{'type' => 'e','content' => [],'attrib' => {'type' => 'integer','value' => '$_GridSize[0]','name' => 'nI'},'name' => 'set'},{'type' => 'e','content' => [],'attrib' => {'type' => 'integer','value' => '$_GridSize[1]','name' => 'nJ'},'name' => 'set'},{'type' => 'e','content' => [],'attrib' => {'type' => 'integer','value' => '$_GridSize[2]','name' => 'nK'},'name' => 'set'},{'type' => 'e','content' => [],'attrib' => {'type' => 'integer','value' => '$_GridSize[3]','name' => 'MaxBlock'},'name' => 'set'},{'type' => 'e','content' => [],'attrib' => {'type' => 'integer','value' => '$_GridSize[4]','name' => 'MaxImplBlock'},'name' => 'set'},{'type' => 'e','content' => [],'attrib' => {'type' => 'integer','value' => '$_nProc and $MaxBlock and $_nProc*$MaxBlock','name' => 'MaxBlockALL'},'name' => 'set'},{'type' => 'e','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! STAND ALONE PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'T','name' => 'UseNewParam'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'T','name' => 'UseNewAxes'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'T','name' => 'DoTimeAccurate'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'T','name' => 'UseCorotation'},'name' => 'parameter'},{'type' => 't','content' => '

#NEWPARAM
T			UseNewParam
T			UseNewAxes
T			DoTimeAccurate
T			UseCorotation

This command can be used to make the standalone code backwards compatible.

If UseNewParam is true, the time frequencies of various commands 
(SAVEPLOT, SAVELOGFILE, STOP etc.) are always read, irrespective of the value 
of DoTimeAccurate and the DoTimeAccurate logical can be set with the TIMEACCURATE command.

If UseNewParam is false, the time frequencies are only read when DoTimeAccurate is true, 
and DoTimeAccurate can be set as the first parameter of the TIMESTEPPING command.

If UseNewAxes is true, the planet\'s rotational and magnetix axes are set by the new
algorithms found in share/Library/src/CON\\_axes, the planet data is set and
stored by share/Library/src/CON\\_planet, and magnetic field information and
mapping is provided by share/Library/src/CON\\_planet_field, and the rotational speed
of the planet is calculated using $v_\\phi=\\Omega \\times r$.

If UseNewAxes is false, the original algorithms in GM/BATSRUS/src/ModCompatibility 
are used. Some of these algorithms are inaccurate, some of them contain bugs,
some of them are inefficient. The algorithms were kept for sake of backwards
compatibility.

The DoTimeAccurate and UseCorotation parameters can be set elsewhere, but their
default values can be set here. This is again useful for backwards compatibility,
since BATSRUS v7.72 and earlier has DoTimeAccurate=F and UseCorotation=F as the
default, while SWMF has the default values DoTimeAccurate=T and UseCorotation=T
(consistent with the assumption that the default behaviour is as realistic as possible).

The default values depend on how the standalone code was installed
(make install STANDALON=???). For STANDALONE=gm and STANDALONE=ih
all the logicals have true default values (consistent with SWMF), 
for STANDALONE=old and STANDALONE=oldtest the default values are false 
(consistent with BATSRUS v7.72 and earlier).
'}],'attrib' => {'if' => '$_IsStandAlone','name' => 'NEWPARAM'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'string','length' => '100','name' => 'StringDescription'},'name' => 'parameter'},{'type' => 't','content' => '

#DESCRIPTION
This is a test run for Jupiter with no rotation.

This command is only used in the stand alone mode.

The StringDescription string can be used to describe the simulation
for which the parameter file is written. The #DESCRIPTION command and
the StringDescription string are saved into the restart file,
which helps in identifying the restart files.

The default value is "Please describe me!", which is self explanatory.
'}],'attrib' => {'if' => '$_IsStandAlone','name' => 'DESCRIPTION'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'F','name' => 'DoEcho'},'name' => 'parameter'},{'type' => 't','content' => '

#ECHO
T                       DoEcho

This command is only used in the stand alone mode.

If the DoEcho variable is true, the input parameters are echoed back.
The default value for DoEcho is .false., but it is a good idea to
set it to true at the beginning of the PARAM.in file.
'}],'attrib' => {'if' => '$_IsStandAlone','name' => 'ECHO'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'integer','default' => '10','min' => '-1','name' => 'DnProgressShort'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'integer','default' => '100','min' => '-1','name' => 'DnProgressLong'},'name' => 'parameter'},{'type' => 't','content' => '
#PROGRESS
10			DnProgressShort
100			DnProgressLong

The frequency of short and long progress reports for BATSRUS in
stand alone mode. These are the defaults. Set -1-s for no progress reports.
'}],'attrib' => {'if' => '$_IsStandAlone','name' => 'PROGRESS'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'T','name' => 'DoTimeAccurate'},'name' => 'parameter'},{'type' => 't','content' => '

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
to accelarate the convergence towards steady state.

The default value depends on how the stand alone code was installed.
See the description of the NEWPARAM command.
'}],'attrib' => {'if' => '$_IsStandAlone','name' => 'TIMEACCURATE'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 't','content' => '

This command is allowed in stand alone mode only for sake of the 
test suite, which contains these commands when the framework is tested.
'}],'attrib' => {'multiple' => 'T','if' => '$_IsStandAlone','name' => 'BEGIN_COMP'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 't','content' => '

This command is allowed in stand alone mode only for sake of the 
test suite, which contains these commands when the framework is tested.
'}],'attrib' => {'multiple' => 'T','if' => '$_IsStandAlone','name' => 'END_COMP'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 't','content' => '

#RUN

This command is only used in stand alone mode.

The #RUN command does not have any parameters. It signals the end
of the current session, and makes BATSRUS execute the session with
the current set of parameters. The parameters for the next session
start after the #RUN command. For the last session there is no
need to use the #RUN command, since the #END command or simply
the end of the PARAM.in file makes BATSRUS execute the last session.
'}],'attrib' => {'if' => '$_IsStandAlone','name' => 'RUN'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 't','content' => '

#END

The #END command signals the end of the included file or the
end of the PARAM.in file. Lines following the #END command are
ignored. It is not required to use the #END command. The end
of the included file or PARAM.in file is equivalent with an 
#END command in the last line.
'}],'attrib' => {'name' => 'END'},'name' => 'command'}],'attrib' => {'name' => 'STAND ALONE MODE'},'name' => 'commandgroup'},{'type' => 'e','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! PLANET COMMANDS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

The planet commands can only be used in stand alone mode and only
when UseNewAxes is set to true (see discussion at the NEWPARAM command).
The commands allow to work with an arbitrary planet.
It is also possible to change some parameters of the planet relative
to the real values.

By default Earth is assumed with its real parameters.
Another planet can be selected with the #PLANET command.
The real planet parameters can be modified and simplified
with the other planet commands listed in this subsection.
These modifier commands cannot preceed the #PLANET command!

'},{'type' => 'e','content' => [{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'default' => 'T','value' => 'EARTH/Earth/earth','name' => 'Earth'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'SATURN/Saturn/saturn','name' => 'Saturn'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'New'},'name' => 'option'}],'attrib' => {'type' => 'string','input' => 'select','name' => 'NamePlanet'},'name' => 'parameter'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'real','min' => '0','name' => 'RadiusPlanet'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','min' => '0','name' => 'MassPlanet'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','min' => '0','name' => 'OmegaPlanet'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','min' => '0','name' => 'TiltRotation'},'name' => 'parameter'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'name' => 'NONE'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'default' => 'T','name' => 'DIPOLE'},'name' => 'option'}],'attrib' => {'type' => 'string','input' => 'select','name' => 'TypeBField'},'name' => 'parameter'}],'attrib' => {'expr' => '$NamePlanet eq \'New\''},'name' => 'if'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'real','max' => '180','min' => '0','name' => 'MagAxisThetaGeo'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','max' => '360','min' => '0','name' => 'MagAxisPhiGeo'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','name' => 'DipoleStrength'},'name' => 'parameter'}],'attrib' => {'expr' => '$TyepBField eq \'DIPOLE\''},'name' => 'if'},{'type' => 'e','content' => [{'type' => 't','content' => '
		PLANET should precede $PlanetCommand
	'}],'attrib' => {'expr' => 'not $PlanetCommand'},'name' => 'rule'},{'type' => 't','content' => '

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
'}],'attrib' => {'if' => '$_IsFirstSession and $_IsStandAlone','name' => 'PLANET'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'T','name' => 'IsRotAxisPrimary'},'name' => 'parameter'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'real','max' => '180','min' => '0','name' => 'RotAxisTheta'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','max' => '360','min' => '0','name' => 'RotAxisPhi'},'name' => 'parameter'}],'attrib' => {'expr' => '$IsRotAxisPrimary'},'name' => 'if'},{'type' => 'e','content' => [],'attrib' => {'type' => 'string','value' => 'ROTATIONAXIS','name' => 'PlanetCommand'},'name' => 'set'},{'type' => 't','content' => '

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
'}],'attrib' => {'if' => '$_IsStandAlone','name' => 'ROTATIONAXIS'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'T','name' => 'UseRotation'},'name' => 'parameter'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'real','name' => 'RotationPeriod'},'name' => 'parameter'}],'attrib' => {'expr' => '$UseRotation'},'name' => 'if'},{'type' => 'e','content' => [],'attrib' => {'type' => 'string','value' => 'MAGNETICAXIS','name' => 'PlanetCommand'},'name' => 'set'},{'type' => 't','content' => '

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
'}],'attrib' => {'if' => '$_IsStandAlone','name' => 'ROTATION'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'T','name' => 'IsMagAxisPrimary'},'name' => 'parameter'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'real','max' => '180','min' => '0','name' => 'MagAxisTheta'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','max' => '360','min' => '0','name' => 'MagAxisPhi'},'name' => 'parameter'}],'attrib' => {'expr' => '$IsMagAxisPrimary'},'name' => 'if'},{'type' => 'e','content' => [],'attrib' => {'type' => 'string','value' => 'MAGNETICAXIS','name' => 'PlanetCommand'},'name' => 'set'},{'type' => 't','content' => '

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
'}],'attrib' => {'if' => '$_IsStandAlone','name' => 'MAGNETICAXIS'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'real','name' => 'DipoleStrength'},'name' => 'parameter'},{'type' => 't','content' => '

#DIPOLE
-3.11e-4		DipoleStrength [Tesla]

The DipoleStrength variable contains the
magnetic equatorial strength of the dipole magnetic field in Tesla.

The default value is the real dipole strength for the planet.
For the Earth the default is taken to be -31100 nT.
The sign is taken to be negative so that the magnetic axis can
point northward as usual.
'}],'attrib' => {'if' => '$_IsStandAlone','name' => 'DIPOLE'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '0.0001','min' => '-1','name' => 'DtUpdateB0'},'name' => 'parameter'},{'type' => 't','content' => '

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
'}],'attrib' => {'if' => '$_IsStandAlone','name' => 'UPDATEB0'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 't','content' => '

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
'}],'attrib' => {'if' => '$_IsStandAlone','name' => 'IDEALAXES'},'name' => 'command'}],'attrib' => {'name' => 'PLANET COMMANDS'},'name' => 'commandgroup'},{'type' => 'e','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!  USER DEFINED INPUT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'F','name' => 'UseUserInnerBcs'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'F','name' => 'UseUserSource'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'F','name' => 'UseUserPerturbation'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'F','name' => 'UseUserOuterBcs'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'F','name' => 'UseUserICs'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'F','name' => 'UseUserSpecifyRefinement'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'F','name' => 'UseUserLogFiles'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'F','name' => 'UseUserWritePlot'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'F','name' => 'UseUserAMR'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'F','name' => 'UseUserEchoInput'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'F','name' => 'UseUserB0'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'F','name' => 'UseUserSetPhysConst'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'F','name' => 'UseUserUpdateStates'},'name' => 'parameter'},{'type' => 't','content' => '

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
'}],'attrib' => {'name' => 'USER_FLAGS'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 't','content' => '

This command signals the beginning of the section of the file which 
is read by the subroutine user\\_read\\_inputs in the user\\_routines.f90 file.
The section ends with the #USERINPUTEND command. There is no XML based parameter
checking in the user section.
'}],'attrib' => {'name' => 'USERINPUTBEGIN'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 't','content' => '

This command signals the end of the section of the file which 
is read by the subroutine user\\_read\\_inputs in the user\\_routines.f90 file.
The section begins with the #USERINPUTBEGIN command. There is no XML based parameter
checking in the user section.
'}],'attrib' => {'name' => 'USERINPUTEND'},'name' => 'command'}],'attrib' => {'name' => 'USER DEFINED INPUT'},'name' => 'commandgroup'},{'type' => 'e','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  TESTING AND TIMING PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'string','length' => '100','name' => 'TestString'},'name' => 'parameter'},{'type' => 't','content' => '
#TEST
read_inputs

! A space separated list of subroutine names. Default is empty string.
!
! Examples:
!   read_inputs  - echo the input parameters following the #TEST line
!   message_count- count messages
!   initial_refinement
!   ...
! Check the subroutines for call setoktest("...",oktest,oktest_me) to
! see the appropriate strings.
'}],'attrib' => {'name' => 'TEST'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'integer','max' => '$nI+2','min' => '-2','name' => 'iTest'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'integer','max' => '$nJ+2','min' => '-2','name' => 'jTest'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'integer','max' => '$nK+2','min' => '-2','name' => 'kTest'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'integer','max' => '$MaxBlock','min' => '1','name' => 'iBlockTest'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'integer','min' => '0','name' => 'iProcTest'},'name' => 'parameter'},{'type' => 't','content' => '
#TESTIJK
1                       iTest           (cell index for testing)
1                       jTest           (cell index for testing)
1                       kTest           (cell index for testing)
1                       BlockTest       (block index for testing)
0                       ProcTest        (processor index for testing)

! The location of test info in terms of indices, block and processor number.
! Note that the user should set #TESTIJK or #TESTXYZ, not both.  If both
! are set, the final one in the session will set the test point.
'}],'attrib' => {'name' => 'TESTIJK'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'real','max' => '$xMax','min' => '$xMin','name' => 'xTest'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','max' => '$yMax','min' => '$yMin','name' => 'yTest'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','max' => '$zMax','min' => '$zMin','name' => 'zTest'},'name' => 'parameter'},{'type' => 't','content' => '
#TESTXYZ
1.5                     xTest           (X coordinate of cell for testing)
-10.5                   yTest           (Y coordinate of cell for testing)
-10.                    zTest           (Z coordinate of cell for testing)

! The location of test info in terms of coordinates.
! Note that the user should set #TESTIJK or #TESTXYZ, not both.  If both
! are set, the final one in the session will set the test point.
'}],'attrib' => {'name' => 'TESTXYZ'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'integer','default' => '-1','min' => '-1','name' => 'nIterTest'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '1e30','min' => '-1','name' => 'TimeTest'},'name' => 'parameter'},{'type' => 't','content' => '

#TESTTIME
-1                      nIterTest       (iteration number to start testing)
10.5                    TimeTest        (time to start testing in seconds)

! The time step and physical time to start testing.
'}],'attrib' => {'name' => 'TESTTIME'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'default' => 'T','value' => '1','name' => 'Rho'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => '2','name' => 'RhoUx'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => '3','name' => 'RhoUy'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => '4','name' => 'RhoUz'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => '5','name' => 'Bx'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => '6','name' => 'By'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => '7','name' => 'Bz'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => '8','name' => 'e'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => '9','name' => 'p'},'name' => 'option'}],'attrib' => {'type' => 'integer','input' => 'select','name' => 'iVarTest'},'name' => 'parameter'},{'type' => 't','content' => '
#TESTVAR
1                       iVarTest

! Index of variable to be tested. Default is rho_="1", ie. density.
'}],'attrib' => {'name' => 'TESTVAR'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'value' => '0','name' => 'all'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'default' => 'T','value' => '1','name' => 'x'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => '2','name' => 'y'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => '3','name' => 'z'},'name' => 'option'}],'attrib' => {'type' => 'integer','input' => 'select','name' => 'iVarTest'},'name' => 'parameter'},{'type' => 't','content' => '
#TESTDIM
1                       iDimTest

! Index of dimension/direction to be tested. Default is X dimension.
'}],'attrib' => {'name' => 'TESTDIM'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'T','name' => 'UseStrict'},'name' => 'parameter'},{'type' => 't','content' => '
#STRICT
T                       UseStrict

! If true then stop when parameters are incompatible. If false, try to
! correct parameters and continue. Default is true, ie. strict mode
'}],'attrib' => {'name' => 'STRICT'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'value' => '-1','name' => 'errors and warnings only'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => '0','name' => 'start and end of sessions'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'default' => 'T','value' => '1','name' => 'normal'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => '10','name' => 'calls on test processor'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => '100','name' => 'calls on all processors'},'name' => 'option'}],'attrib' => {'type' => 'integer','input' => 'select','name' => 'iVarTest'},'name' => 'parameter'},{'type' => 't','content' => '
#VERBOSE
-1                      lVerbose

! Verbosity level controls the amount of output to STDOUT. Default level is 1.
!   lVerbose .le. -1 only warnings and error messages are shown.
!   lVerbose .ge.  0 start and end of sessions is shown.
!   lVerbose .ge.  1 a lot of extra information is given.
!   lVerbose .ge. 10 all calls of set_oktest are shown for the test processor.
!   lVerbose .ge.100 all calls of set_oktest are shown for all processors.
'}],'attrib' => {'name' => 'VERBOSE'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'F','name' => 'DoDebug'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'F','name' => 'DoDebugGhost'},'name' => 'parameter'},{'type' => 't','content' => '
#DEBUG
F                       DoDebug         (use it as if(okdebug.and.oktest)...)
F                       DoDebugGhost    (parameter for show_BLK in library.f90)

! Excessive debug output can be controlled by the global okdebug parameter
'}],'attrib' => {'name' => 'DEBUG'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '7.50','min' => '0','name' => 'CodeVersion'},'name' => 'parameter'},{'type' => 't','content' => '
#CODEVERSION
7.50                    CodeVersion

! Cheks CodeVersion. Prints a WARNING if it differs from the CodeVersion
! defined in ModMain. Used in newer restart header files. 
! Should be given in PARAM.in when reading old restart files, 
! which do not have version info in the header file.
'}],'attrib' => {'if' => '$_IsFirstSession','name' => 'CODEVERSION'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'string','default' => 'MHD','length' => '100','name' => 'NameEquation'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'integer','default' => '8','name' => 'nVar'},'name' => 'parameter'},{'type' => 't','content' => '
#EQUATION
MHD			NameEquation
8			nVar

! Define the equation name and the number of variables.
! If any of these do not agree with the values determined 
! by the code, BATSRUS stops with an error. Used in restart
! header files and can be given in PARAM.in as a check
! and as a description.
'}],'attrib' => {'if' => '$_IsFirstSession','name' => 'EQUATION'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'default' => '$_nByteReal==4','value' => '4','name' => 'single precision (4)'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'default' => '$_nByteReal==8','value' => '8','name' => 'double precision (8)'},'name' => 'option'}],'attrib' => {'type' => 'integer','input' => 'select','name' => 'nByteReal'},'name' => 'parameter'},{'type' => 'e','content' => [{'type' => 't','content' => '
		nByteReal in file must agree with _nByteReal.
	'}],'attrib' => {'expr' => '$nByteReal==$_nByteReal'},'name' => 'rule'},{'type' => 't','content' => '

#PRECISION
8                       nByteReal

! Define the number of bytes in a real number. If it does not agree
! with the value determined by the code, BATSRUS stops with an error.
! This is a check, the internal value is calculated in parallel_setup.
! Used in latest restart header files to check binary compatibility.
! May be given in PARAM.in to enforce a certain precision.
'}],'attrib' => {'if' => '$_IsFirstSession','name' => 'PRECISION'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'integer','default' => '$nI','max' => '$nI','min' => '$nI','name' => 'nI'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'integer','default' => '$nJ','max' => '$nJ','min' => '$nJ','name' => 'nJ'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'integer','default' => '$nK','max' => '$nK','min' => '$nK','name' => 'nK'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'integer','max' => '$MaxBlockALL','min' => '1','name' => 'MinBlockALL'},'name' => 'parameter'},{'type' => 't','content' => '

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
'}],'attrib' => {'if' => '$_IsFirstSession','name' => 'CHECKGRIDSIZE'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 't','content' => '
#BLOCKLEVELSRELOADED

This command means that the restart file contains the information about
the minimum and maximum allowed refinement levels for each block.
This command is only used in the restart header file.
'}],'attrib' => {'name' => 'BLOCKLEVELSRELOADED'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'T','name' => 'UseTiming'},'name' => 'parameter'},{'type' => 'e','content' => [{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'value' => '-3','name' => 'none'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'default' => 'T','value' => '-2','name' => 'final only'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => '-1','name' => 'end of sessions'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'default' => '100','min' => '1','name' => 'every X steps'},'name' => 'optioninput'}],'attrib' => {'type' => 'integer','input' => 'select','name' => 'Frequency'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'integer','default' => '-1','min' => '-1','name' => 'nDepthTiming'},'name' => 'parameter'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'default' => '1','value' => 'cumm','name' => 'cummulative'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'list'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'tree'},'name' => 'option'}],'attrib' => {'type' => 'string','input' => 'select','name' => 'TypeTimingReport'},'name' => 'parameter'}],'attrib' => {'expr' => '$UseTiming'},'name' => 'if'},{'type' => 't','content' => '
#TIMING
T                       UseTiming      (rest of parameters read if true)
-2                      DnTiming       (-3 none, -2 final, -1 each session/AMR)
-1                      nDepthTiming   (-1 for arbitrary depth)
cumm                    TypeTimingReport   (\'cumm\', \'list\', or \'tree\')

! The default values are shown.
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
! \'cumm\' - cummulative list sorted by timings
! \'list\' - list based on caller and sorted by timings
! \'tree\' - tree based on calling sequence
'}],'attrib' => {'name' => 'TIMING'},'name' => 'command'}],'attrib' => {'name' => 'TESTING AND TIMING'},'name' => 'commandgroup'},{'type' => 'e','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! MAIN INITIAL AND BOUNDARY CONDITION PARAMETERS  !!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'type' => 'e','content' => [{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'value' => '1','name' => 'Uniform'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => '2','name' => 'Shock tube'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => '3','name' => 'Heliosphere'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => '5','name' => 'Comet'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => '6','name' => 'Rotation'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => '7','name' => 'Diffusion'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'default' => 'T','value' => '11','name' => 'Earth'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => '12','name' => 'Saturn'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => '13','name' => 'Jupiter'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => '14','name' => 'Venus'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => '21','name' => 'Cylinder'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => '22','name' => 'Sphere'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => '25','name' => 'Arcade'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => '26','name' => 'CME'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => '30','name' => 'Dissipation'},'name' => 'option'}],'attrib' => {'type' => 'integer','input' => 'select','name' => 'iProblem'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'if' => '$iProblem==30','type' => 'string','length' => '20','name' => 'TypeDissipation'},'name' => 'parameter'},{'type' => 't','content' => '
#PROBLEMTYPE
30			iProblem
heat_test1		TypeProblemDiss

! select a problem type which defines defaults for a lot of parameters
!
! Problem type has to be defined as the first item after #TEST..#DEBUG items!
!
!                           iProblem: 1=MHD Uniform Flow
!                                     2=Shock tube
!                                     3=Solar Wind and Inner Heliosphere
!                                     5=Mass-Loaded Comet
!                                     6=Rotation test
!                                     7=Diffusion test
!                                    11=Earth Magnetosphere
!                                    12=Saturn Magnetosphere
!                                    13=Jupiter Magnetosphere
!                                    14=Venus Ionosphere
!                                    21=Conducting Cylinder (2-D)
!                                    22=Conducting Sphere   (3-D)
!                                    25=Arcade
!                                    26=CME
!				     30=Test Dissipative MHD
'}],'attrib' => {'if' => '$_IsFirstSession','required' => 'T','name' => 'PROBLEMTYPE'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'if' => '$_NameComp eq \'GM\'','default' => 'T','value' => 'GSM','name' => 'GeoSolarMagnetic, GSM'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'if' => '$_NameComp ne \'GM\'','default' => 'T','value' => 'HGI','name' => 'HelioGraphicInertial, HGI'},'name' => 'option'}],'attrib' => {'type' => 'string','input' => 'select','name' => 'TypeCoordSystem'},'name' => 'parameter'},{'type' => 't','content' => '

#COORDSYSTEM
GSM			TypeCoordSystem

! TypeCoordSystem defines the coordinate system for the component.
! Currently only one coordinate system is available for GM ("GSM")
! and one for IH or SC ("HGI"). In the near future "GSE" should be also
! an option for GM.
!
! Default is component dependent: "GSM" for GM and "HGI" for IH or SC.
'}],'attrib' => {'if' => '$_IsFirstSession','name' => 'COORDSYSTEM'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'string','default' => 'GM/restartIN','length' => '100','name' => 'NameRestartInDir'},'name' => 'parameter'},{'type' => 't','content' => '

#RESTARTINDIR
GM/restart_n5000	NameRestartInDir

! The NameRestartInDir variable contains the name of the directory
! where restart files are saved relative to the run directory.
! The directory should be inside the subdirectory with the name 
! of the component.
!
! Default value is "GM/restartIN".
'}],'attrib' => {'if' => '$_IsFirstSession','name' => 'RESTARTINDIR'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'integer','default' => '2','min' => '1','name' => 'nRootX'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'integer','default' => '1','min' => '1','name' => 'nRootY'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'integer','default' => '1','min' => '1','name' => 'nRootZ'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '-192.0','name' => 'xMin'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '  64.0','min' => '$xMin','name' => 'xMax'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => ' -64.0','name' => 'yMin'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '  64.0','min' => '$yMin','name' => 'yMax'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => ' -64.0','name' => 'zMin'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '  64.0','min' => '$zMin','name' => 'zMax'},'name' => 'parameter'},{'type' => 't','content' => '
#GRID
2                       nIRoot_D(1)
1                       nJRoot_D(2)
1                       nKRoot_D(3)
-224.                   xMinALL
 32.                    xMaxALL
-64.                    yMinALL
 64.                    yMaxALL
-64.                    zMinALL
 64.                    zMaxALL

! Grid size should always be set.
! nRootX, nRootY, nRootZ define the number of blocks of the base grid, ie.
! the roots of the octree. Each root block must be on a differenet PE.
'}],'attrib' => {'if' => '$_IsFirstSession','required' => 'T','name' => 'GRID'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'name' => 'coupled'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'default' => '$Side ne \'TypeBcEast\'','name' => 'fixed/inflow'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'default' => '$Side eq \'TypeBcEast\'','name' => 'float/outflow'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'heliofloat'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'reflect'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'periodic'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'vary'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'shear'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'linetied'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'raeder'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'arcadetop'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'arcadebot'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'arcadebotcont'},'name' => 'option'}],'attrib' => {'type' => 'string','input' => 'select','name' => '$Side'},'name' => 'parameter'}],'attrib' => {'values' => 'TypeBcEast,TypeBcWest,TypeBcSouth,TypeBcNorth,TypeBcBot,TypeBcTop','name' => 'Side'},'name' => 'foreach'},{'type' => 'e','content' => [{'type' => 't','content' => '
	! East and west BCs must be both periodic or neither
	'}],'attrib' => {'expr' => 'not($TypeBcEast eq \'periodic\' xor $TypeBcWest eq \'periodic\')'},'name' => 'rule'},{'type' => 'e','content' => [{'type' => 't','content' => '
	! South and North BCs must be both periodic or neither
	'}],'attrib' => {'expr' => 'not($TypeBcSouth eq \'periodic\' xor $TypeBcNorth eq \'periodic\')'},'name' => 'rule'},{'type' => 'e','content' => [{'type' => 't','content' => '
	! Bottom and top BCs must be both periodic or neither
	'}],'attrib' => {'expr' => 'not($TypeBcBot eq \'periodic\' xor $TypeBcTop eq \'periodic\')'},'name' => 'rule'},{'type' => 't','content' => '
#OUTERBOUNDARY
outflow                 TypeBcOuter_E(East_)
inflow                  TypeBcOuter_E(West_)
float                   TypeBcOuter_E(South_)
float                   TypeBcOuter_E(North_)
float                   TypeBcOuter_E(Bot_)
float                   TypeBcOuter_E(Top_)

! Default depends on problem type.
! Possible values:
! fixed/inflow  - fixed solarwind values
! fixedB1       - fixed solarwind values without correction for the dipole B0
! float/outflow - zero gradient
! linetied      - float P, rho, and B, reflect all components of U
! raeder        - Jimmy Raeder\'s BC
! reflect       - reflective
! periodic      - periodic
! vary          - time dependent BC (same as fixed for non time_accurate)
! shear         - sheared (intended for shock tube problem only)
! arcadetop     - intended for arcade problem only
! arcadebot     - intended for arcade problem only
! arcadebotcont - intended for arcade problem only
'}],'attrib' => {'name' => 'OUTERBOUNDARY'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 't','content' => '
! Inner boundary types for body 1 and body 2
	'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'name' => 'reflect'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'float'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'fixed'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'default' => 'T','name' => 'ionosphere'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'ionosphereB0/ionosphereb0','name' => 'ionosphereB0'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'ionospherefloat'},'name' => 'option'}],'attrib' => {'type' => 'string','input' => 'select','name' => 'TypeInnerBc'},'name' => 'parameter'},{'type' => 't','content' => '
#INNERBOUNDARY
ionosphere              InnerBCType

! Default is ionosphere for Earth, Saturn, Jupiter, and problem_rotation.
! For all other problems with an inner boundary the default is \'reflect\'.
! If UseIonosphere=.true., velocity is determined by the coupled ionosphere
! model.
!
! Possible values for TypeBcInner are
!
! \'reflect\'     - reflect Vr, reflect Vphi to rotation, float Vtheta,
!                 reflect Br, float Bphi, float Btheta, float rho, float P
! \'float\'       - float Vr, reflect Vphi to rotation, float Vtheta,
!                 float B, float rho, float P
! \'fixed\'       - Vr=0, Vphi=rotation, Vtheta=0
!                 B=B0 (ie B1=0), fix rho, fix P
! \'ionosphere\'  - set V as if ionosphere gave V_iono=0
!                 float B, fix rho, fix P
! \'ionospherefloat\'-set V as if ionosphere gave V_iono=0
!                 float B, float rho, float P
'}],'attrib' => {'name' => 'INNERBOUNDARY'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'F','name' => 'UseExtraBoundary'},'name' => 'parameter'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'string','name' => 'TypeBcExtra'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'F','name' => 'DoFixExtraboundary'},'name' => 'parameter'}],'attrib' => {'expr' => '$UseExtraBoundary'},'name' => 'if'}],'attrib' => {'name' => 'EXTRABOUNDARY'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'integer','default' => '0','max' => '6','min' => '0','name' => 'MaxBoundary'},'name' => 'parameter'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'F','name' => 'DoFixOuterBoundary'},'name' => 'parameter'}],'attrib' => {'expr' => '$MaxBoundary >= 1'},'name' => 'if'},{'type' => 't','content' => '
#FACEOUTERBC
0              MaxBoundary            
F              DoFixOuterBoundary)    !read only for MaxBoundary>=East_(=1).
! if MaxBoundary>=East_(=1) then the outer boundaries with
! the number of boundary being between East_ and MaxBoundary
! are treated using set_BCs.f90 subroutines instead of set_outerBCs.f90 
! if DoFixOuterBoundary==.true., there is no resolution
! change along the outer boundaries with the number of
! of boundary being between East_ and MaxBoundary
'}],'attrib' => {'name' => 'FACEOUTERBC'},'name' => 'command'}],'attrib' => {'name' => 'INITIAL AND BOUNDARY CONDITIONS'},'name' => 'commandgroup'},{'type' => 'e','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! INITIAL TIME AND STEP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'integer','default' => '2000','name' => 'year'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'integer','default' => '3','max' => '12','min' => '1','name' => 'month'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'integer','default' => '21','max' => '31','min' => '1','name' => 'day'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'integer','default' => '0','max' => '23','min' => '0','name' => 'hour'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'integer','default' => '0','max' => '59','min' => '0','name' => 'minute'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'integer','default' => '0','max' => '59','min' => '0','name' => 'second'},'name' => 'parameter'},{'type' => 't','content' => '
#STARTTIME
2000                    StartTime_i(1)=year
3                       StartTime_i(2)=month
21                      StartTime_i(3)=day
10                      StartTime_i(4)=hour
45                      StartTime_i(5)=minute
0                       StartTime_i(6)=second

The #STARTTIME command sets the initial date and time for the
simulation in Greenwich Mean Time (GMT) or Universal Time (UT)
in stand alone mode. 
In the SWMF this command checks start times against the SWMF start time 
and warns if the difference exceeds 1 millisecond.
This time is stored in the BATSRUS restart header file.

The default values are shown above.
This is a date and time when both the rotational and the magnetic axes
have approximately zero tilt towards the Sun.
'}],'attrib' => {'alias' => 'SETREALTIME','if' => '$_IsFirstSession','name' => 'STARTTIME'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '0.0','min' => '0','name' => 'tSimulation'},'name' => 'parameter'},{'type' => 't','content' => '

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
'}],'attrib' => {'if' => '$_IsFirstSession','name' => 'TIMESIMULATION'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'integer','default' => '0','min' => '0','name' => 'nStep'},'name' => 'parameter'},{'type' => 't','content' => '

#NSTEP
100			nStep

! Set nStep for the component. Typically used in the restart.H header file.
! Generally it is not inserted in a PARAM.in file by the user.
!
! The default is nStep=0 as the starting time step with no restart.
'}],'attrib' => {'if' => '$_IsFirstSession','name' => 'NSTEP'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'integer','default' => '-1','min' => '-1','name' => 'nPrevious'},'name' => 'parameter'},{'type' => 't','content' => '

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
'}],'attrib' => {'if' => '$_IsFirstSession','name' => 'NPREVIOUS'},'name' => 'command'}],'attrib' => {'name' => 'INITIAL TIME'},'name' => 'commandgroup'},{'type' => 'e','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  TIME INTEGRATION PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'type' => 'e','content' => [{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'default' => 'T','value' => '1'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => '2'},'name' => 'option'}],'attrib' => {'type' => 'integer','input' => 'select','name' => 'nStage'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '0.8','max' => '1','min' => '0','name' => 'CflExpl'},'name' => 'parameter'},{'type' => 't','content' => '

#TIMESTEPPING
2                       nStage
0.80                    CflExpl

! Parameters for explicit time integration.
! Default is 1 stage and CflExpl=0.8
'}],'attrib' => {'name' => 'TIMESTEPPING'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'F','name' => 'UseDtFixed'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'if' => '$UseDtFixed','type' => 'real','default' => '1.0','min' => '0','name' => 'DtFixedDim'},'name' => 'parameter'},{'type' => 't','content' => '
#FIXEDTIMESTEP
T                       UseDtFixed
10.                     DtFixedDim [sec] (read if UseDtFixed is true)

! Default is UseDtFixed=.false. Effective only if DoTimeAccurate is true.
! If UseDtFixed is true, the time step is fixed to DtFixedDim.
!
! This is useful for debugging explicit schemes.
'}],'attrib' => {'name' => 'FIXEDTIMESTEP'},'name' => 'command'}],'attrib' => {'name' => 'TIME INTEGRATION'},'name' => 'commandgroup'},{'type' => 'e','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! STOPPING CRITERIA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

The commands in this group only work in stand alone mode.

'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'integer','default' => '-1','min' => '-1','name' => 'MaxIteration'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '-1','min' => '-1','name' => 'tSimulationMax'},'name' => 'parameter'},{'type' => 't','content' => '

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
'}],'attrib' => {'if' => '$_IsStandAlone','required' => '$_IsStandAlone','name' => 'STOP'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'T','name' => 'DoCheckStopFile'},'name' => 'parameter'},{'type' => 't','content' => '

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
'}],'attrib' => {'if' => '$_IsStandAlone','name' => 'CHECKSTOPFILE'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '-1','min' => '-1','name' => 'CpuTimeMax'},'name' => 'parameter'},{'type' => 't','content' => '

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
'}],'attrib' => {'if' => '$_IsStandAlone','name' => 'CPUTIMEMAX'},'name' => 'command'}],'attrib' => {'name' => 'STOPPING CRITERIA'},'name' => 'commandgroup'},{'type' => 'e','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  OUTPUT PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'string','default' => 'GM/restartOUT','length' => '100','name' => 'NameRestartOutDir'},'name' => 'parameter'},{'type' => 't','content' => '

#RESTARTOUTDIR
GM/restart_n5000	NameRestartOutDir

! The NameRestartOutDir variable contains the name of the directory
! where restart files are saved relative to the run directory.
! The directory should be inside the subdirectory with the name 
! of the component.
!
! Default value is "GM/restartOUT".
'}],'attrib' => {'name' => 'RESTARTOUTDIR'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'T','name' => 'SaveRestart'},'name' => 'parameter'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'integer','default' => '-1','min' => '-1','name' => 'DnRestart'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '-1','min' => '-1','name' => 'DtRestart'},'name' => 'parameter'}],'attrib' => {'expr' => '$SaveRestart'},'name' => 'if'},{'type' => 't','content' => '
#SAVERESTART
T			saveRestartFile  Rest of parameters read if true
100			DnOutput_i(restart_)
-1.			DtOutput_i(restart_) in seconds. Read if time_accurate!

! Default is save_restartfile=.true. with DnOutput(restart_)=-1, 
! DtOutput(restart_)=-1. This results in the restart file being 
! saved only at the end.  A binary restart file is produced for every 
! block and named as
!
! restartOUT/blkGLOBALBLKNUMBER.rst
!
! In addition the grid is described by
!
! restartOUT/octree.rst
!
! and an ASCII header file is produced with timestep and time info:
!
! restartOUT/restart.H
!
! The restart files are overwritten every time a new restart is done.
'}],'attrib' => {'name' => 'SAVERESTART'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'string','default' => 'GM/IO2','length' => '100','name' => 'NamePlotDir'},'name' => 'parameter'},{'type' => 't','content' => '

The NamePlotDir variable contains the name of the directory
where plot files and logfiles are saved relative to the run directory.
The directory should be inside the subdirectory with the name
of the component.

Default value is "GM/IO2".
'}],'attrib' => {'name' => 'PLOTDIR'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'F','name' => 'DoSaveLogfile'},'name' => 'parameter'},{'type' => 'e','content' => [{'type' => 'e','content' => [{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'default' => 'T','value' => 'MHD','name' => 'MHD vars. dimensional'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'FLX','name' => 'Flux vars. dimensional'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'RAW','name' => 'Raw vars. dimensional'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'VAR','name' => 'Set vars. dimensional'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'default' => 'T','value' => 'mhd','name' => 'MHD vars. scaled'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'flx','name' => 'Flux vars. scaled'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'raw','name' => 'Raw vars. scaled'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'var','name' => 'Set vars. scaled'},'name' => 'option'}],'attrib' => {'type' => 'string','required' => 'T','input' => 'select','name' => 'TypeLogVar'},'name' => 'part'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'exclusive' => 'T','name' => 'none'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'step'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'date'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'time'},'name' => 'option'}],'attrib' => {'multiple' => 'T','type' => 'string','required' => 'F','input' => 'select','name' => 'TypeTime'},'name' => 'part'}],'attrib' => {'type' => 'strings','max' => '4','min' => '1','name' => 'StringLog'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'integer','default' => '1','min' => '-1','name' => 'DnOutput'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '-1','min' => '-1','name' => 'DtOutput'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'if' => '$TypeLogVar =~ /var/i','type' => 'string','length' => '100','name' => 'NameLogVars'},'name' => 'parameter'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'multiple' => 'T','type' => 'real','min' => '$rBody','name' => 'LogRadii'},'name' => 'part'}],'attrib' => {'if' => '($TypeLogVar=~/flx/i or $NameLogVars=~/flx/i)','type' => 'strings','length' => '100','max' => '10','min' => '1','name' => 'StringLogRadii'},'name' => 'parameter'}],'attrib' => {'expr' => '$DoSaveLogfile'},'name' => 'if'},{'type' => 't','content' => '
#SAVELOGFILE
T                       DoSaveLogfile, rest of parameters read if true
VAR step date           StringLog
100                     DnOutput_i(logfile_)
-1.                     DtOutput_i(logfile_) in sec. Read only if time accurate
rho p rhoflx            NameLogVars (variable to write) Read for \'var\' or \'VAR\'
4.0  10.0               rLog  !radii where flx is calc. Read if vars inc. flx.

! Default is save_logfile=.false.
! The logfile can contain averages or point values and other scalar
! quantities.  It is written into an ASCII file named as
!
! IO2/log_timestep.log
!
! The StringLog can contain two groups of information in arbitrary order.
! The first is LogVar which is a single 3 character string that indicates
! the type of variables that are to be writen.  The second group indicates
! the type of time/iteration output format to use.  This second group is
! not required and defaults to something standard for each logvar case.
! Any of the identifiers for the timetype can be includec in arbitrary order.
!
! logvar  = \'mhd\', \'raw\', \'flx\' or \'var\' - unitless output
! logvar  = \'MHD\', \'RAW\', \'FLX\' or \'VAR\' - dimensional output
! timetype = \'none\', \'step\', \'time\', \'date\'
!
! The logvar string is not optional and must be found on the line.
! The timetype is optional - when not specified a logical choice is made
!       by the code
!
! The log_var string defines the variables to print in the log file
! It also controls whether or not the variables will come out in
! dimensional or non-dimensional form by the capatilization of the log_var
! string.
!
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
!
! log_vars is read only when the log_string contains var or VAR.  The choices
! for variables are currently:
!
! Average value on grid: rho rhoUx rhoUy rhoUz Ux Uy Uz Bx By Bz P E
! Value at the test point: rhopnt rhoUxpnt rhoUypnt rhoUxpnt Uxpnt Uypnt Uzpnt
!                          Bxpnt Bypnt Bzpnt B1xpnt B1ypnt B1zpnt
!                          Epnt Ppnt Jxpnt Jypnt Jzpnt
!                          theta1pnt theta2pnt phi1pnt phi2pnt statuspnt
!
! Max or Min on grid:    Pmin Pmax
! Flux values:           Aflx rhoflx Bflx B2flx pvecflx e2dflx
! Other variables:     dt
!
! timetype values mean the following:
!  none  = there will be no indication of time in the logfile (not even an
!                # of steps)
!  step  = # of time steps (n_steps)
!  date  = time is given as an array of 7 integers:  year mo dy hr mn sc msc
!  time  = time is given as a real number - elapsed time since the start of
!          the run.  Units are determined by log_var and unitUSER_t
!
!  these can be listed in any combination in the log_string line
!
! R_log is read only when one of the variables used is a \'flx\' variable.  R_log
! is a list of radii at which to calculate the flux through a sphere.
'}],'attrib' => {'name' => 'SAVELOGFILE'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'integer','default' => '0','min' => '0','name' => 'nSatellite'},'name' => 'parameter'},{'type' => 'e','content' => [{'type' => 'e','content' => [{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'default' => 'T','value' => 'MHD','name' => 'MHD vars. dimensional'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'FUL','name' => 'All vars. dimensional'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'VAR','name' => 'Set vars. dimensional'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'mhd','name' => 'MHD vars. scaled'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'ful','name' => 'All vars. scaled'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'var','name' => 'Set vars. scaled'},'name' => 'option'}],'attrib' => {'type' => 'string','required' => 'T','input' => 'select','name' => 'TypeSatelliteVar'},'name' => 'part'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'default' => 'T','name' => 'file'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'eqn','name' => 'equation'},'name' => 'option'}],'attrib' => {'type' => 'string','required' => 'F','input' => 'select','name' => 'TypeTrajectory'},'name' => 'part'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'exclusive' => 'T','name' => 'none'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'step'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'date'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'time'},'name' => 'option'}],'attrib' => {'multiple' => 'T','type' => 'string','required' => 'F','input' => 'select','name' => 'TypeTime'},'name' => 'part'}],'attrib' => {'type' => 'strings','max' => '5','min' => '1','name' => 'StringSatellite'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'integer','default' => '1','min' => '-1','name' => 'DnOutput'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '-1','min' => '-1','name' => 'DtOutput'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'string','length' => '100','name' => 'NameSatellite'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'if' => '$TypeSatelliteVar =~ /\\bvar\\b/i','type' => 'string','length' => '100','name' => 'NameSatelliteVars'},'name' => 'parameter'}],'attrib' => {'to' => '$nSatellite','from' => '1'},'name' => 'for'},{'type' => 't','content' => '
#SATELLITE
2                       nSatellite
MHD file                StringSatellite (variables and traj type)
100                     DnOutput_i(satellite_)
-1.                     DtOutput_i(satellite_) in sec. ALWAYS READ!
satellite1.dat          Filename or satellite name (Satellite_name(satellite_))
VAR eqn step date       StringSatellite
100                     DnOutput_i(satellite_)
-1.                     DtOutput_i(satellite_) in sec. ALWAYS READ!
satellite2.dat          NameSatellite_i(satellite_)
rho p                   NameSatelliteVars Read if satellitevar=\'var\' or \'VAR\'

! satellite_string can contain the following 3 parts in arbitrary order
!
! satellitevar  = \'mhd\', \'ful\' or \'var\' - unitless output
! satellitevar  = \'MHD\', \'FUL\' or \'VAR\' - dimensional output
! trajectory_type = \'file\' or \'eqn\'
! timetype = \'none\', \'step\', \'time\', \'date\'
!
! satellitevar -> REQUIRED
! trajectory_type -> not required - defaults to \'file\'
! time_type -> not required - a logical default is used
!
! The satellitevar string defines the variables to print in the satellite
! output file.  It also controls whether or not the variables will come out in
! dimensional or non-dimensional form by the capatilization of the
! satellite_vars string.
!
! ALL CAPS  - dimensional
! all lower - dimensionless
!
! \'mhd\' - vars: rho Ux Uy Uz Bx By Bz P Jx Jy Jz
! \'ful\' - vars: rho Ux Uy Uz Bx By Bz P Jx Jy Jz theta1 phi1 theta2 phi2 status
! \'var\' - vars: READ FROM PARAMETER FILE
!
! satellite_vars is read only when the satellite_string is var or VAR.  The
! choices for variables are currently:
!
! rho, rho, rhouy, rhouz, ux, uy, uz
! Bx, By, Bz, B1x, B1y, B1z
! E, P, Jx, Jy, Jz
! theta1,theta2,phi1,phi2,status
!
!
! timetype values mean the following:
!  none  = there will be no indication of time in the logfile (not even an
!                # of steps)
!  step  = # of time steps (n_steps)
!  date  = time is given as an array of 7 integers:  year mo dy hr mn sc msc
!  time  = time is given as a real number - elapsed time since the start of
!          the run.  Units are determined by satellitevar and unitUSER_t
!
!  More than one of these can be listed.  They can be put together in any
!  combination.
'}],'attrib' => {'if' => '$_IsFirstSession','name' => 'SATELLITE'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 't','content' => '
! plot_string must contain the following 3 parts in arbitrary order
...
	'},{'type' => 'e','content' => [],'attrib' => {'type' => 'integer','default' => '0','max' => '100','min' => '0','name' => 'nPlotFile'},'name' => 'parameter'},{'type' => 'e','content' => [{'type' => 'e','content' => [{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'name' => 'TECPLOT','value' => 'tec'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'IDL','value' => 'idl'},'name' => 'option'}],'attrib' => {'type' => 'string','required' => 'T','input' => 'select','name' => 'plotform'},'name' => 'part'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'value' => '3d'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'x=0'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'default' => 'T','value' => 'y=0'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'z=0'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'sph'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'los'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'lin'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'if' => '$plotform =~ /\\bidl\\b/','value' => 'cut'},'name' => 'option'}],'attrib' => {'type' => 'string','required' => 'T','input' => 'select','name' => 'plotarea'},'name' => 'part'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'value' => 'MHD','name' => 'MHD vars. dimensional'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'FUL','name' => 'All vars. dimensional'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'RAW','name' => 'Raw vars. dimensional'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'RAY','name' => 'Ray tracing vars. dim.'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'FLX','name' => 'Flux vars. dimensional'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'SOL','name' => 'Solar vars. dimensional'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'VAR','name' => 'Select dimensional vars.'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'mhd','name' => 'MHD vars. scaled'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'ful','name' => 'All vars. scaled'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'raw','name' => 'Raw vars. scaled'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'ray','name' => 'Ray tracing vars. scaled'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'flx','name' => 'Flux vars. scaled'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'sol','name' => 'Solar vars. scaled'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'var','name' => 'Select scaled vars.'},'name' => 'option'}],'attrib' => {'type' => 'string','required' => 'T','input' => 'select','name' => 'plotvar'},'name' => 'part'}],'attrib' => {'type' => 'strings','max' => '3','min' => '3','name' => 'plotString'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'integer','min' => '-1','name' => 'DnOutput'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','min' => '-1','name' => 'DtOutput'},'name' => 'parameter'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'real','name' => 'xMinCut'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','min' => '$xMinCut','name' => 'xMaxCut'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','name' => 'yMinCut'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','min' => '$yMinCut','name' => 'yMaxCut'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','name' => 'zMinCut'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','min' => '$zMinCut','name' => 'zMaxCut'},'name' => 'parameter'}],'attrib' => {'expr' => '$plotarea =~ /\\bcut\\b/'},'name' => 'if'},{'type' => 'e','content' => [],'attrib' => {'if' => '$plotarea =~ /\\bsph\\b/','type' => 'real','default' => '10','min' => '0','name' => 'radius'},'name' => 'parameter'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '0','name' => 'LosVectorX'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '0.0001','name' => 'LosVectorY'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '1','name' => 'LosVectorZ'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '20','min' => '0','name' => 'xSizeImage'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '20','min' => '0','name' => 'ySizeImage'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '10','name' => 'xOffset'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '10','name' => 'yOffset'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '2.5','min' => '1','name' => 'rOccult'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '0.5','max' => '1','min' => '0','name' => 'MuLimbDarkening'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'integer','default' => '200','min' => '2','name' => 'nPixX'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'integer','default' => '200','min' => '2','name' => 'nPixY'},'name' => 'parameter'}],'attrib' => {'expr' => '$plotarea =~ /\\blos\\b/'},'name' => 'if'},{'type' => 'e','content' => [{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'value' => 'A','name' => 'Advected B'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'default' => 'T','name' => 'B'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'U'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'J'},'name' => 'option'}],'attrib' => {'type' => 'string','input' => 'select','name' => 'NameLine'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'F','name' => 'IsSingleLine'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'integer','default' => '1','max' => '20','min' => '1','name' => 'nLine'},'name' => 'parameter'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'real','name' => 'xStartLine'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','name' => 'yStartLine'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','name' => 'zStartLine'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','name' => 'IsParallel'},'name' => 'parameter'}],'attrib' => {'to' => '$nLine','from' => '1'},'name' => 'for'}],'attrib' => {'expr' => '$plotarea =~ /\\blin\\b/'},'name' => 'if'},{'type' => 'e','content' => [],'attrib' => {'if' => '($plotform=~/\\bidl\\b/ and $plotarea!~/\\b(sph|los|lin)\\b/)','type' => 'real','default' => '-1.0','min' => '-1.0','name' => 'DxPlot'},'name' => 'parameter'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'string','length' => '100','name' => 'plotVars'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'string','length' => '100','name' => 'plotPars'},'name' => 'parameter'}],'attrib' => {'expr' => '$plotvar =~ /\\bvar\\b/i'},'name' => 'if'}],'attrib' => {'to' => '$nPlotFile','from' => '1'},'name' => 'for'},{'type' => 't','content' => '
#SAVEPLOT
6			nPlotfile
3d MHD tec		StringPlot ! 3d MHD data
100			DnOutput
-1.			DtOutput
y=0 VAR idl		plotString ! y=0 cut
-1			DnOutput
100.			DtOutput
2.			PlotDx     ! Read only for format \'idl\'
jx jy jz		PlotVars   ! Read only for content \'var\'
g unitx unitv unitn	PlotPars   ! Read only for content \'var\'
cut ray idl		StringPlot ! ray tracing plot
1			DnOutput
-1.			DtOutput
-10.			PlotRange_ei(x1,3) Read only for area \'cut\'
10.			plotRange_ei(x2,3) Read only for area \'cut\'
-10.			plotRange_ei(y1,3) Read only for area \'cut\'
10.			plotRange_ei(y2,3) Read only for area \'cut\'
-10.			plotRange_ei(z1,3) Read only for area \'cut\'
10.			plotRange_ei(z2,3) Read only for area \'cut\'
1.			plotDx      ! Read only for format \'idl\'
sph flx idl		plotString  ! spherical plot
-1			DnOutput
100.			DtOutput
4.			rPlot - R of spherical cut, Read only for area \'sph\'
los sol idl             StringPlot  ! line of sight plot
-1			DnOutput
100.			DtOutput
1.			LosVector_i(1)
0.			LosVector_i(2)
0.			LosVector_i(3)
30.			xSizeImage
50.			ySizeImage
10.			xOffset
20.			yOffset
5.			rOccult
0.5			MuLimbDarkening
256			nPixX
256			nPixY
lin mhd idl		PlotString  ! field line plot
-1			DnOutput
10.			DtOutput
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

! Default is nplotfile=0
! plot_string must contain the following 3 parts in arbitrary order
!
! plotarea plotvar plotform
!
! plotarea = \'3d\' , \'x=0\', \'y=0\', \'z=0\', \'cut\', \'sph\', \'los\', \'lin\'
! plotvar  = \'mhd\', \'ful\',\'raw\', \'ray\', \'flx\', \'sol\', \'var\' - unitless output
! plotvar  = \'MHD\', \'FUL\',\'RAW\', \'RAY\', \'FLX\', \'SOL\', \'VAR\' - dimensional
! plotform = \'tec\', \'idl\'
!
! NOTES: The plotvar option \'sol\' is only valid for plotarea \'los\'.
!
! The plotarea string defines the 1, 2, or 3D volume of the plotting area:
!
! x=0	- full x=0 plane: xmin=-0.001, xmax=0.001, average for symmetry plane
! y=0	- full y=0 plane: ymin=-0.001, ymax=0.001, average for symmetry plane
! z=0	- full z=0 plane: zmin=-0.001, zmax=0.001, average for symmetry plane
! 3d	- full 3D volume
! cut	- READ PLOTRANGE FROM PARAM.in, only works for plotform=\'idl\'
! sph   - spherical cut at radius R_plot, READ FROM PARAM.in
! los   - line of sight integrated plot
! lin   - one dimensional plot along a field or stream line
!
! The plotvar string defines the plot variables and the equation parameters.
! It also controls whether or not the variables will be plotted in dimensional
! values or as non-dimensional values:
!
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
!
! The plot_string is always followed by the plotting frequency
! DnOutput and for time accurate runs by DtOutput.
!
! Depending on plot_string, further information is read from the parameter file
! in this order:
!
! PlotRange		if plotarea is \'cut\'
! DxPlot		if plotform is \'idl\' and plotarea is not sph, ion, los
! rPlot			if plotarea is \'sph\'
! plotVars		if plotform is \'var\'
! plotPars		if plotform is \'var\'
!
! The plot_range is described by 6 coordinates. If the width in one or two 
! dimensions is less than the smallest cell size within the plotarea, 
! then the plot file will be 2 or 1 dimensional. If the range is thin but
! symmetric about one of the x=0, y=0, or z=0 planes, data will be averaged
! in the postprocessing.
!
! Possible values for plotDx (for IDL files):
!
!  0.5	- fixed resolution (any positive value)
!  0.	- fixed resolution based on the smallest cell in the plotting area
! -1.	- unstructured grid will be produced by PostIDL.exe
!
! rPlot is the radius of the spherical cut for plotarea=\'sph\'
!
! LosVector_i defines the direction of the line of sight integration
! xSizeImage, ySizeImage defines the size of the LOS image
! xOffset, yOffset defines the offset relative to the origin (Sun)
! rOccult defines the minimum distance of the line from the origin (Sun)
! MuLimbDarkening is the limb darkening parameter for the \'wl\' (white light)
!                 and \'pb\' (polarization brightness) plot variables.
!
! The possible values for plot_vars with plotarea \'los\' 
!       are listed in subroutine set_plotvar_los in write_plot_los.f90.
! The possible values for plot_vars for other plot areas
!       are listed in subroutine set_plotvar in write_plot_common.f90.
!
! The possible values for plot_pars 
!       are listed in subroutine set_eqpar in write_plot_common.f90
!
! A plot file is produced by each processor.  This file is ASCII in \'tec\'
! format and can be either binary or ASCII in \'idl\' format as chosen under
! the #SAVEBINARY flag.  The name of the files are
!
! IO2/plotarea_plotvar_plotnumber_timestep_PEnumber.extenstion 
!
! where extension is \'tec\' for the TEC and \'idl\' for the IDL file formats.
! The plotnumber goes from 1 to nplot in the order of the files in PARAM.in.
! After all processors wrote their plot files, processor 0 writes a small 
! ASCII header file named as
!
! IO2/plotarea_plotvar_plotnumber_timestep.headextension
!
! where headextension is:
!           \'T\' for TEC file format
!           \'S\' for TEC and plot_area \'sph\' 
!           \'h\' for IDL file format       
!
! The line of sight integration produces TecPlot and IDL files directly:
!
! IO2/los_plotvar_plotnumber_timestep.extension
!
! where extension is \'dat\' for TecPlot and \'out\' for IDL file formats.
! The IDL output from line of sight integration is always in ASCII format.

'}],'attrib' => {'name' => 'SAVEPLOT'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'T','name' => 'DoSaveBinary'},'name' => 'parameter'},{'type' => 't','content' => '
#SAVEBINARY
T			DoSaveBinary   used only for \'idl\' plot file

! Default is .true. Saves unformatted IO2/*.idl files if true. 
! This is the recommended method, because it is fast and accurate.
! The only advantage of saving IO2/*.idl in formatted text files is
! that it can be processed on another machine or with a different 
! (lower) precision. For example PostIDL.exe may be compiled with 
! single precision to make IO2/*.out files smaller, while BATSRUS.exe is 
! compiled in double precision, to make results more accurate.
'}],'attrib' => {'name' => 'SAVEBINARY'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'F','name' => 'SavePlotsAmr'},'name' => 'parameter'},{'type' => 't','content' => '
#SAVEPLOTSAMR
F			savePlotsAMR to save plots before each AMR

! Default is save_plots_amr=.false.
'}],'attrib' => {'name' => 'SAVEPLOTSAMR'},'name' => 'command'}],'attrib' => {'name' => 'OUTPUT PARAMETERS'},'name' => 'commandgroup'},{'type' => 'e','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  AMR PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'type' => 'e','content' => [{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'default' => '1','name' => 'default'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'all'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'none'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => '3Dbodyfocus'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'spherefocus'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'magnetosphere'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'points'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'helio_init'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'helio_z=4'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'all_then_focus'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'cme'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'points'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'mag_new'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'magnetosphere'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'magneto_fine'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'magneto12'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'magnetosaturn'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'magnetojupiter'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'paleo'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'comet'},'name' => 'option'}],'attrib' => {'type' => 'string','input' => 'select','name' => 'InitialRefineType'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'integer','default' => '4','min' => '0','name' => 'InitialRefineLevel'},'name' => 'parameter'},{'type' => 't','content' => '
#AMRINIT
default			InitialRefineType
4			InitialRefineLevel

! These are the default values for the initial refinement.

! Possible values for InitialRefineType:
! Default depends on problem_type. 
! \'none\'		- Refine no blocks
! \'all\' 		- Refine all blocks
! \'3Dbodyfocus\'		- Refinement focusing on body
! \'spherefocus\'		- Refinement focusing on the orgin, does not require 
!                           a body
! \'points\'      	- Refine around given points
! \'magnetosphere\'	- Refine for generic magnetosphere
! *			- any other value will use default value by ProblemType
'}],'attrib' => {'if' => '$_IsFirstSession','name' => 'AMRINIT'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'integer','default' => '0','min' => '0','name' => 'nRefineLevelIC'},'name' => 'parameter'},{'type' => 't','content' => '
#AMRINITPHYSICS
3			nRefineLevelIC

! Defines number of physics (initial condition) based AMR-s AFTER the 
! geometry based initial AMR-s defined by #AMRINIT were done.
! Only useful if the initial condition has a non-trivial analytic form.
'}],'attrib' => {'if' => '$_IsFirstSession','name' => 'AMRINITPHYSICS'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'integer','default' => '0','min' => '-1','name' => 'minBlockLevel'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'integer','default' => '99','min' => '-1','name' => 'maxBlockLevel'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'F','name' => 'FixBodyLevel'},'name' => 'parameter'},{'type' => 't','content' => '
#AMRLEVELS
0			minBlockLevel
99			maxBlockLevel
F			fixBodyLevel

! Set the minimum/maximum levels that can be affected by AMR.  The usage is as
! follows:
!
! minBlockLevel .ge.0 Cells can be coarsened up to the listed level but not
!                       further.
! minBlockLevel .lt.0 The current grid is ``frozen\'\' for coarsening such that
!                       blocks are not allowed to be coarsened to a size
!                       larger than their current one.
! maxBlockLevel .ge.0 Any cell at a level greater than or equal to
!                       maxBlockLevel is uneffected by AMR (cannot be coarsened
!                       or refined).
! maxBlockLevel .lt.0 The current grid is ``frozen\'\' for refinement such that
!                       blocks are not allowed to be refined to a size
!                       smaller than their current one.
! fixBodyLevel = T    Blocks touching the body cannot be coarsened or refined.
!
! This command has no effect when automatic_refinement is .false.
!
! Note that the user can set either #AMRLEVELS or #AMRRESOLUTION but not
! both.  If both are set, the final one in the session will set the values
! for AMR.
'}],'attrib' => {'name' => 'AMRLEVELS'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '0','min' => '-1','name' => 'minCellDx'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '99999','min' => '-1','name' => 'maxCellDx'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'F','name' => 'FixBodyLevel'},'name' => 'parameter'},{'type' => 't','content' => '
#AMRRESOLUTION
0.			minCellDx
99999.			maxCellDx
F			fixBodyLevel

! Serves the same function as AMRLEVELS. min_block_dx and max_block_dx are
! converted into minBlockLevel and maxBlockLevel when they are read.
! Note that minBlockLevel corresponds to maxCellDx and maxBlockLevel
! corresponds to minCellDx.  See details above.
!
! This command has no effect when automatic_refinement is .false.
!
! Note that the user can set either #AMRLEVELS or #AMRRESOLUTION but not
! both.  If both are set, the final one in the session will set the values
! for AMR.
'}],'attrib' => {'name' => 'AMRRESOLUTION'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'integer','default' => '-1','min' => '-1','name' => 'DnRefine'},'name' => 'parameter'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'F','name' => 'DoAutoRefine'},'name' => 'parameter'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '20','max' => '100','min' => '0','name' => 'percentCoarsen'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '20','max' => '100','min' => '0','name' => 'percentRefine'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'integer','default' => '99999','min' => '1','name' => 'maxTotalBlocks'},'name' => 'parameter'}],'attrib' => {'expr' => '$DoAutoRefine'},'name' => 'if'}],'attrib' => {'expr' => '$DnRefine>0'},'name' => 'if'},{'type' => 't','content' => '
#AMR
2001			dnRefine (frequency in terms of total steps n_step)
T			DoAutoRefine 
0.			percentCoarsen
0.			percentRefine
99999			maxTotalBlocks

! Default for dn_refine is -1, ie. no run time refinement.
'}],'attrib' => {'name' => 'AMR'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'name' => '1'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => '2'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'default' => '1','name' => '3'},'name' => 'option'}],'attrib' => {'type' => 'integer','input' => 'select','name' => 'nRefineCrit'},'name' => 'parameter'},{'type' => 'e','content' => [{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'value' => 'gradt/gradT','name' => 'grad T'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'gradp/gradP','name' => 'grad P'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'gradlogrho','name' => 'grad log(Rho)'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'gradlogP/gradlogp','name' => 'grad log(p)'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'gradE','name' => 'grad E'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'curlV/curlv/curlU/curlu','name' => 'curl U'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'curlB/curlb','name' => 'curl B'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'divU/divu/divV/divv','name' => 'div U'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'divb/divB','name' => 'divB'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'Valfven/vAlfven/valfven','name' => 'vAlfven'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'heliobeta','name' => 'heliospheric beta'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'flux'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'heliocurrentsheet','name' => 'heliospheric current sheet'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'rcurrents/Rcurrents','name' => 'rCurrents'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'transient/Transient','name' => 'Transient'},'name' => 'option'}],'attrib' => {'type' => 'string','input' => 'select','name' => 'TypeRefine'},'name' => 'parameter'},{'type' => 'e','content' => [{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'value' => 'p_dot/P_dot','name' => 'P_dot'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 't_dot/T_dot','name' => 'T_dot'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'default' => 'T','value' => 'rho_dot/Rho_dot','name' => 'Rho_dot'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'RhoU_dot/rhou_dot','name' => 'RhoU_dot'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'Rho_2nd_1/rho_2nd_1','name' => 'Rho_2nd_1'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'Rho_2nd_2/rho_2nd_2','name' => 'Rho_2nd_2'},'name' => 'option'}],'attrib' => {'type' => 'string','input' => 'select','name' => 'TypeTransient'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'F','name' => 'UseSunEarth'},'name' => 'parameter'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'real','name' => 'xEarth'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','name' => 'yEarth'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','name' => 'zEarth'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','name' => 'InvD2Ray'},'name' => 'parameter'}],'attrib' => {'expr' => '$UseSunEarth'},'name' => 'if'}],'attrib' => {'expr' => '$TypeRefine =~ /transient/i'},'name' => 'if'}],'attrib' => {'to' => '$nRefineCrit','from' => '1'},'name' => 'for'},{'type' => 't','content' => '
#AMRCRITERIA
3			nRefineCrit (number of refinement criteria: 1,2 or 3)
gradlogP		RefineCrit_i(1)
divB			RefineCrit_i(2)
Transient		RefineCrit_i(3)
Rho_dot			TypeTransient_I(i) ! Only if \'Transient\' or \'transient\'
T			UseSunEarth 	   ! Only if \'Transient\'
0.00E+00		xEarth		   ! Only if UseSunEarth
2.56E+02 		yEarth		   ! Only if UseSunEarth
0.00E+00		zEarth		   ! Only if UseSunEarth
5.00E-01		InvD2Ray	   ! Only if UseSunEarth

! The default values depend on problem_type. 
! At most three criteria can be given. Possible criteria:
!
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
!
! All the names can also be spelled with all small case letters.
!
! The possible choices for TypeTransient_I 
!
! \'P_dot\' (same as \'p_dot\')
! \'T_dot\' (same as \'t_dot\')
! \'Rho_dot\' (same as \'rho_dot\')
! \'RhoU_dot\' (same as \'rhou_dot\')
! \'B_dot\' (same as \'b_dot\')
! \'Rho_2nd_1\' (same as \'rho_2nd_1\')
! \'Rho_2nd_2\' (same as \'rho_2nd_2\')
! 
! Also, (xEarth,yEarth,zEarth) are the coordinates of the Earth. InvD2Ray is
! a factor that defines how close to the ray Sun-Earth to refine the grid.
! Note that the AMR occurs in a cylinder around the ray.
! Example:: for InvD2Ray = 
!   1 - refine_profile = 0.3679 at distance Rsun/10 from the ray
!   2 - refine_profile = 0.0183 at distance Rsun/10 from the ray
!   3 - refine_profile = 0.0001 at distance Rsun/10 from the ray
'}],'attrib' => {'name' => 'AMRCRITERIA'},'name' => 'command'}],'attrib' => {'name' => 'AMR PARAMETERS'},'name' => 'commandgroup'},{'type' => 'e','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  SCHEME PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'type' => 'e','content' => [{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'default' => 'T','name' => '1'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => '2'},'name' => 'option'}],'attrib' => {'type' => 'integer','input' => 'select','name' => 'nOrder'},'name' => 'parameter'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'value' => 'Sokolov/sokolov/4/AW','name' => 'Sokolov'},'name' => 'option'}],'attrib' => {'type' => 'string','input' => 'select','name' => 'TypeFlux'},'name' => 'parameter'},{'type' => 'e','content' => [{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'default' => 'T','name' => 'minmod'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'beta'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'mc'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'LSG'},'name' => 'option'}],'attrib' => {'type' => 'string','input' => 'select','name' => 'TypeLimiter'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'if' => '$TypeLimiter eq \'beta\'','type' => 'real','default' => '1.2','max' => '2','min' => '1','name' => 'LimiterBeta'},'name' => 'parameter'}],'attrib' => {'expr' => '$nOrder == 2'},'name' => 'if'},{'type' => 't','content' => '
#SCHEME
2			nOrder (1 or 2)
Rusanov			TypeFlux
minmod			TypeLimiter ! Only for nOrder=2
1.2			LimiterBeta ! Only for LimiterType=\'beta\'

! Default values are shown above.
!
! Possible values for TypeFlux:
! \'Sokolov\'     - Sokolov\'s Local Artificial Wind flux
!
! Possible values for TypeLimiter:
! \'minmod\'	- minmod limiter is the most robust 1D limiter
! \'mc\'		- Monotonized Central limiter is sharper but less robust
! \'LSG\'		- Least Squares Gradient: robust but expensive multiD limiter 
! \'beta\'        - Beta limiter
!
! Possible values for LimiterBeta are between 1.0 and 2.0 : 
!  LimiterBeta = 1.0 is the same as the minmod limiter
!  LimiterBeta = 2.0 is the same as the superbee limiter
!  LimiterBeta = 1.2 is the recommended value
'}],'attrib' => {'name' => 'SCHEME'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'T','name' => 'UseNonConservative'},'name' => 'parameter'},{'type' => 't','content' => '
#NONCONSERVATIVE
T		UseNonConservative

! For Earth the default is using non-conservative equations 
! (close to the body).
'}],'attrib' => {'name' => 'NONCONSERVATIVE'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'integer','default' => '1','max' => '3','min' => '0','name' => 'nConservCrit'},'name' => 'parameter'},{'type' => 'e','content' => [{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'default' => 'T','value' => 'r/R/radius/Radius','name' => 'radius'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'p/P','name' => 'p'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'gradp/GradP','name' => 'grad P'},'name' => 'option'}],'attrib' => {'type' => 'string','input' => 'select','name' => 'TypeConservCrit_I'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'if' => '$TypeConservCrit_I =~ /^r|radius$/i','type' => 'real','default' => '2*$rBody','min' => '$rBody','name' => 'rConserv'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'if' => '$TypeConservCrit_I =~ /^p$/i','type' => 'real','default' => '0.05','min' => '0','name' => 'pCoeffConserv'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'if' => '$TypeConservCrit_I =~ /gradp/i','type' => 'real','default' => '0.1','min' => '0','name' => 'GradPCoeffConserv'},'name' => 'parameter'}],'attrib' => {'to' => '$nConservCrit','from' => '1'},'name' => 'for'},{'type' => 't','content' => '
#CONSERVATIVECRITERIA
3		nConservCrit
r		TypeConservCrit_I(1)
6.		rConserv             ! read if TypeConservCrit_I is \'r\'
p		TypeConservCrit_I(2)
0.05		pCoeffConserv	     ! read if TypeConservCrit_I is \'p\'
GradP		TypeConservCrit_I(3)
0.1		GradPCoeffConserv    ! read if TypeConservCrit_I is \'GradP\'

! Select the parts of the grid where the conservative vs. non-conservative
! schemes are applied. The number of criteria is arbitrary, although 
! there is no point applying the same criterion more than once.
! If no criteria is used, the whole domain will use conservative or
! non-conservative equations depending on UseNonConservative set in
! command #NONCONSERVATIVE.
!
! The physics based conservative criteria (\'p\' and \'GradP\')
! select cells which use the non-conservative scheme if ALL of them are true:
!
! \'p\'      - the pressure is smaller than fraction pCoeffConserv of the energy
! \'GradP\'  - the relative gradient of pressure is less than GradPCoeffConserv
!
! The geometry based criteria are applied after the physics based criteria 
! (if any) and they select the non-conservative scheme if ANY of them is true:
!
! \'r\'      - radial distance of the cell is less than rConserv
!
! Default values are nConservCrit = 1 with TypeConservCrit_I(1)=\'r\'
! and rConserv=2*rBody, where rBody has a problem dependent default.
'}],'attrib' => {'name' => 'CONSERVATIVECRITERIA'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'T','name' => 'UseUpdateCheck'},'name' => 'parameter'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '40','max' => '100','min' => '0','name' => 'rhoMin'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '400','min' => '100','name' => 'rhoMax'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '40','max' => '100','min' => '0','name' => 'pMin'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '400','min' => '100','name' => 'pMax'},'name' => 'parameter'}],'attrib' => {'expr' => '$UseUpdateCheck'},'name' => 'if'},{'type' => 't','content' => '
#UPDATECHECK
T			UseUpdateCheck
40.			rhoMin[%]
400.			rhoMax[%]
40.			pMin[%]
400.			pMax[%]

! Default values are shown.  This will adjust the timestep so that
! density and pressure cannot change by more than the given percentages
! in a single timestep.
'}],'attrib' => {'name' => 'UPDATECHECK'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'default' => 'T','name' => '1'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => '2'},'name' => 'option'}],'attrib' => {'type' => 'integer','input' => 'select','name' => 'nOrderProlong'},'name' => 'parameter'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'default' => 'T','value' => 'lr','name' => 'left-right'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'central'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'minmod'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'lr2','name' => 'left-right extrapolate'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'central2','name' => 'central    extrapolate'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => 'minmod2','name' => 'minmod     extrapolate'},'name' => 'option'}],'attrib' => {'if' => '$nOrderProlong==2','type' => 'string','input' => 'select','name' => 'TypeProlong'},'name' => 'parameter'},{'type' => 't','content' => '
#PROLONGATION
2			nOrderProlong (1 or 2 for ghost cells)
lr			TypeProlong  ! Only for nOrderProlong=2

! Default is prolong_order=1. 
! Possible values for prolong_type:
! 1. in message_pass_dir (used if limiter_type is not \'LSG\')
! \'lr\'		- interpolate only with left and right slopes 
! \'central\'	- interpolate only with central difference slope
! \'minmod\' 	- interpolate only with minmod limited slope
! \'lr2\'		- like \'lr\' but extrapolate when necessary
! \'central2\'	- like \'central\' but extrapolate when necessary
! \'minmod2\'	- like \'minmod\' but extrapolate when necessary
! \'lr3\'		- only experimental
!
! 2. in messagepass_all (used if limiter_type is \'LSG\')
! \'lr\',\'lr2\'		- left and right slopes (all interpolation)
! \'central\',\'central2\'	- central differences (all interpolation)
! \'minmod\',\'minmod2\'	- to be implemented
'}],'attrib' => {'name' => 'PROLONGATION'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'default' => 'T','name' => 'm_p_cell FACES ONLY','value' => 'allopt'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'm_p_cell','value' => 'all'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'm_p_dir FACES ONLY','value' => 'opt'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'm_p_dir group by directions','value' => 'dir'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'm_p_dir group by faces     ','value' => 'face'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'name' => 'm_p_dir group by kind and face','value' => 'min'},'name' => 'option'}],'attrib' => {'type' => 'string','input' => 'select','name' => 'TypeMessagePass'},'name' => 'parameter'},{'type' => 't','content' => '
#MESSAGEPASS
allopt			TypeMessagePass

! Default value is shown above.
! Possible values for optimize_message_pass
!
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
!
'}],'attrib' => {'alias' => 'OPTIMIZE','name' => 'MESSAGEPASS'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'T','name' => 'UseDivbSource'},'name' => 'parameter'},{'type' => 'e','content' => [{'type' => 't','content' => '
	! At least one of the options should be true.
	'}],'attrib' => {'expr' => '$UseDivbSource or $UseDivbDiffusion or $UseProjection or $UseConstrainB'},'name' => 'rule'},{'type' => 't','content' => '
#DIVB
T			UseDivbSource

! Default values are shown above.
! At least one of the options should be true.
'}],'attrib' => {'name' => 'DIVB'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'T','name' => 'UseB0Source'},'name' => 'parameter'},{'type' => 't','content' => '
#DIVBSOURCE
T			UseB0Source

! Add extra source terms related to the non-zero divergence and curl of B0.
! Default is true.
'}],'attrib' => {'name' => 'DIVBSOURCE'},'name' => 'command'}],'attrib' => {'name' => 'SCHEME PARAMETERS'},'name' => 'commandgroup'},{'type' => 'e','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  PHYSICS PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '1.6666666667','min' => '1','name' => 'Gamma'},'name' => 'parameter'},{'type' => 't','content' => '
#GAMMA
1.6666666667		g

! Above value is the default.
'}],'attrib' => {'if' => '$_IsFirstSession','name' => 'GAMMA'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '1','min' => '0','name' => 'RhoLeft'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '0','name' => 'UnLeft'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '0','name' => 'Ut1Left'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '0','name' => 'Ut2Left'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '0.75','name' => 'BnLeft'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '1','name' => 'Bt1Left'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '0','name' => 'Bt2Left'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '1','min' => '0','name' => 'pRight'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '0.125','min' => '0','name' => 'RhoRight'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '0','name' => 'UnRight'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '0','name' => 'Ut1Right'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '0','name' => 'Ut2Right'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '0.75','name' => 'BnRight'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '-1','name' => 'Bt1Right'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '0','name' => 'Bt2Right'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '0.1','min' => '0','name' => 'pRight'},'name' => 'parameter'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'default' => 'T','name' => 'no rotation','value' => '0'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => '0.25'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => '0.3333333333333'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => '0.5'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => '1'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => '2'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => '3'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => '4'},'name' => 'option'}],'attrib' => {'type' => 'real','input' => 'select','name' => 'ShockSlope'},'name' => 'parameter'},{'type' => 't','content' => '
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
!
! ShockSlope = 1., 2., 3., 4., 5.      .....
! ShockSlope = 0.5, 0.33333333, 0.25, 0.2, .....
!
! can be used, because these angles can be accurately represented
! on the grid.
'}],'attrib' => {'name' => 'SHOCKTUBE'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '5','min' => '0','name' => 'SwRhoDim'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '181712.175','min' => '0','name' => 'SwTDim'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '-400','max' => '0','name' => 'SwUxDim'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '0','name' => 'SwUyDim'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '0','name' => 'SwUzDim'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '0','name' => 'SwBxDim'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '0','name' => 'SwByDim'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '5','name' => 'SwBzDim'},'name' => 'parameter'},{'type' => 't','content' => '
#SOLARWIND
5.0			SwRhoDim [n/cc]
181712.175		SwTDim [K]
-400.0			SwUxDim [km/s]
0.0			SwUyDim [km/s]
0.0			SwUzDim [km/s]
0.0			SwBxDim [nT]
0.0			SwByDim [nT]
5.0			SwBzDim [nT]

! No default values!
'}],'attrib' => {'if' => '$_IsFirstSession','name' => 'SOLARWIND'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'F','name' => 'UseUpstreamInputFile'},'name' => 'parameter'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'string','length' => '100','name' => 'NameUpstreamFile'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '0','name' => 'SatelliteYPos'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '0','name' => 'SatelliteZPos'},'name' => 'parameter'}],'attrib' => {'expr' => '$UseUpstreamInputFile'},'name' => 'if'},{'type' => 't','content' => '
#UPSTREAM_INPUT_FILE
T			UseUpstreamInputFile (rest of parameters read if true)
IMF.dat                 NameUpstreamFile
0.0                     SatelliteYPos
0.0                     SatelliteZPos

! UseUpstreamInputFile - default is false
! UpstreamFileName     - user specified input file
! Satellite_Y_Pos      - not yet used
! Satellite_Z_Pos      - not yet used
'}],'attrib' => {'name' => 'UPSTREAM_INPUT_FILE'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'F','name' => 'UseBody'},'name' => 'parameter'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '3','min' => '0','name' => 'rBody'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '4','min' => '-1','name' => 'rCurrents'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '1','min' => '0','name' => 'BodyRhoDim'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '10000','min' => '0','name' => 'BodyTDim'},'name' => 'parameter'}],'attrib' => {'expr' => '$UseBody'},'name' => 'if'},{'type' => 't','content' => '
#BODY
T			UseBody (rest of parameters read if true)
3.0			rBody
4.0			rCurrents
1.0			BodyRhoDim (/ccm) density for fixed BC for rho_BLK
10000.0			BodyTDim (K) temperature for fixed BC for P_BLK

! Default values depend on problem_type.
'}],'attrib' => {'alias' => 'MAGNETOSPHERE','if' => '$_IsFirstSession','name' => 'BODY'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'F','name' => 'UseGravity'},'name' => 'parameter'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'default' => 'T','value' => '0','name' => 'central mass'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => '1','name' => 'X direction'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => '2','name' => 'Y direction'},'name' => 'option'},{'type' => 'e','content' => [],'attrib' => {'value' => '3','name' => 'Z direction'},'name' => 'option'}],'attrib' => {'if' => '$UseGravity','type' => 'integer','input' => 'select','name' => 'iDirGravity'},'name' => 'parameter'},{'type' => 't','content' => '
#GRAVITY
T			UseGravity (rest of parameters read if true)
0			GravityDir (0 - central, 1 - X, 2 - Y, 3 - Z direction)

! Default values depend on problem_type.  
'}],'attrib' => {'if' => '$_IsFirstSession','name' => 'GRAVITY'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'F','name' => 'UsePhysicalFAConductance'},'name' => 'parameter'},{'type' => 't','content' => '
#FACONDUCTIVITYMODEL
F			UsePhysicalFAConductance

Default value is shown.
'}],'attrib' => {'name' => 'FACONDUCTIVITYMODEL'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'F','name' => 'UseMassLoading'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'F','name' => 'DoAccelerateMassLoading'},'name' => 'parameter'},{'type' => 't','content' => '
#MASSLOADING
F			UseMassLoading
F			AccelerateMassLoading
'}],'attrib' => {'name' => 'MASSLOADING'},'name' => 'command'}],'attrib' => {'name' => 'PHYSICS PARAMETERS'},'name' => 'commandgroup'},{'type' => 'e','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!! SOLAR PROBLEM TYPES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '2.85E06','min' => '0','name' => 'BodyTDim'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '1.50E8','min' => '0','name' => 'BodyRhoDim'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '25.0','min' => '0','name' => 'qSun'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '1.75','min' => '0','name' => 'tHeat'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '1.0','min' => '0','name' => 'rHeat'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '4.5','min' => '0','name' => 'SigmaHeat'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'F','name' => 'DoInitRope'},'name' => 'parameter'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '0.7','min' => '0','name' => 'CmeA'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '1.2','min' => '0','name' => 'CmeR1'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '1.0','min' => '0','name' => 'CmeR0'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '0.23','min' => '0','name' => 'CmeA1'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '0.0','name' => 'CmeAlpha'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '2.5E-12','min' => '0','name' => 'CmeRho1'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '2.0E-13','min' => '0','name' => 'CmeRho2'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '0.0','max' => '10','min' => '0','name' => 'ModulationRho'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '0.0','max' => '10','min' => '0','name' => 'ModulationP'},'name' => 'parameter'}],'attrib' => {'expr' => '$DoInitRope'},'name' => 'if'},{'type' => 't','content' => '
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

! Default values are shown. Parameters for problem_heliosphere
'}],'attrib' => {'if' => '$_IsFirstSession and $_NameComp ne \'GM\'','name' => 'HELIOSPHERE'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'real','name' => 'HelioDipoleStrength'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '0','max' => '90','min' => '-90','name' => 'HelioDipoleTilt'},'name' => 'parameter'},{'type' => 't','content' => '

#HELIODIPOLE
-3.0                    HelioDipoleStrength [G]
 0.0                    HelioDipoleTilt     [deg]

! Variable HelioDipoleStrength defines the equatorial field strength in Gauss,
! while HelioDipoleTilt is the tilt relative to the ecliptic North 
! (negative sign means towards the planet) in degrees.
!
! Default values are ???
'}],'attrib' => {'if' => '$_NameComp ne \'GM\'','name' => 'HELIODIPOLE'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'T','name' => 'UseInertialFrame'},'name' => 'parameter'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'T','name' => 'UseRotatingBC'},'name' => 'parameter'}],'attrib' => {'expr' => '$UseInertialFrame'},'name' => 'if'},{'type' => 't','content' => '

#HELIOROTATION
T			UseInertialFrame
F			UseRotatingBC (read only if UseInertialFrame is true)

! If UseInertialFrame is false, the heliosphere is modeled in a corotating
! frame. In this frame the inner boundary (the solar surface) is not rotating
! (for now differential rotation is ignored). If UseInertialFrame is true,
! the heliosphere is modeled in an inertial coordinate system.
! In that case UseRotatingBC determines if the inner boundary is rotating
! or the rotation is neglected.
!
! Default values are shown. The #INERTIAL command name is obsolete.
'}],'attrib' => {'alias' => 'INERTIAL','if' => '$_IsFirstSession and $_NameComp ne \'GM\'','name' => 'HELIOROTATION'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'logical','default' => 'F','name' => 'DoSendMHD'},'name' => 'parameter'},{'type' => 't','content' => '

#HELIOTEST
F			DoSendMHD

! If DoSendMHD is true, IH sends the real MHD solution to GM in the coupling.
! If DoSendMHD is false then the values read from the IMF file are sent,
! so there is no real coupling. Mostly used for testing the framework.
!
! Default value is true, ie. real coupling.
'}],'attrib' => {'if' => '$_NameComp ne \'GM\'','name' => 'HELIOTEST'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'default' => 'T','name' => 'Low'},'name' => 'option'}],'attrib' => {'type' => 'string','input' => 'select','name' => 'TypeCme'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '0.7','min' => '0','name' => 'CmeA'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '1.2','min' => '0','name' => 'CmeR1'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '1.0','min' => '0','name' => 'CmeR0'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '0.23','min' => '0','name' => 'CmeA1'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '0.0','name' => 'CmeAlpha'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '2.5E-12','min' => '0','name' => 'CmeRho1'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '2.0E-13','min' => '0','name' => 'CmeRho2'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '1.0','name' => 'CmeB1Dim'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '4.0E5','min' => '0','name' => 'CmeUErupt'},'name' => 'parameter'},{'type' => 't','content' => '
#CME
Low		TypeCme   model type (\'Low\')
0.7		CmeA    [scaled] contraction distance
1.2             CmeR1   [scaled] distance of spheromac from sun center
1.0             CmeR0   [scaled] diameter of spheromac
0.23		CmeA1   [Gauss]  spheromac B field strength
0.0		Cmealpha   [scaled] cme acceleration rate
2.5E-12		CmeRho1 [kg/m^3] density of background corona before contract
2.0E-13		CmeRho2 [kg/m^3] density of background corona after contract 
1.0             CmeB1Dim [Gauss] field strength of dipole-type B field
4.0E5           CmeUErupt  [m/s] cme velocity

! Default values are shown above for B.C. Low\'s CME model
'}],'attrib' => {'if' => '$_IsFirstSession and $_NameComp ne \'GM\'','name' => 'CME'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '1.0E6','min' => '0','name' => 'tArcDim'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '1.0E-12','min' => '0','name' => 'RhoArcDim'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '0.718144','min' => '0','name' => 'bArcDim'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '1.0E6','name' => 'ByArcDim'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '5.0E3','name' => 'UzArcDim'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '0.5','name' => 'Phi0Arc'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '1.3','name' => 'MuArc'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '3','min' => '0','name' => 'ExpArc'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','default' => '0.5','min' => '0','name' => 'WidthArc'},'name' => 'parameter'},{'type' => 't','content' => '
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
'}],'attrib' => {'if' => '$_IsFirstSession and $_NameComp ne \'GM\'','name' => 'ARCADE'},'name' => 'command'}],'attrib' => {'name' => 'SOLAR PROBLEM TYPES'},'name' => 'commandgroup'},{'type' => 'e','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!! COMET PROBLEM TYPE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'real','min' => '0','name' => 'ProdRate'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','min' => '0','name' => 'UrNeutral'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','min' => '0','name' => 'AverageMass'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','min' => '0','name' => 'IonizationRate'},'name' => 'parameter'},{'type' => 'e','content' => [],'attrib' => {'type' => 'real','min' => '0','name' => 'kFriction'},'name' => 'parameter'},{'type' => 't','content' => '
#COMET
1.0E28		ProdRate    - Production rate (#/s)
1.0		UrNeutral   - neutral radial outflow velocity (km/s)
17.0		AverageMass - average particle mass (amu)
1.0E-6		IonizationRate (1/s)
1.7E-9		kFriction - ion-neutral friction rate coefficient (cm^3/s)

! Only used by problem_comet.  Defaults are as shown.
'}],'attrib' => {'if' => '$_IsFirstSession','name' => 'COMET'},'name' => 'command'}],'attrib' => {'name' => 'COMET PROBLEM TYPE'},'name' => 'commandgroup'},{'type' => 'e','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! SCRIPT COMMANDS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'type' => 'e','content' => [{'type' => 'e','content' => [],'attrib' => {'type' => 'string','default' => 'Param/','length' => '100','name' => 'NameIncludeFile'},'name' => 'parameter'},{'type' => 't','content' => '

#INCLUDE
Param/SSS_3000		NameIncludeFile

! Include a library file from Param/ or any file from anywhere else.
'}],'attrib' => {'name' => 'INCLUDE'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 't','content' => '

#RUN

This command is only used in stand alone mode.

Run BATSRUS with the parameters above and then return for the next session
'}],'attrib' => {'if' => '$_IsStandAlone','name' => 'RUN'},'name' => 'command'},{'type' => 'e','content' => [{'type' => 't','content' => '

#END

This command is only used in stand alone mode.

Run the executable with the parameters above and then stop.
In included files #END simply means the end of the included lines.
'}],'attrib' => {'if' => '$_IsStandAlone','name' => 'END'},'name' => 'command'}],'attrib' => {'name' => 'SCRIPT COMMANDS'},'name' => 'commandgroup'},{'type' => 'e','content' => [{'type' => 't','content' => '
	Either command #SOLARWIND or #UPSTREAM_INPUT_FILE must be used!
'}],'attrib' => {'expr' => '($SwRhoDim > 0) or $UseUpstreamInputFile or $_NameComp ne \'GM\''},'name' => 'rule'},{'type' => 'e','content' => [{'type' => 't','content' => '
	Part implicit scheme requires more than 1 implicit block!
'}],'attrib' => {'expr' => '$MaxImplBlock>1 or not $UsePartImplicit or not $MaxImplBlock'},'name' => 'rule'},{'type' => 'e','content' => [{'type' => 't','content' => '
	Full implicit scheme should be used with equal number of 
	explicit and implicit blocks!
'}],'attrib' => {'expr' => '$MaxImplBlock==$MaxBlock or not $UseFullImplicit'},'name' => 'rule'}],'attrib' => {'name' => 'Global Magnetosphere and Inner Heliosphere'},'name' => 'commandList'}];