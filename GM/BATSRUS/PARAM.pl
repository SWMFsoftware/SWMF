#^CFG FILE _FALSE_
$tree = [{'type' => 'e','attrib' => {'name' => 'Global Magnetosphere and Inner Heliosphere'},'content' => [{'type' => 't','content' => '

List of MH (GM, IH and SC) commands used in the PARAM.in file




'},{'type' => 'e','attrib' => {'type' => 'integer','value' => '$_GridSize[0]','name' => 'nI'},'content' => [],'name' => 'set'},{'type' => 'e','attrib' => {'type' => 'integer','value' => '$_GridSize[1]','name' => 'nJ'},'content' => [],'name' => 'set'},{'type' => 'e','attrib' => {'type' => 'integer','value' => '$_GridSize[2]','name' => 'nK'},'content' => [],'name' => 'set'},{'type' => 'e','attrib' => {'type' => 'integer','value' => '$_GridSize[3]','name' => 'MaxBlock'},'content' => [],'name' => 'set'},{'type' => 'e','attrib' => {'type' => 'integer','value' => '$_GridSize[4]','name' => 'MaxImplBlock'},'content' => [],'name' => 'set'},{'type' => 'e','attrib' => {'type' => 'integer','value' => '$_nProc and $MaxBlock and $_nProc*$MaxBlock','name' => 'MaxBlockALL'},'content' => [],'name' => 'set'},{'type' => 'e','attrib' => {'name' => 'STAND ALONE MODE'},'content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! STAND ALONE PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

'},{'type' => 'e','attrib' => {'if' => '$_IsStandAlone','name' => 'NEWPARAM'},'content' => [{'type' => 'e','attrib' => {'type' => 'logical','default' => 'T','name' => 'UseNewParam'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'logical','default' => 'T','name' => 'UseNewAxes'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'logical','default' => 'T','name' => 'DoTimeAccurate'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'logical','default' => 'T','name' => 'UseCorotation'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '

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
'}],'name' => 'command'},{'type' => 'e','attrib' => {'if' => '$_IsStandAlone','name' => 'DESCRIPTION'},'content' => [{'type' => 'e','attrib' => {'type' => 'string','length' => '100','name' => 'StringDescription'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '

#DESCRIPTION
This is a test run for Jupiter with no rotation.

This command is only used in the stand alone mode.

The StringDescription string can be used to describe the simulation
for which the parameter file is written. The #DESCRIPTION command and
the StringDescription string are saved into the restart file,
which helps in identifying the restart files.

The default value is "Please describe me!", which is self explanatory.
'}],'name' => 'command'},{'type' => 'e','attrib' => {'if' => '$_IsStandAlone','name' => 'ECHO'},'content' => [{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'DoEcho'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '

#ECHO
T                       DoEcho

This command is only used in the stand alone mode.

If the DoEcho variable is true, the input parameters are echoed back.
The default value for DoEcho is .false., but it is a good idea to
set it to true at the beginning of the PARAM.in file.
'}],'name' => 'command'},{'type' => 'e','attrib' => {'if' => '$_IsStandAlone','name' => 'PROGRESS'},'content' => [{'type' => 'e','attrib' => {'min' => '-1','type' => 'integer','default' => '10','name' => 'DnProgressShort'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '-1','type' => 'integer','default' => '100','name' => 'DnProgressLong'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '
#PROGRESS
10			DnProgressShort
100			DnProgressLong

The frequency of short and long progress reports for BATSRUS in
stand alone mode. These are the defaults. Set -1-s for no progress reports.
'}],'name' => 'command'},{'type' => 'e','attrib' => {'if' => '$_IsStandAlone','name' => 'TIMEACCURATE'},'content' => [{'type' => 'e','attrib' => {'type' => 'logical','default' => 'T','name' => 'DoTimeAccurate'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '

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
'}],'name' => 'command'},{'type' => 'e','attrib' => {'multiple' => 'T','if' => '$_IsStandAlone','name' => 'BEGIN_COMP'},'content' => [{'type' => 't','content' => '

This command is allowed in stand alone mode only for sake of the 
test suite, which contains these commands when the framework is tested.
'}],'name' => 'command'},{'type' => 'e','attrib' => {'multiple' => 'T','if' => '$_IsStandAlone','name' => 'END_COMP'},'content' => [{'type' => 't','content' => '

This command is allowed in stand alone mode only for sake of the 
test suite, which contains these commands when the framework is tested.
'}],'name' => 'command'},{'type' => 'e','attrib' => {'if' => '$_IsStandAlone','name' => 'RUN'},'content' => [{'type' => 't','content' => '

#RUN

This command is only used in stand alone mode.

The #RUN command does not have any parameters. It signals the end
of the current session, and makes BATSRUS execute the session with
the current set of parameters. The parameters for the next session
start after the #RUN command. For the last session there is no
need to use the #RUN command, since the #END command or simply
the end of the PARAM.in file makes BATSRUS execute the last session.
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'END'},'content' => [{'type' => 't','content' => '

#END

The #END command signals the end of the included file or the
end of the PARAM.in file. Lines following the #END command are
ignored. It is not required to use the #END command. The end
of the included file or PARAM.in file is equivalent with an 
#END command in the last line.
'}],'name' => 'command'}],'name' => 'commandgroup'},{'type' => 'e','attrib' => {'name' => 'PLANET COMMANDS'},'content' => [{'type' => 't','content' => '
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

'},{'type' => 'e','attrib' => {'if' => '$_IsFirstSession and $_IsStandAlone','name' => 'PLANET'},'content' => [{'type' => 'e','attrib' => {'type' => 'string','input' => 'select','name' => 'NamePlanet'},'content' => [{'type' => 'e','attrib' => {'value' => 'EARTH/Earth/earth','default' => 'T','name' => 'Earth'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'SATURN/Saturn/saturn','name' => 'Saturn'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'New'},'content' => [],'name' => 'option'}],'name' => 'parameter'},{'type' => 'e','attrib' => {'expr' => '$NamePlanet eq \'New\''},'content' => [{'type' => 'e','attrib' => {'min' => '0','type' => 'real','name' => 'RadiusPlanet'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '0','type' => 'real','name' => 'MassPlanet'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '0','type' => 'real','name' => 'OmegaPlanet'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '0','type' => 'real','name' => 'TiltRotation'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'string','input' => 'select','name' => 'TypeBField'},'content' => [{'type' => 'e','attrib' => {'name' => 'NONE'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'default' => 'T','name' => 'DIPOLE'},'content' => [],'name' => 'option'}],'name' => 'parameter'}],'name' => 'if'},{'type' => 'e','attrib' => {'expr' => '$TyepBField eq \'DIPOLE\''},'content' => [{'type' => 'e','attrib' => {'max' => '180','min' => '0','type' => 'real','name' => 'MagAxisThetaGeo'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'max' => '360','min' => '0','type' => 'real','name' => 'MagAxisPhiGeo'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','name' => 'DipoleStrength'},'content' => [],'name' => 'parameter'}],'name' => 'if'},{'type' => 'e','attrib' => {'expr' => 'not $PlanetCommand'},'content' => [{'type' => 't','content' => '
		PLANET should precede $PlanetCommand
	'}],'name' => 'rule'},{'type' => 't','content' => '

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
'}],'name' => 'command'},{'type' => 'e','attrib' => {'if' => '$_IsStandAlone','name' => 'ROTATIONAXIS'},'content' => [{'type' => 'e','attrib' => {'type' => 'logical','default' => 'T','name' => 'IsRotAxisPrimary'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'expr' => '$IsRotAxisPrimary'},'content' => [{'type' => 'e','attrib' => {'max' => '180','min' => '0','type' => 'real','name' => 'RotAxisTheta'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'max' => '360','min' => '0','type' => 'real','name' => 'RotAxisPhi'},'content' => [],'name' => 'parameter'}],'name' => 'if'},{'type' => 'e','attrib' => {'type' => 'string','value' => 'ROTATIONAXIS','name' => 'PlanetCommand'},'content' => [],'name' => 'set'},{'type' => 't','content' => '

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
'}],'name' => 'command'},{'type' => 'e','attrib' => {'if' => '$_IsStandAlone','name' => 'ROTATION'},'content' => [{'type' => 'e','attrib' => {'type' => 'logical','default' => 'T','name' => 'UseRotation'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'expr' => '$UseRotation'},'content' => [{'type' => 'e','attrib' => {'type' => 'real','name' => 'RotationPeriod'},'content' => [],'name' => 'parameter'}],'name' => 'if'},{'type' => 'e','attrib' => {'type' => 'string','value' => 'MAGNETICAXIS','name' => 'PlanetCommand'},'content' => [],'name' => 'set'},{'type' => 't','content' => '

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
'}],'name' => 'command'},{'type' => 'e','attrib' => {'if' => '$_IsStandAlone','name' => 'MAGNETICAXIS'},'content' => [{'type' => 'e','attrib' => {'type' => 'logical','default' => 'T','name' => 'IsMagAxisPrimary'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'expr' => '$IsMagAxisPrimary'},'content' => [{'type' => 'e','attrib' => {'max' => '180','min' => '0','type' => 'real','name' => 'MagAxisTheta'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'max' => '360','min' => '0','type' => 'real','name' => 'MagAxisPhi'},'content' => [],'name' => 'parameter'}],'name' => 'if'},{'type' => 'e','attrib' => {'type' => 'string','value' => 'MAGNETICAXIS','name' => 'PlanetCommand'},'content' => [],'name' => 'set'},{'type' => 't','content' => '

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
'}],'name' => 'command'},{'type' => 'e','attrib' => {'if' => '$_IsStandAlone','name' => 'DIPOLE'},'content' => [{'type' => 'e','attrib' => {'type' => 'real','name' => 'DipoleStrength'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '

#DIPOLE
-3.11e-4		DipoleStrength [Tesla]

The DipoleStrength variable contains the
magnetic equatorial strength of the dipole magnetic field in Tesla.

The default value is the real dipole strength for the planet.
For the Earth the default is taken to be -31100 nT.
The sign is taken to be negative so that the magnetic axis can
point northward as usual.
'}],'name' => 'command'},{'type' => 'e','attrib' => {'if' => '$_IsStandAlone','name' => 'UPDATEB0'},'content' => [{'type' => 'e','attrib' => {'min' => '-1','type' => 'real','default' => '0.0001','name' => 'DtUpdateB0'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '

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
'}],'name' => 'command'},{'type' => 'e','attrib' => {'if' => '$_IsStandAlone','name' => 'IDEALAXES'},'content' => [{'type' => 't','content' => '

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
'}],'name' => 'command'}],'name' => 'commandgroup'},{'type' => 'e','attrib' => {'name' => 'USER DEFINED INPUT'},'content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!  USER DEFINED INPUT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

'},{'type' => 'e','attrib' => {'name' => 'USER_FLAGS'},'content' => [{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'UseUserInnerBcs'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'UseUserSource'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'UseUserPerturbation'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'UseUserOuterBcs'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'UseUserICs'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'UseUserSpecifyRefinement'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'UseUserLogFiles'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'UseUserWritePlot'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'UseUserAMR'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'UseUserEchoInput'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'UseUserB0'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'UseUserSetPhysConst'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'UseUserUpdateStates'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '

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
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'USERINPUTBEGIN'},'content' => [{'type' => 't','content' => '

This command signals the beginning of the section of the file which 
is read by the subroutine user\\_read\\_inputs in the user\\_routines.f90 file.
The section ends with the #USERINPUTEND command. There is no XML based parameter
checking in the user section.
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'USERINPUTEND'},'content' => [{'type' => 't','content' => '

This command signals the end of the section of the file which 
is read by the subroutine user\\_read\\_inputs in the user\\_routines.f90 file.
The section begins with the #USERINPUTBEGIN command. There is no XML based parameter
checking in the user section.
'}],'name' => 'command'}],'name' => 'commandgroup'},{'type' => 'e','attrib' => {'name' => 'TESTING AND TIMING'},'content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  TESTING AND TIMING PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'type' => 'e','attrib' => {'name' => 'TEST'},'content' => [{'type' => 'e','attrib' => {'type' => 'string','length' => '100','name' => 'TestString'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '
#TEST
read_inputs

! A space separated list of subroutine names. Default is empty string.
!
! Examples:
!   read_inputs  - echo the input parameters following the #TEST line
!   project_B    - info on projection scheme   
!   implicit     - info on implicit scheme     
!   krylov       - info on the Krylov solver   
!   message_count- count messages
!   initial_refinement
!   ...
! Check the subroutines for call setoktest("...",oktest,oktest_me) to
! see the appropriate strings.
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'TESTIJK'},'content' => [{'type' => 'e','attrib' => {'max' => '$nI+2','min' => '-2','type' => 'integer','name' => 'iTest'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'max' => '$nJ+2','min' => '-2','type' => 'integer','name' => 'jTest'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'max' => '$nK+2','min' => '-2','type' => 'integer','name' => 'kTest'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'max' => '$MaxBlock','min' => '1','type' => 'integer','name' => 'iBlockTest'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '0','type' => 'integer','name' => 'iProcTest'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '
#TESTIJK
1                       iTest           (cell index for testing)
1                       jTest           (cell index for testing)
1                       kTest           (cell index for testing)
1                       iBlockTest      (block index for testing)
0                       iProcTest       (processor index for testing)

! The location of test info in terms of indices, block and processor number.
! Note that the user should set #TESTIJK or #TESTXYZ, not both.  If both
! are set, the final one in the session will set the test point.
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'TESTXYZ'},'content' => [{'type' => 'e','attrib' => {'max' => '$xMax','min' => '$xMin','type' => 'real','name' => 'xTest'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'max' => '$yMax','min' => '$yMin','type' => 'real','name' => 'yTest'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'max' => '$zMax','min' => '$zMin','type' => 'real','name' => 'zTest'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '
#TESTXYZ
1.5                     xTest           (X coordinate of cell for testing)
-10.5                   yTest           (Y coordinate of cell for testing)
-10.                    zTest           (Z coordinate of cell for testing)

! The location of test info in terms of coordinates.
! Note that the user should set #TESTIJK or #TESTXYZ, not both.  If both
! are set, the final one in the session will set the test point.
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'TESTTIME'},'content' => [{'type' => 'e','attrib' => {'min' => '-1','type' => 'integer','default' => '-1','name' => 'nIterTest'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '-1','type' => 'real','default' => '1e30','name' => 'TimeTest'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '

#TESTTIME
-1                      nIterTest       (iteration number to start testing)
10.5                    TimeTest        (time to start testing in seconds)

! The time step and physical time to start testing.
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'TESTVAR'},'content' => [{'type' => 'e','attrib' => {'type' => 'integer','input' => 'select','name' => 'iVarTest'},'content' => [{'type' => 'e','attrib' => {'value' => '1','default' => 'T','name' => 'Rho'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => '2','name' => 'RhoUx'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => '3','name' => 'RhoUy'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => '4','name' => 'RhoUz'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => '5','name' => 'Bx'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => '6','name' => 'By'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => '7','name' => 'Bz'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => '8','name' => 'e'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => '9','name' => 'p'},'content' => [],'name' => 'option'}],'name' => 'parameter'},{'type' => 't','content' => '
#TESTVAR
1                       iVarTest

! Index of variable to be tested. Default is rho_="1", ie. density.
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'TESTDIM'},'content' => [{'type' => 'e','attrib' => {'type' => 'integer','input' => 'select','name' => 'iVarTest'},'content' => [{'type' => 'e','attrib' => {'value' => '0','name' => 'all'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => '1','default' => 'T','name' => 'x'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => '2','name' => 'y'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => '3','name' => 'z'},'content' => [],'name' => 'option'}],'name' => 'parameter'},{'type' => 't','content' => '
#TESTDIM
1                       iDimTest

! Index of dimension/direction to be tested. Default is X dimension.
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'STRICT'},'content' => [{'type' => 'e','attrib' => {'type' => 'logical','default' => 'T','name' => 'UseStrict'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '
#STRICT
T                       UseStrict

! If true then stop when parameters are incompatible. If false, try to
! correct parameters and continue. Default is true, ie. strict mode
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'VERBOSE'},'content' => [{'type' => 'e','attrib' => {'type' => 'integer','input' => 'select','name' => 'lVerbose'},'content' => [{'type' => 'e','attrib' => {'value' => '-1','name' => 'errors and warnings only'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => '0','name' => 'start and end of sessions'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => '1','default' => 'T','name' => 'normal'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => '10','name' => 'calls on test processor'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => '100','name' => 'calls on all processors'},'content' => [],'name' => 'option'}],'name' => 'parameter'},{'type' => 't','content' => '
#VERBOSE
-1                      lVerbose

! Verbosity level controls the amount of output to STDOUT. Default level is 1.
!   lVerbose .le. -1 only warnings and error messages are shown.
!   lVerbose .ge.  0 start and end of sessions is shown.
!   lVerbose .ge.  1 a lot of extra information is given.
!   lVerbose .ge. 10 all calls of set_oktest are shown for the test processor.
!   lVerbose .ge.100 all calls of set_oktest are shown for all processors.
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'DEBUG'},'content' => [{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'DoDebug'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'DoDebugGhost'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '
#DEBUG
F                       DoDebug         (use it as if(okdebug.and.oktest)...)
F                       DoDebugGhost    (parameter for show_BLK in library.f90)

! Excessive debug output can be controlled by the global okdebug parameter
'}],'name' => 'command'},{'type' => 'e','attrib' => {'if' => '$_IsFirstSession','name' => 'CODEVERSION'},'content' => [{'type' => 'e','attrib' => {'type' => 'real','min' => '0','default' => '7.50','name' => 'CodeVersion'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '
#CODEVERSION
7.50                    CodeVersion

! Checks CodeVersion. Prints a WARNING if it differs from the CodeVersion
! defined in ModMain.f90. Used in newer restart header files. 
! Should be given in PARAM.in when reading old restart files, 
! which do not have version info in the header file.
'}],'name' => 'command'},{'type' => 'e','attrib' => {'if' => '$_IsFirstSession','name' => 'EQUATION'},'content' => [{'type' => 'e','attrib' => {'type' => 'string','default' => 'MHD','length' => '100','name' => 'NameEquation'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'integer','default' => '8','name' => 'nVar'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '
#EQUATION
MHD			NameEquation
8			nVar

! Define the equation name and the number of variables.
! If any of these do not agree with the values determined 
! by the code, BATSRUS stops with an error. Used in restart
! header files and can be given in PARAM.in as a check
! and as a description.
'}],'name' => 'command'},{'type' => 'e','attrib' => {'if' => '$_IsFirstSession','name' => 'PRECISION'},'content' => [{'type' => 'e','attrib' => {'type' => 'integer','input' => 'select','name' => 'nByteReal'},'content' => [{'type' => 'e','attrib' => {'value' => '4','default' => '$_nByteReal==4','name' => 'single precision (4)'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => '8','default' => '$_nByteReal==8','name' => 'double precision (8)'},'content' => [],'name' => 'option'}],'name' => 'parameter'},{'type' => 'e','attrib' => {'expr' => '$nByteReal==$_nByteReal'},'content' => [{'type' => 't','content' => '
		nByteReal in file must agree with _nByteReal.
	'}],'name' => 'rule'},{'type' => 't','content' => '

#PRECISION
8                       nByteReal

! Define the number of bytes in a real number. If it does not agree
! with the value determined by the code, BATSRUS stops with an error.
! This is a check, the internal value is calculated in parallel_setup.
! Used in latest restart header files to check binary compatibility.
! May be given in PARAM.in to enforce a certain precision.
'}],'name' => 'command'},{'type' => 'e','attrib' => {'if' => '$_IsFirstSession','name' => 'CHECKGRIDSIZE'},'content' => [{'type' => 'e','attrib' => {'max' => '$nI','min' => '$nI','type' => 'integer','default' => '$nI','name' => 'nI'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'max' => '$nJ','min' => '$nJ','type' => 'integer','default' => '$nJ','name' => 'nJ'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'max' => '$nK','min' => '$nK','type' => 'integer','default' => '$nK','name' => 'nK'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'max' => '$MaxBlockALL','min' => '1','type' => 'integer','name' => 'MinBlockALL'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '

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
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'BLOCKLEVELSRELOADED'},'content' => [{'type' => 't','content' => '
#BLOCKLEVELSRELOADED

This command means that the restart file contains the information about
the minimum and maximum allowed refinement levels for each block.
This command is only used in the restart header file.
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'TIMING'},'content' => [{'type' => 'e','attrib' => {'type' => 'logical','default' => 'T','name' => 'UseTiming'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'expr' => '$UseTiming'},'content' => [{'type' => 'e','attrib' => {'type' => 'integer','input' => 'select','name' => 'Frequency'},'content' => [{'type' => 'e','attrib' => {'value' => '-3','name' => 'none'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => '-2','default' => 'T','name' => 'final only'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => '-1','name' => 'end of sessions'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'min' => '1','default' => '100','name' => 'every X steps'},'content' => [],'name' => 'optioninput'}],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '-1','type' => 'integer','default' => '-1','name' => 'nDepthTiming'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'string','input' => 'select','name' => 'TypeTimingReport'},'content' => [{'type' => 'e','attrib' => {'value' => 'cumm','default' => '1','name' => 'cummulative'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'list'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'tree'},'content' => [],'name' => 'option'}],'name' => 'parameter'}],'name' => 'if'},{'type' => 't','content' => '
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
'}],'name' => 'command'}],'name' => 'commandgroup'},{'type' => 'e','attrib' => {'name' => 'INITIAL AND BOUNDARY CONDITIONS'},'content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! MAIN INITIAL AND BOUNDARY CONDITION PARAMETERS  !!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'type' => 'e','attrib' => {'required' => 'T','if' => '$_IsFirstSession','name' => 'PROBLEMTYPE'},'content' => [{'type' => 'e','attrib' => {'type' => 'integer','input' => 'select','name' => 'iProblem'},'content' => [{'type' => 'e','attrib' => {'value' => '1','name' => 'Uniform'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => '2','name' => 'Shock tube'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => '3','name' => 'Heliosphere'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => '5','name' => 'Comet'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => '6','name' => 'Rotation'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => '7','name' => 'Diffusion'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => '11','default' => 'T','name' => 'Earth'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => '12','name' => 'Saturn'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => '13','name' => 'Jupiter'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => '14','name' => 'Venus'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => '21','name' => 'Cylinder'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => '22','name' => 'Sphere'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => '25','name' => 'Arcade'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => '26','name' => 'CME'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => '30','name' => 'Dissipation'},'content' => [],'name' => 'option'}],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'string','length' => '20','if' => '$iProblem==30','name' => 'TypeDissipation'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '
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
'}],'name' => 'command'},{'type' => 'e','attrib' => {'if' => '$_IsFirstSession','name' => 'COORDSYSTEM'},'content' => [{'type' => 'e','attrib' => {'type' => 'string','input' => 'select','name' => 'TypeCoordSystem'},'content' => [{'type' => 'e','attrib' => {'value' => 'GSM','default' => 'T','if' => '$_NameComp eq \'GM\'','name' => 'GeoSolarMagnetic, GSM'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'HGI','default' => 'T','if' => '$_NameComp ne \'GM\'','name' => 'HelioGraphicInertial, HGI'},'content' => [],'name' => 'option'}],'name' => 'parameter'},{'type' => 't','content' => '

#COORDSYSTEM
GSM			TypeCoordSystem

! TypeCoordSystem defines the coordinate system for the component.
! Currently only one coordinate system is available for GM ("GSM")
! and one for IH or SC ("HGI"). In the near future "GSE" should be also
! an option for GM.
!
! Default is component dependent: "GSM" for GM and "HGI" for IH or SC.
'}],'name' => 'command'},{'type' => 'e','attrib' => {'if' => '$_IsFirstSession','name' => 'RESTARTINDIR'},'content' => [{'type' => 'e','attrib' => {'type' => 'string','default' => 'GM/restartIN','length' => '100','name' => 'NameRestartInDir'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '

#RESTARTINDIR
GM/restart_n5000	NameRestartInDir

! The NameRestartInDir variable contains the name of the directory
! where restart files are saved relative to the run directory.
! The directory should be inside the subdirectory with the name 
! of the component.
!
! Default value is "GM/restartIN".
'}],'name' => 'command'},{'type' => 'e','attrib' => {'if' => '$_IsFirstSession','name' => 'NEWRESTART'},'content' => [{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'DoRestartBFace'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => ' 

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

'}],'name' => 'command'},{'type' => 'e','attrib' => {'required' => 'T','if' => '$_IsFirstSession','name' => 'GRID'},'content' => [{'type' => 'e','attrib' => {'min' => '1','type' => 'integer','default' => '2','name' => 'nRootBlockX'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '1','type' => 'integer','default' => '1','name' => 'nRootBlockY'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '1','type' => 'integer','default' => '1','name' => 'nRootBlockZ'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '-192.0','name' => 'xMin'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '$xMin','type' => 'real','default' => '  64.0','name' => 'xMax'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => ' -64.0','name' => 'yMin'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '$yMin','type' => 'real','default' => '  64.0','name' => 'yMax'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => ' -64.0','name' => 'zMin'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '$zMin','type' => 'real','default' => '  64.0','name' => 'zMax'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '
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
! number of blocks of the base grid, ie. the roots of the octree. 
! By varying these parameters, one can setup a grid which is elongated
! in some direction. The xMin ... zMax parameters define the physical
! size of the grid.
!
! There is no default value, the grid size must always be given.
! The #GRID command should be used before the #SAVEPLOT command.
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'OUTERBOUNDARY'},'content' => [{'type' => 'e','attrib' => {'values' => 'TypeBcEast,TypeBcWest,TypeBcSouth,TypeBcNorth,TypeBcBot,TypeBcTop','name' => 'Side'},'content' => [{'type' => 'e','attrib' => {'type' => 'string','input' => 'select','name' => '$Side'},'content' => [{'type' => 'e','attrib' => {'name' => 'coupled'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'default' => '$Side ne \'TypeBcEast\'','name' => 'fixed/inflow'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'default' => '$Side eq \'TypeBcEast\'','name' => 'float/outflow'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'heliofloat'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'reflect'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'periodic'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'vary'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'shear'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'linetied'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'raeder'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'arcadetop'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'arcadebot'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'arcadebotcont'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'user'},'content' => [],'name' => 'option'}],'name' => 'parameter'}],'name' => 'foreach'},{'type' => 'e','attrib' => {'expr' => 'not($TypeBcEast eq \'periodic\' xor $TypeBcWest eq \'periodic\')'},'content' => [{'type' => 't','content' => '
	! East and west BCs must be both periodic or neither
	'}],'name' => 'rule'},{'type' => 'e','attrib' => {'expr' => 'not($TypeBcSouth eq \'periodic\' xor $TypeBcNorth eq \'periodic\')'},'content' => [{'type' => 't','content' => '
	! South and North BCs must be both periodic or neither
	'}],'name' => 'rule'},{'type' => 'e','attrib' => {'expr' => 'not($TypeBcBot eq \'periodic\' xor $TypeBcTop eq \'periodic\')'},'content' => [{'type' => 't','content' => '
	! Bottom and top BCs must be both periodic or neither
	'}],'name' => 'rule'},{'type' => 't','content' => '
#OUTERBOUNDARY
outflow                 TypeBcEast
inflow                  TypeBcWest
float                   TypeBcSouth
float                   TypeBcNorth
float                   TypeBcBottom
float                   TypeBcTop

! Default depends on problem type.
! Possible values:
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
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'INNERBOUNDARY'},'content' => [{'type' => 't','content' => '
! Inner boundary types for body 1 and body 2
	'},{'type' => 'e','attrib' => {'type' => 'string','input' => 'select','name' => 'TypeBcInner'},'content' => [{'type' => 'e','attrib' => {'name' => 'reflect'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'float'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'fixed'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'default' => 'T','name' => 'ionosphere'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'ionosphereB0/ionosphereb0','name' => 'ionosphereB0'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'ionospherefloat'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'user'},'content' => [],'name' => 'option'}],'name' => 'parameter'},{'type' => 'e','attrib' => {'expr' => '$UseBody2'},'content' => [{'type' => 'e','attrib' => {'type' => 'string','input' => 'select','name' => 'TypeBcBody2'},'content' => [{'type' => 'e','attrib' => {'default' => 'T','name' => 'reflect'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'float'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'fixed'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'ionosphere'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'ionosphereB0'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'ionospherefloat'},'content' => [],'name' => 'option'}],'name' => 'parameter'}],'name' => 'if'},{'type' => 'e','attrib' => {'expr' => 'not($TypeBcBody2 =~ /ionosphere/)'},'content' => [{'type' => 't','content' => '
! Note: for the second body COROTATION AND AN IONOSPHERIC BOUNDARY DO NOT WORK.
	'}],'name' => 'rule'},{'type' => 't','content' => '

#INNERBOUNDARY
ionosphere              TypeBcInner

ionosphere              TypeBcBody2  !read only if UseBody2=.true. 

! This command should appear after the #SECONDBODY command when using 
! two bodies. Note: for the second body COROTATION AND AN IONOSPHERIC 
! BOUNDARY DO NOT WORK.
!
! Default value for TypeBcBody2 is \'reflect\'.


! Possible values for TypeBcInner are:
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
! \'user\'        - user defined
!
! For \'ionosphere\' and \'ionospherefloat\' types and a coupled GM-IE run,
! the velocity at the inner boundary is determined by the ionosphere model.
!
! Default value for TypeBcInner is \'ionosphere\' for problem types
! Earth, Saturn, Jupiter, and rotation.
! For all other problems with an inner boundary the default is \'unknown\',
! so the inner boundary must be set.
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'EXTRABOUNDARY'},'content' => [{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'UseExtraBoundary'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'expr' => '$UseExtraBoundary'},'content' => [{'type' => 'e','attrib' => {'type' => 'string','name' => 'TypeBcExtra'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'DoFixExtraboundary'},'content' => [],'name' => 'parameter'}],'name' => 'if'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'FACEOUTERBC'},'content' => [{'type' => 'e','attrib' => {'max' => '6','min' => '0','type' => 'integer','default' => '0','name' => 'MaxBoundary'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'expr' => '$MaxBoundary >= 1'},'content' => [{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'DoFixOuterBoundary'},'content' => [],'name' => 'parameter'}],'name' => 'if'},{'type' => 't','content' => '
#FACEOUTERBC
0              MaxBoundary            
F              DoFixOuterBoundary)    !read only for MaxBoundary>=East_(=1).
! If MaxBoundary is East_(=1) or more then the outer boundaries with
! the number of boundary being between East_ and MaxBoundary
! are treated using set_BCs.f90 subroutines instead of set_outerBCs.f90 
! if DoFixOuterBoundary is .true., there is no resolution
! change along the outer boundaries with the number of
! of boundary being between East_ and MaxBoundary
'}],'name' => 'command'}],'name' => 'commandgroup'},{'type' => 'e','attrib' => {'name' => 'INITIAL TIME'},'content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! INITIAL TIME AND STEP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'type' => 'e','attrib' => {'alias' => 'SETREALTIME','if' => '$_IsFirstSession','name' => 'STARTTIME'},'content' => [{'type' => 'e','attrib' => {'type' => 'integer','default' => '2000','name' => 'iYear'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'max' => '12','min' => '1','type' => 'integer','default' => '3','name' => 'iMonth'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'max' => '31','min' => '1','type' => 'integer','default' => '21','name' => 'iDay'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'max' => '23','min' => '0','type' => 'integer','default' => '0','name' => 'iHour'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'max' => '59','min' => '0','type' => 'integer','default' => '0','name' => 'iMinute'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'max' => '59','min' => '0','type' => 'integer','default' => '0','name' => 'iSecond'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '
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
'}],'name' => 'command'},{'type' => 'e','attrib' => {'if' => '$_IsFirstSession','name' => 'TIMESIMULATION'},'content' => [{'type' => 'e','attrib' => {'min' => '0','type' => 'real','default' => '0.0','name' => 'tSimulation'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '

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
'}],'name' => 'command'},{'type' => 'e','attrib' => {'if' => '$_IsFirstSession','name' => 'NSTEP'},'content' => [{'type' => 'e','attrib' => {'min' => '0','type' => 'integer','default' => '0','name' => 'nStep'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '

#NSTEP
100			nStep

! Set nStep for the component. Typically used in the restart.H header file.
! Generally it is not inserted in a PARAM.in file by the user.
!
! The default is nStep=0 as the starting time step with no restart.
'}],'name' => 'command'},{'type' => 'e','attrib' => {'if' => '$_IsFirstSession','name' => 'NPREVIOUS'},'content' => [{'type' => 'e','attrib' => {'min' => '-1','type' => 'integer','default' => '-1','name' => 'nPrevious'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '

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
'}],'name' => 'command'}],'name' => 'commandgroup'},{'type' => 'e','attrib' => {'name' => 'TIME INTEGRATION'},'content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  TIME INTEGRATION PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'type' => 'e','attrib' => {'name' => 'TIMESTEPPING'},'content' => [{'type' => 'e','attrib' => {'type' => 'integer','input' => 'select','name' => 'nStage'},'content' => [{'type' => 'e','attrib' => {'value' => '1','default' => 'T'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => '2'},'content' => [],'name' => 'option'}],'name' => 'parameter'},{'type' => 'e','attrib' => {'max' => '1','min' => '0','type' => 'real','default' => '0.8','name' => 'CflExpl'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '

#TIMESTEPPING
2                       nStage
0.80                    CflExpl

! Parameters for explicit time integration.
! Default is 1 stage and CflExpl=0.8
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'FIXEDTIMESTEP'},'content' => [{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'UseDtFixed'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '0','type' => 'real','default' => '1.0','if' => '$UseDtFixed','name' => 'DtFixedDim'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '
#FIXEDTIMESTEP
T                       UseDtFixed
10.                     DtFixedDim [sec] (read if UseDtFixed is true)

! Default is UseDtFixed=.false. Effective only if DoTimeAccurate is true.
! If UseDtFixed is true, the time step is fixed to DtFixedDim.
!
! This is useful for debugging explicit schemes.

! The real application is, however, for implicit and partially
! implicit/local schemes.

'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'PARTLOCAL'},'content' => [{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'UsePartLocal'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '
#PARTLOCAL
T               UsePartLocal

! Default is UsePartLocal=.false. If UsePartLocal is true and the
! run is time accurate, then the blocks selected as "implicit"
! by the criteria defined in #STEPPINGCRITERIA are not used to
! calculate the time step, and all cells are advanced with the
! smaller of the stable and the global time steps.
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'IMPLICIT'},'content' => [{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'UsePointImplicit'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'UsePartImplicit'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'UseFullImplicit'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '0','type' => 'real','default' => '100','if' => '$UsePartImplicit or $UseFullImplicit','name' => 'CflImpl'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'expr' => '$UsePointImplicit + $UsePartImplicit + $UseFullImplicit <= 1'},'content' => [{'type' => 't','content' => '
	At most one of these logicals can be true!
	'}],'name' => 'rule'},{'type' => 't','content' => '

#IMPLICIT
F               UsePointImplicit   
F               UsePartImplicit
F               UseFullImplicit
100.0           CflImpl (read if UsePartImplicit or UseFullImplicit is true)

! Default is false for all logicals. Only one of them can be set to true!
! The CFL number is used in the implicit blocks of the fully or partially
! implicit schemes. Ignored if UseDtFixed is true.
'}],'name' => 'command'}],'name' => 'commandgroup'},{'type' => 'e','attrib' => {'name' => 'IMPLICIT PARAMETERS'},'content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! PARAMETERS FOR FULL AND PART IMPLICIT TIME INTEGRATION !!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'type' => 'e','attrib' => {'alias' => 'STEPPINGCRITERIA','name' => 'IMPLICITCRITERIA'},'content' => [{'type' => 'e','attrib' => {'type' => 'string','input' => 'select','name' => 'TypeImplCrit'},'content' => [{'type' => 'e','attrib' => {'value' => 'dt','default' => 'T','name' => 'Time step'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'r/R','name' => 'Radial distance'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'test','name' => 'Test block'},'content' => [],'name' => 'option'}],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '0','type' => 'real','if' => '$TypeImplCrit eq \'R\'','name' => 'rImplicit'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '

#IMPLICITCRITERIA
r		TypeImplCrit (dt or r or test)
10.0		rImplicit    (only read for TypeImplCrit = r)

! Both #IMPLICITCRITERIA and #STEPPINGCRITERIA are acceptable.
! Only effective if PartImplicit or PartLocal is true in a time accurate run.
! Default value is ImplCritType=\'dt\'.
!
! The options are
!
! if     (TypeImplCrit ==\'dt\'  ) then blocks with DtBLK .gt. DtFixed
! elseif (TypeImplCrit ==\'r\'   ) then blocks with rMinBLK .lt. rImplicit
! elseif (TypeImplCrit ==\'test\') then block iBlockTest on processor iProcTest
!
! are handled with local/implicit scheme. Here DtBlock is the time step
! allowed by the CFL condition for a given block, while rMinBLK is the
! smallest radial distance for all the cells in the block.
! The iBlockTest and iProcTest can be defined in the #TESTIJK command.
! DtFixed must be defined in the #FIXEDTIMESTEP command.
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'IMPLSTEP'},'content' => [{'type' => 'e','attrib' => {'max' => '1','min' => '0','type' => 'real','default' => '1','name' => 'ImplCoeff'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'logical','default' => 'T','name' => 'UseBdf2'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'UseSourceImpl'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '
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
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'IMPLSCHEME'},'content' => [{'type' => 'e','attrib' => {'type' => 'integer','input' => 'select','name' => 'nOrderImpl'},'content' => [{'type' => 'e','attrib' => {'default' => 'T','name' => '1'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => '2'},'content' => [],'name' => 'option'}],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'string','input' => 'select','name' => 'TypeFluxImpl'},'content' => [{'type' => 'e','attrib' => {'value' => 'Roe/roe/1','name' => 'Roe'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'Rusanov/rusanov/2/TVDLF','default' => 'T','name' => 'Rusanov'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'Linde/linde/3/HLLEL','name' => 'Linde'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'Sokolov/sokolov/4/AW','name' => 'Sokolov'},'content' => [],'name' => 'option'}],'name' => 'parameter'},{'type' => 't','content' => '
#IMPLSCHEME
1               nOrderImpl
Rusanov         TypeFluxImpl

! Default values are shown, ie. first order Rusanov scheme.
! This defines the scheme used in the implicit part.
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'NEWTON'},'content' => [{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'UseConservativeImplicit'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'UseNewton'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'expr' => '$UseNewton'},'content' => [{'type' => 'e','attrib' => {'type' => 'logical','default' => 'T','name' => 'UseNewMatrix'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '1','type' => 'integer','default' => '10','name' => 'MaxIterNewton'},'content' => [],'name' => 'parameter'}],'name' => 'if'},{'type' => 't','content' => '
#NEWTON
F		UseConservativeImplicit
T               UseNewton
F               UseNewMatrix  (only read if UseNewton is true)
10              MaxIterNewton (only read if UseNewton is true)

! Default is UseConservativeImplicit=F and UseNewton=F, ie. 
! no conservative fix is used and only one "Newton" iteration is done.
! UseNewMatrix decides whether the Jacobian should be recalculated
! for every Newton iteration. MaxIterNewton is the maximum number
! of Newton iterations before giving up.
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'JACOBIAN'},'content' => [{'type' => 'e','attrib' => {'type' => 'string','input' => 'select','name' => 'TypeJacobian'},'content' => [{'type' => 'e','attrib' => {'value' => 'prec','default' => 'T','name' => 'Preconditioned'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'free','name' => 'No preconditioning'},'content' => [],'name' => 'option'}],'name' => 'parameter'},{'type' => 'e','attrib' => {'max' => '1.e-5','min' => '0','type' => 'real','default' => '$doublePrecision ? 1.e-12 : 1.e-6','name' => 'JacobianEps'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '
#JACOBIAN
prec            TypeJacobian (prec, free)
1.E-12          JacobianEps

! The Jacobian matrix is always calculated with a matrix free approach,
! however it can be preconditioned  (\'prec\'), or not (\'free\')
! Default value is TypeJacobian=\'prec\'.
! JacobianEps contains the machine round off error for numerical derivatives.
! The default value is 1.E-12 for 8 byte reals and 1.E-6 for 4 byte reals.
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'PRECONDITIONER'},'content' => [{'type' => 'e','attrib' => {'type' => 'string','input' => 'select','name' => 'TypePrecondSide'},'content' => [{'type' => 'e','attrib' => {'name' => 'left'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'default' => 'T','name' => 'symmetric'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'right'},'content' => [],'name' => 'option'}],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'string','input' => 'select','name' => 'TypePrecond'},'content' => [{'type' => 'e','attrib' => {'default' => 'T','name' => 'MBILU'},'content' => [],'name' => 'option'}],'name' => 'parameter'},{'type' => 'e','attrib' => {'max' => '1','min' => '0','type' => 'real','default' => '0.5','name' => 'GustafssonPar'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '
#PRECONDITIONER
symmetric       TypePrecondSide (left, symmetric, right)
MBILU           TypePrecond (MBILU)
0.5             GustafssonPar (0. no modification, 1. full modification)

! Default parameters are shown. Right preconditioning does not affect
! the normalization of the residual. The Gustafsson parameter determines
! how much the MBILU preconditioner is modified. The default 0.5 value
! means a relaxed modification.
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'KRYLOV'},'content' => [{'type' => 'e','attrib' => {'type' => 'string','input' => 'select','name' => 'TypeKrylov'},'content' => [{'type' => 'e','attrib' => {'default' => 'T','name' => 'gmres'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'bicgstab'},'content' => [],'name' => 'option'}],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'string','input' => 'select','name' => 'TypeInitKrylov'},'content' => [{'type' => 'e','attrib' => {'value' => 'nul','default' => 'T','name' => '0'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'old','name' => 'previous'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'explicit'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'explicit','name' => 'scaled explicit'},'content' => [],'name' => 'option'}],'name' => 'parameter'},{'type' => 'e','attrib' => {'max' => '0.1','min' => '0','type' => 'real','default' => '0.001','name' => 'ErrorMaxKrylov'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '1','type' => 'integer','default' => '100','name' => 'MaxMatvecKrylov'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '
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
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'KRYLOVSIZE'},'content' => [{'type' => 'e','attrib' => {'min' => '1','type' => 'integer','default' => 'MaxMatvecKrylov','name' => 'nKrylovVector'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '
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
'}],'name' => 'command'}],'name' => 'commandgroup'},{'type' => 'e','attrib' => {'name' => 'STOPPING CRITERIA'},'content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! STOPPING CRITERIA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

The commands in this group only work in stand alone mode.

'},{'type' => 'e','attrib' => {'required' => '$_IsStandAlone','if' => '$_IsStandAlone','name' => 'STOP'},'content' => [{'type' => 'e','attrib' => {'min' => '-1','type' => 'integer','default' => '-1','name' => 'MaxIteration'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '-1','type' => 'real','default' => '-1','name' => 'tSimulationMax'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '

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
'}],'name' => 'command'},{'type' => 'e','attrib' => {'if' => '$_IsStandAlone','name' => 'CHECKSTOPFILE'},'content' => [{'type' => 'e','attrib' => {'type' => 'logical','default' => 'T','name' => 'DoCheckStopFile'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '

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
'}],'name' => 'command'},{'type' => 'e','attrib' => {'if' => '$_IsStandAlone','name' => 'CPUTIMEMAX'},'content' => [{'type' => 'e','attrib' => {'min' => '-1','type' => 'real','default' => '-1','name' => 'CpuTimeMax'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '

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
'}],'name' => 'command'}],'name' => 'commandgroup'},{'type' => 'e','attrib' => {'name' => 'OUTPUT PARAMETERS'},'content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  OUTPUT PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'type' => 'e','attrib' => {'name' => 'RESTARTOUTDIR'},'content' => [{'type' => 'e','attrib' => {'type' => 'string','default' => 'GM/restartOUT','length' => '100','name' => 'NameRestartOutDir'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '

#RESTARTOUTDIR
GM/restart_n5000	NameRestartOutDir

! The NameRestartOutDir variable contains the name of the directory
! where restart files are saved relative to the run directory.
! The directory should be inside the subdirectory with the name 
! of the component.
!
! Default value is "GM/restartOUT".
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'SAVERESTART'},'content' => [{'type' => 'e','attrib' => {'type' => 'logical','default' => 'T','name' => 'DoSaveRestart'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'expr' => '$DoSaveRestart'},'content' => [{'type' => 'e','attrib' => {'min' => '-1','type' => 'integer','default' => '-1','name' => 'DnSaveRestart'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '-1','type' => 'real','default' => '-1','name' => 'DtSaveRestart'},'content' => [],'name' => 'parameter'}],'name' => 'if'},{'type' => 't','content' => '
#SAVERESTART
T			DoSaveRestart Rest of parameters read if true
100			DnSaveRestart
-1.			DtSaveRestart [seconds]

! Default is DoSaveRestart=.true. with DnSaveRestart=-1 and 
! DtSaveRestart=-1. This results in the restart file being 
! saved only at the end.  A binary restart file is produced for every 
! block and named as
!
! RESTARTOUTDIR/blkGLOBALBLKNUMBER.rst
!
! In addition the grid is described by
!
! RESTARTOUTDIR/octree.rst
!
! and an ASCII header file is produced with timestep and time info:
!
! RESTARTOUTDIR/restart.H
!
! The restart files are overwritten every time a new restart is done,
! but one can change the name of the RESTARTOUTDIR with the #RESTARTOUTDIR
! command from session to session. The default directory name is \'restartOUT\'.
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'PLOTDIR'},'content' => [{'type' => 'e','attrib' => {'type' => 'string','default' => 'GM/IO2','length' => '100','name' => 'NamePlotDir'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '

The NamePlotDir variable contains the name of the directory
where plot files and logfiles are saved relative to the run directory.
The directory should be inside the subdirectory with the name
of the component.

Default value is "GM/IO2".
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'SAVELOGFILE'},'content' => [{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'DoSaveLogfile'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'expr' => '$DoSaveLogfile'},'content' => [{'type' => 'e','attrib' => {'max' => '4','min' => '1','type' => 'strings','name' => 'StringLog'},'content' => [{'type' => 'e','attrib' => {'type' => 'string','required' => 'T','input' => 'select','name' => 'TypeLogVar'},'content' => [{'type' => 'e','attrib' => {'value' => 'MHD','default' => 'T','name' => 'MHD vars. dimensional'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'FLX','name' => 'Flux vars. dimensional'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'RAW','name' => 'Raw vars. dimensional'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'VAR','name' => 'Set vars. dimensional'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'mhd','default' => 'T','name' => 'MHD vars. scaled'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'flx','name' => 'Flux vars. scaled'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'raw','name' => 'Raw vars. scaled'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'var','name' => 'Set vars. scaled'},'content' => [],'name' => 'option'}],'name' => 'part'},{'type' => 'e','attrib' => {'type' => 'string','required' => 'F','input' => 'select','multiple' => 'T','name' => 'TypeTime'},'content' => [{'type' => 'e','attrib' => {'exclusive' => 'T','name' => 'none'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'step'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'date'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'time'},'content' => [],'name' => 'option'}],'name' => 'part'}],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '-1','type' => 'integer','default' => '1','name' => 'DnSaveLogfile'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '-1','type' => 'real','default' => '-1','name' => 'DtSaveLogfile'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'string','length' => '100','if' => '$TypeLogVar =~ /var/i','name' => 'NameLogVars'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'max' => '10','min' => '1','type' => 'strings','length' => '100','if' => '($TypeLogVar=~/flx/i or $NameLogVars=~/flx/i)','name' => 'StringLogRadii'},'content' => [{'type' => 'e','attrib' => {'min' => '$rBody','type' => 'real','multiple' => 'T','name' => 'LogRadii'},'content' => [],'name' => 'part'}],'name' => 'parameter'}],'name' => 'if'},{'type' => 't','content' => '
#SAVELOGFILE
T                       DoSaveLogfile, rest of parameters read if true
VAR step date           StringLog
100                     DnSaveLogfile
-1.                     DtSaveLogfile [sec]
rho p rhoflx            NameLogVars (read if StrigLog is \'var\' or \'VAR\')
4.0  10.0               rLog  (radii for the flux. Read if vars include \'flx\')

! Default is DoSaveLogfile=.false.
! The logfile can contain averages or point values and other scalar
! quantities.  It is written into an ASCII file named as
!
! NAMEPLOTDIR/log_TIMESTEP.log
!
! where NAMEPLOTDIR can be defined with the #PLOTDIR command (default is IO2).
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
! Ionosphere values:  cpcpn cpcps                  
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
'}],'name' => 'command'},{'type' => 'e','attrib' => {'if' => '$_IsFirstSession','name' => 'SATELLITE'},'content' => [{'type' => 'e','attrib' => {'min' => '0','type' => 'integer','default' => '0','name' => 'nSatellite'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'from' => '1','to' => '$nSatellite'},'content' => [{'type' => 'e','attrib' => {'max' => '5','min' => '1','type' => 'strings','name' => 'StringSatellite'},'content' => [{'type' => 'e','attrib' => {'type' => 'string','required' => 'T','input' => 'select','name' => 'TypeSatelliteVar'},'content' => [{'type' => 'e','attrib' => {'value' => 'MHD','default' => 'T','name' => 'MHD vars. dimensional'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'FUL','name' => 'All vars. dimensional'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'VAR','name' => 'Set vars. dimensional'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'mhd','name' => 'MHD vars. scaled'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'ful','name' => 'All vars. scaled'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'var','name' => 'Set vars. scaled'},'content' => [],'name' => 'option'}],'name' => 'part'},{'type' => 'e','attrib' => {'type' => 'string','required' => 'F','input' => 'select','name' => 'TypeTrajectory'},'content' => [{'type' => 'e','attrib' => {'default' => 'T','name' => 'file'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'eqn','name' => 'equation'},'content' => [],'name' => 'option'}],'name' => 'part'},{'type' => 'e','attrib' => {'type' => 'string','required' => 'F','input' => 'select','multiple' => 'T','name' => 'TypeTime'},'content' => [{'type' => 'e','attrib' => {'exclusive' => 'T','name' => 'none'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'step'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'date'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'time'},'content' => [],'name' => 'option'}],'name' => 'part'}],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '-1','type' => 'integer','default' => '1','name' => 'DnOutput'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '-1','type' => 'real','default' => '-1','name' => 'DtOutput'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'string','length' => '100','name' => 'NameTrajectoryFile'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'string','length' => '100','if' => '$TypeSatelliteVar =~ /\\bvar\\b/i','name' => 'NameSatelliteVars'},'content' => [],'name' => 'parameter'}],'name' => 'for'},{'type' => 't','content' => '
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
!
! satellitevar   = \'mhd\', \'ful\' or \'var\' (unitless output)
!                 \'MHD\', \'FUL\' or \'VAR\' (dimensional output)
! trajectorytype = \'file\' or \'eqn\'
! timetype       = \'none\', \'step\', \'time\', \'date\'
!
! The \'satellitevar\' part is required, 
! the \'trajectorytype\' part is optional (defaults to \'file\'), and
! the \'timetype\' part is also optional (default depends on satellitevar)
!
! The \'satellitevar\' string defines the variables to print in the satellite
! output file.  It also controls whether or not the variables will come out in
! dimensional or non-dimensional form by the capatilization of the
! satellitevars string: ALL CAPS means dimensional, all lower means 
! dimensionless. 
!
! If \'satellitevar\' is set to \'mhd\', the variables 
! \'rho ux uy uz bx by bz p jx jy jz\' will be saved, while\'ful\' implies
! \'rho ux uy uz bx by bz b1x b1y b1z p jx jy jz\'.
! If satellitevar is set to \'var\' then the list of variables is read 
! from the NameSatelliteVar parameter as a space separated list. 
! The choices for variables are currently:
!
! rho, rho, rhouy, rhouz, ux, uy, uz,
! Bx, By, Bz, B1x, B1y, B1z,
! E, P, Jx, Jy, Jz,
! theta1, theta2, phi1, phi2, status.
!
! If \'trajectorytype\' is \'file\' (default) than the trajectory of the 
! satellite is read from the file given by the NameTrajectoryFile parameter.
! If \'trajectorytype\' is \'eqn\' then the trajectory is defined by an
! equation, which is hard coded in subroutine satellite_trajectory_formula
! in satellites.f90.
!
! The \'timetype\' values mean the following:
!  none  = there will be no indication of time in the logfile 
!          (not even the number of steps),
!  step  = number of time steps (n_steps),
!  date  = time is given as an array of 7 integers:  year mo dy hr mn sc msc,
!  time  = time is given as a real number - elapsed time since the start of
!          the run.  Units are determined by satellitevar and unitUSER_t.
!
!  More than one \'timetype\' can be listed. They can be put together in any
!  combination.
!
! The DnOutput and DtOutput parameters determine the frequency of extracting
! values along the satellite trajectories. 
!
! The extracted satellite information is saved into the files named
!
! PLOTDIR/satellite_NN_TRAJECTORYNAME.sat
!
! where NN is the number of the satellite (e.g. 01), and TRAJECTORYNAME
! is the name of the trajectory file.
!
! The default is nSatellite=0, i.e. no satellite data is saved.
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'SAVEPLOT'},'content' => [{'type' => 'e','attrib' => {'max' => '100','min' => '0','type' => 'integer','default' => '0','name' => 'nPlotFile'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'from' => '1','to' => '$nPlotFile'},'content' => [{'type' => 'e','attrib' => {'max' => '3','min' => '3','type' => 'strings','name' => 'StringPlot'},'content' => [{'type' => 'e','attrib' => {'type' => 'string','required' => 'T','input' => 'select','name' => 'plotform'},'content' => [{'type' => 'e','attrib' => {'value' => 'tec','name' => 'TECPLOT'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'idl','name' => 'IDL'},'content' => [],'name' => 'option'}],'name' => 'part'},{'type' => 'e','attrib' => {'type' => 'string','required' => 'T','input' => 'select','name' => 'plotarea'},'content' => [{'type' => 'e','attrib' => {'value' => '3d'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'x=0'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'y=0','default' => 'T'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'z=0'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'sph'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'los'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'lin'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'cut','if' => '$plotform =~ /\\bidl\\b/'},'content' => [],'name' => 'option'}],'name' => 'part'},{'type' => 'e','attrib' => {'type' => 'string','required' => 'T','input' => 'select','name' => 'plotvar'},'content' => [{'type' => 'e','attrib' => {'value' => 'MHD','name' => 'MHD vars. dimensional'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'FUL','name' => 'All vars. dimensional'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'RAW','name' => 'Raw vars. dimensional'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'RAY','name' => 'Ray tracing vars. dim.'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'FLX','name' => 'Flux vars. dimensional'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'SOL','name' => 'Solar vars. dimensional'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'VAR','name' => 'Select dimensional vars.'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'mhd','name' => 'MHD vars. scaled'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'ful','name' => 'All vars. scaled'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'raw','name' => 'Raw vars. scaled'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'ray','name' => 'Ray tracing vars. scaled'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'flx','name' => 'Flux vars. scaled'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'sol','name' => 'Solar vars. scaled'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'var','name' => 'Select scaled vars.'},'content' => [],'name' => 'option'}],'name' => 'part'}],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '-1','type' => 'integer','name' => 'DnSavePlot'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '-1','type' => 'real','name' => 'DtSavePlot'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'expr' => '$plotarea =~ /\\bcut\\b/'},'content' => [{'type' => 'e','attrib' => {'type' => 'real','name' => 'xMinCut'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '$xMinCut','type' => 'real','name' => 'xMaxCut'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','name' => 'yMinCut'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '$yMinCut','type' => 'real','name' => 'yMaxCut'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','name' => 'zMinCut'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '$zMinCut','type' => 'real','name' => 'zMaxCut'},'content' => [],'name' => 'parameter'}],'name' => 'if'},{'type' => 'e','attrib' => {'min' => '0','type' => 'real','default' => '10','if' => '$plotarea =~ /\\bsph\\b/','name' => 'Radius'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'expr' => '$plotarea =~ /\\blos\\b/'},'content' => [{'type' => 'e','attrib' => {'type' => 'real','default' => '0','name' => 'LosVectorX'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '0.0001','name' => 'LosVectorY'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '1','name' => 'LosVectorZ'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '0','type' => 'real','default' => '20','name' => 'xSizeImage'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '0','type' => 'real','default' => '20','name' => 'ySizeImage'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '10','name' => 'xOffset'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '10','name' => 'yOffset'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '1','type' => 'real','default' => '2.5','name' => 'rOccult'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'max' => '1','min' => '0','type' => 'real','default' => '0.5','name' => 'MuLimbDarkening'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '2','type' => 'integer','default' => '200','name' => 'nPixX'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '2','type' => 'integer','default' => '200','name' => 'nPixY'},'content' => [],'name' => 'parameter'}],'name' => 'if'},{'type' => 'e','attrib' => {'expr' => '$plotarea =~ /\\blin\\b/'},'content' => [{'type' => 'e','attrib' => {'type' => 'string','input' => 'select','name' => 'NameLine'},'content' => [{'type' => 'e','attrib' => {'value' => 'A','name' => 'Advected B'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'default' => 'T','name' => 'B'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'U'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'J'},'content' => [],'name' => 'option'}],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'IsSingleLine'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'max' => '20','min' => '1','type' => 'integer','default' => '1','name' => 'nLine'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'from' => '1','to' => '$nLine'},'content' => [{'type' => 'e','attrib' => {'type' => 'real','name' => 'xStartLine'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','name' => 'yStartLine'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','name' => 'zStartLine'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'logical','name' => 'IsParallel'},'content' => [],'name' => 'parameter'}],'name' => 'for'}],'name' => 'if'},{'type' => 'e','attrib' => {'min' => '-1.0','type' => 'real','default' => '-1.0','if' => '($plotform=~/\\bidl\\b/ and $plotarea!~/\\b(sph|los|lin)\\b/)','name' => 'DxSavePlot'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'expr' => '$plotvar =~ /\\bvar\\b/i'},'content' => [{'type' => 'e','attrib' => {'type' => 'string','length' => '100','name' => 'NameVars'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'string','length' => '100','name' => 'NamePars'},'content' => [],'name' => 'parameter'}],'name' => 'if'}],'name' => 'for'},{'type' => 't','content' => '
#SAVEPLOT
6			nPlotfile
3d MHD tec		StringPlot ! 3d MHD data
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

! Default is nPlotFile=0
! StringPlot must contain the following 3 parts in arbitrary order
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
! DnSavePlot and for time accurate runs by DtSavePlot.
!
! Depending on StringPlot, further information is read from the parameter file
! in this order:
!
! PlotRange		if plotarea is \'cut\'
! DxSavePlot		if plotform is \'idl\' and plotarea is not sph, ion, los
! Radius		if plotarea is \'sph\'
! NameVars		if plotform is \'var\'
! NamePars		if plotform is \'var\'
!
! The PlotRange is described by 6 coordinates. If the width in one or two 
! dimensions is less than the smallest cell size within the plotarea, 
! then the plot file will be 2 or 1 dimensional. If the range is thin but
! symmetric about one of the x=0, y=0, or z=0 planes, data will be averaged
! in the postprocessing.
!
! Possible values for DxSavePlot (for IDL files):
!
!  0.5	- fixed resolution (any positive value)
!  0.	- fixed resolution based on the smallest cell in the plotting area
! -1.	- unstructured grid will be produced by PostIDL.exe
!
! Radius is the radius of the spherical cut for plotarea=\'sph\'
!
! LosVectorX,Y,Z define the direction of the line of sight integration
! xSizeImage, ySizeImage defines the size of the LOS image
! xOffset, yOffset defines the offset relative to the origin (Sun)
! rOccult defines the minimum distance of the line from the origin (Sun)
! MuLimbDarkening is the limb darkening parameter for the \'wl\' (white light)
!                 and \'pb\' (polarization brightness) plot variables.
!
! The possible values for NameVars with plotarea \'los\' 
!       are listed in subroutine set_plotvar_los in write_plot_los.f90.
! The possible values for NameVars for other plot areas
!       are listed in subroutine set_plotvar in write_plot_common.f90.
!
! The possible values for NamePars 
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

'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'SAVEBINARY'},'content' => [{'type' => 'e','attrib' => {'type' => 'logical','default' => 'T','name' => 'DoSaveBinary'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '
#SAVEBINARY
T			DoSaveBinary   used only for \'idl\' plot file

! Default is .true. Saves unformatted IO2/*.idl files if true. 
! This is the recommended method, because it is fast and accurate.
! The only advantage of saving IO2/*.idl in formatted text files is
! that it can be processed on another machine or with a different 
! (lower) precision. For example PostIDL.exe may be compiled with 
! single precision to make IO2/*.out files smaller, while BATSRUS.exe is 
! compiled in double precision, to make results more accurate.
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'SAVEPLOTSAMR'},'content' => [{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'DoSavePlotsAmr'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '
#SAVEPLOTSAMR
F			DoSavePlotsAmr

! Save plots before each AMR. Default is DoSavePlotsAMR=.false.
'}],'name' => 'command'}],'name' => 'commandgroup'},{'type' => 'e','attrib' => {'name' => 'AMR PARAMETERS'},'content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  AMR PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'type' => 'e','attrib' => {'if' => '$_IsFirstSession','name' => 'AMRINIT'},'content' => [{'type' => 'e','attrib' => {'type' => 'string','input' => 'select','name' => 'InitialRefineType'},'content' => [{'type' => 'e','attrib' => {'default' => '1','name' => 'default'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'all'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'none'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => '3Dbodyfocus'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'spherefocus'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'magnetosphere'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'points'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'helio_init'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'helio_z=4'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'all_then_focus'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'cme'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'points'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'mag_new'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'magnetosphere'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'magneto_fine'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'magneto12'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'magnetosaturn'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'magnetojupiter'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'paleo'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'comet'},'content' => [],'name' => 'option'}],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '0','type' => 'integer','default' => '4','name' => 'InitialRefineLevel'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '
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
'}],'name' => 'command'},{'type' => 'e','attrib' => {'if' => '$_IsFirstSession','name' => 'AMRINITPHYSICS'},'content' => [{'type' => 'e','attrib' => {'min' => '0','type' => 'integer','default' => '0','name' => 'nRefineLevelIC'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '
#AMRINITPHYSICS
3			nRefineLevelIC

! Defines number of physics (initial condition) based AMR-s AFTER the 
! geometry based initial AMR-s defined by #AMRINIT were done.
! Only useful if the initial condition has a non-trivial analytic form.
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'AMRLEVELS'},'content' => [{'type' => 'e','attrib' => {'min' => '-1','type' => 'integer','default' => '0','name' => 'MinBlockLevel'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '-1','type' => 'integer','default' => '99','name' => 'MaxBlockLevel'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'DoFixBodyLevel'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '
#AMRLEVELS
0			MinBlockLevel
99			MaxBlockLevel
F			DoFixBodyLevel

! Set the minimum/maximum levels that can be affected by AMR.  The usage is as
! follows:
!
! MinBlockLevel .ge.0 Cells can be coarsened up to the listed level but not
!                       further.
! MinBlockLevel .lt.0 The current grid is ``frozen\'\' for coarsening such that
!                       blocks are not allowed to be coarsened to a size
!                       larger than their current one.
! MaxBlockLevel .ge.0 Any cell at a level greater than or equal to
!                       MaxBlockLevel is uneffected by AMR (cannot be coarsened
!                       or refined).
! MaxBlockLevel .lt.0 The current grid is ``frozen\'\' for refinement such that
!                       blocks are not allowed to be refined to a size
!                       smaller than their current one.
! DoFixBodyLevel = T  Blocks touching the body cannot be coarsened or refined.
!
! This command has no effect when DoAutoRefine is .false. in the #AMR command.
!
! Note that the user can set either #AMRLEVELS or #AMRRESOLUTION but not
! both.  If both are set, the final one in the session will set the values
! for AMR.
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'AMRRESOLUTION'},'content' => [{'type' => 'e','attrib' => {'min' => '-1','type' => 'real','default' => '0','name' => 'DxCellMin'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '-1','type' => 'real','default' => '99999','name' => 'DxCellMax'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'DoFixBodyLevel'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '
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
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'AMR'},'content' => [{'type' => 'e','attrib' => {'min' => '-1','type' => 'integer','default' => '-1','name' => 'DnRefine'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'expr' => '$DnRefine>0'},'content' => [{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'DoAutoRefine'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'expr' => '$DoAutoRefine'},'content' => [{'type' => 'e','attrib' => {'max' => '100','min' => '0','type' => 'real','default' => '20','name' => 'PercentCoarsen'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'max' => '100','min' => '0','type' => 'real','default' => '20','name' => 'PercentRefine'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '1','type' => 'integer','default' => '99999','name' => 'MaxTotalBlocks'},'content' => [],'name' => 'parameter'}],'name' => 'if'}],'name' => 'if'},{'type' => 't','content' => '
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
! the PercentRefine and PercentCoarsen parameters. These per centages
! are approximate only, because the constraints of the block adaptive
! grid may result in more or fewer blocks than prescribed.
! The total number of blocks will not exceed the smaller of the 
! MaxTotalBlocks parameter and the total number of blocks available on all 
! the PE-s (which is determined by the number of PE-s and 
! the MaxBlocks parameter in ModSize.f90).
! 
! Default for DnRefine is -1, ie. no run time refinement.
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'AMRCRITERIA'},'content' => [{'type' => 'e','attrib' => {'type' => 'integer','input' => 'select','name' => 'nRefineCrit'},'content' => [{'type' => 'e','attrib' => {'name' => '1'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => '2'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'default' => '1','name' => '3'},'content' => [],'name' => 'option'}],'name' => 'parameter'},{'type' => 'e','attrib' => {'from' => '1','to' => '$nRefineCrit'},'content' => [{'type' => 'e','attrib' => {'type' => 'string','input' => 'select','name' => 'TypeRefine'},'content' => [{'type' => 'e','attrib' => {'value' => 'gradt/gradT','name' => 'grad T'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'gradp/gradP','name' => 'grad P'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'gradlogrho','name' => 'grad log(Rho)'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'gradlogP/gradlogp','name' => 'grad log(p)'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'gradE','name' => 'grad E'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'curlV/curlv/curlU/curlu','name' => 'curl U'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'curlB/curlb','name' => 'curl B'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'divU/divu/divV/divv','name' => 'div U'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'divb/divB','name' => 'divB'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'Valfven/vAlfven/valfven','name' => 'vAlfven'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'heliobeta','name' => 'heliospheric beta'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'flux'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'heliocurrentsheet','name' => 'heliospheric current sheet'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'rcurrents/Rcurrents','name' => 'rCurrents'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'transient/Transient','name' => 'Transient'},'content' => [],'name' => 'option'}],'name' => 'parameter'},{'type' => 'e','attrib' => {'expr' => '$TypeRefine =~ /transient/i'},'content' => [{'type' => 'e','attrib' => {'type' => 'string','input' => 'select','name' => 'TypeTransient'},'content' => [{'type' => 'e','attrib' => {'value' => 'p_dot/P_dot','name' => 'P_dot'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 't_dot/T_dot','name' => 'T_dot'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'rho_dot/Rho_dot','default' => 'T','name' => 'Rho_dot'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'RhoU_dot/rhou_dot','name' => 'RhoU_dot'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'Rho_2nd_1/rho_2nd_1','name' => 'Rho_2nd_1'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'Rho_2nd_2/rho_2nd_2','name' => 'Rho_2nd_2'},'content' => [],'name' => 'option'}],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'UseSunEarth'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'expr' => '$UseSunEarth'},'content' => [{'type' => 'e','attrib' => {'type' => 'real','name' => 'xEarth'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','name' => 'yEarth'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','name' => 'zEarth'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','name' => 'InvD2Ray'},'content' => [],'name' => 'parameter'}],'name' => 'if'}],'name' => 'if'}],'name' => 'for'},{'type' => 't','content' => '
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
! The possible choices for TypeTransient
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
! Example for InvD2Ray = 
!   1 - refine_profile = 0.3679 at distance Rsun/10 from the ray
!   2 - refine_profile = 0.0183 at distance Rsun/10 from the ray
!   3 - refine_profile = 0.0001 at distance Rsun/10 from the ray
'}],'name' => 'command'}],'name' => 'commandgroup'},{'type' => 'e','attrib' => {'name' => 'SCHEME PARAMETERS'},'content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  SCHEME PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'type' => 'e','attrib' => {'name' => 'SCHEME'},'content' => [{'type' => 'e','attrib' => {'type' => 'integer','input' => 'select','name' => 'nOrder'},'content' => [{'type' => 'e','attrib' => {'default' => 'T','name' => '1'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => '2'},'content' => [],'name' => 'option'}],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'string','input' => 'select','name' => 'TypeFlux'},'content' => [{'type' => 'e','attrib' => {'value' => 'Roe/roe/1','name' => 'Roe'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'Rusanov/rusanov/2/TVDLF','default' => 'T','name' => 'Rusanov'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'Linde/linde/3/HLLEL','name' => 'Linde'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'Sokolov/sokolov/4/AW','name' => 'Sokolov'},'content' => [],'name' => 'option'}],'name' => 'parameter'},{'type' => 'e','attrib' => {'expr' => '$nOrder == 2'},'content' => [{'type' => 'e','attrib' => {'type' => 'string','input' => 'select','name' => 'TypeLimiter'},'content' => [{'type' => 'e','attrib' => {'default' => 'T','name' => 'minmod'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'beta'},'content' => [],'name' => 'option'}],'name' => 'parameter'},{'type' => 'e','attrib' => {'max' => '2','min' => '1','type' => 'real','default' => '1.2','if' => '$TypeLimiter eq \'beta\'','name' => 'LimiterBeta'},'content' => [],'name' => 'parameter'}],'name' => 'if'},{'type' => 't','content' => '
#SCHEME
2			nOrder (1 or 2)
Rusanov			TypeFlux
minmod			TypeLimiter ! Only for nOrder=2
1.2			LimiterBeta ! Only for LimiterType=\'beta\'

! Default values are shown above.
!
! Possible values for TypeFlux:
! \'Rusanov\'     - Rusanov or Lax-Friedrichs flux     
! \'Linde        - Linde\'s HLLEL flux                   
! \'Sokolov\'     - Sokolov\'s Local Artificial Wind flux 
! \'Roe\'         - Roe\'s approximate Riemann flux       
!
! Possible values for TypeLimiter:
! \'minmod\'	- minmod limiter is the most robust 1D limiter
! \'beta\'        - Beta limiter
!
! Possible values for LimiterBeta are between 1.0 and 2.0 : 
!  LimiterBeta = 1.0 is the same as the minmod limiter
!  LimiterBeta = 2.0 is the same as the superbee limiter
!  LimiterBeta = 1.2 is the recommended value
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'NONCONSERVATIVE'},'content' => [{'type' => 'e','attrib' => {'type' => 'logical','default' => 'T','name' => 'UseNonConservative'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '
#NONCONSERVATIVE
T		UseNonConservative

! For Earth the default is using non-conservative equations 
! (close to the body).
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'CONSERVATIVECRITERIA'},'content' => [{'type' => 'e','attrib' => {'max' => '3','min' => '0','type' => 'integer','default' => '1','name' => 'nConservCrit'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'from' => '1','to' => '$nConservCrit'},'content' => [{'type' => 'e','attrib' => {'type' => 'string','input' => 'select','name' => 'TypeConservCrit'},'content' => [{'type' => 'e','attrib' => {'value' => 'r/R/radius/Radius','default' => 'T','name' => 'radius'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'p/P','name' => 'p'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'gradp/GradP','name' => 'grad P'},'content' => [],'name' => 'option'}],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '$rBody','type' => 'real','default' => '2*$rBody','if' => '$TypeConservCrit =~ /^r|radius$/i','name' => 'rConserv'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '0','type' => 'real','default' => '0.05','if' => '$TypeConservCrit =~ /^p$/i','name' => 'pCoeffConserv'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '0','type' => 'real','default' => '0.1','if' => '$TypeConservCrit =~ /gradp/i','name' => 'GradPCoeffConserv'},'content' => [],'name' => 'parameter'}],'name' => 'for'},{'type' => 't','content' => '
#CONSERVATIVECRITERIA
3		nConservCrit
r		TypeConservCrit
6.		rConserv             ! read if TypeConservCrit is \'r\'
p		TypeConservCrit
0.05		pCoeffConserv	     ! read if TypeConservCrit is \'p\'
GradP		TypeConservCrit
0.1		GradPCoeffConserv    ! read if TypeConservCrit is \'GradP\'

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
! Default values are nConservCrit = 1 with TypeConservCrit = \'r\'
! and rConserv=2*rBody, where rBody has a problem dependent default.
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'UPDATECHECK'},'content' => [{'type' => 'e','attrib' => {'type' => 'logical','default' => 'T','name' => 'UseUpdateCheck'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'expr' => '$UseUpdateCheck'},'content' => [{'type' => 'e','attrib' => {'max' => '100','min' => '0','type' => 'real','default' => '40','name' => 'RhoMinPercent'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '100','type' => 'real','default' => '400','name' => 'RhoMaxPercent'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'max' => '100','min' => '0','type' => 'real','default' => '40','name' => 'pMinPercent'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '100','type' => 'real','default' => '400','name' => 'pMaxPercent'},'content' => [],'name' => 'parameter'}],'name' => 'if'},{'type' => 't','content' => '
#UPDATECHECK
T			UseUpdateCheck
40.			RhoMinPercent
400.			RhoMaxPercent
40.			pMinPercent
400.			pMaxPercent

! Default values are shown.  This will adjust the timestep so that
! density and pressure cannot change by more than the given percentages
! in a single timestep.
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'PROLONGATION'},'content' => [{'type' => 'e','attrib' => {'type' => 'integer','input' => 'select','name' => 'nOrderProlong'},'content' => [{'type' => 'e','attrib' => {'default' => 'T','name' => '1'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => '2'},'content' => [],'name' => 'option'}],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'string','input' => 'select','if' => '$nOrderProlong==2','name' => 'TypeProlong'},'content' => [{'type' => 'e','attrib' => {'value' => 'lr','default' => 'T','name' => 'left-right'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'central'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'name' => 'minmod'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'lr2','name' => 'left-right extrapolate'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'central2','name' => 'central    extrapolate'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'minmod2','name' => 'minmod     extrapolate'},'content' => [],'name' => 'option'}],'name' => 'parameter'},{'type' => 't','content' => '
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
'}],'name' => 'command'},{'type' => 'e','attrib' => {'alias' => 'OPTIMIZE','name' => 'MESSAGEPASS'},'content' => [{'type' => 'e','attrib' => {'type' => 'string','input' => 'select','name' => 'TypeMessagePass'},'content' => [{'type' => 'e','attrib' => {'value' => 'allopt','default' => 'T','name' => 'm_p_cell FACES ONLY'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'all','name' => 'm_p_cell'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'opt','name' => 'm_p_dir FACES ONLY'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'dir','name' => 'm_p_dir group by directions'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'face','name' => 'm_p_dir group by faces     '},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'min','name' => 'm_p_dir group by kind and face'},'content' => [],'name' => 'option'}],'name' => 'parameter'},{'type' => 't','content' => '
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
! Constrained transport requires corners, default is \'all\'! 
! Diffusive control requires corners, default is \'all\'!
! Projection uses message_pass_dir for efficiency!
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'BORIS'},'content' => [{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'UseBorisCorrection'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'max' => '1','min' => '0','type' => 'real','default' => '1','if' => '$UseBorisCorrection','name' => 'BorisClightFactor'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '
#BORIS
T			UseBorisCorrection
1.0			BorisClightFactor !Only if UseBorisCorrection is true

! Default is boris_correction=.false.
! Use semi-relativistic MHD equations with speed of light reduced by
! the BorisClightFactor. Set BorisClightFactor=1.0 for true semi-relativistic
! MHD. Gives the same steady state as normal MHD analytically, but there
! can be differences due to discretization errors. 
! You can use either Boris or BorisSimple but not both. 
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'BORISSIMPLE'},'content' => [{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'UseBorisSimple'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'max' => '1','min' => '0','type' => 'real','default' => '1','if' => '$UseBorisSimple','name' => 'BorisClightFactor'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '
#BORISSIMPLE
T			UseBorisSimple
0.05			BorisClightFactor !Only if UseBorisSimple is true

! Default is UseBorisSimple=.false. 
! Use simplified semi-relativistic MHD with speed of light reduced by the
! BorisClightFactor. This is only useful with BorisClightFactor less than 1.
! Should give the same steady state as normal MHD, but there can be a
! difference due to discretization errors.
! You can use either Boris or BorisSimple but not both. 
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'DIVB'},'content' => [{'type' => 'e','attrib' => {'type' => 'logical','default' => 'T','name' => 'UseDivbSource'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'UseDivbDiffusion'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'UseProjection'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'UseConstrainB'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'expr' => '$UseDivbSource or $UseDivbDiffusion or $UseProjection or $UseConstrainB'},'content' => [{'type' => 't','content' => '
	! At least one of the options should be true.
	'}],'name' => 'rule'},{'type' => 'e','attrib' => {'expr' => 'not($UseProjection and ($UseDivbSource or $UseDivbDiffusion or $UseConstrainB))'},'content' => [{'type' => 't','content' => '
	! If UseProjection is true, all others should be false.
	'}],'name' => 'rule'},{'type' => 'e','attrib' => {'expr' => 'not($UseConstrainB and ($UseDivbSource or $UseDivbDiffusion or $UseProjection))'},'content' => [{'type' => 't','content' => '
	! If UseConstrainB is true, all others should be false.
	'}],'name' => 'rule'},{'type' => 't','content' => '
	
#DIVB
T			UseDivbSource
F			UseDivbDiffusion	
F			UseProjection           
F			UseConstrainB           

! Default values are shown above.
! If UseProjection is true, all others should be false.
! If UseConstrainB is true, all others should be false.
! At least one of the options should be true.
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'DIVBSOURCE'},'content' => [{'type' => 'e','attrib' => {'type' => 'logical','default' => 'T','name' => 'UseB0Source'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '
#DIVBSOURCE
T			UseB0Source

! Add extra source terms related to the non-zero divergence and curl of B0.
! Default is true.
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'DIVBDIFFUSION'},'content' => [{'type' => 'e','attrib' => {'max' => '1','min' => '0','type' => 'real','default' => '0.1666667','name' => 'DivbDiffCoeff'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '
#DIVBDIFFUSION
0.1666667		DivbDiffCoeff

! Default value is shown above. 1.0/6.0
! If divb_diffcoeff .gt. 0.5 then cfl .lt. 0.5/DivbDiffCoeff is required!
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'PROJECTION'},'content' => [{'type' => 'e','attrib' => {'type' => 'string','input' => 'select','name' => 'TypeProjectIter'},'content' => [{'type' => 'e','attrib' => {'value' => 'cg','default' => 'T','name' => 'Conjugate Gradients'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'bicgstab','name' => 'BiCGSTAB'},'content' => [],'name' => 'option'}],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'string','input' => 'select','name' => 'TypeProjectStop'},'content' => [{'type' => 'e','attrib' => {'value' => 'rel','default' => 'T','name' => 'Relative norm'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'max','name' => 'Maximum error'},'content' => [],'name' => 'option'}],'name' => 'parameter'},{'type' => 'e','attrib' => {'max' => '1','min' => '0','type' => 'real','default' => '0.1','name' => 'RelativeLimit'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '0','type' => 'real','default' => '0.0','name' => 'AbsoluteLimit'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '1','type' => 'integer','default' => '50','name' => 'MaxMatvec'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '
#PROJECTION
cg			TypeProjectIter:\'cg\' or \'bicgstab\' for iterative scheme
rel			TypeProjectStop:\'rel\' or \'max\' error for stop condition
0.1			RelativeLimit
0.0			AbsoluteLimit 
50			MaxMatvec (upper limit on matrix.vector multipl.)

! Default values are shown above.
!
! For symmetric Laplacian matrix TypeProjectIter=\'cg\' (Conjugate Gradients)
! should be used, as it is faster than BiCGSTAB. In current applications
! the Laplacian matrix is always symmetric.
! 
! The iterative scheme stops when the stopping condition is fulfilled:
!   TypeProjectStop = \'rel\' : 
!        stop if ||div B|| < RelativeLimit*||div B0||
!   TypeProjectStop = \'max\' and RelativeLimit is positive: 
!        stop if max(|div B|) < RelativeLimit*max(|div B0|)
!   TypeProjectStop = \'max\' and RelativeLimit is negative
!        stop if max(|div B|) < AbsoluteLimit
!
!   where ||.|| is the second norm, and B0 is the magnetic
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

'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'CORRECTP'},'content' => [{'type' => 'e','attrib' => {'max' => '1','min' => '0','type' => 'real','default' => '0.01','name' => 'pRatioLow'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'max' => '1','min' => '$pRatioLow','type' => 'real','default' => '0.1','name' => 'pRatioHigh'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '
! Default values are shown. 
!
! The purpose of the correctP subroutine is to remove any discrepancies between
! pressure in the p_BLK variable and the pressure calculated from the
! E_BLK variable. Such discrepancies can be caused by the
! constrained transport scheme and by the projection scheme
! which modify the magnetic energy. The algorithm is the following:
!
! q = eThermal/e
!
!                  q.lt. pRatioLow  E is set to eThermal+(rho*u**2+B**2)/2
! if pRatioLow .lt.q.lt.pRatioHigh  both P and E are modified depending on q
! if pratioHigh.lt.q                P is set to (gamma-1)*(e-(rho*u**2+B**2)/2)
!
! The 2nd case is a linear interpolation between the 2nd and 4th cases.
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'RAYTRACE'},'content' => [{'type' => 'e','attrib' => {'type' => 'logical','default' => 'T','name' => 'UseAccurateIntegral'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'UseAccurateTrace'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'max' => '60','min' => '0.01','type' => 'real','default' => '0.1','name' => 'DtExchangeRay'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '1','type' => 'integer','default' => '1','name' => 'DnRaytrace'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '

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
communitacation. The frequency of the information exchanges 
(in terms of CPU seconds) is given by the DtExchangeRay parameter. 
This is an optimization parameter for speed. Very small values of DtExchangeRay
result in many exchanges with few rays, while very large values result
in infrequent exchages thus some PE-s may become idle (no more work to do).
The optimal value is problem dependent. A typically acceptable value is 
DtExchangeRay = 0.1 seconds (default).

The DnRaytrace parameter contains the minimum number of iterations between
two ray tracings. The default value 1 means that every new step requires
a new trace (since the magnetic field is changing). A larger value implies
that the field does not change significantly in that many time steps.
The ray tracing is always redone if the grid changes due to an AMR.

Default values are UseAccurateIntegral = .true., UseAccurateTrace = .false.,
DtExchangeRay = 0.1 and DnRaytrace=1.
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'IM'},'content' => [{'type' => 'e','attrib' => {'min' => '0','type' => 'real','name' => 'TauCoupleIm'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '

#IM
0.01			TauCoupleIm

! Determine how fast the GM pressure p should be nudged towards the 
! IM pressure pIm. The nudging is done in every GM time step.
! The pIm pressure is updated every time IM->GM coupling occurs.
!
! If TauCoupleIm is less than 1.0, a weighted average is taken: 
!
! p\' = (p + TauCoupleIm*pIm)/(1+TauCoupleIm)
!
! Therefore the smaller TauCoupleIm is the slower the adjustment will be.
! It takes approximately 2/TauCoupleIm time steps to get p close to pIm.
!
! If TauCoupleIm is greater or equal than 1.0, the nudging is based on
! phyical time:
!
! p\' = p + max(1.0, dt/TauCoupleIm)*(pIm - p)
!
! where dt is the (local) time step. It takes about 2*TauCoupleIm seconds
! to get p close to pIm. If the (local) time step dt exceeds TauCoupleIm, 
! p\' = pIm is set in a single step.
!
! The default value is TauCoupleIm=0.01, which corresponds to a rather 
! slow nudging.
'}],'name' => 'command'}],'name' => 'commandgroup'},{'type' => 'e','attrib' => {'name' => 'PHYSICS PARAMETERS'},'content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  PHYSICS PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'type' => 'e','attrib' => {'if' => '$_IsFirstSession','name' => 'GAMMA'},'content' => [{'type' => 'e','attrib' => {'min' => '1','type' => 'real','default' => '1.6666666667','name' => 'Gamma'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '
#GAMMA
1.6666666667		Gamma

! The adiabatic index (ratio of the specific heats for fixed pressure
! and volume. The default value is 5.0/3.0, which is valid for
! monoatomic gas or plasma.
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'SHOCKTUBE'},'content' => [{'type' => 'e','attrib' => {'min' => '0','type' => 'real','default' => '1','name' => 'RhoLeft'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '0','name' => 'UnLeft'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '0','name' => 'Ut1Left'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '0','name' => 'Ut2Left'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '0.75','name' => 'BnLeft'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '1','name' => 'Bt1Left'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '0','name' => 'Bt2Left'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '0','type' => 'real','default' => '1','name' => 'pRight'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '0','type' => 'real','default' => '0.125','name' => 'RhoRight'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '0','name' => 'UnRight'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '0','name' => 'Ut1Right'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '0','name' => 'Ut2Right'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '0.75','name' => 'BnRight'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '-1','name' => 'Bt1Right'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '0','name' => 'Bt2Right'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '0','type' => 'real','default' => '0.1','name' => 'pRight'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','input' => 'select','name' => 'ShockSlope'},'content' => [{'type' => 'e','attrib' => {'value' => '0','default' => 'T','name' => 'no rotation'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => '0.25'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => '0.3333333333333'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => '0.5'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => '1'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => '2'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => '3'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => '4'},'content' => [],'name' => 'option'}],'name' => 'parameter'},{'type' => 't','content' => '
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
'}],'name' => 'command'},{'type' => 'e','attrib' => {'if' => '$_IsFirstSession','name' => 'SOLARWIND'},'content' => [{'type' => 'e','attrib' => {'min' => '0','type' => 'real','default' => '5','name' => 'SwRhoDim'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '0','type' => 'real','default' => '181712.175','name' => 'SwTDim'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'max' => '0','type' => 'real','default' => '-400','name' => 'SwUxDim'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '0','name' => 'SwUyDim'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '0','name' => 'SwUzDim'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '0','name' => 'SwBxDim'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '0','name' => 'SwByDim'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '5','name' => 'SwBzDim'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '
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
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'UPSTREAM_INPUT_FILE'},'content' => [{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'UseUpstreamInputFile'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'expr' => '$UseUpstreamInputFile'},'content' => [{'type' => 'e','attrib' => {'type' => 'string','length' => '100','name' => 'NameUpstreamFile'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '0','name' => 'SatelliteYPos'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '0','name' => 'SatelliteZPos'},'content' => [],'name' => 'parameter'}],'name' => 'if'},{'type' => 't','content' => '
#UPSTREAM_INPUT_FILE
T			UseUpstreamInputFile (rest of parameters read if true)
IMF.dat                 NameUpstreamFile

! Read IMF data from file NameUpstreamFile if UseUpstreamInputFile is true.
! The data file contains all information required for setting the upstream
! boundary conditions. Parameter TypeBcEast should be set to \'vary\' for
! the time dependent boundary condition.
!
! If the #SOLARWIND command is not provided than the first time read from
! the upstream input file will set the normalization of all variables
! in the GM component. Consequently either the #SOLARWIND command or
! the #UPSTREAM_INPUT_FILE command with UseUpstreamInputFile=.true.
! is required by the GM component.
!
! Default is UseUpstreamInputFile = .false.
'}],'name' => 'command'},{'type' => 'e','attrib' => {'alias' => 'MAGNETOSPHERE','if' => '$_IsFirstSession','name' => 'BODY'},'content' => [{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'UseBody'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'expr' => '$UseBody'},'content' => [{'type' => 'e','attrib' => {'min' => '0','type' => 'real','default' => '3','name' => 'rBody'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '-1','type' => 'real','default' => '4','name' => 'rCurrents'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '0','type' => 'real','default' => '1','name' => 'BodyRhoDim'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '0','type' => 'real','default' => '10000','name' => 'BodyTDim'},'content' => [],'name' => 'parameter'}],'name' => 'if'},{'type' => 't','content' => '
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
'}],'name' => 'command'},{'type' => 'e','attrib' => {'if' => '$_IsFirstSession','name' => 'GRAVITY'},'content' => [{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'UseGravity'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'integer','input' => 'select','if' => '$UseGravity','name' => 'iDirGravity'},'content' => [{'type' => 'e','attrib' => {'value' => '0','default' => 'T','name' => 'central mass'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => '1','name' => 'X direction'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => '2','name' => 'Y direction'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => '3','name' => 'Z direction'},'content' => [],'name' => 'option'}],'name' => 'parameter'},{'type' => 't','content' => '
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

'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'MASSLOADING'},'content' => [{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'UseMassLoading'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'DoAccelerateMassLoading'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '
#MASSLOADING
F			UseMassLoading
F			DoAccelerateMassLoading
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'HEATFLUX'},'content' => [{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'UseHeatFlux'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'logical','default' => 'T','name' => 'UseSpitzerForm'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'expr' => 'not $UseSpitzerForm'},'content' => [{'type' => 'e','attrib' => {'type' => 'real','default' => '1.23E-11','name' => 'Kappa0Heat'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '2.5','name' => 'Kappa0Heat'},'content' => [],'name' => 'parameter'}],'name' => 'if'},{'type' => 't','content' => '
#HEATFLUX
T		UseHeatFlux
F		UseSpitzerForm
1.23E-11	Kappa0Heat [W/m/K]	! Only if not UseSpitzerForm
2.50E+00	ExponentHeat [-]	! Only if not UseSpitzerForm
'}],'name' => 'command'},{'type' => 'e','attrib' => {'name' => 'RESISTIVEFLUX'},'content' => [{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'UseResistFlux'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'logical','default' => 'T','name' => 'UseSpitzerForm'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'expr' => 'not $UseSpitzerForm'},'content' => [{'type' => 'e','attrib' => {'type' => 'string','input' => 'select','name' => 'TypeResist'},'content' => [{'type' => 'e','attrib' => {'value' => 'Localized/localized','name' => 'localized'},'content' => [],'name' => 'option'},{'type' => 'e','attrib' => {'value' => 'Constant/constant','default' => 'T','name' => 'constant'},'content' => [],'name' => 'option'}],'name' => 'parameter'},{'type' => 'e','attrib' => {'expr' => '$TypeResist =~ /localized/i'},'content' => [{'type' => 'e','attrib' => {'type' => 'real','default' => '9.69953E+8','name' => 'Eta0Resist'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '150','name' => 'Alpha0Resist'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '0.5','name' => 'yShiftResist'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '0.05','name' => 'TimeInitRise'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '1','name' => 'TimeConstLev'},'content' => [],'name' => 'parameter'}],'name' => 'if'}],'name' => 'if'},{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'UseAnomResist'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'expr' => '$UseAnomResist'},'content' => [{'type' => 'e','attrib' => {'type' => 'real','default' => '1.93991E+09','name' => 'Eta0AnomResist'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '1.93991E+10','name' => 'EtaAnomMaxResist'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '1','name' => 'ThresholdFactorResist'},'content' => [],'name' => 'parameter'}],'name' => 'if'},{'type' => 't','content' => '
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
'}],'name' => 'command'},{'type' => 'e','attrib' => {'if' => '$_IsFirstSession','name' => 'TESTDISSMHD'},'content' => [{'type' => 'e','attrib' => {'type' => 'logical','default' => 'T','name' => 'UseDefaultUnits'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '2.635620E-02','name' => 'Grav0Diss'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '1.640000E-01','name' => 'Beta0Diss'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '1.500000E+06','name' => 'Length0Diss'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '1.159850E+01','name' => 'Time0Diss'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '5.019000E-11','name' => 'Rho0Diss'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '1.000000E+05','name' => 'Tem0Diss'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '6.000000E-01','name' => 'Theta0Diss'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '2.500000E+01','name' => 'Delta0Diss'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '7.000000E+00','name' => 'EpsilonDiss'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '4.500000E+00','name' => 'RhoDifDiss'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '4.000000E-01','name' => 'yShiftDiss'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '5.000000E-01','name' => 'ScaleHeightDiss'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '1.159850E+01','name' => 'ScaleFactorDiss'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '5.000000E-01','name' => 'BZ0iss'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '
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
'}],'name' => 'command'},{'type' => 'e','attrib' => {'if' => '$_IsFirstSession','name' => 'SECONDBODY'},'content' => [{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'UseBody2'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'expr' => '$UseBody2'},'content' => [{'type' => 'e','attrib' => {'min' => '0','type' => 'real','default' => '0.1','name' => 'rBody2'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'max' => '$xMax','min' => '$xMin','type' => 'real','default' => '-40','name' => 'xBody2'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'max' => '$yMax','min' => '$yMin','type' => 'real','default' => '0','name' => 'yBody2'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'max' => '$zMax','min' => '$zMin','type' => 'real','default' => '0','name' => 'zBody2'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '$rBody2','type' => 'real','default' => '1.3*$rBody2','name' => 'rCurrents2'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '0','type' => 'real','default' => '5','name' => 'RhoDimBody2'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '0','type' => 'real','default' => '25000','name' => 'tDimBody2'},'content' => [],'name' => 'parameter'}],'name' => 'if'},{'type' => 't','content' => '

#SECONDBODY
T			UseBody2   ! Rest of the parameters read if .true.
0.01			rBody2 
-40.			xBody2
0.			yBody2
0.			zBody2
0.011                   rCurrents2  !This is unused currently 
5.0			RhoDimBody2 (/ccm) density for fixed BC for rho_BLK
25000.0			TDimBody2 (K) temperature for fixed BC for P_BLK

! Default for UseBody2=.false.   -   All others no defaults!
! This command should appear before the #INNERBOUNDARY command when using
! a second body.
'}],'name' => 'command'},{'type' => 'e','attrib' => {'if' => '$_IsFirstSession','name' => 'DIPOLEBODY2'},'content' => [{'type' => 'e','attrib' => {'type' => 'real','name' => 'BdpDimBody2x'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','name' => 'BdpDimBody2y'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','name' => 'BdpDimBody2z'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '

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

!for now the dipole of the second body can only be aligned with the z-axis
'}],'name' => 'command'}],'name' => 'commandgroup'},{'type' => 'e','attrib' => {'name' => 'SOLAR PROBLEM TYPES'},'content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!! SOLAR PROBLEM TYPES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'type' => 'e','attrib' => {'if' => '$_IsFirstSession and $_NameComp ne \'GM\'','name' => 'HELIOSPHERE'},'content' => [{'type' => 'e','attrib' => {'min' => '0','type' => 'real','default' => '2.85E06','name' => 'BodyTDim'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '0','type' => 'real','default' => '1.50E8','name' => 'BodyRhoDim'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '0','type' => 'real','default' => '25.0','name' => 'qSun'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '0','type' => 'real','default' => '1.75','name' => 'tHeat'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '0','type' => 'real','default' => '1.0','name' => 'rHeat'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '0','type' => 'real','default' => '4.5','name' => 'SigmaHeat'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'DoInitRope'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'expr' => '$DoInitRope'},'content' => [{'type' => 'e','attrib' => {'min' => '0','type' => 'real','default' => '0.7','name' => 'CmeA'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '0','type' => 'real','default' => '1.2','name' => 'CmeR1'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '0','type' => 'real','default' => '1.0','name' => 'CmeR0'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '0','type' => 'real','default' => '0.23','name' => 'CmeA1'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '0.0','name' => 'CmeAlpha'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '0','type' => 'real','default' => '2.5E-12','name' => 'CmeRho1'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '0','type' => 'real','default' => '2.0E-13','name' => 'CmeRho2'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'max' => '10','min' => '0','type' => 'real','default' => '0.0','name' => 'ModulationRho'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'max' => '10','min' => '0','type' => 'real','default' => '0.0','name' => 'ModulationP'},'content' => [],'name' => 'parameter'}],'name' => 'if'},{'type' => 't','content' => '
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
0.0			cRotxGl98 [deg]
0.0			cRotYGl98 [deg]
0.0			cRotZGl98 [deg]

! This command defines the heliosphere parameters with a CME model.
'}],'name' => 'command'},{'type' => 'e','attrib' => {'if' => '$_NameComp ne \'GM\'','name' => 'HELIOUPDATEB0'},'content' => [{'type' => 'e','attrib' => {'min' => '-1','type' => 'real','default' => '0.0001','name' => 'DtUpdateB0'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '

#HELIOUPDATEB0
-1.0			DtUpdateB0 [s]

Set the frequency of updating the B0 field for the solar corona.
A negative value means that the B0 field is not updated.
'}],'name' => 'command'},{'type' => 'e','attrib' => {'if' => '$_NameComp ne \'GM\'','name' => 'HELIODIPOLE'},'content' => [{'type' => 'e','attrib' => {'type' => 'real','name' => 'HelioDipoleStrength'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'max' => '90','min' => '-90','type' => 'real','default' => '0','name' => 'HelioDipoleTilt'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '

#HELIODIPOLE
-3.0                    HelioDipoleStrength [G]
 0.0                    HelioDipoleTilt     [deg]

! Variable HelioDipoleStrength defines the equatorial field strength in Gauss,
! while HelioDipoleTilt is the tilt relative to the ecliptic North 
! (negative sign means towards the planet) in degrees.
!
! Default value is HelioDipoleStrength = 0.0.
'}],'name' => 'command'},{'type' => 'e','attrib' => {'alias' => 'INERTIAL','if' => '$_IsFirstSession and $_NameComp ne \'GM\'','name' => 'HELIOROTATION'},'content' => [{'type' => 'e','attrib' => {'type' => 'logical','default' => 'T','name' => 'UseInertialFrame'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'expr' => '$UseInertialFrame'},'content' => [{'type' => 'e','attrib' => {'type' => 'logical','default' => 'T','name' => 'UseRotatingBC'},'content' => [],'name' => 'parameter'}],'name' => 'if'},{'type' => 't','content' => '

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
'}],'name' => 'command'},{'type' => 'e','attrib' => {'if' => '$_NameComp ne \'GM\'','name' => 'HELIOTEST'},'content' => [{'type' => 'e','attrib' => {'type' => 'logical','default' => 'F','name' => 'DoSendMHD'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '

#HELIOTEST
F			DoSendMHD

! If DoSendMHD is true, IH sends the real MHD solution to GM in the coupling.
! If DoSendMHD is false then the values read from the IMF file are sent,
! so there is no real coupling. Mostly used for testing the framework.
!
! Default value is true, ie. real coupling.
'}],'name' => 'command'},{'type' => 'e','attrib' => {'if' => '$_IsFirstSession and $_NameComp ne \'GM\'','name' => 'CME'},'content' => [{'type' => 'e','attrib' => {'type' => 'string','input' => 'select','name' => 'TypeCme'},'content' => [{'type' => 'e','attrib' => {'default' => 'T','name' => 'Low'},'content' => [],'name' => 'option'}],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '0','type' => 'real','default' => '0.7','name' => 'CmeA'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '0','type' => 'real','default' => '1.2','name' => 'CmeR1'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '0','type' => 'real','default' => '1.0','name' => 'CmeR0'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '0','type' => 'real','default' => '0.23','name' => 'CmeA1'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '0.0','name' => 'CmeAlpha'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '0','type' => 'real','default' => '2.5E-12','name' => 'CmeRho1'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '0','type' => 'real','default' => '2.0E-13','name' => 'CmeRho2'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '1.0','name' => 'CmeB1Dim'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '0','type' => 'real','default' => '4.0E5','name' => 'CmeUErupt'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '
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
'}],'name' => 'command'},{'type' => 'e','attrib' => {'if' => '$_IsFirstSession and $_NameComp ne \'GM\'','name' => 'ARCADE'},'content' => [{'type' => 'e','attrib' => {'min' => '0','type' => 'real','default' => '1.0E6','name' => 'tArcDim'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '0','type' => 'real','default' => '1.0E-12','name' => 'RhoArcDim'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '0','type' => 'real','default' => '0.718144','name' => 'bArcDim'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '1.0E6','name' => 'ByArcDim'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '5.0E3','name' => 'UzArcDim'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '0.5','name' => 'Phi0Arc'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'type' => 'real','default' => '1.3','name' => 'MuArc'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '0','type' => 'real','default' => '3','name' => 'ExpArc'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '0','type' => 'real','default' => '0.5','name' => 'WidthArc'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '
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
'}],'name' => 'command'}],'name' => 'commandgroup'},{'type' => 'e','attrib' => {'name' => 'COMET PROBLEM TYPE'},'content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!! COMET PROBLEM TYPE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'type' => 'e','attrib' => {'if' => '$_IsFirstSession','name' => 'COMET'},'content' => [{'type' => 'e','attrib' => {'min' => '0','type' => 'real','name' => 'ProdRate'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '0','type' => 'real','name' => 'UrNeutral'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '0','type' => 'real','name' => 'AverageMass'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '0','type' => 'real','name' => 'IonizationRate'},'content' => [],'name' => 'parameter'},{'type' => 'e','attrib' => {'min' => '0','type' => 'real','name' => 'kFriction'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '
#COMET
1.0E28		ProdRate    - Production rate (#/s)
1.0		UrNeutral   - neutral radial outflow velocity (km/s)
17.0		AverageMass - average particle mass (amu)
1.0E-6		IonizationRate (1/s)
1.7E-9		kFriction - ion-neutral friction rate coefficient (cm^3/s)

! Only used by problem_comet.  Defaults are as shown.
'}],'name' => 'command'}],'name' => 'commandgroup'},{'type' => 'e','attrib' => {'name' => 'SCRIPT COMMANDS'},'content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! SCRIPT COMMANDS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'type' => 'e','attrib' => {'name' => 'INCLUDE'},'content' => [{'type' => 'e','attrib' => {'type' => 'string','default' => 'Param/','length' => '100','name' => 'NameIncludeFile'},'content' => [],'name' => 'parameter'},{'type' => 't','content' => '

#INCLUDE
Param/SSS_3000		NameIncludeFile

! Include a library file from Param/ or any file from anywhere else.
'}],'name' => 'command'},{'type' => 'e','attrib' => {'if' => '$_IsStandAlone','name' => 'RUN'},'content' => [{'type' => 't','content' => '

#RUN

This command is only used in stand alone mode.

Run BATSRUS with the parameters above and then return for the next session
'}],'name' => 'command'},{'type' => 'e','attrib' => {'if' => '$_IsStandAlone','name' => 'END'},'content' => [{'type' => 't','content' => '

#END

This command is only used in stand alone mode.

Run the executable with the parameters above and then stop.
In included files #END simply means the end of the included lines.
'}],'name' => 'command'}],'name' => 'commandgroup'},{'type' => 'e','attrib' => {'expr' => '($SwRhoDim > 0) or $UseUpstreamInputFile or $_NameComp ne \'GM\''},'content' => [{'type' => 't','content' => '
	Either command #SOLARWIND or #UPSTREAM_INPUT_FILE must be used!
'}],'name' => 'rule'},{'type' => 'e','attrib' => {'expr' => '$MaxImplBlock>1 or not $UsePartImplicit or not $MaxImplBlock'},'content' => [{'type' => 't','content' => '
	Part implicit scheme requires more than 1 implicit block!
'}],'name' => 'rule'},{'type' => 'e','attrib' => {'expr' => '$MaxImplBlock==$MaxBlock or not $UseFullImplicit'},'content' => [{'type' => 't','content' => '
	Full implicit scheme should be used with equal number of 
	explicit and implicit blocks!
'}],'name' => 'rule'}],'name' => 'commandList'}];