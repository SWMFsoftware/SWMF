#^CFG FILE _FALSE_
$tree = [{'type' => 'e','name' => 'commandList','content' => [{'type' => 't','content' => '

List of MH (GM, IH and SC) commands used in the PARAM.in file




'},{'type' => 'e','name' => 'set','content' => [],'attrib' => {'type' => 'integer','name' => 'nI','value' => '$_GridSize[0]'}},{'type' => 'e','name' => 'set','content' => [],'attrib' => {'type' => 'integer','name' => 'nJ','value' => '$_GridSize[1]'}},{'type' => 'e','name' => 'set','content' => [],'attrib' => {'type' => 'integer','name' => 'nK','value' => '$_GridSize[2]'}},{'type' => 'e','name' => 'set','content' => [],'attrib' => {'type' => 'integer','name' => 'MaxBlock','value' => '$_GridSize[3]'}},{'type' => 'e','name' => 'set','content' => [],'attrib' => {'type' => 'integer','name' => 'MaxImplBlock','value' => '$_GridSize[4]'}},{'type' => 'e','name' => 'set','content' => [],'attrib' => {'type' => 'integer','name' => 'MaxBlockALL','value' => '$_nProc and $MaxBlock and $_nProc*$MaxBlock'}},{'type' => 'e','name' => 'commandgroup','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! STAND ALONE PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

'},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseNewParam','default' => 'T'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseNewAxes','default' => 'T'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'DoTimeAccurate','default' => 'T'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseCorotation','default' => 'T'}},{'type' => 't','content' => '

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
'}],'attrib' => {'name' => 'NEWPARAM','if' => '$_IsStandAlone'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'string','name' => 'StringDescription','length' => '100'}},{'type' => 't','content' => '

#DESCRIPTION
This is a test run for Jupiter with no rotation.

This command is only used in the stand alone mode.

The StringDescription string can be used to describe the simulation
for which the parameter file is written. The #DESCRIPTION command and
the StringDescription string are saved into the restart file,
which helps in identifying the restart files.

The default value is "Please describe me!", which is self explanatory.
'}],'attrib' => {'name' => 'DESCRIPTION','if' => '$_IsStandAlone'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'DoEcho','default' => 'F'}},{'type' => 't','content' => '

#ECHO
T                       DoEcho

This command is only used in the stand alone mode.

If the DoEcho variable is true, the input parameters are echoed back.
The default value for DoEcho is .false., but it is a good idea to
set it to true at the beginning of the PARAM.in file.
'}],'attrib' => {'name' => 'ECHO','if' => '$_IsStandAlone'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '-1','type' => 'integer','name' => 'DnProgressShort','default' => '10'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '-1','type' => 'integer','name' => 'DnProgressLong','default' => '100'}},{'type' => 't','content' => '
#PROGRESS
10			DnProgressShort
100			DnProgressLong

The frequency of short and long progress reports for BATSRUS in
stand alone mode. These are the defaults. Set -1-s for no progress reports.
'}],'attrib' => {'name' => 'PROGRESS','if' => '$_IsStandAlone'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'DoTimeAccurate','default' => 'T'}},{'type' => 't','content' => '

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
'}],'attrib' => {'name' => 'TIMEACCURATE','if' => '$_IsStandAlone'}},{'type' => 'e','name' => 'command','content' => [{'type' => 't','content' => '

This command is allowed in stand alone mode only for sake of the 
test suite, which contains these commands when the framework is tested.
'}],'attrib' => {'multiple' => 'T','name' => 'BEGIN_COMP','if' => '$_IsStandAlone'}},{'type' => 'e','name' => 'command','content' => [{'type' => 't','content' => '

This command is allowed in stand alone mode only for sake of the 
test suite, which contains these commands when the framework is tested.
'}],'attrib' => {'multiple' => 'T','name' => 'END_COMP','if' => '$_IsStandAlone'}},{'type' => 'e','name' => 'command','content' => [{'type' => 't','content' => '

#RUN

This command is only used in stand alone mode.

The #RUN command does not have any parameters. It signals the end
of the current session, and makes BATSRUS execute the session with
the current set of parameters. The parameters for the next session
start after the #RUN command. For the last session there is no
need to use the #RUN command, since the #END command or simply
the end of the PARAM.in file makes BATSRUS execute the last session.
'}],'attrib' => {'name' => 'RUN','if' => '$_IsStandAlone'}},{'type' => 'e','name' => 'command','content' => [{'type' => 't','content' => '

#END

The #END command signals the end of the included file or the
end of the PARAM.in file. Lines following the #END command are
ignored. It is not required to use the #END command. The end
of the included file or PARAM.in file is equivalent with an 
#END command in the last line.
'}],'attrib' => {'name' => 'END'}}],'attrib' => {'name' => 'STAND ALONE MODE'}},{'type' => 'e','name' => 'commandgroup','content' => [{'type' => 't','content' => '
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

'},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Earth','value' => 'EARTH/Earth/earth','default' => 'T'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Saturn','value' => 'SATURN/Saturn/saturn'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'New'}}],'attrib' => {'type' => 'string','name' => 'NamePlanet','input' => 'select'}},{'type' => 'e','name' => 'if','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'RadiusPlanet'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'MassPlanet'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'OmegaPlanet'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'TiltRotation'}},{'type' => 'e','name' => 'parameter','content' => [{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'NONE'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'DIPOLE','default' => 'T'}}],'attrib' => {'type' => 'string','name' => 'TypeBField','input' => 'select'}}],'attrib' => {'expr' => '$NamePlanet eq \'New\''}},{'type' => 'e','name' => 'if','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','max' => '180','name' => 'MagAxisThetaGeo'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','max' => '360','name' => 'MagAxisPhiGeo'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'DipoleStrength'}}],'attrib' => {'expr' => '$TyepBField eq \'DIPOLE\''}},{'type' => 'e','name' => 'rule','content' => [{'type' => 't','content' => '
		PLANET should precede $PlanetCommand
	'}],'attrib' => {'expr' => 'not $PlanetCommand'}},{'type' => 't','content' => '

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
'}],'attrib' => {'name' => 'PLANET','if' => '$_IsFirstSession and $_IsStandAlone'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'IsRotAxisPrimary','default' => 'T'}},{'type' => 'e','name' => 'if','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','max' => '180','name' => 'RotAxisTheta'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','max' => '360','name' => 'RotAxisPhi'}}],'attrib' => {'expr' => '$IsRotAxisPrimary'}},{'type' => 'e','name' => 'set','content' => [],'attrib' => {'type' => 'string','name' => 'PlanetCommand','value' => 'ROTATIONAXIS'}},{'type' => 't','content' => '

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
'}],'attrib' => {'name' => 'ROTATIONAXIS','if' => '$_IsStandAlone'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseRotation','default' => 'T'}},{'type' => 'e','name' => 'if','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'RotationPeriod'}}],'attrib' => {'expr' => '$UseRotation'}},{'type' => 'e','name' => 'set','content' => [],'attrib' => {'type' => 'string','name' => 'PlanetCommand','value' => 'MAGNETICAXIS'}},{'type' => 't','content' => '

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
'}],'attrib' => {'name' => 'ROTATION','if' => '$_IsStandAlone'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'IsMagAxisPrimary','default' => 'T'}},{'type' => 'e','name' => 'if','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','max' => '180','name' => 'MagAxisTheta'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','max' => '360','name' => 'MagAxisPhi'}}],'attrib' => {'expr' => '$IsMagAxisPrimary'}},{'type' => 'e','name' => 'set','content' => [],'attrib' => {'type' => 'string','name' => 'PlanetCommand','value' => 'MAGNETICAXIS'}},{'type' => 't','content' => '

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
'}],'attrib' => {'name' => 'MAGNETICAXIS','if' => '$_IsStandAlone'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'DipoleStrength'}},{'type' => 't','content' => '

#DIPOLE
-3.11e-4		DipoleStrength [Tesla]

The DipoleStrength variable contains the
magnetic equatorial strength of the dipole magnetic field in Tesla.

The default value is the real dipole strength for the planet.
For the Earth the default is taken to be -31100 nT.
The sign is taken to be negative so that the magnetic axis can
point northward as usual.
'}],'attrib' => {'name' => 'DIPOLE','if' => '$_IsStandAlone'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '-1','type' => 'real','name' => 'DtUpdateB0','default' => '0.0001'}},{'type' => 't','content' => '

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
'}],'attrib' => {'name' => 'UPDATEB0','if' => '$_IsStandAlone'}},{'type' => 'e','name' => 'command','content' => [{'type' => 't','content' => '

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
'}],'attrib' => {'name' => 'IDEALAXES','if' => '$_IsStandAlone'}}],'attrib' => {'name' => 'PLANET COMMANDS'}},{'type' => 'e','name' => 'commandgroup','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!  USER DEFINED INPUT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

'},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseUserInnerBcs','default' => 'F'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseUserSource','default' => 'F'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseUserPerturbation','default' => 'F'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseUserOuterBcs','default' => 'F'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseUserICs','default' => 'F'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseUserSpecifyRefinement','default' => 'F'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseUserLogFiles','default' => 'F'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseUserWritePlot','default' => 'F'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseUserAMR','default' => 'F'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseUserEchoInput','default' => 'F'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseUserB0','default' => 'F'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseUserSetPhysConst','default' => 'F'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseUserUpdateStates','default' => 'F'}},{'type' => 't','content' => '

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
'}],'attrib' => {'name' => 'USER_FLAGS'}},{'type' => 'e','name' => 'command','content' => [{'type' => 't','content' => '

This command signals the beginning of the section of the file which 
is read by the subroutine user\\_read\\_inputs in the user\\_routines.f90 file.
The section ends with the #USERINPUTEND command. There is no XML based parameter
checking in the user section.
'}],'attrib' => {'name' => 'USERINPUTBEGIN'}},{'type' => 'e','name' => 'command','content' => [{'type' => 't','content' => '

This command signals the end of the section of the file which 
is read by the subroutine user\\_read\\_inputs in the user\\_routines.f90 file.
The section begins with the #USERINPUTBEGIN command. There is no XML based parameter
checking in the user section.
'}],'attrib' => {'name' => 'USERINPUTEND'}}],'attrib' => {'name' => 'USER DEFINED INPUT'}},{'type' => 'e','name' => 'commandgroup','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  TESTING AND TIMING PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'string','name' => 'TestString','length' => '100'}},{'type' => 't','content' => '
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
'}],'attrib' => {'name' => 'TEST'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '-2','type' => 'integer','max' => '$nI+2','name' => 'iTest'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '-2','type' => 'integer','max' => '$nJ+2','name' => 'jTest'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '-2','type' => 'integer','max' => '$nK+2','name' => 'kTest'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '1','type' => 'integer','max' => '$MaxBlock','name' => 'iBlockTest'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'integer','name' => 'iProcTest'}},{'type' => 't','content' => '
#TESTIJK
1                       iTest           (cell index for testing)
1                       jTest           (cell index for testing)
1                       kTest           (cell index for testing)
1                       iBlockTest      (block index for testing)
0                       iProcTest       (processor index for testing)

! The location of test info in terms of indices, block and processor number.
! Note that the user should set #TESTIJK or #TESTXYZ, not both.  If both
! are set, the final one in the session will set the test point.
'}],'attrib' => {'name' => 'TESTIJK'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '$xMin','type' => 'real','max' => '$xMax','name' => 'xTest'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '$yMin','type' => 'real','max' => '$yMax','name' => 'yTest'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '$zMin','type' => 'real','max' => '$zMax','name' => 'zTest'}},{'type' => 't','content' => '
#TESTXYZ
1.5                     xTest           (X coordinate of cell for testing)
-10.5                   yTest           (Y coordinate of cell for testing)
-10.                    zTest           (Z coordinate of cell for testing)

! The location of test info in terms of coordinates.
! Note that the user should set #TESTIJK or #TESTXYZ, not both.  If both
! are set, the final one in the session will set the test point.
'}],'attrib' => {'name' => 'TESTXYZ'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '-1','type' => 'integer','name' => 'nIterTest','default' => '-1'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '-1','type' => 'real','name' => 'TimeTest','default' => '1e30'}},{'type' => 't','content' => '

#TESTTIME
-1                      nIterTest       (iteration number to start testing)
10.5                    TimeTest        (time to start testing in seconds)

! The time step and physical time to start testing.
'}],'attrib' => {'name' => 'TESTTIME'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Rho','value' => '1','default' => 'T'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'RhoUx','value' => '2'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'RhoUy','value' => '3'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'RhoUz','value' => '4'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Bx','value' => '5'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'By','value' => '6'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Bz','value' => '7'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'e','value' => '8'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'p','value' => '9'}}],'attrib' => {'type' => 'integer','name' => 'iVarTest','input' => 'select'}},{'type' => 't','content' => '
#TESTVAR
1                       iVarTest

! Index of variable to be tested. Default is rho_="1", ie. density.
'}],'attrib' => {'name' => 'TESTVAR'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'all','value' => '0'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'x','value' => '1','default' => 'T'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'y','value' => '2'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'z','value' => '3'}}],'attrib' => {'type' => 'integer','name' => 'iVarTest','input' => 'select'}},{'type' => 't','content' => '
#TESTDIM
1                       iDimTest

! Index of dimension/direction to be tested. Default is X dimension.
'}],'attrib' => {'name' => 'TESTDIM'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseStrict','default' => 'T'}},{'type' => 't','content' => '
#STRICT
T                       UseStrict

! If true then stop when parameters are incompatible. If false, try to
! correct parameters and continue. Default is true, ie. strict mode
'}],'attrib' => {'name' => 'STRICT'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'errors and warnings only','value' => '-1'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'start and end of sessions','value' => '0'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'normal','value' => '1','default' => 'T'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'calls on test processor','value' => '10'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'calls on all processors','value' => '100'}}],'attrib' => {'type' => 'integer','name' => 'lVerbose','input' => 'select'}},{'type' => 't','content' => '
#VERBOSE
-1                      lVerbose

! Verbosity level controls the amount of output to STDOUT. Default level is 1.
!   lVerbose .le. -1 only warnings and error messages are shown.
!   lVerbose .ge.  0 start and end of sessions is shown.
!   lVerbose .ge.  1 a lot of extra information is given.
!   lVerbose .ge. 10 all calls of set_oktest are shown for the test processor.
!   lVerbose .ge.100 all calls of set_oktest are shown for all processors.
'}],'attrib' => {'name' => 'VERBOSE'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'DoDebug','default' => 'F'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'DoDebugGhost','default' => 'F'}},{'type' => 't','content' => '
#DEBUG
F                       DoDebug         (use it as if(okdebug.and.oktest)...)
F                       DoDebugGhost    (parameter for show_BLK in library.f90)

! Excessive debug output can be controlled by the global okdebug parameter
'}],'attrib' => {'name' => 'DEBUG'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'CodeVersion','default' => '7.50'}},{'type' => 't','content' => '
#CODEVERSION
7.50                    CodeVersion

! Checks CodeVersion. Prints a WARNING if it differs from the CodeVersion
! defined in ModMain.f90. Used in newer restart header files. 
! Should be given in PARAM.in when reading old restart files, 
! which do not have version info in the header file.
'}],'attrib' => {'name' => 'CODEVERSION','if' => '$_IsFirstSession'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'string','name' => 'NameEquation','length' => '100','default' => 'MHD'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'integer','name' => 'nVar','default' => '8'}},{'type' => 't','content' => '
#EQUATION
MHD			NameEquation
8			nVar

! Define the equation name and the number of variables.
! If any of these do not agree with the values determined 
! by the code, BATSRUS stops with an error. Used in restart
! header files and can be given in PARAM.in as a check
! and as a description.
'}],'attrib' => {'name' => 'EQUATION','if' => '$_IsFirstSession'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'single precision (4)','value' => '4','default' => '$_nByteReal==4'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'double precision (8)','value' => '8','default' => '$_nByteReal==8'}}],'attrib' => {'type' => 'integer','name' => 'nByteReal','input' => 'select'}},{'type' => 'e','name' => 'rule','content' => [{'type' => 't','content' => '
		nByteReal in file must agree with _nByteReal.
	'}],'attrib' => {'expr' => '$nByteReal==$_nByteReal'}},{'type' => 't','content' => '

#PRECISION
8                       nByteReal

! Define the number of bytes in a real number. If it does not agree
! with the value determined by the code, BATSRUS stops with an error.
! This is a check, the internal value is calculated in parallel_setup.
! Used in latest restart header files to check binary compatibility.
! May be given in PARAM.in to enforce a certain precision.
'}],'attrib' => {'name' => 'PRECISION','if' => '$_IsFirstSession'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '$nI','type' => 'integer','max' => '$nI','name' => 'nI','default' => '$nI'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '$nJ','type' => 'integer','max' => '$nJ','name' => 'nJ','default' => '$nJ'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '$nK','type' => 'integer','max' => '$nK','name' => 'nK','default' => '$nK'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '1','type' => 'integer','max' => '$MaxBlockALL','name' => 'MinBlockALL'}},{'type' => 't','content' => '

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
'}],'attrib' => {'name' => 'CHECKGRIDSIZE','if' => '$_IsFirstSession'}},{'type' => 'e','name' => 'command','content' => [{'type' => 't','content' => '
#BLOCKLEVELSRELOADED

This command means that the restart file contains the information about
the minimum and maximum allowed refinement levels for each block.
This command is only used in the restart header file.
'}],'attrib' => {'name' => 'BLOCKLEVELSRELOADED'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseTiming','default' => 'T'}},{'type' => 'e','name' => 'if','content' => [{'type' => 'e','name' => 'parameter','content' => [{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'none','value' => '-3'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'final only','value' => '-2','default' => 'T'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'end of sessions','value' => '-1'}},{'type' => 'e','name' => 'optioninput','content' => [],'attrib' => {'min' => '1','name' => 'every X steps','default' => '100'}}],'attrib' => {'type' => 'integer','name' => 'Frequency','input' => 'select'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '-1','type' => 'integer','name' => 'nDepthTiming','default' => '-1'}},{'type' => 'e','name' => 'parameter','content' => [{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'cummulative','value' => 'cumm','default' => '1'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'list'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'tree'}}],'attrib' => {'type' => 'string','name' => 'TypeTimingReport','input' => 'select'}}],'attrib' => {'expr' => '$UseTiming'}},{'type' => 't','content' => '
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
'}],'attrib' => {'name' => 'TIMING'}}],'attrib' => {'name' => 'TESTING AND TIMING'}},{'type' => 'e','name' => 'commandgroup','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! MAIN INITIAL AND BOUNDARY CONDITION PARAMETERS  !!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Uniform','value' => '1'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Shock tube','value' => '2'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Heliosphere','value' => '3'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Comet','value' => '5'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Rotation','value' => '6'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Diffusion','value' => '7'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Earth','value' => '11','default' => 'T'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Saturn','value' => '12'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Jupiter','value' => '13'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Venus','value' => '14'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Cylinder','value' => '21'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Sphere','value' => '22'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Arcade','value' => '25'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'CME','value' => '26'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Dissipation','value' => '30'}}],'attrib' => {'type' => 'integer','name' => 'iProblem','input' => 'select'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'string','name' => 'TypeDissipation','if' => '$iProblem==30','length' => '20'}},{'type' => 't','content' => '
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
'}],'attrib' => {'name' => 'PROBLEMTYPE','if' => '$_IsFirstSession','required' => 'T'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'GeoSolarMagnetic, GSM','if' => '$_NameComp eq \'GM\'','value' => 'GSM','default' => 'T'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'HelioGraphicInertial, HGI','if' => '$_NameComp ne \'GM\'','value' => 'HGI','default' => 'T'}}],'attrib' => {'type' => 'string','name' => 'TypeCoordSystem','input' => 'select'}},{'type' => 't','content' => '

#COORDSYSTEM
GSM			TypeCoordSystem

! TypeCoordSystem defines the coordinate system for the component.
! Currently only one coordinate system is available for GM ("GSM")
! and one for IH or SC ("HGI"). In the near future "GSE" should be also
! an option for GM.
!
! Default is component dependent: "GSM" for GM and "HGI" for IH or SC.
'}],'attrib' => {'name' => 'COORDSYSTEM','if' => '$_IsFirstSession'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'string','name' => 'NameRestartInDir','length' => '100','default' => 'GM/restartIN'}},{'type' => 't','content' => '

#RESTARTINDIR
GM/restart_n5000	NameRestartInDir

! The NameRestartInDir variable contains the name of the directory
! where restart files are saved relative to the run directory.
! The directory should be inside the subdirectory with the name 
! of the component.
!
! Default value is "GM/restartIN".
'}],'attrib' => {'name' => 'RESTARTINDIR','if' => '$_IsFirstSession'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'DoRestartBFace','default' => 'F'}},{'type' => 't','content' => ' 

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

'}],'attrib' => {'name' => 'NEWRESTART','if' => '$_IsFirstSession'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '1','type' => 'integer','name' => 'nRootBlockX','default' => '2'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '1','type' => 'integer','name' => 'nRootBlockY','default' => '1'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '1','type' => 'integer','name' => 'nRootBlockZ','default' => '1'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'xMin','default' => '-192.0'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '$xMin','type' => 'real','name' => 'xMax','default' => '  64.0'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'yMin','default' => ' -64.0'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '$yMin','type' => 'real','name' => 'yMax','default' => '  64.0'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'zMin','default' => ' -64.0'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '$zMin','type' => 'real','name' => 'zMax','default' => '  64.0'}},{'type' => 't','content' => '
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
'}],'attrib' => {'name' => 'GRID','if' => '$_IsFirstSession','required' => 'T'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'foreach','content' => [{'type' => 'e','name' => 'parameter','content' => [{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'coupled'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'fixed/inflow','default' => '$Side ne \'TypeBcEast\''}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'float/outflow','default' => '$Side eq \'TypeBcEast\''}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'heliofloat'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'reflect'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'periodic'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'vary'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'shear'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'linetied'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'raeder'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'arcadetop'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'arcadebot'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'arcadebotcont'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'user'}}],'attrib' => {'type' => 'string','name' => '$Side','input' => 'select'}}],'attrib' => {'name' => 'Side','values' => 'TypeBcEast,TypeBcWest,TypeBcSouth,TypeBcNorth,TypeBcBot,TypeBcTop'}},{'type' => 'e','name' => 'rule','content' => [{'type' => 't','content' => '
	! East and west BCs must be both periodic or neither
	'}],'attrib' => {'expr' => 'not($TypeBcEast eq \'periodic\' xor $TypeBcWest eq \'periodic\')'}},{'type' => 'e','name' => 'rule','content' => [{'type' => 't','content' => '
	! South and North BCs must be both periodic or neither
	'}],'attrib' => {'expr' => 'not($TypeBcSouth eq \'periodic\' xor $TypeBcNorth eq \'periodic\')'}},{'type' => 'e','name' => 'rule','content' => [{'type' => 't','content' => '
	! Bottom and top BCs must be both periodic or neither
	'}],'attrib' => {'expr' => 'not($TypeBcBot eq \'periodic\' xor $TypeBcTop eq \'periodic\')'}},{'type' => 't','content' => '
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
'}],'attrib' => {'name' => 'OUTERBOUNDARY'}},{'type' => 'e','name' => 'command','content' => [{'type' => 't','content' => '
! Inner boundary types for body 1 and body 2
	'},{'type' => 'e','name' => 'parameter','content' => [{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'reflect'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'float'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'fixed'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'ionosphere','default' => 'T'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'ionosphereB0','value' => 'ionosphereB0/ionosphereb0'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'ionospherefloat'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'coronatoih'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'user'}}],'attrib' => {'type' => 'string','name' => 'TypeBcInner','input' => 'select'}},{'type' => 'e','name' => 'if','content' => [{'type' => 'e','name' => 'parameter','content' => [{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'reflect','default' => 'T'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'float'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'fixed'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'ionosphere'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'ionosphereB0'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'ionospherefloat'}}],'attrib' => {'type' => 'string','name' => 'TypeBcBody2','input' => 'select'}}],'attrib' => {'expr' => '$UseBody2'}},{'type' => 'e','name' => 'rule','content' => [{'type' => 't','content' => '
! Note: for the second body COROTATION AND AN IONOSPHERIC BOUNDARY DO NOT WORK.
	'}],'attrib' => {'expr' => 'not($TypeBcBody2 =~ /ionosphere/)'}},{'type' => 't','content' => '

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
! \'coronatoih\'  - IH component obtains inner boundary from the SC component
! \'user\'        - user defined
!
! For \'ionosphere\' and \'ionospherefloat\' types and a coupled GM-IE run,
! the velocity at the inner boundary is determined by the ionosphere model.
!
! Default value for TypeBcInner is \'ionosphere\' for problem types
! Earth, Saturn, Jupiter, and rotation.
! For all other problems with an inner boundary the default is \'unknown\',
! so the inner boundary must be set.
'}],'attrib' => {'name' => 'INNERBOUNDARY'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseExtraBoundary','default' => 'F'}},{'type' => 'e','name' => 'if','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'string','name' => 'TypeBcExtra'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'DoFixExtraboundary','default' => 'F'}}],'attrib' => {'expr' => '$UseExtraBoundary'}}],'attrib' => {'name' => 'EXTRABOUNDARY'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'integer','max' => '6','name' => 'MaxBoundary','default' => '0'}},{'type' => 'e','name' => 'if','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'DoFixOuterBoundary','default' => 'F'}}],'attrib' => {'expr' => '$MaxBoundary >= 1'}},{'type' => 't','content' => '
#FACEOUTERBC
0              MaxBoundary            
F              DoFixOuterBoundary)    !read only for MaxBoundary>=East_(=1).
! If MaxBoundary is East_(=1) or more then the outer boundaries with
! the number of boundary being between East_ and MaxBoundary
! are treated using set_BCs.f90 subroutines instead of set_outerBCs.f90 
! if DoFixOuterBoundary is .true., there is no resolution
! change along the outer boundaries with the number of
! of boundary being between East_ and MaxBoundary
'}],'attrib' => {'name' => 'FACEOUTERBC'}}],'attrib' => {'name' => 'INITIAL AND BOUNDARY CONDITIONS'}},{'type' => 'e','name' => 'commandgroup','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! INITIAL TIME AND STEP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'integer','name' => 'iYear','default' => '2000'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '1','type' => 'integer','max' => '12','name' => 'iMonth','default' => '3'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '1','type' => 'integer','max' => '31','name' => 'iDay','default' => '21'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'integer','max' => '23','name' => 'iHour','default' => '0'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'integer','max' => '59','name' => 'iMinute','default' => '0'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'integer','max' => '59','name' => 'iSecond','default' => '0'}},{'type' => 't','content' => '
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
'}],'attrib' => {'alias' => 'SETREALTIME','name' => 'STARTTIME','if' => '$_IsFirstSession'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'tSimulation','default' => '0.0'}},{'type' => 't','content' => '

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
'}],'attrib' => {'name' => 'TIMESIMULATION','if' => '$_IsFirstSession'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'integer','name' => 'nStep','default' => '0'}},{'type' => 't','content' => '

#NSTEP
100			nStep

! Set nStep for the component. Typically used in the restart.H header file.
! Generally it is not inserted in a PARAM.in file by the user.
!
! The default is nStep=0 as the starting time step with no restart.
'}],'attrib' => {'name' => 'NSTEP','if' => '$_IsFirstSession'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '-1','type' => 'integer','name' => 'nPrevious','default' => '-1'}},{'type' => 't','content' => '

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
'}],'attrib' => {'name' => 'NPREVIOUS','if' => '$_IsFirstSession'}}],'attrib' => {'name' => 'INITIAL TIME'}},{'type' => 'e','name' => 'commandgroup','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  TIME INTEGRATION PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [{'type' => 'e','name' => 'option','content' => [],'attrib' => {'value' => '1','default' => 'T'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'value' => '2'}}],'attrib' => {'type' => 'integer','name' => 'nStage','input' => 'select'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','max' => '1','name' => 'CflExpl','default' => '0.8'}},{'type' => 't','content' => '

#TIMESTEPPING
2                       nStage
0.80                    CflExpl

! Parameters for explicit time integration.
! Default is 1 stage and CflExpl=0.8
'}],'attrib' => {'name' => 'TIMESTEPPING'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseDtFixed','default' => 'F'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'DtFixedDim','if' => '$UseDtFixed','default' => '1.0'}},{'type' => 't','content' => '
#FIXEDTIMESTEP
T                       UseDtFixed
10.                     DtFixedDim [sec] (read if UseDtFixed is true)

! Default is UseDtFixed=.false. Effective only if DoTimeAccurate is true.
! If UseDtFixed is true, the time step is fixed to DtFixedDim.
!
! This is useful for debugging explicit schemes.

! The real application is, however, for implicit and partially
! implicit/local schemes.

'}],'attrib' => {'name' => 'FIXEDTIMESTEP'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UsePartLocal','default' => 'F'}},{'type' => 't','content' => '
#PARTLOCAL
T               UsePartLocal

! Default is UsePartLocal=.false. If UsePartLocal is true and the
! run is time accurate, then the blocks selected as "implicit"
! by the criteria defined in #STEPPINGCRITERIA are not used to
! calculate the time step, and all cells are advanced with the
! smaller of the stable and the global time steps.
'}],'attrib' => {'name' => 'PARTLOCAL'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UsePointImplicit','default' => 'F'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UsePartImplicit','default' => 'F'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseFullImplicit','default' => 'F'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'CflImpl','if' => '$UsePartImplicit or $UseFullImplicit','default' => '100'}},{'type' => 'e','name' => 'rule','content' => [{'type' => 't','content' => '
	At most one of these logicals can be true!
	'}],'attrib' => {'expr' => '$UsePointImplicit + $UsePartImplicit + $UseFullImplicit <= 1'}},{'type' => 't','content' => '

#IMPLICIT
F               UsePointImplicit   
F               UsePartImplicit
F               UseFullImplicit
100.0           CflImpl (read if UsePartImplicit or UseFullImplicit is true)

! Default is false for all logicals. Only one of them can be set to true!
! The CFL number is used in the implicit blocks of the fully or partially
! implicit schemes. Ignored if UseDtFixed is true.
'}],'attrib' => {'name' => 'IMPLICIT'}}],'attrib' => {'name' => 'TIME INTEGRATION'}},{'type' => 'e','name' => 'commandgroup','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! PARAMETERS FOR FULL AND PART IMPLICIT TIME INTEGRATION !!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Time step','value' => 'dt','default' => 'T'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Radial distance','value' => 'r/R'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Test block','value' => 'test'}}],'attrib' => {'type' => 'string','name' => 'TypeImplCrit','input' => 'select'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'rImplicit','if' => '$TypeImplCrit eq \'R\''}},{'type' => 't','content' => '

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
'}],'attrib' => {'alias' => 'STEPPINGCRITERIA','name' => 'IMPLICITCRITERIA'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','max' => '1','name' => 'ImplCoeff','default' => '1'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseBdf2','default' => 'T'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseSourceImpl','default' => 'F'}},{'type' => 't','content' => '
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
'}],'attrib' => {'name' => 'IMPLSTEP'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => '1','default' => 'T'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => '2'}}],'attrib' => {'type' => 'integer','name' => 'nOrderImpl','input' => 'select'}},{'type' => 'e','name' => 'parameter','content' => [{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Roe','value' => 'Roe/roe/1'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Rusanov','value' => 'Rusanov/rusanov/2/TVDLF','default' => 'T'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Linde','value' => 'Linde/linde/3/HLLEL'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Sokolov','value' => 'Sokolov/sokolov/4/AW'}}],'attrib' => {'type' => 'string','name' => 'TypeFluxImpl','input' => 'select'}},{'type' => 't','content' => '
#IMPLSCHEME
1               nOrderImpl
Rusanov         TypeFluxImpl

! Default values are shown, ie. first order Rusanov scheme.
! This defines the scheme used in the implicit part.
'}],'attrib' => {'name' => 'IMPLSCHEME'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseConservativeImplicit','default' => 'F'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseNewton','default' => 'F'}},{'type' => 'e','name' => 'if','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseNewMatrix','default' => 'T'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '1','type' => 'integer','name' => 'MaxIterNewton','default' => '10'}}],'attrib' => {'expr' => '$UseNewton'}},{'type' => 't','content' => '
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
'}],'attrib' => {'name' => 'NEWTON'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Preconditioned','value' => 'prec','default' => 'T'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'No preconditioning','value' => 'free'}}],'attrib' => {'type' => 'string','name' => 'TypeJacobian','input' => 'select'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','max' => '1.e-5','name' => 'JacobianEps','default' => '$doublePrecision ? 1.e-12 : 1.e-6'}},{'type' => 't','content' => '
#JACOBIAN
prec            TypeJacobian (prec, free)
1.E-12          JacobianEps

! The Jacobian matrix is always calculated with a matrix free approach,
! however it can be preconditioned  (\'prec\'), or not (\'free\')
! Default value is TypeJacobian=\'prec\'.
! JacobianEps contains the machine round off error for numerical derivatives.
! The default value is 1.E-12 for 8 byte reals and 1.E-6 for 4 byte reals.
'}],'attrib' => {'name' => 'JACOBIAN'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'left'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'symmetric','default' => 'T'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'right'}}],'attrib' => {'type' => 'string','name' => 'TypePrecondSide','input' => 'select'}},{'type' => 'e','name' => 'parameter','content' => [{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'MBILU','default' => 'T'}}],'attrib' => {'type' => 'string','name' => 'TypePrecond','input' => 'select'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','max' => '1','name' => 'GustafssonPar','default' => '0.5'}},{'type' => 't','content' => '
#PRECONDITIONER
symmetric       TypePrecondSide (left, symmetric, right)
MBILU           TypePrecond (MBILU)
0.5             GustafssonPar (0. no modification, 1. full modification)

! Default parameters are shown. Right preconditioning does not affect
! the normalization of the residual. The Gustafsson parameter determines
! how much the MBILU preconditioner is modified. The default 0.5 value
! means a relaxed modification.
'}],'attrib' => {'name' => 'PRECONDITIONER'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'gmres','default' => 'T'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'bicgstab'}}],'attrib' => {'type' => 'string','name' => 'TypeKrylov','input' => 'select'}},{'type' => 'e','name' => 'parameter','content' => [{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => '0','value' => 'nul','default' => 'T'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'previous','value' => 'old'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'explicit'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'scaled explicit','value' => 'explicit'}}],'attrib' => {'type' => 'string','name' => 'TypeInitKrylov','input' => 'select'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','max' => '0.1','name' => 'ErrorMaxKrylov','default' => '0.001'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '1','type' => 'integer','name' => 'MaxMatvecKrylov','default' => '100'}},{'type' => 't','content' => '
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
'}],'attrib' => {'name' => 'KRYLOV'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '1','type' => 'integer','name' => 'nKrylovVector','default' => 'MaxMatvecKrylov'}},{'type' => 't','content' => '
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
'}],'attrib' => {'name' => 'KRYLOVSIZE'}}],'attrib' => {'name' => 'IMPLICIT PARAMETERS'}},{'type' => 'e','name' => 'commandgroup','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! STOPPING CRITERIA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

The commands in this group only work in stand alone mode.

'},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '-1','type' => 'integer','name' => 'MaxIteration','default' => '-1'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '-1','type' => 'real','name' => 'tSimulationMax','default' => '-1'}},{'type' => 't','content' => '

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
'}],'attrib' => {'name' => 'STOP','if' => '$_IsStandAlone','required' => '$_IsStandAlone'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'DoCheckStopFile','default' => 'T'}},{'type' => 't','content' => '

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
'}],'attrib' => {'name' => 'CHECKSTOPFILE','if' => '$_IsStandAlone'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '-1','type' => 'real','name' => 'CpuTimeMax','default' => '-1'}},{'type' => 't','content' => '

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
'}],'attrib' => {'name' => 'CPUTIMEMAX','if' => '$_IsStandAlone'}}],'attrib' => {'name' => 'STOPPING CRITERIA'}},{'type' => 'e','name' => 'commandgroup','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  OUTPUT PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'string','name' => 'NameRestartOutDir','length' => '100','default' => 'GM/restartOUT'}},{'type' => 't','content' => '

#RESTARTOUTDIR
GM/restart_n5000	NameRestartOutDir

! The NameRestartOutDir variable contains the name of the directory
! where restart files are saved relative to the run directory.
! The directory should be inside the subdirectory with the name 
! of the component.
!
! Default value is "GM/restartOUT".
'}],'attrib' => {'name' => 'RESTARTOUTDIR'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'DoSaveRestart','default' => 'T'}},{'type' => 'e','name' => 'if','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '-1','type' => 'integer','name' => 'DnSaveRestart','default' => '-1'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '-1','type' => 'real','name' => 'DtSaveRestart','default' => '-1'}}],'attrib' => {'expr' => '$DoSaveRestart'}},{'type' => 't','content' => '
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
'}],'attrib' => {'name' => 'SAVERESTART'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'string','name' => 'NamePlotDir','length' => '100','default' => 'GM/IO2'}},{'type' => 't','content' => '

The NamePlotDir variable contains the name of the directory
where plot files and logfiles are saved relative to the run directory.
The directory should be inside the subdirectory with the name
of the component.

Default value is "GM/IO2".
'}],'attrib' => {'name' => 'PLOTDIR'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'DoSaveLogfile','default' => 'F'}},{'type' => 'e','name' => 'if','content' => [{'type' => 'e','name' => 'parameter','content' => [{'type' => 'e','name' => 'part','content' => [{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'MHD vars. dimensional','value' => 'MHD','default' => 'T'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Flux vars. dimensional','value' => 'FLX'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Raw vars. dimensional','value' => 'RAW'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Set vars. dimensional','value' => 'VAR'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'MHD vars. scaled','value' => 'mhd','default' => 'T'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Flux vars. scaled','value' => 'flx'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Raw vars. scaled','value' => 'raw'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Set vars. scaled','value' => 'var'}}],'attrib' => {'type' => 'string','name' => 'TypeLogVar','required' => 'T','input' => 'select'}},{'type' => 'e','name' => 'part','content' => [{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'none','exclusive' => 'T'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'step'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'date'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'time'}}],'attrib' => {'type' => 'string','multiple' => 'T','name' => 'TypeTime','required' => 'F','input' => 'select'}}],'attrib' => {'min' => '1','type' => 'strings','max' => '4','name' => 'StringLog'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '-1','type' => 'integer','name' => 'DnSaveLogfile','default' => '1'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '-1','type' => 'real','name' => 'DtSaveLogfile','default' => '-1'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'string','name' => 'NameLogVars','if' => '$TypeLogVar =~ /var/i','length' => '100'}},{'type' => 'e','name' => 'parameter','content' => [{'type' => 'e','name' => 'part','content' => [],'attrib' => {'min' => '$rBody','type' => 'real','multiple' => 'T','name' => 'LogRadii'}}],'attrib' => {'min' => '1','type' => 'strings','max' => '10','name' => 'StringLogRadii','if' => '($TypeLogVar=~/flx/i or $NameLogVars=~/flx/i)','length' => '100'}}],'attrib' => {'expr' => '$DoSaveLogfile'}},{'type' => 't','content' => '
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
'}],'attrib' => {'name' => 'SAVELOGFILE'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'integer','name' => 'nSatellite','default' => '0'}},{'type' => 'e','name' => 'for','content' => [{'type' => 'e','name' => 'parameter','content' => [{'type' => 'e','name' => 'part','content' => [{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'MHD vars. dimensional','value' => 'MHD','default' => 'T'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'All vars. dimensional','value' => 'FUL'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Set vars. dimensional','value' => 'VAR'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'MHD vars. scaled','value' => 'mhd'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'All vars. scaled','value' => 'ful'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Set vars. scaled','value' => 'var'}}],'attrib' => {'type' => 'string','name' => 'TypeSatelliteVar','required' => 'T','input' => 'select'}},{'type' => 'e','name' => 'part','content' => [{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'file','default' => 'T'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'equation','value' => 'eqn'}}],'attrib' => {'type' => 'string','name' => 'TypeTrajectory','required' => 'F','input' => 'select'}},{'type' => 'e','name' => 'part','content' => [{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'none','exclusive' => 'T'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'step'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'date'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'time'}}],'attrib' => {'type' => 'string','multiple' => 'T','name' => 'TypeTime','required' => 'F','input' => 'select'}}],'attrib' => {'min' => '1','type' => 'strings','max' => '5','name' => 'StringSatellite'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '-1','type' => 'integer','name' => 'DnOutput','default' => '1'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '-1','type' => 'real','name' => 'DtOutput','default' => '-1'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'string','name' => 'NameTrajectoryFile','length' => '100'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'string','name' => 'NameSatelliteVars','if' => '$TypeSatelliteVar =~ /\\bvar\\b/i','length' => '100'}}],'attrib' => {'from' => '1','to' => '$nSatellite'}},{'type' => 't','content' => '
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
'}],'attrib' => {'name' => 'SATELLITE','if' => '$_IsFirstSession'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'integer','max' => '100','name' => 'nPlotFile','default' => '0'}},{'type' => 'e','name' => 'for','content' => [{'type' => 'e','name' => 'parameter','content' => [{'type' => 'e','name' => 'part','content' => [{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'TECPLOT','value' => 'tec'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'IDL','value' => 'idl'}}],'attrib' => {'type' => 'string','name' => 'plotform','required' => 'T','input' => 'select'}},{'type' => 'e','name' => 'part','content' => [{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => '3D','value' => '3d/3d_'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'value' => 'x=0'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'value' => 'y=0','default' => 'T'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'value' => 'z=0'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'value' => 'sph'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'value' => 'los'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'value' => 'lin'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'if' => '$plotform =~ /\\bidl\\b/','value' => 'cut'}}],'attrib' => {'type' => 'string','name' => 'plotarea','required' => 'T','input' => 'select'}},{'type' => 'e','name' => 'part','content' => [{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'MHD vars. dimensional','value' => 'MHD'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'All vars. dimensional','value' => 'FUL'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Raw vars. dimensional','value' => 'RAW'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Ray tracing vars. dim.','value' => 'RAY'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Flux vars. dimensional','value' => 'FLX'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Solar vars. dimensional','value' => 'SOL'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Position vars. dimensional','if' => '$plotarea eq \'lin\'','value' => 'POS'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Select dimensional vars.','value' => 'VAR'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'MHD vars. scaled','value' => 'mhd'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'All vars. scaled','value' => 'ful'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Raw vars. scaled','value' => 'raw'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Ray tracing vars. scaled','value' => 'ray'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Flux vars. scaled','value' => 'flx'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Solar vars. scaled','value' => 'sol'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Position vars. scaled','if' => '$plotarea eq \'lin\'','value' => 'pos'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Select scaled vars.','value' => 'var'}}],'attrib' => {'type' => 'string','name' => 'plotvar','required' => 'T','input' => 'select'}}],'attrib' => {'min' => '3','type' => 'strings','max' => '3','name' => 'StringPlot'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '-1','type' => 'integer','name' => 'DnSavePlot'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '-1','type' => 'real','name' => 'DtSavePlot'}},{'type' => 'e','name' => 'if','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'xMinCut'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '$xMinCut','type' => 'real','name' => 'xMaxCut'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'yMinCut'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '$yMinCut','type' => 'real','name' => 'yMaxCut'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'zMinCut'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '$zMinCut','type' => 'real','name' => 'zMaxCut'}}],'attrib' => {'expr' => '$plotarea =~ /\\bcut\\b/'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'Radius','if' => '$plotarea =~ /\\bsph\\b/','default' => '10'}},{'type' => 'e','name' => 'if','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'LosVectorX','default' => '0'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'LosVectorY','default' => '0.0001'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'LosVectorZ','default' => '1'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'xSizeImage','default' => '20'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'ySizeImage','default' => '20'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'xOffset','default' => '10'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'yOffset','default' => '10'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '1','type' => 'real','name' => 'rOccult','default' => '2.5'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','max' => '1','name' => 'MuLimbDarkening','default' => '0.5'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '2','type' => 'integer','name' => 'nPixX','default' => '200'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '2','type' => 'integer','name' => 'nPixY','default' => '200'}}],'attrib' => {'expr' => '$plotarea =~ /\\blos\\b/'}},{'type' => 'e','name' => 'if','content' => [{'type' => 'e','name' => 'parameter','content' => [{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Advected B','value' => 'A'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'B','default' => 'T'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'U'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'J'}}],'attrib' => {'type' => 'string','name' => 'NameLine','input' => 'select'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'IsSingleLine','default' => 'F'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '1','type' => 'integer','max' => '20','name' => 'nLine','default' => '1'}},{'type' => 'e','name' => 'for','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'xStartLine'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'yStartLine'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'zStartLine'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'IsParallel'}}],'attrib' => {'from' => '1','to' => '$nLine'}}],'attrib' => {'expr' => '$plotarea =~ /\\blin\\b/'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '-1.0','type' => 'real','name' => 'DxSavePlot','if' => '($plotform=~/\\bidl\\b/ and $plotarea!~/\\b(sph|los|lin)\\b/)','default' => '-1.0'}},{'type' => 'e','name' => 'if','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'string','name' => 'NameVars','length' => '100'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'string','name' => 'NamePars','length' => '100'}}],'attrib' => {'expr' => '$plotvar =~ /\\bvar\\b/i'}}],'attrib' => {'from' => '1','to' => '$nPlotFile'}},{'type' => 't','content' => '
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

'}],'attrib' => {'name' => 'SAVEPLOT'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'DoSaveBinary','default' => 'T'}},{'type' => 't','content' => '
#SAVEBINARY
T			DoSaveBinary   used only for \'idl\' plot file

! Default is .true. Saves unformatted IO2/*.idl files if true. 
! This is the recommended method, because it is fast and accurate.
! The only advantage of saving IO2/*.idl in formatted text files is
! that it can be processed on another machine or with a different 
! (lower) precision. For example PostIDL.exe may be compiled with 
! single precision to make IO2/*.out files smaller, while BATSRUS.exe is 
! compiled in double precision, to make results more accurate.
'}],'attrib' => {'name' => 'SAVEBINARY'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'DoSavePlotsAmr','default' => 'F'}},{'type' => 't','content' => '
#SAVEPLOTSAMR
F			DoSavePlotsAmr

! Save plots before each AMR. Default is DoSavePlotsAMR=.false.
'}],'attrib' => {'name' => 'SAVEPLOTSAMR'}}],'attrib' => {'name' => 'OUTPUT PARAMETERS'}},{'type' => 'e','name' => 'commandgroup','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  AMR PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'default','default' => '1'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'all'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'none'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => '3Dbodyfocus'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'spherefocus'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'magnetosphere'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'points'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'coupledhelio'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'helio_init'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'helio_z=4'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'all_then_focus'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'cme'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'points'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'mag_new'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'magnetosphere'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'magneto_fine'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'magneto12'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'magnetosaturn'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'magnetojupiter'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'paleo'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'comet'}}],'attrib' => {'type' => 'string','name' => 'InitialRefineType','input' => 'select'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'integer','name' => 'InitialRefineLevel','default' => '4'}},{'type' => 't','content' => '
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
'}],'attrib' => {'name' => 'AMRINIT','if' => '$_IsFirstSession'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'integer','name' => 'nRefineLevelIC','default' => '0'}},{'type' => 't','content' => '
#AMRINITPHYSICS
3			nRefineLevelIC

! Defines number of physics (initial condition) based AMR-s AFTER the 
! geometry based initial AMR-s defined by #AMRINIT were done.
! Only useful if the initial condition has a non-trivial analytic form.
'}],'attrib' => {'name' => 'AMRINITPHYSICS','if' => '$_IsFirstSession'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '-1','type' => 'integer','name' => 'MinBlockLevel','default' => '0'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '-1','type' => 'integer','name' => 'MaxBlockLevel','default' => '99'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'DoFixBodyLevel','default' => 'F'}},{'type' => 't','content' => '
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
'}],'attrib' => {'name' => 'AMRLEVELS'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '-1','type' => 'real','name' => 'DxCellMin','default' => '0'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '-1','type' => 'real','name' => 'DxCellMax','default' => '99999'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'DoFixBodyLevel','default' => 'F'}},{'type' => 't','content' => '
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
'}],'attrib' => {'name' => 'AMRRESOLUTION'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '-1','type' => 'integer','name' => 'DnRefine','default' => '-1'}},{'type' => 'e','name' => 'if','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'DoAutoRefine','default' => 'F'}},{'type' => 'e','name' => 'if','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','max' => '100','name' => 'PercentCoarsen','default' => '20'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','max' => '100','name' => 'PercentRefine','default' => '20'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '1','type' => 'integer','name' => 'MaxTotalBlocks','default' => '99999'}}],'attrib' => {'expr' => '$DoAutoRefine'}}],'attrib' => {'expr' => '$DnRefine>0'}},{'type' => 't','content' => '
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
'}],'attrib' => {'name' => 'AMR'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => '1'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => '2'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => '3','default' => '1'}}],'attrib' => {'type' => 'integer','name' => 'nRefineCrit','input' => 'select'}},{'type' => 'e','name' => 'for','content' => [{'type' => 'e','name' => 'parameter','content' => [{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'grad T','value' => 'gradt/gradT'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'grad P','value' => 'gradp/gradP'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'grad log(Rho)','value' => 'gradlogrho'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'grad log(p)','value' => 'gradlogP/gradlogp'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'grad E','value' => 'gradE'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'curl U','value' => 'curlV/curlv/curlU/curlu'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'curl B','value' => 'curlB/curlb'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'div U','value' => 'divU/divu/divV/divv'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'divB','value' => 'divb/divB'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'vAlfven','value' => 'Valfven/vAlfven/valfven'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'heliospheric beta','value' => 'heliobeta'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'flux'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'heliospheric current sheet','value' => 'heliocurrentsheet'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'rCurrents','value' => 'rcurrents/Rcurrents'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Transient','value' => 'transient/Transient'}}],'attrib' => {'type' => 'string','name' => 'TypeRefine','input' => 'select'}},{'type' => 'e','name' => 'if','content' => [{'type' => 'e','name' => 'parameter','content' => [{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'P_dot','value' => 'p_dot/P_dot'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'T_dot','value' => 't_dot/T_dot'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Rho_dot','value' => 'rho_dot/Rho_dot','default' => 'T'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'RhoU_dot','value' => 'RhoU_dot/rhou_dot'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Rho_2nd_1','value' => 'Rho_2nd_1/rho_2nd_1'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Rho_2nd_2','value' => 'Rho_2nd_2/rho_2nd_2'}}],'attrib' => {'type' => 'string','name' => 'TypeTransient','input' => 'select'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseSunEarth','default' => 'F'}},{'type' => 'e','name' => 'if','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'xEarth'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'yEarth'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'zEarth'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'InvD2Ray'}}],'attrib' => {'expr' => '$UseSunEarth'}}],'attrib' => {'expr' => '$TypeRefine =~ /transient/i'}}],'attrib' => {'from' => '1','to' => '$nRefineCrit'}},{'type' => 't','content' => '
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
'}],'attrib' => {'name' => 'AMRCRITERIA'}}],'attrib' => {'name' => 'AMR PARAMETERS'}},{'type' => 'e','name' => 'commandgroup','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  SCHEME PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => '1','default' => 'T'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => '2'}}],'attrib' => {'type' => 'integer','name' => 'nOrder','input' => 'select'}},{'type' => 'e','name' => 'parameter','content' => [{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Roe','value' => 'Roe/roe/1'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Rusanov','value' => 'Rusanov/rusanov/2/TVDLF','default' => 'T'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Linde','value' => 'Linde/linde/3/HLLEL'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Sokolov','value' => 'Sokolov/sokolov/4/AW'}}],'attrib' => {'type' => 'string','name' => 'TypeFlux','input' => 'select'}},{'type' => 'e','name' => 'if','content' => [{'type' => 'e','name' => 'parameter','content' => [{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'minmod','default' => 'T'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'beta'}}],'attrib' => {'type' => 'string','name' => 'TypeLimiter','input' => 'select'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '1','type' => 'real','max' => '2','name' => 'LimiterBeta','if' => '$TypeLimiter eq \'beta\'','default' => '1.2'}}],'attrib' => {'expr' => '$nOrder == 2'}},{'type' => 't','content' => '
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
'}],'attrib' => {'name' => 'SCHEME'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseNonConservative','default' => 'T'}},{'type' => 't','content' => '
#NONCONSERVATIVE
T		UseNonConservative

! For Earth the default is using non-conservative equations 
! (close to the body).
'}],'attrib' => {'name' => 'NONCONSERVATIVE'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'integer','max' => '3','name' => 'nConservCrit','default' => '1'}},{'type' => 'e','name' => 'for','content' => [{'type' => 'e','name' => 'parameter','content' => [{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'radius','value' => 'r/R/radius/Radius','default' => 'T'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'p','value' => 'p/P'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'grad P','value' => 'gradp/GradP'}}],'attrib' => {'type' => 'string','name' => 'TypeConservCrit','input' => 'select'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '$rBody','type' => 'real','name' => 'rConserv','if' => '$TypeConservCrit =~ /^r|radius$/i','default' => '2*$rBody'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'pCoeffConserv','if' => '$TypeConservCrit =~ /^p$/i','default' => '0.05'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'GradPCoeffConserv','if' => '$TypeConservCrit =~ /gradp/i','default' => '0.1'}}],'attrib' => {'from' => '1','to' => '$nConservCrit'}},{'type' => 't','content' => '
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
'}],'attrib' => {'name' => 'CONSERVATIVECRITERIA'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseUpdateCheck','default' => 'T'}},{'type' => 'e','name' => 'if','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','max' => '100','name' => 'RhoMinPercent','default' => '40'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '100','type' => 'real','name' => 'RhoMaxPercent','default' => '400'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','max' => '100','name' => 'pMinPercent','default' => '40'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '100','type' => 'real','name' => 'pMaxPercent','default' => '400'}}],'attrib' => {'expr' => '$UseUpdateCheck'}},{'type' => 't','content' => '
#UPDATECHECK
T			UseUpdateCheck
40.			RhoMinPercent
400.			RhoMaxPercent
40.			pMinPercent
400.			pMaxPercent

! Default values are shown.  This will adjust the timestep so that
! density and pressure cannot change by more than the given percentages
! in a single timestep.
'}],'attrib' => {'name' => 'UPDATECHECK'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => '1','default' => 'T'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => '2'}}],'attrib' => {'type' => 'integer','name' => 'nOrderProlong','input' => 'select'}},{'type' => 'e','name' => 'parameter','content' => [{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'left-right','value' => 'lr','default' => 'T'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'central'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'minmod'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'left-right extrapolate','value' => 'lr2'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'central    extrapolate','value' => 'central2'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'minmod     extrapolate','value' => 'minmod2'}}],'attrib' => {'type' => 'string','name' => 'TypeProlong','if' => '$nOrderProlong==2','input' => 'select'}},{'type' => 't','content' => '
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
'}],'attrib' => {'name' => 'PROLONGATION'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'm_p_cell FACES ONLY','value' => 'allopt','default' => 'T'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'm_p_cell','value' => 'all'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'm_p_dir FACES ONLY','value' => 'opt'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'm_p_dir group by directions','value' => 'dir'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'm_p_dir group by faces     ','value' => 'face'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'm_p_dir group by kind and face','value' => 'min'}}],'attrib' => {'type' => 'string','name' => 'TypeMessagePass','input' => 'select'}},{'type' => 't','content' => '
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
'}],'attrib' => {'alias' => 'OPTIMIZE','name' => 'MESSAGEPASS'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseBorisCorrection','default' => 'F'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','max' => '1','name' => 'BorisClightFactor','if' => '$UseBorisCorrection','default' => '1'}},{'type' => 't','content' => '
#BORIS
T			UseBorisCorrection
1.0			BorisClightFactor !Only if UseBorisCorrection is true

! Default is boris_correction=.false.
! Use semi-relativistic MHD equations with speed of light reduced by
! the BorisClightFactor. Set BorisClightFactor=1.0 for true semi-relativistic
! MHD. Gives the same steady state as normal MHD analytically, but there
! can be differences due to discretization errors. 
! You can use either Boris or BorisSimple but not both. 
'}],'attrib' => {'name' => 'BORIS'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseBorisSimple','default' => 'F'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','max' => '1','name' => 'BorisClightFactor','if' => '$UseBorisSimple','default' => '1'}},{'type' => 't','content' => '
#BORISSIMPLE
T			UseBorisSimple
0.05			BorisClightFactor !Only if UseBorisSimple is true

! Default is UseBorisSimple=.false. 
! Use simplified semi-relativistic MHD with speed of light reduced by the
! BorisClightFactor. This is only useful with BorisClightFactor less than 1.
! Should give the same steady state as normal MHD, but there can be a
! difference due to discretization errors.
! You can use either Boris or BorisSimple but not both. 
'}],'attrib' => {'name' => 'BORISSIMPLE'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseDivbSource','default' => 'T'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseDivbDiffusion','default' => 'F'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseProjection','default' => 'F'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseConstrainB','default' => 'F'}},{'type' => 'e','name' => 'rule','content' => [{'type' => 't','content' => '
	! At least one of the options should be true.
	'}],'attrib' => {'expr' => '$UseDivbSource or $UseDivbDiffusion or $UseProjection or $UseConstrainB'}},{'type' => 'e','name' => 'rule','content' => [{'type' => 't','content' => '
	! If UseProjection is true, all others should be false.
	'}],'attrib' => {'expr' => 'not($UseProjection and ($UseDivbSource or $UseDivbDiffusion or $UseConstrainB))'}},{'type' => 'e','name' => 'rule','content' => [{'type' => 't','content' => '
	! If UseConstrainB is true, all others should be false.
	'}],'attrib' => {'expr' => 'not($UseConstrainB and ($UseDivbSource or $UseDivbDiffusion or $UseProjection))'}},{'type' => 't','content' => '
	
#DIVB
T			UseDivbSource
F			UseDivbDiffusion	
F			UseProjection           
F			UseConstrainB           

! Default values are shown above.
! If UseProjection is true, all others should be false.
! If UseConstrainB is true, all others should be false.
! At least one of the options should be true.
'}],'attrib' => {'name' => 'DIVB'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseB0Source','default' => 'T'}},{'type' => 't','content' => '
#DIVBSOURCE
T			UseB0Source

! Add extra source terms related to the non-zero divergence and curl of B0.
! Default is true.
'}],'attrib' => {'name' => 'DIVBSOURCE'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','max' => '1','name' => 'DivbDiffCoeff','default' => '0.1666667'}},{'type' => 't','content' => '
#DIVBDIFFUSION
0.1666667		DivbDiffCoeff

! Default value is shown above. 1.0/6.0
! If divb_diffcoeff .gt. 0.5 then cfl .lt. 0.5/DivbDiffCoeff is required!
'}],'attrib' => {'name' => 'DIVBDIFFUSION'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Conjugate Gradients','value' => 'cg','default' => 'T'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'BiCGSTAB','value' => 'bicgstab'}}],'attrib' => {'type' => 'string','name' => 'TypeProjectIter','input' => 'select'}},{'type' => 'e','name' => 'parameter','content' => [{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Relative norm','value' => 'rel','default' => 'T'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Maximum error','value' => 'max'}}],'attrib' => {'type' => 'string','name' => 'TypeProjectStop','input' => 'select'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','max' => '1','name' => 'RelativeLimit','default' => '0.1'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'AbsoluteLimit','default' => '0.0'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '1','type' => 'integer','name' => 'MaxMatvec','default' => '50'}},{'type' => 't','content' => '
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

'}],'attrib' => {'name' => 'PROJECTION'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','max' => '1','name' => 'pRatioLow','default' => '0.01'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '$pRatioLow','type' => 'real','max' => '1','name' => 'pRatioHigh','default' => '0.1'}},{'type' => 't','content' => '
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
'}],'attrib' => {'name' => 'CORRECTP'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseAccurateIntegral','default' => 'T'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseAccurateTrace','default' => 'F'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0.01','type' => 'real','max' => '60','name' => 'DtExchangeRay','default' => '0.1'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '1','type' => 'integer','name' => 'DnRaytrace','default' => '1'}},{'type' => 't','content' => '

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
'}],'attrib' => {'name' => 'RAYTRACE'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'TauCoupleIm'}},{'type' => 't','content' => '

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
'}],'attrib' => {'name' => 'IM'}}],'attrib' => {'name' => 'SCHEME PARAMETERS'}},{'type' => 'e','name' => 'commandgroup','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  PHYSICS PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '1','type' => 'real','name' => 'Gamma','default' => '1.6666666667'}},{'type' => 't','content' => '
#GAMMA
1.6666666667		Gamma

! The adiabatic index (ratio of the specific heats for fixed pressure
! and volume. The default value is 5.0/3.0, which is valid for
! monoatomic gas or plasma.
'}],'attrib' => {'name' => 'GAMMA','if' => '$_IsFirstSession'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'RhoLeft','default' => '1'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'UnLeft','default' => '0'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'Ut1Left','default' => '0'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'Ut2Left','default' => '0'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'BnLeft','default' => '0.75'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'Bt1Left','default' => '1'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'Bt2Left','default' => '0'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'pRight','default' => '1'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'RhoRight','default' => '0.125'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'UnRight','default' => '0'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'Ut1Right','default' => '0'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'Ut2Right','default' => '0'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'BnRight','default' => '0.75'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'Bt1Right','default' => '-1'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'Bt2Right','default' => '0'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'pRight','default' => '0.1'}},{'type' => 'e','name' => 'parameter','content' => [{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'no rotation','value' => '0','default' => 'T'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'value' => '0.25'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'value' => '0.3333333333333'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'value' => '0.5'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'value' => '1'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'value' => '2'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'value' => '3'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'value' => '4'}}],'attrib' => {'type' => 'real','name' => 'ShockSlope','input' => 'select'}},{'type' => 't','content' => '
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
'}],'attrib' => {'name' => 'SHOCKTUBE'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '-1','type' => 'real','name' => 'SwRhoDim','default' => '5'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '-1','type' => 'real','name' => 'SwTDim','default' => '181712.175'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','max' => '0','name' => 'SwUxDim','default' => '-400'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'SwUyDim','default' => '0'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'SwUzDim','default' => '0'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'SwBxDim','default' => '0'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'SwByDim','default' => '0'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'SwBzDim','default' => '5'}},{'type' => 't','content' => '
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
'}],'attrib' => {'name' => 'SOLARWIND','if' => '$_IsFirstSession'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseUpstreamInputFile','default' => 'F'}},{'type' => 'e','name' => 'if','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'string','name' => 'NameUpstreamFile','length' => '100'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'SatelliteYPos','default' => '0'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'SatelliteZPos','default' => '0'}}],'attrib' => {'expr' => '$UseUpstreamInputFile'}},{'type' => 't','content' => '
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
'}],'attrib' => {'name' => 'UPSTREAM_INPUT_FILE'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseBody','default' => 'F'}},{'type' => 'e','name' => 'if','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'rBody','default' => '3'}},{'type' => 'e','name' => 'if','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '-1','type' => 'real','name' => 'rCurrents','default' => '4'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'BodyRhoDim','default' => '1'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'BodyTDim','default' => '10000'}}],'attrib' => {'expr' => '$_NameComp eq \'GM\''}}],'attrib' => {'expr' => '$UseBody'}},{'type' => 't','content' => '
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
'}],'attrib' => {'alias' => 'MAGNETOSPHERE','name' => 'BODY','if' => '$_IsFirstSession'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseGravity','default' => 'F'}},{'type' => 'e','name' => 'parameter','content' => [{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'central mass','value' => '0','default' => 'T'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'X direction','value' => '1'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Y direction','value' => '2'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Z direction','value' => '3'}}],'attrib' => {'type' => 'integer','name' => 'iDirGravity','if' => '$UseGravity','input' => 'select'}},{'type' => 't','content' => '
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

'}],'attrib' => {'name' => 'GRAVITY','if' => '$_IsFirstSession'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseMassLoading','default' => 'F'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'DoAccelerateMassLoading','default' => 'F'}},{'type' => 't','content' => '
#MASSLOADING
F			UseMassLoading
F			DoAccelerateMassLoading
'}],'attrib' => {'name' => 'MASSLOADING'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseHeatFlux','default' => 'F'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseSpitzerForm','default' => 'T'}},{'type' => 'e','name' => 'if','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'Kappa0Heat','default' => '1.23E-11'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'Kappa0Heat','default' => '2.5'}}],'attrib' => {'expr' => 'not $UseSpitzerForm'}},{'type' => 't','content' => '
#HEATFLUX
T		UseHeatFlux
F		UseSpitzerForm
1.23E-11	Kappa0Heat [W/m/K]	! Only if not UseSpitzerForm
2.50E+00	ExponentHeat [-]	! Only if not UseSpitzerForm
'}],'attrib' => {'name' => 'HEATFLUX'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseResistFlux','default' => 'F'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseSpitzerForm','default' => 'T'}},{'type' => 'e','name' => 'if','content' => [{'type' => 'e','name' => 'parameter','content' => [{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'localized','value' => 'Localized/localized'}},{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'constant','value' => 'Constant/constant','default' => 'T'}}],'attrib' => {'type' => 'string','name' => 'TypeResist','input' => 'select'}},{'type' => 'e','name' => 'if','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'Eta0Resist','default' => '9.69953E+8'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'Alpha0Resist','default' => '150'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'yShiftResist','default' => '0.5'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'TimeInitRise','default' => '0.05'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'TimeConstLev','default' => '1'}}],'attrib' => {'expr' => '$TypeResist =~ /localized/i'}}],'attrib' => {'expr' => 'not $UseSpitzerForm'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseAnomResist','default' => 'F'}},{'type' => 'e','name' => 'if','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'Eta0AnomResist','default' => '1.93991E+09'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'EtaAnomMaxResist','default' => '1.93991E+10'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'ThresholdFactorResist','default' => '1'}}],'attrib' => {'expr' => '$UseAnomResist'}},{'type' => 't','content' => '
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
'}],'attrib' => {'name' => 'RESISTIVEFLUX'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseDefaultUnits','default' => 'T'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'Grav0Diss','default' => '2.635620E-02'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'Beta0Diss','default' => '1.640000E-01'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'Length0Diss','default' => '1.500000E+06'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'Time0Diss','default' => '1.159850E+01'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'Rho0Diss','default' => '5.019000E-11'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'Tem0Diss','default' => '1.000000E+05'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'Theta0Diss','default' => '6.000000E-01'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'Delta0Diss','default' => '2.500000E+01'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'EpsilonDiss','default' => '7.000000E+00'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'RhoDifDiss','default' => '4.500000E+00'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'yShiftDiss','default' => '4.000000E-01'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'ScaleHeightDiss','default' => '5.000000E-01'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'ScaleFactorDiss','default' => '1.159850E+01'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'BZ0iss','default' => '5.000000E-01'}},{'type' => 't','content' => '
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
'}],'attrib' => {'name' => 'TESTDISSMHD','if' => '$_IsFirstSession'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseBody2','default' => 'F'}},{'type' => 'e','name' => 'if','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'rBody2','default' => '0.1'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '$xMin','type' => 'real','max' => '$xMax','name' => 'xBody2','default' => '-40'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '$yMin','type' => 'real','max' => '$yMax','name' => 'yBody2','default' => '0'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '$zMin','type' => 'real','max' => '$zMax','name' => 'zBody2','default' => '0'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '$rBody2','type' => 'real','name' => 'rCurrents2','default' => '1.3*$rBody2'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'RhoDimBody2','default' => '5'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'tDimBody2','default' => '25000'}}],'attrib' => {'expr' => '$UseBody2'}},{'type' => 't','content' => '

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
'}],'attrib' => {'name' => 'SECONDBODY','if' => '$_IsFirstSession'}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'BdpDimBody2x'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'BdpDimBody2y'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'BdpDimBody2z'}},{'type' => 't','content' => '

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
'}],'attrib' => {'name' => 'DIPOLEBODY2','if' => '$_IsFirstSession'}}],'attrib' => {'name' => 'PHYSICS PARAMETERS'}},{'type' => 'e','name' => 'commandgroup','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!! SOLAR PROBLEM TYPES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'BodyTDim','default' => '2.85E06'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'BodyRhoDim','default' => '1.50E8'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'qSun','default' => '25.0'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'tHeat','default' => '1.75'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'rHeat','default' => '1.0'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'SigmaHeat','default' => '4.5'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'DoInitRope','default' => 'F'}},{'type' => 'e','name' => 'if','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'CmeA','default' => '0.7'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'CmeR1','default' => '1.2'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'CmeR0','default' => '1.0'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'CmeA1','default' => '0.23'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'CmeAlpha','default' => '0.0'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'CmeRho1','default' => '2.5E-12'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'CmeRho2','default' => '2.0E-13'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','max' => '10','name' => 'ModulationRho','default' => '0.0'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','max' => '10','name' => 'ModulationP','default' => '0.0'}}],'attrib' => {'expr' => '$DoInitRope'}},{'type' => 't','content' => '
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
'}],'attrib' => {'name' => 'HELIOSPHERE','if' => '$_IsFirstSession and $_NameComp ne \'GM\''}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '-1','type' => 'real','name' => 'DtUpdateB0','default' => '0.0001'}},{'type' => 't','content' => '

#HELIOUPDATEB0
-1.0			DtUpdateB0 [s]

Set the frequency of updating the B0 field for the solar corona.
A negative value means that the B0 field is not updated.
'}],'attrib' => {'name' => 'HELIOUPDATEB0','if' => '$_NameComp ne \'GM\''}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'HelioDipoleStrength'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '-90','type' => 'real','max' => '90','name' => 'HelioDipoleTilt','default' => '0'}},{'type' => 't','content' => '

#HELIODIPOLE
-3.0                    HelioDipoleStrength [G]
 0.0                    HelioDipoleTilt     [deg]

! Variable HelioDipoleStrength defines the equatorial field strength in Gauss,
! while HelioDipoleTilt is the tilt relative to the ecliptic North 
! (negative sign means towards the planet) in degrees.
!
! Default value is HelioDipoleStrength = 0.0.
'}],'attrib' => {'name' => 'HELIODIPOLE','if' => '$_NameComp ne \'GM\''}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseInertialFrame','default' => 'T'}},{'type' => 'e','name' => 'if','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'UseRotatingBC','default' => 'T'}}],'attrib' => {'expr' => '$UseInertialFrame'}},{'type' => 't','content' => '

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
'}],'attrib' => {'alias' => 'INERTIAL','name' => 'HELIOROTATION','if' => '$_IsFirstSession and $_NameComp ne \'GM\''}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'logical','name' => 'DoSendMHD','default' => 'F'}},{'type' => 't','content' => '

#HELIOTEST
F			DoSendMHD

! If DoSendMHD is true, IH sends the real MHD solution to GM in the coupling.
! If DoSendMHD is false then the values read from the IMF file are sent,
! so there is no real coupling. Mostly used for testing the framework.
!
! Default value is true, ie. real coupling.
'}],'attrib' => {'name' => 'HELIOTEST','if' => '$_NameComp ne \'GM\''}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [{'type' => 'e','name' => 'option','content' => [],'attrib' => {'name' => 'Low','default' => 'T'}}],'attrib' => {'type' => 'string','name' => 'TypeCme','input' => 'select'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'CmeA','default' => '0.7'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'CmeR1','default' => '1.2'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'CmeR0','default' => '1.0'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'CmeA1','default' => '0.23'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'CmeAlpha','default' => '0.0'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'CmeRho1','default' => '2.5E-12'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'CmeRho2','default' => '2.0E-13'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'CmeB1Dim','default' => '1.0'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'CmeUErupt','default' => '4.0E5'}},{'type' => 't','content' => '
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
'}],'attrib' => {'name' => 'CME','if' => '$_IsFirstSession and $_NameComp ne \'GM\''}},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'tArcDim','default' => '1.0E6'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'RhoArcDim','default' => '1.0E-12'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'bArcDim','default' => '0.718144'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'ByArcDim','default' => '1.0E6'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'UzArcDim','default' => '5.0E3'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'Phi0Arc','default' => '0.5'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'real','name' => 'MuArc','default' => '1.3'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'ExpArc','default' => '3'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'WidthArc','default' => '0.5'}},{'type' => 't','content' => '
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
'}],'attrib' => {'name' => 'ARCADE','if' => '$_IsFirstSession and $_NameComp ne \'GM\''}}],'attrib' => {'name' => 'SOLAR PROBLEM TYPES'}},{'type' => 'e','name' => 'commandgroup','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!! COMET PROBLEM TYPE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'ProdRate'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'UrNeutral'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'AverageMass'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'IonizationRate'}},{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'min' => '0','type' => 'real','name' => 'kFriction'}},{'type' => 't','content' => '
#COMET
1.0E28		ProdRate    - Production rate (#/s)
1.0		UrNeutral   - neutral radial outflow velocity (km/s)
17.0		AverageMass - average particle mass (amu)
1.0E-6		IonizationRate (1/s)
1.7E-9		kFriction - ion-neutral friction rate coefficient (cm^3/s)

! Only used by problem_comet.  Defaults are as shown.
'}],'attrib' => {'name' => 'COMET','if' => '$_IsFirstSession'}}],'attrib' => {'name' => 'COMET PROBLEM TYPE'}},{'type' => 'e','name' => 'commandgroup','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! SCRIPT COMMANDS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'type' => 'e','name' => 'command','content' => [{'type' => 'e','name' => 'parameter','content' => [],'attrib' => {'type' => 'string','name' => 'NameIncludeFile','length' => '100','default' => 'Param/'}},{'type' => 't','content' => '

#INCLUDE
Param/SSS_3000		NameIncludeFile

! Include a library file from Param/ or any file from anywhere else.
'}],'attrib' => {'name' => 'INCLUDE'}},{'type' => 'e','name' => 'command','content' => [{'type' => 't','content' => '

#RUN

This command is only used in stand alone mode.

Run BATSRUS with the parameters above and then return for the next session
'}],'attrib' => {'name' => 'RUN','if' => '$_IsStandAlone'}},{'type' => 'e','name' => 'command','content' => [{'type' => 't','content' => '

#END

This command is only used in stand alone mode.

Run the executable with the parameters above and then stop.
In included files #END simply means the end of the included lines.
'}],'attrib' => {'name' => 'END','if' => '$_IsStandAlone'}}],'attrib' => {'name' => 'SCRIPT COMMANDS'}},{'type' => 'e','name' => 'rule','content' => [{'type' => 't','content' => '
	Either command #SOLARWIND or #UPSTREAM_INPUT_FILE must be used!
'}],'attrib' => {'expr' => '($SwRhoDim > 0) or $UseUpstreamInputFile or $_NameComp ne \'GM\''}},{'type' => 'e','name' => 'rule','content' => [{'type' => 't','content' => '
	Part implicit scheme requires more than 1 implicit block!
'}],'attrib' => {'expr' => '$MaxImplBlock>1 or not $UsePartImplicit or not $MaxImplBlock'}},{'type' => 'e','name' => 'rule','content' => [{'type' => 't','content' => '
	Full implicit scheme should be used with equal number of 
	explicit and implicit blocks!
'}],'attrib' => {'expr' => '$MaxImplBlock==$MaxBlock or not $UseFullImplicit'}}],'attrib' => {'name' => 'BATSRUS: GM, SC and IH Components'}}];