#^CFG FILE _FALSE_
$tree = [{'name' => 'commandList','type' => 'e','attrib' => {'name' => 'Global Magnetosphere and Inner Heliosphere'},'content' => [{'type' => 't','content' => '

List of MH (GM, IH and SC) commands used in the PARAM.in file




'},{'name' => 'set','type' => 'e','attrib' => {'name' => 'nI','type' => 'integer','value' => '$_GridSize[0]'},'content' => []},{'name' => 'set','type' => 'e','attrib' => {'name' => 'nJ','type' => 'integer','value' => '$_GridSize[1]'},'content' => []},{'name' => 'set','type' => 'e','attrib' => {'name' => 'nK','type' => 'integer','value' => '$_GridSize[2]'},'content' => []},{'name' => 'set','type' => 'e','attrib' => {'name' => 'MaxBlock','type' => 'integer','value' => '$_GridSize[3]'},'content' => []},{'name' => 'set','type' => 'e','attrib' => {'name' => 'MaxImplBlock','type' => 'integer','value' => '$_GridSize[4]'},'content' => []},{'name' => 'set','type' => 'e','attrib' => {'name' => 'MaxBlockALL','type' => 'integer','value' => '$_nProc and $MaxBlock and $_nProc*$MaxBlock'},'content' => []},{'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'STAND ALONE MODE'},'content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! STAND ALONE PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

'},{'name' => 'command','type' => 'e','attrib' => {'name' => 'NEWPARAM','if' => '$_IsStandAlone'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseNewParam','type' => 'logical','default' => 'T'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseNewAxes','type' => 'logical','default' => 'T'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoTimeAccurate','type' => 'logical','default' => 'T'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseCorotation','type' => 'logical','default' => 'T'},'content' => []},{'type' => 't','content' => '

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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'DESCRIPTION','if' => '$_IsStandAlone'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'StringDescription','type' => 'string','length' => '100'},'content' => []},{'type' => 't','content' => '

#DESCRIPTION
This is a test run for Jupiter with no rotation.

This command is only used in the stand alone mode.

The StringDescription string can be used to describe the simulation
for which the parameter file is written. The #DESCRIPTION command and
the StringDescription string are saved into the restart file,
which helps in identifying the restart files.

The default value is "Please describe me!", which is self explanatory.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'ECHO','if' => '$_IsStandAlone'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoEcho','type' => 'logical','default' => 'F'},'content' => []},{'type' => 't','content' => '

#ECHO
T                       DoEcho

This command is only used in the stand alone mode.

If the DoEcho variable is true, the input parameters are echoed back.
The default value for DoEcho is .false., but it is a good idea to
set it to true at the beginning of the PARAM.in file.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'PROGRESS','if' => '$_IsStandAlone'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DnProgressShort','min' => '-1','type' => 'integer','default' => '10'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DnProgressLong','min' => '-1','type' => 'integer','default' => '100'},'content' => []},{'type' => 't','content' => '
#PROGRESS
10			DnProgressShort
100			DnProgressLong

The frequency of short and long progress reports for BATSRUS in
stand alone mode. These are the defaults. Set -1-s for no progress reports.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'TIMEACCURATE','if' => '$_IsStandAlone'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoTimeAccurate','type' => 'logical','default' => 'T'},'content' => []},{'type' => 't','content' => '

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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'BEGIN_COMP','multiple' => 'T','if' => '$_IsStandAlone'},'content' => [{'type' => 't','content' => '

This command is allowed in stand alone mode only for sake of the 
test suite, which contains these commands when the framework is tested.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'END_COMP','multiple' => 'T','if' => '$_IsStandAlone'},'content' => [{'type' => 't','content' => '

This command is allowed in stand alone mode only for sake of the 
test suite, which contains these commands when the framework is tested.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'RUN','if' => '$_IsStandAlone'},'content' => [{'type' => 't','content' => '

#RUN

This command is only used in stand alone mode.

The #RUN command does not have any parameters. It signals the end
of the current session, and makes BATSRUS execute the session with
the current set of parameters. The parameters for the next session
start after the #RUN command. For the last session there is no
need to use the #RUN command, since the #END command or simply
the end of the PARAM.in file makes BATSRUS execute the last session.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'END'},'content' => [{'type' => 't','content' => '

#END

The #END command signals the end of the included file or the
end of the PARAM.in file. Lines following the #END command are
ignored. It is not required to use the #END command. The end
of the included file or PARAM.in file is equivalent with an 
#END command in the last line.
'}]}]},{'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'PLANET COMMANDS'},'content' => [{'type' => 't','content' => '
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

'},{'name' => 'command','type' => 'e','attrib' => {'name' => 'PLANET','if' => '$_IsFirstSession and $_IsStandAlone'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'NamePlanet','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'Earth','value' => 'EARTH/Earth/earth','default' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Saturn','value' => 'SATURN/Saturn/saturn'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'New'},'content' => []}]},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$NamePlanet eq \'New\''},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'RadiusPlanet','min' => '0','type' => 'real'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'MassPlanet','min' => '0','type' => 'real'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'OmegaPlanet','min' => '0','type' => 'real'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TiltRotation','min' => '0','type' => 'real'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypeBField','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'NONE'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'DIPOLE','default' => 'T'},'content' => []}]}]},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$TyepBField eq \'DIPOLE\''},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'MagAxisThetaGeo','min' => '0','type' => 'real','max' => '180'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'MagAxisPhiGeo','min' => '0','type' => 'real','max' => '360'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DipoleStrength','type' => 'real'},'content' => []}]},{'name' => 'rule','type' => 'e','attrib' => {'expr' => 'not $PlanetCommand'},'content' => [{'type' => 't','content' => '
		PLANET should precede $PlanetCommand
	'}]},{'type' => 't','content' => '

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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'ROTATIONAXIS','if' => '$_IsStandAlone'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'IsRotAxisPrimary','type' => 'logical','default' => 'T'},'content' => []},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$IsRotAxisPrimary'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'RotAxisTheta','min' => '0','type' => 'real','max' => '180'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'RotAxisPhi','min' => '0','type' => 'real','max' => '360'},'content' => []}]},{'name' => 'set','type' => 'e','attrib' => {'name' => 'PlanetCommand','type' => 'string','value' => 'ROTATIONAXIS'},'content' => []},{'type' => 't','content' => '

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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'ROTATION','if' => '$_IsStandAlone'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseRotation','type' => 'logical','default' => 'T'},'content' => []},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$UseRotation'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'RotationPeriod','type' => 'real'},'content' => []}]},{'name' => 'set','type' => 'e','attrib' => {'name' => 'PlanetCommand','type' => 'string','value' => 'MAGNETICAXIS'},'content' => []},{'type' => 't','content' => '

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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'MAGNETICAXIS','if' => '$_IsStandAlone'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'IsMagAxisPrimary','type' => 'logical','default' => 'T'},'content' => []},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$IsMagAxisPrimary'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'MagAxisTheta','min' => '0','type' => 'real','max' => '180'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'MagAxisPhi','min' => '0','type' => 'real','max' => '360'},'content' => []}]},{'name' => 'set','type' => 'e','attrib' => {'name' => 'PlanetCommand','type' => 'string','value' => 'MAGNETICAXIS'},'content' => []},{'type' => 't','content' => '

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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'DIPOLE','if' => '$_IsStandAlone'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DipoleStrength','type' => 'real'},'content' => []},{'type' => 't','content' => '

#DIPOLE
-3.11e-4		DipoleStrength [Tesla]

The DipoleStrength variable contains the
magnetic equatorial strength of the dipole magnetic field in Tesla.

The default value is the real dipole strength for the planet.
For the Earth the default is taken to be -31100 nT.
The sign is taken to be negative so that the magnetic axis can
point northward as usual.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'UPDATEB0','if' => '$_IsStandAlone'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DtUpdateB0','min' => '-1','type' => 'real','default' => '0.0001'},'content' => []},{'type' => 't','content' => '

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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'IDEALAXES','if' => '$_IsStandAlone'},'content' => [{'type' => 't','content' => '

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
'}]}]},{'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'USER DEFINED INPUT'},'content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!  USER DEFINED INPUT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

'},{'name' => 'command','type' => 'e','attrib' => {'name' => 'USER_FLAGS'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseUserInnerBcs','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseUserSource','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseUserPerturbation','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseUserOuterBcs','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseUserICs','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseUserSpecifyRefinement','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseUserLogFiles','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseUserWritePlot','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseUserAMR','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseUserEchoInput','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseUserB0','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseUserSetPhysConst','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseUserUpdateStates','type' => 'logical','default' => 'F'},'content' => []},{'type' => 't','content' => '

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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'USERINPUTBEGIN'},'content' => [{'type' => 't','content' => '

This command signals the beginning of the section of the file which 
is read by the subroutine user\\_read\\_inputs in the user\\_routines.f90 file.
The section ends with the #USERINPUTEND command. There is no XML based parameter
checking in the user section.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'USERINPUTEND'},'content' => [{'type' => 't','content' => '

This command signals the end of the section of the file which 
is read by the subroutine user\\_read\\_inputs in the user\\_routines.f90 file.
The section begins with the #USERINPUTBEGIN command. There is no XML based parameter
checking in the user section.
'}]}]},{'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'TESTING AND TIMING'},'content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  TESTING AND TIMING PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'name' => 'command','type' => 'e','attrib' => {'name' => 'TEST'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TestString','type' => 'string','length' => '100'},'content' => []},{'type' => 't','content' => '
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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'TESTIJK'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'iTest','min' => '-2','type' => 'integer','max' => '$nI+2'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'jTest','min' => '-2','type' => 'integer','max' => '$nJ+2'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'kTest','min' => '-2','type' => 'integer','max' => '$nK+2'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'iBlockTest','min' => '1','type' => 'integer','max' => '$MaxBlock'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'iProcTest','min' => '0','type' => 'integer'},'content' => []},{'type' => 't','content' => '
#TESTIJK
1                       iTest           (cell index for testing)
1                       jTest           (cell index for testing)
1                       kTest           (cell index for testing)
1                       BlockTest       (block index for testing)
0                       ProcTest        (processor index for testing)

! The location of test info in terms of indices, block and processor number.
! Note that the user should set #TESTIJK or #TESTXYZ, not both.  If both
! are set, the final one in the session will set the test point.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'TESTXYZ'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'xTest','min' => '$xMin','type' => 'real','max' => '$xMax'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'yTest','min' => '$yMin','type' => 'real','max' => '$yMax'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'zTest','min' => '$zMin','type' => 'real','max' => '$zMax'},'content' => []},{'type' => 't','content' => '
#TESTXYZ
1.5                     xTest           (X coordinate of cell for testing)
-10.5                   yTest           (Y coordinate of cell for testing)
-10.                    zTest           (Z coordinate of cell for testing)

! The location of test info in terms of coordinates.
! Note that the user should set #TESTIJK or #TESTXYZ, not both.  If both
! are set, the final one in the session will set the test point.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'TESTTIME'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nIterTest','min' => '-1','type' => 'integer','default' => '-1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TimeTest','min' => '-1','type' => 'real','default' => '1e30'},'content' => []},{'type' => 't','content' => '

#TESTTIME
-1                      nIterTest       (iteration number to start testing)
10.5                    TimeTest        (time to start testing in seconds)

! The time step and physical time to start testing.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'TESTVAR'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'iVarTest','type' => 'integer','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'Rho','value' => '1','default' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'RhoUx','value' => '2'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'RhoUy','value' => '3'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'RhoUz','value' => '4'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Bx','value' => '5'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'By','value' => '6'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Bz','value' => '7'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'e','value' => '8'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'p','value' => '9'},'content' => []}]},{'type' => 't','content' => '
#TESTVAR
1                       iVarTest

! Index of variable to be tested. Default is rho_="1", ie. density.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'TESTDIM'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'iVarTest','type' => 'integer','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'all','value' => '0'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'x','value' => '1','default' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'y','value' => '2'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'z','value' => '3'},'content' => []}]},{'type' => 't','content' => '
#TESTDIM
1                       iDimTest

! Index of dimension/direction to be tested. Default is X dimension.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'STRICT'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseStrict','type' => 'logical','default' => 'T'},'content' => []},{'type' => 't','content' => '
#STRICT
T                       UseStrict

! If true then stop when parameters are incompatible. If false, try to
! correct parameters and continue. Default is true, ie. strict mode
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'VERBOSE'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'iVarTest','type' => 'integer','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'errors and warnings only','value' => '-1'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'start and end of sessions','value' => '0'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'normal','value' => '1','default' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'calls on test processor','value' => '10'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'calls on all processors','value' => '100'},'content' => []}]},{'type' => 't','content' => '
#VERBOSE
-1                      lVerbose

! Verbosity level controls the amount of output to STDOUT. Default level is 1.
!   lVerbose .le. -1 only warnings and error messages are shown.
!   lVerbose .ge.  0 start and end of sessions is shown.
!   lVerbose .ge.  1 a lot of extra information is given.
!   lVerbose .ge. 10 all calls of set_oktest are shown for the test processor.
!   lVerbose .ge.100 all calls of set_oktest are shown for all processors.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'DEBUG'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoDebug','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoDebugGhost','type' => 'logical','default' => 'F'},'content' => []},{'type' => 't','content' => '
#DEBUG
F                       DoDebug         (use it as if(okdebug.and.oktest)...)
F                       DoDebugGhost    (parameter for show_BLK in library.f90)

! Excessive debug output can be controlled by the global okdebug parameter
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'CODEVERSION','if' => '$_IsFirstSession'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CodeVersion','min' => '0','type' => 'real','default' => '7.50'},'content' => []},{'type' => 't','content' => '
#CODEVERSION
7.50                    CodeVersion

! Cheks CodeVersion. Prints a WARNING if it differs from the CodeVersion
! defined in ModMain. Used in newer restart header files. 
! Should be given in PARAM.in when reading old restart files, 
! which do not have version info in the header file.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'EQUATION','if' => '$_IsFirstSession'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'NameEquation','type' => 'string','default' => 'MHD','length' => '100'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nVar','type' => 'integer','default' => '8'},'content' => []},{'type' => 't','content' => '
#EQUATION
MHD			NameEquation
8			nVar

! Define the equation name and the number of variables.
! If any of these do not agree with the values determined 
! by the code, BATSRUS stops with an error. Used in restart
! header files and can be given in PARAM.in as a check
! and as a description.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'PRECISION','if' => '$_IsFirstSession'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nByteReal','type' => 'integer','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'single precision (4)','value' => '4','default' => '$_nByteReal==4'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'double precision (8)','value' => '8','default' => '$_nByteReal==8'},'content' => []}]},{'name' => 'rule','type' => 'e','attrib' => {'expr' => '$nByteReal==$_nByteReal'},'content' => [{'type' => 't','content' => '
		nByteReal in file must agree with _nByteReal.
	'}]},{'type' => 't','content' => '

#PRECISION
8                       nByteReal

! Define the number of bytes in a real number. If it does not agree
! with the value determined by the code, BATSRUS stops with an error.
! This is a check, the internal value is calculated in parallel_setup.
! Used in latest restart header files to check binary compatibility.
! May be given in PARAM.in to enforce a certain precision.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'CHECKGRIDSIZE','if' => '$_IsFirstSession'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nI','min' => '$nI','type' => 'integer','max' => '$nI','default' => '$nI'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nJ','min' => '$nJ','type' => 'integer','max' => '$nJ','default' => '$nJ'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nK','min' => '$nK','type' => 'integer','max' => '$nK','default' => '$nK'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'MinBlockALL','min' => '1','type' => 'integer','max' => '$MaxBlockALL'},'content' => []},{'type' => 't','content' => '

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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'BLOCKLEVELSRELOADED'},'content' => [{'type' => 't','content' => '
#BLOCKLEVELSRELOADED

This command means that the restart file contains the information about
the minimum and maximum allowed refinement levels for each block.
This command is only used in the restart header file.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'TIMING'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseTiming','type' => 'logical','default' => 'T'},'content' => []},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$UseTiming'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Frequency','type' => 'integer','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'none','value' => '-3'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'final only','value' => '-2','default' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'end of sessions','value' => '-1'},'content' => []},{'name' => 'optioninput','type' => 'e','attrib' => {'name' => 'every X steps','min' => '1','default' => '100'},'content' => []}]},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nDepthTiming','min' => '-1','type' => 'integer','default' => '-1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypeTimingReport','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'cummulative','value' => 'cumm','default' => '1'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'list'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'tree'},'content' => []}]}]},{'type' => 't','content' => '
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
'}]}]},{'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'INITIAL AND BOUNDARY CONDITIONS'},'content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! MAIN INITIAL AND BOUNDARY CONDITION PARAMETERS  !!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'name' => 'command','type' => 'e','attrib' => {'required' => 'T','name' => 'PROBLEMTYPE','if' => '$_IsFirstSession'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'iProblem','type' => 'integer','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'Uniform','value' => '1'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Shock tube','value' => '2'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Heliosphere','value' => '3'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Comet','value' => '5'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Rotation','value' => '6'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Diffusion','value' => '7'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Earth','value' => '11','default' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Saturn','value' => '12'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Jupiter','value' => '13'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Venus','value' => '14'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Cylinder','value' => '21'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Sphere','value' => '22'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Arcade','value' => '25'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'CME','value' => '26'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Dissipation','value' => '30'},'content' => []}]},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypeDissipation','if' => '$iProblem==30','type' => 'string','length' => '20'},'content' => []},{'type' => 't','content' => '
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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'COORDSYSTEM','if' => '$_IsFirstSession'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypeCoordSystem','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'GeoSolarMagnetic, GSM','if' => '$_NameComp eq \'GM\'','value' => 'GSM','default' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'HelioGraphicInertial, HGI','if' => '$_NameComp ne \'GM\'','value' => 'HGI','default' => 'T'},'content' => []}]},{'type' => 't','content' => '

#COORDSYSTEM
GSM			TypeCoordSystem

! TypeCoordSystem defines the coordinate system for the component.
! Currently only one coordinate system is available for GM ("GSM")
! and one for IH or SC ("HGI"). In the near future "GSE" should be also
! an option for GM.
!
! Default is component dependent: "GSM" for GM and "HGI" for IH or SC.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'RESTARTINDIR','if' => '$_IsFirstSession'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'NameRestartInDir','type' => 'string','default' => 'GM/restartIN','length' => '100'},'content' => []},{'type' => 't','content' => '

#RESTARTINDIR
GM/restart_n5000	NameRestartInDir

! The NameRestartInDir variable contains the name of the directory
! where restart files are saved relative to the run directory.
! The directory should be inside the subdirectory with the name 
! of the component.
!
! Default value is "GM/restartIN".
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'NEWRESTART','if' => '$_IsFirstSession'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'IsRestartBFace','type' => 'logical','default' => 'F'},'content' => []},{'type' => 't','content' => '
! The restartIN/restart.H file always contains the #NEWRESTART command.
! This command is really used only in the restart headerfile.  Generally
! it is not inserted in a PARAM.in file by the user.
!
! Other than setting RestartBFace (used by the Constrained Transport scheme)
! The #NEWRESTART command also sets the following global variables:
!   DoRestart=.true.          read restart files
!   DoRestartGhost=.false.    no ghost cells are saved into restart file
!   DoRestartReals=.true.     only real numbers are saved in blk*.rst files
'}]},{'name' => 'command','type' => 'e','attrib' => {'required' => 'T','name' => 'GRID','if' => '$_IsFirstSession'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nRootX','min' => '1','type' => 'integer','default' => '2'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nRootY','min' => '1','type' => 'integer','default' => '1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nRootZ','min' => '1','type' => 'integer','default' => '1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'xMin','type' => 'real','default' => '-192.0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'xMax','min' => '$xMin','type' => 'real','default' => '  64.0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'yMin','type' => 'real','default' => ' -64.0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'yMax','min' => '$yMin','type' => 'real','default' => '  64.0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'zMin','type' => 'real','default' => ' -64.0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'zMax','min' => '$zMin','type' => 'real','default' => '  64.0'},'content' => []},{'type' => 't','content' => '
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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'OUTERBOUNDARY'},'content' => [{'name' => 'foreach','type' => 'e','attrib' => {'values' => 'TypeBcEast,TypeBcWest,TypeBcSouth,TypeBcNorth,TypeBcBot,TypeBcTop','name' => 'Side'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => '$Side','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'coupled'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'fixed/inflow','default' => '$Side ne \'TypeBcEast\''},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'float/outflow','default' => '$Side eq \'TypeBcEast\''},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'reflect'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'periodic'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'vary'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'shear'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'linetied'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'raeder'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'arcadetop'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'arcadebot'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'arcadebotcont'},'content' => []}]}]},{'name' => 'rule','type' => 'e','attrib' => {'expr' => 'not($TypeBcEast eq \'periodic\' xor $TypeBcWest eq \'periodic\')'},'content' => [{'type' => 't','content' => '
	! East and west BCs must be both periodic or neither
	'}]},{'name' => 'rule','type' => 'e','attrib' => {'expr' => 'not($TypeBcSouth eq \'periodic\' xor $TypeBcNorth eq \'periodic\')'},'content' => [{'type' => 't','content' => '
	! South and North BCs must be both periodic or neither
	'}]},{'name' => 'rule','type' => 'e','attrib' => {'expr' => 'not($TypeBcBot eq \'periodic\' xor $TypeBcTop eq \'periodic\')'},'content' => [{'type' => 't','content' => '
	! Bottom and top BCs must be both periodic or neither
	'}]},{'type' => 't','content' => '
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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'INNERBOUNDARY'},'content' => [{'type' => 't','content' => '
! Inner boundary types for body 1 and body 2
	'},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypeInnerBc','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'reflect'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'float'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'fixed'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'ionosphere','default' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'ionosphereB0','value' => 'ionosphereB0/ionosphereb0'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'ionospherefloat'},'content' => []}]},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$UseBody2'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypeInnerBcBody2','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'reflect','default' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'float'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'fixed'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'ionosphere'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'ionosphereB0'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'ionospherefloat'},'content' => []}]}]},{'name' => 'rule','type' => 'e','attrib' => {'expr' => 'not($TypeInnerBcBody2 =~ /ionosphere/)'},'content' => [{'type' => 't','content' => '
! Note: for the second body COROTATION AND AN IONOSPHERIC BOUNDARY DO NOT WORK.
	'}]},{'type' => 't','content' => '

#INNERBOUNDARY
ionosphere              InnerBCType

ionosphere              InnerBCTypeBody2  !read only if UseBody2=.true. 

!This command should appear after the #SECONDBODY command if using 2 bodies
! Note:  for the second body COROTATION AND AN IONOSPHERIC BOUNDARY DO NOT
!        WORK.
! Default boundary for the second body is reflect.


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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'EXTRABOUNDARY'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseExtraBoundary','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$UseExtraBoundary'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypeBcExtra','type' => 'string'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoFixExtraboundary','type' => 'logical','default' => 'F'},'content' => []}]}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'FACEOUTERBC'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'MaxBoundary','min' => '0','type' => 'integer','max' => '6','default' => '0'},'content' => []},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$MaxBoundary >= 1'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoFixOuterBoundary','type' => 'logical','default' => 'F'},'content' => []}]},{'type' => 't','content' => '
#FACEOUTERBC
0              MaxBoundary            
F              DoFixOuterBoundary)    !read only for MaxBoundary>=East_(=1).
! if MaxBoundary>=East_(=1) then the outer boundaries with
! the number of boundary being between East_ and MaxBoundary
! are treated using set_BCs.f90 subroutines instead of set_outerBCs.f90 
! if DoFixOuterBoundary==.true., there is no resolution
! change along the outer boundaries with the number of
! of boundary being between East_ and MaxBoundary
'}]}]},{'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'INITIAL TIME'},'content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! INITIAL TIME AND STEP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'name' => 'command','type' => 'e','attrib' => {'name' => 'STARTTIME','if' => '$_IsFirstSession','alias' => 'SETREALTIME'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'year','type' => 'integer','default' => '2000'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'month','min' => '1','type' => 'integer','max' => '12','default' => '3'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'day','min' => '1','type' => 'integer','max' => '31','default' => '21'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'hour','min' => '0','type' => 'integer','max' => '23','default' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'minute','min' => '0','type' => 'integer','max' => '59','default' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'second','min' => '0','type' => 'integer','max' => '59','default' => '0'},'content' => []},{'type' => 't','content' => '
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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'TIMESIMULATION','if' => '$_IsFirstSession'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'tSimulation','min' => '0','type' => 'real','default' => '0.0'},'content' => []},{'type' => 't','content' => '

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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'NSTEP','if' => '$_IsFirstSession'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nStep','min' => '0','type' => 'integer','default' => '0'},'content' => []},{'type' => 't','content' => '

#NSTEP
100			nStep

! Set nStep for the component. Typically used in the restart.H header file.
! Generally it is not inserted in a PARAM.in file by the user.
!
! The default is nStep=0 as the starting time step with no restart.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'NPREVIOUS','if' => '$_IsFirstSession'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nPrevious','min' => '-1','type' => 'integer','default' => '-1'},'content' => []},{'type' => 't','content' => '

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
'}]}]},{'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'TIME INTEGRATION'},'content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  TIME INTEGRATION PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'name' => 'command','type' => 'e','attrib' => {'name' => 'TIMESTEPPING'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nStage','type' => 'integer','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'value' => '1','default' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'value' => '2'},'content' => []}]},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CflExpl','min' => '0','type' => 'real','max' => '1','default' => '0.8'},'content' => []},{'type' => 't','content' => '

#TIMESTEPPING
2                       nStage
0.80                    CflExpl

! Parameters for explicit time integration.
! Default is 1 stage and CflExpl=0.8
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'FIXEDTIMESTEP'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseDtFixed','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DtFixedDim','min' => '0','if' => '$UseDtFixed','type' => 'real','default' => '1.0'},'content' => []},{'type' => 't','content' => '
#FIXEDTIMESTEP
T                       UseDtFixed
10.                     DtFixedDim [sec] (read if UseDtFixed is true)

! Default is UseDtFixed=.false. Effective only if DoTimeAccurate is true.
! If UseDtFixed is true, the time step is fixed to DtFixedDim.
!
! This is useful for debugging explicit schemes.

! The real application is, however, for implicit and partially
! implicit/local schemes.

'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'PARTLOCAL'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UsePartLocal','type' => 'logical','default' => 'F'},'content' => []},{'type' => 't','content' => '
#PARTLOCAL
T               UsePartLocal

! Default is UsePartLocal=.false. If UsePartLocal is true and the
! run is time accurate, then the blocks selected as "implicit"
! by the criteria defined in #STEPPINGCRITERIA are not used to
! calculate the time step, and all cells are advanced with the
! smaller of the stable and the global time steps.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'IMPLICIT'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UsePointImplicit','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UsePartImplicit','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseFullImplicit','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CflImpl','min' => '0','if' => '$UsePartImplicit or $UseFullImplicit','type' => 'real','default' => '100'},'content' => []},{'name' => 'rule','type' => 'e','attrib' => {'expr' => '$UsePointImplicit + $UsePartImplicit + $UseFullImplicit <= 1'},'content' => [{'type' => 't','content' => '
	At most one of these logicals can be true!
	'}]},{'type' => 't','content' => '

#IMPLICIT
F               UsePointImplicit   
F               UsePartImplicit
F               UseFullImplicit
100.0           CflImpl (read if UsePartImplicit or UseFullImplicit is true)

! Default is false for all logicals. Only one of them can be set to true!
! The CFL number is used in the implicit blocks of the fully or partially
! implicit schemes. Ignored if UseDtFixed is true.
'}]}]},{'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'IMPLICIT PARAMETERS'},'content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! PARAMETERS FOR FULL AND PART IMPLICIT TIME INTEGRATION !!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'name' => 'command','type' => 'e','attrib' => {'name' => 'IMPLICITCRITERIA','alias' => 'STEPPINGCRITERIA'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypeImplCrit','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'Time step','value' => 'dt','default' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Radial distance','value' => 'r/R'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Test block','value' => 'test'},'content' => []}]},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'rImplicit','min' => '0','if' => '$TypeImplCrit eq \'R\'','type' => 'real'},'content' => []},{'type' => 't','content' => '
! Both #IMPLICITCRITERIA and #STEPPINGCRITERIA are acceptable.
! Only effective if PartImplicit or PartLocal is true in a time accurate run.
! Default value is ImplCritType=\'dt\'.
!
! The options are
!
! If     (TypeImplCrit ==\'dt\'  ) then blocks with DtBLK .gt. DtFixed
! ElseIf (TypeImplCrit ==\'R\'   ) then blocks with rMinBLK .lt. rImplicit
! ElseIf (TypeImplCrit ==\'test\') then block iBlockTest on processor iProcTest
!
! are handled with local/implicit scheme.
! DtFixed must be defined in #FIXEDTIMESTEP
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'IMPLSTEP'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'ImplCoeff','min' => '0','type' => 'real','max' => '1','default' => '1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseBdf2','type' => 'logical','default' => 'T'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseSourceImpl','type' => 'logical','default' => 'F'},'content' => []},{'type' => 't','content' => '
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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'IMPLSCHEME'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nOrderImpl','type' => 'integer','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => '1','default' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => '2'},'content' => []}]},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypeFluxImpl','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'Roe','value' => 'Roe/roe/1'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Rusanov','value' => 'Rusanov/rusanov/2/TVDLF','default' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Linde','value' => 'Linde/linde/3/HLLEL'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Sokolov','value' => 'Sokolov/sokolov/4/AW'},'content' => []}]},{'type' => 't','content' => '
#IMPLSCHEME
1               nOrderImpl
Rusanov         TypeFluxImpl

! Default values are shown, ie. first order Rusanov scheme.
! This defines the scheme used in the implicit part.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'NEWTON'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseConservativeImplicit','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseNewton','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$UseNewton'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseNewMatrix','type' => 'logical','default' => 'T'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'MaxIterNewton','min' => '1','type' => 'integer','default' => '10'},'content' => []}]},{'type' => 't','content' => '
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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'JACOBIAN'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypeJacobian','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'Preconditioned','value' => 'prec','default' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'No preconditioning','value' => 'free'},'content' => []}]},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'JacobianEps','min' => '0','type' => 'real','max' => '1.e-5','default' => '$doublePrecision ? 1.e-12 : 1.e-6'},'content' => []},{'type' => 't','content' => '
#JACOBIAN
prec            TypeJacobian (prec, free)
1.E-12          JacobianEps

! The Jacobian matrix is always calculated with a matrix free approach,
! however it can be preconditioned  (\'prec\'), or not (\'free\')
! Default value is TypeJacobian=\'prec\'.
! JacobianEps contains the machine round off error for numerical derivatives.
! The default value is 1.E-12 for 8 byte reals and 1.E-6 for 4 byte reals.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'PRECONDITIONER'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypePrecondSide','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'left'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'symmetric','default' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'right'},'content' => []}]},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypePrecond','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'MBILU','default' => 'T'},'content' => []}]},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'GustafssonPar','min' => '0','type' => 'real','max' => '1','default' => '0.5'},'content' => []},{'type' => 't','content' => '
#PRECONDITIONER
symmetric       TypePrecondSide (left, symmetric, right)
MBILU           TypePrecond (MBILU)
0.5             GustafssonPar (0. no modification, 1. full modification)

! Default parameters are shown. Right preconditioning does not affect
! the normalization of the residual. The Gustafsson parameter determines
! how much the MBILU preconditioner is modified. The default 0.5 value
! means a relaxed modification.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'KRYLOV'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypeKrylov','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'gmres','default' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'bicgstab'},'content' => []}]},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypeInitKrylov','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => '0','value' => 'nul','default' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'previous','value' => 'old'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'explicit'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'scaled explicit','value' => 'explicit'},'content' => []}]},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'ErrorMaxKrylov','min' => '0','type' => 'real','max' => '0.1','default' => '0.001'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'MaxMatvecKrylov','min' => '1','type' => 'integer','default' => '100'},'content' => []},{'type' => 't','content' => '
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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'KRYLOVSIZE'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nKrylovVector','min' => '1','type' => 'integer','default' => 'MaxMatvecKrylov'},'content' => []},{'type' => 't','content' => '
#KRYLOVSIZE
10		nKrylovVector

! The number of Krylov vectors only matters for GMRES (TypeKrylov=\'gmres\').
! If GMRES does not converge within nKrylovVector iterations, it needs
! a restart, which usually degrade its convergence rate and robustness.
! So nKrylovVector should exceed the number of iterations, on the other
! hand it should not exceed the maximum number of iterations MaxMatvecKrylov.
! On the other hand the dynamically allocated memory is also proportional 
! to nKrylovVector. The default is nKrylovVector=MaxMatvecKrylov (in #KRYLOV)
! which can be overwritten by #KRYLOVSIZE after the #KRYLOV command (if any).
'}]}]},{'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'STOPPING CRITERIA'},'content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! STOPPING CRITERIA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

The commands in this group only work in stand alone mode.

'},{'name' => 'command','type' => 'e','attrib' => {'required' => '$_IsStandAlone','name' => 'STOP','if' => '$_IsStandAlone'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'MaxIteration','min' => '-1','type' => 'integer','default' => '-1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'tSimulationMax','min' => '-1','type' => 'real','default' => '-1'},'content' => []},{'type' => 't','content' => '

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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'CHECKSTOPFILE','if' => '$_IsStandAlone'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoCheckStopFile','type' => 'logical','default' => 'T'},'content' => []},{'type' => 't','content' => '

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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'CPUTIMEMAX','if' => '$_IsStandAlone'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CpuTimeMax','min' => '-1','type' => 'real','default' => '-1'},'content' => []},{'type' => 't','content' => '

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
'}]}]},{'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'OUTPUT PARAMETERS'},'content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  OUTPUT PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'name' => 'command','type' => 'e','attrib' => {'name' => 'RESTARTOUTDIR'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'NameRestartOutDir','type' => 'string','default' => 'GM/restartOUT','length' => '100'},'content' => []},{'type' => 't','content' => '

#RESTARTOUTDIR
GM/restart_n5000	NameRestartOutDir

! The NameRestartOutDir variable contains the name of the directory
! where restart files are saved relative to the run directory.
! The directory should be inside the subdirectory with the name 
! of the component.
!
! Default value is "GM/restartOUT".
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'SAVERESTART'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'SaveRestart','type' => 'logical','default' => 'T'},'content' => []},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$SaveRestart'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DnRestart','min' => '-1','type' => 'integer','default' => '-1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DtRestart','min' => '-1','type' => 'real','default' => '-1'},'content' => []}]},{'type' => 't','content' => '
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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'PLOTDIR'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'NamePlotDir','type' => 'string','default' => 'GM/IO2','length' => '100'},'content' => []},{'type' => 't','content' => '

The NamePlotDir variable contains the name of the directory
where plot files and logfiles are saved relative to the run directory.
The directory should be inside the subdirectory with the name
of the component.

Default value is "GM/IO2".
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'SAVELOGFILE'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoSaveLogfile','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$DoSaveLogfile'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'StringLog','min' => '1','type' => 'strings','max' => '4'},'content' => [{'name' => 'part','type' => 'e','attrib' => {'required' => 'T','name' => 'TypeLogVar','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'MHD vars. dimensional','value' => 'MHD','default' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Flux vars. dimensional','value' => 'FLX'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Raw vars. dimensional','value' => 'RAW'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Set vars. dimensional','value' => 'VAR'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'MHD vars. scaled','value' => 'mhd','default' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Flux vars. scaled','value' => 'flx'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Raw vars. scaled','value' => 'raw'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Set vars. scaled','value' => 'var'},'content' => []}]},{'name' => 'part','type' => 'e','attrib' => {'required' => 'F','name' => 'TypeTime','multiple' => 'T','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'none','exclusive' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'step'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'date'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'time'},'content' => []}]}]},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DnOutput','min' => '-1','type' => 'integer','default' => '1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DtOutput','min' => '-1','type' => 'real','default' => '-1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'NameLogVars','if' => '$TypeLogVar =~ /var/i','type' => 'string','length' => '100'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'StringLogRadii','min' => '1','if' => '($TypeLogVar=~/flx/i or $NameLogVars=~/flx/i)','type' => 'strings','max' => '10','length' => '100'},'content' => [{'name' => 'part','type' => 'e','attrib' => {'name' => 'LogRadii','multiple' => 'T','min' => '$rBody','type' => 'real'},'content' => []}]}]},{'type' => 't','content' => '
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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'SATELLITE','if' => '$_IsFirstSession'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nSatellite','min' => '0','type' => 'integer','default' => '0'},'content' => []},{'name' => 'for','type' => 'e','attrib' => {'from' => '1','to' => '$nSatellite'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'StringSatellite','min' => '1','type' => 'strings','max' => '5'},'content' => [{'name' => 'part','type' => 'e','attrib' => {'required' => 'T','name' => 'TypeSatelliteVar','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'MHD vars. dimensional','value' => 'MHD','default' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'All vars. dimensional','value' => 'FUL'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Set vars. dimensional','value' => 'VAR'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'MHD vars. scaled','value' => 'mhd'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'All vars. scaled','value' => 'ful'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Set vars. scaled','value' => 'var'},'content' => []}]},{'name' => 'part','type' => 'e','attrib' => {'required' => 'F','name' => 'TypeTrajectory','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'file','default' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'equation','value' => 'eqn'},'content' => []}]},{'name' => 'part','type' => 'e','attrib' => {'required' => 'F','name' => 'TypeTime','multiple' => 'T','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'none','exclusive' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'step'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'date'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'time'},'content' => []}]}]},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DnOutput','min' => '-1','type' => 'integer','default' => '1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DtOutput','min' => '-1','type' => 'real','default' => '-1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'NameSatellite','type' => 'string','length' => '100'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'NameSatelliteVars','if' => '$TypeSatelliteVar =~ /\\bvar\\b/i','type' => 'string','length' => '100'},'content' => []}]},{'type' => 't','content' => '
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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'SAVEPLOT'},'content' => [{'type' => 't','content' => '
! plot_string must contain the following 3 parts in arbitrary order
...
	'},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nPlotFile','min' => '0','type' => 'integer','max' => '100','default' => '0'},'content' => []},{'name' => 'for','type' => 'e','attrib' => {'from' => '1','to' => '$nPlotFile'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'plotString','min' => '3','type' => 'strings','max' => '3'},'content' => [{'name' => 'part','type' => 'e','attrib' => {'required' => 'T','name' => 'plotform','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'TECPLOT','value' => 'tec'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'IDL','value' => 'idl'},'content' => []}]},{'name' => 'part','type' => 'e','attrib' => {'required' => 'T','name' => 'plotarea','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'value' => '3d'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'value' => 'x=0'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'value' => 'y=0','default' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'value' => 'z=0'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'value' => 'sph'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'value' => 'los'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'value' => 'lin'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'if' => '$plotform =~ /\\bidl\\b/','value' => 'cut'},'content' => []}]},{'name' => 'part','type' => 'e','attrib' => {'required' => 'T','name' => 'plotvar','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'MHD vars. dimensional','value' => 'MHD'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'All vars. dimensional','value' => 'FUL'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Raw vars. dimensional','value' => 'RAW'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Ray tracing vars. dim.','value' => 'RAY'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Flux vars. dimensional','value' => 'FLX'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Solar vars. dimensional','value' => 'SOL'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Select dimensional vars.','value' => 'VAR'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'MHD vars. scaled','value' => 'mhd'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'All vars. scaled','value' => 'ful'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Raw vars. scaled','value' => 'raw'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Ray tracing vars. scaled','value' => 'ray'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Flux vars. scaled','value' => 'flx'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Solar vars. scaled','value' => 'sol'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Select scaled vars.','value' => 'var'},'content' => []}]}]},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DnOutput','min' => '-1','type' => 'integer'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DtOutput','min' => '-1','type' => 'real'},'content' => []},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$plotarea =~ /\\bcut\\b/'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'xMinCut','type' => 'real'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'xMaxCut','min' => '$xMinCut','type' => 'real'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'yMinCut','type' => 'real'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'yMaxCut','min' => '$yMinCut','type' => 'real'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'zMinCut','type' => 'real'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'zMaxCut','min' => '$zMinCut','type' => 'real'},'content' => []}]},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'radius','min' => '0','if' => '$plotarea =~ /\\bsph\\b/','type' => 'real','default' => '10'},'content' => []},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$plotarea =~ /\\blos\\b/'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'LosVectorX','type' => 'real','default' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'LosVectorY','type' => 'real','default' => '0.0001'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'LosVectorZ','type' => 'real','default' => '1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'xSizeImage','min' => '0','type' => 'real','default' => '20'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'ySizeImage','min' => '0','type' => 'real','default' => '20'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'xOffset','type' => 'real','default' => '10'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'yOffset','type' => 'real','default' => '10'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'rOccult','min' => '1','type' => 'real','default' => '2.5'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'MuLimbDarkening','min' => '0','type' => 'real','max' => '1','default' => '0.5'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nPixX','min' => '2','type' => 'integer','default' => '200'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nPixY','min' => '2','type' => 'integer','default' => '200'},'content' => []}]},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$plotarea =~ /\\blin\\b/'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'NameLine','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'Advected B','value' => 'A'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'B','default' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'U'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'J'},'content' => []}]},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'IsSingleLine','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nLine','min' => '1','type' => 'integer','max' => '20','default' => '1'},'content' => []},{'name' => 'for','type' => 'e','attrib' => {'from' => '1','to' => '$nLine'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'xStartLine','type' => 'real'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'yStartLine','type' => 'real'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'zStartLine','type' => 'real'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'IsParallel','type' => 'logical'},'content' => []}]}]},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DxPlot','min' => '-1.0','if' => '($plotform=~/\\bidl\\b/ and $plotarea!~/\\b(sph|los|lin)\\b/)','type' => 'real','default' => '-1.0'},'content' => []},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$plotvar =~ /\\bvar\\b/i'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'plotVars','type' => 'string','length' => '100'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'plotPars','type' => 'string','length' => '100'},'content' => []}]}]},{'type' => 't','content' => '
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

'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'SAVEBINARY'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoSaveBinary','type' => 'logical','default' => 'T'},'content' => []},{'type' => 't','content' => '
#SAVEBINARY
T			DoSaveBinary   used only for \'idl\' plot file

! Default is .true. Saves unformatted IO2/*.idl files if true. 
! This is the recommended method, because it is fast and accurate.
! The only advantage of saving IO2/*.idl in formatted text files is
! that it can be processed on another machine or with a different 
! (lower) precision. For example PostIDL.exe may be compiled with 
! single precision to make IO2/*.out files smaller, while BATSRUS.exe is 
! compiled in double precision, to make results more accurate.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'SAVEPLOTSAMR'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'SavePlotsAmr','type' => 'logical','default' => 'F'},'content' => []},{'type' => 't','content' => '
#SAVEPLOTSAMR
F			savePlotsAMR to save plots before each AMR

! Default is save_plots_amr=.false.
'}]}]},{'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'AMR PARAMETERS'},'content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  AMR PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'name' => 'command','type' => 'e','attrib' => {'name' => 'AMRINIT','if' => '$_IsFirstSession'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'InitialRefineType','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'default','default' => '1'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'all'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'none'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => '3Dbodyfocus'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'spherefocus'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'magnetosphere'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'points'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'helio_init'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'helio_z=4'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'all_then_focus'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'cme'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'points'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'mag_new'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'magnetosphere'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'magneto_fine'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'magneto12'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'magnetosaturn'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'magnetojupiter'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'paleo'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'comet'},'content' => []}]},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'InitialRefineLevel','min' => '0','type' => 'integer','default' => '4'},'content' => []},{'type' => 't','content' => '
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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'AMRINITPHYSICS','if' => '$_IsFirstSession'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nRefineLevelIC','min' => '0','type' => 'integer','default' => '0'},'content' => []},{'type' => 't','content' => '
#AMRINITPHYSICS
3			nRefineLevelIC

! Defines number of physics (initial condition) based AMR-s AFTER the 
! geometry based initial AMR-s defined by #AMRINIT were done.
! Only useful if the initial condition has a non-trivial analytic form.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'AMRLEVELS'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'minBlockLevel','min' => '-1','type' => 'integer','default' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'maxBlockLevel','min' => '-1','type' => 'integer','default' => '99'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'FixBodyLevel','type' => 'logical','default' => 'F'},'content' => []},{'type' => 't','content' => '
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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'AMRRESOLUTION'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'minCellDx','min' => '-1','type' => 'real','default' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'maxCellDx','min' => '-1','type' => 'real','default' => '99999'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'FixBodyLevel','type' => 'logical','default' => 'F'},'content' => []},{'type' => 't','content' => '
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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'AMR'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DnRefine','min' => '-1','type' => 'integer','default' => '-1'},'content' => []},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$DnRefine>0'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoAutoRefine','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$DoAutoRefine'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'percentCoarsen','min' => '0','type' => 'real','max' => '100','default' => '20'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'percentRefine','min' => '0','type' => 'real','max' => '100','default' => '20'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'maxTotalBlocks','min' => '1','type' => 'integer','default' => '99999'},'content' => []}]}]},{'type' => 't','content' => '
#AMR
2001			dnRefine (frequency in terms of total steps n_step)
T			DoAutoRefine 
0.			percentCoarsen
0.			percentRefine
99999			maxTotalBlocks

! Default for dn_refine is -1, ie. no run time refinement.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'AMRCRITERIA'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nRefineCrit','type' => 'integer','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => '1'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => '2'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => '3','default' => '1'},'content' => []}]},{'name' => 'for','type' => 'e','attrib' => {'from' => '1','to' => '$nRefineCrit'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypeRefine','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'grad T','value' => 'gradt/gradT'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'grad P','value' => 'gradp/gradP'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'grad log(Rho)','value' => 'gradlogrho'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'grad log(p)','value' => 'gradlogP/gradlogp'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'grad E','value' => 'gradE'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'curl U','value' => 'curlV/curlv/curlU/curlu'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'curl B','value' => 'curlB/curlb'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'div U','value' => 'divU/divu/divV/divv'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'divB','value' => 'divb/divB'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'vAlfven','value' => 'Valfven/vAlfven/valfven'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'heliospheric beta','value' => 'heliobeta'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'flux'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'heliospheric current sheet','value' => 'heliocurrentsheet'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'rCurrents','value' => 'rcurrents/Rcurrents'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Transient','value' => 'transient/Transient'},'content' => []}]},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$TypeRefine =~ /transient/i'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypeTransient','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'P_dot','value' => 'p_dot/P_dot'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'T_dot','value' => 't_dot/T_dot'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Rho_dot','value' => 'rho_dot/Rho_dot','default' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'RhoU_dot','value' => 'RhoU_dot/rhou_dot'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Rho_2nd_1','value' => 'Rho_2nd_1/rho_2nd_1'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Rho_2nd_2','value' => 'Rho_2nd_2/rho_2nd_2'},'content' => []}]},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseSunEarth','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$UseSunEarth'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'xEarth','type' => 'real'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'yEarth','type' => 'real'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'zEarth','type' => 'real'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'InvD2Ray','type' => 'real'},'content' => []}]}]}]},{'type' => 't','content' => '
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
'}]}]},{'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'SCHEME PARAMETERS'},'content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  SCHEME PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'name' => 'command','type' => 'e','attrib' => {'name' => 'SCHEME'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nOrder','type' => 'integer','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => '1','default' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => '2'},'content' => []}]},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypeFlux','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'Roe','value' => 'Roe/roe/1'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Rusanov','value' => 'Rusanov/rusanov/2/TVDLF','default' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Linde','value' => 'Linde/linde/3/HLLEL'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Sokolov','value' => 'Sokolov/sokolov/4/AW'},'content' => []}]},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$nOrder == 2'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypeLimiter','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'minmod','default' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'beta'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'mc'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'LSG'},'content' => []}]},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'LimiterBeta','min' => '1','if' => '$TypeLimiter eq \'beta\'','type' => 'real','max' => '2','default' => '1.2'},'content' => []}]},{'type' => 't','content' => '
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
! \'mc\'		- Monotonized Central limiter is sharper but less robust
! \'LSG\'		- Least Squares Gradient: robust but expensive multiD limiter 
! \'beta\'        - Beta limiter
!
! Possible values for LimiterBeta are between 1.0 and 2.0 : 
!  LimiterBeta = 1.0 is the same as the minmod limiter
!  LimiterBeta = 2.0 is the same as the superbee limiter
!  LimiterBeta = 1.2 is the recommended value
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'NONCONSERVATIVE'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseNonConservative','type' => 'logical','default' => 'T'},'content' => []},{'type' => 't','content' => '
#NONCONSERVATIVE
T		UseNonConservative

! For Earth the default is using non-conservative equations 
! (close to the body).
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'CONSERVATIVECRITERIA'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nConservCrit','min' => '0','type' => 'integer','max' => '3','default' => '1'},'content' => []},{'name' => 'for','type' => 'e','attrib' => {'from' => '1','to' => '$nConservCrit'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypeConservCrit_I','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'radius','value' => 'r/R/radius/Radius','default' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'p','value' => 'p/P'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'grad P','value' => 'gradp/GradP'},'content' => []}]},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'rConserv','min' => '$rBody','if' => '$TypeConservCrit_I =~ /^r|radius$/i','type' => 'real','default' => '2*$rBody'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'pCoeffConserv','min' => '0','if' => '$TypeConservCrit_I =~ /^p$/i','type' => 'real','default' => '0.05'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'GradPCoeffConserv','min' => '0','if' => '$TypeConservCrit_I =~ /gradp/i','type' => 'real','default' => '0.1'},'content' => []}]},{'type' => 't','content' => '
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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'UPDATECHECK'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseUpdateCheck','type' => 'logical','default' => 'T'},'content' => []},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$UseUpdateCheck'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'rhoMin','min' => '0','type' => 'real','max' => '100','default' => '40'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'rhoMax','min' => '100','type' => 'real','default' => '400'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'pMin','min' => '0','type' => 'real','max' => '100','default' => '40'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'pMax','min' => '100','type' => 'real','default' => '400'},'content' => []}]},{'type' => 't','content' => '
#UPDATECHECK
T			UseUpdateCheck
40.			rhoMin[%]
400.			rhoMax[%]
40.			pMin[%]
400.			pMax[%]

! Default values are shown.  This will adjust the timestep so that
! density and pressure cannot change by more than the given percentages
! in a single timestep.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'PROLONGATION'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nOrderProlong','type' => 'integer','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => '1','default' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => '2'},'content' => []}]},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypeProlong','if' => '$nOrderProlong==2','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'left-right','value' => 'lr','default' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'central'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'minmod'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'left-right extrapolate','value' => 'lr2'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'central    extrapolate','value' => 'central2'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'minmod     extrapolate','value' => 'minmod2'},'content' => []}]},{'type' => 't','content' => '
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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'MESSAGEPASS','alias' => 'OPTIMIZE'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypeMessagePass','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'm_p_cell FACES ONLY','value' => 'allopt','default' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'm_p_cell','value' => 'all'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'm_p_dir FACES ONLY','value' => 'opt'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'm_p_dir group by directions','value' => 'dir'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'm_p_dir group by faces     ','value' => 'face'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'm_p_dir group by kind and face','value' => 'min'},'content' => []}]},{'type' => 't','content' => '
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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'BORIS'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseBorisCorrection','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'BorisClightFactor','min' => '0','if' => '$UseBorisCorrection','type' => 'real','max' => '1','default' => '1'},'content' => []},{'type' => 't','content' => '
#BORIS
T			UseBorisCorrection
1.0			BorisClightFactor !Only if UseBorisCorrection is true

! Default is boris_correction=.false.
! Use semi-relativistic MHD equations with speed of light reduced by
! the BorisClightFactor. Set BorisClightFactor=1.0 for true semi-relativistic
! MHD. Gives the same steady state as normal MHD analytically, but there
! can be differences due to discretization errors. 
! You can use either Boris or BorisSimple but not both. 
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'BORISSIMPLE'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseBorisSimple','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'BorisClightFactor','min' => '0','if' => '$UseBorisSimple','type' => 'real','max' => '1','default' => '1'},'content' => []},{'type' => 't','content' => '
#BORISSIMPLE
T			UseBorisSimple
0.05			BorisClightFactor !Only if UseBorisSimple is true

! Default is UseBorisSimple=.false. 
! Use simplified semi-relativistic MHD with speed of light reduced by the
! BorisClightFactor. This is only useful with BorisClightFactor less than 1.
! Should give the same steady state as normal MHD, but there can be a
! difference due to discretization errors.
! You can use either Boris or BorisSimple but not both. 
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'DIVB'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseDivbSource','type' => 'logical','default' => 'T'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseDivbDiffusion','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseProjection','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseConstrainB','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'rule','type' => 'e','attrib' => {'expr' => '$UseDivbSource or $UseDivbDiffusion or $UseProjection or $UseConstrainB'},'content' => [{'type' => 't','content' => '
	! At least one of the options should be true.
	'}]},{'name' => 'rule','type' => 'e','attrib' => {'expr' => 'not($UseProjection and ($UseDivbSource or $UseDivbDiffusion or $UseConstrainB))'},'content' => [{'type' => 't','content' => '
	! If UseProjection is true, all others should be false.
	'}]},{'name' => 'rule','type' => 'e','attrib' => {'expr' => 'not($UseConstrainB and ($UseDivbSource or $UseDivbDiffusion or $UseProjection))'},'content' => [{'type' => 't','content' => '
	! If UseConstrainB is true, all others should be false.
	'}]},{'type' => 't','content' => '
	
#DIVB
T			UseDivbSource
F			UseDivbDiffusion	
F			UseProjection           
F			UseConstrainB           

! Default values are shown above.
! If UseProjection is true, all others should be false.
! If UseConstrainB is true, all others should be false.
! At least one of the options should be true.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'DIVBSOURCE'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseB0Source','type' => 'logical','default' => 'T'},'content' => []},{'type' => 't','content' => '
#DIVBSOURCE
T			UseB0Source

! Add extra source terms related to the non-zero divergence and curl of B0.
! Default is true.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'DIVBDIFFUSION'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DivbDiffCoeff','min' => '0','type' => 'real','max' => '1','default' => '0.1666667'},'content' => []},{'type' => 't','content' => '
#DIVBDIFFUSION
0.1666667		DivbDiffCoeff

! Default value is shown above. 1.0/6.0
! If divb_diffcoeff .gt. 0.5 then cfl .lt. 0.5/DivbDiffCoeff is required!
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'PROJECTION'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypeProjectIter','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'Conjugate Gradients','value' => 'cg','default' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'BiCGSTAB','value' => 'bicgstab'},'content' => []}]},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypeProjectStop','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'Relative norm','value' => 'rel','default' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Maximum error','value' => 'max'},'content' => []}]},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'RelativeLimit','min' => '0','type' => 'real','max' => '1','default' => '0.1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'AbsoluteLimit','min' => '0','type' => 'real','default' => '0.0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'MaxMatvec','min' => '1','type' => 'integer','default' => '50'},'content' => []},{'type' => 't','content' => '
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

'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'CORRECTP'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'pRatioLow','min' => '0','type' => 'real','max' => '1','default' => '0.01'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'pRatioHigh','min' => '$pRatioLow','type' => 'real','max' => '1','default' => '0.1'},'content' => []},{'type' => 't','content' => '
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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'RAYTRACE'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseAccurateIntegral','type' => 'logical','default' => 'T'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseAccurateTrace','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DtExchangeRay','min' => '0.01','type' => 'real','max' => '60','default' => '0.1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DnRaytrace','min' => '1','type' => 'integer','default' => '1'},'content' => []},{'type' => 't','content' => '

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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'IM'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TauCoupleIm','min' => '0','type' => 'real'},'content' => []},{'type' => 't','content' => '

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
'}]}]},{'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'PHYSICS PARAMETERS'},'content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  PHYSICS PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'name' => 'command','type' => 'e','attrib' => {'name' => 'GAMMA','if' => '$_IsFirstSession'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Gamma','min' => '1','type' => 'real','default' => '1.6666666667'},'content' => []},{'type' => 't','content' => '
#GAMMA
1.6666666667		g

! Above value is the default.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'SHOCKTUBE'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'RhoLeft','min' => '0','type' => 'real','default' => '1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UnLeft','type' => 'real','default' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Ut1Left','type' => 'real','default' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Ut2Left','type' => 'real','default' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'BnLeft','type' => 'real','default' => '0.75'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Bt1Left','type' => 'real','default' => '1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Bt2Left','type' => 'real','default' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'pRight','min' => '0','type' => 'real','default' => '1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'RhoRight','min' => '0','type' => 'real','default' => '0.125'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UnRight','type' => 'real','default' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Ut1Right','type' => 'real','default' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Ut2Right','type' => 'real','default' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'BnRight','type' => 'real','default' => '0.75'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Bt1Right','type' => 'real','default' => '-1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Bt2Right','type' => 'real','default' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'pRight','min' => '0','type' => 'real','default' => '0.1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'ShockSlope','type' => 'real','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'no rotation','value' => '0','default' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'value' => '0.25'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'value' => '0.3333333333333'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'value' => '0.5'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'value' => '1'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'value' => '2'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'value' => '3'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'value' => '4'},'content' => []}]},{'type' => 't','content' => '
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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'SOLARWIND','if' => '$_IsFirstSession'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'SwRhoDim','min' => '0','type' => 'real','default' => '5'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'SwTDim','min' => '0','type' => 'real','default' => '181712.175'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'SwUxDim','type' => 'real','max' => '0','default' => '-400'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'SwUyDim','type' => 'real','default' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'SwUzDim','type' => 'real','default' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'SwBxDim','type' => 'real','default' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'SwByDim','type' => 'real','default' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'SwBzDim','type' => 'real','default' => '5'},'content' => []},{'type' => 't','content' => '
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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'UPSTREAM_INPUT_FILE'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseUpstreamInputFile','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$UseUpstreamInputFile'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'NameUpstreamFile','type' => 'string','length' => '100'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'SatelliteYPos','type' => 'real','default' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'SatelliteZPos','type' => 'real','default' => '0'},'content' => []}]},{'type' => 't','content' => '
#UPSTREAM_INPUT_FILE
T			UseUpstreamInputFile (rest of parameters read if true)
IMF.dat                 NameUpstreamFile
0.0                     SatelliteYPos
0.0                     SatelliteZPos

! UseUpstreamInputFile - default is false
! UpstreamFileName     - user specified input file
! Satellite_Y_Pos      - not yet used
! Satellite_Z_Pos      - not yet used
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'BODY','if' => '$_IsFirstSession','alias' => 'MAGNETOSPHERE'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseBody','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$UseBody'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'rBody','min' => '0','type' => 'real','default' => '3'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'rCurrents','min' => '-1','type' => 'real','default' => '4'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'BodyRhoDim','min' => '0','type' => 'real','default' => '1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'BodyTDim','min' => '0','type' => 'real','default' => '10000'},'content' => []}]},{'type' => 't','content' => '
#BODY
T			UseBody (rest of parameters read if true)
3.0			rBody
4.0			rCurrents
1.0			BodyRhoDim (/ccm) density for fixed BC for rho_BLK
10000.0			BodyTDim (K) temperature for fixed BC for P_BLK

! Default values depend on problem_type.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'GRAVITY','if' => '$_IsFirstSession'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseGravity','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'iDirGravity','if' => '$UseGravity','type' => 'integer','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'central mass','value' => '0','default' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'X direction','value' => '1'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Y direction','value' => '2'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Z direction','value' => '3'},'content' => []}]},{'type' => 't','content' => '
#GRAVITY
T			UseGravity (rest of parameters read if true)
0			GravityDir (0 - central, 1 - X, 2 - Y, 3 - Z direction)

! Default values depend on problem_type.  

! When a second body is used the gravity direction for the second body
!  is independent of the GravityDir value.  Gravity due to the second body
!  is radially inward toward the second body.

'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'FACONDUCTIVITYMODEL'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UsePhysicalFAConductance','type' => 'logical','default' => 'F'},'content' => []},{'type' => 't','content' => '
#FACONDUCTIVITYMODEL
F			UsePhysicalFAConductance

Default value is shown.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'MASSLOADING'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseMassLoading','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoAccelerateMassLoading','type' => 'logical','default' => 'F'},'content' => []},{'type' => 't','content' => '
#MASSLOADING
F			UseMassLoading
F			AccelerateMassLoading
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'HEATFLUX'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseHeatFlux','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseSpitzerForm','type' => 'logical','default' => 'T'},'content' => []},{'name' => 'if','type' => 'e','attrib' => {'expr' => 'not $UseSpitzerForm'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Kappa0Heat','type' => 'real','default' => '1.23E-11'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Kappa0Heat','type' => 'real','default' => '2.5'},'content' => []}]},{'type' => 't','content' => '
#HEATFLUX
T		UseHeatFlux
F		UseSpitzerForm
1.23E-11	Kappa0Heat [W/m/K]	! Only if not UseSpitzerForm
2.50E+00	ExponentHeat [-]	! Only if not UseSpitzerForm
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'RESISTIVEFLUX'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseResistFlux','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseSpitzerForm','type' => 'logical','default' => 'T'},'content' => []},{'name' => 'if','type' => 'e','attrib' => {'expr' => 'not $UseSpitzerForm'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypeResist','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'localized','value' => 'Localized/localized'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'constant','value' => 'Constant/constant','default' => 'T'},'content' => []}]},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$TypeResist =~ /localized/i'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Eta0Resist','type' => 'real','default' => '9.69953E+8'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Alpha0Resist','type' => 'real','default' => '150'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'yShiftResist','type' => 'real','default' => '0.5'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TimeInitRise','type' => 'real','default' => '0.05'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TimeConstLev','type' => 'real','default' => '1'},'content' => []}]}]},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseAnomResist','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$UseAnomResist'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Eta0AnomResist','type' => 'real','default' => '1.93991E+09'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'EtaAnomMaxResist','type' => 'real','default' => '1.93991E+10'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'ThresholdFactorResist','type' => 'real','default' => '1'},'content' => []}]},{'type' => 't','content' => '
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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'TESTDISSMHD','if' => '$_IsFirstSession'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseDefaultUnits','type' => 'logical','default' => 'T'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Grav0Diss','type' => 'real','default' => '2.635620E-02'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Beta0Diss','type' => 'real','default' => '1.640000E-01'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Length0Diss','type' => 'real','default' => '1.500000E+06'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Time0Diss','type' => 'real','default' => '1.159850E+01'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Rho0Diss','type' => 'real','default' => '5.019000E-11'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Tem0Diss','type' => 'real','default' => '1.000000E+05'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Theta0Diss','type' => 'real','default' => '6.000000E-01'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Delta0Diss','type' => 'real','default' => '2.500000E+01'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'EpsilonDiss','type' => 'real','default' => '7.000000E+00'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'RhoDifDiss','type' => 'real','default' => '4.500000E+00'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'yShiftDiss','type' => 'real','default' => '4.000000E-01'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'ScaleHeightDiss','type' => 'real','default' => '5.000000E-01'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'ScaleFactorDiss','type' => 'real','default' => '1.159850E+01'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'BZ0iss','type' => 'real','default' => '5.000000E-01'},'content' => []},{'type' => 't','content' => '
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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'SECONDBODY','if' => '$_IsFirstSession'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseBody2','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$UseBody2'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'rBody2','min' => '0','type' => 'real','default' => '0.1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'xBody2','min' => '$xMin','type' => 'real','max' => '$xMax','default' => '-40'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'yBody2','min' => '$yMin','type' => 'real','max' => '$yMax','default' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'zBody2','min' => '$zMin','type' => 'real','max' => '$zMax','default' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'rCurrents2','min' => '$rBody2','type' => 'real','default' => '1.3*$rBody2'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'RhoDimBody2','min' => '0','type' => 'real','default' => '5'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'tDimBody2','min' => '0','type' => 'real','default' => '25000'},'content' => []}]},{'type' => 't','content' => '

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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'DIPOLEBODY2','if' => '$_IsFirstSession'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'BdpDimBody2x','type' => 'real'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'BdpDimBody2y','type' => 'real'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'BdpDimBody2z','type' => 'real'},'content' => []},{'type' => 't','content' => '

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
'}]}]},{'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'SOLAR PROBLEM TYPES'},'content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!! SOLAR PROBLEM TYPES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'name' => 'command','type' => 'e','attrib' => {'name' => 'HELIOSPHERE','if' => '$_IsFirstSession and $_NameComp ne \'GM\''},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'BodyTDim','min' => '0','type' => 'real','default' => '2.85E06'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'BodyRhoDim','min' => '0','type' => 'real','default' => '1.50E8'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'qSun','min' => '0','type' => 'real','default' => '25.0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'tHeat','min' => '0','type' => 'real','default' => '1.75'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'rHeat','min' => '0','type' => 'real','default' => '1.0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'SigmaHeat','min' => '0','type' => 'real','default' => '4.5'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoInitRope','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$DoInitRope'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CmeA','min' => '0','type' => 'real','default' => '0.7'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CmeR1','min' => '0','type' => 'real','default' => '1.2'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CmeR0','min' => '0','type' => 'real','default' => '1.0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CmeA1','min' => '0','type' => 'real','default' => '0.23'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CmeAlpha','type' => 'real','default' => '0.0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CmeRho1','min' => '0','type' => 'real','default' => '2.5E-12'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CmeRho2','min' => '0','type' => 'real','default' => '2.0E-13'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'ModulationRho','min' => '0','type' => 'real','max' => '10','default' => '0.0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'ModulationP','min' => '0','type' => 'real','max' => '10','default' => '0.0'},'content' => []}]},{'type' => 't','content' => '
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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'HELIODIPOLE','if' => '$_NameComp ne \'GM\''},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'HelioDipoleStrength','type' => 'real'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'HelioDipoleTilt','min' => '-90','type' => 'real','max' => '90','default' => '0'},'content' => []},{'type' => 't','content' => '

#HELIODIPOLE
-3.0                    HelioDipoleStrength [G]
 0.0                    HelioDipoleTilt     [deg]

! Variable HelioDipoleStrength defines the equatorial field strength in Gauss,
! while HelioDipoleTilt is the tilt relative to the ecliptic North 
! (negative sign means towards the planet) in degrees.
!
! Default values are ???
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'HELIOROTATION','if' => '$_IsFirstSession and $_NameComp ne \'GM\'','alias' => 'INERTIAL'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseInertialFrame','type' => 'logical','default' => 'T'},'content' => []},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$UseInertialFrame'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseRotatingBC','type' => 'logical','default' => 'T'},'content' => []}]},{'type' => 't','content' => '

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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'HELIOTEST','if' => '$_NameComp ne \'GM\''},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoSendMHD','type' => 'logical','default' => 'F'},'content' => []},{'type' => 't','content' => '

#HELIOTEST
F			DoSendMHD

! If DoSendMHD is true, IH sends the real MHD solution to GM in the coupling.
! If DoSendMHD is false then the values read from the IMF file are sent,
! so there is no real coupling. Mostly used for testing the framework.
!
! Default value is true, ie. real coupling.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'CME','if' => '$_IsFirstSession and $_NameComp ne \'GM\''},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypeCme','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'Low','default' => 'T'},'content' => []}]},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CmeA','min' => '0','type' => 'real','default' => '0.7'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CmeR1','min' => '0','type' => 'real','default' => '1.2'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CmeR0','min' => '0','type' => 'real','default' => '1.0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CmeA1','min' => '0','type' => 'real','default' => '0.23'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CmeAlpha','type' => 'real','default' => '0.0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CmeRho1','min' => '0','type' => 'real','default' => '2.5E-12'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CmeRho2','min' => '0','type' => 'real','default' => '2.0E-13'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CmeB1Dim','type' => 'real','default' => '1.0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CmeUErupt','min' => '0','type' => 'real','default' => '4.0E5'},'content' => []},{'type' => 't','content' => '
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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'ARCADE','if' => '$_IsFirstSession and $_NameComp ne \'GM\''},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'tArcDim','min' => '0','type' => 'real','default' => '1.0E6'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'RhoArcDim','min' => '0','type' => 'real','default' => '1.0E-12'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'bArcDim','min' => '0','type' => 'real','default' => '0.718144'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'ByArcDim','type' => 'real','default' => '1.0E6'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UzArcDim','type' => 'real','default' => '5.0E3'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Phi0Arc','type' => 'real','default' => '0.5'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'MuArc','type' => 'real','default' => '1.3'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'ExpArc','min' => '0','type' => 'real','default' => '3'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'WidthArc','min' => '0','type' => 'real','default' => '0.5'},'content' => []},{'type' => 't','content' => '
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
'}]}]},{'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'COMET PROBLEM TYPE'},'content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!! COMET PROBLEM TYPE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'name' => 'command','type' => 'e','attrib' => {'name' => 'COMET','if' => '$_IsFirstSession'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'ProdRate','min' => '0','type' => 'real'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UrNeutral','min' => '0','type' => 'real'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'AverageMass','min' => '0','type' => 'real'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'IonizationRate','min' => '0','type' => 'real'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'kFriction','min' => '0','type' => 'real'},'content' => []},{'type' => 't','content' => '
#COMET
1.0E28		ProdRate    - Production rate (#/s)
1.0		UrNeutral   - neutral radial outflow velocity (km/s)
17.0		AverageMass - average particle mass (amu)
1.0E-6		IonizationRate (1/s)
1.7E-9		kFriction - ion-neutral friction rate coefficient (cm^3/s)

! Only used by problem_comet.  Defaults are as shown.
'}]}]},{'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'SCRIPT COMMANDS'},'content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! SCRIPT COMMANDS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'name' => 'command','type' => 'e','attrib' => {'name' => 'INCLUDE'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'NameIncludeFile','type' => 'string','default' => 'Param/','length' => '100'},'content' => []},{'type' => 't','content' => '

#INCLUDE
Param/SSS_3000		NameIncludeFile

! Include a library file from Param/ or any file from anywhere else.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'RUN','if' => '$_IsStandAlone'},'content' => [{'type' => 't','content' => '

#RUN

This command is only used in stand alone mode.

Run BATSRUS with the parameters above and then return for the next session
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'END','if' => '$_IsStandAlone'},'content' => [{'type' => 't','content' => '

#END

This command is only used in stand alone mode.

Run the executable with the parameters above and then stop.
In included files #END simply means the end of the included lines.
'}]}]},{'name' => 'rule','type' => 'e','attrib' => {'expr' => '($SwRhoDim > 0) or $UseUpstreamInputFile'},'content' => [{'type' => 't','content' => '
	Either command #SOLARWIND or #UPSTREAM_INPUT_FILE must be used!
'}]},{'name' => 'rule','type' => 'e','attrib' => {'expr' => '$MaxImplBlock>1 or not $UsePartImplicit or not $MaxImplBlock'},'content' => [{'type' => 't','content' => '
	Part implicit scheme requires more than 1 implicit block!
'}]},{'name' => 'rule','type' => 'e','attrib' => {'expr' => '$MaxImplBlock==$MaxBlock or not $UseFullImplicit'},'content' => [{'type' => 't','content' => '
	Full implicit scheme should be used with equal number of 
	explicit and implicit blocks!
'}]}]}];