#^CFG FILE _FALSE_
$tree = [{'content' => [{'content' => '

List of MH (GM, IH and SC) commands used in the PARAM.in file


','type' => 't'},{'content' => [],'attrib' => {'value' => '$_GridSize[0]','type' => 'integer','name' => 'nI'},'name' => 'set','type' => 'e'},{'content' => [],'attrib' => {'value' => '$_GridSize[1]','type' => 'integer','name' => 'nJ'},'name' => 'set','type' => 'e'},{'content' => [],'attrib' => {'value' => '$_GridSize[2]','type' => 'integer','name' => 'nK'},'name' => 'set','type' => 'e'},{'content' => [],'attrib' => {'value' => '$_GridSize[3]','type' => 'integer','name' => 'MaxBlock'},'name' => 'set','type' => 'e'},{'content' => [],'attrib' => {'value' => '$_GridSize[4]','type' => 'integer','name' => 'MaxImplBlock'},'name' => 'set','type' => 'e'},{'content' => [],'attrib' => {'value' => '$_nProc and $MaxBlock and $_nProc*$MaxBlock','type' => 'integer','name' => 'MaxBlockALL'},'name' => 'set','type' => 'e'},{'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! STAND ALONE PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

','type' => 't'},{'content' => [{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseNewParam','default' => 'T'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseNewAxes','default' => 'T'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'logical','name' => 'DoTimeAccurate','default' => 'T'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseCorotation','default' => 'T'},'name' => 'parameter','type' => 'e'},{'content' => '

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
algorithms found in\\\\
share/Library/src/CON\\_axes, the planet data is set and
stored by share/Library/src/CON\\_planet, and magnetic field information and
mapping is provided by share/Library/src/CON\\_planet_field, and the rotational speed
of the planet is calculated using $v_{\\phi}=\\Omega \\times r$.

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
','type' => 't'}],'attrib' => {'name' => 'NEWPARAM','if' => '$_IsStandAlone'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [{'content' => [],'attrib' => {'name' => 'SC','default' => 'T','if' => '$_NameComp eq \'SC\''},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'IH','default' => 'T','if' => '$_NameComp eq \'IH\''},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'GM','default' => 'T','if' => '$_NameComp eq \'GM\''},'name' => 'option','type' => 'e'}],'attrib' => {'input' => 'select','type' => 'string','name' => 'NameComp'},'name' => 'parameter','type' => 'e'},{'content' => '

#COMPONENT
GM			NameComp

This command is only used in the stand alone mode.

The NameComp variable contains the two-character component ID
for the component which BATSRUS is representing.
If NameComp does not agree with the value of the NameThisComp
variable, BATSRUS stops with an error message.
This command is saved into the restart header file for consistency check.

There is no default value: if the command is not given, the component ID is not checked.
','type' => 't'}],'attrib' => {'name' => 'COMPONENT'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'length' => '100','type' => 'string','name' => 'StringDescription'},'name' => 'parameter','type' => 'e'},{'content' => '

#DESCRIPTION
This is a test run for Jupiter with no rotation.

This command is only used in the stand alone mode.

The StringDescription string can be used to describe the simulation
for which the parameter file is written. The #DESCRIPTION command and
the StringDescription string are saved into the restart file,
which helps in identifying the restart files.

The default value is "Please describe me!", which is self explanatory.
','type' => 't'}],'attrib' => {'name' => 'DESCRIPTION','if' => '$_IsStandAlone'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'logical','name' => 'DoEcho','default' => 'F'},'name' => 'parameter','type' => 'e'},{'content' => '

#ECHO
T                       DoEcho

This command is only used in the stand alone mode.

If the DoEcho variable is true, the input parameters are echoed back.
The default value for DoEcho is .false., but it is a good idea to
set it to true at the beginning of the PARAM.in file.
','type' => 't'}],'attrib' => {'name' => 'ECHO','if' => '$_IsStandAlone'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'integer','name' => 'DnProgressShort','default' => '10','min' => '-1'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'integer','name' => 'DnProgressLong','default' => '100','min' => '-1'},'name' => 'parameter','type' => 'e'},{'content' => '
#PROGRESS
10			DnProgressShort
100			DnProgressLong

The frequency of short and long progress reports for BATSRUS in
stand alone mode. These are the defaults. Set -1-s for no progress reports.
','type' => 't'}],'attrib' => {'name' => 'PROGRESS','if' => '$_IsStandAlone'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'logical','name' => 'DoTimeAccurate','default' => 'T'},'name' => 'parameter','type' => 'e'},{'content' => '

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
','type' => 't'}],'attrib' => {'name' => 'TIMEACCURATE','if' => '$_IsStandAlone'},'name' => 'command','type' => 'e'},{'content' => [{'content' => '

This command is allowed in stand alone mode only for sake of the 
test suite, which contains these commands when the framework is tested.
','type' => 't'}],'attrib' => {'name' => 'BEGIN_COMP','multiple' => 'T','if' => '$_IsStandAlone'},'name' => 'command','type' => 'e'},{'content' => [{'content' => '

This command is allowed in stand alone mode only for sake of the 
test suite, which contains these commands when the framework is tested.
','type' => 't'}],'attrib' => {'name' => 'END_COMP','multiple' => 'T','if' => '$_IsStandAlone'},'name' => 'command','type' => 'e'},{'content' => [{'content' => '

#RUN

This command is only used in stand alone mode.

The #RUN command does not have any parameters. It signals the end
of the current session, and makes BATSRUS execute the session with
the current set of parameters. The parameters for the next session
start after the #RUN command. For the last session there is no
need to use the #RUN command, since the #END command or simply
the end of the PARAM.in file makes BATSRUS execute the last session.
','type' => 't'}],'attrib' => {'name' => 'RUN','if' => '$_IsStandAlone'},'name' => 'command','type' => 'e'},{'content' => [{'content' => '

#END

The #END command signals the end of the included file or the
end of the PARAM.in file. Lines following the #END command are
ignored. It is not required to use the #END command. The end
of the included file or PARAM.in file is equivalent with an 
#END command in the last line.
','type' => 't'}],'attrib' => {'name' => 'END'},'name' => 'command','type' => 'e'}],'attrib' => {'name' => 'STAND ALONE MODE'},'name' => 'commandgroup','type' => 'e'},{'content' => [{'content' => '
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

','type' => 't'},{'content' => [{'content' => [{'content' => [],'attrib' => {'value' => 'EARTH/Earth/earth','name' => 'Earth','default' => 'T'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'SATURN/Saturn/saturn','name' => 'Saturn'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'New'},'name' => 'option','type' => 'e'}],'attrib' => {'input' => 'select','type' => 'string','name' => 'NamePlanet'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'real','name' => 'RadiusPlanet','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'MassPlanet','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'OmegaPlanet','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'TiltRotation','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'name' => 'NONE'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'DIPOLE','default' => 'T'},'name' => 'option','type' => 'e'}],'attrib' => {'input' => 'select','type' => 'string','name' => 'TypeBField'},'name' => 'parameter','type' => 'e'}],'attrib' => {'expr' => '$NamePlanet eq \'New\''},'name' => 'if','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'max' => '180','type' => 'real','name' => 'MagAxisThetaGeo','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '360','type' => 'real','name' => 'MagAxisPhiGeo','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'DipoleStrength'},'name' => 'parameter','type' => 'e'}],'attrib' => {'expr' => '$TyepBField eq \'DIPOLE\''},'name' => 'if','type' => 'e'},{'content' => [{'content' => '
		PLANET should precede $PlanetCommand
	','type' => 't'}],'attrib' => {'expr' => 'not $PlanetCommand'},'name' => 'rule','type' => 'e'},{'content' => '

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
','type' => 't'}],'attrib' => {'name' => 'PLANET','if' => '$_IsFirstSession and $_IsStandAlone'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'logical','name' => 'IsRotAxisPrimary','default' => 'T'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'max' => '180','type' => 'real','name' => 'RotAxisTheta','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '360','type' => 'real','name' => 'RotAxisPhi','min' => '0'},'name' => 'parameter','type' => 'e'}],'attrib' => {'expr' => '$IsRotAxisPrimary'},'name' => 'if','type' => 'e'},{'content' => [],'attrib' => {'value' => 'ROTATIONAXIS','type' => 'string','name' => 'PlanetCommand'},'name' => 'set','type' => 'e'},{'content' => '

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
','type' => 't'}],'attrib' => {'name' => 'ROTATIONAXIS','if' => '$_IsStandAlone'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseRotation','default' => 'T'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'real','name' => 'RotationPeriod'},'name' => 'parameter','type' => 'e'}],'attrib' => {'expr' => '$UseRotation'},'name' => 'if','type' => 'e'},{'content' => [],'attrib' => {'value' => 'MAGNETICAXIS','type' => 'string','name' => 'PlanetCommand'},'name' => 'set','type' => 'e'},{'content' => '

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
','type' => 't'}],'attrib' => {'name' => 'ROTATION','if' => '$_IsStandAlone'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'logical','name' => 'IsMagAxisPrimary','default' => 'T'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'max' => '180','type' => 'real','name' => 'MagAxisTheta','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '360','type' => 'real','name' => 'MagAxisPhi','min' => '0'},'name' => 'parameter','type' => 'e'}],'attrib' => {'expr' => '$IsMagAxisPrimary'},'name' => 'if','type' => 'e'},{'content' => [],'attrib' => {'value' => 'MAGNETICAXIS','type' => 'string','name' => 'PlanetCommand'},'name' => 'set','type' => 'e'},{'content' => '

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
','type' => 't'}],'attrib' => {'name' => 'MAGNETICAXIS','if' => '$_IsStandAlone'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'real','name' => 'DipoleStrength'},'name' => 'parameter','type' => 'e'},{'content' => '

#DIPOLE
-3.11e-4		DipoleStrength [Tesla]

The DipoleStrength variable contains the
magnetic equatorial strength of the dipole magnetic field in Tesla.

The default value is the real dipole strength for the planet.
For the Earth the default is taken to be -31100 nT.
The sign is taken to be negative so that the magnetic axis can
point northward as usual.
','type' => 't'}],'attrib' => {'name' => 'DIPOLE','if' => '$_IsStandAlone'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'real','name' => 'DtUpdateB0','default' => '0.0001','min' => '-1'},'name' => 'parameter','type' => 'e'},{'content' => '

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
','type' => 't'}],'attrib' => {'name' => 'UPDATEB0','if' => '$_IsStandAlone'},'name' => 'command','type' => 'e'},{'content' => [{'content' => '

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
','type' => 't'}],'attrib' => {'name' => 'IDEALAXES','if' => '$_IsStandAlone'},'name' => 'command','type' => 'e'}],'attrib' => {'name' => 'PLANET COMMANDS'},'name' => 'commandgroup','type' => 'e'},{'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!  USER DEFINED INPUT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

','type' => 't'},{'content' => [{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseUserInnerBcs','default' => 'F'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseUserSource','default' => 'F'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseUserPerturbation','default' => 'F'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseUserOuterBcs','default' => 'F'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseUserICs','default' => 'F'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseUserSpecifyRefinement','default' => 'F'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseUserLogFiles','default' => 'F'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseUserWritePlot','default' => 'F'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseUserAMR','default' => 'F'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseUserEchoInput','default' => 'F'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseUserB0','default' => 'F'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseUserSetPhysConst','default' => 'F'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseUserUpdateStates','default' => 'F'},'name' => 'parameter','type' => 'e'},{'content' => '

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
','type' => 't'}],'attrib' => {'name' => 'USER_FLAGS'},'name' => 'command','type' => 'e'},{'content' => [{'content' => '

This command signals the beginning of the section of the file which 
is read by the subroutine user\\_read\\_inputs in the user\\_routines.f90 file.
The section ends with the #USERINPUTEND command. There is no XML based parameter
checking in the user section.
','type' => 't'}],'attrib' => {'name' => 'USERINPUTBEGIN'},'name' => 'command','type' => 'e'},{'content' => [{'content' => '

This command signals the end of the section of the file which 
is read by the subroutine user\\_read\\_inputs in the user\\_routines.f90 file.
The section begins with the #USERINPUTBEGIN command. There is no XML based parameter
checking in the user section.
','type' => 't'}],'attrib' => {'name' => 'USERINPUTEND'},'name' => 'command','type' => 'e'}],'attrib' => {'name' => 'USER DEFINED INPUT'},'name' => 'commandgroup','type' => 'e'},{'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  TESTING AND TIMING PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
','type' => 't'},{'content' => [{'content' => [],'attrib' => {'length' => '100','type' => 'string','name' => 'TestString'},'name' => 'parameter','type' => 'e'},{'content' => '
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
','type' => 't'}],'attrib' => {'name' => 'TEST'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'max' => '$nI+2','type' => 'integer','name' => 'iTest','min' => '-2'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '$nJ+2','type' => 'integer','name' => 'jTest','min' => '-2'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '$nK+2','type' => 'integer','name' => 'kTest','min' => '-2'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '$MaxBlock','type' => 'integer','name' => 'iBlockTest','min' => '1'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'integer','name' => 'iProcTest','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => '
#TESTIJK
1                       iTest           (cell index for testing)
1                       jTest           (cell index for testing)
1                       kTest           (cell index for testing)
1                       iBlockTest      (block index for testing)
0                       iProcTest       (processor index for testing)

! The location of test info in terms of indices, block and processor number.
! Note that the user should set #TESTIJK or #TESTXYZ, not both.  If both
! are set, the final one in the session will set the test point.
','type' => 't'}],'attrib' => {'name' => 'TESTIJK'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'max' => '$xMax','type' => 'real','name' => 'xTest','min' => '$xMin'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '$yMax','type' => 'real','name' => 'yTest','min' => '$yMin'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '$zMax','type' => 'real','name' => 'zTest','min' => '$zMin'},'name' => 'parameter','type' => 'e'},{'content' => '
#TESTXYZ
1.5                     xTest           (X coordinate of cell for testing)
-10.5                   yTest           (Y coordinate of cell for testing)
-10.                    zTest           (Z coordinate of cell for testing)

! The location of test info in terms of coordinates.
! Note that the user should set #TESTIJK or #TESTXYZ, not both.  If both
! are set, the final one in the session will set the test point.
','type' => 't'}],'attrib' => {'name' => 'TESTXYZ'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'integer','name' => 'nIterTest','default' => '-1','min' => '-1'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'TimeTest','default' => '1e30','min' => '-1'},'name' => 'parameter','type' => 'e'},{'content' => '

#TESTTIME
-1                      nIterTest       (iteration number to start testing)
10.5                    TimeTest        (time to start testing in seconds)

! The time step and physical time to start testing.
','type' => 't'}],'attrib' => {'name' => 'TESTTIME'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [{'content' => [],'attrib' => {'value' => '1','name' => 'Rho','default' => 'T'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => '2','name' => 'RhoUx'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => '3','name' => 'RhoUy'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => '4','name' => 'RhoUz'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => '5','name' => 'Bx'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => '6','name' => 'By'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => '7','name' => 'Bz'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => '8','name' => 'e'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => '9','name' => 'p'},'name' => 'option','type' => 'e'}],'attrib' => {'input' => 'select','type' => 'integer','name' => 'iVarTest'},'name' => 'parameter','type' => 'e'},{'content' => '
#TESTVAR
1                       iVarTest

! Index of variable to be tested. Default is rho_="1", ie. density.
','type' => 't'}],'attrib' => {'name' => 'TESTVAR'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [{'content' => [],'attrib' => {'value' => '0','name' => 'all'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => '1','name' => 'x','default' => 'T'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => '2','name' => 'y'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => '3','name' => 'z'},'name' => 'option','type' => 'e'}],'attrib' => {'input' => 'select','type' => 'integer','name' => 'iVarTest'},'name' => 'parameter','type' => 'e'},{'content' => '
#TESTDIM
1                       iDimTest

! Index of dimension/direction to be tested. Default is X dimension.
','type' => 't'}],'attrib' => {'name' => 'TESTDIM'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseStrict','default' => 'T'},'name' => 'parameter','type' => 'e'},{'content' => '
#STRICT
T                       UseStrict

! If true then stop when parameters are incompatible. If false, try to
! correct parameters and continue. Default is true, ie. strict mode
','type' => 't'}],'attrib' => {'name' => 'STRICT'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [{'content' => [],'attrib' => {'value' => '-1','name' => 'errors and warnings only'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => '0','name' => 'start and end of sessions'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => '1','name' => 'normal','default' => 'T'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => '10','name' => 'calls on test processor'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => '100','name' => 'calls on all processors'},'name' => 'option','type' => 'e'}],'attrib' => {'input' => 'select','type' => 'integer','name' => 'lVerbose'},'name' => 'parameter','type' => 'e'},{'content' => '
#VERBOSE
-1                      lVerbose

! Verbosity level controls the amount of output to STDOUT. Default level is 1.
!   lVerbose .le. -1 only warnings and error messages are shown.
!   lVerbose .ge.  0 start and end of sessions is shown.
!   lVerbose .ge.  1 a lot of extra information is given.
!   lVerbose .ge. 10 all calls of set_oktest are shown for the test processor.
!   lVerbose .ge.100 all calls of set_oktest are shown for all processors.
','type' => 't'}],'attrib' => {'name' => 'VERBOSE'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'logical','name' => 'DoDebug','default' => 'F'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'logical','name' => 'DoDebugGhost','default' => 'F'},'name' => 'parameter','type' => 'e'},{'content' => '
#DEBUG
F                       DoDebug         (use it as if(okdebug.and.oktest)...)
F                       DoDebugGhost    (parameter for show_BLK in library.f90)

! Excessive debug output can be controlled by the global okdebug parameter
','type' => 't'}],'attrib' => {'name' => 'DEBUG'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'real','name' => 'CodeVersion','default' => '7.50','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => '
#CODEVERSION
7.50                    CodeVersion

! Checks CodeVersion. Prints a WARNING if it differs from the CodeVersion
! defined in ModMain.f90. Used in newer restart header files. 
! Should be given in PARAM.in when reading old restart files, 
! which do not have version info in the header file.
','type' => 't'}],'attrib' => {'name' => 'CODEVERSION','if' => '$_IsFirstSession'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'length' => '100','type' => 'string','name' => 'NameEquation','default' => 'MHD'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'integer','name' => 'nVar','default' => '8'},'name' => 'parameter','type' => 'e'},{'content' => '
#EQUATION
MHD			NameEquation
8			nVar

! Define the equation name and the number of variables.
! If any of these do not agree with the values determined 
! by the code, BATSRUS stops with an error. Used in restart
! header files and can be given in PARAM.in as a check
! and as a description.
','type' => 't'}],'attrib' => {'name' => 'EQUATION','if' => '$_IsFirstSession'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [{'content' => [],'attrib' => {'value' => '4','name' => 'single precision (4)','default' => '$_nByteReal==4'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => '8','name' => 'double precision (8)','default' => '$_nByteReal==8'},'name' => 'option','type' => 'e'}],'attrib' => {'input' => 'select','type' => 'integer','name' => 'nByteReal'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => '
		nByteReal in file must agree with _nByteReal.
	','type' => 't'}],'attrib' => {'expr' => '$nByteReal==$_nByteReal'},'name' => 'rule','type' => 'e'},{'content' => '

#PRECISION
8                       nByteReal

! Define the number of bytes in a real number. If it does not agree
! with the value determined by the code, BATSRUS stops with an error.
! This is a check, the internal value is calculated in parallel_setup.
! Used in latest restart header files to check binary compatibility.
! May be given in PARAM.in to enforce a certain precision.
','type' => 't'}],'attrib' => {'name' => 'PRECISION','if' => '$_IsFirstSession'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'max' => '$nI','type' => 'integer','name' => 'nI','default' => '$nI','min' => '$nI'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '$nJ','type' => 'integer','name' => 'nJ','default' => '$nJ','min' => '$nJ'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '$nK','type' => 'integer','name' => 'nK','default' => '$nK','min' => '$nK'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '$MaxBlockALL','type' => 'integer','name' => 'MinBlockALL','min' => '1'},'name' => 'parameter','type' => 'e'},{'content' => '

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
','type' => 't'}],'attrib' => {'name' => 'CHECKGRIDSIZE','if' => '$_IsFirstSession'},'name' => 'command','type' => 'e'},{'content' => [{'content' => '
#BLOCKLEVELSRELOADED

This command means that the restart file contains the information about
the minimum and maximum allowed refinement levels for each block.
This command is only used in the restart header file.
','type' => 't'}],'attrib' => {'name' => 'BLOCKLEVELSRELOADED'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseTiming','default' => 'T'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [{'content' => [],'attrib' => {'value' => '-3','name' => 'none'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => '-2','name' => 'final only','default' => 'T'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => '-1','name' => 'end of sessions'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'every X steps','default' => '100','min' => '1'},'name' => 'optioninput','type' => 'e'}],'attrib' => {'input' => 'select','type' => 'integer','name' => 'Frequency'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'integer','name' => 'nDepthTiming','default' => '-1','min' => '-1'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'value' => 'cumm','name' => 'cummulative','default' => '1'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'list'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'tree'},'name' => 'option','type' => 'e'}],'attrib' => {'input' => 'select','type' => 'string','name' => 'TypeTimingReport'},'name' => 'parameter','type' => 'e'}],'attrib' => {'expr' => '$UseTiming'},'name' => 'if','type' => 'e'},{'content' => '
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
','type' => 't'}],'attrib' => {'name' => 'TIMING'},'name' => 'command','type' => 'e'}],'attrib' => {'name' => 'TESTING AND TIMING'},'name' => 'commandgroup','type' => 'e'},{'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! MAIN INITIAL AND BOUNDARY CONDITION PARAMETERS  !!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
','type' => 't'},{'content' => [{'content' => [{'content' => [],'attrib' => {'value' => '1','name' => 'Uniform'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => '2','name' => 'Shock tube'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => '3','name' => 'Heliosphere'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => '5','name' => 'Comet'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => '6','name' => 'Rotation'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => '7','name' => 'Diffusion'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => '11','name' => 'Earth','default' => 'T'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => '12','name' => 'Saturn'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => '13','name' => 'Jupiter'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => '14','name' => 'Venus'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => '21','name' => 'Cylinder'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => '22','name' => 'Sphere'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => '25','name' => 'Arcade'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => '26','name' => 'CME'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => '30','name' => 'Dissipation'},'name' => 'option','type' => 'e'}],'attrib' => {'input' => 'select','type' => 'integer','name' => 'iProblem'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'length' => '20','type' => 'string','name' => 'TypeDissipation','if' => '$iProblem==30'},'name' => 'parameter','type' => 'e'},{'content' => '
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
','type' => 't'}],'attrib' => {'name' => 'PROBLEMTYPE','if' => '$_IsFirstSession','required' => 'T'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [{'content' => [],'attrib' => {'value' => 'GSM','name' => 'GeoSolarMagnetic, GSM','default' => 'T','if' => '$_NameComp eq \'GM\''},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'HGI','name' => 'HelioGraphicInertial, HGI','default' => 'T','if' => '$_NameComp ne \'GM\''},'name' => 'option','type' => 'e'}],'attrib' => {'input' => 'select','type' => 'string','name' => 'TypeCoordSystem'},'name' => 'parameter','type' => 'e'},{'content' => '

#COORDSYSTEM
GSM			TypeCoordSystem

! TypeCoordSystem defines the coordinate system for the component.
! Currently only one coordinate system is available for GM ("GSM")
! and one for IH or SC ("HGI"). In the near future "GSE" should be also
! an option for GM.
!
! Default is component dependent: "GSM" for GM and "HGI" for IH or SC.
','type' => 't'}],'attrib' => {'name' => 'COORDSYSTEM','if' => '$_IsFirstSession'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'length' => '100','type' => 'string','name' => 'NameRestartInDir','default' => 'GM/restartIN'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => '
		Restart input directory $NameRestartInDir must exist
	','type' => 't'}],'attrib' => {'expr' => '-d $NameRestartInDir'},'name' => 'rule','type' => 'e'},{'content' => '

#RESTARTINDIR
GM/restart_n5000	NameRestartInDir

! The NameRestartInDir variable contains the name of the directory
! where restart files are saved relative to the run directory.
! The directory should be inside the subdirectory with the name 
! of the component.
!
! Default value is "GM/restartIN".
','type' => 't'}],'attrib' => {'name' => 'RESTARTINDIR','if' => '$_IsFirstSession'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'logical','name' => 'DoRestartBFace','default' => 'F'},'name' => 'parameter','type' => 'e'},{'content' => ' 

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

','type' => 't'}],'attrib' => {'name' => 'NEWRESTART','if' => '$_IsFirstSession'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'integer','name' => 'nRootBlockX','default' => '2','min' => '1'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'integer','name' => 'nRootBlockY','default' => '1','min' => '1'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'integer','name' => 'nRootBlockZ','default' => '1','min' => '1'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'xMin','default' => '-192.0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'xMax','default' => '  64.0','min' => '$xMin'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'yMin','default' => ' -64.0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'yMax','default' => '  64.0','min' => '$yMin'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'zMin','default' => ' -64.0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'zMax','default' => '  64.0','min' => '$zMin'},'name' => 'parameter','type' => 'e'},{'content' => '
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
','type' => 't'}],'attrib' => {'name' => 'GRID','if' => '$_IsFirstSession','required' => 'T'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [{'content' => [{'content' => [],'attrib' => {'name' => 'coupled'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'fixed/inflow','default' => '$Side ne \'TypeBcEast\''},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'float/outflow','default' => '$Side eq \'TypeBcEast\''},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'heliofloat'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'reflect'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'periodic'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'vary'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'shear'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'linetied'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'raeder'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'arcadetop'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'arcadebot'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'arcadebotcont'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'user'},'name' => 'option','type' => 'e'}],'attrib' => {'input' => 'select','type' => 'string','name' => '$Side'},'name' => 'parameter','type' => 'e'}],'attrib' => {'values' => 'TypeBcEast,TypeBcWest,TypeBcSouth,TypeBcNorth,TypeBcBot,TypeBcTop','name' => 'Side'},'name' => 'foreach','type' => 'e'},{'content' => [{'content' => '
		East and west BCs must be both periodic or neither
	','type' => 't'}],'attrib' => {'expr' => 'not($TypeBcEast eq \'periodic\' xor $TypeBcWest eq \'periodic\')'},'name' => 'rule','type' => 'e'},{'content' => [{'content' => '
		South and North BCs must be both periodic or neither
	','type' => 't'}],'attrib' => {'expr' => 'not($TypeBcSouth eq \'periodic\' xor $TypeBcNorth eq \'periodic\')'},'name' => 'rule','type' => 'e'},{'content' => [{'content' => '
		Bottom and top BCs must be both periodic or neither
	','type' => 't'}],'attrib' => {'expr' => 'not($TypeBcBot eq \'periodic\' xor $TypeBcTop eq \'periodic\')'},'name' => 'rule','type' => 'e'},{'content' => '
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
','type' => 't'}],'attrib' => {'name' => 'OUTERBOUNDARY'},'name' => 'command','type' => 'e'},{'content' => [{'content' => '
! Inner boundary types for body 1 and body 2
	','type' => 't'},{'content' => [{'content' => [],'attrib' => {'name' => 'reflect'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'float'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'fixed'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'ionosphere','default' => 'T'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'ionosphereB0/ionosphereb0','name' => 'ionosphereB0'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'ionospherefloat'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'coronatoih'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'user'},'name' => 'option','type' => 'e'}],'attrib' => {'input' => 'select','type' => 'string','name' => 'TypeBcInner'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [{'content' => [],'attrib' => {'name' => 'reflect','default' => 'T'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'float'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'fixed'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'ionosphere'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'ionosphereB0'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'ionospherefloat'},'name' => 'option','type' => 'e'}],'attrib' => {'input' => 'select','type' => 'string','name' => 'TypeBcBody2'},'name' => 'parameter','type' => 'e'}],'attrib' => {'expr' => '$UseBody2'},'name' => 'if','type' => 'e'},{'content' => [{'content' => ' 
	For the second body COROTATION AND AN IONOSPHERIC BOUNDARY DO NOT WORK.
	','type' => 't'}],'attrib' => {'expr' => 'not($TypeBcBody2 =~ /ionosphere/)'},'name' => 'rule','type' => 'e'},{'content' => '

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
','type' => 't'}],'attrib' => {'name' => 'INNERBOUNDARY'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseExtraBoundary','default' => 'F'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'string','name' => 'TypeBcExtra'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'logical','name' => 'DoFixExtraboundary','default' => 'F'},'name' => 'parameter','type' => 'e'}],'attrib' => {'expr' => '$UseExtraBoundary'},'name' => 'if','type' => 'e'}],'attrib' => {'name' => 'EXTRABOUNDARY'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'max' => '6','type' => 'integer','name' => 'MaxBoundary','default' => '0','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'logical','name' => 'DoFixOuterBoundary','default' => 'F'},'name' => 'parameter','type' => 'e'}],'attrib' => {'expr' => '$MaxBoundary >= 1'},'name' => 'if','type' => 'e'},{'content' => '
#FACEOUTERBC
0              MaxBoundary            
F              DoFixOuterBoundary)    !read only for MaxBoundary>=East_(=1).
! If MaxBoundary is East_(=1) or more then the outer boundaries with
! the number of boundary being between East_ and MaxBoundary
! are treated using set_BCs.f90 subroutines instead of set_outerBCs.f90 
! if DoFixOuterBoundary is .true., there is no resolution
! change along the outer boundaries with the number of
! of boundary being between East_ and MaxBoundary
','type' => 't'}],'attrib' => {'name' => 'FACEOUTERBC'},'name' => 'command','type' => 'e'}],'attrib' => {'name' => 'INITIAL AND BOUNDARY CONDITIONS'},'name' => 'commandgroup','type' => 'e'},{'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! INITIAL TIME AND STEP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
','type' => 't'},{'content' => [{'content' => [],'attrib' => {'type' => 'integer','name' => 'iYear','default' => '2000'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '12','type' => 'integer','name' => 'iMonth','default' => '3','min' => '1'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '31','type' => 'integer','name' => 'iDay','default' => '21','min' => '1'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '23','type' => 'integer','name' => 'iHour','default' => '0','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '59','type' => 'integer','name' => 'iMinute','default' => '0','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '59','type' => 'integer','name' => 'iSecond','default' => '0','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => '
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
','type' => 't'}],'attrib' => {'alias' => 'SETREALTIME','name' => 'STARTTIME','if' => '$_IsFirstSession'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'real','name' => 'tSimulation','default' => '0.0','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => '

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
','type' => 't'}],'attrib' => {'name' => 'TIMESIMULATION','if' => '$_IsFirstSession'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'integer','name' => 'nStep','default' => '0','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => '

#NSTEP
100			nStep

! Set nStep for the component. Typically used in the restart.H header file.
! Generally it is not inserted in a PARAM.in file by the user.
!
! The default is nStep=0 as the starting time step with no restart.
','type' => 't'}],'attrib' => {'name' => 'NSTEP','if' => '$_IsFirstSession'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'integer','name' => 'nPrevious','default' => '-1','min' => '-1'},'name' => 'parameter','type' => 'e'},{'content' => '

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
','type' => 't'}],'attrib' => {'name' => 'NPREVIOUS','if' => '$_IsFirstSession'},'name' => 'command','type' => 'e'}],'attrib' => {'name' => 'INITIAL TIME'},'name' => 'commandgroup','type' => 'e'},{'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  TIME INTEGRATION PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
','type' => 't'},{'content' => [{'content' => [{'content' => [],'attrib' => {'value' => '1','default' => 'T'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => '2'},'name' => 'option','type' => 'e'}],'attrib' => {'input' => 'select','type' => 'integer','name' => 'nStage'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '1','type' => 'real','name' => 'CflExpl','default' => '0.8','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => '

#TIMESTEPPING
2                       nStage
0.80                    CflExpl

! Parameters for explicit time integration.
! Default is 1 stage and CflExpl=0.8
','type' => 't'}],'attrib' => {'name' => 'TIMESTEPPING'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseDtFixed','default' => 'F'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'DtFixedDim','default' => '1.0','min' => '0','if' => '$UseDtFixed'},'name' => 'parameter','type' => 'e'},{'content' => '
#FIXEDTIMESTEP
T                       UseDtFixed
10.                     DtFixedDim [sec] (read if UseDtFixed is true)

! Default is UseDtFixed=.false. Effective only if DoTimeAccurate is true.
! If UseDtFixed is true, the time step is fixed to DtFixedDim.
!
! This is useful for debugging explicit schemes.

! The real application is, however, for implicit and partially
! implicit/local schemes.

','type' => 't'}],'attrib' => {'name' => 'FIXEDTIMESTEP'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'logical','name' => 'UsePartLocal','default' => 'F'},'name' => 'parameter','type' => 'e'},{'content' => '
#PARTLOCAL
T               UsePartLocal

! Default is UsePartLocal=.false. If UsePartLocal is true and the
! run is time accurate, then the blocks selected as "implicit"
! by the criteria defined in #STEPPINGCRITERIA are not used to
! calculate the time step, and all cells are advanced with the
! smaller of the stable and the global time steps.
','type' => 't'}],'attrib' => {'name' => 'PARTLOCAL'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'logical','name' => 'UsePointImplicit','default' => 'F'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'logical','name' => 'UsePartImplicit','default' => 'F'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseFullImplicit','default' => 'F'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'CflImpl','default' => '100','min' => '0','if' => '$UsePartImplicit or $UseFullImplicit'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => '
		At most one of these logicals can be true!
	','type' => 't'}],'attrib' => {'expr' => '$UsePointImplicit + $UsePartImplicit + $UseFullImplicit <= 1'},'name' => 'rule','type' => 'e'},{'content' => '

#IMPLICIT
F               UsePointImplicit   
F               UsePartImplicit
F               UseFullImplicit
100.0           CflImpl (read if UsePartImplicit or UseFullImplicit is true)

! Default is false for all logicals. Only one of them can be set to true!
! The CFL number is used in the implicit blocks of the fully or partially
! implicit schemes. Ignored if UseDtFixed is true.
','type' => 't'}],'attrib' => {'name' => 'IMPLICIT'},'name' => 'command','type' => 'e'}],'attrib' => {'name' => 'TIME INTEGRATION'},'name' => 'commandgroup','type' => 'e'},{'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! PARAMETERS FOR FULL AND PART IMPLICIT TIME INTEGRATION !!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
','type' => 't'},{'content' => [{'content' => [{'content' => [],'attrib' => {'value' => 'dt','name' => 'Time step','default' => 'T'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'r/R','name' => 'Radial distance'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'test','name' => 'Test block'},'name' => 'option','type' => 'e'}],'attrib' => {'input' => 'select','type' => 'string','name' => 'TypeImplCrit'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'rImplicit','min' => '0','if' => '$TypeImplCrit eq \'R\''},'name' => 'parameter','type' => 'e'},{'content' => '

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
','type' => 't'}],'attrib' => {'alias' => 'STEPPINGCRITERIA','name' => 'IMPLICITCRITERIA'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'max' => '1','type' => 'real','name' => 'ImplCoeff','default' => '1','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseBdf2','default' => 'T'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseSourceImpl','default' => 'F'},'name' => 'parameter','type' => 'e'},{'content' => '
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
','type' => 't'}],'attrib' => {'name' => 'IMPLSTEP'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [{'content' => [],'attrib' => {'name' => '1','default' => 'T'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => '2'},'name' => 'option','type' => 'e'}],'attrib' => {'input' => 'select','type' => 'integer','name' => 'nOrderImpl'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'value' => 'Roe/roe/1','name' => 'Roe'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'Rusanov/rusanov/2/TVDLF','name' => 'Rusanov','default' => 'T'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'Linde/linde/3/HLLEL','name' => 'Linde'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'Sokolov/sokolov/4/AW','name' => 'Sokolov'},'name' => 'option','type' => 'e'}],'attrib' => {'input' => 'select','type' => 'string','name' => 'TypeFluxImpl'},'name' => 'parameter','type' => 'e'},{'content' => '
#IMPLSCHEME
1               nOrderImpl
Rusanov         TypeFluxImpl

! Default values are shown, ie. first order Rusanov scheme.
! This defines the scheme used in the implicit part.
','type' => 't'}],'attrib' => {'name' => 'IMPLSCHEME'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseConservativeImplicit','default' => 'F'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseNewton','default' => 'F'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseNewMatrix','default' => 'T'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'integer','name' => 'MaxIterNewton','default' => '10','min' => '1'},'name' => 'parameter','type' => 'e'}],'attrib' => {'expr' => '$UseNewton'},'name' => 'if','type' => 'e'},{'content' => '
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
','type' => 't'}],'attrib' => {'name' => 'NEWTON'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [{'content' => [],'attrib' => {'value' => 'prec','name' => 'Preconditioned','default' => 'T'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'free','name' => 'No preconditioning'},'name' => 'option','type' => 'e'}],'attrib' => {'input' => 'select','type' => 'string','name' => 'TypeJacobian'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '1.e-5','type' => 'real','name' => 'JacobianEps','default' => '$doublePrecision ? 1.e-12 : 1.e-6','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => '
#JACOBIAN
prec            TypeJacobian (prec, free)
1.E-12          JacobianEps

! The Jacobian matrix is always calculated with a matrix free approach,
! however it can be preconditioned  (\'prec\'), or not (\'free\')
! Default value is TypeJacobian=\'prec\'.
! JacobianEps contains the machine round off error for numerical derivatives.
! The default value is 1.E-12 for 8 byte reals and 1.E-6 for 4 byte reals.
','type' => 't'}],'attrib' => {'name' => 'JACOBIAN'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [{'content' => [],'attrib' => {'name' => 'left'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'symmetric','default' => 'T'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'right'},'name' => 'option','type' => 'e'}],'attrib' => {'input' => 'select','type' => 'string','name' => 'TypePrecondSide'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'name' => 'MBILU','default' => 'T'},'name' => 'option','type' => 'e'}],'attrib' => {'input' => 'select','type' => 'string','name' => 'TypePrecond'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '1','type' => 'real','name' => 'GustafssonPar','default' => '0.5','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => '
#PRECONDITIONER
symmetric       TypePrecondSide (left, symmetric, right)
MBILU           TypePrecond (MBILU)
0.5             GustafssonPar (0. no modification, 1. full modification)

! Default parameters are shown. Right preconditioning does not affect
! the normalization of the residual. The Gustafsson parameter determines
! how much the MBILU preconditioner is modified. The default 0.5 value
! means a relaxed modification.
','type' => 't'}],'attrib' => {'name' => 'PRECONDITIONER'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [{'content' => [],'attrib' => {'name' => 'gmres','default' => 'T'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'bicgstab'},'name' => 'option','type' => 'e'}],'attrib' => {'input' => 'select','type' => 'string','name' => 'TypeKrylov'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'value' => 'nul','name' => '0','default' => 'T'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'old','name' => 'previous'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'explicit'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'explicit','name' => 'scaled explicit'},'name' => 'option','type' => 'e'}],'attrib' => {'input' => 'select','type' => 'string','name' => 'TypeInitKrylov'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '0.1','type' => 'real','name' => 'ErrorMaxKrylov','default' => '0.001','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'integer','name' => 'MaxMatvecKrylov','default' => '100','min' => '1'},'name' => 'parameter','type' => 'e'},{'content' => '
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
','type' => 't'}],'attrib' => {'name' => 'KRYLOV'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'integer','name' => 'nKrylovVector','default' => 'MaxMatvecKrylov','min' => '1'},'name' => 'parameter','type' => 'e'},{'content' => '
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
','type' => 't'}],'attrib' => {'name' => 'KRYLOVSIZE'},'name' => 'command','type' => 'e'}],'attrib' => {'name' => 'IMPLICIT PARAMETERS'},'name' => 'commandgroup','type' => 'e'},{'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! STOPPING CRITERIA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

The commands in this group only work in stand alone mode.

','type' => 't'},{'content' => [{'content' => [],'attrib' => {'type' => 'integer','name' => 'MaxIteration','default' => '-1','min' => '-1'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'tSimulationMax','default' => '-1','min' => '-1'},'name' => 'parameter','type' => 'e'},{'content' => '

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
','type' => 't'}],'attrib' => {'name' => 'STOP','if' => '$_IsStandAlone','required' => '$_IsStandAlone'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'logical','name' => 'DoCheckStopFile','default' => 'T'},'name' => 'parameter','type' => 'e'},{'content' => '

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
','type' => 't'}],'attrib' => {'name' => 'CHECKSTOPFILE','if' => '$_IsStandAlone'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'real','name' => 'CpuTimeMax','default' => '-1','min' => '-1'},'name' => 'parameter','type' => 'e'},{'content' => '

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
','type' => 't'}],'attrib' => {'name' => 'CPUTIMEMAX','if' => '$_IsStandAlone'},'name' => 'command','type' => 'e'}],'attrib' => {'name' => 'STOPPING CRITERIA'},'name' => 'commandgroup','type' => 'e'},{'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  OUTPUT PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
','type' => 't'},{'content' => [{'content' => [],'attrib' => {'length' => '100','type' => 'string','name' => 'NameRestartOutDir','default' => 'GM/restartOUT'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => '
		Restart output directory $NameRestartOutDir must exist
	','type' => 't'}],'attrib' => {'expr' => '-d $NameRestartOutDir'},'name' => 'rule','type' => 'e'},{'content' => '

#RESTARTOUTDIR
GM/restart_n5000	NameRestartOutDir

! The NameRestartOutDir variable contains the name of the directory
! where restart files are saved relative to the run directory.
! The directory should be inside the subdirectory with the name 
! of the component.
!
! Default value is "GM/restartOUT".
','type' => 't'}],'attrib' => {'name' => 'RESTARTOUTDIR'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'logical','name' => 'DoSaveRestart','default' => 'T'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'integer','name' => 'DnSaveRestart','default' => '-1','min' => '-1'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'DtSaveRestart','default' => '-1','min' => '-1'},'name' => 'parameter','type' => 'e'}],'attrib' => {'expr' => '$DoSaveRestart'},'name' => 'if','type' => 'e'},{'content' => '
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
','type' => 't'}],'attrib' => {'name' => 'SAVERESTART'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'length' => '100','type' => 'string','name' => 'NamePlotDir','default' => 'GM/IO2'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => '
		Plot directory $NamePlotDir must exist
	','type' => 't'}],'attrib' => {'expr' => '-d $NamePlotDir'},'name' => 'rule','type' => 'e'},{'content' => '

The NamePlotDir variable contains the name of the directory
where plot files and logfiles are saved relative to the run directory.
The directory should be inside the subdirectory with the name
of the component.

Default value is "GM/IO2".
','type' => 't'}],'attrib' => {'name' => 'PLOTDIR'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'logical','name' => 'DoSaveLogfile','default' => 'F'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [{'content' => [{'content' => [],'attrib' => {'value' => 'MHD','name' => 'MHD vars. dimensional','default' => 'T'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'FLX','name' => 'Flux vars. dimensional'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'RAW','name' => 'Raw vars. dimensional'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'VAR','name' => 'Set vars. dimensional'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'mhd','name' => 'MHD vars. scaled','default' => 'T'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'flx','name' => 'Flux vars. scaled'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'raw','name' => 'Raw vars. scaled'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'var','name' => 'Set vars. scaled'},'name' => 'option','type' => 'e'}],'attrib' => {'input' => 'select','type' => 'string','name' => 'TypeLogVar','required' => 'T'},'name' => 'part','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'exclusive' => 'T','name' => 'none'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'step'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'date'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'time'},'name' => 'option','type' => 'e'}],'attrib' => {'input' => 'select','type' => 'string','name' => 'TypeTime','multiple' => 'T','required' => 'F'},'name' => 'part','type' => 'e'}],'attrib' => {'max' => '4','type' => 'strings','name' => 'StringLog','min' => '1'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'integer','name' => 'DnSaveLogfile','default' => '1','min' => '-1'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'DtSaveLogfile','default' => '-1','min' => '-1'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'length' => '100','type' => 'string','name' => 'NameLogVars','if' => '$TypeLogVar =~ /var/i'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'real','name' => 'LogRadii','multiple' => 'T','min' => '$rBody'},'name' => 'part','type' => 'e'}],'attrib' => {'length' => '100','max' => '10','type' => 'strings','name' => 'StringLogRadii','min' => '1','if' => '($TypeLogVar=~/flx/i or $NameLogVars=~/flx/i)'},'name' => 'parameter','type' => 'e'}],'attrib' => {'expr' => '$DoSaveLogfile'},'name' => 'if','type' => 'e'},{'content' => '
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
','type' => 't'}],'attrib' => {'name' => 'SAVELOGFILE'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'integer','name' => 'nSatellite','default' => '0','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [{'content' => [{'content' => [],'attrib' => {'value' => 'MHD','name' => 'MHD vars. dimensional','default' => 'T'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'FUL','name' => 'All vars. dimensional'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'VAR','name' => 'Set vars. dimensional'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'mhd','name' => 'MHD vars. scaled'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'ful','name' => 'All vars. scaled'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'var','name' => 'Set vars. scaled'},'name' => 'option','type' => 'e'}],'attrib' => {'input' => 'select','type' => 'string','name' => 'TypeSatelliteVar','required' => 'T'},'name' => 'part','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'name' => 'file','default' => 'T'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'eqn','name' => 'equation'},'name' => 'option','type' => 'e'}],'attrib' => {'input' => 'select','type' => 'string','name' => 'TypeTrajectory','required' => 'F'},'name' => 'part','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'exclusive' => 'T','name' => 'none'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'step'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'date'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'time'},'name' => 'option','type' => 'e'}],'attrib' => {'input' => 'select','type' => 'string','name' => 'TypeTime','multiple' => 'T','required' => 'F'},'name' => 'part','type' => 'e'}],'attrib' => {'max' => '5','type' => 'strings','name' => 'StringSatellite','min' => '1'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'integer','name' => 'DnOutput','default' => '1','min' => '-1'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'DtOutput','default' => '-1','min' => '-1'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'length' => '100','type' => 'string','name' => 'NameTrajectoryFile'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => '
			Trajectory file $NameTrajectoryFile must exist
		','type' => 't'}],'attrib' => {'expr' => '-f $NameTrajectoryFile'},'name' => 'rule','type' => 'e'},{'content' => [],'attrib' => {'length' => '100','type' => 'string','name' => 'NameSatelliteVars','if' => '$TypeSatelliteVar =~ /\\bvar\\b/i'},'name' => 'parameter','type' => 'e'}],'attrib' => {'to' => '$nSatellite','from' => '1'},'name' => 'for','type' => 'e'},{'content' => '
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
','type' => 't'}],'attrib' => {'name' => 'SATELLITE','if' => '$_IsFirstSession'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'max' => '100','type' => 'integer','name' => 'nPlotFile','default' => '0','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [{'content' => [{'content' => [],'attrib' => {'value' => 'tec','name' => 'TECPLOT'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'idl','name' => 'IDL'},'name' => 'option','type' => 'e'}],'attrib' => {'input' => 'select','type' => 'string','name' => 'plotform','required' => 'T'},'name' => 'part','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'value' => '3d/3d_','name' => '3D'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'x=0'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'y=0','default' => 'T'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'z=0'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'sph'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'los'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'lin'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'cut','if' => '$plotform =~ /\\bidl\\b/'},'name' => 'option','type' => 'e'}],'attrib' => {'input' => 'select','type' => 'string','name' => 'plotarea','required' => 'T'},'name' => 'part','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'value' => 'MHD','name' => 'MHD vars. dimensional'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'FUL','name' => 'All vars. dimensional'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'RAW','name' => 'Raw vars. dimensional'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'RAY','name' => 'Ray tracing vars. dim.'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'FLX','name' => 'Flux vars. dimensional'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'SOL','name' => 'Solar vars. dimensional'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'POS','name' => 'Position vars. dimensional','if' => '$plotarea eq \'lin\''},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'VAR','name' => 'Select dimensional vars.'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'mhd','name' => 'MHD vars. scaled'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'ful','name' => 'All vars. scaled'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'raw','name' => 'Raw vars. scaled'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'ray','name' => 'Ray tracing vars. scaled'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'flx','name' => 'Flux vars. scaled'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'sol','name' => 'Solar vars. scaled'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'pos','name' => 'Position vars. scaled','if' => '$plotarea eq \'lin\''},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'var','name' => 'Select scaled vars.'},'name' => 'option','type' => 'e'}],'attrib' => {'input' => 'select','type' => 'string','name' => 'plotvar','required' => 'T'},'name' => 'part','type' => 'e'}],'attrib' => {'max' => '3','type' => 'strings','name' => 'StringPlot','min' => '3'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'integer','name' => 'DnSavePlot','min' => '-1'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'DtSavePlot','min' => '-1'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'real','name' => 'xMinCut'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'xMaxCut','min' => '$xMinCut'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'yMinCut'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'yMaxCut','min' => '$yMinCut'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'zMinCut'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'zMaxCut','min' => '$zMinCut'},'name' => 'parameter','type' => 'e'}],'attrib' => {'expr' => '$plotarea =~ /\\bcut\\b/'},'name' => 'if','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'Radius','default' => '10','min' => '0','if' => '$plotarea =~ /\\bsph\\b/'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'real','name' => 'LosVectorX','default' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'LosVectorY','default' => '0.0001'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'LosVectorZ','default' => '1'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'xSizeImage','default' => '20','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'ySizeImage','default' => '20','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'xOffset','default' => '10'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'yOffset','default' => '10'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'rOccult','default' => '2.5','min' => '1'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '1','type' => 'real','name' => 'MuLimbDarkening','default' => '0.5','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'integer','name' => 'nPixX','default' => '200','min' => '2'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'integer','name' => 'nPixY','default' => '200','min' => '2'},'name' => 'parameter','type' => 'e'}],'attrib' => {'expr' => '$plotarea =~ /\\blos\\b/'},'name' => 'if','type' => 'e'},{'content' => [{'content' => [{'content' => [],'attrib' => {'value' => 'A','name' => 'Advected B'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'B','default' => 'T'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'U'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'J'},'name' => 'option','type' => 'e'}],'attrib' => {'input' => 'select','type' => 'string','name' => 'NameLine'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'logical','name' => 'IsSingleLine','default' => 'F'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '20','type' => 'integer','name' => 'nLine','default' => '1','min' => '1'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'real','name' => 'xStartLine'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'yStartLine'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'zStartLine'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'logical','name' => 'IsParallel'},'name' => 'parameter','type' => 'e'}],'attrib' => {'to' => '$nLine','from' => '1'},'name' => 'for','type' => 'e'}],'attrib' => {'expr' => '$plotarea =~ /\\blin\\b/'},'name' => 'if','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'DxSavePlot','default' => '-1.0','min' => '-1.0','if' => '($plotform=~/\\bidl\\b/ and $plotarea!~/\\b(sph|los|lin)\\b/)'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'length' => '100','type' => 'string','name' => 'NameVars'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'length' => '100','type' => 'string','name' => 'NamePars'},'name' => 'parameter','type' => 'e'}],'attrib' => {'expr' => '$plotvar =~ /\\bvar\\b/i'},'name' => 'if','type' => 'e'}],'attrib' => {'to' => '$nPlotFile','from' => '1'},'name' => 'for','type' => 'e'},{'content' => '
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

','type' => 't'}],'attrib' => {'name' => 'SAVEPLOT'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'logical','name' => 'DoSaveBinary','default' => 'T'},'name' => 'parameter','type' => 'e'},{'content' => '
#SAVEBINARY
T			DoSaveBinary   used only for \'idl\' plot file

! Default is .true. Saves unformatted IO2/*.idl files if true. 
! This is the recommended method, because it is fast and accurate.
! The only advantage of saving IO2/*.idl in formatted text files is
! that it can be processed on another machine or with a different 
! (lower) precision. For example PostIDL.exe may be compiled with 
! single precision to make IO2/*.out files smaller, while BATSRUS.exe is 
! compiled in double precision, to make results more accurate.
','type' => 't'}],'attrib' => {'name' => 'SAVEBINARY'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'logical','name' => 'DoSavePlotsAmr','default' => 'F'},'name' => 'parameter','type' => 'e'},{'content' => '
#SAVEPLOTSAMR
F			DoSavePlotsAmr

! Save plots before each AMR. Default is DoSavePlotsAMR=.false.
','type' => 't'}],'attrib' => {'name' => 'SAVEPLOTSAMR'},'name' => 'command','type' => 'e'}],'attrib' => {'name' => 'OUTPUT PARAMETERS'},'name' => 'commandgroup','type' => 'e'},{'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  AMR PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
','type' => 't'},{'content' => [{'content' => [{'content' => [],'attrib' => {'name' => 'default','default' => '1'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'all'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'none'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => '3Dbodyfocus'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'spherefocus'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'magnetosphere'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'points'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'coupledhelio'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'helio_init'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'helio_z=4'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'all_then_focus'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'cme'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'points'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'mag_new'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'magnetosphere'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'magneto_fine'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'magneto12'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'magnetosaturn'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'magnetojupiter'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'paleo'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'comet'},'name' => 'option','type' => 'e'}],'attrib' => {'input' => 'select','type' => 'string','name' => 'InitialRefineType'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'integer','name' => 'InitialRefineLevel','default' => '4','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => '
#AMRINIT
default			TypeRefineInit
4			nRefineLevelInit

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
','type' => 't'}],'attrib' => {'name' => 'AMRINIT','if' => '$_IsFirstSession'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'integer','name' => 'nRefineLevelIC','default' => '0','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => '
#AMRINITPHYSICS
3			nRefineLevelIC

! Defines number of physics (initial condition) based AMR-s AFTER the 
! geometry based initial AMR-s defined by #AMRINIT were done.
! Only useful if the initial condition has a non-trivial analytic form.
','type' => 't'}],'attrib' => {'name' => 'AMRINITPHYSICS','if' => '$_IsFirstSession'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'length' => '100','type' => 'string','name' => 'Resolution','if' => '$_command =~ /RESOLUTION/'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'integer','name' => 'nLevel','min' => '0','if' => '$_command =~ /LEVEL/'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'name' => 'all','default' => 'T'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'box'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'brick'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'sphere'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'sphere0'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'shell'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'shell0'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'cylinderx'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'cylinderx0'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'cylindery'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'cylindery0'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'cylinderz'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'cylinderz0'},'name' => 'option','type' => 'e'}],'attrib' => {'input' => 'select','type' => 'string','name' => 'NameArea'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'real','name' => 'xMinBox'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'yMinBox'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'zMinBox'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'xMaxBox'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'yMaxBox'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'zMaxBox'},'name' => 'parameter','type' => 'e'}],'attrib' => {'expr' => '$NameArea =~ /box/'},'name' => 'if','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'real','name' => 'xCenterBrick'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'yCenterBrick'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'zCenterBrick'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'xSizeBrick','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'ySizeBrick','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'zSizeBrick','min' => '0'},'name' => 'parameter','type' => 'e'}],'attrib' => {'expr' => '$NameArea =~ /brick/i'},'name' => 'if','type' => 'e'},{'content' => [{'content' => [{'content' => [],'attrib' => {'type' => 'real','name' => 'xCenterShell'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'yCenterShell'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'zCenterShell'},'name' => 'parameter','type' => 'e'}],'attrib' => {'expr' => '$NameArea !~ /0/i'},'name' => 'if','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'r1Shell','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'r2Shell','min' => '0'},'name' => 'parameter','type' => 'e'}],'attrib' => {'expr' => '$NameArea =~ /shell/i'},'name' => 'if','type' => 'e'},{'content' => [{'content' => [{'content' => [],'attrib' => {'type' => 'real','name' => 'xCenterSphere'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'yCenterSphere'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'zCenterSphere'},'name' => 'parameter','type' => 'e'}],'attrib' => {'expr' => '$NameArea !~ /sphere0/i'},'name' => 'if','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'rSphere','min' => '0'},'name' => 'parameter','type' => 'e'}],'attrib' => {'expr' => '$NameArea =~ /sphere/i'},'name' => 'if','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'real','name' => 'x1Cylinder'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'yCylinder','if' => '$NameArea!~/0/'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'zCylinder','if' => '$NameArea!~/0/'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'x2Cylinder'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'rCylinder','min' => '0'},'name' => 'parameter','type' => 'e'}],'attrib' => {'expr' => '$NameArea =~ /cylinderx/i'},'name' => 'if','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'real','name' => 'xCylinder','if' => '$NameArea!~/0/'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'y1Cylinder'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'zCylinder','if' => '$NameArea!~/0/'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'y2Cylinder'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'rCylinder','min' => '0'},'name' => 'parameter','type' => 'e'}],'attrib' => {'expr' => '$NameArea =~ /cylindery/i'},'name' => 'if','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'real','name' => 'xCylinder','if' => '$NameArea!~/0/'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'yCylinder','if' => '$NameArea!~/0/'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'z1Cylinder'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'z2Cylinder'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'rCylinder','min' => '0'},'name' => 'parameter','type' => 'e'}],'attrib' => {'expr' => '$NameArea =~ /cylinderz/i'},'name' => 'if','type' => 'e'},{'content' => '
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
-48.0			xCenterBrick
  0.0			yCenterBrick
 -8.0			zCenterBrick
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
cylinderx0		NameArea
-60.0			x1Cylinder
-10.0			x2Cylinder
 20.0			rCylinder

#GRIDRESOLUTION
1/8			Resolution
cylindery		NameArea
-60.0			xCylinder
-20.0			y1Cylinder
  0.0			zCylinder
 20.0			y2Cylinder
 10.0			rCylinder

The #GRIDRESOLUTION and #GRIDLEVEL commands allow to set the grid resolution
or refinement level, respectively, in a given area. The Resolution parameter
refers to the size of the cell in the X direction (Dx) and it can be given
as a real number (e.g. 0.125) as the inverse of the real number (e.g. 1/8).
The nLevel parameter is an integer with level 0 meaning no refinement relative
to the root block, while level N is a refinement by 2 to the power N.

The next parameter NameArea defines the shape of the area. 
All computational blocks that intersect the area and have a coarser
resolution than the resolution set for the area are selected for refinement.
There are the following basic shapes: 
\'all\', \'box\', \'brick\', \'sphere\', \'shell\', \'cylinderx\', \'cylindery\' 
and \'cylinderz\'.
The area \'all\' refers to the whole computational domain, and it can be
used to set the overall minimum resolution. The area \'box\' is a box
aligned with the X, Y and Z axes, and it is given with the coordinates
of two diagonally opposite corners. The area \'brick\' has the same shape
as \'box\', but it is defines with the center of the brick and the 
size of the brick. The area \'sphere\' is a sphere around an arbitrary point,
which is defined with the center point and the radius of the sphere.
The area \'shell\' consists of the volume between two concentric spherical
surfaces, which is given with the center point and the two radii.
The area \'cylinderx\' is a cylinder with an axis parallel with the X axis,
and it is given with the 3 coordinates of one end of the axis of the cylinder,
the X coordinate of the other end of the axis and the radius of the cylinder.
The areas \'cylindery\' and \'cylinderz\' are cylinders parallel with the 
Y and Z axes, respectively, and are defined analogously as \'cylinderx\'.
The areas \'sphere0\' and \'shell0\' are the same as areas \'sphere\' and \'shell\'
except that their center is in the origin, so only the radii need to be
given. The areas \'cylinderx0\', \'cylindery0\' and \'cylinderz0\' are 
cylinders with their axes lying on the X, Y and Z axes respectively,
so only two coordinates and a radius are required.


','type' => 't'}],'attrib' => {'alias' => 'GRIDLEVEL','name' => 'GRIDRESOLUTION','multiple' => 'T'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'integer','name' => 'MinBlockLevel','default' => '0','min' => '-1'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'integer','name' => 'MaxBlockLevel','default' => '99','min' => '-1'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'logical','name' => 'DoFixBodyLevel','default' => 'F'},'name' => 'parameter','type' => 'e'},{'content' => '
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
','type' => 't'}],'attrib' => {'name' => 'AMRLEVELS'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'real','name' => 'DxCellMin','default' => '0','min' => '-1'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'DxCellMax','default' => '99999','min' => '-1'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'logical','name' => 'DoFixBodyLevel','default' => 'F'},'name' => 'parameter','type' => 'e'},{'content' => '
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
','type' => 't'}],'attrib' => {'name' => 'AMRRESOLUTION'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'integer','name' => 'DnRefine','default' => '-1','min' => '-1'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'logical','name' => 'DoAutoRefine','default' => 'F'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'max' => '100','type' => 'real','name' => 'PercentCoarsen','default' => '20','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '100','type' => 'real','name' => 'PercentRefine','default' => '20','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'integer','name' => 'MaxTotalBlocks','default' => '99999','min' => '1'},'name' => 'parameter','type' => 'e'}],'attrib' => {'expr' => '$DoAutoRefine'},'name' => 'if','type' => 'e'}],'attrib' => {'expr' => '$DnRefine>0'},'name' => 'if','type' => 'e'},{'content' => '
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
','type' => 't'}],'attrib' => {'name' => 'AMR'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [{'content' => [],'attrib' => {'name' => '1'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => '2'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => '3','default' => '1'},'name' => 'option','type' => 'e'}],'attrib' => {'input' => 'select','type' => 'integer','name' => 'nRefineCrit'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [{'content' => [],'attrib' => {'value' => 'gradt/gradT','name' => 'grad T'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'gradp/gradP','name' => 'grad P'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'gradlogrho','name' => 'grad log(Rho)'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'gradlogP/gradlogp','name' => 'grad log(p)'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'gradE','name' => 'grad E'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'curlV/curlv/curlU/curlu','name' => 'curl U'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'curlB/curlb','name' => 'curl B'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'divU/divu/divV/divv','name' => 'div U'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'divb/divB','name' => 'divB'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'Valfven/vAlfven/valfven','name' => 'vAlfven'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'heliobeta','name' => 'heliospheric beta'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'flux'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'heliocurrentsheet','name' => 'heliospheric current sheet'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'rcurrents/Rcurrents','name' => 'rCurrents'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'transient/Transient','name' => 'Transient'},'name' => 'option','type' => 'e'}],'attrib' => {'input' => 'select','type' => 'string','name' => 'TypeRefine'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [{'content' => [],'attrib' => {'value' => 'p_dot/P_dot','name' => 'P_dot'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 't_dot/T_dot','name' => 'T_dot'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'rho_dot/Rho_dot','name' => 'Rho_dot','default' => 'T'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'RhoU_dot/rhou_dot','name' => 'RhoU_dot'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'Rho_2nd_1/rho_2nd_1','name' => 'Rho_2nd_1'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'Rho_2nd_2/rho_2nd_2','name' => 'Rho_2nd_2'},'name' => 'option','type' => 'e'}],'attrib' => {'input' => 'select','type' => 'string','name' => 'TypeTransient'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseSunEarth','default' => 'F'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'real','name' => 'xEarth'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'yEarth'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'zEarth'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'InvD2Ray'},'name' => 'parameter','type' => 'e'}],'attrib' => {'expr' => '$UseSunEarth'},'name' => 'if','type' => 'e'}],'attrib' => {'expr' => '$TypeRefine =~ /transient/i'},'name' => 'if','type' => 'e'}],'attrib' => {'to' => '$nRefineCrit','from' => '1'},'name' => 'for','type' => 'e'},{'content' => '
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
','type' => 't'}],'attrib' => {'name' => 'AMRCRITERIA'},'name' => 'command','type' => 'e'}],'attrib' => {'name' => 'AMR PARAMETERS'},'name' => 'commandgroup','type' => 'e'},{'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  SCHEME PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
','type' => 't'},{'content' => [{'content' => [{'content' => [],'attrib' => {'name' => '1','default' => 'T'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => '2'},'name' => 'option','type' => 'e'}],'attrib' => {'input' => 'select','type' => 'integer','name' => 'nOrder'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'value' => 'Roe/roe/1','name' => 'Roe'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'Rusanov/rusanov/2/TVDLF','name' => 'Rusanov','default' => 'T'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'Linde/linde/3/HLLEL','name' => 'Linde'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'Sokolov/sokolov/4/AW','name' => 'Sokolov'},'name' => 'option','type' => 'e'}],'attrib' => {'input' => 'select','type' => 'string','name' => 'TypeFlux'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [{'content' => [],'attrib' => {'name' => 'minmod','default' => 'T'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'beta'},'name' => 'option','type' => 'e'}],'attrib' => {'input' => 'select','type' => 'string','name' => 'TypeLimiter'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '2','type' => 'real','name' => 'LimiterBeta','default' => '1.2','min' => '1','if' => '$TypeLimiter eq \'beta\''},'name' => 'parameter','type' => 'e'}],'attrib' => {'expr' => '$nOrder == 2'},'name' => 'if','type' => 'e'},{'content' => '
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
','type' => 't'}],'attrib' => {'name' => 'SCHEME'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseNonConservative','default' => 'T'},'name' => 'parameter','type' => 'e'},{'content' => '
#NONCONSERVATIVE
T		UseNonConservative

! For Earth the default is using non-conservative equations 
! (close to the body).
','type' => 't'}],'attrib' => {'name' => 'NONCONSERVATIVE'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'max' => '3','type' => 'integer','name' => 'nConservCrit','default' => '1','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [{'content' => [],'attrib' => {'value' => 'r/R/radius/Radius','name' => 'radius','default' => 'T'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'parabola/paraboloid','name' => 'parabola'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'p/P','name' => 'p'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'gradp/GradP','name' => 'grad P'},'name' => 'option','type' => 'e'}],'attrib' => {'input' => 'select','type' => 'string','name' => 'TypeConservCrit'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'rConserv','default' => '2*$rBody','min' => '$rBody','if' => '$TypeConservCrit =~ /^r|radius$/i'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'xParabolaConserv','default' => '6','min' => '0','if' => '$TypeConservCrit =~ /^parabol/i'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'yParabolaConserv','default' => '36','min' => '0','if' => '$TypeConservCrit =~ /^parabol/i'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'pCoeffConserv','default' => '0.05','min' => '0','if' => '$TypeConservCrit =~ /^p$/i'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'GradPCoeffConserv','default' => '0.1','min' => '0','if' => '$TypeConservCrit =~ /gradp/i'},'name' => 'parameter','type' => 'e'}],'attrib' => {'to' => '$nConservCrit','from' => '1'},'name' => 'for','type' => 'e'},{'content' => '
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
! \'r\'        - radial distance of the cell is less than rConserv
! \'parabola\' - x less than xParabolaConserv + (y**2+z**2)/yParabolaConserv
!
! Default values are nConservCrit = 1 with TypeConservCrit = \'r\'
! and rConserv=2*rBody, where rBody has a problem dependent default.
','type' => 't'}],'attrib' => {'name' => 'CONSERVATIVECRITERIA'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseUpdateCheck','default' => 'T'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'max' => '100','type' => 'real','name' => 'RhoMinPercent','default' => '40','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'RhoMaxPercent','default' => '400','min' => '100'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '100','type' => 'real','name' => 'pMinPercent','default' => '40','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'pMaxPercent','default' => '400','min' => '100'},'name' => 'parameter','type' => 'e'}],'attrib' => {'expr' => '$UseUpdateCheck'},'name' => 'if','type' => 'e'},{'content' => '
#UPDATECHECK
T			UseUpdateCheck
40.			RhoMinPercent
400.			RhoMaxPercent
40.			pMinPercent
400.			pMaxPercent

! Default values are shown.  This will adjust the timestep so that
! density and pressure cannot change by more than the given percentages
! in a single timestep.
','type' => 't'}],'attrib' => {'name' => 'UPDATECHECK'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [{'content' => [],'attrib' => {'name' => '1','default' => 'T'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => '2'},'name' => 'option','type' => 'e'}],'attrib' => {'input' => 'select','type' => 'integer','name' => 'nOrderProlong'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'value' => 'lr','name' => 'left-right','default' => 'T'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'central'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => 'minmod'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'lr2','name' => 'left-right extrapolate'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'central2','name' => 'central    extrapolate'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'minmod2','name' => 'minmod     extrapolate'},'name' => 'option','type' => 'e'}],'attrib' => {'input' => 'select','type' => 'string','name' => 'TypeProlong','if' => '$nOrderProlong==2'},'name' => 'parameter','type' => 'e'},{'content' => '
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
','type' => 't'}],'attrib' => {'name' => 'PROLONGATION'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [{'content' => [],'attrib' => {'value' => 'allopt','name' => 'm_p_cell FACES ONLY','default' => 'T'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'all','name' => 'm_p_cell'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'opt','name' => 'm_p_dir FACES ONLY'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'dir','name' => 'm_p_dir group by directions'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'face','name' => 'm_p_dir group by faces     '},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'min','name' => 'm_p_dir group by kind and face'},'name' => 'option','type' => 'e'}],'attrib' => {'input' => 'select','type' => 'string','name' => 'TypeMessagePass'},'name' => 'parameter','type' => 'e'},{'content' => '
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
','type' => 't'}],'attrib' => {'alias' => 'OPTIMIZE','name' => 'MESSAGEPASS'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseBorisCorrection','default' => 'F'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '1','type' => 'real','name' => 'BorisClightFactor','default' => '1','min' => '0','if' => '$UseBorisCorrection'},'name' => 'parameter','type' => 'e'},{'content' => '
#BORIS
T			UseBorisCorrection
1.0			BorisClightFactor !Only if UseBorisCorrection is true

! Default is boris_correction=.false.
! Use semi-relativistic MHD equations with speed of light reduced by
! the BorisClightFactor. Set BorisClightFactor=1.0 for true semi-relativistic
! MHD. Gives the same steady state as normal MHD analytically, but there
! can be differences due to discretization errors. 
! You can use either Boris or BorisSimple but not both. 
','type' => 't'}],'attrib' => {'name' => 'BORIS'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseBorisSimple','default' => 'F'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '1','type' => 'real','name' => 'BorisClightFactor','default' => '1','min' => '0','if' => '$UseBorisSimple'},'name' => 'parameter','type' => 'e'},{'content' => '
#BORISSIMPLE
T			UseBorisSimple
0.05			BorisClightFactor !Only if UseBorisSimple is true

! Default is UseBorisSimple=.false. 
! Use simplified semi-relativistic MHD with speed of light reduced by the
! BorisClightFactor. This is only useful with BorisClightFactor less than 1.
! Should give the same steady state as normal MHD, but there can be a
! difference due to discretization errors.
! You can use either Boris or BorisSimple but not both. 
','type' => 't'}],'attrib' => {'name' => 'BORISSIMPLE'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseDivbSource','default' => 'T'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseDivbDiffusion','default' => 'F'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseProjection','default' => 'F'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseConstrainB','default' => 'F'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => '
		At least one of the options should be true.
	','type' => 't'}],'attrib' => {'expr' => '$UseDivbSource or $UseDivbDiffusion or $UseProjection or $UseConstrainB'},'name' => 'rule','type' => 'e'},{'content' => [{'content' => '
		If UseProjection is true, all others should be false.
	','type' => 't'}],'attrib' => {'expr' => 'not($UseProjection and ($UseDivbSource or $UseDivbDiffusion or $UseConstrainB))'},'name' => 'rule','type' => 'e'},{'content' => [{'content' => '
		If UseConstrainB is true, all others should be false.
	','type' => 't'}],'attrib' => {'expr' => 'not($UseConstrainB and ($UseDivbSource or $UseDivbDiffusion or $UseProjection))'},'name' => 'rule','type' => 'e'},{'content' => '
	
#DIVB
T			UseDivbSource
F			UseDivbDiffusion	
F			UseProjection           
F			UseConstrainB           

! Default values are shown above.
! If UseProjection is true, all others should be false.
! If UseConstrainB is true, all others should be false.
! At least one of the options should be true.
','type' => 't'}],'attrib' => {'name' => 'DIVB'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseB0Source','default' => 'T'},'name' => 'parameter','type' => 'e'},{'content' => '
#DIVBSOURCE
T			UseB0Source

! Add extra source terms related to the non-zero divergence and curl of B0.
! Default is true.
','type' => 't'}],'attrib' => {'name' => 'DIVBSOURCE'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'max' => '1','type' => 'real','name' => 'DivbDiffCoeff','default' => '0.1666667','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => '
#DIVBDIFFUSION
0.1666667		DivbDiffCoeff

! Default value is shown above. 1.0/6.0
! If divb_diffcoeff .gt. 0.5 then cfl .lt. 0.5/DivbDiffCoeff is required!
','type' => 't'}],'attrib' => {'name' => 'DIVBDIFFUSION'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [{'content' => [],'attrib' => {'value' => 'cg','name' => 'Conjugate Gradients','default' => 'T'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'bicgstab','name' => 'BiCGSTAB'},'name' => 'option','type' => 'e'}],'attrib' => {'input' => 'select','type' => 'string','name' => 'TypeProjectIter'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'value' => 'rel','name' => 'Relative norm','default' => 'T'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'max','name' => 'Maximum error'},'name' => 'option','type' => 'e'}],'attrib' => {'input' => 'select','type' => 'string','name' => 'TypeProjectStop'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '1','type' => 'real','name' => 'RelativeLimit','default' => '0.1','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'AbsoluteLimit','default' => '0.0','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'integer','name' => 'MaxMatvec','default' => '50','min' => '1'},'name' => 'parameter','type' => 'e'},{'content' => '
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

','type' => 't'}],'attrib' => {'name' => 'PROJECTION'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'max' => '1','type' => 'real','name' => 'pRatioLow','default' => '0.01','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '1','type' => 'real','name' => 'pRatioHigh','default' => '0.1','min' => '$pRatioLow'},'name' => 'parameter','type' => 'e'},{'content' => '
#CORRECTP
0.01			pRatioLow
0.1			pRatioHigh

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
','type' => 't'}],'attrib' => {'name' => 'CORRECTP'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseAccurateIntegral','default' => 'T'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseAccurateTrace','default' => 'F'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '60','type' => 'real','name' => 'DtExchangeRay','default' => '0.1','min' => '0.01'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'integer','name' => 'DnRaytrace','default' => '1','min' => '1'},'name' => 'parameter','type' => 'e'},{'content' => '

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
','type' => 't'}],'attrib' => {'name' => 'RAYTRACE'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'real','name' => 'TauCoupleIm','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => '

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
','type' => 't'}],'attrib' => {'name' => 'IM'},'name' => 'command','type' => 'e'}],'attrib' => {'name' => 'SCHEME PARAMETERS'},'name' => 'commandgroup','type' => 'e'},{'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  PHYSICS PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
','type' => 't'},{'content' => [{'content' => [],'attrib' => {'type' => 'real','name' => 'Gamma','default' => '1.6666666667','min' => '1'},'name' => 'parameter','type' => 'e'},{'content' => '
#GAMMA
1.6666666667		Gamma

! The adiabatic index (ratio of the specific heats for fixed pressure
! and volume. The default value is 5.0/3.0, which is valid for
! monoatomic gas or plasma.
','type' => 't'}],'attrib' => {'name' => 'GAMMA','if' => '$_IsFirstSession'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'real','name' => 'RhoLeft','default' => '1','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'UnLeft','default' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'Ut1Left','default' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'Ut2Left','default' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'BnLeft','default' => '0.75'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'Bt1Left','default' => '1'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'Bt2Left','default' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'pRight','default' => '1','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'RhoRight','default' => '0.125','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'UnRight','default' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'Ut1Right','default' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'Ut2Right','default' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'BnRight','default' => '0.75'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'Bt1Right','default' => '-1'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'Bt2Right','default' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'pRight','default' => '0.1','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'value' => '0','name' => 'no rotation','default' => 'T'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => '0.25'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => '0.3333333333333'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => '0.5'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => '1'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => '2'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => '3'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => '4'},'name' => 'option','type' => 'e'}],'attrib' => {'input' => 'select','type' => 'real','name' => 'ShockSlope'},'name' => 'parameter','type' => 'e'},{'content' => '
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
','type' => 't'}],'attrib' => {'name' => 'SHOCKTUBE'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'real','name' => 'SwRhoDim','default' => '5','min' => '-1'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'SwTDim','default' => '181712.175','min' => '-1'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '0','type' => 'real','name' => 'SwUxDim','default' => '-400'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'SwUyDim','default' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'SwUzDim','default' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'SwBxDim','default' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'SwByDim','default' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'SwBzDim','default' => '5'},'name' => 'parameter','type' => 'e'},{'content' => '
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
','type' => 't'}],'attrib' => {'name' => 'SOLARWIND','if' => '$_IsFirstSession'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseUpstreamInputFile','default' => 'F'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'length' => '100','type' => 'string','name' => 'NameUpstreamFile'},'name' => 'parameter','type' => 'e'}],'attrib' => {'expr' => '$UseUpstreamInputFile'},'name' => 'if','type' => 'e'},{'content' => [{'content' => '
		Upstream file $NameUpstreamFile must exist
	','type' => 't'}],'attrib' => {'expr' => '-f $NameUpstreamFile'},'name' => 'rule','type' => 'e'},{'content' => '
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
','type' => 't'}],'attrib' => {'name' => 'UPSTREAM_INPUT_FILE'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseBody','default' => 'F'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'real','name' => 'rBody','default' => '3','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'real','name' => 'rCurrents','default' => '4','min' => '-1'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'BodyRhoDim','default' => '1','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'BodyTDim','default' => '10000','min' => '0'},'name' => 'parameter','type' => 'e'}],'attrib' => {'expr' => '$_NameComp eq \'GM\''},'name' => 'if','type' => 'e'}],'attrib' => {'expr' => '$UseBody'},'name' => 'if','type' => 'e'},{'content' => '
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
','type' => 't'}],'attrib' => {'alias' => 'MAGNETOSPHERE','name' => 'BODY','if' => '$_IsFirstSession'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseGravity','default' => 'F'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'value' => '0','name' => 'central mass','default' => 'T'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => '1','name' => 'X direction'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => '2','name' => 'Y direction'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => '3','name' => 'Z direction'},'name' => 'option','type' => 'e'}],'attrib' => {'input' => 'select','type' => 'integer','name' => 'iDirGravity','if' => '$UseGravity'},'name' => 'parameter','type' => 'e'},{'content' => '
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

','type' => 't'}],'attrib' => {'name' => 'GRAVITY','if' => '$_IsFirstSession'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseMassLoading','default' => 'F'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'logical','name' => 'DoAccelerateMassLoading','default' => 'F'},'name' => 'parameter','type' => 'e'},{'content' => '
#MASSLOADING
F			UseMassLoading
F			DoAccelerateMassLoading
','type' => 't'}],'attrib' => {'name' => 'MASSLOADING'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseHeatFlux','default' => 'F'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseSpitzerForm','default' => 'T'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'real','name' => 'Kappa0Heat','default' => '1.23E-11'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'Kappa0Heat','default' => '2.5'},'name' => 'parameter','type' => 'e'}],'attrib' => {'expr' => 'not $UseSpitzerForm'},'name' => 'if','type' => 'e'},{'content' => '
#HEATFLUX
T		UseHeatFlux
F		UseSpitzerForm
1.23E-11	Kappa0Heat [W/m/K]	! Only if not UseSpitzerForm
2.50E+00	ExponentHeat [-]	! Only if not UseSpitzerForm
','type' => 't'}],'attrib' => {'name' => 'HEATFLUX'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseResistFlux','default' => 'F'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseSpitzerForm','default' => 'T'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [{'content' => [],'attrib' => {'value' => 'Localized/localized','name' => 'localized'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'value' => 'Constant/constant','name' => 'constant','default' => 'T'},'name' => 'option','type' => 'e'}],'attrib' => {'input' => 'select','type' => 'string','name' => 'TypeResist'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'real','name' => 'Eta0Resist','default' => '9.69953E+8'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'Alpha0Resist','default' => '150'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'yShiftResist','default' => '0.5'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'TimeInitRise','default' => '0.05'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'TimeConstLev','default' => '1'},'name' => 'parameter','type' => 'e'}],'attrib' => {'expr' => '$TypeResist =~ /localized/i'},'name' => 'if','type' => 'e'}],'attrib' => {'expr' => 'not $UseSpitzerForm'},'name' => 'if','type' => 'e'},{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseAnomResist','default' => 'F'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'real','name' => 'Eta0AnomResist','default' => '1.93991E+09'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'EtaAnomMaxResist','default' => '1.93991E+10'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'ThresholdFactorResist','default' => '1'},'name' => 'parameter','type' => 'e'}],'attrib' => {'expr' => '$UseAnomResist'},'name' => 'if','type' => 'e'},{'content' => '
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
','type' => 't'}],'attrib' => {'name' => 'RESISTIVEFLUX'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseDefaultUnits','default' => 'T'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'Grav0Diss','default' => '2.635620E-02'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'Beta0Diss','default' => '1.640000E-01'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'Length0Diss','default' => '1.500000E+06'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'Time0Diss','default' => '1.159850E+01'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'Rho0Diss','default' => '5.019000E-11'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'Tem0Diss','default' => '1.000000E+05'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'Theta0Diss','default' => '6.000000E-01'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'Delta0Diss','default' => '2.500000E+01'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'EpsilonDiss','default' => '7.000000E+00'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'RhoDifDiss','default' => '4.500000E+00'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'yShiftDiss','default' => '4.000000E-01'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'ScaleHeightDiss','default' => '5.000000E-01'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'ScaleFactorDiss','default' => '1.159850E+01'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'BZ0iss','default' => '5.000000E-01'},'name' => 'parameter','type' => 'e'},{'content' => '
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
','type' => 't'}],'attrib' => {'name' => 'TESTDISSMHD','if' => '$_IsFirstSession'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseBody2','default' => 'F'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'real','name' => 'rBody2','default' => '0.1','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '$xMax','type' => 'real','name' => 'xBody2','default' => '-40','min' => '$xMin'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '$yMax','type' => 'real','name' => 'yBody2','default' => '0','min' => '$yMin'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '$zMax','type' => 'real','name' => 'zBody2','default' => '0','min' => '$zMin'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'rCurrents2','default' => '1.3*$rBody2','min' => '$rBody2'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'RhoDimBody2','default' => '5','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'tDimBody2','default' => '25000','min' => '0'},'name' => 'parameter','type' => 'e'}],'attrib' => {'expr' => '$UseBody2'},'name' => 'if','type' => 'e'},{'content' => '

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
','type' => 't'}],'attrib' => {'name' => 'SECONDBODY','if' => '$_IsFirstSession'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'real','name' => 'BdpDimBody2x'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'BdpDimBody2y'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'BdpDimBody2z'},'name' => 'parameter','type' => 'e'},{'content' => '

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
','type' => 't'}],'attrib' => {'name' => 'DIPOLEBODY2','if' => '$_IsFirstSession'},'name' => 'command','type' => 'e'}],'attrib' => {'name' => 'PHYSICS PARAMETERS'},'name' => 'commandgroup','type' => 'e'},{'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!! SOLAR PROBLEM TYPES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
','type' => 't'},{'content' => [{'content' => [],'attrib' => {'type' => 'real','name' => 'BodyTDim','default' => '2.85E06','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'BodyRhoDim','default' => '1.50E8','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'qSun','default' => '25.0','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'tHeat','default' => '1.75','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'rHeat','default' => '1.0','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'SigmaHeat','default' => '4.5','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'logical','name' => 'DoInitRope','default' => 'F'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'real','name' => 'CmeA','default' => '0.7','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'CmeR1','default' => '1.2','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'CmeR0','default' => '1.0','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'CmeA1','default' => '0.23'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'CmeAlpha','default' => '0.0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'CmeRho1','default' => '2.5E-12','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'CmeRho2','default' => '2.0E-13','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '10','type' => 'real','name' => 'ModulationRho','default' => '0.0','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '10','type' => 'real','name' => 'ModulationP','default' => '0.0','min' => '0'},'name' => 'parameter','type' => 'e'}],'attrib' => {'expr' => '$DoInitRope'},'name' => 'if','type' => 'e'},{'content' => '
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

This command defines the heliosphere parameters with a CME model.
The coronal eruptive event generator is based on the
Gibson and Low (GL) analytical solution prescribing a
three-dimensional twisted magnetic flux rope in
hydrostatic equilibrium in the presence of gravity.
The GL solution is described in the Astrophysical
Journal, volume 493, page 460.
This flux rope is formed by applying a mathematical
stretching operation to axisymmetric speromak flux
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
rope being bodily expelled from the conona.  Eruption energy
increases with flux rope size, field strength, stretching
deformation and the buoyancy of the flux rope.
','type' => 't'}],'attrib' => {'name' => 'HELIOSPHERE','if' => '$_IsFirstSession and $_NameComp ne \'GM\''},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'real','name' => 'DtUpdateB0','default' => '0.0001','min' => '-1'},'name' => 'parameter','type' => 'e'},{'content' => '

#HELIOUPDATEB0
-1.0			DtUpdateB0 [s]

Set the frequency of updating the B0 field for the solar corona.
A negative value means that the B0 field is not updated.
','type' => 't'}],'attrib' => {'name' => 'HELIOUPDATEB0','if' => '$_NameComp ne \'GM\''},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'real','name' => 'HelioDipoleStrength'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '90','type' => 'real','name' => 'HelioDipoleTilt','default' => '0','min' => '-90'},'name' => 'parameter','type' => 'e'},{'content' => '

#HELIODIPOLE
-3.0                    HelioDipoleStrength [G]
 0.0                    HelioDipoleTilt     [deg]

! Variable HelioDipoleStrength defines the equatorial field strength in Gauss,
! while HelioDipoleTilt is the tilt relative to the ecliptic North 
! (negative sign means towards the planet) in degrees.
!
! Default value is HelioDipoleStrength = 0.0.
','type' => 't'}],'attrib' => {'name' => 'HELIODIPOLE','if' => '$_NameComp ne \'GM\''},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseInertialFrame','default' => 'T'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'logical','name' => 'UseRotatingBC','default' => 'T'},'name' => 'parameter','type' => 'e'}],'attrib' => {'expr' => '$UseInertialFrame'},'name' => 'if','type' => 'e'},{'content' => '

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
','type' => 't'}],'attrib' => {'alias' => 'INERTIAL','name' => 'HELIOROTATION','if' => '$_IsFirstSession and $_NameComp ne \'GM\''},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'logical','name' => 'DoSendMHD','default' => 'F'},'name' => 'parameter','type' => 'e'},{'content' => '

#HELIOTEST
F			DoSendMHD

! If DoSendMHD is true, IH sends the real MHD solution to GM in the coupling.
! If DoSendMHD is false then the values read from the IMF file are sent,
! so there is no real coupling. Mostly used for testing the framework.
!
! Default value is true, ie. real coupling.
','type' => 't'}],'attrib' => {'name' => 'HELIOTEST','if' => '$_NameComp ne \'GM\''},'name' => 'command','type' => 'e'},{'content' => [{'content' => [{'content' => [],'attrib' => {'name' => 'Low','default' => 'T'},'name' => 'option','type' => 'e'}],'attrib' => {'input' => 'select','type' => 'string','name' => 'TypeCme'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'CmeA','default' => '0.7','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'CmeR1','default' => '1.2','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'CmeR0','default' => '1.0','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'CmeA1','default' => '0.23'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'CmeAlpha','default' => '0.0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'CmeRho1','default' => '2.5E-12','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'CmeRho2','default' => '2.0E-13','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'CmeB1Dim','default' => '1.0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'CmeUErupt','default' => '4.0E5','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => '
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
stretching operation to axisymmetric speromak flux
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
rope being bodily expelled from the conona.  Eruption energy
increases with flux rope size, field strength, stretching
deformation and the buoyancy of the flux rope.
Default values are shown above for the GL flux rope CME model.
','type' => 't'}],'attrib' => {'name' => 'CME','if' => '$_IsFirstSession and $_NameComp ne \'GM\''},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'type' => 'real','name' => 'tArcDim','default' => '1.0E6','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'RhoArcDim','default' => '1.0E-12','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'bArcDim','default' => '0.718144','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'ByArcDim','default' => '1.0E6'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'UzArcDim','default' => '5.0E3'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'Phi0Arc','default' => '0.5'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'MuArc','default' => '1.3'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'ExpArc','default' => '3','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'WidthArc','default' => '0.5','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => '
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
','type' => 't'}],'attrib' => {'name' => 'ARCADE','if' => '$_IsFirstSession and $_NameComp ne \'GM\''},'name' => 'command','type' => 'e'}],'attrib' => {'name' => 'SOLAR PROBLEM TYPES'},'name' => 'commandgroup','type' => 'e'},{'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!! COMET PROBLEM TYPE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
','type' => 't'},{'content' => [{'content' => [],'attrib' => {'type' => 'real','name' => 'ProdRate','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'UrNeutral','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'AverageMass','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'IonizationRate','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'type' => 'real','name' => 'kFriction','min' => '0'},'name' => 'parameter','type' => 'e'},{'content' => '
#COMET
1.0E28		ProdRate    - Production rate (#/s)
1.0		UrNeutral   - neutral radial outflow velocity (km/s)
17.0		AverageMass - average particle mass (amu)
1.0E-6		IonizationRate (1/s)
1.7E-9		kFriction - ion-neutral friction rate coefficient (cm^3/s)

! Only used by problem_comet.  Defaults are as shown.
','type' => 't'}],'attrib' => {'name' => 'COMET','if' => '$_IsFirstSession'},'name' => 'command','type' => 'e'}],'attrib' => {'name' => 'COMET PROBLEM TYPE'},'name' => 'commandgroup','type' => 'e'},{'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! SCRIPT COMMANDS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
','type' => 't'},{'content' => [{'content' => [],'attrib' => {'length' => '100','type' => 'string','name' => 'NameIncludeFile','default' => 'Param/'},'name' => 'parameter','type' => 'e'},{'content' => '

#INCLUDE
Param/SSS_3000		NameIncludeFile

! Include a library file from Param/ or any file from anywhere else.
','type' => 't'}],'attrib' => {'name' => 'INCLUDE'},'name' => 'command','type' => 'e'}],'attrib' => {'name' => 'SCRIPT COMMANDS'},'name' => 'commandgroup','type' => 'e'},{'content' => [{'content' => '
	Either command #SOLARWIND or #UPSTREAM_INPUT_FILE must be used!
','type' => 't'}],'attrib' => {'expr' => '($SwRhoDim > 0) or $UseUpstreamInputFile or $_NameComp ne \'GM\''},'name' => 'rule','type' => 'e'},{'content' => [{'content' => '
	Part implicit scheme requires more than 1 implicit block!
','type' => 't'}],'attrib' => {'expr' => '$MaxImplBlock>1 or not $UsePartImplicit or not $MaxImplBlock'},'name' => 'rule','type' => 'e'},{'content' => [{'content' => '
	Full implicit scheme should be used with equal number of 
	explicit and implicit blocks!
','type' => 't'}],'attrib' => {'expr' => '$MaxImplBlock==$MaxBlock or not $UseFullImplicit'},'name' => 'rule','type' => 'e'}],'attrib' => {'name' => 'BATSRUS: GM, SC and IH Components'},'name' => 'commandList','type' => 'e'}];