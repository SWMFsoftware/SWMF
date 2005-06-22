#^CFG FILE _FALSE_
$tree = [{'name' => 'commandList','attrib' => {'name' => 'BATSRUS: GM, SC and IH Components'},'type' => 'e','content' => [{'type' => 't','content' => '

List of MH (GM, IH and SC) commands used in the PARAM.in file


'},{'name' => 'set','attrib' => {'name' => 'nI','type' => 'integer','value' => '$_GridSize[0]'},'type' => 'e','content' => []},{'name' => 'set','attrib' => {'name' => 'nJ','type' => 'integer','value' => '$_GridSize[1]'},'type' => 'e','content' => []},{'name' => 'set','attrib' => {'name' => 'nK','type' => 'integer','value' => '$_GridSize[2]'},'type' => 'e','content' => []},{'name' => 'set','attrib' => {'name' => 'MaxBlock','type' => 'integer','value' => '$_GridSize[3]'},'type' => 'e','content' => []},{'name' => 'set','attrib' => {'name' => 'MaxImplBlock','type' => 'integer','value' => '$_GridSize[4]'},'type' => 'e','content' => []},{'name' => 'set','attrib' => {'name' => 'MaxBlockALL','type' => 'integer','value' => '$_nProc and $MaxBlock and $_nProc*$MaxBlock'},'type' => 'e','content' => []},{'name' => 'set','attrib' => {'name' => 'NameRestartOutDir','type' => 'string','value' => '$_NameComp/restartOUT'},'type' => 'e','content' => []},{'name' => 'set','attrib' => {'name' => 'NamePlotDir','type' => 'string','value' => '$_NameComp/IO2'},'type' => 'e','content' => []},{'name' => 'commandgroup','attrib' => {'name' => 'STAND ALONE MODE'},'type' => 'e','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! STAND ALONE PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

'},{'name' => 'command','attrib' => {'name' => 'NEWPARAM','if' => '$_IsStandAlone'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => 'T','name' => 'UseNewParam','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => 'T','name' => 'UseNewAxes','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => 'T','name' => 'DoTimeAccurate','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => 'T','name' => 'UseRotatingBc','type' => 'logical'},'type' => 'e','content' => []},{'type' => 't','content' => '

#NEWPARAM
T			UseNewParam
T			UseNewAxes
T			DoTimeAccurate
T			UseRotatingBc

This command can be used to make the standalone code backwards compatible.

If UseNewParam is true, the time frequencies of various commands 
(SAVEPLOT, SAVELOGFILE, STOP etc.) are always read, irrespective of the value 
of DoTimeAccurate and the DoTimeAccurate logical can be set with the TIMEACCURATE command.

If UseNewParam is false, the time frequencies are only read when 
DoTimeAccurate is true, and DoTimeAccurate can be set as the first 
parameter of the TIMESTEPPING command.

If UseNewAxes is true, the planet\'s rotational and magnetic axes 
are set by the new algorithms found in
share/Library/src/CON\\_axes, the planet data is set and
stored by share/Library/src/CON\\_planet, and magnetic field information and
mapping is provided by share/Library/src/CON\\_planet_field, 
and the rotational speed
of the planet is calculated using $v_{\\phi}=\\Omega \\times r$.

If UseNewAxes is false, the original algorithms in 
GM/BATSRUS/src/ModCompatibility are used. 
Some of these algorithms are inaccurate, some of them contain bugs and
some of them are inefficient. The algorithms were kept for the sake 
of backwards compatibility.

The DoTimeAccurate and UseRotatingBc parameters can be set elsewhere, but their
default values can be set here. This is again useful for backwards 
compatibility,
since BATSRUS v7.72 and earlier has DoTimeAccurate=F and UseRotatingBc=F as the
default, while SWMF has the default values DoTimeAccurate=T and UseRotatingBc=T
(consistent with the assumption that the default behavior is as realistic as 
possible).

The default values depend on how the standalone code was installed
(make install STANDALONE=???). For STANDALONE=gm and STANDALONE=ih
all the logicals have true default values, while for 
STANDALONE=sc all logicals are true except for UseRotatingBc 
(consistent with the SWMF).
For STANDALONE=old and STANDALONE=oldtest the default values are false 
(consistent with BATSRUS v7.72 and earlier).
'}]},{'name' => 'command','attrib' => {'name' => 'COMPONENT'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'input' => 'select','name' => 'NameComp','type' => 'string'},'type' => 'e','content' => [{'name' => 'option','attrib' => {'default' => 'T','name' => 'SC','if' => '$_NameComp eq \'SC\''},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'default' => 'T','name' => 'IH','if' => '$_NameComp eq \'IH\''},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'default' => 'T','name' => 'GM','if' => '$_NameComp eq \'GM\''},'type' => 'e','content' => []}]},{'type' => 't','content' => '

#COMPONENT
GM			NameComp

This command is only used in the stand alone mode.

The NameComp variable contains the two-character component ID
for the component which BATSRUS is representing.
If NameComp does not agree with the value of the NameThisComp
variable, BATSRUS stops with an error message.
This command is saved into the restart header file for consistency check.

There is no default value: if the command is not given, the component ID is not checked.
'}]},{'name' => 'command','attrib' => {'name' => 'DESCRIPTION','if' => '$_IsStandAlone'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'length' => '100','name' => 'StringDescription','type' => 'string'},'type' => 'e','content' => []},{'type' => 't','content' => '

#DESCRIPTION
This is a test run for Jupiter with no rotation.

This command is only used in the stand alone mode.

The StringDescription string can be used to describe the simulation
for which the parameter file is written. The #DESCRIPTION command and
the StringDescription string are saved into the restart file,
which helps in identifying the restart files.

The default value is ``Please describe me!", which is self explanatory.
'}]},{'name' => 'command','attrib' => {'name' => 'ECHO','if' => '$_IsStandAlone'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'DoEcho','type' => 'logical'},'type' => 'e','content' => []},{'type' => 't','content' => '

#ECHO
T                       DoEcho

This command is only used in the stand alone mode.

If the DoEcho variable is true, the input parameters are echoed back.
The default value for DoEcho is .false., but it is a good idea to
set it to true at the beginning of the PARAM.in file.
'}]},{'name' => 'command','attrib' => {'name' => 'PROGRESS','if' => '$_IsStandAlone'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => '10','name' => 'DnProgressShort','min' => '-1','type' => 'integer'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '100','name' => 'DnProgressLong','min' => '-1','type' => 'integer'},'type' => 'e','content' => []},{'type' => 't','content' => '
#PROGRESS
10			DnProgressShort
100			DnProgressLong

The frequency of short and long progress reports for BATSRUS in
stand alone mode. These are the defaults. Set -1-s for no progress reports.
'}]},{'name' => 'command','attrib' => {'name' => 'TIMEACCURATE','if' => '$_IsStandAlone'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => 'T','name' => 'DoTimeAccurate','type' => 'logical'},'type' => 'e','content' => []},{'type' => 't','content' => '

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
'}]},{'name' => 'command','attrib' => {'multiple' => 'T','name' => 'BEGIN_COMP','if' => '$_IsStandAlone'},'type' => 'e','content' => [{'type' => 't','content' => '

This command is allowed in stand alone mode only for the sake of the 
test suite, which contains these commands when the framework is tested.
'}]},{'name' => 'command','attrib' => {'multiple' => 'T','name' => 'END_COMP','if' => '$_IsStandAlone'},'type' => 'e','content' => [{'type' => 't','content' => '

This command is allowed in stand alone mode only for the sake of the 
test suite, which contains these commands when the framework is tested.
'}]},{'name' => 'command','attrib' => {'name' => 'RUN','if' => '$_IsStandAlone'},'type' => 'e','content' => [{'type' => 't','content' => '

#RUN

This command is only used in stand alone mode.

The #RUN command does not have any parameters. It signals the end
of the current session, and makes BATSRUS execute the session with
the current set of parameters. The parameters for the next session
start after the #RUN command. For the last session there is no
need to use the #RUN command, since the #END command or simply
the end of the PARAM.in file makes BATSRUS execute the last session.
'}]},{'name' => 'command','attrib' => {'name' => 'END'},'type' => 'e','content' => [{'type' => 't','content' => '

#END

The #END command signals the end of the included file or the
end of the PARAM.in file. Lines following the #END command are
ignored. It is not required to use the #END command. The end
of the included file or PARAM.in file is equivalent with an 
#END command in the last line.
'}]}]},{'name' => 'commandgroup','attrib' => {'name' => 'PLANET COMMANDS'},'type' => 'e','content' => [{'type' => 't','content' => '
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
These modified commands cannot precede the #PLANET command!

'},{'name' => 'command','attrib' => {'name' => 'PLANET','if' => '$_IsFirstSession and $_IsStandAlone'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'input' => 'select','name' => 'NamePlanet','type' => 'string'},'type' => 'e','content' => [{'name' => 'option','attrib' => {'default' => 'T','name' => 'Earth','value' => 'EARTH/Earth/earth'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'Saturn','value' => 'SATURN/Saturn/saturn'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'New'},'type' => 'e','content' => []}]},{'name' => 'if','attrib' => {'expr' => '$NamePlanet eq \'New\''},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'name' => 'RadiusPlanet','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'name' => 'MassPlanet','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'name' => 'OmegaPlanet','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'name' => 'TiltRotation','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'input' => 'select','name' => 'TypeBField','type' => 'string'},'type' => 'e','content' => [{'name' => 'option','attrib' => {'name' => 'NONE'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'default' => 'T','name' => 'DIPOLE'},'type' => 'e','content' => []}]}]},{'name' => 'if','attrib' => {'expr' => '$TyepBField eq \'DIPOLE\''},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'name' => 'MagAxisThetaGeo','max' => '180','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'name' => 'MagAxisPhiGeo','max' => '360','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'name' => 'DipoleStrength','type' => 'real'},'type' => 'e','content' => []}]},{'name' => 'rule','attrib' => {'expr' => 'not $PlanetCommand'},'type' => 'e','content' => [{'type' => 't','content' => '
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
'}]},{'name' => 'command','attrib' => {'name' => 'ROTATIONAXIS','if' => '$_IsStandAlone'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => 'T','name' => 'IsRotAxisPrimary','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'if','attrib' => {'expr' => '$IsRotAxisPrimary'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'name' => 'RotAxisTheta','max' => '180','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'name' => 'RotAxisPhi','max' => '360','min' => '0','type' => 'real'},'type' => 'e','content' => []}]},{'name' => 'set','attrib' => {'name' => 'PlanetCommand','type' => 'string','value' => 'ROTATIONAXIS'},'type' => 'e','content' => []},{'type' => 't','content' => '

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
'}]},{'name' => 'command','attrib' => {'name' => 'ROTATION','if' => '$_IsStandAlone'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => 'T','name' => 'UseRotation','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'if','attrib' => {'expr' => '$UseRotation'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'name' => 'RotationPeriod','type' => 'real'},'type' => 'e','content' => []}]},{'name' => 'set','attrib' => {'name' => 'PlanetCommand','type' => 'string','value' => 'MAGNETICAXIS'},'type' => 'e','content' => []},{'type' => 't','content' => '

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
'}]},{'name' => 'command','attrib' => {'name' => 'MAGNETICAXIS','if' => '$_IsStandAlone'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => 'T','name' => 'IsMagAxisPrimary','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'if','attrib' => {'expr' => '$IsMagAxisPrimary'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'name' => 'MagAxisTheta','max' => '180','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'name' => 'MagAxisPhi','max' => '360','min' => '0','type' => 'real'},'type' => 'e','content' => []}]},{'name' => 'set','attrib' => {'name' => 'PlanetCommand','type' => 'string','value' => 'MAGNETICAXIS'},'type' => 'e','content' => []},{'type' => 't','content' => '

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
'}]},{'name' => 'command','attrib' => {'name' => 'DIPOLE','if' => '$_IsStandAlone'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'name' => 'DipoleStrength','type' => 'real'},'type' => 'e','content' => []},{'type' => 't','content' => '

#DIPOLE
-3.11e-5		DipoleStrength [Tesla]

The DipoleStrength variable contains the
magnetic equatorial strength of the dipole magnetic field in Tesla.

The default value is the real dipole strength for the planet.
For the Earth the default is taken to be -31100 nT.
The sign is taken to be negative so that the magnetic axis can
point northward as usual.
'}]},{'name' => 'command','attrib' => {'name' => 'UPDATEB0','if' => '$_IsStandAlone'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => '0.0001','name' => 'DtUpdateB0','min' => '-1','type' => 'real'},'type' => 'e','content' => []},{'type' => 't','content' => '

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
'}]},{'name' => 'command','attrib' => {'name' => 'IDEALAXES','if' => '$_IsStandAlone'},'type' => 'e','content' => [{'type' => 't','content' => '

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
'}]}]},{'name' => 'commandgroup','attrib' => {'name' => 'USER DEFINED INPUT'},'type' => 'e','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!  USER DEFINED INPUT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

'},{'name' => 'command','attrib' => {'name' => 'USER_FLAGS'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'UseUserInnerBcs','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'UseUserSource','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'UseUserPerturbation','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'UseUserOuterBcs','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'UseUserICs','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'UseUserSpecifyRefinement','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'UseUserLogFiles','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'UseUserWritePlot','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'UseUserAMR','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'UseUserEchoInput','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'UseUserB0','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'UseUserSetPhysConst','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'UseUserUpdateStates','type' => 'logical'},'type' => 'e','content' => []},{'type' => 't','content' => '

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
'}]},{'name' => 'command','attrib' => {'name' => 'USERINPUTBEGIN'},'type' => 'e','content' => [{'type' => 't','content' => '

This command signals the beginning of the section of the file which 
is read by the subroutine user\\_read\\_inputs in the user\\_routines.f90 file.
The section ends with the #USERINPUTEND command. There is no XML based parameter
checking in the user section.
'}]},{'name' => 'command','attrib' => {'name' => 'USERINPUTEND'},'type' => 'e','content' => [{'type' => 't','content' => '

This command signals the end of the section of the file which 
is read by the subroutine user\\_read\\_inputs in the user\\_routines.f90 file.
The section begins with the #USERINPUTBEGIN command. There is no XML based parameter
checking in the user section.
'}]}]},{'name' => 'commandgroup','attrib' => {'name' => 'TESTING AND TIMING'},'type' => 'e','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  TESTING AND TIMING PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'name' => 'command','attrib' => {'name' => 'TEST'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'length' => '100','name' => 'TestString','type' => 'string'},'type' => 'e','content' => []},{'type' => 't','content' => '
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
'}]},{'name' => 'command','attrib' => {'name' => 'TESTIJK'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'name' => 'iTest','max' => '$nI+2','min' => '-2','type' => 'integer'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'name' => 'jTest','max' => '$nJ+2','min' => '-2','type' => 'integer'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'name' => 'kTest','max' => '$nK+2','min' => '-2','type' => 'integer'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'name' => 'iBlockTest','max' => '$MaxBlock','min' => '1','type' => 'integer'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'name' => 'iProcTest','min' => '0','type' => 'integer'},'type' => 'e','content' => []},{'type' => 't','content' => '
#TESTIJK
1                       iTest           (cell index for testing)
1                       jTest           (cell index for testing)
1                       kTest           (cell index for testing)
1                       iBlockTest      (block index for testing)
0                       iProcTest       (processor index for testing)

! The location of test info in terms of indices, block and processor number.
! Note that the user should set #TESTIJK or #TESTXYZ, not both.  If both
! are set, the final one in the session will set the test point.
'}]},{'name' => 'command','attrib' => {'name' => 'TESTXYZ'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'name' => 'xTest','max' => '$xMax','min' => '$xMin','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'name' => 'yTest','max' => '$yMax','min' => '$yMin','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'name' => 'zTest','max' => '$zMax','min' => '$zMin','type' => 'real'},'type' => 'e','content' => []},{'type' => 't','content' => '
#TESTXYZ
1.5                     xTest           (X coordinate of cell for testing)
-10.5                   yTest           (Y coordinate of cell for testing)
-10.                    zTest           (Z coordinate of cell for testing)

! The location of test info in terms of coordinates.
! Note that the user should set #TESTIJK or #TESTXYZ, not both.  If both
! are set, the final one in the session will set the test point.
'}]},{'name' => 'command','attrib' => {'name' => 'TESTTIME'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => '-1','name' => 'nIterTest','min' => '-1','type' => 'integer'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '1e30','name' => 'TimeTest','min' => '-1','type' => 'real'},'type' => 'e','content' => []},{'type' => 't','content' => '

#TESTTIME
-1                      nIterTest       (iteration number to start testing)
10.5                    TimeTest        (time to start testing in seconds)

! The time step and physical time to start testing.
'}]},{'name' => 'command','attrib' => {'name' => 'TESTVAR'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'input' => 'select','name' => 'iVarTest','type' => 'integer'},'type' => 'e','content' => [{'name' => 'option','attrib' => {'default' => 'T','name' => 'Rho','value' => '1'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'RhoUx','value' => '2'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'RhoUy','value' => '3'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'RhoUz','value' => '4'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'Bx','value' => '5'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'By','value' => '6'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'Bz','value' => '7'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'e','value' => '8'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'p','value' => '9'},'type' => 'e','content' => []}]},{'type' => 't','content' => '
#TESTVAR
1                       iVarTest

! Index of variable to be tested. Default is rho_="1", i.e. density.
'}]},{'name' => 'command','attrib' => {'name' => 'TESTDIM'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'input' => 'select','name' => 'iVarTest','type' => 'integer'},'type' => 'e','content' => [{'name' => 'option','attrib' => {'name' => 'all','value' => '0'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'default' => 'T','name' => 'x','value' => '1'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'y','value' => '2'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'z','value' => '3'},'type' => 'e','content' => []}]},{'type' => 't','content' => '
#TESTDIM
1                       iDimTest

! Index of dimension/direction to be tested. Default is X dimension.
'}]},{'name' => 'command','attrib' => {'name' => 'STRICT'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => 'T','name' => 'UseStrict','type' => 'logical'},'type' => 'e','content' => []},{'type' => 't','content' => '
#STRICT
T                       UseStrict

! If true then stop when parameters are incompatible. If false, try to
! correct parameters and continue. Default is true, i.e. strict mode
'}]},{'name' => 'command','attrib' => {'name' => 'VERBOSE'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'input' => 'select','name' => 'lVerbose','type' => 'integer'},'type' => 'e','content' => [{'name' => 'option','attrib' => {'name' => 'errors and warnings only','value' => '-1'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'start and end of sessions','value' => '0'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'default' => 'T','name' => 'normal','value' => '1'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'calls on test processor','value' => '10'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'calls on all processors','value' => '100'},'type' => 'e','content' => []}]},{'type' => 't','content' => '
#VERBOSE
-1                      lVerbose

! Verbosity level controls the amount of output to STDOUT. Default level is 1.
!\\\\
!   lVerbose $\\leq -1$ only warnings and error messages are shown.\\\\
!   lVerbose $\\geq  0$ start and end of sessions is shown.\\\\
!   lVerbose $\\leq  1$ a lot of extra information is given.\\\\
!   lVerbose $\\leq 10$ all calls of set_oktest are shown for the test processor.\\\\
!   lVerbose $\\leq 100$ all calls of set_oktest are shown for all processors.\\\\
'}]},{'name' => 'command','attrib' => {'name' => 'DEBUG'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'DoDebug','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'DoDebugGhost','type' => 'logical'},'type' => 'e','content' => []},{'type' => 't','content' => '
#DEBUG
F                       DoDebug         (use it as if(okdebug.and.oktest)...)
F                       DoDebugGhost    (parameter for show_BLK in library.f90)

! Excessive debug output can be controlled by the global okdebug parameter
'}]},{'name' => 'command','attrib' => {'name' => 'CODEVERSION','if' => '$_IsFirstSession'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => '7.50','name' => 'CodeVersion','type' => 'real','min' => '0'},'type' => 'e','content' => []},{'type' => 't','content' => '
#CODEVERSION
7.50                    CodeVersion

! Checks CodeVersion. Prints a WARNING if it differs from the CodeVersion
! defined in ModMain.f90. Used in newer restart header files. 
! Should be given in PARAM.in when reading old restart files, 
! which do not have version info in the header file.
'}]},{'name' => 'command','attrib' => {'name' => 'EQUATION','if' => '$_IsFirstSession'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => 'MHD','length' => '100','name' => 'NameEquation','type' => 'string'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '8','name' => 'nVar','type' => 'integer'},'type' => 'e','content' => []},{'type' => 't','content' => '
#EQUATION
MHD			NameEquation
8			nVar

! Define the equation name and the number of variables.
! If any of these do not agree with the values determined 
! by the code, BATSRUS stops with an error. Used in restart
! header files and can be given in PARAM.in as a check
! and as a description.
'}]},{'name' => 'command','attrib' => {'name' => 'PRECISION','if' => '$_IsFirstSession'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'input' => 'select','name' => 'nByteReal','type' => 'integer'},'type' => 'e','content' => [{'name' => 'option','attrib' => {'default' => '$_nByteReal==4','name' => 'single precision (4)','value' => '4'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'default' => '$_nByteReal==8','name' => 'double precision (8)','value' => '8'},'type' => 'e','content' => []}]},{'name' => 'rule','attrib' => {'expr' => '$nByteReal==$_nByteReal'},'type' => 'e','content' => [{'type' => 't','content' => '
		nByteReal in file must agree with _nByteReal.
	'}]},{'type' => 't','content' => '

#PRECISION
8                       nByteReal

! Define the number of bytes in a real number. If it does not agree
! with the value determined by the code, BATSRUS stops with an error.
! This is used in latest restart header files to check binary 
! compatibility between the restart file and the compiled code
! The command may also be given in PARAM.in to enforce a certain precision.
'}]},{'name' => 'command','attrib' => {'name' => 'CHECKGRIDSIZE','if' => '$_IsFirstSession'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => '$nI','name' => 'nI','max' => '$nI','min' => '$nI','type' => 'integer'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '$nJ','name' => 'nJ','max' => '$nJ','min' => '$nJ','type' => 'integer'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '$nK','name' => 'nK','max' => '$nK','min' => '$nK','type' => 'integer'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'name' => 'MinBlockALL','max' => '$MaxBlockALL','min' => '1','type' => 'integer'},'type' => 'e','content' => []},{'type' => 't','content' => '

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
'}]},{'name' => 'command','attrib' => {'name' => 'BLOCKLEVELSRELOADED'},'type' => 'e','content' => [{'type' => 't','content' => '
#BLOCKLEVELSRELOADED

This command means that the restart file contains the information about
the minimum and maximum allowed refinement levels for each block.
This command is only used in the restart header file.
'}]},{'name' => 'command','attrib' => {'name' => 'TIMING','if' => '$_IsStandAlone'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => 'T','name' => 'UseTiming','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'if','attrib' => {'expr' => '$UseTiming'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'input' => 'select','name' => 'Frequency','type' => 'integer'},'type' => 'e','content' => [{'name' => 'option','attrib' => {'name' => 'none','value' => '-3'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'default' => 'T','name' => 'final only','value' => '-2'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'end of sessions','value' => '-1'},'type' => 'e','content' => []},{'name' => 'optioninput','attrib' => {'default' => '100','name' => 'every X steps','min' => '1'},'type' => 'e','content' => []}]},{'name' => 'parameter','attrib' => {'default' => '-1','name' => 'nDepthTiming','min' => '-1','type' => 'integer'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'input' => 'select','name' => 'TypeTimingReport','type' => 'string'},'type' => 'e','content' => [{'name' => 'option','attrib' => {'default' => '1','name' => 'cumulative','value' => 'cumu'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'list'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'tree'},'type' => 'e','content' => []}]}]},{'type' => 't','content' => '
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
'}]}]},{'name' => 'commandgroup','attrib' => {'name' => 'INITIAL AND BOUNDARY CONDITIONS'},'type' => 'e','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! MAIN INITIAL AND BOUNDARY CONDITION PARAMETERS  !!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'name' => 'command','attrib' => {'required' => 'T','name' => 'PROBLEMTYPE','if' => '$_IsFirstSession'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'input' => 'select','name' => 'iProblem','type' => 'integer'},'type' => 'e','content' => [{'name' => 'option','attrib' => {'name' => 'Uniform','value' => '1'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'Shock tube','value' => '2'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'Heliosphere','value' => '3'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'Comet','value' => '5'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'Rotation','value' => '6'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'Diffusion','value' => '7'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'default' => 'T','name' => 'Earth','value' => '11'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'Saturn','value' => '12'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'Jupiter','value' => '13'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'Venus','value' => '14'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'Cylinder','value' => '21'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'Sphere','value' => '22'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'Arcade','value' => '25'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'CME','value' => '26'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'Dissipation','value' => '30'},'type' => 'e','content' => []}]},{'name' => 'parameter','attrib' => {'length' => '20','name' => 'TypeDissipation','if' => '$iProblem==30','type' => 'string'},'type' => 'e','content' => []},{'type' => 't','content' => '
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
'}]},{'name' => 'command','attrib' => {'alias' => 'COORDINATESYSTEM','multiple' => 'T','name' => 'COORDSYSTEM'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'input' => 'select','name' => 'TypeCoordSystem','type' => 'string'},'type' => 'e','content' => [{'name' => 'option','attrib' => {'default' => 'T','name' => 'GSM','if' => '$_NameComp eq \'GM\''},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'default' => 'T','name' => 'HGI','if' => '$_NameComp eq \'IH\''},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'HGC','if' => '$_NameComp eq \'IH\''},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'default' => 'T','name' => 'HGR','if' => '$_NameComp eq \'SC\''},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'HGI','if' => '$_NameComp eq \'SC\''},'type' => 'e','content' => []}]},{'type' => 't','content' => '

#COORDSYSTEM
GSM			TypeCoordSystem

! TypeCoordSystem defines the coordinate system for the component.
! Currently only one coordinate system is available for GM ("GSM")
! and two for IH ("HGI" or "HGC") and two for SC ("HGR" or "HGI"). 
! In the future "GSE" should be also an option for GM.
! The coordinate systems are defined in share/Library/src/CON_axes.
!
! Default is component dependent: "GSM" for GM, "HGI" for IH, and "HGR" for SC.
'}]},{'name' => 'command','attrib' => {'name' => 'RESTARTINDIR','if' => '$_IsFirstSession'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => 'GM/restartIN','length' => '100','name' => 'NameRestartInDir','type' => 'string'},'type' => 'e','content' => []},{'name' => 'rule','attrib' => {'expr' => '-d $NameRestartInDir'},'type' => 'e','content' => [{'type' => 't','content' => '
		Restart input directory $NameRestartInDir must exist!
	'}]},{'type' => 't','content' => '

#RESTARTINDIR
GM/restart_n5000	NameRestartInDir

! The NameRestartInDir variable contains the name of the directory
! where restart files are saved relative to the run directory.
! The directory should be inside the subdirectory with the name 
! of the component.
!
! Default value is "GM/restartIN".
'}]},{'name' => 'command','attrib' => {'name' => 'NEWRESTART','if' => '$_IsFirstSession'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'DoRestartBFace','type' => 'logical'},'type' => 'e','content' => []},{'type' => 't','content' => ' 

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

'}]},{'name' => 'command','attrib' => {'required' => 'T','name' => 'GRID','if' => '$_IsFirstSession'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => '2','name' => 'nRootBlockX','min' => '1','type' => 'integer'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '1','name' => 'nRootBlockY','min' => '1','type' => 'integer'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '1','name' => 'nRootBlockZ','min' => '1','type' => 'integer'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '-192.0','name' => 'xMin','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '  64.0','name' => 'xMax','min' => '$xMin','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => ' -64.0','name' => 'yMin','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '  64.0','name' => 'yMax','min' => '$yMin','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => ' -64.0','name' => 'zMin','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '  64.0','name' => 'zMax','min' => '$zMin','type' => 'real'},'type' => 'e','content' => []},{'type' => 't','content' => '
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
'}]},{'name' => 'command','attrib' => {'name' => 'OUTERBOUNDARY'},'type' => 'e','content' => [{'name' => 'foreach','attrib' => {'values' => 'TypeBcEast,TypeBcWest,TypeBcSouth,TypeBcNorth,TypeBcBot,TypeBcTop','name' => 'Side'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'input' => 'select','name' => '$Side','type' => 'string'},'type' => 'e','content' => [{'name' => 'option','attrib' => {'name' => 'coupled'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'default' => '$Side ne \'TypeBcEast\'','name' => 'fixed/inflow'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'default' => '$Side eq \'TypeBcEast\'','name' => 'float/outflow'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'heliofloat'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'reflect'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'periodic'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'vary'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'shear'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'linetied'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'raeder'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'arcadetop'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'arcadebot'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'arcadebotcont'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'user'},'type' => 'e','content' => []}]}]},{'name' => 'rule','attrib' => {'expr' => 'not($TypeBcEast eq \'periodic\' xor $TypeBcWest eq \'periodic\')'},'type' => 'e','content' => [{'type' => 't','content' => '
		East and west BCs must be both periodic or neither
	'}]},{'name' => 'rule','attrib' => {'expr' => 'not($TypeBcSouth eq \'periodic\' xor $TypeBcNorth eq \'periodic\')'},'type' => 'e','content' => [{'type' => 't','content' => '
		South and North BCs must be both periodic or neither
	'}]},{'name' => 'rule','attrib' => {'expr' => 'not($TypeBcBot eq \'periodic\' xor $TypeBcTop eq \'periodic\')'},'type' => 'e','content' => [{'type' => 't','content' => '
		Bottom and top BCs must be both periodic or neither
	'}]},{'type' => 't','content' => '
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
'}]},{'name' => 'command','attrib' => {'name' => 'INNERBOUNDARY'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'input' => 'select','name' => 'TypeBcInner','type' => 'string'},'type' => 'e','content' => [{'name' => 'option','attrib' => {'name' => 'reflect'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'float'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'fixed'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'default' => 'T','name' => 'ionosphere'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'ionosphereB0','value' => 'ionosphereB0/ionosphereb0'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'ionospherefloat'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'coronatoih'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'user'},'type' => 'e','content' => []}]},{'name' => 'if','attrib' => {'expr' => '$UseBody2'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'input' => 'select','name' => 'TypeBcBody2','type' => 'string'},'type' => 'e','content' => [{'name' => 'option','attrib' => {'default' => 'T','name' => 'reflect'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'float'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'fixed'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'ionosphere'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'ionosphereB0'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'ionospherefloat'},'type' => 'e','content' => []}]}]},{'name' => 'rule','attrib' => {'expr' => 'not($TypeBcBody2 =~ /ionosphere/)'},'type' => 'e','content' => [{'type' => 't','content' => ' 
	For the second body COROTATION AND AN IONOSPHERIC BOUNDARY DO NOT WORK.
	'}]},{'type' => 't','content' => '

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
'}]},{'name' => 'command','attrib' => {'name' => 'EXTRABOUNDARY'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'UseExtraBoundary','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'if','attrib' => {'expr' => '$UseExtraBoundary'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'name' => 'TypeBcExtra','type' => 'string'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'DoFixExtraboundary','type' => 'logical'},'type' => 'e','content' => []}]}]},{'name' => 'command','attrib' => {'name' => 'FACEOUTERBC'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => '0','name' => 'MaxBoundary','max' => '6','min' => '0','type' => 'integer'},'type' => 'e','content' => []},{'name' => 'if','attrib' => {'expr' => '$MaxBoundary >= 1'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'DoFixOuterBoundary','type' => 'logical'},'type' => 'e','content' => []}]},{'type' => 't','content' => '
#FACEOUTERBC
0              MaxBoundary            
F              DoFixOuterBoundary)    !read only for MaxBoundary>=East_(=1).
! If MaxBoundary is East_(=1) or more then the outer boundaries with
! the number of boundary being between East_ and MaxBoundary
! are treated using set_BCs.f90 subroutines instead of set_outerBCs.f90 
! if DoFixOuterBoundary is .true., there is no resolution
! change along the outer boundaries with the number of
! of boundary being between East_ and MaxBoundary
'}]}]},{'name' => 'commandgroup','attrib' => {'name' => 'INITIAL TIME'},'type' => 'e','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! INITIAL TIME AND STEP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'name' => 'command','attrib' => {'alias' => 'SETREALTIME','name' => 'STARTTIME','if' => '$_IsFirstSession'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => '2000','name' => 'iYear','type' => 'integer'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '3','name' => 'iMonth','max' => '12','min' => '1','type' => 'integer'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '21','name' => 'iDay','max' => '31','min' => '1','type' => 'integer'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '0','name' => 'iHour','max' => '23','min' => '0','type' => 'integer'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '0','name' => 'iMinute','max' => '59','min' => '0','type' => 'integer'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '0','name' => 'iSecond','max' => '59','min' => '0','type' => 'integer'},'type' => 'e','content' => []},{'type' => 't','content' => '
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
'}]},{'name' => 'command','attrib' => {'name' => 'TIMESIMULATION','if' => '$_IsFirstSession'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => '0.0','name' => 'tSimulation','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'type' => 't','content' => '

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
'}]},{'name' => 'command','attrib' => {'name' => 'NSTEP','if' => '$_IsFirstSession'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => '0','name' => 'nStep','min' => '0','type' => 'integer'},'type' => 'e','content' => []},{'type' => 't','content' => '

#NSTEP
100			nStep

! Set nStep for the component. Typically used in the restart.H header file.
! Generally it is not inserted in a PARAM.in file by the user.
!
! The default is nStep=0 as the starting time step with no restart.
'}]},{'name' => 'command','attrib' => {'name' => 'NPREVIOUS','if' => '$_IsFirstSession'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => '-1','name' => 'nPrevious','min' => '-1','type' => 'integer'},'type' => 'e','content' => []},{'type' => 't','content' => '

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
'}]}]},{'name' => 'commandgroup','attrib' => {'name' => 'TIME INTEGRATION'},'type' => 'e','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  TIME INTEGRATION PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'name' => 'command','attrib' => {'name' => 'TIMESTEPPING'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'input' => 'select','name' => 'nStage','type' => 'integer'},'type' => 'e','content' => [{'name' => 'option','attrib' => {'default' => 'T','value' => '1'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'value' => '2'},'type' => 'e','content' => []}]},{'name' => 'parameter','attrib' => {'default' => '0.8','name' => 'CflExpl','max' => '1','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'type' => 't','content' => '

#TIMESTEPPING
2                       nStage
0.80                    CflExpl

! Parameters for explicit time integration.
! Default is 1 stage and CflExpl=0.8
'}]},{'name' => 'command','attrib' => {'name' => 'FIXEDTIMESTEP'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'UseDtFixed','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '1.0','name' => 'DtFixedDim','if' => '$UseDtFixed','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'type' => 't','content' => '
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

'}]},{'name' => 'command','attrib' => {'name' => 'PARTSTEADY'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'UsePartSteady','type' => 'logical'},'type' => 'e','content' => []},{'type' => 't','content' => '

! If UsePartSteady is true, the partially steady state algorithm is used.
! Only blocks which are changing or next to changing blocks are evolved.
! This can speed up the calculation if part of the domain is in a 
! numerical steady state. Default value is .false.
'}]},{'name' => 'command','attrib' => {'name' => 'PARTLOCAL'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'UsePartLocal','type' => 'logical'},'type' => 'e','content' => []},{'type' => 't','content' => '
#PARTLOCAL
T               UsePartLocal

! Default is UsePartLocal=.false. If UsePartLocal is true and the
! run is time accurate, then the blocks selected as "implicit"
! by the criteria defined in #STEPPINGCRITERIA are not used to
! calculate the time step, and all cells are advanced with the
! smaller of the stable and the global time steps.
'}]},{'name' => 'command','attrib' => {'name' => 'IMPLICIT'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'UsePointImplicit','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'UsePartImplicit','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'UseFullImplicit','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '100','name' => 'CflImpl','if' => '$UsePartImplicit or $UseFullImplicit','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'rule','attrib' => {'expr' => '$UsePointImplicit + $UsePartImplicit + $UseFullImplicit <= 1'},'type' => 'e','content' => [{'type' => 't','content' => '
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
'}]}]},{'name' => 'commandgroup','attrib' => {'name' => 'IMPLICIT PARAMETERS'},'type' => 'e','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! PARAMETERS FOR FULL AND PART IMPLICIT TIME INTEGRATION !!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'name' => 'command','attrib' => {'alias' => 'STEPPINGCRITERIA','name' => 'IMPLICITCRITERIA'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'input' => 'select','name' => 'TypeImplCrit','type' => 'string'},'type' => 'e','content' => [{'name' => 'option','attrib' => {'default' => 'T','name' => 'Time step','value' => 'dt'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'Radial distance','value' => 'r/R'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'Test block','value' => 'test'},'type' => 'e','content' => []}]},{'name' => 'parameter','attrib' => {'name' => 'rImplicit','if' => '$TypeImplCrit eq \'R\'','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'type' => 't','content' => '

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
'}]},{'name' => 'command','attrib' => {'alias' => 'IMPLICITSTEP','name' => 'IMPLSTEP'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => '1','name' => 'ImplCoeff','max' => '1','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => 'T','name' => 'UseBdf2','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'UseSourceImpl','type' => 'logical'},'type' => 'e','content' => []},{'type' => 't','content' => '
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
'}]},{'name' => 'command','attrib' => {'alias' => 'IMPLICITSCHEME','name' => 'IMPLSCHEME'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'input' => 'select','name' => 'nOrderImpl','type' => 'integer'},'type' => 'e','content' => [{'name' => 'option','attrib' => {'default' => 'T','name' => '1'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => '2'},'type' => 'e','content' => []}]},{'name' => 'parameter','attrib' => {'input' => 'select','name' => 'TypeFluxImpl','type' => 'string'},'type' => 'e','content' => [{'name' => 'option','attrib' => {'name' => 'Roe','value' => 'Roe/roe/1'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'default' => 'T','name' => 'Rusanov','value' => 'Rusanov/rusanov/2/TVDLF'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'Linde','value' => 'Linde/linde/3/HLLEL'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'Sokolov','value' => 'Sokolov/sokolov/4/AW'},'type' => 'e','content' => []}]},{'type' => 't','content' => '
#IMPLSCHEME
1               nOrderImpl
Rusanov         TypeFluxImpl

! This command defines the scheme used in the implicit part (\'left hand side\').
! The default order is first order. The default scheme is the same as the
! scheme selected for the explicit part. 
'}]},{'name' => 'command','attrib' => {'alias' => 'IMPLICITCHECK','name' => 'IMPLCHECK'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => '0.3','name' => 'RejectStepLevel','max' => '0.9','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '0.5','name' => 'RejectStepFactor','max' => '0.9','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '0.6','name' => 'ReduceStepLevel','max' => '0.9','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '0.9','name' => 'ReduceStepFactor','max' => '1','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '0.8','name' => 'IncreaseStepLevel','max' => '1','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '1.05','name' => 'IncreaseStepFactor','max' => '2','min' => '1','type' => 'real'},'type' => 'e','content' => []},{'type' => 't','content' => '

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
'}]},{'name' => 'command','attrib' => {'name' => 'NEWTON'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'UseConservativeImplicit','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'UseNewton','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'if','attrib' => {'expr' => '$UseNewton'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => 'T','name' => 'UseNewMatrix','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '10','name' => 'MaxIterNewton','min' => '1','type' => 'integer'},'type' => 'e','content' => []}]},{'type' => 't','content' => '
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
'}]},{'name' => 'command','attrib' => {'name' => 'JACOBIAN'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'input' => 'select','name' => 'TypeJacobian','type' => 'string'},'type' => 'e','content' => [{'name' => 'option','attrib' => {'default' => 'T','name' => 'Preconditioned','value' => 'prec'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'No preconditioning','value' => 'free'},'type' => 'e','content' => []}]},{'name' => 'parameter','attrib' => {'default' => '$doublePrecision ? 1.e-12 : 1.e-6','name' => 'JacobianEps','max' => '1.e-5','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'type' => 't','content' => '
#JACOBIAN
prec            TypeJacobian (prec, free)
1.E-12          JacobianEps

! The Jacobian matrix is always calculated with a matrix free approach,
! however it can be preconditioned  (\'prec\'), or not (\'free\').  The
! Default value is TypeJacobian=\'prec\'.
! JacobianEps contains the machine round off error for numerical derivatives.
! The default value is 1.E-12 for 8 byte reals and 1.E-6 for 4 byte reals.
'}]},{'name' => 'command','attrib' => {'name' => 'PRECONDITIONER'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'input' => 'select','name' => 'TypePrecondSide','type' => 'string'},'type' => 'e','content' => [{'name' => 'option','attrib' => {'name' => 'left'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'default' => 'T','name' => 'symmetric'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'right'},'type' => 'e','content' => []}]},{'name' => 'parameter','attrib' => {'input' => 'select','name' => 'TypePrecond','type' => 'string'},'type' => 'e','content' => [{'name' => 'option','attrib' => {'default' => 'T','name' => 'MBILU'},'type' => 'e','content' => []}]},{'name' => 'parameter','attrib' => {'default' => '0.5','name' => 'GustafssonPar','max' => '1','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'type' => 't','content' => '
#PRECONDITIONER
symmetric       TypePrecondSide (left, symmetric, right)
MBILU           TypePrecond (MBILU)
0.5             GustafssonPar (0. no modification, 1. full modification)

! Default parameters are shown. Right preconditioning does not affect
! the normalization of the residual. The Gustafsson parameter determines
! how much the MBILU preconditioner is modified. The default 0.5 value
! means a relaxed modification.
'}]},{'name' => 'command','attrib' => {'name' => 'KRYLOV'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'input' => 'select','name' => 'TypeKrylov','type' => 'string'},'type' => 'e','content' => [{'name' => 'option','attrib' => {'default' => 'T','name' => 'gmres'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'bicgstab'},'type' => 'e','content' => []}]},{'name' => 'parameter','attrib' => {'input' => 'select','name' => 'TypeInitKrylov','type' => 'string'},'type' => 'e','content' => [{'name' => 'option','attrib' => {'default' => 'T','name' => '0','value' => 'nul'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'previous','value' => 'old'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'explicit'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'scaled explicit','value' => 'explicit'},'type' => 'e','content' => []}]},{'name' => 'parameter','attrib' => {'default' => '0.001','name' => 'ErrorMaxKrylov','max' => '0.1','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '100','name' => 'MaxMatvecKrylov','min' => '1','type' => 'integer'},'type' => 'e','content' => []},{'type' => 't','content' => '
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
'}]},{'name' => 'command','attrib' => {'name' => 'KRYLOVSIZE'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => 'MaxMatvecKrylov','name' => 'nKrylovVector','min' => '1','type' => 'integer'},'type' => 'e','content' => []},{'type' => 't','content' => '
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
'}]}]},{'name' => 'commandgroup','attrib' => {'name' => 'STOPPING CRITERIA'},'type' => 'e','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! STOPPING CRITERIA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

The commands in this group only work in stand alone mode.

'},{'name' => 'command','attrib' => {'required' => '$_IsStandAlone','name' => 'STOP','if' => '$_IsStandAlone'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => '-1','name' => 'MaxIteration','min' => '-1','type' => 'integer'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '-1','name' => 'tSimulationMax','min' => '-1','type' => 'real'},'type' => 'e','content' => []},{'type' => 't','content' => '

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
'}]},{'name' => 'command','attrib' => {'name' => 'CHECKSTOPFILE','if' => '$_IsStandAlone'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => 'T','name' => 'DoCheckStopFile','type' => 'logical'},'type' => 'e','content' => []},{'type' => 't','content' => '

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
'}]},{'name' => 'command','attrib' => {'name' => 'CPUTIMEMAX','if' => '$_IsStandAlone'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => '-1','name' => 'CpuTimeMax','min' => '-1','type' => 'real'},'type' => 'e','content' => []},{'type' => 't','content' => '

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
'}]}]},{'name' => 'commandgroup','attrib' => {'name' => 'OUTPUT PARAMETERS'},'type' => 'e','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  OUTPUT PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'name' => 'command','attrib' => {'name' => 'RESTARTOUTDIR'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => 'GM/restartOUT','length' => '100','name' => 'NameRestartOutDir','type' => 'string'},'type' => 'e','content' => []},{'name' => 'rule','attrib' => {'expr' => '-d $NameRestartOutDir'},'type' => 'e','content' => [{'type' => 't','content' => '
		Restart output directory $NameRestartOutDir must exist
	'}]},{'type' => 't','content' => '

#RESTARTOUTDIR
GM/restart_n5000	NameRestartOutDir

! The NameRestartOutDir variable contains the name of the directory
! where restart files are saved relative to the run directory.
! The directory should be inside the subdirectory with the name 
! of the component.
!
! Default value is "GM/restartOUT".
'}]},{'name' => 'command','attrib' => {'name' => 'SAVERESTART'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => 'T','name' => 'DoSaveRestart','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'if','attrib' => {'expr' => '$DoSaveRestart'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => '-1','name' => 'DnSaveRestart','min' => '-1','type' => 'integer'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '-1','name' => 'DtSaveRestart','min' => '-1','type' => 'real'},'type' => 'e','content' => []}]},{'type' => 't','content' => '
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
'}]},{'name' => 'command','attrib' => {'name' => 'PLOTDIR'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => 'GM/IO2','length' => '100','name' => 'NamePlotDir','type' => 'string'},'type' => 'e','content' => []},{'name' => 'rule','attrib' => {'expr' => '-d $NamePlotDir'},'type' => 'e','content' => [{'type' => 't','content' => '
		Plot directory $NamePlotDir must exist
	'}]},{'type' => 't','content' => '

The NamePlotDir variable contains the name of the directory
where plot files and logfiles are saved relative to the run directory.
The directory should be inside the subdirectory with the name
of the component.

Default value is "GM/IO2".
'}]},{'name' => 'command','attrib' => {'name' => 'SAVELOGFILE'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'DoSaveLogfile','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'if','attrib' => {'expr' => '$DoSaveLogfile'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'name' => 'StringLog','max' => '4','min' => '1','type' => 'strings'},'type' => 'e','content' => [{'name' => 'part','attrib' => {'required' => 'T','input' => 'select','name' => 'TypeLogVar','type' => 'string'},'type' => 'e','content' => [{'name' => 'option','attrib' => {'default' => 'T','name' => 'MHD vars. dimensional','value' => 'MHD'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'Flux vars. dimensional','value' => 'FLX'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'Raw vars. dimensional','value' => 'RAW'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'Set vars. dimensional','value' => 'VAR'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'default' => 'T','name' => 'MHD vars. scaled','value' => 'mhd'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'Flux vars. scaled','value' => 'flx'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'Raw vars. scaled','value' => 'raw'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'Set vars. scaled','value' => 'var'},'type' => 'e','content' => []}]},{'name' => 'part','attrib' => {'required' => 'F','input' => 'select','multiple' => 'T','name' => 'TypeTime','type' => 'string'},'type' => 'e','content' => [{'name' => 'option','attrib' => {'exclusive' => 'T','name' => 'none'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'step'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'date'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'time'},'type' => 'e','content' => []}]}]},{'name' => 'parameter','attrib' => {'default' => '1','name' => 'DnSaveLogfile','min' => '-1','type' => 'integer'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '-1','name' => 'DtSaveLogfile','min' => '-1','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'length' => '100','name' => 'NameLogVars','if' => '$TypeLogVar =~ /var/i','type' => 'string'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'length' => '100','name' => 'StringLogRadii','if' => '($TypeLogVar=~/flx/i or $NameLogVars=~/flx/i)','max' => '10','min' => '1','type' => 'strings'},'type' => 'e','content' => [{'name' => 'part','attrib' => {'multiple' => 'T','name' => 'LogRadii','min' => '$rBody','type' => 'real'},'type' => 'e','content' => []}]}]},{'type' => 't','content' => '
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
'}]},{'name' => 'command','attrib' => {'name' => 'SATELLITE','if' => '$_IsFirstSession'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => '0','name' => 'nSatellite','min' => '0','type' => 'integer'},'type' => 'e','content' => []},{'name' => 'for','attrib' => {'from' => '1','to' => '$nSatellite'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'name' => 'StringSatellite','max' => '5','min' => '1','type' => 'strings'},'type' => 'e','content' => [{'name' => 'part','attrib' => {'required' => 'T','input' => 'select','name' => 'TypeSatelliteVar','type' => 'string'},'type' => 'e','content' => [{'name' => 'option','attrib' => {'default' => 'T','name' => 'MHD vars. dimensional','value' => 'MHD'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'All vars. dimensional','value' => 'FUL'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'Set vars. dimensional','value' => 'VAR'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'MHD vars. scaled','value' => 'mhd'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'All vars. scaled','value' => 'ful'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'Set vars. scaled','value' => 'var'},'type' => 'e','content' => []}]},{'name' => 'part','attrib' => {'required' => 'F','input' => 'select','name' => 'TypeTrajectory','type' => 'string'},'type' => 'e','content' => [{'name' => 'option','attrib' => {'default' => 'T','name' => 'file'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'equation','value' => 'eqn'},'type' => 'e','content' => []}]},{'name' => 'part','attrib' => {'required' => 'F','input' => 'select','multiple' => 'T','name' => 'TypeTime','type' => 'string'},'type' => 'e','content' => [{'name' => 'option','attrib' => {'exclusive' => 'T','name' => 'none'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'step'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'date'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'time'},'type' => 'e','content' => []}]}]},{'name' => 'parameter','attrib' => {'default' => '1','name' => 'DnOutput','min' => '-1','type' => 'integer'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '-1','name' => 'DtOutput','min' => '-1','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'length' => '100','name' => 'NameTrajectoryFile','type' => 'string'},'type' => 'e','content' => []},{'name' => 'rule','attrib' => {'expr' => '-f $NameTrajectoryFile'},'type' => 'e','content' => [{'type' => 't','content' => '
			Trajectory file $NameTrajectoryFile must exist
		'}]},{'name' => 'parameter','attrib' => {'length' => '100','name' => 'NameSatelliteVars','if' => '$TypeSatelliteVar =~ /\\bvar\\b/i','type' => 'string'},'type' => 'e','content' => []}]},{'type' => 't','content' => '
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


'}]},{'name' => 'command','attrib' => {'name' => 'SAVEPLOT'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => '0','name' => 'nPlotFile','max' => '100','min' => '0','type' => 'integer'},'type' => 'e','content' => []},{'name' => 'for','attrib' => {'from' => '1','to' => '$nPlotFile'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'name' => 'StringPlot','max' => '3','min' => '3','type' => 'strings'},'type' => 'e','content' => [{'name' => 'part','attrib' => {'required' => 'T','input' => 'select','name' => 'plotform','type' => 'string'},'type' => 'e','content' => [{'name' => 'option','attrib' => {'name' => 'TECPLOT','value' => 'tec'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'IDL','value' => 'idl'},'type' => 'e','content' => []}]},{'name' => 'part','attrib' => {'required' => 'T','input' => 'select','name' => 'plotarea','type' => 'string'},'type' => 'e','content' => [{'name' => 'option','attrib' => {'name' => '3D','value' => '3d/3d_'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'value' => 'x=0'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'default' => 'T','value' => 'y=0'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'value' => 'z=0'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'value' => 'sph'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'value' => 'los'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'value' => 'lin'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'if' => '$plotform =~ /\\bidl\\b/','value' => 'cut'},'type' => 'e','content' => []}]},{'name' => 'part','attrib' => {'required' => 'T','input' => 'select','name' => 'plotvar','type' => 'string'},'type' => 'e','content' => [{'name' => 'option','attrib' => {'name' => 'MHD vars. dimensional','value' => 'MHD'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'All vars. dimensional','value' => 'FUL'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'Raw vars. dimensional','value' => 'RAW'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'Ray tracing vars. dim.','value' => 'RAY'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'Flux vars. dimensional','value' => 'FLX'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'Solar vars. dimensional','value' => 'SOL'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'Position vars. dimensional','if' => '$plotarea eq \'lin\'','value' => 'POS'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'Select dimensional vars.','value' => 'VAR'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'MHD vars. scaled','value' => 'mhd'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'All vars. scaled','value' => 'ful'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'Raw vars. scaled','value' => 'raw'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'Ray tracing vars. scaled','value' => 'ray'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'Flux vars. scaled','value' => 'flx'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'Solar vars. scaled','value' => 'sol'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'Position vars. scaled','if' => '$plotarea eq \'lin\'','value' => 'pos'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'Select scaled vars.','value' => 'var'},'type' => 'e','content' => []}]}]},{'name' => 'parameter','attrib' => {'name' => 'DnSavePlot','min' => '-1','type' => 'integer'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'name' => 'DtSavePlot','min' => '-1','type' => 'real'},'type' => 'e','content' => []},{'name' => 'if','attrib' => {'expr' => '$plotarea =~ /\\bcut\\b/'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'name' => 'xMinCut','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'name' => 'xMaxCut','min' => '$xMinCut','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'name' => 'yMinCut','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'name' => 'yMaxCut','min' => '$yMinCut','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'name' => 'zMinCut','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'name' => 'zMaxCut','min' => '$zMinCut','type' => 'real'},'type' => 'e','content' => []}]},{'name' => 'parameter','attrib' => {'default' => '10','name' => 'Radius','if' => '$plotarea =~ /\\bsph\\b/','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'if','attrib' => {'expr' => '$plotarea =~ /\\blos\\b/'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => '0','name' => 'LosVectorX','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '0.0001','name' => 'LosVectorY','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '1','name' => 'LosVectorZ','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '20','name' => 'xSizeImage','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '20','name' => 'ySizeImage','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '10','name' => 'xOffset','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '10','name' => 'yOffset','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '2.5','name' => 'rOccult','min' => '1','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '0.5','name' => 'MuLimbDarkening','max' => '1','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '200','name' => 'nPixX','min' => '2','type' => 'integer'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '200','name' => 'nPixY','min' => '2','type' => 'integer'},'type' => 'e','content' => []}]},{'name' => 'if','attrib' => {'expr' => '$plotarea =~ /\\blin\\b/'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'input' => 'select','name' => 'NameLine','type' => 'string'},'type' => 'e','content' => [{'name' => 'option','attrib' => {'name' => 'Advected B','value' => 'A'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'default' => 'T','name' => 'B'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'U'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'J'},'type' => 'e','content' => []}]},{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'IsSingleLine','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '1','name' => 'nLine','max' => '20','min' => '1','type' => 'integer'},'type' => 'e','content' => []},{'name' => 'for','attrib' => {'from' => '1','to' => '$nLine'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'name' => 'xStartLine','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'name' => 'yStartLine','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'name' => 'zStartLine','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'name' => 'IsParallel','type' => 'logical'},'type' => 'e','content' => []}]}]},{'name' => 'parameter','attrib' => {'default' => '-1.0','name' => 'DxSavePlot','if' => '($plotform=~/\\bidl\\b/ and $plotarea!~/\\b(sph|los|lin)\\b/)','min' => '-1.0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'if','attrib' => {'expr' => '$plotvar =~ /\\bvar\\b/i'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'length' => '100','name' => 'NameVars','type' => 'string'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'length' => '100','name' => 'NamePars','type' => 'string'},'type' => 'e','content' => []}]}]},{'type' => 't','content' => '
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

! Default is nPlotFile=0. \\\\
!
! \\noindent
! StringPlot must contain the following 3 parts in arbitrary order
!\\begin{verbatim}
! plotarea plotvar plotform
!
! plotarea = \'3d\' , \'x=0\', \'y=0\', \'z=0\', \'cut\', \'sph\', \'los\', \'lin\'
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
! PlotRange		if plotarea is \'cut\'
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
! For Tecplot (tec) file, plot range is read but only 1 dimension is used.  
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

'}]},{'name' => 'command','attrib' => {'name' => 'SAVEBINARY'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => 'T','name' => 'DoSaveBinary','type' => 'logical'},'type' => 'e','content' => []},{'type' => 't','content' => '
#SAVEBINARY
T			DoSaveBinary   used only for \'idl\' plot file

! Default is .true. Saves unformatted IO2/*.idl files if true. 
! This is the recommended method, because it is fast and accurate.
! The only advantage of saving IO2/*.idl in formatted text files is
! that it can be processed on another machine or with a different 
! (lower) precision. For example PostIDL.exe may be compiled with 
! single precision to make IO2/*.out files smaller, while BATSRUS.exe is 
! compiled in double precision, to make results more accurate.
'}]},{'name' => 'command','attrib' => {'name' => 'SAVEPLOTSAMR'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'DoSavePlotsAmr','type' => 'logical'},'type' => 'e','content' => []},{'type' => 't','content' => '
#SAVEPLOTSAMR
F			DoSavePlotsAmr

! Save plots before each AMR. Default is DoSavePlotsAMR=.false.
'}]}]},{'name' => 'commandgroup','attrib' => {'name' => 'AMR PARAMETERS'},'type' => 'e','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  AMR PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'name' => 'command','attrib' => {'name' => 'AMRINIT','if' => '$_IsFirstSession'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'input' => 'select','name' => 'InitialRefineType','type' => 'string'},'type' => 'e','content' => [{'name' => 'option','attrib' => {'default' => '1','name' => 'default'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'all'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'none'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => '3Dbodyfocus'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'spherefocus'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'magnetosphere'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'points'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'coupledhelio'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'helio_init'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'helio_z=4'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'all_then_focus'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'cme'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'points'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'mag_new'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'magnetosphere'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'magneto_fine'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'magneto12'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'magnetosaturn'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'magnetojupiter'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'paleo'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'comet'},'type' => 'e','content' => []}]},{'name' => 'parameter','attrib' => {'default' => '4','name' => 'InitialRefineLevel','min' => '0','type' => 'integer'},'type' => 'e','content' => []},{'type' => 't','content' => '
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
'}]},{'name' => 'command','attrib' => {'name' => 'AMRINITPHYSICS','if' => '$_IsFirstSession'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => '0','name' => 'nRefineLevelIC','min' => '0','type' => 'integer'},'type' => 'e','content' => []},{'type' => 't','content' => '
#AMRINITPHYSICS
3			nRefineLevelIC

! Defines number of physics (initial condition) based AMR-s AFTER the 
! geometry based initial AMR-s defined by #AMRINIT were done.
! Only useful if the initial condition has a non-trivial analytic form.
'}]},{'name' => 'command','attrib' => {'alias' => 'GRIDLEVEL','multiple' => 'T','name' => 'GRIDRESOLUTION'},'type' => 'e','content' => [{'name' => 'set','attrib' => {'name' => 'RotateArea','value' => '0'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'name' => 'Resolution','if' => '$_command =~ /RESOLUTION/','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'name' => 'nLevel','if' => '$_command =~ /LEVEL/','min' => '0','type' => 'integer'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'name' => 'StringArea','max' => '2','min' => '1','type' => 'strings'},'type' => 'e','content' => [{'name' => 'part','attrib' => {'required' => 'T','input' => 'select','name' => 'NameArea','type' => 'string'},'type' => 'e','content' => [{'name' => 'option','attrib' => {'default' => 'T','name' => 'initial','value' => 'init/initial'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'all'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'box'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'brick'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'brick0'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'sphere'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'sphere0'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'shell'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'shell0'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'cylinderx'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'cylinderx0'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'cylindery'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'cylindery0'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'cylinderz'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'cylinderz0'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'ringx'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'ringx0'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'ringy'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'ringy0'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'ringz'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'ringz0'},'type' => 'e','content' => []}]},{'name' => 'part','attrib' => {'required' => 'F','input' => 'select','name' => 'RotateArea','type' => 'string'},'type' => 'e','content' => [{'name' => 'option','attrib' => {'name' => 'rotated'},'type' => 'e','content' => []}]}]},{'name' => 'if','attrib' => {'expr' => '$NameArea !~ /box|all|init|0/'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => '0','name' => 'xCenter','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '0','name' => 'yCenter','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '0','name' => 'zCenter','type' => 'real'},'type' => 'e','content' => []}]},{'name' => 'if','attrib' => {'expr' => '$NameArea =~ /box/'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'name' => 'xMinBox','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'name' => 'yMinBox','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'name' => 'zMinBox','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'name' => 'xMaxBox','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'name' => 'yMaxBox','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'name' => 'zMaxBox','type' => 'real'},'type' => 'e','content' => []}]},{'name' => 'if','attrib' => {'expr' => '$NameArea =~ /brick/i'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'name' => 'xSizeBrick','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'name' => 'ySizeBrick','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'name' => 'zSizeBrick','min' => '0','type' => 'real'},'type' => 'e','content' => []}]},{'name' => 'if','attrib' => {'expr' => '$NameArea =~ /sphere/i'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'name' => 'Radius','min' => '0','type' => 'real'},'type' => 'e','content' => []}]},{'name' => 'if','attrib' => {'expr' => '$NameArea =~ /shell/i'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'name' => 'Radius1','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'name' => 'Radius2','min' => '0','type' => 'real'},'type' => 'e','content' => []}]},{'name' => 'if','attrib' => {'expr' => '$NameArea =~ /cylinder/i'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'name' => 'LengthCylinder','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'name' => 'Radius','min' => '0','type' => 'real'},'type' => 'e','content' => []}]},{'name' => 'if','attrib' => {'expr' => '$NameArea =~ /ring/i'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'name' => 'HeightRing','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'name' => 'Radius1','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'name' => 'Radius2','min' => '0','type' => 'real'},'type' => 'e','content' => []}]},{'name' => 'if','attrib' => {'expr' => '$RotateArea =~ /rotated/i'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => '0','name' => 'xRotate','max' => '360','min' => '-360','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '0','name' => 'yRotate','max' => '360','min' => '-360','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '0','name' => 'zRotate','max' => '360','min' => '-360','type' => 'real'},'type' => 'e','content' => []}]},{'type' => 't','content' => '

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
'}]},{'name' => 'command','attrib' => {'name' => 'AMRLEVELS'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => '0','name' => 'MinBlockLevel','min' => '-1','type' => 'integer'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '99','name' => 'MaxBlockLevel','min' => '-1','type' => 'integer'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'DoFixBodyLevel','type' => 'logical'},'type' => 'e','content' => []},{'type' => 't','content' => '
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
'}]},{'name' => 'command','attrib' => {'name' => 'AMRRESOLUTION'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => '0','name' => 'DxCellMin','min' => '-1','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '99999','name' => 'DxCellMax','min' => '-1','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'DoFixBodyLevel','type' => 'logical'},'type' => 'e','content' => []},{'type' => 't','content' => '
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
'}]},{'name' => 'command','attrib' => {'name' => 'AMR'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => '-1','name' => 'DnRefine','min' => '-1','type' => 'integer'},'type' => 'e','content' => []},{'name' => 'if','attrib' => {'expr' => '$DnRefine>0'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'DoAutoRefine','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'if','attrib' => {'expr' => '$DoAutoRefine'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => '20','name' => 'PercentCoarsen','max' => '100','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '20','name' => 'PercentRefine','max' => '100','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '99999','name' => 'MaxTotalBlocks','min' => '1','type' => 'integer'},'type' => 'e','content' => []}]}]},{'type' => 't','content' => '
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
'}]},{'name' => 'command','attrib' => {'name' => 'AMRCRITERIA'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'input' => 'select','name' => 'nRefineCrit','type' => 'integer'},'type' => 'e','content' => [{'name' => 'option','attrib' => {'name' => '1'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => '2'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'default' => '1','name' => '3'},'type' => 'e','content' => []}]},{'name' => 'for','attrib' => {'from' => '1','to' => '$nRefineCrit'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'input' => 'select','name' => 'TypeRefine','type' => 'string'},'type' => 'e','content' => [{'name' => 'option','attrib' => {'name' => 'grad T','value' => 'gradt/gradT'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'grad P','value' => 'gradp/gradP'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'grad log(Rho)','value' => 'gradlogrho'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'grad log(p)','value' => 'gradlogP/gradlogp'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'grad E','value' => 'gradE'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'curl U','value' => 'curlV/curlv/curlU/curlu'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'curl B','value' => 'curlB/curlb'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'div U','value' => 'divU/divu/divV/divv'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'divB','value' => 'divb/divB'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'vAlfven','value' => 'Valfven/vAlfven/valfven'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'heliospheric beta','value' => 'heliobeta'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'flux'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'heliospheric current sheet','value' => 'heliocurrentsheet'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'rCurrents','value' => 'rcurrents/Rcurrents'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'Transient','value' => 'transient/Transient'},'type' => 'e','content' => []}]},{'name' => 'if','attrib' => {'expr' => '$TypeRefine =~ /transient/i'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'input' => 'select','name' => 'TypeTransient','type' => 'string'},'type' => 'e','content' => [{'name' => 'option','attrib' => {'name' => 'P_dot','value' => 'p_dot/P_dot'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'T_dot','value' => 't_dot/T_dot'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'default' => 'T','name' => 'Rho_dot','value' => 'rho_dot/Rho_dot'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'RhoU_dot','value' => 'RhoU_dot/rhou_dot'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'Rho_2nd_1','value' => 'Rho_2nd_1/rho_2nd_1'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'Rho_2nd_2','value' => 'Rho_2nd_2/rho_2nd_2'},'type' => 'e','content' => []}]},{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'UseSunEarth','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'if','attrib' => {'expr' => '$UseSunEarth'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'name' => 'xEarth','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'name' => 'yEarth','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'name' => 'zEarth','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'name' => 'InvD2Ray','type' => 'real'},'type' => 'e','content' => []}]}]}]},{'type' => 't','content' => '
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
'}]}]},{'name' => 'commandgroup','attrib' => {'name' => 'SCHEME PARAMETERS'},'type' => 'e','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  SCHEME PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'name' => 'command','attrib' => {'name' => 'SCHEME'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'input' => 'select','name' => 'nOrder','type' => 'integer'},'type' => 'e','content' => [{'name' => 'option','attrib' => {'default' => 'T','name' => '1'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => '2'},'type' => 'e','content' => []}]},{'name' => 'parameter','attrib' => {'input' => 'select','name' => 'TypeFlux','type' => 'string'},'type' => 'e','content' => [{'name' => 'option','attrib' => {'name' => 'Roe','value' => 'Roe/roe/1'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'default' => 'T','name' => 'Rusanov','value' => 'Rusanov/rusanov/2/TVDLF'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'Linde','value' => 'Linde/linde/3/HLLEL'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'Sokolov','value' => 'Sokolov/sokolov/4/AW'},'type' => 'e','content' => []}]},{'name' => 'if','attrib' => {'expr' => '$nOrder == 2'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'input' => 'select','name' => 'TypeLimiter','type' => 'string'},'type' => 'e','content' => [{'name' => 'option','attrib' => {'default' => 'T','name' => 'minmod'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'mc'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'beta'},'type' => 'e','content' => []}]},{'name' => 'parameter','attrib' => {'default' => '1.2','name' => 'LimiterBeta','if' => '$TypeLimiter ne \'minmod\'','max' => '2','min' => '1','type' => 'real'},'type' => 'e','content' => []}]},{'type' => 't','content' => '
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
'}]},{'name' => 'command','attrib' => {'name' => 'NONCONSERVATIVE'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => 'T','name' => 'UseNonConservative','type' => 'logical'},'type' => 'e','content' => []},{'type' => 't','content' => '
#NONCONSERVATIVE
T		UseNonConservative

! For Earth the default is using non-conservative equations 
! (close to the body).
'}]},{'name' => 'command','attrib' => {'name' => 'CONSERVATIVECRITERIA'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => '1','name' => 'nConservCrit','max' => '3','min' => '0','type' => 'integer'},'type' => 'e','content' => []},{'name' => 'for','attrib' => {'from' => '1','to' => '$nConservCrit'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'input' => 'select','name' => 'TypeConservCrit','type' => 'string'},'type' => 'e','content' => [{'name' => 'option','attrib' => {'default' => 'T','name' => 'radius','value' => 'r/R/radius/Radius'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'parabola','value' => 'parabola/paraboloid'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'p','value' => 'p/P'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'grad P','value' => 'gradp/GradP'},'type' => 'e','content' => []}]},{'name' => 'parameter','attrib' => {'default' => '6','name' => 'rConserv','if' => '$TypeConservCrit =~ /^r|radius$/i','min' => '$rBody','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '6','name' => 'xParabolaConserv','if' => '$TypeConservCrit =~ /^parabol/i','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '36','name' => 'yParabolaConserv','if' => '$TypeConservCrit =~ /^parabol/i','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '0.05','name' => 'pCoeffConserv','if' => '$TypeConservCrit =~ /^p$/i','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '0.1','name' => 'GradPCoeffConserv','if' => '$TypeConservCrit =~ /gradp/i','min' => '0','type' => 'real'},'type' => 'e','content' => []}]},{'type' => 't','content' => '

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
'}]},{'name' => 'command','attrib' => {'name' => 'UPDATECHECK'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => 'T','name' => 'UseUpdateCheck','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'if','attrib' => {'expr' => '$UseUpdateCheck'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => '40','name' => 'RhoMinPercent','max' => '100','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '400','name' => 'RhoMaxPercent','min' => '100','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '40','name' => 'pMinPercent','max' => '100','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '400','name' => 'pMaxPercent','min' => '100','type' => 'real'},'type' => 'e','content' => []}]},{'type' => 't','content' => '
#UPDATECHECK
T			UseUpdateCheck
40.			RhoMinPercent
400.			RhoMaxPercent
40.			pMinPercent
400.			pMaxPercent

! Default values are shown.  This will adjust the timestep so that
! density and pressure cannot change by more than the given percentages
! in a single timestep.
'}]},{'name' => 'command','attrib' => {'name' => 'PROLONGATION'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'input' => 'select','name' => 'nOrderProlong','type' => 'integer'},'type' => 'e','content' => [{'name' => 'option','attrib' => {'default' => 'T','name' => '1'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => '2'},'type' => 'e','content' => []}]},{'name' => 'parameter','attrib' => {'input' => 'select','name' => 'TypeProlong','if' => '$nOrderProlong==2','type' => 'string'},'type' => 'e','content' => [{'name' => 'option','attrib' => {'default' => 'T','name' => 'left-right','value' => 'lr'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'central'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'minmod'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'left-right extrapolate','value' => 'lr2'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'central    extrapolate','value' => 'central2'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'minmod     extrapolate','value' => 'minmod2'},'type' => 'e','content' => []}]},{'type' => 't','content' => '
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
'}]},{'name' => 'command','attrib' => {'alias' => 'OPTIMIZE','name' => 'MESSAGEPASS'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'input' => 'select','name' => 'TypeMessagePass','type' => 'string'},'type' => 'e','content' => [{'name' => 'option','attrib' => {'default' => 'T','name' => 'm_p_cell FACES ONLY','value' => 'allopt'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'm_p_cell','value' => 'all'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'm_p_dir FACES ONLY','value' => 'opt'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'm_p_dir group by directions','value' => 'dir'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'm_p_dir group by faces     ','value' => 'face'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'm_p_dir group by kind and face','value' => 'min'},'type' => 'e','content' => []}]},{'type' => 't','content' => '
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
'}]},{'name' => 'command','attrib' => {'name' => 'TVDRESCHANGE'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'UseTvdAtReschange','type' => 'logical'},'type' => 'e','content' => []},{'type' => 't','content' => '
#TVDRESCHANGE
T		UseTvdAtResChange

! For UseTvdAtResChange=T a second order TVD limited scheme is used 
! at the resolution changes. 
!
! Default value is false, which results in first order 
! prolongation and restriction operators at the resolution changes.
'}]},{'name' => 'command','attrib' => {'name' => 'BORIS'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'UseBorisCorrection','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '1','name' => 'BorisClightFactor','if' => '$UseBorisCorrection','max' => '1','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'type' => 't','content' => '
#BORIS
T			UseBorisCorrection
1.0			BorisClightFactor !Only if UseBorisCorrection is true

! Default is boris_correction=.false.
! Use semi-relativistic MHD equations with speed of light reduced by
! the BorisClightFactor. Set BorisClightFactor=1.0 for true semi-relativistic
! MHD. Gives the same steady state as normal MHD analytically, but there
! can be differences due to discretization errors. 
! You can use either Boris or BorisSimple but not both. 
'}]},{'name' => 'command','attrib' => {'name' => 'BORISSIMPLE'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'UseBorisSimple','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '1','name' => 'BorisClightFactor','if' => '$UseBorisSimple','max' => '1','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'type' => 't','content' => '
#BORISSIMPLE
T			UseBorisSimple
0.05			BorisClightFactor !Only if UseBorisSimple is true

! Default is UseBorisSimple=.false. 
! Use simplified semi-relativistic MHD with speed of light reduced by the
! BorisClightFactor. This is only useful with BorisClightFactor less than 1.
! Should give the same steady state as normal MHD, but there can be a
! difference due to discretization errors.
! You can use either Boris or BorisSimple but not both. 
'}]},{'name' => 'command','attrib' => {'name' => 'DIVB'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => 'T','name' => 'UseDivbSource','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'UseDivbDiffusion','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'UseProjection','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'UseConstrainB','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'rule','attrib' => {'expr' => '$UseDivbSource or $UseDivbDiffusion or $UseProjection or $UseConstrainB'},'type' => 'e','content' => [{'type' => 't','content' => '
		At least one of the options should be true.
	'}]},{'name' => 'rule','attrib' => {'expr' => 'not($UseProjection and ($UseDivbSource or $UseDivbDiffusion or $UseConstrainB))'},'type' => 'e','content' => [{'type' => 't','content' => '
		If UseProjection is true, all others should be false.
	'}]},{'name' => 'rule','attrib' => {'expr' => 'not($UseConstrainB and ($UseDivbSource or $UseDivbDiffusion or $UseProjection))'},'type' => 'e','content' => [{'type' => 't','content' => '
		If UseConstrainB is true, all others should be false.
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
'}]},{'name' => 'command','attrib' => {'name' => 'DIVBSOURCE'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => 'T','name' => 'UseB0Source','type' => 'logical'},'type' => 'e','content' => []},{'type' => 't','content' => '
#DIVBSOURCE
T			UseB0Source

! Add extra source terms related to the non-zero divergence and curl of B0.
! Default is true.
'}]},{'name' => 'command','attrib' => {'name' => 'DIVBDIFFUSION'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => '0.1666667','name' => 'DivbDiffCoeff','max' => '1','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'type' => 't','content' => '
#DIVBDIFFUSION
0.1666667		DivbDiffCoeff

! Default value is shown above. 1.0/6.0
! If divb_diffcoeff .gt. 0.5 then cfl .lt. 0.5/DivbDiffCoeff is required!
'}]},{'name' => 'command','attrib' => {'name' => 'PROJECTION'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'input' => 'select','name' => 'TypeProjectIter','type' => 'string'},'type' => 'e','content' => [{'name' => 'option','attrib' => {'default' => 'T','name' => 'Conjugate Gradients','value' => 'cg'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'BiCGSTAB','value' => 'bicgstab'},'type' => 'e','content' => []}]},{'name' => 'parameter','attrib' => {'input' => 'select','name' => 'TypeProjectStop','type' => 'string'},'type' => 'e','content' => [{'name' => 'option','attrib' => {'default' => 'T','name' => 'Relative norm','value' => 'rel'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'Maximum error','value' => 'max'},'type' => 'e','content' => []}]},{'name' => 'parameter','attrib' => {'default' => '0.1','name' => 'RelativeLimit','max' => '1','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '0.0','name' => 'AbsoluteLimit','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '50','name' => 'MaxMatvec','min' => '1','type' => 'integer'},'type' => 'e','content' => []},{'type' => 't','content' => '
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


'}]},{'name' => 'command','attrib' => {'name' => 'CORRECTP'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => '0.01','name' => 'pRatioLow','max' => '1','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '0.1','name' => 'pRatioHigh','max' => '1','min' => '$pRatioLow','type' => 'real'},'type' => 'e','content' => []},{'type' => 't','content' => '
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
'}]},{'name' => 'command','attrib' => {'name' => 'RAYTRACE'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => 'T','name' => 'UseAccurateIntegral','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'UseAccurateTrace','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '0.1','name' => 'DtExchangeRay','max' => '60','min' => '0.01','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '1','name' => 'DnRaytrace','min' => '1','type' => 'integer'},'type' => 'e','content' => []},{'type' => 't','content' => '

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
'}]},{'name' => 'command','attrib' => {'name' => 'IM'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'name' => 'TauCoupleIm','min' => '1','type' => 'real'},'type' => 'e','content' => []},{'type' => 't','content' => '

#IM
20.0			TauCoupleIm

Same as command IMCOUPLING, except it does not read the two logicals.
See description for command IMCOUPLING.

The default value is TauCoupleIm=20.0, which corresponds to typical nudging.
'}]},{'name' => 'command','attrib' => {'name' => 'IMCOUPLING'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'name' => 'TauCoupleIm','min' => '1','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => 'T','name' => 'DoCoupleImPressure','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'DoCoupleImDensity','type' => 'logical'},'type' => 'e','content' => []},{'type' => 't','content' => '

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
'}]}]},{'name' => 'commandgroup','attrib' => {'name' => 'PHYSICS PARAMETERS'},'type' => 'e','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  PHYSICS PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'name' => 'command','attrib' => {'name' => 'GAMMA','if' => '$_IsFirstSession'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => '1.6666666667','name' => 'Gamma','min' => '1','type' => 'real'},'type' => 'e','content' => []},{'type' => 't','content' => '
#GAMMA
1.6666666667		Gamma

! The adiabatic index (ratio of the specific heats for fixed pressure
! and volume. The default value is 5.0/3.0, which is valid for
! monoatomic gas or plasma.
'}]},{'name' => 'command','attrib' => {'name' => 'SHOCKTUBE'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => '1','name' => 'RhoLeft','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '0','name' => 'UnLeft','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '0','name' => 'Ut1Left','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '0','name' => 'Ut2Left','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '0.75','name' => 'BnLeft','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '1','name' => 'Bt1Left','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '0','name' => 'Bt2Left','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '1','name' => 'pRight','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '0.125','name' => 'RhoRight','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '0','name' => 'UnRight','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '0','name' => 'Ut1Right','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '0','name' => 'Ut2Right','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '0.75','name' => 'BnRight','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '-1','name' => 'Bt1Right','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '0','name' => 'Bt2Right','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '0.1','name' => 'pRight','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'input' => 'select','name' => 'ShockSlope','type' => 'real'},'type' => 'e','content' => [{'name' => 'option','attrib' => {'default' => 'T','name' => 'no rotation','value' => '0'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'value' => '0.25'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'value' => '0.3333333333333'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'value' => '0.5'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'value' => '1'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'value' => '2'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'value' => '3'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'value' => '4'},'type' => 'e','content' => []}]},{'type' => 't','content' => '
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
'}]},{'name' => 'command','attrib' => {'name' => 'SOLARWIND','if' => '$_IsFirstSession'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => '5','name' => 'SwRhoDim','min' => '-1','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '181712.175','name' => 'SwTDim','min' => '-1','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '-400','name' => 'SwUxDim','max' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '0','name' => 'SwUyDim','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '0','name' => 'SwUzDim','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '0','name' => 'SwBxDim','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '0','name' => 'SwByDim','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '5','name' => 'SwBzDim','type' => 'real'},'type' => 'e','content' => []},{'type' => 't','content' => '
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
'}]},{'name' => 'command','attrib' => {'name' => 'UPSTREAM_INPUT_FILE'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'UseUpstreamInputFile','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'if','attrib' => {'expr' => '$UseUpstreamInputFile'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'length' => '100','name' => 'NameUpstreamFile','type' => 'string'},'type' => 'e','content' => []}]},{'name' => 'rule','attrib' => {'expr' => '-f $NameUpstreamFile'},'type' => 'e','content' => [{'type' => 't','content' => '
		Upstream file $NameUpstreamFile must exist
	'}]},{'type' => 't','content' => '
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
'}]},{'name' => 'command','attrib' => {'alias' => 'MAGNETOSPHERE','name' => 'BODY','if' => '$_IsFirstSession'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'UseBody','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'if','attrib' => {'expr' => '$UseBody'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => '3','name' => 'rBody','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'if','attrib' => {'expr' => '$_NameComp eq \'GM\''},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => '4','name' => 'rCurrents','min' => '-1','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '1','name' => 'BodyRhoDim','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '10000','name' => 'BodyTDim','min' => '0','type' => 'real'},'type' => 'e','content' => []}]}]},{'type' => 't','content' => '
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
'}]},{'name' => 'command','attrib' => {'name' => 'GRAVITY','if' => '$_IsFirstSession'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'UseGravity','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'input' => 'select','name' => 'iDirGravity','if' => '$UseGravity','type' => 'integer'},'type' => 'e','content' => [{'name' => 'option','attrib' => {'default' => 'T','name' => 'central mass','value' => '0'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'X direction','value' => '1'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'Y direction','value' => '2'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'name' => 'Z direction','value' => '3'},'type' => 'e','content' => []}]},{'type' => 't','content' => '
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

'}]},{'name' => 'command','attrib' => {'name' => 'MASSLOADING'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'UseMassLoading','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'DoAccelerateMassLoading','type' => 'logical'},'type' => 'e','content' => []},{'type' => 't','content' => '
#MASSLOADING
F			UseMassLoading
F			DoAccelerateMassLoading
'}]},{'name' => 'command','attrib' => {'name' => 'HEATFLUX'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'UseHeatFlux','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => 'T','name' => 'UseSpitzerForm','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'if','attrib' => {'expr' => 'not $UseSpitzerForm'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => '1.23E-11','name' => 'Kappa0Heat','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '2.5','name' => 'Kappa0Heat','type' => 'real'},'type' => 'e','content' => []}]},{'type' => 't','content' => '
#HEATFLUX
T		UseHeatFlux
F		UseSpitzerForm
1.23E-11	Kappa0Heat [W/m/K]	! Only if not UseSpitzerForm
2.50E+00	ExponentHeat [-]	! Only if not UseSpitzerForm
'}]},{'name' => 'command','attrib' => {'name' => 'RESISTIVEFLUX'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'UseResistFlux','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => 'T','name' => 'UseSpitzerForm','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'if','attrib' => {'expr' => 'not $UseSpitzerForm'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'input' => 'select','name' => 'TypeResist','type' => 'string'},'type' => 'e','content' => [{'name' => 'option','attrib' => {'name' => 'localized','value' => 'Localized/localized'},'type' => 'e','content' => []},{'name' => 'option','attrib' => {'default' => 'T','name' => 'constant','value' => 'Constant/constant'},'type' => 'e','content' => []}]},{'name' => 'parameter','attrib' => {'default' => '9.69953E+8','name' => 'Eta0Resist','type' => 'real'},'type' => 'e','content' => []},{'name' => 'if','attrib' => {'expr' => '$TypeResist =~ /localized/i'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => '150','name' => 'Alpha0Resist','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '0.5','name' => 'yShiftResist','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '0.05','name' => 'TimeInitRise','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '1','name' => 'TimeConstLev','type' => 'real'},'type' => 'e','content' => []}]}]},{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'UseAnomResist','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'if','attrib' => {'expr' => '$UseAnomResist'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => '1.93991E+09','name' => 'Eta0AnomResist','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '1.93991E+10','name' => 'EtaAnomMaxResist','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '1','name' => 'ThresholdFactorResist','type' => 'real'},'type' => 'e','content' => []}]},{'type' => 't','content' => '
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
'}]},{'name' => 'command','attrib' => {'name' => 'TESTDISSMHD','if' => '$_IsFirstSession'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => 'T','name' => 'UseDefaultUnits','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '2.635620E-02','name' => 'Grav0Diss','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '1.640000E-01','name' => 'Beta0Diss','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '1.500000E+06','name' => 'Length0Diss','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '1.159850E+01','name' => 'Time0Diss','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '5.019000E-11','name' => 'Rho0Diss','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '1.000000E+05','name' => 'Tem0Diss','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '6.000000E-01','name' => 'Theta0Diss','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '2.500000E+01','name' => 'Delta0Diss','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '7.000000E+00','name' => 'EpsilonDiss','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '4.500000E+00','name' => 'RhoDifDiss','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '4.000000E-01','name' => 'yShiftDiss','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '5.000000E-01','name' => 'ScaleHeightDiss','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '1.159850E+01','name' => 'ScaleFactorDiss','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '5.000000E-01','name' => 'BZ0iss','type' => 'real'},'type' => 'e','content' => []},{'type' => 't','content' => '
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
'}]},{'name' => 'command','attrib' => {'name' => 'SECONDBODY','if' => '$_IsFirstSession'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'UseBody2','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'if','attrib' => {'expr' => '$UseBody2'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => '0.1','name' => 'rBody2','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '-40','name' => 'xBody2','max' => '$xMax','min' => '$xMin','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '0','name' => 'yBody2','max' => '$yMax','min' => '$yMin','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '0','name' => 'zBody2','max' => '$zMax','min' => '$zMin','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '1.3*$rBody2','name' => 'rCurrents2','min' => '$rBody2','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '5','name' => 'RhoDimBody2','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '25000','name' => 'tDimBody2','min' => '0','type' => 'real'},'type' => 'e','content' => []}]},{'type' => 't','content' => '

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
'}]},{'name' => 'command','attrib' => {'name' => 'DIPOLEBODY2','if' => '$_IsFirstSession'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'name' => 'BdpDimBody2x','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'name' => 'BdpDimBody2y','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'name' => 'BdpDimBody2z','type' => 'real'},'type' => 'e','content' => []},{'type' => 't','content' => '

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
'}]}]},{'name' => 'commandgroup','attrib' => {'name' => 'SOLAR PROBLEM TYPES'},'type' => 'e','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!! SOLAR PROBLEM TYPES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'name' => 'command','attrib' => {'name' => 'HELIOSPHERE','if' => '$_IsFirstSession and $_NameComp ne \'GM\''},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => '2.85E06','name' => 'BodyTDim','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '1.50E8','name' => 'BodyRhoDim','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '25.0','name' => 'qSun','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '1.75','name' => 'tHeat','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '1.0','name' => 'rHeat','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '4.5','name' => 'SigmaHeat','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'DoInitRope','type' => 'logical'},'type' => 'e','content' => []},{'name' => 'if','attrib' => {'expr' => '$DoInitRope'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => '0.7','name' => 'CmeA','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '1.2','name' => 'CmeR1','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '1.0','name' => 'CmeR0','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '0.23','name' => 'CmeA1','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '0.0','name' => 'CmeAlpha','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '2.5E-12','name' => 'CmeRho1','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '2.0E-13','name' => 'CmeRho2','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '0.0','name' => 'ModulationRho','max' => '10','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '0.0','name' => 'ModulationP','max' => '10','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '0.0','name' => 'OrientationGL98','max' => '360','min' => '-360','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '0.0','name' => 'LatitudeGL98','max' => '90','min' => '-90','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '0.0','name' => 'LongitudeGL98','max' => '360','min' => '-360','type' => 'real'},'type' => 'e','content' => []}]},{'type' => 't','content' => '

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
'}]},{'name' => 'command','attrib' => {'name' => 'HELIOUPDATEB0','if' => '$_NameComp ne \'GM\''},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => '0.0001','name' => 'DtUpdateB0','min' => '-1','type' => 'real'},'type' => 'e','content' => []},{'type' => 't','content' => '

#HELIOUPDATEB0
-1.0			DtUpdateB0 [s]

Set the frequency of updating the B0 field for the solar corona.
A negative value means that the B0 field is not updated.
'}]},{'name' => 'command','attrib' => {'name' => 'HELIODIPOLE','if' => '$_NameComp ne \'GM\''},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'name' => 'HelioDipoleStrength','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '0','name' => 'HelioDipoleTilt','max' => '90','min' => '-90','type' => 'real'},'type' => 'e','content' => []},{'type' => 't','content' => '

#HELIODIPOLE
-3.0                    HelioDipoleStrength [G]
 0.0                    HelioDipoleTilt     [deg]

! Variable HelioDipoleStrength defines the equatorial field strength in Gauss,
! while HelioDipoleTilt is the tilt relative to the ecliptic North 
! (negative sign means towards the planet) in degrees.
!
! Default value is HelioDipoleStrength = 0.0.
'}]},{'name' => 'command','attrib' => {'name' => 'HELIOTEST','if' => '$_NameComp ne \'GM\''},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => 'F','name' => 'DoSendMHD','type' => 'logical'},'type' => 'e','content' => []},{'type' => 't','content' => '

#HELIOTEST
F			DoSendMHD

! If DoSendMHD is true, IH sends the real MHD solution to GM in the coupling.
! If DoSendMHD is false then the values read from the IMF file are sent,
! so there is no real coupling. Mostly used for testing the framework.
!
! Default value is true, i.e. real coupling.
'}]},{'name' => 'command','attrib' => {'name' => 'HELIOBUFFERGRID','if' => '$_IsFirstSession and $_NameComp eq \'IH\''},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => '19','name' => 'rBuffMin','min' => '1','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '21','name' => 'rBuffMax','min' => '1','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '45','name' => 'nThetaBuff','min' => '18','type' => 'integer'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '90','name' => 'nPhiBuff','min' => '36','type' => 'integer'},'type' => 'e','content' => []},{'type' => 't','content' => '
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
'}]},{'name' => 'command','attrib' => {'name' => 'CME','if' => '$_IsFirstSession and $_NameComp ne \'GM\''},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'input' => 'select','name' => 'TypeCme','type' => 'string'},'type' => 'e','content' => [{'name' => 'option','attrib' => {'default' => 'T','name' => 'Low'},'type' => 'e','content' => []}]},{'name' => 'parameter','attrib' => {'default' => '0.7','name' => 'CmeA','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '1.2','name' => 'CmeR1','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '1.0','name' => 'CmeR0','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '0.23','name' => 'CmeA1','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '0.0','name' => 'CmeAlpha','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '2.5E-12','name' => 'CmeRho1','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '2.0E-13','name' => 'CmeRho2','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '1.0','name' => 'CmeB1Dim','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '4.0E5','name' => 'CmeUErupt','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'type' => 't','content' => '
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
'}]},{'name' => 'command','attrib' => {'name' => 'ARCADE','if' => '$_IsFirstSession and $_NameComp ne \'GM\''},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => '1.0E6','name' => 'tArcDim','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '1.0E-12','name' => 'RhoArcDim','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '0.718144','name' => 'bArcDim','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '1.0E6','name' => 'ByArcDim','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '5.0E3','name' => 'UzArcDim','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '0.5','name' => 'Phi0Arc','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '1.3','name' => 'MuArc','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '3','name' => 'ExpArc','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'default' => '0.5','name' => 'WidthArc','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'type' => 't','content' => '
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
'}]}]},{'name' => 'commandgroup','attrib' => {'name' => 'COMET PROBLEM TYPE'},'type' => 'e','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!! COMET PROBLEM TYPE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'name' => 'command','attrib' => {'name' => 'COMET','if' => '$_IsFirstSession'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'name' => 'ProdRate','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'name' => 'UrNeutral','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'name' => 'AverageMass','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'name' => 'IonizationRate','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'name' => 'parameter','attrib' => {'name' => 'kFriction','min' => '0','type' => 'real'},'type' => 'e','content' => []},{'type' => 't','content' => '
#COMET
1.0E28		ProdRate    - Production rate (#/s)
1.0		UrNeutral   - neutral radial outflow velocity (km/s)
17.0		AverageMass - average particle mass (amu)
1.0E-6		IonizationRate (1/s)
1.7E-9		kFriction - ion-neutral friction rate coefficient (cm^3/s)

! Only used by problem_comet.  Defaults are as shown.
'}]}]},{'name' => 'commandgroup','attrib' => {'name' => 'SCRIPT COMMANDS'},'type' => 'e','content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! SCRIPT COMMANDS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'name' => 'command','attrib' => {'name' => 'INCLUDE'},'type' => 'e','content' => [{'name' => 'parameter','attrib' => {'default' => 'Param/','length' => '100','name' => 'NameIncludeFile','type' => 'string'},'type' => 'e','content' => []},{'type' => 't','content' => '

#INCLUDE
Param/SSS_3000		NameIncludeFile

! Include a library file from Param/ or any file from anywhere else.
'}]}]},{'name' => 'rule','attrib' => {'expr' => '($SwRhoDim > 0) or $UseUpstreamInputFile or $_NameComp ne \'GM\''},'type' => 'e','content' => [{'type' => 't','content' => '
	Either command #SOLARWIND or #UPSTREAM_INPUT_FILE must be used!
'}]},{'name' => 'rule','attrib' => {'expr' => '$MaxImplBlock>1 or not $UsePartImplicit or not $MaxImplBlock'},'type' => 'e','content' => [{'type' => 't','content' => '
	Part implicit scheme requires more than 1 implicit block!
'}]},{'name' => 'rule','attrib' => {'expr' => '$MaxImplBlock==$MaxBlock or not $UseFullImplicit'},'type' => 'e','content' => [{'type' => 't','content' => '
	Full implicit scheme should be used with equal number of 
	explicit and implicit blocks!
'}]},{'name' => 'rule','attrib' => {'expr' => '-d $NameRestartOutDir or not $_IsFirstSession'},'type' => 'e','content' => [{'type' => 't','content' => '
	Output restart directory $NameRestartOutDir should exist!
'}]},{'name' => 'rule','attrib' => {'expr' => '-d $NamePlotDir or not $_IsFirstSession'},'type' => 'e','content' => [{'type' => 't','content' => '
	Plot directory $NamePlotDir should exist!
'}]}]}];