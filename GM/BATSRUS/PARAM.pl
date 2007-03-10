#^CFG FILE _FALSE_
$tree = [{'name' => 'commandList','type' => 'e','attrib' => {'name' => 'BATSRUS: GM, SC and IH Components'},'content' => [{'type' => 't','content' => '

List of MH (GM, IH and SC) commands used in the PARAM.in file


'},{'name' => 'set','type' => 'e','attrib' => {'name' => 'nI','type' => 'integer','value' => '$_GridSize[0]'},'content' => []},{'name' => 'set','type' => 'e','attrib' => {'name' => 'nJ','type' => 'integer','value' => '$_GridSize[1]'},'content' => []},{'name' => 'set','type' => 'e','attrib' => {'name' => 'nK','type' => 'integer','value' => '$_GridSize[2]'},'content' => []},{'name' => 'set','type' => 'e','attrib' => {'name' => 'MaxBlock','type' => 'integer','value' => '$_GridSize[3]'},'content' => []},{'name' => 'set','type' => 'e','attrib' => {'name' => 'MaxImplBlock','type' => 'integer','value' => '$_GridSize[4]'},'content' => []},{'name' => 'set','type' => 'e','attrib' => {'name' => 'MaxBlockALL','type' => 'integer','value' => '$_nProc and $MaxBlock and $_nProc*$MaxBlock'},'content' => []},{'name' => 'set','type' => 'e','attrib' => {'name' => 'NameRestartOutDir','type' => 'string','value' => '$_NameComp/restartOUT'},'content' => []},{'name' => 'set','type' => 'e','attrib' => {'name' => 'NamePlotDir','type' => 'string','value' => '$_NameComp/IO2'},'content' => []},{'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'STAND ALONE MODE'},'content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! STAND ALONE PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

'},{'name' => 'command','type' => 'e','attrib' => {'name' => 'COMPONENT'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'NameComp','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'SC','default' => 'T','if' => '$_NameComp eq \'SC\''},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'IH','default' => 'T','if' => '$_NameComp eq \'IH\''},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'GM','default' => 'T','if' => '$_NameComp eq \'GM\''},'content' => []}]},{'type' => 't','content' => '

#COMPONENT
GM			NameComp

This command is only used in the stand alone mode.

The NameComp variable contains the two-character component ID
for the component which BATSRUS is representing.
If NameComp does not agree with the value of the NameThisComp
variable, BATSRUS stops with an error message.
This command is saved into the restart header file for consistency check.

There is no default value: if the command is not given, the component ID is not checked.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'DESCRIPTION','if' => '$_IsStandAlone'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'StringDescription','length' => '100','type' => 'string'},'content' => []},{'type' => 't','content' => '

#DESCRIPTION
This is a test run for Jupiter with no rotation.

This command is only used in the stand alone mode.

The StringDescription string can be used to describe the simulation
for which the parameter file is written. The #DESCRIPTION command and
the StringDescription string are saved into the restart file,
which helps in identifying the restart files.

The default value is ``Please describe me!", which is self explanatory.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'ECHO','if' => '$_IsStandAlone'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoEcho','type' => 'logical','default' => 'F'},'content' => []},{'type' => 't','content' => '

#ECHO
T                       DoEcho

This command is only used in the stand alone mode.

If the DoEcho variable is true, the input parameters are echoed back.
The default value for DoEcho is .false., but it is a good idea to
set it to true at the beginning of the PARAM.in file.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'PROGRESS','if' => '$_IsStandAlone'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DnProgressShort','type' => 'integer','default' => '10','min' => '-1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DnProgressLong','type' => 'integer','default' => '100','min' => '-1'},'content' => []},{'type' => 't','content' => '
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
to accelerate the convergence towards steady state.

The default value depends on how the stand alone code was installed.
See the description of the NEWPARAM command.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'BEGIN_COMP','multiple' => 'T','if' => '$_IsStandAlone'},'content' => [{'type' => 't','content' => '

This command is allowed in stand alone mode only for the sake of the 
test suite, which contains these commands when the framework is tested.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'END_COMP','multiple' => 'T','if' => '$_IsStandAlone'},'content' => [{'type' => 't','content' => '

This command is allowed in stand alone mode only for the sake of the 
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
'}]}]},{'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'PLANET PARAETERS'},'content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! PLANET PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

The planet commands can only be used in stand alone mode.
The commands allow to work with an arbitrary planet.
It is also possible to change some parameters of the planet relative
to the real values.

By default Earth is assumed with its real parameters.
Another planet (moon, comet) can be selected with the #PLANET
(#MOON, #COMET) command. The real planet parameters can be modified 
and simplified with the other planet commands listed in this subsection.
These modified commands cannot precede the #PLANET command!

'},{'name' => 'command','type' => 'e','attrib' => {'name' => 'PLANET','if' => '$_IsFirstSession and $_IsStandAlone','alias' => 'MOON,COMET'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'NamePlanet','type' => 'string','case' => 'upper','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'MERCURY'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'VENUS'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'EARTH','default' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'MARS'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'JUPITER'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'SATURN'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'URANUS'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'NEPTUNE'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'PLUTO'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'MOON'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'IO'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'EUROPA'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'TITAN'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'NEW'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'HALLEY'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'COMET1P'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'BORRELLY'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'COMET19P'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'CHURYUMOVGERASIMENKO'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'COMET67P'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'HALEBOPP'},'content' => []}]},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$NamePlanet eq \'NEW\''},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'RadiusPlanet','type' => 'real','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'MassPlanet','type' => 'real','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'OmegaPlanet','type' => 'real','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TiltRotation','type' => 'real','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypeBField','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'NONE'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'DIPOLE','default' => 'T'},'content' => []}]}]},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$TyepBField eq \'DIPOLE\''},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'MagAxisThetaGeo','type' => 'real','min' => '0','max' => '180'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'MagAxisPhiGeo','type' => 'real','min' => '0','max' => '360'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DipoleStrength','type' => 'real'},'content' => []}]},{'name' => 'rule','type' => 'e','attrib' => {'expr' => 'not $PlanetCommand'},'content' => [{'type' => 't','content' => '
		PLANET should precede $PlanetCommand
	'}]},{'type' => 't','content' => '

#PLANET
NEW			NamePlanet (rest of parameters read for unknown planet)
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

The default value is NamePlanet="Earth". Although many other planets
and some of the moons are recognized, some of the parameters, 
like the equinox time are not yet properly set.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'ROTATIONAXIS','if' => '$_IsStandAlone'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'IsRotAxisPrimary','type' => 'logical','default' => 'T'},'content' => []},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$IsRotAxisPrimary'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'RotAxisTheta','type' => 'real','min' => '0','max' => '180'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'RotAxisPhi','type' => 'real','min' => '0','max' => '360'},'content' => []}]},{'name' => 'set','type' => 'e','attrib' => {'name' => 'PlanetCommand','type' => 'string','value' => 'ROTATIONAXIS'},'content' => []},{'type' => 't','content' => '

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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'MAGNETICAXIS','if' => '$_IsStandAlone'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'IsMagAxisPrimary','type' => 'logical','default' => 'T'},'content' => []},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$IsMagAxisPrimary'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'MagAxisTheta','type' => 'real','min' => '0','max' => '180'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'MagAxisPhi','type' => 'real','min' => '0','max' => '360'},'content' => []}]},{'name' => 'set','type' => 'e','attrib' => {'name' => 'PlanetCommand','type' => 'string','value' => 'MAGNETICAXIS'},'content' => []},{'type' => 't','content' => '

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
-3.11e-5		DipoleStrength [Tesla]

The DipoleStrength variable contains the
magnetic equatorial strength of the dipole magnetic field in Tesla.

The default value is the real dipole strength for the planet.
For the Earth the default is taken to be -31100 nT.
The sign is taken to be negative so that the magnetic axis can
point northward as usual.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'UPDATEB0','if' => '$_IsStandAlone'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DtUpdateB0','type' => 'real','default' => '0.0001','min' => '-1'},'content' => []},{'type' => 't','content' => '

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
it is identical with the commands:

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

'},{'name' => 'command','type' => 'e','attrib' => {'name' => 'USERFLAGS','alias' => 'USER_FLAGS'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseUserInnerBcs','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseUserSource','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseUserPerturbation','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseUserOuterBcs','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseUserICs','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseUserSpecifyRefinement','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseUserLogFiles','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseUserWritePlot','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseUserAMR','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseUserEchoInput','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseUserB0','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseUserInitSession','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseUserUpdateStates','type' => 'logical','default' => 'F'},'content' => []},{'type' => 't','content' => '

#USERFLAGS
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
F                       UseUserInitSession
F                       UseUserUpdateStates

This command controls the use of user defined routines in ModUser.f90.
For each flag that is set, an associated routine will be called in 
the user module.  Default is .false. for all flags.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'USERINPUTBEGIN'},'content' => [{'type' => 't','content' => '

#USERINPUTBEGIN

This command signals the beginning of the section of the file which 
is read by the subroutine user\\_read\\_inputs in the ModUser.f90 file.
The section ends with the #USERINPUTEND command. 
There is no XML based parameter checking in the user section.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'USERINPUTEND'},'content' => [{'type' => 't','content' => '

#USERINPUTEND

This command signals the end of the section of the file which 
is read by the subroutine user\\_read\\_inputs in the ModUser.f90 file.
The section begins with the #USERINPUTBEGIN command. 
There is no XML based parameter checking in the user section.
'}]}]},{'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'TESTING AND TIMING'},'content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  TESTING AND TIMING PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'name' => 'command','type' => 'e','attrib' => {'name' => 'TEST'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TestString','length' => '100','type' => 'string'},'content' => []},{'type' => 't','content' => '
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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'TESTIJK'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'iTest','type' => 'integer','min' => '-2','max' => '$nI+2'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'jTest','type' => 'integer','min' => '-2','max' => '$nJ+2'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'kTest','type' => 'integer','min' => '-2','max' => '$nK+2'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'iBlockTest','type' => 'integer','min' => '1','max' => '$MaxBlock'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'iProcTest','type' => 'integer','min' => '0'},'content' => []},{'type' => 't','content' => '
#TESTIJK
1                       iTest           (cell index for testing)
1                       jTest           (cell index for testing)
1                       kTest           (cell index for testing)
1                       iBlockTest      (block index for testing)
0                       iProcTest       (processor index for testing)

! The location of test info in terms of indices, block and processor number.
! Note that the user should set #TESTIJK or #TESTXYZ, not both.  If both
! are set, the final one in the session will set the test point.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'TESTXYZ'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'xTest','type' => 'real','min' => '$xMin','max' => '$xMax'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'yTest','type' => 'real','min' => '$yMin','max' => '$yMax'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'zTest','type' => 'real','min' => '$zMin','max' => '$zMax'},'content' => []},{'type' => 't','content' => '
#TESTXYZ
1.5                     xTest           (X coordinate of cell for testing)
-10.5                   yTest           (Y coordinate of cell for testing)
-10.                    zTest           (Z coordinate of cell for testing)

! The location of test info in terms of coordinates.
! Note that the user should set #TESTIJK or #TESTXYZ, not both.  If both
! are set, the final one in the session will set the test point.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'TESTTIME'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nIterTest','type' => 'integer','default' => '-1','min' => '-1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TimeTest','type' => 'real','default' => '1e30','min' => '-1'},'content' => []},{'type' => 't','content' => '

#TESTTIME
-1                      nIterTest       (iteration number to start testing)
10.5                    TimeTest        (time to start testing in seconds)

! The time step and physical time to start testing.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'TESTVAR'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'iVarTest','type' => 'integer','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'Rho','default' => 'T','value' => '1'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'RhoUx','value' => '2'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'RhoUy','value' => '3'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'RhoUz','value' => '4'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Bx','value' => '5'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'By','value' => '6'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Bz','value' => '7'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'e','value' => '8'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'p','value' => '9'},'content' => []}]},{'type' => 't','content' => '
#TESTVAR
1                       iVarTest

! Index of variable to be tested. Default is rho_="1", i.e. density.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'TESTDIM'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'iVarTest','type' => 'integer','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'all','value' => '0'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'x','default' => 'T','value' => '1'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'y','value' => '2'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'z','value' => '3'},'content' => []}]},{'type' => 't','content' => '
#TESTDIM
1                       iDimTest

! Index of dimension/direction to be tested. Default is X dimension.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'STRICT','multiple' => 'T'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseStrict','type' => 'logical','default' => 'T'},'content' => []},{'type' => 't','content' => '
#STRICT
T                       UseStrict

! If true then stop when parameters are incompatible. If false, try to
! correct parameters and continue. Default is true, i.e. strict mode
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'VERBOSE'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'lVerbose','type' => 'integer','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'errors and warnings only','value' => '-1'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'start and end of sessions','value' => '0'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'normal','default' => 'T','value' => '1'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'calls on test processor','value' => '10'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'calls on all processors','value' => '100'},'content' => []}]},{'type' => 't','content' => '
#VERBOSE
-1                      lVerbose

! Verbosity level controls the amount of output to STDOUT. Default level is 1.
!\\\\
!   lVerbose $\\leq -1$ only warnings and error messages are shown.\\\\
!   lVerbose $\\geq  0$ start and end of sessions is shown.\\\\
!   lVerbose $\\leq  1$ a lot of extra information is given.\\\\
!   lVerbose $\\leq 10$ all calls of set_oktest are shown for the test processor.\\\\
!   lVerbose $\\leq 100$ all calls of set_oktest are shown for all processors.\\\\
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'DEBUG'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoDebug','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoDebugGhost','type' => 'logical','default' => 'F'},'content' => []},{'type' => 't','content' => '
#DEBUG
F                       DoDebug         (use it as if(okdebug.and.oktest)...)
F                       DoDebugGhost    (parameter for show_BLK in library.f90)

! Excessive debug output can be controlled by the global okdebug parameter
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'CODEVERSION','if' => '$_IsFirstSession'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CodeVersion','type' => 'real','default' => '7.50','min' => '0'},'content' => []},{'type' => 't','content' => '
#CODEVERSION
7.50                    CodeVersion

! Checks CodeVersion. Prints a WARNING if it differs from the CodeVersion
! defined in ModMain.f90. Used in newer restart header files. 
! Should be given in PARAM.in when reading old restart files, 
! which do not have version info in the header file.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'USERMODULE','if' => '$_IsFirstSession'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'NameUserModule','length' => '100','type' => 'string','default' => 'EMPTY'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'VersionUserModule','type' => 'real','default' => '1.0','min' => '0'},'content' => []},{'type' => 't','content' => '

#USERMODULE
TEST PROBLEM Smith
1.3			VersionUserModule

Checks the selected user module. If the name or the version number
differs from that of the compiled user module, a warning is written,
and the code stops in strict mode (see #STRICT command). 
This command and its parameters are written into the restart header file too,
so the user module is checked when a restart is done. 
There are no default values. If the command is not present, the user 
module is not checked.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'EQUATION','if' => '$_IsFirstSession'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'NameEquation','length' => '100','type' => 'string','default' => 'MHD'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nVar','type' => 'integer','default' => '8'},'content' => []},{'type' => 't','content' => '
#EQUATION
MHD			NameEquation
8			nVar

! Define the equation name and the number of variables.
! If any of these do not agree with the values determined 
! by the code, BATSRUS stops with an error. Used in restart
! header files and can be given in PARAM.in as a check
! and as a description.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'PRECISION','if' => '$_IsFirstSession'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nByteReal','type' => 'integer','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'single precision (4)','default' => '$_nByteReal==4','value' => '4'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'double precision (8)','default' => '$_nByteReal==8','value' => '8'},'content' => []}]},{'name' => 'rule','type' => 'e','attrib' => {'expr' => '$nByteReal==$_nByteReal or not $UseStrict'},'content' => [{'type' => 't','content' => '
		nByteReal in file must agree with _nByteReal in strict mode.
	'}]},{'type' => 't','content' => '

#PRECISION
8                       nByteReal

Define the number of bytes in a real number. If it does not agree
with the value determined by the code, BATSRUS stops with an error
unless the strict mode is switched off.
This is used in restart header files to store (and check) the precision
of the restart files. It is now possible to read restart files with
a precision that differs from the precision the code is compiled with,
but strict mode has to be switched off with the #STRICT command.
The #PRECISION command may also be used to enforce a certain precision.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'CHECKGRIDSIZE','if' => '$_IsFirstSession'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nI','type' => 'integer','default' => '$nI','min' => '$nI','max' => '$nI'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nJ','type' => 'integer','default' => '$nJ','min' => '$nJ','max' => '$nJ'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nK','type' => 'integer','default' => '$nK','min' => '$nK','max' => '$nK'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'MinBlockALL','type' => 'integer','min' => '1','max' => '$MaxBlockALL'},'content' => []},{'type' => 't','content' => '

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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'TIMING','if' => '$_IsStandAlone'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseTiming','type' => 'logical','default' => 'T'},'content' => []},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$UseTiming'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Frequency','type' => 'integer','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'none','value' => '-3'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'final only','default' => 'T','value' => '-2'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'end of sessions','value' => '-1'},'content' => []},{'name' => 'optioninput','type' => 'e','attrib' => {'name' => 'every X steps','default' => '100','min' => '1'},'content' => []}]},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nDepthTiming','type' => 'integer','default' => '-1','min' => '-1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypeTimingReport','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'cumulative','default' => '1','value' => 'cumu'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'list'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'tree'},'content' => []}]}]},{'type' => 't','content' => '
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
'}]}]},{'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'INITIAL AND BOUNDARY CONDITIONS'},'content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! MAIN INITIAL AND BOUNDARY CONDITION PARAMETERS  !!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'name' => 'command','type' => 'e','attrib' => {'name' => 'NORMALIZATION','if' => '$_IsFirstSession'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypeNormalization','type' => 'string','case' => 'upper','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'SI'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'PLANETARY','default' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'SOLARWIND'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'NONE'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'READ'},'content' => []}]},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$TypeNormalization =~ /READ/'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'No2SiUnitX','type' => 'real','default' => '1','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'No2SiUnitU','type' => 'real','default' => '1','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'No2SiUnitRho','type' => 'real','default' => '1','min' => '0'},'content' => []}]},{'type' => 't','content' => '
#NORMALIZATION
READ			TypeNormalization
1000.0			No2SiUnitX   (only read if TypeNormalization=READ)
1000.0			No2SiUnitU   (only read if TypeNormalization=READ)
1.0e-6			No2SiUnitRho (only read if TypeNormalization=READ)

This command determines what units are used internally in BATSRUS.
The units are normalized so that several physical constants become
unity (e.g. the permeability of vacuum), so the equations are simpler
in the code. The normalization also helps to keep the various
quantities within reasonable ranges. For example density of space
plasma is very small in SI units, so it is better to use some normalization,
like amu/cm$^3$. Also note that distances and positions (like grid size,
grid resolution, plotting resolution, radius of the inner body etc) 
are always read in normalized units from the PARAM.in file.
Other quantities are read in I/O units (see the #IOUNITS command).

The normalization of the distance, velocity and density are 
determined by the TypeNormalization parameter. The normalization
of all other quantities are derived from these three values.
It is important to note that the normalization of number density
(defined as the density normalization divided by the proton mass)
is usually not consistent with the inverse cube of the normalization 
of distance. 

Possible values for TypeNormalization are NONE, PLANETARY, SOLARWIND, 
and READ. 

If TypeNormalization="NONE" then the distance, velocity and density
units are the SI units, i.e. meter, meter/sec, and kg/m$^3$.
Note that the magnetic field and the temperature are still normalized
differently from SI units so that the Alfven speed is $B/\\sqrt{\\rho}$
and the ion temperature is simply $p/(\\rho/AverageIonMass)$,
where the AverageIonMass is given relative to the mass of proton.

If TypeNormalization="PLANETARY" then the distance unit is the
radius of the central body (moon, planet, or the Sun). If there is
no central body, the length normalization is 1km. The velocity unit
is km/s, and the density unit is amu/cm3.

If TypeNormalization="SOLARWIND" then the distance unit is the 
same as for PLANETARY units, but the velocity and density are normalized to
the density and the sound speed of the solar wind. This normalization
is very impractical, because it depends on the solar wind values
that are variable, and may not even make sense (e.g. for a shock tube test
or a Tokamak problem). This normalization is only kept for sake of 
backwards compatibility.

Finally TypeNormalization="READ" reads the three basic normalization units from
the PARAM.in file as shown in the example. This allows arbitrary normalization.

The restart header file saves the normalization with TypeNormalization="READ"
and the actual values of the distance, velocity and density normalization
factors. This avoids the problem of continuing the run with inconsistent
normalization (e.g. if the SOLARWIND normalization is used and the
solar wind parameters have been changed). It also allows other programs
to read the data saved in the restart files and convert them to appropriate
units. 

The default normalization is PLANETARY.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'IOUNITS','if' => '$_IsFirstSession'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypeIoUnit','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'SI'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'PLANETARY','default' => '$_NameComp eq \'GM\''},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'HELIOSPHERIC','default' => '$_NameComp ne \'GM\''},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'NONE'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'USER'},'content' => []}]},{'type' => 't','content' => '
#IOUNITS
PLANETARY			TypeIoUnit

This command determines the physical units of various parameters read from the
PARAM.in file and written out into log files and plot files (if they are
dimensional. The units are determined by the TypeIoUnit string.
Note that distances and positions are always read in normalized units from 
PARAM.in but they are written out in I/O units. 
In most cases the two coincides.

Also note that the I/O units are NOT necessarily physically consistent units. 
For example one cannot divide distance with time and compare it with
the velocity because they may be in inconsistent units. 
One needs to convert into some consistent units before the various quantities 
can be combined.

If TypeIoUnits="SI" the input and output values are taken in SI units
(m, s, kg, etc).

The PLANETARY units use the radius of the planet for distance, 
seconds for time, amu/cm$^3$ for mass density, cm$^-3$ for number denisty, 
km/s for speed, nPa for pressure, nT for magnetic field, micro Amper/m$^2$ 
for current density, mV/m for electric field, nT/planet radius for div B,
and degrees for angles. For any other quantity SI units are used.
If there is no planet (see the #PLANET command) then the distance
unit is 1 km.

The HELIOSPHERIC units use the solar radius for distance, seconds for time,
km/s for velocity, degrees for angle, and
CGS units for mass density, number density, pressure,
magnetic field, momentum, energy density, current, and div B.

When TypeIoUnit="NONE" the input and output units are the same as
the normalized units (see the #NORMALIZATION command).

Finally when TypeIoUnit="USER", the user can modify the I/O units 
(Io2Si_V) and the names of the units (NameTecUnit_V and NameIdlUnit_V)
in the subroutine user_io_units of the user module.
Initially the values are set to SI units.

The #IOUNITS command and the value of TypeIoUnits is saved into the 
restart header file so that one continues with the same I/O units 
after restart.

The default is "PLANETARY" unit if BATSRUS is used as the GM component
and "HELIOSPHERIC" otherwise (SC or IH). 
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'RESTARTINDIR','if' => '$_IsFirstSession'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'NameRestartInDir','length' => '100','type' => 'string','default' => 'GM/restartIN'},'content' => []},{'name' => 'rule','type' => 'e','attrib' => {'expr' => '-d $NameRestartInDir'},'content' => [{'type' => 't','content' => '
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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'RESTARTINFILE','if' => '$_IsFirstSession'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypeRestartInFile','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'block','default' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'one'},'content' => []}]},{'type' => 't','content' => '

#RESTARTINFILE
one			TypeRestartInFile

This command is normally saved in the restart header file, and it describes
how the restart data was saved: 
into separate files for each grid block (\'block\') 
or into a single direct access file (\'one\').
The default value is \'block\'.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'NEWRESTART','if' => '$_IsFirstSession'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoRestartBFace','type' => 'logical','default' => 'F'},'content' => []},{'type' => 't','content' => ' 

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

'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'OUTERBOUNDARY'},'content' => [{'name' => 'foreach','type' => 'e','attrib' => {'values' => 'TypeBcEast,TypeBcWest,TypeBcSouth,TypeBcNorth,TypeBcBot,TypeBcTop','name' => 'Side'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => '$Side','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'coupled'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'vary/inflow','default' => '$Side eq \'TypeBcWest\''},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'float/outflow','default' => '$Side ne \'TypeBcWest\''},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'heliofloat'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'reflect'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'periodic'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'fixed'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'shear'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'linetied'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'raeder'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'arcadetop'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'arcadebot'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'arcadebotcont'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'user'},'content' => []}]}]},{'name' => 'rule','type' => 'e','attrib' => {'expr' => 'not($TypeBcEast eq \'periodic\' xor $TypeBcWest eq \'periodic\')'},'content' => [{'type' => 't','content' => '
		East and west BCs must be both periodic or neither
	'}]},{'name' => 'rule','type' => 'e','attrib' => {'expr' => 'not($TypeBcSouth eq \'periodic\' xor $TypeBcNorth eq \'periodic\')'},'content' => [{'type' => 't','content' => '
		South and North BCs must be both periodic or neither
	'}]},{'name' => 'rule','type' => 'e','attrib' => {'expr' => 'not($TypeBcBot eq \'periodic\' xor $TypeBcTop eq \'periodic\')'},'content' => [{'type' => 't','content' => '
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
! fixed         - fixed solarwind values
! fixedB1       - fixed solarwind values without correction for the dipole B0
! float/outflow - zero gradient
! heliofloat    - floating for the SC component (requires #FACEOUTERBC)
! linetied      - float P, rho, and B, reflect all components of U
! raeder        - Jimmy Raeder\'s BC
! reflect       - reflective
! periodic      - periodic
! vary/inflow   - time dependent BC (same as fixed for non time_accurate)
! shear         - sheared (intended for shock tube problem only)
! arcadetop     - intended for arcade problem only
! arcadebot     - intended for arcade problem only
! arcadebotcont - intended for arcade problem only
! user          - user defined
!\\end{verbatim}
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'INNERBOUNDARY'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypeBcInner','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'reflect'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'float'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'fixed'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'ionosphere','default' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'ionosphereB0','value' => 'ionosphereB0/ionosphereb0'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'ionospherefloat'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'coronatoih'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'user'},'content' => []}]},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$UseBody2'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypeBcBody2','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'reflect','default' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'float'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'fixed'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'ionosphere'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'ionosphereB0'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'ionospherefloat'},'content' => []}]}]},{'name' => 'rule','type' => 'e','attrib' => {'expr' => 'not($TypeBcBody2 =~ /ionosphere/)'},'content' => [{'type' => 't','content' => ' 
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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'EXTRABOUNDARY'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseExtraBoundary','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$UseExtraBoundary'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypeExtraBoundary','type' => 'string'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoFixExtraBoundary','type' => 'logical','default' => 'F'},'content' => []}]},{'type' => 't','content' => '

#EXTRABOUNDARY
T		UseExtraBoundary
user		TypeExtraBoundary
F		DoFixExtraboundary

It UseExtraBoundary is true, the user can define an extra boundary
condition in the user files. The TypeExtraBoundary can be used
to select a certain type. The DoFixExtraboundary controls
if resolution change is allowed.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'FACEOUTERBC'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'MaxBoundary','type' => 'integer','default' => '0','min' => '0','max' => '6'},'content' => []},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$MaxBoundary >= 1'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoFixOuterBoundary','type' => 'logical','default' => 'F'},'content' => []}]},{'type' => 't','content' => '
#FACEOUTERBC
0              MaxBoundary            
F              DoFixOuterBoundary)    !read only for MaxBoundary>=East_(=1).

If MaxBoundary is East_(=1) or more then the outer boundaries indexed
between East_ and MaxBoundary are treated using set_BCs.f90 subroutines 
(face values are defined) instead of set_outerBCs.f90 (ghost cell values
are defined).
If DoFixOuterBoundary is .true., there is no resolution
change along the outer boundaries indexed from East_ to MaxBoundary.
'}]}]},{'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'GRID GEOMETRY'},'content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!! GRID GEOMETRY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

'},{'name' => 'command','type' => 'e','attrib' => {'name' => 'GRID','required' => 'T','if' => '$_IsFirstSession'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nRootBlockX','type' => 'integer','default' => '2','min' => '1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nRootBlockY','type' => 'integer','default' => '1','min' => '1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nRootBlockZ','type' => 'integer','default' => '1','min' => '1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'xMin','type' => 'real','default' => '-192.0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'xMax','type' => 'real','default' => '  64.0','min' => '$xMin'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'yMin','type' => 'real','default' => ' -64.0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'yMax','type' => 'real','default' => '  64.0','min' => '$yMin'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'zMin','type' => 'real','default' => ' -64.0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'zMax','type' => 'real','default' => '  64.0','min' => '$zMin'},'content' => []},{'type' => 't','content' => '
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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'COORDSYSTEM','multiple' => 'T','alias' => 'COORDINATESYSTEM'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypeCoordSystem','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'GSM','default' => 'T','if' => '$_NameComp eq \'GM\''},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'HGI','default' => 'T','if' => '$_NameComp eq \'IH\''},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'HGC','if' => '$_NameComp eq \'IH\''},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'HGR','default' => 'T','if' => '$_NameComp eq \'SC\''},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'HGI','if' => '$_NameComp eq \'SC\''},'content' => []}]},{'type' => 't','content' => '

#COORDSYSTEM
GSM			TypeCoordSystem

! TypeCoordSystem defines the coordinate system for the component.
! Currently only one coordinate system is available for GM ("GSM")
! and two for IH ("HGI" or "HGC") and two for SC ("HGR" or "HGI"). 
! In the future "GSE" should be also an option for GM.
! The coordinate systems are defined in share/Library/src/CON_axes.
!
! Default is component dependent: "GSM" for GM, "HGI" for IH, and "HGR" for SC.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'GRIDGEOMETRY','if' => '$_IsFirstSession','alias' => 'COVARIANTGEOMETRY'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypeGeometry','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'cartesian','default' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'cylindrical'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'spherical'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'spherical_lnr'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'axial_torus'},'content' => []}]},{'type' => 't','content' => '
#GRIDGEOMETRY
spherical_lnr			TypeGeometry
		
This command determines the geometry of the grid. Possible values
are cartesian, cylindrical, spherical, spherical_lnr (spherical
with a logarithmic radial coordinate, and axial (axially symmetric
toroidal grid).

BATSRUS uses a Cartesian geometry by default. Setting this command
allows using more complicated grids. Note that setting "cartesian"
is not identical with not setting the command at all, because
the code will use the general formulas. The results should be essentially
identical though.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'LIMITRADIUS','if' => '$_IsFirstSession','alias' => 'LIMITGENCOORD1'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'rMin','type' => 'real','default' => '0','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'rMax','type' => 'real','default' => '1000','min' => '$rMin'},'content' => []},{'type' => 't','content' => '
#LIMITRADIUS
 10.0			rMin
100.0			rMax

This command allows limiting the first coordinate (usually the radius) 
of the grid when general coordinates are used. One can exclude the 
origin of a spherical grid, for example.

The default is to use the whole radial extent allowed by the #GRID command.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'TORUSSIZE'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'rTorusLarge','type' => 'real','default' => '6.0','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'rTorusSmall','type' => 'real','default' => '0.5','min' => '0'},'content' => []},{'type' => 't','content' => '
#TORUSSIZE
6.0		rTorusLarge
0.5		rTorusSmall

This command can sets the large and small radii for a toroidal grid
(see the #GRID command). The default values are shown.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'VERTEXBASEDGRID','if' => '$_IsFirstSession and $_IsStandAlone'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseVertexBasedGrid','type' => 'logical','default' => 'F'},'content' => []},{'type' => 't','content' => '

#VERTEXBASEDGRID
F                         UseVertexBasedGrid

For a vertex-based logically cartesian (spherical, cylindircal) grid 
(UseVertexBasedGrid=.true.) the node coordinates are defined
in terms of an arbitrary pointwide transformation of nodes of an 
original cartesian (spherical,cylindrical) block adaptive grid.
Advantage: the possiblity to use the arbitrary transformation.
Disadvantages: the cell center coordinates can not be definied unambigously
and the difference of the state variables across the face does not evaluate
the gradient in the direction, normal to this face (stricly speaking).
Cell-centered grids are used if UseVertexBasedGrid=.false. (default value)
Advantage: for some particular geometries (spherical, cylindrical) the 
control volumes are the Voronoy cells (any face is perpendicular to the line
connecting the centers of the neighboring cells). 
Disadvantages: even in these particular cases it is not easy to properly 
define the face area vectors at the resolution change. More general 
cell-centered grid either is not logically cartesian, or does not consist of 
the Voronoy cells only.

'}]}]},{'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'INITIAL TIME'},'content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! INITIAL TIME AND STEP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'name' => 'command','type' => 'e','attrib' => {'name' => 'STARTTIME','if' => '$_IsFirstSession','alias' => 'SETREALTIME'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'iYear','type' => 'integer','default' => '2000'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'iMonth','type' => 'integer','default' => '3','min' => '1','max' => '12'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'iDay','type' => 'integer','default' => '21','min' => '1','max' => '31'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'iHour','type' => 'integer','default' => '0','min' => '0','max' => '23'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'iMinute','type' => 'integer','default' => '0','min' => '0','max' => '59'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'iSecond','type' => 'integer','default' => '0','min' => '0','max' => '59'},'content' => []},{'type' => 't','content' => '
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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'TIMESIMULATION','if' => '$_IsFirstSession'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'tSimulation','type' => 'real','default' => '0.0','min' => '0'},'content' => []},{'type' => 't','content' => '

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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'NSTEP','if' => '$_IsFirstSession'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nStep','type' => 'integer','default' => '0','min' => '0'},'content' => []},{'type' => 't','content' => '

#NSTEP
100			nStep

! Set nStep for the component. Typically used in the restart.H header file.
! Generally it is not inserted in a PARAM.in file by the user.
!
! The default is nStep=0 as the starting time step with no restart.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'NPREVIOUS','if' => '$_IsFirstSession'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nPrevious','type' => 'integer','default' => '-1','min' => '-1'},'content' => []},{'type' => 't','content' => '

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
'},{'name' => 'command','type' => 'e','attrib' => {'name' => 'TIMESTEPPING'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nStage','type' => 'integer','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'value' => '1'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'default' => 'T','value' => '2'},'content' => []}]},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CflExpl','type' => 'real','default' => '0.8','min' => '0','max' => '1'},'content' => []},{'type' => 't','content' => '

#TIMESTEPPING
2                       nStage (1 or 2)
0.80                    CflExpl

Parameters for time integration. For explicit time integration nStage
is the number of stages in the Runge-Kutta scheme. For implicit time
stepping nStage=2 corresponds to the BDF2 scheme that
uses the previous time step to make the scheme 2nd order accurate in time.
For explicit time stepping the CPU time is proportional to the number of 
stages. In time accurate runs the 1-stage explicit time stepping scheme 
may work reasonably well with second order spatial discretization, 
especially if the time step is limited to a very small value. 
Using a one stage scheme can speed up the code by a factor of two with 
little compromise in accuracy.

For local time stepping (steady state mode) one should always use the 
2-stage scheme with 2-nd order spatial accuracy to avoid oscillations
(or make the CFL number small, about 0.4). 

For implicit scheme the second order BDF2 scheme is more accurate 
but not more expensive than the first order backward Euler scheme, 
so it is a good idea to use nStage=nOrder. 

To achieve consistency between the spatial and temporal orders of accuracy,
the #SCHEME command always sets nStage to be the same as nOrder. 
The #TIMESTEPPING command can be used to reset nStage if required.

Default is 2-stage scheme with CflExpl=0.8
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'FIXEDTIMESTEP'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseDtFixed','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DtFixedDim','type' => 'real','default' => '1.0','min' => '0','if' => '$UseDtFixed'},'content' => []},{'type' => 't','content' => '
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

'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'PARTSTEADY'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UsePartSteady','type' => 'logical','default' => 'F'},'content' => []},{'type' => 't','content' => '

! If UsePartSteady is true, the partially steady state algorithm is used.
! Only blocks which are changing or next to changing blocks are evolved.
! This scheme can speed up the calculation if part 
! of the domain is in a numerical steady state. 
! In steady state runs the code stops when a full steady state is
! achieved. The conditions for checking the numerical steady state are set 
! by the #PARTSTEADYCRITERIA command.
! Default value is UsePartSteady = .false.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'PARTSTEADYCRITERIA','alias' => 'STEADYCRITERIA'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'MinCheckVar','type' => 'integer','default' => '1','min' => '1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'MaxCheckVar','type' => 'integer','default' => '8','min' => '$MinCheckVar'},'content' => []},{'name' => 'for','type' => 'e','attrib' => {'from' => '$MinCheckVar','to' => '$MaxCheckVar'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'RelativeEps','type' => 'real','default' => '0.001','min' => '0','max' => '1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'AbsoluteEps','type' => 'real','default' => '0.0001','min' => '0'},'content' => []}]},{'type' => 't','content' => '
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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'POINTIMPLICIT'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UsePointImplicit','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'BetaPointImplicit','type' => 'real','default' => '0.5','min' => '0.5','if' => '$UsePointImplicit','max' => '1'},'content' => []},{'type' => 't','content' => '
#POINTIMPLICIT
T		UsePointImplicit
0.5		BetaPointImplicit (read if UsePointImplicit is true)

Switches on or off the point implicit scheme. The BetaPointImplicit
parameter (in the 0.5 to 1.0 range) determines the order of accuracy 
for a 2-stage scheme. If BetaPointImplicit=0.5 the point implicit scheme 
is second order accurate in time when used in a 2-stage scheme. 
Larger values may be more robust, but only first order accurate in time.
For a 1-stage scheme (and in the 1st stage of a 2-stage scheme) the
BetaPointImplicit parameter is ignored, and the coefficient is set to 1.

The default values are UsePointImplicit=false and BetaPointImplicit=0.5.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'PARTLOCAL'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UsePartLocal','type' => 'logical','default' => 'F'},'content' => []},{'type' => 't','content' => '
#PARTLOCAL
T               UsePartLocal

! Default is UsePartLocal=.false. If UsePartLocal is true and the
! run is time accurate, then the blocks selected as "implicit"
! by the criteria defined in #STEPPINGCRITERIA are not used to
! calculate the time step, and all cells are advanced with the
! smaller of the stable and the global time steps.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'IMPLICIT'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UsePointImplicit','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UsePartImplicit','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseFullImplicit','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CflImpl','type' => 'real','default' => '100','min' => '0','if' => '$UsePartImplicit or $UseFullImplicit'},'content' => []},{'name' => 'rule','type' => 'e','attrib' => {'expr' => '$UsePointImplicit + $UsePartImplicit + $UseFullImplicit <= 1'},'content' => [{'type' => 't','content' => '
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
'},{'name' => 'command','type' => 'e','attrib' => {'name' => 'IMPLICITCRITERIA','alias' => 'STEPPINGCRITERIA'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypeImplCrit','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'Time step','default' => 'T','value' => 'dt'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Radial distance','value' => 'r/R'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Test block','value' => 'test'},'content' => []}]},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'rImplicit','type' => 'real','min' => '0','if' => '$TypeImplCrit eq \'R\''},'content' => []},{'type' => 't','content' => '

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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'PARTIMPL','alias' => 'PARTIMPLICIT'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UsePartImplicit2','type' => 'logical','default' => 'F'},'content' => []},{'type' => 't','content' => '

#PARTIMPLICIT
T		UsePartImplicit2

If UsePartImplicit2 is set to true, the explicit scheme is executed in all
blocks before the implicit scheme is applied in the implicit blocks. This way
the fluxes at the explicit/implicit interface are second order accurate,
and the overall part implicit scheme will be fully second order in time.
When this switch is false, the explicit/implicit interface fluxes are only
first order accurate in time.
A potential drawback of the second order scheme is that the explicit scheme 
may crash in the implicit blocks. This could be avoided with a more 
sophisticated implementation. There may also be a slight speed penalty, 
because the explicit scheme is applied in more blocks. 

The default is UsePartImplicit2 = false for now, which is safe and 
backward compatible.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'IMPLSTEP','alias' => 'IMPLICITSTEP'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'ImplCoeff','type' => 'real','default' => '1','min' => '0','max' => '1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseBdf2','type' => 'logical','default' => 'T'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseSourceImpl','type' => 'logical','default' => 'F'},'content' => []},{'type' => 't','content' => '

#IMPLSTEP
0.5			ImplCoeff
F			UseBdf2
F			UseSourceImpl

The ImplCoeff is the beta coefficient in front of the implicit terms
for the two-level implicit scheme.
The UseBdf2 parameter decides if the 3 level BDF2 scheme is used or 
a 2 level scheme.
UseSourceImpl true means that the preconditioner should take point
source terms into account. 

For steady state run the default is the backward Euler scheme, which 
corresponds to ImplCoeff=1.0 and UsedBdf2=F.
For second order time accurate run the default is UseBdf2=T, since
BDF2 is a 3 level second order in time and stable implicit scheme.
In both cases the default value for UseSourceImpl is false.

The default values can be overwritten with #IMPLSTEP, but only
after the #TIMESTEPPING command!
For example one could use the 2-level trapezoid scheme with
ImplCoeff=0.5 and UseBDF2=F as shown in the example above. 
The BDF2 scheme determines the coefficient for the implicit terms itself, 
but ImplCoeff is still used in the first time step and after AMR-s, when
the code switches back to the two-level scheme for one time step.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'IMPLSCHEME','alias' => 'IMPLICITSCHEME'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nOrderImpl','type' => 'integer','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => '1','default' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => '2'},'content' => []}]},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypeFluxImpl','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'Roe','value' => 'Roe/roe/1'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Rusanov','default' => 'T','value' => 'Rusanov/rusanov/2/TVDLF'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Linde','value' => 'Linde/linde/3/HLLEL'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Sokolov','value' => 'Sokolov/sokolov/4/AW'},'content' => []}]},{'type' => 't','content' => '
#IMPLSCHEME
1               nOrderImpl
Rusanov         TypeFluxImpl

! This command defines the scheme used in the implicit part (\'left hand side\').
! The default order is first order. The default scheme is the same as the
! scheme selected for the explicit part. 
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'IMPLCHECK','alias' => 'IMPLICITCHECK'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'RejectStepLevel','type' => 'real','default' => '0.3','min' => '0','max' => '0.9'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'RejectStepFactor','type' => 'real','default' => '0.5','min' => '0','max' => '0.9'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'ReduceStepLevel','type' => 'real','default' => '0.6','min' => '0','max' => '0.9'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'ReduceStepFactor','type' => 'real','default' => '0.9','min' => '0','max' => '1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'IncreaseStepLevel','type' => 'real','default' => '0.8','min' => '0','max' => '1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'IncreaseStepFactor','type' => 'real','default' => '1.05','min' => '1','max' => '2'},'content' => []},{'type' => 't','content' => '

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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'NEWTON'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseConservativeImplicit','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseNewton','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$UseNewton'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseNewMatrix','type' => 'logical','default' => 'T'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'MaxIterNewton','type' => 'integer','default' => '10','min' => '1'},'content' => []}]},{'type' => 't','content' => '
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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'JACOBIAN'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypeJacobian','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'Preconditioned','default' => 'T','value' => 'prec'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'No preconditioning','value' => 'free'},'content' => []}]},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'JacobianEps','type' => 'real','default' => '$doublePrecision ? 1.e-12 : 1.e-6','min' => '0','max' => '1.e-5'},'content' => []},{'type' => 't','content' => '
#JACOBIAN
prec            TypeJacobian (prec, free)
1.E-12          JacobianEps

! The Jacobian matrix is always calculated with a matrix free approach,
! however it can be preconditioned  (\'prec\'), or not (\'free\').  The
! Default value is TypeJacobian=\'prec\'.
! JacobianEps contains the machine round off error for numerical derivatives.
! The default value is 1.E-12 for 8 byte reals and 1.E-6 for 4 byte reals.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'PRECONDITIONER'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypePrecondSide','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'left'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'symmetric','default' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'right'},'content' => []}]},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypePrecond','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'MBILU','default' => 'T'},'content' => []}]},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'GustafssonPar','type' => 'real','default' => '0.5','min' => '0','max' => '1'},'content' => []},{'type' => 't','content' => '
#PRECONDITIONER
symmetric       TypePrecondSide (left, symmetric, right)
MBILU           TypePrecond (MBILU)
0.5             GustafssonPar (0. no modification, 1. full modification)

! Default parameters are shown. Right preconditioning does not affect
! the normalization of the residual. The Gustafsson parameter determines
! how much the MBILU preconditioner is modified. The default 0.5 value
! means a relaxed modification.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'KRYLOV'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypeKrylov','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'gmres','default' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'bicgstab'},'content' => []}]},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypeInitKrylov','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => '0','default' => 'T','value' => 'nul'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'previous','value' => 'old'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'explicit'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'scaled explicit','value' => 'explicit'},'content' => []}]},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'ErrorMaxKrylov','type' => 'real','default' => '0.001','min' => '0','max' => '0.1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'MaxMatvecKrylov','type' => 'integer','default' => '100','min' => '1'},'content' => []},{'type' => 't','content' => '
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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'KRYLOVSIZE'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nKrylovVector','type' => 'integer','default' => 'MaxMatvecKrylov','min' => '1'},'content' => []},{'type' => 't','content' => '
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
'}]}]},{'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'STOPPING CRITERIA'},'content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! STOPPING CRITERIA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

The commands in this group only work in stand alone mode.

'},{'name' => 'command','type' => 'e','attrib' => {'name' => 'STOP','required' => '$_IsStandAlone','if' => '$_IsStandAlone'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'MaxIteration','type' => 'integer','default' => '-1','min' => '-1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'tSimulationMax','type' => 'real','default' => '-1','min' => '-1'},'content' => []},{'type' => 't','content' => '

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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'CPUTIMEMAX','if' => '$_IsStandAlone'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CpuTimeMax','type' => 'real','default' => '-1','min' => '-1'},'content' => []},{'type' => 't','content' => '

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
'},{'name' => 'command','type' => 'e','attrib' => {'name' => 'RESTARTOUTDIR'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'NameRestartOutDir','length' => '100','type' => 'string','default' => 'GM/restartOUT'},'content' => []},{'name' => 'rule','type' => 'e','attrib' => {'expr' => '-d $NameRestartOutDir'},'content' => [{'type' => 't','content' => '
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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'RESTARTOUTFILE'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypeRestartOutFile','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'block'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'one','default' => 'T'},'content' => []}]},{'type' => 't','content' => '

#RESTARTOUTFILE
one			TypeRestartOutFile

This command determines if the restart information is saved as an individual 
file for each block (block) or into a single direct access file containing all
blocks (one).  The default value is \'one\'.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'SAVERESTART'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoSaveRestart','type' => 'logical','default' => 'T'},'content' => []},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$DoSaveRestart'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DnSaveRestart','type' => 'integer','default' => '-1','min' => '-1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DtSaveRestart','type' => 'real','default' => '-1','min' => '-1'},'content' => []}]},{'type' => 't','content' => '
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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'PLOTDIR'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'NamePlotDir','length' => '100','type' => 'string','default' => 'GM/IO2'},'content' => []},{'name' => 'rule','type' => 'e','attrib' => {'expr' => '-d $NamePlotDir'},'content' => [{'type' => 't','content' => '
		Plot directory $NamePlotDir must exist
	'}]},{'type' => 't','content' => '

The NamePlotDir variable contains the name of the directory
where plot files and logfiles are saved relative to the run directory.
The directory should be inside the subdirectory with the name
of the component.

Default value is "GM/IO2".
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'SAVELOGFILE'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoSaveLogfile','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$DoSaveLogfile'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'StringLog','type' => 'strings','min' => '1','max' => '4'},'content' => [{'name' => 'part','type' => 'e','attrib' => {'name' => 'TypeLogVar','type' => 'string','required' => 'T','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'MHD vars. dimensional','default' => 'T','value' => 'MHD'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Flux vars. dimensional','value' => 'FLX'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Raw vars. dimensional','value' => 'RAW'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Set vars. dimensional','value' => 'VAR'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'MHD vars. scaled','default' => 'T','value' => 'mhd'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Flux vars. scaled','value' => 'flx'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Raw vars. scaled','value' => 'raw'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Set vars. scaled','value' => 'var'},'content' => []}]},{'name' => 'part','type' => 'e','attrib' => {'name' => 'TypeTime','multiple' => 'T','type' => 'string','required' => 'F','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'none','exclusive' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'step'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'date'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'time'},'content' => []}]}]},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DnSaveLogfile','type' => 'integer','default' => '1','min' => '-1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DtSaveLogfile','type' => 'real','default' => '-1','min' => '-1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'NameLogVars','length' => '100','type' => 'string','if' => '$TypeLogVar =~ /var/i'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'StringLogRadii','length' => '100','type' => 'strings','min' => '1','if' => '($TypeLogVar=~/flx/i or $NameLogVars=~/flx/i)','max' => '10'},'content' => [{'name' => 'part','type' => 'e','attrib' => {'name' => 'LogRadii','multiple' => 'T','type' => 'real','min' => '$rBody'},'content' => []}]}]},{'type' => 't','content' => '
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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'SATELLITE','if' => '$_IsFirstSession'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nSatellite','type' => 'integer','default' => '0','min' => '0'},'content' => []},{'name' => 'for','type' => 'e','attrib' => {'from' => '1','to' => '$nSatellite'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'StringSatellite','type' => 'strings','min' => '1','max' => '5'},'content' => [{'name' => 'part','type' => 'e','attrib' => {'name' => 'TypeSatelliteVar','type' => 'string','required' => 'T','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'MHD vars. dimensional','default' => 'T','value' => 'MHD'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'All vars. dimensional','value' => 'FUL'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Set vars. dimensional','value' => 'VAR'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'MHD vars. scaled','value' => 'mhd'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'All vars. scaled','value' => 'ful'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Set vars. scaled','value' => 'var'},'content' => []}]},{'name' => 'part','type' => 'e','attrib' => {'name' => 'TypeTrajectory','type' => 'string','required' => 'F','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'file','default' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'equation','value' => 'eqn'},'content' => []}]},{'name' => 'part','type' => 'e','attrib' => {'name' => 'TypeTime','multiple' => 'T','type' => 'string','required' => 'F','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'none','exclusive' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'step'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'date'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'time'},'content' => []}]}]},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DnOutput','type' => 'integer','default' => '1','min' => '-1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DtOutput','type' => 'real','default' => '-1','min' => '-1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'NameTrajectoryFile','length' => '100','type' => 'string'},'content' => []},{'name' => 'rule','type' => 'e','attrib' => {'expr' => '-f $NameTrajectoryFile'},'content' => [{'type' => 't','content' => '
			Trajectory file $NameTrajectoryFile must exist
		'}]},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'NameSatelliteVars','length' => '100','type' => 'string','if' => '$TypeSatelliteVar =~ /\\bvar\\b/i'},'content' => []}]},{'type' => 't','content' => '
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

! The file containing the satellite trajectory should include data in the 
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
'}]},{'type' => 't','content' => '

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


'},{'name' => 'command','type' => 'e','attrib' => {'name' => 'STEADYSTATESATELLITE','if' => '$_IsFirstSession'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'SatelliteTimeStart','type' => 'real','default' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'SatelliteTimeEnd','type' => 'real','default' => '0'},'content' => []},{'type' => 't','content' => '
#STEADYSTATESATELLITE
-86400.0                     SatelliteTimeStart [sec]
86400.0                      SatelliteTimeEnd   [sec]
-3600.0                      SatelliteTimeStart [sec]
3600.0                       SatelliteTimeEnd   [sec]

! In the non-time-accurate mode the numerical simulation result converges
! to a steady-state solution. In the course of this simulation mode, the 
! progress in the iteration number is not associated with an increase in the physical 
! time, and the ultimate solution is a  "snapshot" of the parameter distribution at the 
! time instant set by the  #STARTTIME command. Since time does not run, a satellite 
! position cannot be determined in terms of the simulation time. Instead, the parameters 
! along a cut of the satellite trajectory can be saved on file for a given iteration number. 
! The trajectory points can be naturally parameterized by time, so that the cut can be 
! specified with the choice of the start time, end time, and time interval.
! 
! The command #STEADYSTATESATELLITE is required for a steady-state simulation.
! For each of the satellites, the SatelliteTimeStart is a real value that sets the start of 
! trajectory cut,  while SatelliteTimeEnd sets the end of the trajectory cut. Both are in
! seconds with respect to to the time given in #STARTTIME. A negative value means the time 
! prior to  the #STARTTIME. 

! The DtOutput from the #SATELLITE command specifies the frequency of the points along the 
! satellite trajectory for the non-time-accurate mode, while DnOutput keeps to control the
! iteration number at which the data at the trajectory cut are written to the satellite 
! output file.
!
! For more than one satellite (two satellites in the above given example), the start and 
! end times should be set for all of them. 
'}]},{'type' => 't','content' => '

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

'},{'name' => 'command','type' => 'e','attrib' => {'name' => 'SAVEPLOT'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nPlotFile','type' => 'integer','default' => '0','min' => '0','max' => '100'},'content' => []},{'name' => 'for','type' => 'e','attrib' => {'from' => '1','to' => '$nPlotFile'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'StringPlot','type' => 'strings','min' => '3','max' => '3'},'content' => [{'name' => 'part','type' => 'e','attrib' => {'name' => 'plotform','type' => 'string','required' => 'T','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'TECPLOT','value' => 'tec'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'IDL','value' => 'idl'},'content' => []}]},{'name' => 'part','type' => 'e','attrib' => {'name' => 'plotarea','type' => 'string','required' => 'T','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => '3D','value' => '3d/3d_'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'value' => 'x=0'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'default' => 'T','value' => 'y=0'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'value' => 'z=0'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'value' => 'sph'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'value' => 'los'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'value' => 'lin'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'value' => 'cut'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'value' => 'slc','if' => '$plotform =~ /\\btec\\b/'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'value' => 'dpl','if' => '$plotform =~ /\\btec\\b/'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'value' => 'blk','if' => '$plotform =~ /\\btec\\b/'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => '1D','value' => '1d/1d_','if' => '$plotform =~ /\\btec\\b/'},'content' => []}]},{'name' => 'part','type' => 'e','attrib' => {'name' => 'plotvar','type' => 'string','required' => 'T','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'MHD vars. dimensional','value' => 'MHD'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'All vars. dimensional','value' => 'FUL'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Raw vars. dimensional','value' => 'RAW'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Ray tracing vars. dim.','value' => 'RAY'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Flux vars. dimensional','value' => 'FLX'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Solar vars. dimensional','value' => 'SOL'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Position vars. dimensional','value' => 'POS','if' => '$plotarea eq \'lin\''},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Select dimensional vars.','value' => 'VAR'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Basic Block vars. dimensional','value' => 'BBK'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'MHD vars. scaled','value' => 'mhd'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'All vars. scaled','value' => 'ful'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Raw vars. scaled','value' => 'raw'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Ray tracing vars. scaled','value' => 'ray'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Flux vars. scaled','value' => 'flx'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Solar vars. scaled','value' => 'sol'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Position vars. scaled','value' => 'pos','if' => '$plotarea eq \'lin\''},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Select scaled vars.','value' => 'var'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Basic Block vars. scaled','value' => 'bbk'},'content' => []}]}]},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DnSavePlot','type' => 'integer','min' => '-1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DtSavePlot','type' => 'real','min' => '-1'},'content' => []},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$plotarea =~ /\\bdpl|cut|slc\\b/'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'xMinCut','type' => 'real'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'xMaxCut','type' => 'real','min' => '$xMinCut'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'yMinCut','type' => 'real'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'yMaxCut','type' => 'real','min' => '$yMinCut'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'zMinCut','type' => 'real'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'zMaxCut','type' => 'real','min' => '$zMinCut'},'content' => []}]},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$plotarea =~ /\\bslc\\b/'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'xPoint','type' => 'real'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'yPoint','type' => 'real'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'zPoint','type' => 'real'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'xNormal','type' => 'real'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'yNormal','type' => 'real'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'zNormal','type' => 'real'},'content' => []}]},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$plotarea =~ /\\bblk\\b/'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'xPoint','type' => 'real'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'yPoint','type' => 'real'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'zPoint','type' => 'real'},'content' => []}]},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Radius','type' => 'real','default' => '10','min' => '0','if' => '$plotarea =~ /\\bsph\\b/'},'content' => []},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$plotarea =~ /\\blos\\b/'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'ObsPosX','type' => 'real','default' => '-215'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'ObsPosY','type' => 'real','default' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'ObsPosZ','type' => 'real','default' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'OffsetAngle','type' => 'real','default' => '0','min' => '-89','max' => '89'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'rSizeImage','type' => 'real','default' => '32','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'xOffset','type' => 'real','default' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'yOffset','type' => 'real','default' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'rOccult','type' => 'real','default' => '2','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'MuLimbDarkening','type' => 'real','default' => '0.5','min' => '0','max' => '1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nPix','type' => 'integer','default' => '200','min' => '2'},'content' => []}]},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$plotarea =~ /\\blin\\b/'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'NameLine','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'Advected B','value' => 'A'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'B','default' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'U'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'J'},'content' => []}]},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'IsSingleLine','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nLine','type' => 'integer','default' => '1','min' => '1','max' => '20'},'content' => []},{'name' => 'for','type' => 'e','attrib' => {'from' => '1','to' => '$nLine'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'xStartLine','type' => 'real'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'yStartLine','type' => 'real'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'zStartLine','type' => 'real'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'IsParallel','type' => 'logical'},'content' => []}]}]},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DxSavePlot','type' => 'real','default' => '-1.0','min' => '-1.0','if' => '($plotform=~/\\bidl\\b/ and $plotarea!~/\\b(sph|los|lin)\\b/)'},'content' => []},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$plotvar =~ /\\bvar\\b/i'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'NameVars','length' => '100','type' => 'string'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'NamePars','length' => '100','type' => 'string'},'content' => []}]}]},{'type' => 't','content' => '
#SAVEPLOT
10			nPlotfile
3d  MHD tec		StringPlot ! 3d plot with MHD data
100			DnSavePlot
-1.			DtSavePlot
y=0 VAR idl		StringPlot ! y=0 plane plot with listed variables
-1			DnSavePlot
100.			DtSavePlot
0.25			DxSavePlot ! resolution (for IDL plots)
jx jy jz		NameVars
g rbody			NamePars
cut ray idl		StringPlot ! 3D cut IDL (ONLY!) plot with raytrace info
1			DnSavePlot
-1.			DtSavePlot
-10.			xMinCut    
10.			xMaxCut    
-10.			yMinCut    
10.			yMaxCut    
-10.			zMinCut    
10.			zMaxCut    
-1.			DxSavePlot ! unstructured grid (for IDL plots)
sph flx idl		StringPlot ! spherical plot
-1			DnSavePlot
100.			DtSavePlot
4.			Radius     ! of spherical surface
los sol idl             StringPlot ! line of sight plot
-1			DnSavePlot
100.			DtSavePlot
-215.			ObsPosX
0.			ObsPosY
0.			ObsPosZ
0.                      OffsetAngle
32.			rSizeImage
0.			xOffset
0.			yOffset
3.			rOccult
0.5			MuLimbDarkening
300			nPix
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
dpl MHD tec		StringPlot  ! dipole slice Tecplot (ONLY!) plot
-1			DnSavePlot
10.			DtSavePlot
-10.			xMinCut
 10.			xMaxCut
-10.			yMinCut
 10.			yMaxCut
-10.			zMinCut
 10.			zMaxCut
slc MHD tec		StringPlot  ! general slice Tecplot (ONLY!) plot
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
blk MHD tec		StringPlot  ! general block Tecplot (ONLY!) plot
-1			DnSavePlot
10.			DtSavePlot
 5.			xPoint
 1.			yPoint
 1.			zPoint
1d  BBK tec		StringPlot ! 1d plot with BLK data, Tecplot (ONLY!)
1000			DnSavePlot
-1.			DtSavePlot

Default is nPlotFile=0. \\\\

\\noindent
StringPlot must contain the following 3 parts in arbitrary order
\\begin{verbatim}
plotarea plotvar plotform

plotarea = \'3d\' , \'x=0\', \'y=0\', \'z=0\', \'cut\', \'dpl\', \'slc\', \'sph\', \'los\', \'lin\', \'blk\', \'1d\'
plotvar  = \'mhd\', \'ful\', \'raw\', \'ray\', \'flx\', \'sol\', \'pos\', \'var\', \'bbk\' - normalized
plotvar  = \'MHD\', \'FUL\', \'RAW\', \'RAY\', \'FLX\', \'SOL\', \'pos\', \'VAR\', \'BBK\' - dimensional
plotform = \'tec\', \'idl\'
\\end{verbatim}
NOTES: The plotvar option \'sol\' is only valid for plotarea \'los\' and
       the plotvar option \'pos\' is only valid for plotarea \'lin\'.\\\\

\\noindent
The plotarea string defines the 1, 2, or 3D volume of the plotting area:
\\begin{verbatim}
1d    - full 1D, one point per block for all blocks
x=0   - full x=0 plane: average for symmetry plane
y=0   - full y=0 plane: average for symmetry plane
z=0   - full z=0 plane: average for symmetry plane
3d    - full 3D volume
cut   - 3D, 2D or 1D rectangular cut (IDL)/ a 2D rectangular cut (Tecplot)
dpl   - cut at dipole \'equator\', uses PLOTRANGE to clip plot
slc   - 2D slice defined with a point and normal, uses PLOTRANGE to clip plot
sph   - spherical surface cut at the given radius
los   - line of sight integrated plot
lin   - one dimensional plot along a field or stream line
blk   - 3D single block cell centered data, block specified point location
\\end{verbatim}
The plotvar string defines the plot variables and the equation parameters.
It also controls whether or not the variables will be plotted in dimensional
values or as non-dimensional values:
\\begin{verbatim}
 ALL CAPS  - dimensional
 all lower - dimensionless

 \'mhd\' - vars: rho Ux Uy Uz E Bx By Bz P Jx Jy Jz
         pars: g eta
 \'ful\' - vars: rho Ux Uy Uz E Bx By Bz B1x B1y B1z P Jx Jy Jz
         pars: g eta
 \'raw\' - vars: rho rhoUx rhoUy rhoUz E Bx By Bz P b1x b1y b1z divb
         pars: g eta
 \'ray\' - vars: bx by bz theta1 phi1 theta2 phi2 status blk
         pars: R_ray
 \'flx\' - vars: rho rhoUr Br jr pvecr
         pars: g eta
 \'var\' - vars: READ FROM PARAMETER FILE
         pars: READ FROM PARAMETER FILE
 \'sol\' - vars: wl pb
         pars: mu
 \'bbk\' - vars: dx pe blk blkall
         pars: 
\\end{verbatim}
The plot_string is always followed by the plotting frequencies
DnSavePlot and DtSavePlot.\\\\

\\noindent
Depending on StringPlot, further information is read from the parameter file
in this order:
\\begin{verbatim}
 PlotRange		if plotarea is \'cut\', \'dpl\', or \'slc\'
 Point			if plotarea is \'slc\', or \'blk\'
 Normal			if plotarea is \'slc\'
 DxSavePlot		if plotform is \'idl\' and plotarea is not sph, ion, los
 Radius                 if plotarea is \'sph\'
 NameVars		if plotform is \'var\'
 NamePars		if plotform is \'var\'
\\end{verbatim}
The PlotRange is described by 6 coordinates. 
For IDL plots, if the width in one or two 
dimensions is less than the smallest cell size within the plotarea, 
then the plot file will be 2 or 1 dimensional, respectively.
If the range is thin but symmetric about one of the x=0, y=0, or z=0 planes, 
data will be averaged in the postprocessing.\\\\

For Tecplot (tec) file and \'cut\', Plotrange is read but 
only 1 dimension is used.  
Cuts are entire x, y, or z = constant planes (2D only, 1D or 3D cuts are not
implemented.  For x=constant, for example, the y and z ranges 
do not matter as long at they are "wider" than the x range.  The slice will be 
located at the average of the two x ranges.  So, for example to save a plot in
a x=-5 constant plane cut in tec. The following would work for the plot range:
\\begin{verbatim}
 -5.01			xMinCut
 -4.99			xMaxCut
 -10.			yMinCut
  10.			yMaxCut
 -10.			zMinCut
  10.			zMaxCut
\\end{verbatim}
The \'dpl\' and \'slc\' Tecplot plots use Plotrange like the IDL plots, and will
clip the cut plane when it exits the defined box.

\\noindent
Point is described by the coordinate of any point on the cut plane, often the origin,
or inside of a 3D block. Normal is the coordinate of a vector normal to the plane.  
If the normal in any given coordinate direction is less than 0.01, 
then no cuts are computed for cell edges parallel to that coordinate direction.
For example, the following would result in only computing cuts on cell 
edges parallel to the Z axis.
\\begin{verbatim}
  0.0			xNormal
  0.0			yNormal
  1.0			zNormal
\\end{verbatim}

\\noindent
Possible values for DxSavePlot (for IDL files):
\\begin{verbatim}
  any positive value	- fixed resolution
  0.	- fixed resolution based on the smallest cell in the plotting area
 -1.	- unstructured grid will be produced by PostIDL.exe
\\end{verbatim}
Radius is the radius of the spherical cut for plotarea=\'sph\'

The line-of-sight (plotarea \'los\') plots calculate integrals along the
lines of sight of some quantity and create a 2D Cartesian square
shaped grid of the integrated values. Only the circle enclosed in the
square is actually calculated and the corners are filled in with
zeros.  The image plane always contains the origin of the
computational domain (usually the center of the Sun).  By default the
image plane is orthogonal to the observers position relative to the
origin. The image plane can be rotated around the Z axis with an
offset angle. By default the center of the image is the observer
projected onto the image plane, but the center of the image can be
offset.  Since the central object (the Sun) contains extremely large
values, an occultational disk is used to block the lines of sight
going through the Sun.  The variables which control the direction of
the lines of sight and the grid position and resolution are the
following:
\\begin{verbatim}
 ObsPosX,ObsPosY,ObsPosZ- the position of the observer
 OffsetAngle            - the angle between the image plane normal direction
                          and the Sun\'s direction from the observer position
 rSizeImage             - the radius of the LOS image
 xOffset, yOffset       - offset relative to the observer projected onto 
                          the image plane 
 rOccult                - the radius of the occulting disk
 MuLimbDarkening        - the limb darkening parameter for the \'wl\' 
                          (white light) and \'pb\' (polarization brightness) 
                          plot variables.
 nPix                   - the number of pixels in each direction
\\end{verbatim}

\\noindent
The possible values for NameVars with plotarea \'los\' 
are listed in subroutine set_plotvar_los in write_plot_los.f90. \\\\

\\noindent
The possible values for NameVars for other plot areas
are listed in subroutine set_plotvar in write_plot_common.f90.\\\\

\\noindent
The possible values for NamePars are listed in subroutine 
set_eqpar in write_plot_common.f90\\\\

A plot file is produced by each processor.  This file is ASCII in \'tec\'
format and can be either binary or ASCII in \'idl\' format as chosen under
the #SAVEBINARY flag.  The name of the files are
\\begin{verbatim}
 IO2/plotarea_plotvar_plotnumber_timeinfo_PEnumber.extension
\\end{verbatim}
where extension is \'tec\' for the TEC and \'idl\' for the IDL file formats.
The plotnumber goes from 1 to nplot in the order of the files in PARAM.in.
The \'timeinfo\' contains simulation time as hours-minutes-seconds
(for time accurate runs only), and time step number.
Spherical plot area \'sph\' creates two files per processor starting with
\'spN\' and \'spS\' for the northern and southern hemispheres, respectively.  

After all processors wrote their plot files, processor 0 writes a small 
ASCII header file named as
\\begin{verbatim}
 IO2/plotarea_plotvar_plotnumber_timestep.headextension
\\end{verbatim}
where headextension is:
\\begin{verbatim}
           \'T\' for TEC file format
           \'S\' for TEC and plot_area \'sph\' 
           \'h\' for IDL file format       
\\end{verbatim}

\\noindent
The line of sight integration produces TecPlot and IDL files directly:
\\begin{verbatim}
 IO2/los_plotvar_plotnumber_timestep.extension
\\end{verbatim}
where extension is \'dat\' for TecPlot and \'out\' for IDL file formats.
The IDL output from line of sight integration is always in ASCII format.

'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'SAVEBINARY'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoSaveBinary','type' => 'logical','default' => 'T'},'content' => []},{'type' => 't','content' => '
#SAVEBINARY
T			DoSaveBinary   used only for \'idl\' plot file

! Default is .true. Saves unformatted IO2/*.idl files if true. 
! This is the recommended method, because it is fast and accurate.
! The only advantage of saving IO2/*.idl in formatted text files is
! that it can be processed on another machine or with a different 
! (lower) precision. For example PostIDL.exe may be compiled with 
! single precision to make IO2/*.out files smaller, while BATSRUS.exe is 
! compiled in double precision to make results more accurate.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'SAVEPLOTSAMR'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoSavePlotsAmr','type' => 'logical','default' => 'F'},'content' => []},{'type' => 't','content' => '
#SAVEPLOTSAMR
F			DoSavePlotsAmr

! Save plots before each AMR. Default is DoSavePlotsAMR=.false.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'FLUSH'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoFlush','type' => 'logical','default' => 'T'},'content' => []},{'type' => 't','content' => '

#FLUSH
F			DoFlush

If the DoFlush variable is true, the output is flushed when
subroutine ModUtility::flush_unit is called. This is used in the 
log and satellite files. The flush is useful to see the output immediately, 
and to avoid truncated files when the code crashes,
but on some systems the flush may be very slow. 

The default is to flush the output, i.e. DoFlush=T.
'}]}]},{'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'AMR PARAMETERS'},'content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  AMR PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'name' => 'command','type' => 'e','attrib' => {'name' => 'AMRINIT','if' => '$_IsFirstSession'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'InitialRefineType','length' => '100','type' => 'string','default' => 'none'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'InitialRefineLevel','type' => 'integer','default' => '4','min' => '0'},'content' => []},{'type' => 't','content' => '
#AMRINIT
3Dbodyfocus		TypeRefineInit
4			nRefineLevelInit

This command defines the geometry based refinement type and
the number of initial refinement levels. The exact meaning
of the various refinement types can be figured out from the
source code in specify_refinement.f90 and/or the 
subroutine user_specify_derinement in the user module. 

This is an obscure way of defining the grid, because the meaning
of the various names keep changing from version to version.
It is much better to use the #GRIDRESOLUTION and #GRIDLEVEL commands
that define the grid structure in geometric terms.

The default is TypeRefineInit=\'none\' that results in no refinement
of the base grid.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'AMRINITPHYSICS','if' => '$_IsFirstSession'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nRefineLevelIC','type' => 'integer','default' => '0','min' => '0'},'content' => []},{'type' => 't','content' => '
#AMRINITPHYSICS
3			nRefineLevelIC

! Defines number of physics (initial condition) based AMR-s AFTER the 
! geometry based initial AMR-s defined by #AMRINIT were done.
! Only useful if the initial condition has a non-trivial analytic form.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'GRIDRESOLUTION','multiple' => 'T','alias' => 'GRIDLEVEL'},'content' => [{'name' => 'set','type' => 'e','attrib' => {'name' => 'RotateArea','value' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Resolution','type' => 'real','min' => '0','if' => '$_command =~ /RESOLUTION/'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nLevel','type' => 'integer','min' => '0','if' => '$_command =~ /LEVEL/'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'StringArea','type' => 'strings','min' => '1','max' => '2'},'content' => [{'name' => 'part','type' => 'e','attrib' => {'name' => 'NameArea','type' => 'string','required' => 'T','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'initial','default' => 'T','value' => 'init/initial'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'all'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'box'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'brick'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'brick0'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'sphere'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'sphere0'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'shell'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'shell0'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'cylinderx'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'cylinderx0'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'cylindery'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'cylindery0'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'cylinderz'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'cylinderz0'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'ringx'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'ringx0'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'ringy'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'ringy0'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'ringz'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'ringz0'},'content' => []}]},{'name' => 'part','type' => 'e','attrib' => {'name' => 'RotateArea','type' => 'string','required' => 'F','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'rotated'},'content' => []}]}]},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$NameArea !~ /box|all|init|0/'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'xCenter','type' => 'real','default' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'yCenter','type' => 'real','default' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'zCenter','type' => 'real','default' => '0'},'content' => []}]},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$NameArea =~ /box/'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'xMinBox','type' => 'real'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'yMinBox','type' => 'real'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'zMinBox','type' => 'real'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'xMaxBox','type' => 'real'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'yMaxBox','type' => 'real'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'zMaxBox','type' => 'real'},'content' => []}]},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$NameArea =~ /brick/i'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'xSizeBrick','type' => 'real','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'ySizeBrick','type' => 'real','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'zSizeBrick','type' => 'real','min' => '0'},'content' => []}]},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$NameArea =~ /sphere/i'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Radius','type' => 'real','min' => '0'},'content' => []}]},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$NameArea =~ /shell/i'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Radius1','type' => 'real','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Radius2','type' => 'real','min' => '0'},'content' => []}]},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$NameArea =~ /cylinder/i'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'LengthCylinder','type' => 'real','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Radius','type' => 'real','min' => '0'},'content' => []}]},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$NameArea =~ /ring/i'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'HeightRing','type' => 'real','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Radius1','type' => 'real','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Radius2','type' => 'real','min' => '0'},'content' => []}]},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$RotateArea =~ /rotated/i'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'xRotate','type' => 'real','default' => '0','min' => '-360','max' => '360'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'yRotate','type' => 'real','default' => '0','min' => '-360','max' => '360'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'zRotate','type' => 'real','default' => '0','min' => '-360','max' => '360'},'content' => []}]},{'type' => 't','content' => '

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

If NameArea is set to \'initial\' or \'init\', the highest initial resolution
or level is set by the command. This is similar to the #AMRINIT command,
but there is no \'TypeRefineInit\' parameter, because all refinement areas 
are defined with the additional #GRIDLEVEL and #GRIDRESOLUTION commands.
Do not use #AMRINIT if the #GRIDLEVEL or #GRIDRESOLUTION command is used 
with the \'initial\' area! The #AMRINIT command may be combined with
the #GRIDLEVEL and #GRIDRESOLUTION commands (these add additional areas
of refinement to the areas defined by TypeRefineInit), but in general 
the use of #AMRINIT is not recommended, because the refinement definitions
are defined in the source code which is difficult to understand and may
change from version to version.

For other values of NameArea, the command specifies the shape of the area. 
where the blocks are to be refined. If the desired grid resolution is finer
than the initial resolution, then initially the grid will be refined
to the initial resolution only, but the area will be further refined 
in subsequent pre-specified adaptive mesh refinements (AMRs) during the run 
(see the #AMR command). Once the resolution reaches the
desired level, the AMR-s will not do further refinement. If a grid block
is covered by more than one areas, the area with the finest resolution
determines the desired grid resolution.

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

If the area name contains the number \'0\', the center is taken to be at the 
origin. If the word \'rotated\' is added, the area can be rotated by 3 angles 
around the X, Y and Z axes in this order.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'AMRLEVELS'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'MinBlockLevel','type' => 'integer','default' => '0','min' => '-1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'MaxBlockLevel','type' => 'integer','default' => '99','min' => '-1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoFixBodyLevel','type' => 'logical','default' => 'F'},'content' => []},{'type' => 't','content' => '
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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'AMRRESOLUTION'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DxCellMin','type' => 'real','default' => '0','min' => '-1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DxCellMax','type' => 'real','default' => '99999','min' => '-1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoFixBodyLevel','type' => 'logical','default' => 'F'},'content' => []},{'type' => 't','content' => '
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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'AMR'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DnRefine','type' => 'integer','default' => '-1','min' => '-1'},'content' => []},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$DnRefine>0'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoAutoRefine','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$DoAutoRefine'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'PercentCoarsen','type' => 'real','default' => '20','min' => '0','max' => '100'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'PercentRefine','type' => 'real','default' => '20','min' => '0','max' => '100'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'MaxTotalBlocks','type' => 'integer','default' => '99999','min' => '1'},'content' => []}]}]},{'type' => 't','content' => '
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
! based on the TypeRefineInit parameter given in the #AMRINIT command
! and/or the areas and resolutions defined by the 
! #GRIDLEVEL and #GRIDRESOLUTION commands.
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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'AMRCRITERIA'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nRefineCrit','type' => 'integer','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => '0'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => '1'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => '2'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => '3','default' => 'T'},'content' => []}]},{'name' => 'for','type' => 'e','attrib' => {'from' => '1','to' => '$nRefineCrit'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypeRefine','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'grad T','value' => 'gradt/gradT'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'grad P','value' => 'gradp/gradP'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'grad log(Rho)','value' => 'gradlogrho'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'grad log(p)','value' => 'gradlogP/gradlogp'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'grad E','value' => 'gradE'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'curl U','value' => 'curlV/curlv/curlU/curlu'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'curl B','value' => 'curlB/curlb'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'div U','value' => 'divU/divu/divV/divv'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'divB','value' => 'divb/divB'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'vAlfven','value' => 'Valfven/vAlfven/valfven'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'heliospheric beta','value' => 'heliobeta'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'flux'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'heliospheric current sheet','value' => 'heliocurrentsheet'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'rCurrents','value' => 'rcurrents/Rcurrents'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Transient','value' => 'transient/Transient'},'content' => []}]},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$TypeRefine =~ /transient/i'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypeTransient','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'P_dot','value' => 'p_dot/P_dot'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'T_dot','value' => 't_dot/T_dot'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Rho_dot','value' => 'rho_dot/Rho_dot','default' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'RhoU_dot','value' => 'RhoU_dot/rhou_dot'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Rho_2nd_1','value' => 'Rho_2nd_1/rho_2nd_1'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Rho_2nd_2','value' => 'Rho_2nd_2/rho_2nd_2'},'content' => []}]},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseSunEarth','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$UseSunEarth'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'xEarth','type' => 'real'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'yEarth','type' => 'real'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'zEarth','type' => 'real'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'InvD2Ray','type' => 'real'},'content' => []}]}]}]},{'type' => 't','content' => '
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

The default values depend on problem_type.

At most three criteria can be given. If nRefineCrit is set to zero,
the blocks are not ordered. This can be used to refine or coarsen
all the blocks limited by the minimum and maximum levels only
(see commands #AMRLEVELS and #AMRRESOLUTION). If nRefineCrit is 1, 2, or 3
then the criteria can be chosen from the following list:
\\begin{verbatim}
  \'gradT\'		- gradient of temperature
  \'gradP\'		- gradient of pressure
  \'gradlogrho\'		- gradient of log(rho)
  \'gradlogP\'		- gradient of log(P)
  \'gradE\'		- gradient of electric field magnitude
  \'curlV\',\'curlU\' 	- magnitude of curl of velocity
  \'curlB\'		- magnitude of current
  \'divU\', \'divV\'	- divergence of velocity
  \'divB\'		- div B
  \'vAlfven\',\'Valfven\'	- Alfven speed
  \'heliobeta\' 		- special function for heliosphere $R^2 B^2/rho$
  \'flux\'		- radial mass flux
  \'heliocurrentsheet\'	- refinement in the currentsheet of the heliosphere
  \'Rcurrents\'		- refinement near Rcurrents value
\\end{verbatim}
All the names can also be spelled with all small case letters.\\\\

\\noindent
The possible choices for TypeTransient:
\\begin{verbatim}
  \'P_dot\' (same as \'p_dot\')
  \'T_dot\' (same as \'t_dot\')
  \'Rho_dot\' (same as \'rho_dot\')
  \'RhoU_dot\' (same as \'rhou_dot\')
  \'B_dot\' (same as \'b_dot\')
  \'Rho_2nd_1\' (same as \'rho_2nd_1\')
  \'Rho_2nd_2\' (same as \'rho_2nd_2\')
\\end{verbatim}

Also, (xEarth,yEarth,zEarth) are the coordinates of the Earth. InvD2Ray is
a factor that defines how close to the ray Sun-Earth to refine the grid.
Note that the AMR occurs in a cylinder around the ray.
Example for InvD2Ray =
\\begin{verbatim}
   1 - refine_profile = 0.3679 at distance Rsun/10 from the ray
   2 - refine_profile = 0.0183 at distance Rsun/10 from the ray
   3 - refine_profile = 0.0001 at distance Rsun/10 from the ray
\\end{verbatim}
'}]}]},{'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'SCHEME PARAMETERS'},'content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  SCHEME PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'name' => 'command','type' => 'e','attrib' => {'name' => 'SCHEME'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nOrder','type' => 'integer','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => '1'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => '2','default' => 'T'},'content' => []}]},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypeFlux','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'Roe','value' => 'Roe/roe/1'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Rusanov','default' => 'T','value' => 'Rusanov/rusanov/2/TVDLF'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Linde','value' => 'Linde/linde/3/HLLEL'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Sokolov','value' => 'Sokolov/sokolov/4/AW'},'content' => []}]},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$nOrder == 2'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypeLimiter','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'minmod','default' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'mc'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'mc3'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'beta'},'content' => []}]},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'LimiterBeta','type' => 'real','default' => '1.2','min' => '1','if' => '$TypeLimiter ne \'minmod\'','max' => '2'},'content' => []}]},{'type' => 't','content' => '
#SCHEME
2			nOrder (1 or 2)
Rusanov			TypeFlux
mc3			TypeLimiter ! Only for nOrder=2
1.2			LimiterBeta ! Only if TypeLimiter is NOT \'minmod\'

The nOrder parameter determines the spatial and temporal accuracy of
the scheme. The spatially first order scheme uses a one-stage time integration,
while the spatially second order scheme either uses a two-stage Runge-Kutta 
explicit or a BDF2 three-level implicit time discretization. 

\\noindent
Possible values for TypeFlux:
\\begin{verbatim}
 \'Rusanov\'     - Rusanov or Lax-Friedrichs flux     
 \'Linde        - Linde\'s HLLEL flux                   
 \'Sokolov\'     - Sokolov\'s Local Artificial Wind flux 
 \'Roe\'         - Roe\'s approximate Riemann flux       
\\end{verbatim}
Possible values for TypeLimiter for spatially second order schemes:
\\begin{verbatim}
 \'minmod\'      - minmod limiter is the most robust and diffusive limiter
 \'mc\'          - monotonized central limiter with a beta parameter
 \'mc3\'         - third order (in the linear sense) limiter with beta parameter
 \'beta\'        - beta limiter is less robust than the mc limiter for 
                 the same beta value
\\end{verbatim}
Possible values for LimiterBeta (for limiters othen than minmod)
are between 1.0 and 2.0 : 
\\begin{verbatim}
  LimiterBeta = 1.0 is the same as the minmod limiter
  LimiterBeta = 2.0 for the beta limiter is the same as the superbee limiter
  LimiterBeta = 1.5 is a typical value for the mc/mc3 limiters
  LimiterBeta = 1.2 is the recommended value for the beta limiter
\\end{verbatim}
The default is the second order Rusanov scheme with the minmod limiter.

'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'NONCONSERVATIVE'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseNonConservative','type' => 'logical','default' => 'T'},'content' => []},{'type' => 't','content' => '
#NONCONSERVATIVE
T		UseNonConservative

For Earth the default is using non-conservative equations 
(close to the body).
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'CONSERVATIVECRITERIA'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nConservCrit','type' => 'integer','default' => '1','min' => '0','max' => '3'},'content' => []},{'name' => 'for','type' => 'e','attrib' => {'from' => '1','to' => '$nConservCrit'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypeConservCrit','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'radius','default' => 'T','value' => 'r/R/radius/Radius'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'parabola','value' => 'parabola/paraboloid'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'p','value' => 'p/P'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'grad P','value' => 'gradp/GradP'},'content' => []}]},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'rConserv','type' => 'real','default' => '6','min' => '$rBody','if' => '$TypeConservCrit =~ /^r|radius$/i'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'xParabolaConserv','type' => 'real','default' => '6','min' => '0','if' => '$TypeConservCrit =~ /^parabol/i'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'yParabolaConserv','type' => 'real','default' => '36','min' => '0','if' => '$TypeConservCrit =~ /^parabol/i'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'pCoeffConserv','type' => 'real','default' => '0.05','min' => '0','if' => '$TypeConservCrit =~ /^p$/i'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'GradPCoeffConserv','type' => 'real','default' => '0.1','min' => '0','if' => '$TypeConservCrit =~ /gradp/i'},'content' => []}]},{'type' => 't','content' => '

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

Select the parts of the grid where the conservative vs. non-conservative
schemes are applied. The number of criteria is arbitrary, although 
there is no point applying the same criterion more than once.

If no criteria is used, the whole domain will use conservative or
non-conservative equations depending on UseNonConservative set in
command #NONCONSERVATIVE.

The physics based conservative criteria (\'p\' and \'GradP\')
select cells which use the non-conservative scheme if ALL of them are true:
\\begin{verbatim}
 \'p\'      - the pressure is smaller than fraction pCoeffConserv of the energy
 \'GradP\'  - the relative gradient of pressure is less than GradPCoeffConserv
\\end{verbatim}
 The geometry based criteria are applied after the physics based criteria 
 (if any) and they select the non-conservative scheme if ANY of them is true:
\\begin{verbatim}
 \'r\'        - radial distance of the cell is less than rConserv
 \'parabola\' - x less than xParabolaConserv - (y**2+z**2)/yParabolaConserv
\\end{verbatim}
Default values are nConservCrit = 1 with TypeConservCrit = \'r\'
and rConserv=2*rBody, where rBody has a problem dependent default.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'UPDATECHECK'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseUpdateCheck','type' => 'logical','default' => 'T'},'content' => []},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$UseUpdateCheck'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'RhoMinPercent','type' => 'real','default' => '40','min' => '0','max' => '100'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'RhoMaxPercent','type' => 'real','default' => '400','min' => '100'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'pMinPercent','type' => 'real','default' => '40','min' => '0','max' => '100'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'pMaxPercent','type' => 'real','default' => '400','min' => '100'},'content' => []}]},{'type' => 't','content' => '
#UPDATECHECK
T			UseUpdateCheck
40.			RhoMinPercent
400.			RhoMaxPercent
40.			pMinPercent
400.			pMaxPercent

Default values are shown.  This will adjust the timestep so that
density and pressure cannot change by more than the given percentages
in a single timestep.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'MULTISPECIES'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoReplaceDensity','type' => 'logical','default' => 'T'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'SpeciesPercentCheck','type' => 'real','default' => '1.0','min' => '0','max' => '100'},'content' => []},{'type' => 't','content' => '
#MULTISPECIES
T			DoReplaceDensity
1.0			SpeciesPercentCheck

This command is only useful for multispecies equations.
If the DoReplaceDensity is true, the total density is replaced with
the sum of the species densities. The SpeciesPercentCheck parameter
determines if a certain species density should or should not be
checked for large changes. If SpeciesPercentCheck is 0, all species
are checked, if it is 1, then only species with densities reaching
or exceeding 1 per cent are checked for large changes (see the
#UPDATECHECK command).

Default values are shown.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'PROLONGATION'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nOrderProlong','type' => 'integer','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => '1','default' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => '2'},'content' => []}]},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypeProlong','type' => 'string','input' => 'select','if' => '$nOrderProlong==2'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'left-right','default' => 'T','value' => 'lr'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'central'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'minmod'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'left-right extrapolate','value' => 'lr2'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'central    extrapolate','value' => 'central2'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'minmod     extrapolate','value' => 'minmod2'},'content' => []}]},{'type' => 't','content' => '
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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'MESSAGEPASS','alias' => 'OPTIMIZE'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypeMessagePass','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'm_p_cell FACES ONLY','default' => 'T','value' => 'allopt'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'm_p_cell','value' => 'all'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'm_p_dir FACES ONLY','value' => 'opt'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'm_p_dir group by directions','value' => 'dir'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'm_p_dir group by faces     ','value' => 'face'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'm_p_dir group by kind and face','value' => 'min'},'content' => []}]},{'type' => 't','content' => '
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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'RESCHANGE','alias' => 'RESOLUTIONCHANGE'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseAccurateReschange','type' => 'logical','default' => 'F'},'content' => []},{'type' => 't','content' => '
#RESCHANGE
T		UseAccurateResChange

! For UseAccurateResChange=T a second order accurate, upwind and oscillation 
! free scheme is used at the resolution changes. 
!
! Default value is false, which results in first order 
! prolongation and restriction operators at the resolution changes.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'TVDRESCHANGE'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseTvdReschange','type' => 'logical','default' => 'F'},'content' => []},{'type' => 't','content' => '
#TVDRESCHANGE
T		UseTvdResChange

! For UseTvdResChange=T an almost second order and partially downwinded 
! TVD limited scheme is used at the resolution changes. 
!
! Default value is false, which results in first order 
! prolongation and restriction operators at the resolution changes.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'BORIS'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseBorisCorrection','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'BorisClightFactor','type' => 'real','default' => '1','min' => '0','if' => '$UseBorisCorrection','max' => '1'},'content' => []},{'type' => 't','content' => '
#BORIS
T			UseBorisCorrection
1.0			BorisClightFactor !Only if UseBorisCorrection is true

! Default is boris_correction=.false.
! Use semi-relativistic MHD equations with speed of light reduced by
! the BorisClightFactor. Set BorisClightFactor=1.0 for true semi-relativistic
! MHD. Gives the same steady state as normal MHD analytically, but there
! can be differences due to discretization errors. 
! You can use either Boris or BorisSimple but not both. 
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'BORISSIMPLE'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseBorisSimple','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'BorisClightFactor','type' => 'real','default' => '1','min' => '0','if' => '$UseBorisSimple','max' => '1'},'content' => []},{'type' => 't','content' => '
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
		At least one of the options should be true.
	'}]},{'name' => 'rule','type' => 'e','attrib' => {'expr' => 'not($UseProjection and ($UseDivbSource or $UseDivbDiffusion or $UseConstrainB))'},'content' => [{'type' => 't','content' => '
		If UseProjection is true, all others should be false.
	'}]},{'name' => 'rule','type' => 'e','attrib' => {'expr' => 'not($UseConstrainB and ($UseDivbSource or $UseDivbDiffusion or $UseProjection))'},'content' => [{'type' => 't','content' => '
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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'DIVBSOURCE'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseB0Source','type' => 'logical','default' => 'T'},'content' => []},{'type' => 't','content' => '
#DIVBSOURCE
T			UseB0Source

! Add extra source terms related to the non-zero divergence and curl of B0.
! Default is true.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'PROJECTION'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypeProjectIter','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'Conjugate Gradients','default' => 'T','value' => 'cg'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'BiCGSTAB','value' => 'bicgstab'},'content' => []}]},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypeProjectStop','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'Relative norm','default' => 'T','value' => 'rel'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Maximum error','value' => 'max'},'content' => []}]},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'RelativeLimit','type' => 'real','default' => '0.1','min' => '0','max' => '1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'AbsoluteLimit','type' => 'real','default' => '0.0','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'MaxMatvec','type' => 'integer','default' => '50','min' => '1'},'content' => []},{'type' => 't','content' => '
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


'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'CORRECTP'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'pRatioLow','type' => 'real','default' => '0.01','min' => '0','max' => '1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'pRatioHigh','type' => 'real','default' => '0.1','min' => '$pRatioLow','max' => '1'},'content' => []},{'type' => 't','content' => '
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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'RAYTRACE'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseAccurateIntegral','type' => 'logical','default' => 'T'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseAccurateTrace','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DtExchangeRay','type' => 'real','default' => '0.1','min' => '0.01','max' => '60'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DnRaytrace','type' => 'integer','default' => '1','min' => '1'},'content' => []},{'type' => 't','content' => '

#RAYTRACE
T			UseAccurateIntegral
T			UseAccurateTrace
0.1			DtExchangeRay [sec]
1			DnRaytrace

Raytracing (field-line tracing) is needed to couple the GM with the IM or RB 
components. It can also be used to create plot files with open-closed 
field line information. There are two algorithms implemented for integrating 
rays and for tracing rays.

The UseAccurateIntegral parameter is kept for backwards compatibility
only, it is ignored. The field line integrals are always calculated with 
the accurate algorithm, which follows the lines all the way. 

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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'IM'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TauCoupleIm','type' => 'real','min' => '1'},'content' => []},{'type' => 't','content' => '

#IM
20.0			TauCoupleIm

Same as command IMCOUPLING, except it does not read the two logicals.
See description for command IMCOUPLING.

The default value is TauCoupleIm=20.0, which corresponds to typical nudging.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'IMCOUPLING'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TauCoupleIm','type' => 'real','min' => '1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoCoupleImPressure','type' => 'logical','default' => 'T'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoCoupleImDensity','type' => 'logical','default' => 'F'},'content' => []},{'type' => 't','content' => '

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
'}]}]},{'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'PHYSICS PARAMETERS'},'content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  PHYSICS PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'name' => 'command','type' => 'e','attrib' => {'name' => 'GAMMA','if' => '$_IsFirstSession'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Gamma','type' => 'real','default' => '1.6666666667','min' => '1'},'content' => []},{'type' => 't','content' => '
#GAMMA
1.6666666667		Gamma

The adiabatic index (ratio of the specific heats for fixed pressure
and volume. The default value is 5.0/3.0, which is valid for
monoatomic gas or plasma.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'PLASMA'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'AverageIonMass','type' => 'real','default' => '1','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'AverageIonCharge','type' => 'real','default' => '1','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'ElectronTemperatureRatio','type' => 'real','default' => '0','min' => '0'},'content' => []},{'type' => 't','content' => '
#PLASMA
1.0		AverageIonMass [amu]
1.0		AverageIonCharge [e]
0.0		ElectronTemperatureRatio

The AverageIonMass determines the average mass of ions (and strongly
coupled neutrals) in atomic mass units (amu).
This parameter is only used if there are no species (UseMultiSpecies=.false.
in ModEquation). The number density is n=rho/AverageIonMass.
For a pure hydrogene plasma AverageIonMass=1.0, while
for a mix of 90 per cent hydrogene and 10 per cent helium AverageIonMass=1.4.

The AverageIonCharge determines the average charge of ions (and strongly
coupled neutrals) in electron charge units.
The electron density is ne=n*AverageIonCharge. 
For a fully ionized hydrogene plasma AverageIonCharge=1.0,
for a fully ionized helium plasma AverageIonCharge=2.0, while
for a 10 per cent ionized hydrogene plasma AverageIonCharge=0.1.

The ElectronTemperatureRatio determines the ratio of electron
and ion temperatures. The ion temperature Te = T * ElectronTemperatureRatio
where T is the ion temperature. The total pressure p = n*k*T + ne*k*Te,
so T = p/(n*k+ne*k*ElectronTemperatureRatio). If the electrons and ions are
in temperature equilibrium, ElectronTemperatureRatio=1.0.

In a real plasma all these values can vary in space and time,
but in a single fluid/species MHD description using these constants is
the best one can do. In multispecies MHD the number density can
be determined accurately as n = sum(RhoSpecies_V/(ProtonMass*MassSpecies_V)).

The default values are shown above, which corresponds to
a fully ionized hydrogen plasma where the electron temperature is
negligible relative to the ion temperature. This default is 
backwards compatible with previous versions of the code.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'SHOCKTUBE','if' => '$_IsFirstSession'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'RhoLeft','type' => 'real','default' => '1','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UnLeft','type' => 'real','default' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Ut1Left','type' => 'real','default' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Ut2Left','type' => 'real','default' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'BnLeft','type' => 'real','default' => '0.75'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Bt1Left','type' => 'real','default' => '1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Bt2Left','type' => 'real','default' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'pRight','type' => 'real','default' => '1','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'RhoRight','type' => 'real','default' => '0.125','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UnRight','type' => 'real','default' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Ut1Right','type' => 'real','default' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Ut2Right','type' => 'real','default' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'BnRight','type' => 'real','default' => '0.75'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Bt1Right','type' => 'real','default' => '-1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Bt2Right','type' => 'real','default' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'pRight','type' => 'real','default' => '0.1','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'ShockSlope','type' => 'real','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'no rotation','default' => 'T','value' => '0'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'value' => '0.25'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'value' => '0.3333333333333'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'value' => '0.5'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'value' => '1'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'value' => '2'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'value' => '3'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'value' => '4'},'content' => []}]},{'type' => 't','content' => '
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

If the #SHOCKTUBE command is used in the first session, a 
two state initial condition is set.
The shock is rotated if ShockSlope is not 0, and the tangent of 
the rotation angle is ShockSlope. 
When the shock is rotated, it is best used in combination
with sheared outer boundaries, but then only
\\begin{verbatim}
  ShockSlope = 1., 2., 3., 4., 5.      .....
  ShockSlope = 0.5, 0.33333333, 0.25, 0.2, .....
\\end{verbatim}
can be used, because these angles can be accurately represented
on the grid. The default initial condition is a constant state
with the solar wind parameters.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'SOLARWIND','if' => '$_IsFirstSession'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'SwRhoDim','type' => 'real','default' => '5','min' => '-1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'SwTDim','type' => 'real','default' => '181712.175','min' => '-1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'SwUxDim','type' => 'real','default' => '-400','max' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'SwUyDim','type' => 'real','default' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'SwUzDim','type' => 'real','default' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'SwBxDim','type' => 'real','default' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'SwByDim','type' => 'real','default' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'SwBzDim','type' => 'real','default' => '5'},'content' => []},{'type' => 't','content' => '
#SOLARWIND
5.0			SwRhoDim [n/cc]
181712.175		SwTDim [K]
-400.0			SwUxDim [km/s]
0.0			SwUyDim [km/s]
0.0			SwUzDim [km/s]
0.0			SwBxDim [nT]
0.0			SwByDim [nT]
5.0			SwBzDim [nT]

This command defines the solar wind parameters for the GM component.
The default values are all 0.0-s.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'SOLARWINDFILE','alias' => 'UPSTREAM_INPUT_FILE'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseSolarWindFile','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$UseSolarWindFile'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'NameSolarWindFile','length' => '100','type' => 'string'},'content' => []}]},{'name' => 'rule','type' => 'e','attrib' => {'expr' => '-f $NameSolarWindFile'},'content' => [{'type' => 't','content' => '
		Solar wind file $NameSolarWindFile must exist
	'}]},{'type' => 't','content' => '
#SOLARWINDFILE
T			UseSolarWindFile (rest of parameters read if true)
IMF.dat                 NameSolarWindFile

Default is UseSolarWindFile = .false.

Read IMF data from file NameSolarWindFile if UseSolarWindFile is true.
The data file contains all information required for setting the upstream
boundary conditions. Parameter TypeBcWest should be set to \'vary\' for
the time dependent boundary condition.

If the #SOLARWIND command is not provided then the first time read from
the solar wind file will set the normalization of all variables
in the GM component. Consequently either the #SOLARWIND command or
the #SOLARWINDFILE command with UseSolarWindFile=.true.
is required by the GM component.

The input files are strutured similar to the PARAM.in file.  There are 
{\\tt #commands} that can be  inserted as well as the data.
The file containing the upstream conditions should include data in the 
following order:
\\begin{verbatim}
yr mn dy hr min sec msec bx by bz vx vy vz dens temp
\\end{verbatim}
The units of the variables should be:
\\begin{verbatim}
Magnetic field (b)     nT
Velocity (v)           km/s
Number Density (dens)  cm^-3
Temperature (Temp)     K
\\end{verbatim}

The input files  can have the following optional commands at the beginning
\\begin{verbatim}
#COOR
GSM          The coordinate system of the data (GSM or GSE)

#PLANE       The input data represents values on a tilted plane
30.0         Angle to rotate in the XY plane
30.0         Angle to rotate in the XZ plane

#POSITION    Origin for the plane rotation (see #PLANE)
30.0         Y location
30.0         Z location

#ZEROBX
T            T means Bx is ignored and set to zero

#TIMEDELAY
3600.0       A real number in seconds by which to delay the input
\\end{verbatim}

Finally, the data should be preceded by a {\\tt #START}.  The beginning of 
a typical solar wind input file might look like:
\\begin{verbatim}
#COOR
GSM

#ZEROBX
T

#START
 2004  6  24   0   0  58   0  2.9  -3.1 - 3.7  -300.0  0.0  0.0  5.3  2.00E+04
 2004  6  24   0   1  58   0  3.0  -3.2 - 3.6  -305.0  0.0  0.0  5.4  2.01E+04
\\end{verbatim}

The maximum number of lines of data allowed in the input file is 50,000.  
However, this can be modified by changing the variable Max_Upstream_Npts 
in the file GM/BATSRUS/get_solar_wind_point.f90.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'BODY','if' => '$_IsFirstSession','alias' => 'MAGNETOSPHERE'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseBody','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$UseBody'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'rBody','type' => 'real','default' => '3','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'rCurrents','type' => 'real','default' => '4','min' => '-1','if' => '$_NameComp eq \'GM\''},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'BodyRhoDim','type' => 'real','default' => '1','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'BodyTDim','type' => 'real','default' => '1000','min' => '0'},'content' => []}]},{'type' => 't','content' => '
#BODY
T			UseBody (rest of parameters read if true)
3.0			rBody (user units)
4.0			rCurrents (only read for GM component)
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
! The default values depend on the component: 
! For GM rBody=3, rCurrents=4, BodyRhoDim = 5amu/cc and BodyTDim = 25000 K. 
! For SC/IH: rBody=1, BodyRhoDim = 1.5e8 amu/cc, BodyTDim = 2.8e6 K.
! The rBody and rCurrents are given in planetary/solar radii units.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'GRAVITY','if' => '$_IsFirstSession'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseGravity','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'iDirGravity','type' => 'integer','input' => 'select','if' => '$UseGravity'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'central mass','default' => 'T','value' => '0'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'X direction','value' => '1'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Y direction','value' => '2'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'Z direction','value' => '3'},'content' => []}]},{'type' => 't','content' => '
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

'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'HALLRESISTIVITY'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseHallResist','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$UseHallResist'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'HallFactor','type' => 'real','default' => '1','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'HallCmaxFactor','type' => 'real','default' => '1','min' => '0','max' => '1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'HallHyperFactor','type' => 'real','default' => '0','min' => '0','max' => '1'},'content' => []}]},{'type' => 't','content' => '

#HALLRESISTIVITY
T		UseHallResist (rest of parameters read only if true)
3e12		HallFactor
0.5		HallCmaxFactor
0.2		HallHyperFactor

If UseHallResist is true the Hall resistivity is used.
The off-diagonal Hall elements of the resistivity tensor
are multiplied by HallFactor. If HallFactor is 1 then the 
physical Hall resistivity is used.

If HallCmaxFactor is 1.0 the maximum propagation speed takes into
account the full whistler wave speed. If it is 0, the wave speed
is not midified. For values betwen 1 and 0 a fraction of the whistler
wave speed is added. The full speed is needed for the stability
of the explicit scheme (unless the whistler speed is very small
and/or the diagonal part of the resistivity tensor is dominant).
For the implicit scheme it is better to have HallCmaxFactor=dtexpl/dt,
otherwise there is a significant phase error.

If HallHyperFactor is larger than zero, a 2nd order hyper resistivity
is added: proportional to Delta $x^2 d^4 B/\\Delta x^4$. The 0.1 to 0.2
is a good range to suppress the island formation in the GEM challenge.

Default is UseHallResist false.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'HALLREGION'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'NameHallRegion','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'all','default' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'shell'},'content' => []}]},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$NameHallRegion eq \'shell\''},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'R1Hall','type' => 'real','default' => '0','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'R2Hall','type' => 'real','default' => '0','min' => '$R1Hall'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'HallWidth','type' => 'real','default' => '0','min' => '0'},'content' => []}]},{'type' => 't','content' => '

#HALLREGION
shell			NameHallRegion
 3.5			R1Hall
20.0			R2Hall
 2.0			HallWidth

The NameHallRegion parameter determines the region where the Hall effect
is taken into account. For value "all" the Hall effect is used everywhere
in the computational domain. For "shell" the region is between radii
R1Hall and R2Hall. The edge of the shell is smoothed to avoid a sharp
discontinuity. The width of the smoothing is given by HallWidth.

'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'MASSLOADING'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseMassLoading','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoAccelerateMassLoading','type' => 'logical','default' => 'F'},'content' => []},{'type' => 't','content' => '
#MASSLOADING
F			UseMassLoading
F			DoAccelerateMassLoading
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'HEATFLUX'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseHeatFlux','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseSpitzerForm','type' => 'logical','default' => 'T'},'content' => []},{'name' => 'if','type' => 'e','attrib' => {'expr' => 'not $UseSpitzerForm'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Kappa0Heat','type' => 'real','default' => '1.23E-11'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Kappa0Heat','type' => 'real','default' => '2.5'},'content' => []}]},{'type' => 't','content' => '
#HEATFLUX
T		UseHeatFlux
F		UseSpitzerForm
1.23E-11	Kappa0Heat [W/m/K]	! Only if not UseSpitzerForm
2.50E+00	ExponentHeat [-]	! Only if not UseSpitzerForm
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'RESISTIVITY'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseResistivity','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$UseResistivity'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypeResistivity','type' => 'string','case' => 'lower','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'constant','default' => 'T'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'spitzer'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'anomalous'},'content' => []},{'name' => 'option','type' => 'e','attrib' => {'name' => 'user'},'content' => []}]},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$TypeResistivity =~ /constant/'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Eta0Si','type' => 'real','default' => '1.0E+11','min' => '0.0'},'content' => []}]},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$TypeResistivity =~ /spitzer/'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CoulombLogarithm','type' => 'real','default' => '20'},'content' => []}]},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$TypeResistivity =~ /anomalous/'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Eta0Si','type' => 'real','default' => '1.0E+9','min' => '0.0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Eta0AnomSi','type' => 'real','default' => '2E+09'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'EtaMaxAnomSi','type' => 'real','default' => '2E+10'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'jCritAnomSi','type' => 'real','default' => '1.0E-9','min' => '0'},'content' => []}]}]},{'type' => 't','content' => '
#RESISTIVITY
T		UseResistivity (rest of parameters read only if set to true)
anomalous	TypeResistivity
1.0E+9		Eta0Si       [m2/s] (read for constant and anomalous only)
2.0E+9		Eta0AnomSi   [m2/s] 
2.0E+10		EtaMaxAnomSi [m2/s]
1.0E-9		jCritAnomSi  [A/m2]

If UseResistitivy is false, no resistivity is included.
If UseResistivity is true, then one can select a constant resistivity,
the classical Spitzer resistivity with a CoulombLogarithm parameter 
(typical values are in the range of 10 to 30), 
anomalous resistivity with a critical current, or a user defined resistivity.
For TypeResistivity=\'Spitzer\' the resistivity is very low in space plasma.
For TypeResistivity=\'constant\' the resistivity is uniformly set to Eta0Si.
For TypeResistivity=\'anomalous\' the anomalous resistivity is 
Eta0Si + Eta0AnomSi*(j/jCritAnomSi-1) limited by 0 and EtaMaxAnomSi.
Here j is the absolute value of the current density in SI units.

The default is UseResistivity=.false., i.e. ideal MHD.
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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'SECONDBODY','if' => '$_IsFirstSession'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseBody2','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$UseBody2'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'rBody2','type' => 'real','default' => '0.1','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'xBody2','type' => 'real','default' => '-40','min' => '$xMin','max' => '$xMax'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'yBody2','type' => 'real','default' => '0','min' => '$yMin','max' => '$yMax'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'zBody2','type' => 'real','default' => '0','min' => '$zMin','max' => '$zMax'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'rCurrents2','type' => 'real','default' => '1.3*$rBody2','min' => '$rBody2'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'RhoDimBody2','type' => 'real','default' => '5','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'tDimBody2','type' => 'real','default' => '25000','min' => '0'},'content' => []}]},{'type' => 't','content' => '

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

! For now the dipole of the second body can only be aligned with the z-axis
'}]}]},{'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'SOLAR PROBLEM TYPES'},'content' => [{'type' => 't','content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!! SOLAR PROBLEM TYPES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'},{'name' => 'command','type' => 'e','attrib' => {'name' => 'HELIOSPHERE','if' => '$_IsFirstSession and $_NameComp ne \'GM\''},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'BodyTDim','type' => 'real','default' => '2.85E06','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'BodyRhoDim','type' => 'real','default' => '1.50E8','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'qSun','type' => 'real','default' => '25.0','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'tHeat','type' => 'real','default' => '1.75','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'rHeat','type' => 'real','default' => '1.0','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'SigmaHeat','type' => 'real','default' => '4.5','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoInitRope','type' => 'logical','default' => 'F'},'content' => []},{'name' => 'if','type' => 'e','attrib' => {'expr' => '$DoInitRope'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CmeA','type' => 'real','default' => '0.7','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CmeR1','type' => 'real','default' => '1.2','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CmeR0','type' => 'real','default' => '1.0','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CmeA1','type' => 'real','default' => '0.23'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CmeAlpha','type' => 'real','default' => '0.0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CmeRho1','type' => 'real','default' => '2.5E-12','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CmeRho2','type' => 'real','default' => '2.0E-13','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'ModulationRho','type' => 'real','default' => '0.0','min' => '0','max' => '10'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'ModulationP','type' => 'real','default' => '0.0','min' => '0','max' => '10'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'OrientationGL98','type' => 'real','default' => '0.0','min' => '-360','max' => '360'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'LatitudeGL98','type' => 'real','default' => '0.0','min' => '-90','max' => '90'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'LongitudeGL98','type' => 'real','default' => '0.0','min' => '-360','max' => '360'},'content' => []}]},{'type' => 't','content' => '

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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'HELIOUPDATEB0','if' => '$_NameComp ne \'GM\''},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DtUpdateB0','type' => 'real','default' => '0.0001','min' => '-1'},'content' => []},{'type' => 't','content' => '

#HELIOUPDATEB0
-1.0			DtUpdateB0 [s]

Set the frequency of updating the B0 field for the solar corona.
A negative value means that the B0 field is not updated.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'HELIODIPOLE','if' => '$_NameComp ne \'GM\''},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'HelioDipoleStrength','type' => 'real'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'HelioDipoleTilt','type' => 'real','default' => '0','min' => '-90','max' => '90'},'content' => []},{'type' => 't','content' => '

#HELIODIPOLE
-3.0                    HelioDipoleStrength [G]
 0.0                    HelioDipoleTilt     [deg]

! Variable HelioDipoleStrength defines the equatorial field strength in Gauss,
! while HelioDipoleTilt is the tilt relative to the ecliptic North 
! (negative sign means towards the planet) in degrees.
!
! Default value is HelioDipoleStrength = 0.0.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'HELIOTEST','if' => '$_NameComp ne \'GM\''},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoSendMHD','type' => 'logical','default' => 'F'},'content' => []},{'type' => 't','content' => '

#HELIOTEST
F			DoSendMHD

! If DoSendMHD is true, IH sends the real MHD solution to GM in the coupling.
! If DoSendMHD is false then the values read from the IMF file are sent,
! so there is no real coupling. Mostly used for testing the framework.
!
! Default value is true, i.e. real coupling.
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'HELIOBUFFERGRID','if' => '$_IsFirstSession and $_NameComp eq \'IH\''},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'rBuffMin','type' => 'real','default' => '19','min' => '1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'rBuffMax','type' => 'real','default' => '21','min' => '1'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nThetaBuff','type' => 'integer','default' => '45','min' => '18'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nPhiBuff','type' => 'integer','default' => '90','min' => '36'},'content' => []},{'type' => 't','content' => '
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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'CME','if' => '$_IsFirstSession and $_NameComp ne \'GM\''},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypeCme','type' => 'string','input' => 'select'},'content' => [{'name' => 'option','type' => 'e','attrib' => {'name' => 'Low','default' => 'T'},'content' => []}]},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CmeA','type' => 'real','default' => '0.7','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CmeR1','type' => 'real','default' => '1.2','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CmeR0','type' => 'real','default' => '1.0','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CmeA1','type' => 'real','default' => '0.23'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CmeAlpha','type' => 'real','default' => '0.0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CmeRho1','type' => 'real','default' => '2.5E-12','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CmeRho2','type' => 'real','default' => '2.0E-13','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CmeB1Dim','type' => 'real','default' => '1.0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CmeUErupt','type' => 'real','default' => '4.0E5','min' => '0'},'content' => []},{'type' => 't','content' => '
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
'}]},{'name' => 'command','type' => 'e','attrib' => {'name' => 'ARCADE','if' => '$_IsFirstSession and $_NameComp ne \'GM\''},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'tArcDim','type' => 'real','default' => '1.0E6','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'RhoArcDim','type' => 'real','default' => '1.0E-12','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'bArcDim','type' => 'real','default' => '0.718144','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'ByArcDim','type' => 'real','default' => '1.0E6'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UzArcDim','type' => 'real','default' => '5.0E3'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Phi0Arc','type' => 'real','default' => '0.5'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'MuArc','type' => 'real','default' => '1.3'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'ExpArc','type' => 'real','default' => '3','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'WidthArc','type' => 'real','default' => '0.5','min' => '0'},'content' => []},{'type' => 't','content' => '
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
'},{'name' => 'command','type' => 'e','attrib' => {'name' => 'COMET','if' => '$_IsFirstSession'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'ProdRate','type' => 'real','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UrNeutral','type' => 'real','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'AverageMass','type' => 'real','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'IonizationRate','type' => 'real','min' => '0'},'content' => []},{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'kFriction','type' => 'real','min' => '0'},'content' => []},{'type' => 't','content' => '
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
'},{'name' => 'command','type' => 'e','attrib' => {'name' => 'INCLUDE'},'content' => [{'name' => 'parameter','type' => 'e','attrib' => {'name' => 'NameIncludeFile','length' => '100','type' => 'string','default' => 'Param/'},'content' => []},{'type' => 't','content' => '

#INCLUDE
Param/SSS_3000		NameIncludeFile

! Include a library file from Param/ or any file from anywhere else.
'}]}]},{'name' => 'rule','type' => 'e','attrib' => {'expr' => '($SwRhoDim > 0) or $UseSolarWindFile or $_NameComp ne \'GM\''},'content' => [{'type' => 't','content' => '
	Either command #SOLARWIND or #SOLARWINDFILE must be used!
'}]},{'name' => 'rule','type' => 'e','attrib' => {'expr' => '$MaxImplBlock>1 or not $UsePartImplicit or not $MaxImplBlock'},'content' => [{'type' => 't','content' => '
	Part implicit scheme requires more than 1 implicit block!
'}]},{'name' => 'rule','type' => 'e','attrib' => {'expr' => '$MaxImplBlock==$MaxBlock or not $UseFullImplicit'},'content' => [{'type' => 't','content' => '
	Full implicit scheme should be used with equal number of 
	explicit and implicit blocks!
'}]},{'name' => 'rule','type' => 'e','attrib' => {'expr' => '-d $NameRestartOutDir or not $_IsFirstSession'},'content' => [{'type' => 't','content' => '
	Output restart directory $NameRestartOutDir should exist!
'}]},{'name' => 'rule','type' => 'e','attrib' => {'expr' => '-d $NamePlotDir or not $_IsFirstSession'},'content' => [{'type' => 't','content' => '
	Plot directory $NamePlotDir should exist!
'}]}]}];