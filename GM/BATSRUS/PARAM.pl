#^CFG FILE _FALSE_
$tree = [{'content' => [{'content' => '

List of MH (GM, IH and SC) commands used in the PARAM.in file




','type' => 't'},{'content' => [],'name' => 'set','type' => 'e','attrib' => {'name' => 'nI','value' => '$_GridSize[0]','type' => 'integer'}},{'content' => [],'name' => 'set','type' => 'e','attrib' => {'name' => 'nJ','value' => '$_GridSize[1]','type' => 'integer'}},{'content' => [],'name' => 'set','type' => 'e','attrib' => {'name' => 'nK','value' => '$_GridSize[2]','type' => 'integer'}},{'content' => [],'name' => 'set','type' => 'e','attrib' => {'name' => 'MaxBlock','value' => '$_GridSize[3]','type' => 'integer'}},{'content' => [],'name' => 'set','type' => 'e','attrib' => {'name' => 'MaxImplBlock','value' => '$_GridSize[4]','type' => 'integer'}},{'content' => [],'name' => 'set','type' => 'e','attrib' => {'name' => 'MaxBlockALL','value' => '$_nProc and $MaxBlock and $_nProc*$MaxBlock','type' => 'integer'}},{'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! STAND ALONE PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

','type' => 't'},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseNewParam','type' => 'logical','default' => 'T'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseNewAxes','type' => 'logical','default' => 'T'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoTimeAccurate','type' => 'logical','default' => 'T'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseCorotation','type' => 'logical','default' => 'T'}},{'content' => '

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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_IsStandAlone','name' => 'NEWPARAM'}},{'content' => [{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'if' => '$_NameComp eq \'SC\'','name' => 'SC','default' => 'T'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'if' => '$_NameComp eq \'IH\'','name' => 'IH','default' => 'T'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'if' => '$_NameComp eq \'GM\'','name' => 'GM','default' => 'T'}}],'name' => 'parameter','type' => 'e','attrib' => {'input' => 'options','name' => 'NameComp','type' => 'string'}},{'content' => '

#COMPONENT
GM			NameComp

This command is only used in the stand alone mode.

The NameComp variable contains the two-character component ID
for the component which BATSRUS is representing.
If NameComp does not agree with the value of the NameThisComp
variable, BATSRUS stops with an error message.
This command is saved into the restart header file for consistency check.

There is no default value: if the command is not given, the component ID is not checked.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'COMPONENT'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'length' => '100','name' => 'StringDescription','type' => 'string'}},{'content' => '

#DESCRIPTION
This is a test run for Jupiter with no rotation.

This command is only used in the stand alone mode.

The StringDescription string can be used to describe the simulation
for which the parameter file is written. The #DESCRIPTION command and
the StringDescription string are saved into the restart file,
which helps in identifying the restart files.

The default value is "Please describe me!", which is self explanatory.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_IsStandAlone','name' => 'DESCRIPTION'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoEcho','type' => 'logical','default' => 'F'}},{'content' => '

#ECHO
T                       DoEcho

This command is only used in the stand alone mode.

If the DoEcho variable is true, the input parameters are echoed back.
The default value for DoEcho is .false., but it is a good idea to
set it to true at the beginning of the PARAM.in file.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_IsStandAlone','name' => 'ECHO'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DnProgressShort','type' => 'integer','default' => '10','min' => '-1'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DnProgressLong','type' => 'integer','default' => '100','min' => '-1'}},{'content' => '
#PROGRESS
10			DnProgressShort
100			DnProgressLong

The frequency of short and long progress reports for BATSRUS in
stand alone mode. These are the defaults. Set -1-s for no progress reports.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_IsStandAlone','name' => 'PROGRESS'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoTimeAccurate','type' => 'logical','default' => 'T'}},{'content' => '

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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_IsStandAlone','name' => 'TIMEACCURATE'}},{'content' => [{'content' => '

This command is allowed in stand alone mode only for sake of the 
test suite, which contains these commands when the framework is tested.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_IsStandAlone','multiple' => 'T','name' => 'BEGIN_COMP'}},{'content' => [{'content' => '

This command is allowed in stand alone mode only for sake of the 
test suite, which contains these commands when the framework is tested.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_IsStandAlone','multiple' => 'T','name' => 'END_COMP'}},{'content' => [{'content' => '

#RUN

This command is only used in stand alone mode.

The #RUN command does not have any parameters. It signals the end
of the current session, and makes BATSRUS execute the session with
the current set of parameters. The parameters for the next session
start after the #RUN command. For the last session there is no
need to use the #RUN command, since the #END command or simply
the end of the PARAM.in file makes BATSRUS execute the last session.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_IsStandAlone','name' => 'RUN'}},{'content' => [{'content' => '

#END

The #END command signals the end of the included file or the
end of the PARAM.in file. Lines following the #END command are
ignored. It is not required to use the #END command. The end
of the included file or PARAM.in file is equivalent with an 
#END command in the last line.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'END'}}],'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'STAND ALONE MODE'}},{'content' => [{'content' => '
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

','type' => 't'},{'content' => [{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Earth','value' => 'EARTH/Earth/earth','default' => 'T'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Saturn','value' => 'SATURN/Saturn/saturn'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'New'}}],'name' => 'parameter','type' => 'e','attrib' => {'input' => 'select','name' => 'NamePlanet','type' => 'string'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'RadiusPlanet','type' => 'real','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'MassPlanet','type' => 'real','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'OmegaPlanet','type' => 'real','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TiltRotation','type' => 'real','min' => '0'}},{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'NONE'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'DIPOLE','default' => 'T'}}],'name' => 'parameter','type' => 'e','attrib' => {'input' => 'select','name' => 'TypeBField','type' => 'string'}}],'name' => 'if','type' => 'e','attrib' => {'expr' => '$NamePlanet eq \'New\''}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '180','name' => 'MagAxisThetaGeo','type' => 'real','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '360','name' => 'MagAxisPhiGeo','type' => 'real','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DipoleStrength','type' => 'real'}}],'name' => 'if','type' => 'e','attrib' => {'expr' => '$TyepBField eq \'DIPOLE\''}},{'content' => [{'content' => '
		PLANET should precede $PlanetCommand
	','type' => 't'}],'name' => 'rule','type' => 'e','attrib' => {'expr' => 'not $PlanetCommand'}},{'content' => '

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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_IsFirstSession and $_IsStandAlone','name' => 'PLANET'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'IsRotAxisPrimary','type' => 'logical','default' => 'T'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '180','name' => 'RotAxisTheta','type' => 'real','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '360','name' => 'RotAxisPhi','type' => 'real','min' => '0'}}],'name' => 'if','type' => 'e','attrib' => {'expr' => '$IsRotAxisPrimary'}},{'content' => [],'name' => 'set','type' => 'e','attrib' => {'name' => 'PlanetCommand','value' => 'ROTATIONAXIS','type' => 'string'}},{'content' => '

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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_IsStandAlone','name' => 'ROTATIONAXIS'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseRotation','type' => 'logical','default' => 'T'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'RotationPeriod','type' => 'real'}}],'name' => 'if','type' => 'e','attrib' => {'expr' => '$UseRotation'}},{'content' => [],'name' => 'set','type' => 'e','attrib' => {'name' => 'PlanetCommand','value' => 'MAGNETICAXIS','type' => 'string'}},{'content' => '

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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_IsStandAlone','name' => 'ROTATION'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'IsMagAxisPrimary','type' => 'logical','default' => 'T'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '180','name' => 'MagAxisTheta','type' => 'real','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '360','name' => 'MagAxisPhi','type' => 'real','min' => '0'}}],'name' => 'if','type' => 'e','attrib' => {'expr' => '$IsMagAxisPrimary'}},{'content' => [],'name' => 'set','type' => 'e','attrib' => {'name' => 'PlanetCommand','value' => 'MAGNETICAXIS','type' => 'string'}},{'content' => '

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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_IsStandAlone','name' => 'MAGNETICAXIS'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DipoleStrength','type' => 'real'}},{'content' => '

#DIPOLE
-3.11e-4		DipoleStrength [Tesla]

The DipoleStrength variable contains the
magnetic equatorial strength of the dipole magnetic field in Tesla.

The default value is the real dipole strength for the planet.
For the Earth the default is taken to be -31100 nT.
The sign is taken to be negative so that the magnetic axis can
point northward as usual.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_IsStandAlone','name' => 'DIPOLE'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DtUpdateB0','type' => 'real','default' => '0.0001','min' => '-1'}},{'content' => '

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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_IsStandAlone','name' => 'UPDATEB0'}},{'content' => [{'content' => '

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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_IsStandAlone','name' => 'IDEALAXES'}}],'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'PLANET COMMANDS'}},{'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!  USER DEFINED INPUT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

','type' => 't'},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseUserInnerBcs','type' => 'logical','default' => 'F'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseUserSource','type' => 'logical','default' => 'F'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseUserPerturbation','type' => 'logical','default' => 'F'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseUserOuterBcs','type' => 'logical','default' => 'F'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseUserICs','type' => 'logical','default' => 'F'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseUserSpecifyRefinement','type' => 'logical','default' => 'F'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseUserLogFiles','type' => 'logical','default' => 'F'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseUserWritePlot','type' => 'logical','default' => 'F'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseUserAMR','type' => 'logical','default' => 'F'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseUserEchoInput','type' => 'logical','default' => 'F'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseUserB0','type' => 'logical','default' => 'F'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseUserSetPhysConst','type' => 'logical','default' => 'F'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseUserUpdateStates','type' => 'logical','default' => 'F'}},{'content' => '

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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'USER_FLAGS'}},{'content' => [{'content' => '

This command signals the beginning of the section of the file which 
is read by the subroutine user\\_read\\_inputs in the user\\_routines.f90 file.
The section ends with the #USERINPUTEND command. There is no XML based parameter
checking in the user section.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'USERINPUTBEGIN'}},{'content' => [{'content' => '

This command signals the end of the section of the file which 
is read by the subroutine user\\_read\\_inputs in the user\\_routines.f90 file.
The section begins with the #USERINPUTBEGIN command. There is no XML based parameter
checking in the user section.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'USERINPUTEND'}}],'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'USER DEFINED INPUT'}},{'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  TESTING AND TIMING PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
','type' => 't'},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'length' => '100','name' => 'TestString','type' => 'string'}},{'content' => '
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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'TEST'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '$nI+2','name' => 'iTest','type' => 'integer','min' => '-2'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '$nJ+2','name' => 'jTest','type' => 'integer','min' => '-2'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '$nK+2','name' => 'kTest','type' => 'integer','min' => '-2'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '$MaxBlock','name' => 'iBlockTest','type' => 'integer','min' => '1'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'iProcTest','type' => 'integer','min' => '0'}},{'content' => '
#TESTIJK
1                       iTest           (cell index for testing)
1                       jTest           (cell index for testing)
1                       kTest           (cell index for testing)
1                       iBlockTest      (block index for testing)
0                       iProcTest       (processor index for testing)

! The location of test info in terms of indices, block and processor number.
! Note that the user should set #TESTIJK or #TESTXYZ, not both.  If both
! are set, the final one in the session will set the test point.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'TESTIJK'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '$xMax','name' => 'xTest','type' => 'real','min' => '$xMin'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '$yMax','name' => 'yTest','type' => 'real','min' => '$yMin'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '$zMax','name' => 'zTest','type' => 'real','min' => '$zMin'}},{'content' => '
#TESTXYZ
1.5                     xTest           (X coordinate of cell for testing)
-10.5                   yTest           (Y coordinate of cell for testing)
-10.                    zTest           (Z coordinate of cell for testing)

! The location of test info in terms of coordinates.
! Note that the user should set #TESTIJK or #TESTXYZ, not both.  If both
! are set, the final one in the session will set the test point.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'TESTXYZ'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nIterTest','type' => 'integer','default' => '-1','min' => '-1'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TimeTest','type' => 'real','default' => '1e30','min' => '-1'}},{'content' => '

#TESTTIME
-1                      nIterTest       (iteration number to start testing)
10.5                    TimeTest        (time to start testing in seconds)

! The time step and physical time to start testing.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'TESTTIME'}},{'content' => [{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Rho','value' => '1','default' => 'T'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'RhoUx','value' => '2'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'RhoUy','value' => '3'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'RhoUz','value' => '4'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Bx','value' => '5'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'By','value' => '6'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Bz','value' => '7'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'e','value' => '8'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'p','value' => '9'}}],'name' => 'parameter','type' => 'e','attrib' => {'input' => 'select','name' => 'iVarTest','type' => 'integer'}},{'content' => '
#TESTVAR
1                       iVarTest

! Index of variable to be tested. Default is rho_="1", ie. density.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'TESTVAR'}},{'content' => [{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'all','value' => '0'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'x','value' => '1','default' => 'T'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'y','value' => '2'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'z','value' => '3'}}],'name' => 'parameter','type' => 'e','attrib' => {'input' => 'select','name' => 'iVarTest','type' => 'integer'}},{'content' => '
#TESTDIM
1                       iDimTest

! Index of dimension/direction to be tested. Default is X dimension.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'TESTDIM'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseStrict','type' => 'logical','default' => 'T'}},{'content' => '
#STRICT
T                       UseStrict

! If true then stop when parameters are incompatible. If false, try to
! correct parameters and continue. Default is true, ie. strict mode
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'STRICT'}},{'content' => [{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'errors and warnings only','value' => '-1'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'start and end of sessions','value' => '0'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'normal','value' => '1','default' => 'T'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'calls on test processor','value' => '10'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'calls on all processors','value' => '100'}}],'name' => 'parameter','type' => 'e','attrib' => {'input' => 'select','name' => 'lVerbose','type' => 'integer'}},{'content' => '
#VERBOSE
-1                      lVerbose

! Verbosity level controls the amount of output to STDOUT. Default level is 1.
!   lVerbose .le. -1 only warnings and error messages are shown.
!   lVerbose .ge.  0 start and end of sessions is shown.
!   lVerbose .ge.  1 a lot of extra information is given.
!   lVerbose .ge. 10 all calls of set_oktest are shown for the test processor.
!   lVerbose .ge.100 all calls of set_oktest are shown for all processors.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'VERBOSE'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoDebug','type' => 'logical','default' => 'F'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoDebugGhost','type' => 'logical','default' => 'F'}},{'content' => '
#DEBUG
F                       DoDebug         (use it as if(okdebug.and.oktest)...)
F                       DoDebugGhost    (parameter for show_BLK in library.f90)

! Excessive debug output can be controlled by the global okdebug parameter
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'DEBUG'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CodeVersion','type' => 'real','default' => '7.50','min' => '0'}},{'content' => '
#CODEVERSION
7.50                    CodeVersion

! Checks CodeVersion. Prints a WARNING if it differs from the CodeVersion
! defined in ModMain.f90. Used in newer restart header files. 
! Should be given in PARAM.in when reading old restart files, 
! which do not have version info in the header file.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_IsFirstSession','name' => 'CODEVERSION'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'length' => '100','name' => 'NameEquation','type' => 'string','default' => 'MHD'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nVar','type' => 'integer','default' => '8'}},{'content' => '
#EQUATION
MHD			NameEquation
8			nVar

! Define the equation name and the number of variables.
! If any of these do not agree with the values determined 
! by the code, BATSRUS stops with an error. Used in restart
! header files and can be given in PARAM.in as a check
! and as a description.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_IsFirstSession','name' => 'EQUATION'}},{'content' => [{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'single precision (4)','value' => '4','default' => '$_nByteReal==4'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'double precision (8)','value' => '8','default' => '$_nByteReal==8'}}],'name' => 'parameter','type' => 'e','attrib' => {'input' => 'select','name' => 'nByteReal','type' => 'integer'}},{'content' => [{'content' => '
		nByteReal in file must agree with _nByteReal.
	','type' => 't'}],'name' => 'rule','type' => 'e','attrib' => {'expr' => '$nByteReal==$_nByteReal'}},{'content' => '

#PRECISION
8                       nByteReal

! Define the number of bytes in a real number. If it does not agree
! with the value determined by the code, BATSRUS stops with an error.
! This is a check, the internal value is calculated in parallel_setup.
! Used in latest restart header files to check binary compatibility.
! May be given in PARAM.in to enforce a certain precision.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_IsFirstSession','name' => 'PRECISION'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '$nI','name' => 'nI','type' => 'integer','default' => '$nI','min' => '$nI'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '$nJ','name' => 'nJ','type' => 'integer','default' => '$nJ','min' => '$nJ'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '$nK','name' => 'nK','type' => 'integer','default' => '$nK','min' => '$nK'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '$MaxBlockALL','name' => 'MinBlockALL','type' => 'integer','min' => '1'}},{'content' => '

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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_IsFirstSession','name' => 'CHECKGRIDSIZE'}},{'content' => [{'content' => '
#BLOCKLEVELSRELOADED

This command means that the restart file contains the information about
the minimum and maximum allowed refinement levels for each block.
This command is only used in the restart header file.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'BLOCKLEVELSRELOADED'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseTiming','type' => 'logical','default' => 'T'}},{'content' => [{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'none','value' => '-3'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'final only','value' => '-2','default' => 'T'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'end of sessions','value' => '-1'}},{'content' => [],'name' => 'optioninput','type' => 'e','attrib' => {'name' => 'every X steps','default' => '100','min' => '1'}}],'name' => 'parameter','type' => 'e','attrib' => {'input' => 'select','name' => 'Frequency','type' => 'integer'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nDepthTiming','type' => 'integer','default' => '-1','min' => '-1'}},{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'cummulative','value' => 'cumm','default' => '1'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'list'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'tree'}}],'name' => 'parameter','type' => 'e','attrib' => {'input' => 'select','name' => 'TypeTimingReport','type' => 'string'}}],'name' => 'if','type' => 'e','attrib' => {'expr' => '$UseTiming'}},{'content' => '
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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'TIMING'}}],'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'TESTING AND TIMING'}},{'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! MAIN INITIAL AND BOUNDARY CONDITION PARAMETERS  !!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
','type' => 't'},{'content' => [{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Uniform','value' => '1'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Shock tube','value' => '2'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Heliosphere','value' => '3'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Comet','value' => '5'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Rotation','value' => '6'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Diffusion','value' => '7'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Earth','value' => '11','default' => 'T'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Saturn','value' => '12'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Jupiter','value' => '13'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Venus','value' => '14'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Cylinder','value' => '21'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Sphere','value' => '22'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Arcade','value' => '25'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'CME','value' => '26'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Dissipation','value' => '30'}}],'name' => 'parameter','type' => 'e','attrib' => {'input' => 'select','name' => 'iProblem','type' => 'integer'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'if' => '$iProblem==30','length' => '20','name' => 'TypeDissipation','type' => 'string'}},{'content' => '
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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_IsFirstSession','name' => 'PROBLEMTYPE','required' => 'T'}},{'content' => [{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'if' => '$_NameComp eq \'GM\'','name' => 'GeoSolarMagnetic, GSM','value' => 'GSM','default' => 'T'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'if' => '$_NameComp ne \'GM\'','name' => 'HelioGraphicInertial, HGI','value' => 'HGI','default' => 'T'}}],'name' => 'parameter','type' => 'e','attrib' => {'input' => 'select','name' => 'TypeCoordSystem','type' => 'string'}},{'content' => '

#COORDSYSTEM
GSM			TypeCoordSystem

! TypeCoordSystem defines the coordinate system for the component.
! Currently only one coordinate system is available for GM ("GSM")
! and one for IH or SC ("HGI"). In the near future "GSE" should be also
! an option for GM.
!
! Default is component dependent: "GSM" for GM and "HGI" for IH or SC.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_IsFirstSession','name' => 'COORDSYSTEM'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'length' => '100','name' => 'NameRestartInDir','type' => 'string','default' => 'GM/restartIN'}},{'content' => '

#RESTARTINDIR
GM/restart_n5000	NameRestartInDir

! The NameRestartInDir variable contains the name of the directory
! where restart files are saved relative to the run directory.
! The directory should be inside the subdirectory with the name 
! of the component.
!
! Default value is "GM/restartIN".
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_IsFirstSession','name' => 'RESTARTINDIR'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoRestartBFace','type' => 'logical','default' => 'F'}},{'content' => ' 

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

','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_IsFirstSession','name' => 'NEWRESTART'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nRootBlockX','type' => 'integer','default' => '2','min' => '1'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nRootBlockY','type' => 'integer','default' => '1','min' => '1'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nRootBlockZ','type' => 'integer','default' => '1','min' => '1'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'xMin','type' => 'real','default' => '-192.0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'xMax','type' => 'real','default' => '  64.0','min' => '$xMin'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'yMin','type' => 'real','default' => ' -64.0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'yMax','type' => 'real','default' => '  64.0','min' => '$yMin'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'zMin','type' => 'real','default' => ' -64.0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'zMax','type' => 'real','default' => '  64.0','min' => '$zMin'}},{'content' => '
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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_IsFirstSession','name' => 'GRID','required' => 'T'}},{'content' => [{'content' => [{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'coupled'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'fixed/inflow','default' => '$Side ne \'TypeBcEast\''}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'float/outflow','default' => '$Side eq \'TypeBcEast\''}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'heliofloat'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'reflect'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'periodic'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'vary'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'shear'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'linetied'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'raeder'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'arcadetop'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'arcadebot'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'arcadebotcont'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'user'}}],'name' => 'parameter','type' => 'e','attrib' => {'input' => 'select','name' => '$Side','type' => 'string'}}],'name' => 'foreach','type' => 'e','attrib' => {'values' => 'TypeBcEast,TypeBcWest,TypeBcSouth,TypeBcNorth,TypeBcBot,TypeBcTop','name' => 'Side'}},{'content' => [{'content' => '
	! East and west BCs must be both periodic or neither
	','type' => 't'}],'name' => 'rule','type' => 'e','attrib' => {'expr' => 'not($TypeBcEast eq \'periodic\' xor $TypeBcWest eq \'periodic\')'}},{'content' => [{'content' => '
	! South and North BCs must be both periodic or neither
	','type' => 't'}],'name' => 'rule','type' => 'e','attrib' => {'expr' => 'not($TypeBcSouth eq \'periodic\' xor $TypeBcNorth eq \'periodic\')'}},{'content' => [{'content' => '
	! Bottom and top BCs must be both periodic or neither
	','type' => 't'}],'name' => 'rule','type' => 'e','attrib' => {'expr' => 'not($TypeBcBot eq \'periodic\' xor $TypeBcTop eq \'periodic\')'}},{'content' => '
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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'OUTERBOUNDARY'}},{'content' => [{'content' => '
! Inner boundary types for body 1 and body 2
	','type' => 't'},{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'reflect'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'float'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'fixed'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'ionosphere','default' => 'T'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'ionosphereB0','value' => 'ionosphereB0/ionosphereb0'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'ionospherefloat'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'coronatoih'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'user'}}],'name' => 'parameter','type' => 'e','attrib' => {'input' => 'select','name' => 'TypeBcInner','type' => 'string'}},{'content' => [{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'reflect','default' => 'T'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'float'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'fixed'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'ionosphere'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'ionosphereB0'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'ionospherefloat'}}],'name' => 'parameter','type' => 'e','attrib' => {'input' => 'select','name' => 'TypeBcBody2','type' => 'string'}}],'name' => 'if','type' => 'e','attrib' => {'expr' => '$UseBody2'}},{'content' => [{'content' => '
! Note: for the second body COROTATION AND AN IONOSPHERIC BOUNDARY DO NOT WORK.
	','type' => 't'}],'name' => 'rule','type' => 'e','attrib' => {'expr' => 'not($TypeBcBody2 =~ /ionosphere/)'}},{'content' => '

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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'INNERBOUNDARY'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseExtraBoundary','type' => 'logical','default' => 'F'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TypeBcExtra','type' => 'string'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoFixExtraboundary','type' => 'logical','default' => 'F'}}],'name' => 'if','type' => 'e','attrib' => {'expr' => '$UseExtraBoundary'}}],'name' => 'command','type' => 'e','attrib' => {'name' => 'EXTRABOUNDARY'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '6','name' => 'MaxBoundary','type' => 'integer','default' => '0','min' => '0'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoFixOuterBoundary','type' => 'logical','default' => 'F'}}],'name' => 'if','type' => 'e','attrib' => {'expr' => '$MaxBoundary >= 1'}},{'content' => '
#FACEOUTERBC
0              MaxBoundary            
F              DoFixOuterBoundary)    !read only for MaxBoundary>=East_(=1).
! If MaxBoundary is East_(=1) or more then the outer boundaries with
! the number of boundary being between East_ and MaxBoundary
! are treated using set_BCs.f90 subroutines instead of set_outerBCs.f90 
! if DoFixOuterBoundary is .true., there is no resolution
! change along the outer boundaries with the number of
! of boundary being between East_ and MaxBoundary
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'FACEOUTERBC'}}],'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'INITIAL AND BOUNDARY CONDITIONS'}},{'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! INITIAL TIME AND STEP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
','type' => 't'},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'iYear','type' => 'integer','default' => '2000'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '12','name' => 'iMonth','type' => 'integer','default' => '3','min' => '1'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '31','name' => 'iDay','type' => 'integer','default' => '21','min' => '1'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '23','name' => 'iHour','type' => 'integer','default' => '0','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '59','name' => 'iMinute','type' => 'integer','default' => '0','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '59','name' => 'iSecond','type' => 'integer','default' => '0','min' => '0'}},{'content' => '
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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_IsFirstSession','name' => 'STARTTIME','alias' => 'SETREALTIME'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'tSimulation','type' => 'real','default' => '0.0','min' => '0'}},{'content' => '

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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_IsFirstSession','name' => 'TIMESIMULATION'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nStep','type' => 'integer','default' => '0','min' => '0'}},{'content' => '

#NSTEP
100			nStep

! Set nStep for the component. Typically used in the restart.H header file.
! Generally it is not inserted in a PARAM.in file by the user.
!
! The default is nStep=0 as the starting time step with no restart.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_IsFirstSession','name' => 'NSTEP'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nPrevious','type' => 'integer','default' => '-1','min' => '-1'}},{'content' => '

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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_IsFirstSession','name' => 'NPREVIOUS'}}],'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'INITIAL TIME'}},{'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  TIME INTEGRATION PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
','type' => 't'},{'content' => [{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'value' => '1','default' => 'T'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'value' => '2'}}],'name' => 'parameter','type' => 'e','attrib' => {'input' => 'select','name' => 'nStage','type' => 'integer'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '1','name' => 'CflExpl','type' => 'real','default' => '0.8','min' => '0'}},{'content' => '

#TIMESTEPPING
2                       nStage
0.80                    CflExpl

! Parameters for explicit time integration.
! Default is 1 stage and CflExpl=0.8
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'TIMESTEPPING'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseDtFixed','type' => 'logical','default' => 'F'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'if' => '$UseDtFixed','name' => 'DtFixedDim','type' => 'real','default' => '1.0','min' => '0'}},{'content' => '
#FIXEDTIMESTEP
T                       UseDtFixed
10.                     DtFixedDim [sec] (read if UseDtFixed is true)

! Default is UseDtFixed=.false. Effective only if DoTimeAccurate is true.
! If UseDtFixed is true, the time step is fixed to DtFixedDim.
!
! This is useful for debugging explicit schemes.

! The real application is, however, for implicit and partially
! implicit/local schemes.

','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'FIXEDTIMESTEP'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UsePartLocal','type' => 'logical','default' => 'F'}},{'content' => '
#PARTLOCAL
T               UsePartLocal

! Default is UsePartLocal=.false. If UsePartLocal is true and the
! run is time accurate, then the blocks selected as "implicit"
! by the criteria defined in #STEPPINGCRITERIA are not used to
! calculate the time step, and all cells are advanced with the
! smaller of the stable and the global time steps.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'PARTLOCAL'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UsePointImplicit','type' => 'logical','default' => 'F'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UsePartImplicit','type' => 'logical','default' => 'F'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseFullImplicit','type' => 'logical','default' => 'F'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'if' => '$UsePartImplicit or $UseFullImplicit','name' => 'CflImpl','type' => 'real','default' => '100','min' => '0'}},{'content' => [{'content' => '
	At most one of these logicals can be true!
	','type' => 't'}],'name' => 'rule','type' => 'e','attrib' => {'expr' => '$UsePointImplicit + $UsePartImplicit + $UseFullImplicit <= 1'}},{'content' => '

#IMPLICIT
F               UsePointImplicit   
F               UsePartImplicit
F               UseFullImplicit
100.0           CflImpl (read if UsePartImplicit or UseFullImplicit is true)

! Default is false for all logicals. Only one of them can be set to true!
! The CFL number is used in the implicit blocks of the fully or partially
! implicit schemes. Ignored if UseDtFixed is true.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'IMPLICIT'}}],'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'TIME INTEGRATION'}},{'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! PARAMETERS FOR FULL AND PART IMPLICIT TIME INTEGRATION !!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
','type' => 't'},{'content' => [{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Time step','value' => 'dt','default' => 'T'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Radial distance','value' => 'r/R'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Test block','value' => 'test'}}],'name' => 'parameter','type' => 'e','attrib' => {'input' => 'select','name' => 'TypeImplCrit','type' => 'string'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'if' => '$TypeImplCrit eq \'R\'','name' => 'rImplicit','type' => 'real','min' => '0'}},{'content' => '

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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'IMPLICITCRITERIA','alias' => 'STEPPINGCRITERIA'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '1','name' => 'ImplCoeff','type' => 'real','default' => '1','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseBdf2','type' => 'logical','default' => 'T'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseSourceImpl','type' => 'logical','default' => 'F'}},{'content' => '
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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'IMPLSTEP'}},{'content' => [{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => '1','default' => 'T'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => '2'}}],'name' => 'parameter','type' => 'e','attrib' => {'input' => 'select','name' => 'nOrderImpl','type' => 'integer'}},{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Roe','value' => 'Roe/roe/1'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Rusanov','value' => 'Rusanov/rusanov/2/TVDLF','default' => 'T'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Linde','value' => 'Linde/linde/3/HLLEL'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Sokolov','value' => 'Sokolov/sokolov/4/AW'}}],'name' => 'parameter','type' => 'e','attrib' => {'input' => 'select','name' => 'TypeFluxImpl','type' => 'string'}},{'content' => '
#IMPLSCHEME
1               nOrderImpl
Rusanov         TypeFluxImpl

! Default values are shown, ie. first order Rusanov scheme.
! This defines the scheme used in the implicit part.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'IMPLSCHEME'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseConservativeImplicit','type' => 'logical','default' => 'F'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseNewton','type' => 'logical','default' => 'F'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseNewMatrix','type' => 'logical','default' => 'T'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'MaxIterNewton','type' => 'integer','default' => '10','min' => '1'}}],'name' => 'if','type' => 'e','attrib' => {'expr' => '$UseNewton'}},{'content' => '
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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'NEWTON'}},{'content' => [{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Preconditioned','value' => 'prec','default' => 'T'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'No preconditioning','value' => 'free'}}],'name' => 'parameter','type' => 'e','attrib' => {'input' => 'select','name' => 'TypeJacobian','type' => 'string'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '1.e-5','name' => 'JacobianEps','type' => 'real','default' => '$doublePrecision ? 1.e-12 : 1.e-6','min' => '0'}},{'content' => '
#JACOBIAN
prec            TypeJacobian (prec, free)
1.E-12          JacobianEps

! The Jacobian matrix is always calculated with a matrix free approach,
! however it can be preconditioned  (\'prec\'), or not (\'free\')
! Default value is TypeJacobian=\'prec\'.
! JacobianEps contains the machine round off error for numerical derivatives.
! The default value is 1.E-12 for 8 byte reals and 1.E-6 for 4 byte reals.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'JACOBIAN'}},{'content' => [{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'left'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'symmetric','default' => 'T'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'right'}}],'name' => 'parameter','type' => 'e','attrib' => {'input' => 'select','name' => 'TypePrecondSide','type' => 'string'}},{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'MBILU','default' => 'T'}}],'name' => 'parameter','type' => 'e','attrib' => {'input' => 'select','name' => 'TypePrecond','type' => 'string'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '1','name' => 'GustafssonPar','type' => 'real','default' => '0.5','min' => '0'}},{'content' => '
#PRECONDITIONER
symmetric       TypePrecondSide (left, symmetric, right)
MBILU           TypePrecond (MBILU)
0.5             GustafssonPar (0. no modification, 1. full modification)

! Default parameters are shown. Right preconditioning does not affect
! the normalization of the residual. The Gustafsson parameter determines
! how much the MBILU preconditioner is modified. The default 0.5 value
! means a relaxed modification.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'PRECONDITIONER'}},{'content' => [{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'gmres','default' => 'T'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'bicgstab'}}],'name' => 'parameter','type' => 'e','attrib' => {'input' => 'select','name' => 'TypeKrylov','type' => 'string'}},{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => '0','value' => 'nul','default' => 'T'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'previous','value' => 'old'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'explicit'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'scaled explicit','value' => 'explicit'}}],'name' => 'parameter','type' => 'e','attrib' => {'input' => 'select','name' => 'TypeInitKrylov','type' => 'string'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '0.1','name' => 'ErrorMaxKrylov','type' => 'real','default' => '0.001','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'MaxMatvecKrylov','type' => 'integer','default' => '100','min' => '1'}},{'content' => '
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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'KRYLOV'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nKrylovVector','type' => 'integer','default' => 'MaxMatvecKrylov','min' => '1'}},{'content' => '
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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'KRYLOVSIZE'}}],'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'IMPLICIT PARAMETERS'}},{'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! STOPPING CRITERIA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

The commands in this group only work in stand alone mode.

','type' => 't'},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'MaxIteration','type' => 'integer','default' => '-1','min' => '-1'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'tSimulationMax','type' => 'real','default' => '-1','min' => '-1'}},{'content' => '

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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_IsStandAlone','name' => 'STOP','required' => '$_IsStandAlone'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoCheckStopFile','type' => 'logical','default' => 'T'}},{'content' => '

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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_IsStandAlone','name' => 'CHECKSTOPFILE'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CpuTimeMax','type' => 'real','default' => '-1','min' => '-1'}},{'content' => '

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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_IsStandAlone','name' => 'CPUTIMEMAX'}}],'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'STOPPING CRITERIA'}},{'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  OUTPUT PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
','type' => 't'},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'length' => '100','name' => 'NameRestartOutDir','type' => 'string','default' => 'GM/restartOUT'}},{'content' => '

#RESTARTOUTDIR
GM/restart_n5000	NameRestartOutDir

! The NameRestartOutDir variable contains the name of the directory
! where restart files are saved relative to the run directory.
! The directory should be inside the subdirectory with the name 
! of the component.
!
! Default value is "GM/restartOUT".
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'RESTARTOUTDIR'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoSaveRestart','type' => 'logical','default' => 'T'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DnSaveRestart','type' => 'integer','default' => '-1','min' => '-1'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DtSaveRestart','type' => 'real','default' => '-1','min' => '-1'}}],'name' => 'if','type' => 'e','attrib' => {'expr' => '$DoSaveRestart'}},{'content' => '
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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'SAVERESTART'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'length' => '100','name' => 'NamePlotDir','type' => 'string','default' => 'GM/IO2'}},{'content' => '

The NamePlotDir variable contains the name of the directory
where plot files and logfiles are saved relative to the run directory.
The directory should be inside the subdirectory with the name
of the component.

Default value is "GM/IO2".
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'PLOTDIR'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoSaveLogfile','type' => 'logical','default' => 'F'}},{'content' => [{'content' => [{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'MHD vars. dimensional','value' => 'MHD','default' => 'T'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Flux vars. dimensional','value' => 'FLX'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Raw vars. dimensional','value' => 'RAW'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Set vars. dimensional','value' => 'VAR'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'MHD vars. scaled','value' => 'mhd','default' => 'T'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Flux vars. scaled','value' => 'flx'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Raw vars. scaled','value' => 'raw'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Set vars. scaled','value' => 'var'}}],'name' => 'part','type' => 'e','attrib' => {'input' => 'select','name' => 'TypeLogVar','required' => 'T','type' => 'string'}},{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'exclusive' => 'T','name' => 'none'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'step'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'date'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'time'}}],'name' => 'part','type' => 'e','attrib' => {'input' => 'select','multiple' => 'T','name' => 'TypeTime','required' => 'F','type' => 'string'}}],'name' => 'parameter','type' => 'e','attrib' => {'max' => '4','name' => 'StringLog','type' => 'strings','min' => '1'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DnSaveLogfile','type' => 'integer','default' => '1','min' => '-1'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DtSaveLogfile','type' => 'real','default' => '-1','min' => '-1'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'if' => '$TypeLogVar =~ /var/i','length' => '100','name' => 'NameLogVars','type' => 'string'}},{'content' => [{'content' => [],'name' => 'part','type' => 'e','attrib' => {'multiple' => 'T','name' => 'LogRadii','type' => 'real','min' => '$rBody'}}],'name' => 'parameter','type' => 'e','attrib' => {'if' => '($TypeLogVar=~/flx/i or $NameLogVars=~/flx/i)','length' => '100','max' => '10','name' => 'StringLogRadii','type' => 'strings','min' => '1'}}],'name' => 'if','type' => 'e','attrib' => {'expr' => '$DoSaveLogfile'}},{'content' => '
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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'SAVELOGFILE'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nSatellite','type' => 'integer','default' => '0','min' => '0'}},{'content' => [{'content' => [{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'MHD vars. dimensional','value' => 'MHD','default' => 'T'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'All vars. dimensional','value' => 'FUL'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Set vars. dimensional','value' => 'VAR'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'MHD vars. scaled','value' => 'mhd'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'All vars. scaled','value' => 'ful'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Set vars. scaled','value' => 'var'}}],'name' => 'part','type' => 'e','attrib' => {'input' => 'select','name' => 'TypeSatelliteVar','required' => 'T','type' => 'string'}},{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'file','default' => 'T'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'equation','value' => 'eqn'}}],'name' => 'part','type' => 'e','attrib' => {'input' => 'select','name' => 'TypeTrajectory','required' => 'F','type' => 'string'}},{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'exclusive' => 'T','name' => 'none'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'step'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'date'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'time'}}],'name' => 'part','type' => 'e','attrib' => {'input' => 'select','multiple' => 'T','name' => 'TypeTime','required' => 'F','type' => 'string'}}],'name' => 'parameter','type' => 'e','attrib' => {'max' => '5','name' => 'StringSatellite','type' => 'strings','min' => '1'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DnOutput','type' => 'integer','default' => '1','min' => '-1'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DtOutput','type' => 'real','default' => '-1','min' => '-1'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'length' => '100','name' => 'NameTrajectoryFile','type' => 'string'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'if' => '$TypeSatelliteVar =~ /\\bvar\\b/i','length' => '100','name' => 'NameSatelliteVars','type' => 'string'}}],'name' => 'for','type' => 'e','attrib' => {'to' => '$nSatellite','from' => '1'}},{'content' => '
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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_IsFirstSession','name' => 'SATELLITE'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '100','name' => 'nPlotFile','type' => 'integer','default' => '0','min' => '0'}},{'content' => [{'content' => [{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'TECPLOT','value' => 'tec'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'IDL','value' => 'idl'}}],'name' => 'part','type' => 'e','attrib' => {'input' => 'select','name' => 'plotform','required' => 'T','type' => 'string'}},{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => '3D','value' => '3d/3d_'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'value' => 'x=0'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'value' => 'y=0','default' => 'T'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'value' => 'z=0'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'value' => 'sph'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'value' => 'los'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'value' => 'lin'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'if' => '$plotform =~ /\\bidl\\b/','value' => 'cut'}}],'name' => 'part','type' => 'e','attrib' => {'input' => 'select','name' => 'plotarea','required' => 'T','type' => 'string'}},{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'MHD vars. dimensional','value' => 'MHD'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'All vars. dimensional','value' => 'FUL'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Raw vars. dimensional','value' => 'RAW'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Ray tracing vars. dim.','value' => 'RAY'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Flux vars. dimensional','value' => 'FLX'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Solar vars. dimensional','value' => 'SOL'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'if' => '$plotarea eq \'lin\'','name' => 'Position vars. dimensional','value' => 'POS'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Select dimensional vars.','value' => 'VAR'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'MHD vars. scaled','value' => 'mhd'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'All vars. scaled','value' => 'ful'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Raw vars. scaled','value' => 'raw'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Ray tracing vars. scaled','value' => 'ray'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Flux vars. scaled','value' => 'flx'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Solar vars. scaled','value' => 'sol'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'if' => '$plotarea eq \'lin\'','name' => 'Position vars. scaled','value' => 'pos'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Select scaled vars.','value' => 'var'}}],'name' => 'part','type' => 'e','attrib' => {'input' => 'select','name' => 'plotvar','required' => 'T','type' => 'string'}}],'name' => 'parameter','type' => 'e','attrib' => {'max' => '3','name' => 'StringPlot','type' => 'strings','min' => '3'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DnSavePlot','type' => 'integer','min' => '-1'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DtSavePlot','type' => 'real','min' => '-1'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'xMinCut','type' => 'real'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'xMaxCut','type' => 'real','min' => '$xMinCut'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'yMinCut','type' => 'real'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'yMaxCut','type' => 'real','min' => '$yMinCut'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'zMinCut','type' => 'real'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'zMaxCut','type' => 'real','min' => '$zMinCut'}}],'name' => 'if','type' => 'e','attrib' => {'expr' => '$plotarea =~ /\\bcut\\b/'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'if' => '$plotarea =~ /\\bsph\\b/','name' => 'Radius','type' => 'real','default' => '10','min' => '0'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'LosVectorX','type' => 'real','default' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'LosVectorY','type' => 'real','default' => '0.0001'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'LosVectorZ','type' => 'real','default' => '1'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'xSizeImage','type' => 'real','default' => '20','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'ySizeImage','type' => 'real','default' => '20','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'xOffset','type' => 'real','default' => '10'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'yOffset','type' => 'real','default' => '10'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'rOccult','type' => 'real','default' => '2.5','min' => '1'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '1','name' => 'MuLimbDarkening','type' => 'real','default' => '0.5','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nPixX','type' => 'integer','default' => '200','min' => '2'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nPixY','type' => 'integer','default' => '200','min' => '2'}}],'name' => 'if','type' => 'e','attrib' => {'expr' => '$plotarea =~ /\\blos\\b/'}},{'content' => [{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Advected B','value' => 'A'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'B','default' => 'T'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'U'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'J'}}],'name' => 'parameter','type' => 'e','attrib' => {'input' => 'select','name' => 'NameLine','type' => 'string'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'IsSingleLine','type' => 'logical','default' => 'F'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '20','name' => 'nLine','type' => 'integer','default' => '1','min' => '1'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'xStartLine','type' => 'real'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'yStartLine','type' => 'real'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'zStartLine','type' => 'real'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'IsParallel','type' => 'logical'}}],'name' => 'for','type' => 'e','attrib' => {'to' => '$nLine','from' => '1'}}],'name' => 'if','type' => 'e','attrib' => {'expr' => '$plotarea =~ /\\blin\\b/'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'if' => '($plotform=~/\\bidl\\b/ and $plotarea!~/\\b(sph|los|lin)\\b/)','name' => 'DxSavePlot','type' => 'real','default' => '-1.0','min' => '-1.0'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'length' => '100','name' => 'NameVars','type' => 'string'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'length' => '100','name' => 'NamePars','type' => 'string'}}],'name' => 'if','type' => 'e','attrib' => {'expr' => '$plotvar =~ /\\bvar\\b/i'}}],'name' => 'for','type' => 'e','attrib' => {'to' => '$nPlotFile','from' => '1'}},{'content' => '
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

','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'SAVEPLOT'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoSaveBinary','type' => 'logical','default' => 'T'}},{'content' => '
#SAVEBINARY
T			DoSaveBinary   used only for \'idl\' plot file

! Default is .true. Saves unformatted IO2/*.idl files if true. 
! This is the recommended method, because it is fast and accurate.
! The only advantage of saving IO2/*.idl in formatted text files is
! that it can be processed on another machine or with a different 
! (lower) precision. For example PostIDL.exe may be compiled with 
! single precision to make IO2/*.out files smaller, while BATSRUS.exe is 
! compiled in double precision, to make results more accurate.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'SAVEBINARY'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoSavePlotsAmr','type' => 'logical','default' => 'F'}},{'content' => '
#SAVEPLOTSAMR
F			DoSavePlotsAmr

! Save plots before each AMR. Default is DoSavePlotsAMR=.false.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'SAVEPLOTSAMR'}}],'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'OUTPUT PARAMETERS'}},{'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  AMR PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
','type' => 't'},{'content' => [{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'default','default' => '1'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'all'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'none'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => '3Dbodyfocus'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'spherefocus'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'magnetosphere'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'points'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'coupledhelio'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'helio_init'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'helio_z=4'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'all_then_focus'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'cme'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'points'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'mag_new'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'magnetosphere'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'magneto_fine'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'magneto12'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'magnetosaturn'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'magnetojupiter'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'paleo'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'comet'}}],'name' => 'parameter','type' => 'e','attrib' => {'input' => 'select','name' => 'InitialRefineType','type' => 'string'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'InitialRefineLevel','type' => 'integer','default' => '4','min' => '0'}},{'content' => '
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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_IsFirstSession','name' => 'AMRINIT'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nRefineLevelIC','type' => 'integer','default' => '0','min' => '0'}},{'content' => '
#AMRINITPHYSICS
3			nRefineLevelIC

! Defines number of physics (initial condition) based AMR-s AFTER the 
! geometry based initial AMR-s defined by #AMRINIT were done.
! Only useful if the initial condition has a non-trivial analytic form.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_IsFirstSession','name' => 'AMRINITPHYSICS'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'MinBlockLevel','type' => 'integer','default' => '0','min' => '-1'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'MaxBlockLevel','type' => 'integer','default' => '99','min' => '-1'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoFixBodyLevel','type' => 'logical','default' => 'F'}},{'content' => '
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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'AMRLEVELS'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DxCellMin','type' => 'real','default' => '0','min' => '-1'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DxCellMax','type' => 'real','default' => '99999','min' => '-1'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoFixBodyLevel','type' => 'logical','default' => 'F'}},{'content' => '
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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'AMRRESOLUTION'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DnRefine','type' => 'integer','default' => '-1','min' => '-1'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoAutoRefine','type' => 'logical','default' => 'F'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '100','name' => 'PercentCoarsen','type' => 'real','default' => '20','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '100','name' => 'PercentRefine','type' => 'real','default' => '20','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'MaxTotalBlocks','type' => 'integer','default' => '99999','min' => '1'}}],'name' => 'if','type' => 'e','attrib' => {'expr' => '$DoAutoRefine'}}],'name' => 'if','type' => 'e','attrib' => {'expr' => '$DnRefine>0'}},{'content' => '
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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'AMR'}},{'content' => [{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => '1'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => '2'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => '3','default' => '1'}}],'name' => 'parameter','type' => 'e','attrib' => {'input' => 'select','name' => 'nRefineCrit','type' => 'integer'}},{'content' => [{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'grad T','value' => 'gradt/gradT'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'grad P','value' => 'gradp/gradP'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'grad log(Rho)','value' => 'gradlogrho'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'grad log(p)','value' => 'gradlogP/gradlogp'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'grad E','value' => 'gradE'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'curl U','value' => 'curlV/curlv/curlU/curlu'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'curl B','value' => 'curlB/curlb'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'div U','value' => 'divU/divu/divV/divv'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'divB','value' => 'divb/divB'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'vAlfven','value' => 'Valfven/vAlfven/valfven'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'heliospheric beta','value' => 'heliobeta'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'flux'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'heliospheric current sheet','value' => 'heliocurrentsheet'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'rCurrents','value' => 'rcurrents/Rcurrents'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Transient','value' => 'transient/Transient'}}],'name' => 'parameter','type' => 'e','attrib' => {'input' => 'select','name' => 'TypeRefine','type' => 'string'}},{'content' => [{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'P_dot','value' => 'p_dot/P_dot'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'T_dot','value' => 't_dot/T_dot'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Rho_dot','value' => 'rho_dot/Rho_dot','default' => 'T'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'RhoU_dot','value' => 'RhoU_dot/rhou_dot'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Rho_2nd_1','value' => 'Rho_2nd_1/rho_2nd_1'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Rho_2nd_2','value' => 'Rho_2nd_2/rho_2nd_2'}}],'name' => 'parameter','type' => 'e','attrib' => {'input' => 'select','name' => 'TypeTransient','type' => 'string'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseSunEarth','type' => 'logical','default' => 'F'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'xEarth','type' => 'real'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'yEarth','type' => 'real'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'zEarth','type' => 'real'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'InvD2Ray','type' => 'real'}}],'name' => 'if','type' => 'e','attrib' => {'expr' => '$UseSunEarth'}}],'name' => 'if','type' => 'e','attrib' => {'expr' => '$TypeRefine =~ /transient/i'}}],'name' => 'for','type' => 'e','attrib' => {'to' => '$nRefineCrit','from' => '1'}},{'content' => '
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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'AMRCRITERIA'}}],'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'AMR PARAMETERS'}},{'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  SCHEME PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
','type' => 't'},{'content' => [{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => '1','default' => 'T'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => '2'}}],'name' => 'parameter','type' => 'e','attrib' => {'input' => 'select','name' => 'nOrder','type' => 'integer'}},{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Roe','value' => 'Roe/roe/1'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Rusanov','value' => 'Rusanov/rusanov/2/TVDLF','default' => 'T'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Linde','value' => 'Linde/linde/3/HLLEL'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Sokolov','value' => 'Sokolov/sokolov/4/AW'}}],'name' => 'parameter','type' => 'e','attrib' => {'input' => 'select','name' => 'TypeFlux','type' => 'string'}},{'content' => [{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'minmod','default' => 'T'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'beta'}}],'name' => 'parameter','type' => 'e','attrib' => {'input' => 'select','name' => 'TypeLimiter','type' => 'string'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'if' => '$TypeLimiter eq \'beta\'','max' => '2','name' => 'LimiterBeta','type' => 'real','default' => '1.2','min' => '1'}}],'name' => 'if','type' => 'e','attrib' => {'expr' => '$nOrder == 2'}},{'content' => '
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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'SCHEME'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseNonConservative','type' => 'logical','default' => 'T'}},{'content' => '
#NONCONSERVATIVE
T		UseNonConservative

! For Earth the default is using non-conservative equations 
! (close to the body).
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'NONCONSERVATIVE'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '3','name' => 'nConservCrit','type' => 'integer','default' => '1','min' => '0'}},{'content' => [{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'radius','value' => 'r/R/radius/Radius','default' => 'T'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'p','value' => 'p/P'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'grad P','value' => 'gradp/GradP'}}],'name' => 'parameter','type' => 'e','attrib' => {'input' => 'select','name' => 'TypeConservCrit','type' => 'string'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'if' => '$TypeConservCrit =~ /^r|radius$/i','name' => 'rConserv','type' => 'real','default' => '2*$rBody','min' => '$rBody'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'if' => '$TypeConservCrit =~ /^p$/i','name' => 'pCoeffConserv','type' => 'real','default' => '0.05','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'if' => '$TypeConservCrit =~ /gradp/i','name' => 'GradPCoeffConserv','type' => 'real','default' => '0.1','min' => '0'}}],'name' => 'for','type' => 'e','attrib' => {'to' => '$nConservCrit','from' => '1'}},{'content' => '
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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'CONSERVATIVECRITERIA'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseUpdateCheck','type' => 'logical','default' => 'T'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '100','name' => 'RhoMinPercent','type' => 'real','default' => '40','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'RhoMaxPercent','type' => 'real','default' => '400','min' => '100'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '100','name' => 'pMinPercent','type' => 'real','default' => '40','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'pMaxPercent','type' => 'real','default' => '400','min' => '100'}}],'name' => 'if','type' => 'e','attrib' => {'expr' => '$UseUpdateCheck'}},{'content' => '
#UPDATECHECK
T			UseUpdateCheck
40.			RhoMinPercent
400.			RhoMaxPercent
40.			pMinPercent
400.			pMaxPercent

! Default values are shown.  This will adjust the timestep so that
! density and pressure cannot change by more than the given percentages
! in a single timestep.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'UPDATECHECK'}},{'content' => [{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => '1','default' => 'T'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => '2'}}],'name' => 'parameter','type' => 'e','attrib' => {'input' => 'select','name' => 'nOrderProlong','type' => 'integer'}},{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'left-right','value' => 'lr','default' => 'T'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'central'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'minmod'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'left-right extrapolate','value' => 'lr2'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'central    extrapolate','value' => 'central2'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'minmod     extrapolate','value' => 'minmod2'}}],'name' => 'parameter','type' => 'e','attrib' => {'if' => '$nOrderProlong==2','input' => 'select','name' => 'TypeProlong','type' => 'string'}},{'content' => '
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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'PROLONGATION'}},{'content' => [{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'm_p_cell FACES ONLY','value' => 'allopt','default' => 'T'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'm_p_cell','value' => 'all'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'm_p_dir FACES ONLY','value' => 'opt'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'm_p_dir group by directions','value' => 'dir'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'm_p_dir group by faces     ','value' => 'face'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'm_p_dir group by kind and face','value' => 'min'}}],'name' => 'parameter','type' => 'e','attrib' => {'input' => 'select','name' => 'TypeMessagePass','type' => 'string'}},{'content' => '
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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'MESSAGEPASS','alias' => 'OPTIMIZE'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseBorisCorrection','type' => 'logical','default' => 'F'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'if' => '$UseBorisCorrection','max' => '1','name' => 'BorisClightFactor','type' => 'real','default' => '1','min' => '0'}},{'content' => '
#BORIS
T			UseBorisCorrection
1.0			BorisClightFactor !Only if UseBorisCorrection is true

! Default is boris_correction=.false.
! Use semi-relativistic MHD equations with speed of light reduced by
! the BorisClightFactor. Set BorisClightFactor=1.0 for true semi-relativistic
! MHD. Gives the same steady state as normal MHD analytically, but there
! can be differences due to discretization errors. 
! You can use either Boris or BorisSimple but not both. 
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'BORIS'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseBorisSimple','type' => 'logical','default' => 'F'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'if' => '$UseBorisSimple','max' => '1','name' => 'BorisClightFactor','type' => 'real','default' => '1','min' => '0'}},{'content' => '
#BORISSIMPLE
T			UseBorisSimple
0.05			BorisClightFactor !Only if UseBorisSimple is true

! Default is UseBorisSimple=.false. 
! Use simplified semi-relativistic MHD with speed of light reduced by the
! BorisClightFactor. This is only useful with BorisClightFactor less than 1.
! Should give the same steady state as normal MHD, but there can be a
! difference due to discretization errors.
! You can use either Boris or BorisSimple but not both. 
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'BORISSIMPLE'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseDivbSource','type' => 'logical','default' => 'T'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseDivbDiffusion','type' => 'logical','default' => 'F'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseProjection','type' => 'logical','default' => 'F'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseConstrainB','type' => 'logical','default' => 'F'}},{'content' => [{'content' => '
	! At least one of the options should be true.
	','type' => 't'}],'name' => 'rule','type' => 'e','attrib' => {'expr' => '$UseDivbSource or $UseDivbDiffusion or $UseProjection or $UseConstrainB'}},{'content' => [{'content' => '
	! If UseProjection is true, all others should be false.
	','type' => 't'}],'name' => 'rule','type' => 'e','attrib' => {'expr' => 'not($UseProjection and ($UseDivbSource or $UseDivbDiffusion or $UseConstrainB))'}},{'content' => [{'content' => '
	! If UseConstrainB is true, all others should be false.
	','type' => 't'}],'name' => 'rule','type' => 'e','attrib' => {'expr' => 'not($UseConstrainB and ($UseDivbSource or $UseDivbDiffusion or $UseProjection))'}},{'content' => '
	
#DIVB
T			UseDivbSource
F			UseDivbDiffusion	
F			UseProjection           
F			UseConstrainB           

! Default values are shown above.
! If UseProjection is true, all others should be false.
! If UseConstrainB is true, all others should be false.
! At least one of the options should be true.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'DIVB'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseB0Source','type' => 'logical','default' => 'T'}},{'content' => '
#DIVBSOURCE
T			UseB0Source

! Add extra source terms related to the non-zero divergence and curl of B0.
! Default is true.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'DIVBSOURCE'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '1','name' => 'DivbDiffCoeff','type' => 'real','default' => '0.1666667','min' => '0'}},{'content' => '
#DIVBDIFFUSION
0.1666667		DivbDiffCoeff

! Default value is shown above. 1.0/6.0
! If divb_diffcoeff .gt. 0.5 then cfl .lt. 0.5/DivbDiffCoeff is required!
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'DIVBDIFFUSION'}},{'content' => [{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Conjugate Gradients','value' => 'cg','default' => 'T'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'BiCGSTAB','value' => 'bicgstab'}}],'name' => 'parameter','type' => 'e','attrib' => {'input' => 'select','name' => 'TypeProjectIter','type' => 'string'}},{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Relative norm','value' => 'rel','default' => 'T'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Maximum error','value' => 'max'}}],'name' => 'parameter','type' => 'e','attrib' => {'input' => 'select','name' => 'TypeProjectStop','type' => 'string'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '1','name' => 'RelativeLimit','type' => 'real','default' => '0.1','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'AbsoluteLimit','type' => 'real','default' => '0.0','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'MaxMatvec','type' => 'integer','default' => '50','min' => '1'}},{'content' => '
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

','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'PROJECTION'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '1','name' => 'pRatioLow','type' => 'real','default' => '0.01','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '1','name' => 'pRatioHigh','type' => 'real','default' => '0.1','min' => '$pRatioLow'}},{'content' => '
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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'CORRECTP'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseAccurateIntegral','type' => 'logical','default' => 'T'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseAccurateTrace','type' => 'logical','default' => 'F'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '60','name' => 'DtExchangeRay','type' => 'real','default' => '0.1','min' => '0.01'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DnRaytrace','type' => 'integer','default' => '1','min' => '1'}},{'content' => '

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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'RAYTRACE'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TauCoupleIm','type' => 'real','min' => '0'}},{'content' => '

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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'IM'}}],'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'SCHEME PARAMETERS'}},{'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!  PHYSICS PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
','type' => 't'},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Gamma','type' => 'real','default' => '1.6666666667','min' => '1'}},{'content' => '
#GAMMA
1.6666666667		Gamma

! The adiabatic index (ratio of the specific heats for fixed pressure
! and volume. The default value is 5.0/3.0, which is valid for
! monoatomic gas or plasma.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_IsFirstSession','name' => 'GAMMA'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'RhoLeft','type' => 'real','default' => '1','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UnLeft','type' => 'real','default' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Ut1Left','type' => 'real','default' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Ut2Left','type' => 'real','default' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'BnLeft','type' => 'real','default' => '0.75'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Bt1Left','type' => 'real','default' => '1'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Bt2Left','type' => 'real','default' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'pRight','type' => 'real','default' => '1','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'RhoRight','type' => 'real','default' => '0.125','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UnRight','type' => 'real','default' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Ut1Right','type' => 'real','default' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Ut2Right','type' => 'real','default' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'BnRight','type' => 'real','default' => '0.75'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Bt1Right','type' => 'real','default' => '-1'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Bt2Right','type' => 'real','default' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'pRight','type' => 'real','default' => '0.1','min' => '0'}},{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'no rotation','value' => '0','default' => 'T'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'value' => '0.25'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'value' => '0.3333333333333'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'value' => '0.5'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'value' => '1'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'value' => '2'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'value' => '3'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'value' => '4'}}],'name' => 'parameter','type' => 'e','attrib' => {'input' => 'select','name' => 'ShockSlope','type' => 'real'}},{'content' => '
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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'SHOCKTUBE'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'SwRhoDim','type' => 'real','default' => '5','min' => '-1'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'SwTDim','type' => 'real','default' => '181712.175','min' => '-1'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '0','name' => 'SwUxDim','type' => 'real','default' => '-400'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'SwUyDim','type' => 'real','default' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'SwUzDim','type' => 'real','default' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'SwBxDim','type' => 'real','default' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'SwByDim','type' => 'real','default' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'SwBzDim','type' => 'real','default' => '5'}},{'content' => '
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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_IsFirstSession','name' => 'SOLARWIND'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseUpstreamInputFile','type' => 'logical','default' => 'F'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'length' => '100','name' => 'NameUpstreamFile','type' => 'string'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'SatelliteYPos','type' => 'real','default' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'SatelliteZPos','type' => 'real','default' => '0'}}],'name' => 'if','type' => 'e','attrib' => {'expr' => '$UseUpstreamInputFile'}},{'content' => '
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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'UPSTREAM_INPUT_FILE'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseBody','type' => 'logical','default' => 'F'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'rBody','type' => 'real','default' => '3','min' => '0'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'rCurrents','type' => 'real','default' => '4','min' => '-1'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'BodyRhoDim','type' => 'real','default' => '1','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'BodyTDim','type' => 'real','default' => '10000','min' => '0'}}],'name' => 'if','type' => 'e','attrib' => {'expr' => '$_NameComp eq \'GM\''}}],'name' => 'if','type' => 'e','attrib' => {'expr' => '$UseBody'}},{'content' => '
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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_IsFirstSession','name' => 'BODY','alias' => 'MAGNETOSPHERE'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseGravity','type' => 'logical','default' => 'F'}},{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'central mass','value' => '0','default' => 'T'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'X direction','value' => '1'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Y direction','value' => '2'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Z direction','value' => '3'}}],'name' => 'parameter','type' => 'e','attrib' => {'if' => '$UseGravity','input' => 'select','name' => 'iDirGravity','type' => 'integer'}},{'content' => '
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

','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_IsFirstSession','name' => 'GRAVITY'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseMassLoading','type' => 'logical','default' => 'F'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoAccelerateMassLoading','type' => 'logical','default' => 'F'}},{'content' => '
#MASSLOADING
F			UseMassLoading
F			DoAccelerateMassLoading
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'MASSLOADING'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseHeatFlux','type' => 'logical','default' => 'F'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseSpitzerForm','type' => 'logical','default' => 'T'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Kappa0Heat','type' => 'real','default' => '1.23E-11'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Kappa0Heat','type' => 'real','default' => '2.5'}}],'name' => 'if','type' => 'e','attrib' => {'expr' => 'not $UseSpitzerForm'}},{'content' => '
#HEATFLUX
T		UseHeatFlux
F		UseSpitzerForm
1.23E-11	Kappa0Heat [W/m/K]	! Only if not UseSpitzerForm
2.50E+00	ExponentHeat [-]	! Only if not UseSpitzerForm
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'HEATFLUX'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseResistFlux','type' => 'logical','default' => 'F'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseSpitzerForm','type' => 'logical','default' => 'T'}},{'content' => [{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'localized','value' => 'Localized/localized'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'constant','value' => 'Constant/constant','default' => 'T'}}],'name' => 'parameter','type' => 'e','attrib' => {'input' => 'select','name' => 'TypeResist','type' => 'string'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Eta0Resist','type' => 'real','default' => '9.69953E+8'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Alpha0Resist','type' => 'real','default' => '150'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'yShiftResist','type' => 'real','default' => '0.5'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TimeInitRise','type' => 'real','default' => '0.05'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TimeConstLev','type' => 'real','default' => '1'}}],'name' => 'if','type' => 'e','attrib' => {'expr' => '$TypeResist =~ /localized/i'}}],'name' => 'if','type' => 'e','attrib' => {'expr' => 'not $UseSpitzerForm'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseAnomResist','type' => 'logical','default' => 'F'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Eta0AnomResist','type' => 'real','default' => '1.93991E+09'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'EtaAnomMaxResist','type' => 'real','default' => '1.93991E+10'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'ThresholdFactorResist','type' => 'real','default' => '1'}}],'name' => 'if','type' => 'e','attrib' => {'expr' => '$UseAnomResist'}},{'content' => '
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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'RESISTIVEFLUX'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseDefaultUnits','type' => 'logical','default' => 'T'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Grav0Diss','type' => 'real','default' => '2.635620E-02'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Beta0Diss','type' => 'real','default' => '1.640000E-01'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Length0Diss','type' => 'real','default' => '1.500000E+06'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Time0Diss','type' => 'real','default' => '1.159850E+01'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Rho0Diss','type' => 'real','default' => '5.019000E-11'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Tem0Diss','type' => 'real','default' => '1.000000E+05'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Theta0Diss','type' => 'real','default' => '6.000000E-01'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Delta0Diss','type' => 'real','default' => '2.500000E+01'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'EpsilonDiss','type' => 'real','default' => '7.000000E+00'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'RhoDifDiss','type' => 'real','default' => '4.500000E+00'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'yShiftDiss','type' => 'real','default' => '4.000000E-01'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'ScaleHeightDiss','type' => 'real','default' => '5.000000E-01'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'ScaleFactorDiss','type' => 'real','default' => '1.159850E+01'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'BZ0iss','type' => 'real','default' => '5.000000E-01'}},{'content' => '
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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_IsFirstSession','name' => 'TESTDISSMHD'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseBody2','type' => 'logical','default' => 'F'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'rBody2','type' => 'real','default' => '0.1','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '$xMax','name' => 'xBody2','type' => 'real','default' => '-40','min' => '$xMin'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '$yMax','name' => 'yBody2','type' => 'real','default' => '0','min' => '$yMin'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '$zMax','name' => 'zBody2','type' => 'real','default' => '0','min' => '$zMin'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'rCurrents2','type' => 'real','default' => '1.3*$rBody2','min' => '$rBody2'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'RhoDimBody2','type' => 'real','default' => '5','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'tDimBody2','type' => 'real','default' => '25000','min' => '0'}}],'name' => 'if','type' => 'e','attrib' => {'expr' => '$UseBody2'}},{'content' => '

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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_IsFirstSession','name' => 'SECONDBODY'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'BdpDimBody2x','type' => 'real'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'BdpDimBody2y','type' => 'real'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'BdpDimBody2z','type' => 'real'}},{'content' => '

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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_IsFirstSession','name' => 'DIPOLEBODY2'}}],'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'PHYSICS PARAMETERS'}},{'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!! SOLAR PROBLEM TYPES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
','type' => 't'},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'BodyTDim','type' => 'real','default' => '2.85E06','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'BodyRhoDim','type' => 'real','default' => '1.50E8','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'qSun','type' => 'real','default' => '25.0','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'tHeat','type' => 'real','default' => '1.75','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'rHeat','type' => 'real','default' => '1.0','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'SigmaHeat','type' => 'real','default' => '4.5','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoInitRope','type' => 'logical','default' => 'F'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CmeA','type' => 'real','default' => '0.7','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CmeR1','type' => 'real','default' => '1.2','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CmeR0','type' => 'real','default' => '1.0','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CmeA1','type' => 'real','default' => '0.23','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CmeAlpha','type' => 'real','default' => '0.0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CmeRho1','type' => 'real','default' => '2.5E-12','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CmeRho2','type' => 'real','default' => '2.0E-13','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '10','name' => 'ModulationRho','type' => 'real','default' => '0.0','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '10','name' => 'ModulationP','type' => 'real','default' => '0.0','min' => '0'}}],'name' => 'if','type' => 'e','attrib' => {'expr' => '$DoInitRope'}},{'content' => '
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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_IsFirstSession and $_NameComp ne \'GM\'','name' => 'HELIOSPHERE'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DtUpdateB0','type' => 'real','default' => '0.0001','min' => '-1'}},{'content' => '

#HELIOUPDATEB0
-1.0			DtUpdateB0 [s]

Set the frequency of updating the B0 field for the solar corona.
A negative value means that the B0 field is not updated.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_NameComp ne \'GM\'','name' => 'HELIOUPDATEB0'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'HelioDipoleStrength','type' => 'real'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '90','name' => 'HelioDipoleTilt','type' => 'real','default' => '0','min' => '-90'}},{'content' => '

#HELIODIPOLE
-3.0                    HelioDipoleStrength [G]
 0.0                    HelioDipoleTilt     [deg]

! Variable HelioDipoleStrength defines the equatorial field strength in Gauss,
! while HelioDipoleTilt is the tilt relative to the ecliptic North 
! (negative sign means towards the planet) in degrees.
!
! Default value is HelioDipoleStrength = 0.0.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_NameComp ne \'GM\'','name' => 'HELIODIPOLE'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseInertialFrame','type' => 'logical','default' => 'T'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseRotatingBC','type' => 'logical','default' => 'T'}}],'name' => 'if','type' => 'e','attrib' => {'expr' => '$UseInertialFrame'}},{'content' => '

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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_IsFirstSession and $_NameComp ne \'GM\'','name' => 'HELIOROTATION','alias' => 'INERTIAL'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoSendMHD','type' => 'logical','default' => 'F'}},{'content' => '

#HELIOTEST
F			DoSendMHD

! If DoSendMHD is true, IH sends the real MHD solution to GM in the coupling.
! If DoSendMHD is false then the values read from the IMF file are sent,
! so there is no real coupling. Mostly used for testing the framework.
!
! Default value is true, ie. real coupling.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_NameComp ne \'GM\'','name' => 'HELIOTEST'}},{'content' => [{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => 'Low','default' => 'T'}}],'name' => 'parameter','type' => 'e','attrib' => {'input' => 'select','name' => 'TypeCme','type' => 'string'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CmeA','type' => 'real','default' => '0.7','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CmeR1','type' => 'real','default' => '1.2','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CmeR0','type' => 'real','default' => '1.0','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CmeA1','type' => 'real','default' => '0.23','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CmeAlpha','type' => 'real','default' => '0.0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CmeRho1','type' => 'real','default' => '2.5E-12','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CmeRho2','type' => 'real','default' => '2.0E-13','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CmeB1Dim','type' => 'real','default' => '1.0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CmeUErupt','type' => 'real','default' => '4.0E5','min' => '0'}},{'content' => '
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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_IsFirstSession and $_NameComp ne \'GM\'','name' => 'CME'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'tArcDim','type' => 'real','default' => '1.0E6','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'RhoArcDim','type' => 'real','default' => '1.0E-12','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'bArcDim','type' => 'real','default' => '0.718144','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'ByArcDim','type' => 'real','default' => '1.0E6'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UzArcDim','type' => 'real','default' => '5.0E3'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Phi0Arc','type' => 'real','default' => '0.5'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'MuArc','type' => 'real','default' => '1.3'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'ExpArc','type' => 'real','default' => '3','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'WidthArc','type' => 'real','default' => '0.5','min' => '0'}},{'content' => '
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
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_IsFirstSession and $_NameComp ne \'GM\'','name' => 'ARCADE'}}],'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'SOLAR PROBLEM TYPES'}},{'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!! COMET PROBLEM TYPE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
','type' => 't'},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'ProdRate','type' => 'real','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UrNeutral','type' => 'real','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'AverageMass','type' => 'real','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'IonizationRate','type' => 'real','min' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'kFriction','type' => 'real','min' => '0'}},{'content' => '
#COMET
1.0E28		ProdRate    - Production rate (#/s)
1.0		UrNeutral   - neutral radial outflow velocity (km/s)
17.0		AverageMass - average particle mass (amu)
1.0E-6		IonizationRate (1/s)
1.7E-9		kFriction - ion-neutral friction rate coefficient (cm^3/s)

! Only used by problem_comet.  Defaults are as shown.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_IsFirstSession','name' => 'COMET'}}],'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'COMET PROBLEM TYPE'}},{'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! SCRIPT COMMANDS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
','type' => 't'},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'length' => '100','name' => 'NameIncludeFile','type' => 'string','default' => 'Param/'}},{'content' => '

#INCLUDE
Param/SSS_3000		NameIncludeFile

! Include a library file from Param/ or any file from anywhere else.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'INCLUDE'}},{'content' => [{'content' => '

#RUN

This command is only used in stand alone mode.

Run BATSRUS with the parameters above and then return for the next session
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_IsStandAlone','name' => 'RUN'}},{'content' => [{'content' => '

#END

This command is only used in stand alone mode.

Run the executable with the parameters above and then stop.
In included files #END simply means the end of the included lines.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'if' => '$_IsStandAlone','name' => 'END'}}],'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'SCRIPT COMMANDS'}},{'content' => [{'content' => '
	Either command #SOLARWIND or #UPSTREAM_INPUT_FILE must be used!
','type' => 't'}],'name' => 'rule','type' => 'e','attrib' => {'expr' => '($SwRhoDim > 0) or $UseUpstreamInputFile or $_NameComp ne \'GM\''}},{'content' => [{'content' => '
	Part implicit scheme requires more than 1 implicit block!
','type' => 't'}],'name' => 'rule','type' => 'e','attrib' => {'expr' => '$MaxImplBlock>1 or not $UsePartImplicit or not $MaxImplBlock'}},{'content' => [{'content' => '
	Full implicit scheme should be used with equal number of 
	explicit and implicit blocks!
','type' => 't'}],'name' => 'rule','type' => 'e','attrib' => {'expr' => '$MaxImplBlock==$MaxBlock or not $UseFullImplicit'}}],'name' => 'commandList','type' => 'e','attrib' => {'name' => 'BATSRUS: GM, SC and IH Components'}}];