#^CFG FILE _FALSE_
$tree = [{'content' => [{'content' => '

GITM is typically run with a {\\tt UAM.in} file, which resides in the
directory that you are running from. In the framework it obtains 
parameters from the PARAM.in file.

','type' => 't'},{'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

In order to run GITM, a starting time and an ending time must be
specified.  These are specified using the following commands:

','type' => 't'},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'iYear','type' => 'integer','default' => '2000'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '12','name' => 'iMonth','type' => 'integer','min' => '1','default' => '3'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '31','name' => 'iDay','type' => 'integer','min' => '1','default' => '21'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '23','name' => 'iHour','type' => 'integer','min' => '0','default' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '59','name' => 'iMinute','type' => 'integer','min' => '0','default' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '59','name' => 'iSecond','type' => 'integer','min' => '0','default' => '0'}},{'content' => '

#TIMESTART
1999			iYear
03			iMonth
18			iDay
00			iHour
00			iMinute
00			iSecond

This command is only used in the stand alone mode.

The #STARTTIME command sets the initial date and time for the simulation.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'TIMESTART','if' => '$_IsStandAlone','alias' => 'STARTTIME'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'iYear','type' => 'integer','default' => '2000'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '12','name' => 'iMonth','type' => 'integer','min' => '1','default' => '3'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '31','name' => 'iDay','type' => 'integer','min' => '1','default' => '21'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '23','name' => 'iHour','type' => 'integer','min' => '0','default' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '59','name' => 'iMinute','type' => 'integer','min' => '0','default' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '59','name' => 'iSecond','type' => 'integer','min' => '0','default' => '0'}},{'content' => '

#TIMEEND
1999			iYear
03			iMonth
25			iDay
00			iHour
00			iMinute
00			iDay

This command is only used in the stand alone mode.

The #TIMEEND command sets the final date and time for the simulation.
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'TIMEEND','if' => '$_IsStandAlone','alias' => 'ENDTIME'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DoRestart','type' => 'logical','default' => 'F'}},{'content' => [{'content' => '
                File UA/RestartIN/header.rst should exist!
	','type' => 't'}],'name' => 'rule','type' => 'e','attrib' => {'expr' => '-f \'UA/RestartIN/header.rst\' or not $DoRestart'}},{'content' => '

#RESTART
F			DoRestart

There are two commands that are typically not input by a user, but are
specified in the restart header file that is read in the exact same
way as the input file.  It is possible to set these variables, though.

','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'RESTART'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'iStep','type' => 'integer','default' => ''}},{'content' => '

#ISTEP
1

This sets the current iStep to the read in value instead of starting at
iteration 1.

','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'ISTEP'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'tsimulation','type' => 'real','default' => ''}},{'content' => '

#TSIMULATION
0.0

This offsets the current start time by the read in amount.  It is simply
added to the #STARTTIME input time.

','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'TSIMULATION'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'CPUTimeMax','type' => 'real, seconds','default' => ''}},{'content' => '

#CPUTIMEMAX
7200.0

When you are running on a queue-based system, you can use this command
to set the exact amount of time that the code should run, then stop
and write a restart file.  It is good to give a small buffer so the
code has time to write files before the queue stops.  This buffer time
is quite dependent upon the system.  On fast I/O machines, I typically
give a buffer of only a couple of minutes.  On some systems, I
sometimes give a full half hour, just to make absolutely sure the code
will write all of the correct files and exit before the queue is up.

','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'CPUTIMEMAX'}}],'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'Time Variables'}},{'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

','type' => 't'},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseMSIS','type' => 'logical','default' => ''}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseIRI','type' => 'logical','default' => ''}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TempMin','type' => 'real','default' => ''}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'TempMax','type' => 'real','default' => ''}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Number Density','type' => 'real','min' => '0.01','default' => ''}}],'name' => 'for','type' => 'e','attrib' => {'to' => '3','from' => '1'}}],'name' => 'if','type' => 'e','attrib' => {'expr' => '$UseMSIS'}},{'content' => '

#INITIAL
F			UseMSIS
T			UseIRI
200.0			TBottom
1200.0			TTop
5.0e17			NDensity1
7.0e18			NDensity2
3.0e19			NDensity3

On the Earth, empirical models exist which can be used to derive a
background atmosphere and ionosphere.  These are MSIS (thermosphere
and IRI (ionosphere).  If MSIS is used, then the all of the species
densities are set using MSIS.  There are 2 species which MSIS does not
include: [NO] and [N($^2$D)].  We have made up some formula for
setting these two species.  The neutral temperature is also set using
MSIS if UseMSIS = T.

If UseMSIS = F, then GITM reads in the temperature at the bottom
and top of the atmosphere (for the initial condition), and the number
density at the bottom of the atmosphere for all of the major species
(nSpecies, which is set in ModPlanet.f90).

It UseIRI = T, the number densities of the ion species are set
by IRI.  If it is .false., then the initial densities are set to some
very, very small value and the ionosphere is grown out of the chemistry
with the neutral atmosphere.

The variables TempMin, etc., are only read in if {\\tt UseMSIS = F}.

','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'INITIAL'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseApex','type' => 'logical','default' => 'F'}},{'content' => '

#APEX
T			UseApex

A model of the magnetic field of the Earth can also be used.  This
variable sets whether to use a realistic magnetic field (T) or a
dipole (F). In the current framework only the dipole works.

The default value is false.

','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'APEX'}}],'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'Initial Conditions'}},{'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

','type' => 't'},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'f107','type' => 'real','min' => '65','default' => ''}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'f107a','type' => 'real','min' => '65','default' => ''}},{'content' => '

#F107
150.0			f107
150.0			f107a

The $f10.7$ is a proxy for how bright the Sun is in a given set of
wavelengths.  For a low value (70), the temperature of the atmosphere
will be low, and the ionospheric density will be small.  For a high
value (300), the temperature will be above 1500 K, and ionospheric
electron densities will be well above $10^{12} /m^3$.  This is used in
the routine {\\tt calc_euv.f90} to determine the solar flux at the top
of the atmosphere.

','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'F107'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'HemisphericPower','type' => 'real','min' => '0.0','default' => ''}},{'content' => '

#HPI
10.0		HPI

The hemispheric power index (HPI) describes how much power is in the
hemispherically summed precipitating electrons (in Gigawatts).  This
is a number that is typically in the 1-10 range, but during active
times can reach 300.  This is used in the {\\tt get_potential.f90}
routine to get the auroral inputs at the top of the atmosphere.

','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'HPI'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '9.3','name' => 'kp','type' => 'real','min' => '0.0','default' => ''}},{'content' => '

#KP
1.0			kp

KP is a 3 hour index that summarizes the general activity level of the
magnetosphere.  It has a range from 0-9.  Currently, KP is not used in
GITM.

','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'KP'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'bx','type' => 'real','default' => ''}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'by','type' => 'real','default' => ''}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'bz','type' => 'real','default' => ''}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'vx','type' => 'real','default' => ''}},{'content' => '

#SOLARWIND
0.0			Bx
0.0			By
-2.0			Bz
400.0			Vx

The interplanetary magnetic field (IMF) and solar wind velocity are
used by a number of empirical models of the ionospheric potential.
The IMF components typically range between $\\pm 10 nT$ at Earth.  The
fields have reached values as high as 75 $nT$.  The $B_z$ component is
typically the most geoeffective at Earth, such that negative $B_z$ can
cause large ionospheric potentials.  The velocity typically ranges
between 350-600 $km/s$, although it has been known to go upwards of
1700 $km/s$.

','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'SOLARWIND'}}],'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'Indices'}},{'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Conversely, you can input time-dependent values of the solar wind and
IMF, HPI, Kp, f10.7, etc.  There are currently three methods for
inputing these quantities.

It is quite easy to incorporate other methods.  These three methods
are located in the {\\tt srcIO} directory with appropriate file names.
You can simply copy one of the files, rename the subroutine, modify it
to read in the appropriate data, add it to the {\\tt Makefile}, add a
flag in the {\\tt set_inputs.f90} (in both the {\\tt src} and {\\tt
srcIO} directories), then compile it and debug it.

','type' => 't'},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'cAMIEFileNorth','type' => 'string','default' => 'none','length' => '80'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'cAMIEFileSouth','type' => 'string','default' => 'none','length' => '80'}},{'content' => '

#AMIEFILES
b19980504n
b19980504s

Instead of using empirical models of the ionospheric potential and
auroral precipitation, you can use model results from the assimilative
mapping of ionospheric electrodynamics (AMIE) technique.  If you do,
you have to specify a Northern hemisphere and Southern hemisphere
file.

','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'AMIEFILES'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'cFileName','type' => 'string','default' => '','length' => '80'}},{'content' => '

#MHD_INDICES
imf.dat

The first method only inputs the solar wind velocity, density,
temperature and IMF.


','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'MHD_INDICES'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'cFileName','type' => 'string','length' => '80'}},{'content' => '

#NGDC_INDICES
spidr.dat

The second method takes data from the NOAA SPIDR interface.  You can
download almost all of the parameters in this format.

','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'NGDC_INDICES'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'cFileName','type' => 'string','default' => '','length' => '80'}},{'content' => '

#NOAAHPI_INDICES
power_1998.txt

The third method only accepts HPI data from the NOAA satellites.

','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'NOAAHPI_INDICES'}}],'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'Index Files'}},{'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

','type' => 't'},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'iDebugLevel','type' => 'integer','default' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'iDebugProc','type' => 'integer','default' => '0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DtReport','type' => 'real','default' => '10.0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseBarriers','type' => 'logical','default' => 'F'}},{'content' => '

#DEBUG
1			iDebugLevel
8			iDebugProc
60.0			DtReport
F			UseBarriers

Sometimes debugging can be a real pain.  This command makes it
slightly easier by allowing you to output more stuff.  The {\\tt
iDebugLevel} variable controls the amount of information output, with
0 outputting only a time-step and a message when output files are
written, and 10 being a torrent of so much information you can\'t read
it all.  You can also choose which CPU is outputting the information -
remember that MPI counts from 0 (not from 1, as most people do).  The
{\\tt DtReport} variabile says how often the time-report is given.  The
{\\tt UseBarriers} variabile is supposed to stop the code fairly often
to make sure all of the processors are on the same page, but there is
a bug in this is the {\\tt \\#SATELLITES} are used (don\'t ask).

','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'DEBUG'}}],'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'Information'}},{'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

I have tried to make the code quite versatile.  This means that most
fundemental parameters have flags so you can turn them off and on.
Most of the time, these should be left on, since all of the being {\\tt
T} means that you are running the ``physically correct\'\' condition.
But, if you want to turn something off to experiment with the physics,
you can do this.

Some of these are not really physically consistent yet.  For example,
the variable {\\tt UseDiffusion} turns off both the Eddy and Molecular
diffusion in the neutral density calculations, but leaves the Eddy
diffusion on in the temperature equation.  Also, if you turn off {\\tt
UseConduction}, the Eddy diffusion in the temperature equation is
turned off.  So, things need to be fixed a little bit.  Most of the
options really only turn off one thing, though.

This is for the neutral temperature equations and {\\bf not} for the
electron and ion equations.

','type' => 't'},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseSolarHeating','type' => 'logical','default' => 'T'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseJouleHeating','type' => 'logical','default' => 'T'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseAuroralHeating','type' => 'logical','default' => 'T'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseNOCooling  ','type' => 'logical','default' => 'T'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseOCooling   ','type' => 'logical','default' => 'T'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseConduction ','type' => 'logical','default' => 'T'}},{'content' => '

#THERMO
T		 	UseSolarHeating
T			UseJouleHeating
T		 	UseAuroralHeating
T		 	UseNOCooling
T		 	UseOCooling
T		 	UseConduction

','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'THERMO'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseDiffusion','type' => 'logical','default' => ''}},{'content' => '

#DIFFUSION
T

This only applies to the neutral densities, and includes both Eddy and
Molecular diffusion.  It should be seperated shortly.

','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'DIFFUSION'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UsePressureGradient','type' => 'logical','default' => 'T'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseIonDrag      ','type' => 'logical','default' => 'T'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseViscosity    ','type' => 'logical','default' => 'T'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseCoriolis     ','type' => 'logical','default' => 'T'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseGravity      ','type' => 'logical','default' => 'T'}},{'content' => '

#FORCING
T		UsePressureGradient
T		UseIonDrag
T		UseViscosity
T		UseCoriolis
T		UseGravity

The {\\tt UsePressureGradient} variable is ignored in this version of
GITM, since pressure solved for self-consistently within the solver.
Everything else works as a source term (in {\\tt calc_sources.f90},
except if {\\tt UseGravity = F}, gravity is zeroed in {\\tt
initialize.f90}.

','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'FORCING'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseExB             ','type' => 'logical','default' => 'T'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseIonPressureGradient','type' => 'logical','default' => 'T'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseIonGravity      ','type' => 'logical','default' => 'T'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseNeutralDrag     ','type' => 'logical','default' => 'T'}},{'content' => '

#IONFORCING
T		UseExB
T		UseIonPressure
T		UseGravity
T		UseNeutralDrag



All of these variables are used within {\\tt calc_ion_v.f90}.

','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'IONFORCING'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseIonChemistry','type' => 'logical','default' => ''}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseIonAdvection','type' => 'logical','default' => ''}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseNeutralChemistry','type' => 'logical','default' => ''}},{'content' => '

#CHEMISTRY
T		UseIonChemistry
T		UseIonAdvection	
T		UseNeutralChemistry

You can turn off the chemisty and the ion advection with these terms.


','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'CHEMISTRY'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DtPotential','type' => 'real, seconds','default' => '60.0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DtAurora','type' => 'real, seconds','default' => '60.0'}},{'content' => '

#ELECTRODYNAMICS
60.0
60.0

The electric potential and aurora are two of the most expensive
routines to run.  In addition, they typically don\'t change on a
global-scale on more than a 1-minute cadence.  So, you can set these
values to something on the order of 60 seconds.  If you are using
higher temporal resolution IMF parameters, you can set them as low as
you want.
  
','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'ELECTRODYNAMICS'}}],'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'The Control of Nature'}},{'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

','type' => 't'},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nBlocksLon','type' => 'integer','default' => '1'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nBlocksLat','type' => 'integer','default' => '1'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '90.0','name' => 'LatStart','type' => 'real','min' => '-91.0','default' => '-90.0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '91.0','name' => 'LatEnd  ','type' => 'real','min' => '-90.0','default' => '90.0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '360.0','name' => 'LonStart','type' => 'real','min' => '-180.0','default' => '180.0'}},{'content' => '

#GRID
8			nBlocksLon
4			nBlocksLat
-90.0			LatStart
90.0			LatEnd
180.0			LonStart

If LatStart and LatEnd are set to less than -90 and
greater than 90, respectively, then GITM does a whole
sphere.  If not, it models between the two.
If you want to do 1-D, set nLons=1, nLats=1 in
ModSize.f90, then recompile, then set LatStart
and LonStart to the point on the Globe you want
to model.

','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'GRID'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '90.0','name' => 'ConcentrationLatitude','type' => 'real','min' => '-90.0','default' => '65.0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '1.0','name' => 'StretchingPercentage','type' => 'real','min' => '0.0','default' => '0.0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'max' => '10.0','name' => 'StretchingFactor','type' => 'real','min' => '0.0','default' => ''}},{'content' => '

#STRETCH
65.0  			ConcentrationLatitude
0.0   			StretchingPercentage
1.0  			StretchingFactor

The stretched grid is concentrated around the ConcentrationLatitude.
The stretching is controlled by StretchingPercentage: 0 means no
stretching, 1.0 means a lot. The StretchingFactor provides further control:
greater than 1 means stretch less, less than 1 means stretch more.

','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'STRETCH'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'AltMin','type' => 'real','min' => '0.0','default' => '95.0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'AltMax','type' => 'real','min' => '$AltMin','default' => '600.0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'UseStretchedAltitude','type' => 'logical','default' => 'T'}},{'content' => '

#ALTITUDE
95.0			AltMin (km)
600.0			AltMax (km)
T			Stretched grid in altitude

','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'ALTITUDE'}}],'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'Controlling the Grid'}},{'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

','type' => 't'},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DtRestart','type' => 'real','default' => '3600.0'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nOutputTypes','type' => 'integer','default' => '1'}},{'content' => [{'content' => [{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => '3DALL','default' => 'T'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => '1DALL'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => '2DGEO'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => '2DMAG'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => '3DNEUTRAL'}},{'content' => [],'name' => 'option','type' => 'e','attrib' => {'name' => '3DION'}}],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'Outputtype','type' => 'string','input' => 'select'}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DtPlot','type' => 'real','default' => '3600.0'}}],'name' => 'for','type' => 'e','attrib' => {'to' => '$nOutputTypes','from' => '1'}},{'content' => '

#SAVEPLOTS
1800.0
1
3DALL
900.0

','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'SAVEPLOTS'}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'nSats ','type' => 'integer - max = 20','default' => ''}},{'content' => [{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'SatFile1','type' => 'string','default' => ''}},{'content' => [],'name' => 'parameter','type' => 'e','attrib' => {'name' => 'DtPlot1','type' => 'real','default' => ''}}],'name' => 'for','type' => 'e','attrib' => {'to' => '$nSats','from' => '1'}},{'content' => '

#SATELLITES
2
guvi.2002041623.in
15.0
stfd.fpi.in
60.0

','type' => 't'}],'name' => 'command','type' => 'e','attrib' => {'name' => 'SATELLITES'}}],'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'Output'}},{'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

','type' => 't'},{'content' => [],'name' => 'command','type' => 'e','attrib' => {'name' => 'END','if' => '$_IsStandAlone'}}],'name' => 'commandgroup','type' => 'e','attrib' => {'name' => 'This is the end'}},{'content' => [{'content' => '
        Directory UA/RestartOUT should exist!
','type' => 't'}],'name' => 'rule','type' => 'e','attrib' => {'expr' => '-d \'UA/RestartOUT\' or not $_IsFirstSession'}},{'content' => [{'content' => '
        Output directory UA/data should exist!
','type' => 't'}],'name' => 'rule','type' => 'e','attrib' => {'expr' => '-d \'UA/data\' or not $_IsFirstSession'}}],'name' => 'commandList','type' => 'e','attrib' => {'name' => 'Global Ionosphere Thermosphere Model: UA Component'}}];