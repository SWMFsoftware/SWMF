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

','type' => 't'},{'content' => [{'content' => [],'attrib' => {'default' => '2000','type' => 'integer','name' => 'iYear'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '12','default' => '3','min' => '1','type' => 'integer','name' => 'iMonth'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '31','default' => '21','min' => '1','type' => 'integer','name' => 'iDay'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '23','default' => '0','min' => '0','type' => 'integer','name' => 'iHour'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '59','default' => '0','min' => '0','type' => 'integer','name' => 'iMinute'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '59','default' => '0','min' => '0','type' => 'integer','name' => 'iSecond'},'name' => 'parameter','type' => 'e'},{'content' => '

#TIMESTART
1999			iYear
03			iMonth
18			iDay
00			iHour
00			iMinute
00			iSecond

This command is only used in the stand alone mode.

The #STARTTIME command sets the initial date and time for the simulation.
','type' => 't'}],'attrib' => {'if' => '$_IsStandAlone','alias' => 'STARTTIME','name' => 'TIMESTART'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'default' => '2000','type' => 'integer','name' => 'iYear'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '12','default' => '3','min' => '1','type' => 'integer','name' => 'iMonth'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '31','default' => '21','min' => '1','type' => 'integer','name' => 'iDay'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '23','default' => '0','min' => '0','type' => 'integer','name' => 'iHour'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '59','default' => '0','min' => '0','type' => 'integer','name' => 'iMinute'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '59','default' => '0','min' => '0','type' => 'integer','name' => 'iSecond'},'name' => 'parameter','type' => 'e'},{'content' => '

#TIMEEND
1999			iYear
03			iMonth
25			iDay
00			iHour
00			iMinute
00			iDay

This command is only used in the stand alone mode.

The #TIMEEND command sets the final date and time for the simulation.
','type' => 't'}],'attrib' => {'if' => '$_IsStandAlone','alias' => 'ENDTIME','name' => 'TIMEEND'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'default' => 'F','type' => 'logical','name' => 'DoRestart'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => '
                File UA/restartIN/header.rst should exist!
	','type' => 't'}],'attrib' => {'expr' => '-f \'UA/restartIN/header.rst\' or not $DoRestart'},'name' => 'rule','type' => 'e'},{'content' => '

#RESTART
F			DoRestart

There are two commands that are typically not input by a user, but are
specified in the restart header file that is read in the exact same
way as the input file.  It is possible to set these variables, though.

','type' => 't'}],'attrib' => {'name' => 'RESTART'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'default' => '','type' => 'integer','name' => 'iStep'},'name' => 'parameter','type' => 'e'},{'content' => '

#ISTEP
1

This sets the current iStep to the read in value instead of starting at
iteration 1.

','type' => 't'}],'attrib' => {'name' => 'ISTEP'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'default' => '','type' => 'real','name' => 'tsimulation'},'name' => 'parameter','type' => 'e'},{'content' => '

#TSIMULATION
0.0

This offsets the current start time by the read in amount.  It is simply
added to the #STARTTIME input time.

','type' => 't'}],'attrib' => {'name' => 'TSIMULATION'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'default' => '','type' => 'real, seconds','name' => 'CPUTimeMax'},'name' => 'parameter','type' => 'e'},{'content' => '

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

','type' => 't'}],'attrib' => {'name' => 'CPUTIMEMAX'},'name' => 'command','type' => 'e'}],'attrib' => {'name' => 'Time Variables'},'name' => 'commandgroup','type' => 'e'},{'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

','type' => 't'},{'content' => [{'content' => [],'attrib' => {'default' => '','type' => 'logical','name' => 'UseMSIS'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'default' => '','type' => 'logical','name' => 'UseIRI'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'default' => '','type' => 'real','name' => 'TempMin'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'default' => '','type' => 'real','name' => 'TempMax'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'default' => '','min' => '0.01','type' => 'real','name' => 'Number Density'},'name' => 'parameter','type' => 'e'}],'attrib' => {'from' => '1','to' => '3'},'name' => 'for','type' => 'e'}],'attrib' => {'expr' => '$UseMSIS'},'name' => 'if','type' => 'e'},{'content' => '

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

','type' => 't'}],'attrib' => {'name' => 'INITIAL'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'default' => 'F','type' => 'logical','name' => 'UseApex'},'name' => 'parameter','type' => 'e'},{'content' => '

#APEX
T			UseApex

A model of the magnetic field of the Earth can also be used.  This
variable sets whether to use a realistic magnetic field (T) or a
dipole (F). In the current framework only the dipole works.

The default value is false.

','type' => 't'}],'attrib' => {'name' => 'APEX'},'name' => 'command','type' => 'e'}],'attrib' => {'name' => 'Initial Conditions'},'name' => 'commandgroup','type' => 'e'},{'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

','type' => 't'},{'content' => [{'content' => [],'attrib' => {'default' => '','min' => '65','type' => 'real','name' => 'f107'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'default' => '','min' => '65','type' => 'real','name' => 'f107a'},'name' => 'parameter','type' => 'e'},{'content' => '

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

','type' => 't'}],'attrib' => {'name' => 'F107'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'default' => '','min' => '0.0','type' => 'real','name' => 'HemisphericPower'},'name' => 'parameter','type' => 'e'},{'content' => '

#HPI
10.0		HPI

The hemispheric power index (HPI) describes how much power is in the
hemispherically summed precipitating electrons (in Gigawatts).  This
is a number that is typically in the 1-10 range, but during active
times can reach 300.  This is used in the {\\tt get_potential.f90}
routine to get the auroral inputs at the top of the atmosphere.

','type' => 't'}],'attrib' => {'name' => 'HPI'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'max' => '9.3','default' => '','min' => '0.0','type' => 'real','name' => 'kp'},'name' => 'parameter','type' => 'e'},{'content' => '

#KP
1.0			kp

KP is a 3 hour index that summarizes the general activity level of the
magnetosphere.  It has a range from 0-9.  Currently, KP is not used in
GITM.

','type' => 't'}],'attrib' => {'name' => 'KP'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'default' => '','type' => 'real','name' => 'bx'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'default' => '','type' => 'real','name' => 'by'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'default' => '','type' => 'real','name' => 'bz'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'default' => '','type' => 'real','name' => 'vx'},'name' => 'parameter','type' => 'e'},{'content' => '

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

','type' => 't'}],'attrib' => {'name' => 'SOLARWIND'},'name' => 'command','type' => 'e'}],'attrib' => {'name' => 'Indices'},'name' => 'commandgroup','type' => 'e'},{'content' => [{'content' => '
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

','type' => 't'},{'content' => [{'content' => [],'attrib' => {'default' => 'none','length' => '80','type' => 'string','name' => 'cAMIEFileNorth'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'default' => 'none','length' => '80','type' => 'string','name' => 'cAMIEFileSouth'},'name' => 'parameter','type' => 'e'},{'content' => '

#AMIEFILES
b19980504n
b19980504s

Instead of using empirical models of the ionospheric potential and
auroral precipitation, you can use model results from the assimilative
mapping of ionospheric electrodynamics (AMIE) technique.  If you do,
you have to specify a Northern hemisphere and Southern hemisphere
file.

','type' => 't'}],'attrib' => {'name' => 'AMIEFILES'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'default' => '','length' => '80','type' => 'string','name' => 'cFileName'},'name' => 'parameter','type' => 'e'},{'content' => '

#MHD_INDICES
imf.dat

The first method only inputs the solar wind velocity, density,
temperature and IMF.


','type' => 't'}],'attrib' => {'name' => 'MHD_INDICES'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'length' => '80','type' => 'string','name' => 'cFileName'},'name' => 'parameter','type' => 'e'},{'content' => '

#NGDC_INDICES
spidr.dat

The second method takes data from the NOAA SPIDR interface.  You can
download almost all of the parameters in this format.

','type' => 't'}],'attrib' => {'name' => 'NGDC_INDICES'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'default' => '','length' => '80','type' => 'string','name' => 'cFileName'},'name' => 'parameter','type' => 'e'},{'content' => '

#NOAAHPI_INDICES
power_1998.txt

The third method only accepts HPI data from the NOAA satellites.

','type' => 't'}],'attrib' => {'name' => 'NOAAHPI_INDICES'},'name' => 'command','type' => 'e'}],'attrib' => {'name' => 'Index Files'},'name' => 'commandgroup','type' => 'e'},{'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

','type' => 't'},{'content' => [{'content' => [],'attrib' => {'default' => '0','type' => 'integer','name' => 'iDebugLevel'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'default' => '0','type' => 'integer','name' => 'iDebugProc'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'default' => '10.0','type' => 'real','name' => 'DtReport'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'default' => 'F','type' => 'logical','name' => 'UseBarriers'},'name' => 'parameter','type' => 'e'},{'content' => '

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
{\\tt DtReport} variable says how often the time-report is given.  The
{\\tt UseBarriers} variable is supposed to stop the code fairly often
to make sure all of the processors are on the same page, but there is
a bug in this is the {\\tt \\#SATELLITES} are used (don\'t ask).

','type' => 't'}],'attrib' => {'name' => 'DEBUG'},'name' => 'command','type' => 'e'}],'attrib' => {'name' => 'Information'},'name' => 'commandgroup','type' => 'e'},{'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

The GITM code development has been aimed toward making the code quite versatile.  
This means that most
fundamental parameters have flags so you can turn them off and on.
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

','type' => 't'},{'content' => [{'content' => [],'attrib' => {'default' => 'T','type' => 'logical','name' => 'UseSolarHeating'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'default' => 'T','type' => 'logical','name' => 'UseJouleHeating'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'default' => 'T','type' => 'logical','name' => 'UseAuroralHeating'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'default' => 'T','type' => 'logical','name' => 'UseNOCooling  '},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'default' => 'T','type' => 'logical','name' => 'UseOCooling   '},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'default' => 'T','type' => 'logical','name' => 'UseConduction '},'name' => 'parameter','type' => 'e'},{'content' => '

#THERMO
T		 	UseSolarHeating
T			UseJouleHeating
T		 	UseAuroralHeating
T		 	UseNOCooling
T		 	UseOCooling
T		 	UseConduction

','type' => 't'}],'attrib' => {'name' => 'THERMO'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'default' => '','type' => 'logical','name' => 'UseDiffusion'},'name' => 'parameter','type' => 'e'},{'content' => '

#DIFFUSION
T

This only applies to the neutral densities, and includes both Eddy and
Molecular diffusion.  It should be separated shortly.

','type' => 't'}],'attrib' => {'name' => 'DIFFUSION'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'default' => 'T','type' => 'logical','name' => 'UsePressureGradient'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'default' => 'T','type' => 'logical','name' => 'UseIonDrag      '},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'default' => 'T','type' => 'logical','name' => 'UseViscosity    '},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'default' => 'T','type' => 'logical','name' => 'UseCoriolis     '},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'default' => 'T','type' => 'logical','name' => 'UseGravity      '},'name' => 'parameter','type' => 'e'},{'content' => '

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

','type' => 't'}],'attrib' => {'name' => 'FORCING'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'default' => 'T','type' => 'logical','name' => 'UseExB             '},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'default' => 'T','type' => 'logical','name' => 'UseIonPressureGradient'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'default' => 'T','type' => 'logical','name' => 'UseIonGravity      '},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'default' => 'T','type' => 'logical','name' => 'UseNeutralDrag     '},'name' => 'parameter','type' => 'e'},{'content' => '

#IONFORCING
T		UseExB
T		UseIonPressure
T		UseGravity
T		UseNeutralDrag



All of these variables are used within {\\tt calc_ion_v.f90}.

','type' => 't'}],'attrib' => {'name' => 'IONFORCING'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'default' => '','type' => 'logical','name' => 'UseIonChemistry'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'default' => '','type' => 'logical','name' => 'UseIonAdvection'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'default' => '','type' => 'logical','name' => 'UseNeutralChemistry'},'name' => 'parameter','type' => 'e'},{'content' => '

#CHEMISTRY
T		UseIonChemistry
T		UseIonAdvection	
T		UseNeutralChemistry

You can turn off the chemistry and the ion advection with these terms.


','type' => 't'}],'attrib' => {'name' => 'CHEMISTRY'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'default' => '60.0','type' => 'real, seconds','name' => 'DtPotential'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'default' => '60.0','type' => 'real, seconds','name' => 'DtAurora'},'name' => 'parameter','type' => 'e'},{'content' => '

#ELECTRODYNAMICS
60.0
60.0

The electric potential and aurora are two of the most expensive
routines to run.  In addition, they typically don\'t change on a
global-scale on more than a 1-minute cadence.  So, you can set these
values to something on the order of 60 seconds.  If you are using
higher temporal resolution IMF parameters, you can set them as low as
you want.
  
','type' => 't'}],'attrib' => {'name' => 'ELECTRODYNAMICS'},'name' => 'command','type' => 'e'}],'attrib' => {'name' => 'The Control of Nature'},'name' => 'commandgroup','type' => 'e'},{'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

','type' => 't'},{'content' => [{'content' => [],'attrib' => {'default' => '1','type' => 'integer','name' => 'nBlocksLon'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'default' => '1','type' => 'integer','name' => 'nBlocksLat'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '90.0','default' => '-90.0','min' => '-91.0','type' => 'real','name' => 'LatStart'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '91.0','default' => '90.0','min' => '-90.0','type' => 'real','name' => 'LatEnd  '},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '360.0','default' => '180.0','min' => '-180.0','type' => 'real','name' => 'LonStart'},'name' => 'parameter','type' => 'e'},{'content' => '

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

','type' => 't'}],'attrib' => {'name' => 'GRID'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'max' => '90.0','default' => '65.0','min' => '-90.0','type' => 'real','name' => 'ConcentrationLatitude'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '1.0','default' => '0.0','min' => '0.0','type' => 'real','name' => 'StretchingPercentage'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'max' => '10.0','default' => '','min' => '0.0','type' => 'real','name' => 'StretchingFactor'},'name' => 'parameter','type' => 'e'},{'content' => '

#STRETCH
65.0  			ConcentrationLatitude
0.0   			StretchingPercentage
1.0  			StretchingFactor

The stretched grid is concentrated around the ConcentrationLatitude.
The stretching is controlled by StretchingPercentage: 0 means no
stretching, 1.0 means a lot. The StretchingFactor provides further control:
greater than 1 means stretch less, less than 1 means stretch more.

','type' => 't'}],'attrib' => {'name' => 'STRETCH'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'default' => '95.0','min' => '0.0','type' => 'real','name' => 'AltMin'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'default' => '600.0','min' => '$AltMin','type' => 'real','name' => 'AltMax'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'default' => 'T','type' => 'logical','name' => 'UseStretchedAltitude'},'name' => 'parameter','type' => 'e'},{'content' => '

#ALTITUDE
95.0			AltMin (km)
600.0			AltMax (km)
T			Stretched grid in altitude

','type' => 't'}],'attrib' => {'name' => 'ALTITUDE'},'name' => 'command','type' => 'e'}],'attrib' => {'name' => 'Controlling the Grid'},'name' => 'commandgroup','type' => 'e'},{'content' => [{'content' => '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

','type' => 't'},{'content' => [{'content' => [],'attrib' => {'default' => '3600.0','type' => 'real','name' => 'DtRestart'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'default' => '1','type' => 'integer','name' => 'nOutputTypes'},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [{'content' => [],'attrib' => {'default' => 'T','name' => '3DALL'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => '1DALL'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => '2DGEO'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => '2DMAG'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => '3DNEUTRAL'},'name' => 'option','type' => 'e'},{'content' => [],'attrib' => {'name' => '3DION'},'name' => 'option','type' => 'e'}],'attrib' => {'type' => 'string','name' => 'Outputtype','input' => 'select'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'default' => '3600.0','type' => 'real','name' => 'DtPlot'},'name' => 'parameter','type' => 'e'}],'attrib' => {'from' => '1','to' => '$nOutputTypes'},'name' => 'for','type' => 'e'},{'content' => '

#SAVEPLOTS
1800.0
1
3DALL
900.0

','type' => 't'}],'attrib' => {'name' => 'SAVEPLOTS'},'name' => 'command','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'default' => '','type' => 'integer - max = 20','name' => 'nSats '},'name' => 'parameter','type' => 'e'},{'content' => [{'content' => [],'attrib' => {'default' => '','type' => 'string','name' => 'SatFile1'},'name' => 'parameter','type' => 'e'},{'content' => [],'attrib' => {'default' => '','type' => 'real','name' => 'DtPlot1'},'name' => 'parameter','type' => 'e'}],'attrib' => {'from' => '1','to' => '$nSats'},'name' => 'for','type' => 'e'},{'content' => '

#SATELLITES
2
guvi.2002041623.in
15.0
stfd.fpi.in
60.0

','type' => 't'}],'attrib' => {'name' => 'SATELLITES'},'name' => 'command','type' => 'e'}],'attrib' => {'name' => 'Output'},'name' => 'commandgroup','type' => 'e'},{'content' => [{'content' => '
        Directory UA/restartOUT should exist!
','type' => 't'}],'attrib' => {'expr' => '-d \'UA/restartOUT\' or not $_IsFirstSession'},'name' => 'rule','type' => 'e'},{'content' => [{'content' => '
        Output directory UA/data should exist!
','type' => 't'}],'attrib' => {'expr' => '-d \'UA/data\' or not $_IsFirstSession'},'name' => 'rule','type' => 'e'}],'attrib' => {'name' => 'Global Ionosphere Thermosphere Model 2: UA Component'},'name' => 'commandList','type' => 'e'}];