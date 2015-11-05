INSTALLATION
============

Julia:
  1) install julia version >= 0.4

  2) in your home directory put the following line into
     .juliarc.jl  (create the file if it does not exist)

     push!(LOAD_PATH, pwd())

     This includes the current directory in which julia is run
     to the LOAD_PATH. (The list of places where Julia looks for modules to
     load)

  3) a)
     OSX users: Because OSX is stupid, create the following symlink in order
     to be able to launch julia properly from the command line

     sudo ln -s /Applications/Julia-x.x.x.app/Contents/Resources/julia/bin/julia /usr/bin/julia

     where x.x.x is your version of the Julia install.

     Linux users: You are much cooler than OSX users. No such step necessary
     for you.

     b)
     start julia from command line by typing:
     julia

     then execute the following commands one by one from within the Julia
     REPL to install necessary packages

     Pkg.add("DataFrames")
     Pkg.add("HDF5")
     Pkg.add("JLD")
     Pkg.add("PyPlot")
     Pkg.update()


Spice:
  1) download and unzip the cspice library. (tested for version N0065)
  2) done

--------------------------------------------------------------------------------

CONFIGURATION
=============
Before you can run the line of sight integration tool you need to go through
a couple of configuration steps. All configurations can be done through the
Config.jl script which is found at:

AMPS/utility/LOS/Config.jl

the general use is to call

julia Config.jl --option

The settings of the configuration process are saved to
AMPS/utility/LOS/.userSettings.conf

Running:
julia Config.jl --auto
will take you through the necessary configuration steps one by one.
The other options below can also be used after or before the --auto run in
order to overwrite/change settings done in --auto.


the following options are available ( * ) are mandatory setups)

--tmpdir <path to temporary directory> ( * )
  This directory is used to store files necessary for the LOS tool. It can be
  picked freely. In this directory the Config.jl file will create two subdirs
  'lib' and 'input'
  The 'lib' directory will contain shared libraries for spice and user defined
  c functions. (see below)
  The 'input' directory will be populated with the AMPS data files and the
  triangulated surface mesh files.
  It is not necessary to manually move any files into this directories, they
  will be updated on later steps by the Config.jl file.


--spicelib <path to cspice/lib/> ( * )
  Full path to the spice directory which contains the files
  cspice.a and csupport.a.
  Those files will then be copied into the 'lib' folder and compiled into a
  shared library (spice.so on linux, spice.dylib on OSX)


--kernelfile <full path to spice metafile>
  Full path and file name to a spice metafile that contains the full list
  of spice kernels to be loaded for the calculations
  --> the spice routine will call furnsh(metafile) on this file.


--clib <full path to custom c function definition>
  custom user file containing a function definition according to:

  void
  ColumnIntegrationFactor (double minSize,
                           double maxSize,
                           double r,
                           double * result){
  result[0]=666.0;

  // minSize = minimum size of dust particle
  // maxSize = maximum size of dust particle
  // r       = distance from observer along LOS
  // result  = provide a custom factor to be multiplied with the column density.
  // e.g. used for brightness calculation of dust grains
  }


--noclib
  do not use any shared c library
  If you never defined --clib this is unnecessary. It is only used to remove
  a c library that has been used in previous runs.


--meshfile        .ply file of the body surface mesh")
--meshfileshadow  .ply file of the body surface mesh for shadow calc.")
--docheckshadow   yes or no if shadow calculation is needed")
--datafile        .full path to h5 AMPS output file")
--clean           remove 'lib' and 'input' dirs in tmpfile")
--help            show this message"


RUN
===
Before you can use the LOS tool you need to modify the AMPS data with the
prepareAmpsData.jl script. This will create a new HDF5 (.h5) file with only
a subset of data from the original file.


7) Prepare AMPS data by executing:

julia prepareAmpsData.jl FileName

this will create a .h5 file which is later used
to load the data
