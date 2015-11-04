INSTALLATION
------------

Julia:
  1) install julia version >= 0.4

  2) in your home directory put the following into
     .juliarc.jl  (create the file if necessary)
     push!(LOAD_PATH, pwd())

  3) start julia from command line by typing
     julia

     then execute the following commands one by one from within
     julia to install necessary packages

     Pkg.add("DataFrames")
     Pkg.add("HDF5")
     Pkg.add("JLD")
     Pkg.add("PyPlot")
     Pkg.update()

SPICE:
  
  1) download and unzip the cspice library






7) Prepare AMPS data by executing:

julia prepareAmpsData.jl FileName

this will create a .h5 file which is later used
to load the data




:
