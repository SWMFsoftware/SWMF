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


4) download and unzip the cspice library

5) go to cspice/lib/ and compile a shared library:

6) execut the following commands

ar -x cspice.a
ar -x csupport.a
gcc -shared -fPIC -lm *.o -o spice.so

7) Prepare AMPS data by executing:

julia prepareAmpsData.jl FileName

this will create a .h5 file which is later used
to load the data




:
