To install and run iPIC3D, you need: 
   cmake (at least version 2.8), MPI and HDF5 (optional).

# create a build directory
  mkdir build && cd build
  cmake /path/to/root/where/CMakeLists.txt/located

# compile, if successful, you will find an executable named iPIC3D in build directory
  make

# run a test case: copy an inputfile named as testXXX.inp from /inputfiles to build directory
# make sure you create an folder for output as specified in the input file
# make sure no_of_proc = XLEN x YLEN x ZLEN as specified in the input file
  mpiexec -n no_of_proc ./iPIC3D  inputfilename.inp
