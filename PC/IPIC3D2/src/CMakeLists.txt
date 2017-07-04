cmake_minimum_required(VERSION 2.8.8) 
# compiler set in ../cmake/cmake_template.cmake.XeonPhi
#message ("for Xeon Phi:")
#message ("cmake .. -DCMAKE_TOOLCHAIN_FILE=../cmake/cmake_template.cmake.XeonPhi")
#message ("for Xeon:")
#message ("cmake .. -DCMAKE_TOOLCHAIN_FILE=../cmake/cmake_template.cmake.Xeon")
#
# Project declaration
#

project(iPic3D)

#
# Set compiler flags per system
#
if (${CMAKE_SYSTEM_PROCESSOR} STREQUAL "k1om") ## Xeon Phi
   option(IPIC_XEONPHI "ipic xeon phi standard compile flags" on)
   if(IPIC_XEONPHI) 
     set(CMAKE_CXX_FLAGS "-O3 -openmp -fno-exceptions -fp-model fast=2 -vec-report -mmic")
   else()
     set(CMAKE_CXX_FLAGS "-mmic")
   endif()
   #set(CMAKE_CXX_FLAGS "$(CMAKE_CXX_FLAGS) -mmic")
   #set(CMAKE_CXX_FLAGS "-O3 -openmp -fno-exceptions -vec-report -mmic")
   #set(CMAKE_CXX_FLAGS "-openmp -g -mmic") # set flags for Xeon Phi, totalview
elseif (${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64") ## Xeon
   option(IPIC_XEON "optimisation with icpc"  OFF)
   if(IPIC_XEON) 
     #reporting: -g -vec-report
     #optimization: -O3 -xHost -fno-exceptions
     set(CMAKE_CXX_FLAGS "-openmp -fno-exceptions -O3 -xHost -vec-report")
     set(CMAKE_CXX_COMPILER "icpc")  
   endif()
   #set(CMAKE_CXX_FLAGS "-O3 -openmp -g -xHost -fno-exceptions -vec-report")
   set(CMAKE_CXX_FLAGS  "-O3 ") # -DPRINTPCL for printing particle subcycles  -DNO_HDF5 to disable output
else()
   set(CMAKE_CXX_FLAGS "-O3")
endif()

#
# Find third class libraries
#
if (${CMAKE_SYSTEM_PROCESSOR} STREQUAL "k1om") ## Xeon Phi
   set(CMAKE_PREFIX_PATH /opt/hdf5/1.8.10-patch1-mic)			
   set(VARIOUS_LIB /opt/intel/lib/mic)
   set (EXTRA_LIBS ${VARIOUS_LIB}/libimf.so ${VARIOUS_LIB}/libsvml.so ${VARIOUS_LIB}/libirng.so ${VARIOUS_LIB}/libintlc.so ${MPELIB})
else()
   find_package(MPI)
   set(EXTRA_LIBS "")   
endif()

SITE_NAME(myhostname)
if (${myhostname} MATCHES "(.*)pdc(.*)")
   set(HDF5_INCLUDE_DIRS /opt/cray/hdf5/1.8.7/cray/73/include)
elseif(${myhostname} MATCHES "localhost")
   set(CMAKE_CXX_COMPILER mpicc)  
   set(CMAKE_CXX_FLAGS "-O2")
   set(CMAKE_PREFIX_PATH /Users/Ivy/Software/install/hdf5-1.8.14/)
   find_package(HDF5 COMPONENTS HL C )
else()
   find_package(HDF5 COMPONENTS HL C )   
endif()


#
# include and lib directories
#

# include_directories: files there are accessible to the project
include_directories(
        include
        ${HDF5_INCLUDE_DIRS}
        ${HDF5_HL_INCLUDE_DIRS}
        ${MPI_INCLUDE_PATH}
)

link_directories(
        ${HDF5_LIBRARY_DIRS}
	${HDF5_hdf5_hl_LIBRARY_DIRS}
        ${MPI_LIBRARY_DIRS}
	${EXTRA_LIBS}
)

#
# Header file list
#

file(
        GLOB
        inc_files
        include/*.h
)

#
# Source file list
#

file(
        GLOB
        src_files
        ConfigFile/src/*.cpp
        PSKOutput3D/*.cpp
        bc/*.cpp
        communication/*.cpp
        fields/*.cpp
        grids/*.cpp
        inputoutput/*.cpp
        mathlib/*.cpp
        mpidata/*.cpp
        particles/*.cpp
        performances/*.cpp
        processtopology/*.cpp
        solvers/*.cpp
        utility/*.cpp
        main/*.cpp
)

#
# Macro definitions
#

set(TEST_B $ENV{BATSRUS})
if(DEFINED TEST_B)
  add_definitions( -DBATSRUS )
  message(" WARNING: BATSRUS flag is active.")
else(DEFINED TEST_B)
  message(" INFO: BATSRUS is not active.")
endif(DEFINED TEST_B)

#
# Executable declaration
#

# Particle solver
add_executable(
        iPic3D
        iPic3D.cpp
)

option(STATIC_LINK "Choose static link"  on)
if(STATIC_LINK) 
	add_library(
	    iPic3Dlib      #name of the library
		STATIC	   #type of the library
	    ${inc_files}    # stuff to build the library
	    ${src_files}
	)
else()
	add_library(
	    iPic3Dlib      #name of the library
		SHARED	   #type of the library
	    ${inc_files}    # stuff to build the library
	    ${src_files}
	)
endif()	


#
# Link external libraries
#
target_link_libraries(
         iPic3Dlib
         ${HDF5_LIBRARIES}
         ${HDF5_HL_LIBRARIES}
         ${MPI_LIBRARIES}
	 ${EXTRA_LIBS}
)

target_link_libraries(
         iPic3D
         iPic3Dlib
         #pmi
)

## to save the executable in the folder where the CMakeLists.txt file is, i.e. CMAKE_CURRENT_SOURCE_DIR
set_target_properties(iPic3D PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

## debug releases have a _d appended to the executable
set_target_properties(iPic3D PROPERTIES DEBUG_POSTFIX "_d")


message("Which system am I compiling for:")
message("MYHOSTNAME is ${myhostname}")
message("CMAKE_SYSTEM_PROCESSOR is ${CMAKE_SYSTEM_PROCESSOR}")
message("Compiler & compiler flags:")
message("CMAKE_CXX_COMPILER is ${CMAKE_CXX_COMPILER}")
message("CMAKE_CXX_FLAGS is ${CMAKE_CXX_FLAGS}")
message("HDF5_INCLUDE_DIRS is ${HDF5_INCLUDE_DIRS}")
message("HDF5_LIBRARIES"   is ${HDF5_LIBRARIES})
message("HDF5_HL_LIBRARIES is ${HDF5_HL_LIBRARIES}")
message("MPI_LIBRARIES is ${MPI_LIBRARIES}")
