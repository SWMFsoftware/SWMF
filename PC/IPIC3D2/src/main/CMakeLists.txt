cmake_minimum_required(VERSION 2.8.8)

#
# Project declaration
#

project(iPic3Dlib)

#
# Set exec path
#

set(LIBRARY_OUTPUT_PATH lib)
set(IPIC_HOME ${CMAKE_CURRENT_SOURCE_DIR}/..)

#
# Find third class libraries
#

option(IPIC_HDF5_OUTPUT "Write output" ON)
if(IPIC_HDF5_OUTPUT)
find_package(HDF5 COMPONENTS HL C REQUIRED)
endif()
find_package(MPI REQUIRED)

#
# Header file list
#

file(
        GLOB
        inc_files
        ${IPIC_HOME}/include/*.h
)

#
# Source file list
#

file(
        GLOB
        src_files
        ${IPIC_HOME}/ConfigFile/src/*.cpp
        ${IPIC_HOME}/bc/*.cpp
        ${IPIC_HOME}/communication/*.cpp
        ${IPIC_HOME}/fields/*.cpp
        ${IPIC_HOME}/grids/*.cpp
        ${IPIC_HOME}/mathlib/*.cpp
        ${IPIC_HOME}/mpidata/*.cpp
        ${IPIC_HOME}/particles/*.cpp
        ${IPIC_HOME}/performances/*.cpp
        ${IPIC_HOME}/processtopology/*.cpp
        ${IPIC_HOME}/solvers/*.cpp
        ${IPIC_HOME}/utility/*.cpp
        ${IPIC_HOME}/main/*.cpp
)
# currently HDF5 is the only output format
if(IPIC_HDF5_OUTPUT)
  file(
        GLOB
        hdf_files
        ${IPIC_HOME}/PSKOutput3D/*.cpp
        ${IPIC_HOME}/inputoutput/*.cpp
  )
endif()

#
# Compilation options
#

if(NOT IPIC_HDF5_OUTPUT)
  add_definitions(-DNO_HDF5)
endif()
option(IPIC_PHDF5 "Output is done using Parallel HDF5 without the H5hut library")
if(IPIC_PHDF5)
  add_definitions(-DPHDF5)
  #add_definitions(-DPARALLEL_IO)
endif()
option(IPIC_H5HUT_OUTPUT "Output is done using Parallel HDF5 with the H5hut library")
if(IPIC_H5HUT_OUTPUT)
  add_subdirectory(${IPIC_HOME}/H5hut-io H5hut-io)
  set(H5HUTIO_INCLUDE_DIRS ${IPIC_HOME}/H5hut-io/include)
  set(H5HUTIO_LIBRARIES H5hut-io)
  set(H5HUT_INCLUDE_DIRS $ENV{H5HUT_DIR}/include)

  add_definitions(-DPHDF5)
  add_definitions(-DUSEH5HUT)
  add_definitions(-DPARALLEL_IO) # needed to correct for a bug in H5hut
endif()

option(IPIC_BATSRUS "Compile library with coupling code for BATS-R-US" OFF)
if(IPIC_BATSRUS)
  add_definitions( -DBATSRUS )
endif()

option(IPIC_XEONPHI "Compile options for Xeon Phi" OFF)
if(IPIC_XEONPHI)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
endif()

#
# include and lib directories
#

include_directories(
        ${IPIC_HOME}/include
        ${H5HUT_INCLUDE_DIRS}
        ${H5HUTIO_INCLUDE_DIRS}
        ${HDF5_INCLUDE_DIRS}
        ${MPI_INCLUDE_PATH}
)

link_directories(
        ${H5HUTIO_LIBRARY_DIRS}
)

#
# Create the library
#

add_library(
        iPic3Dlib
        STATIC
        ${inc_files}
        ${src_files}
        ${hdf_files}
)

#
# Link external libraries
#

target_link_libraries(
         iPic3Dlib
         ${H5HUTIO_LIBRARIES}
         #${HDF5_LIBRARIES}
)
