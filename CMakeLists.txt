# CMakeLists.txt for libeicsmear.
#
# To build using cmake, create a build directory, navigate
# to it and run cmake. e.g. assuming we start in the eicsmear
# source directory:
#  mkdir build
#  cd build
#  cmake .. (or whatever the path to the source directory is)
#
# For PYTHIA6-related code to be built the user must
# provide the following argument to cmake:
#  -DPYTHIA6_LIBDIR=<path to directory containing libPythia6>
#
# You can specify an install directory via
#  -DCMAKE_INSTALL_PREFIX=<path>
# then do
#  make install
# to install the headers and libraries under that location.
# There is no "make uninstall" but (on Unix-like systems)
# you can do
#  xargs rm < install_manifest.txt
# from the cmake build directory.
#
# Run everything from this single file because the sources
# from each subdirectory have to be combined in a single
# shared library, rather than being one library per directory.

cmake_minimum_required(VERSION 3.1)

project(EICSMEAR)



# Add this project's module path to the default path.
# Needed to locate this project's FindROOT.cmake.
set(
   CMAKE_MODULE_PATH
   ${CMAKE_MODULE_PATH}
   ${CMAKE_SOURCE_DIR}/cmake/Modules/
)

# This module sets up ROOT information
# It defines:
# ROOT_FOUND             If the ROOT is found
# ROOT_INCLUDE_DIR       PATH to the include directory
# ROOT_INCLUDE_DIRS      PATH to the include directories (not cached)
# ROOT_LIBRARIES         Most common libraries
# ROOT_<name>_LIBRARY    Full path to the library <name>
# ROOT_LIBRARY_DIR       PATH to the library directory
# ROOT_ETC_DIR           PATH to the etc directory
# ROOT_DEFINITIONS       Compiler definitions
# ROOT_CXX_FLAGS         Compiler flags to used by client packages
# ROOT_C_FLAGS           Compiler flags to used by client packages
# ROOT_EXE_LINKER_FLAGS  Linker flags to used by client packages
#
# Updated by K. Smith (ksmith37@nd.edu) to properly handle
#  dependencies in ROOT_GENERATE_DICTIONARY

find_package(ROOT)
# ROOT6 needs C++11 standard
string(SUBSTRING ${ROOT_VERSION} 0 1 ROOT_MAJOR)

# message ( " ROOT_VERSION = " ${ROOT_VERSION})
# message ( " ROOT_MAJOR = " ${ROOT_MAJOR})
IF(${ROOT_MAJOR} GREATER 5)
  if(NOT DEFINED CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD 11)
  endif()
  
  if( ${CMAKE_CXX_STANDARD} LESS 11 )
    message ( FATAL_ERROR " ROOT_VERSION = " ${ROOT_VERSION} " needs CMAKE_CXX_STANDARD >= 11 " )
  endif()

ENDIF()

## Need to have this base for rootcint
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/include/eicsmear)

# List all main eRHIC headers EXCEPT PYTHIA-related files.
# We don't need to prefix "${CMAKE_SOURCE_DIR}/include/" before all
# the header names for the main library (cmake figures it out) BUT
# we do for it to work nicely with the ROOT dictionary creation.
# Otherwise the dictionary-generation code doesn't seem to locate
# the header files properly.
SET(
   ERHIC_HEADERS
   eicsmear/erhic/BeamParticles.h
   eicsmear/erhic/DisKinematics.h
   eicsmear/erhic/EventBase.h
   eicsmear/erhic/EventDis.h
   eicsmear/erhic/EventDjangoh.h
   eicsmear/erhic/EventDpmjet.h
   eicsmear/erhic/EventFactory.h
   eicsmear/erhic/EventGmcTrans.h
   eicsmear/erhic/EventMC.h
   eicsmear/erhic/EventMCFilterABC.h
   eicsmear/erhic/EventMilou.h
   eicsmear/erhic/EventPepsi.h
   eicsmear/erhic/EventPythia.h
   eicsmear/erhic/EventRapgap.h
   eicsmear/erhic/File.h
   eicsmear/erhic/Forester.h
   eicsmear/erhic/Kinematics.h
   eicsmear/erhic/Particle.h
   eicsmear/erhic/ParticleIdentifier.h
   eicsmear/erhic/ParticleMC.h
   eicsmear/erhic/Pid.h
   eicsmear/erhic/VirtualEvent.h
   eicsmear/erhic/VirtualParticle.h
   )

INSTALL(DIRECTORY ${CMAKE_SOURCE_DIR}/include/eicsmear/erhic DESTINATION include/eicsmear )
INSTALL(FILES ${CMAKE_SOURCE_DIR}/include/eicsmear/functions.h DESTINATION include/eicsmear )
#include_directories(${CMAKE_SOURCE_DIR}/include/eicsmear/erhic)

# Smearing-related headers.
# These are referenced by a different LinkDef to ERHIC_HEADERS so
# we list the separately.
SET(
   SMEAR_HEADERS
   eicsmear/smear/Acceptance.h
   eicsmear/smear/Bremsstrahlung.h
   eicsmear/smear/Detector.h
   eicsmear/smear/Device.h
   eicsmear/smear/Distributor.h
   eicsmear/smear/EventS.h
   eicsmear/smear/EventSmear.h
   eicsmear/smear/EventDisFactory.h
   eicsmear/smear/FormulaString.h
   eicsmear/smear/functions.h
   eicsmear/smear/ParticleID.h
   eicsmear/smear/ParticleMCS.h
   eicsmear/smear/PerfectID.h
   eicsmear/smear/PlanarTracker.h
   eicsmear/smear/RadialTracker.h
   eicsmear/smear/Smear.h
   eicsmear/smear/Smearer.h
   eicsmear/smear/Tracker.h
)
INSTALL(DIRECTORY ${CMAKE_SOURCE_DIR}/include/eicsmear/smear DESTINATION include/eicsmear )

# List all main eRHIC source files EXCEPT PYTHIA-related files.
# This includes smearing code as the separate LinkDefs have no role here.
SET(
   ERHIC_SOURCE
   src/functions.cxx
   src/erhic/BeamParticles.cxx
   src/erhic/BuildTree.cxx
   src/erhic/EventDis.cxx
   src/erhic/EventDjangoh.cxx
   src/erhic/EventDpmjet.cxx
   src/erhic/EventFactory.cxx
   src/erhic/EventGmcTrans.cxx
   src/erhic/EventMC.cxx
   src/erhic/EventMilou.cxx
   src/erhic/EventPepsi.cxx
   src/erhic/EventPythia.cxx
   src/erhic/EventRapgap.cxx
   src/erhic/File.cxx
   src/erhic/Forester.cxx
   src/erhic/Kinematics.cxx
   src/erhic/ParticleIdentifier.cxx
   src/erhic/ParticleMC.cxx
   src/erhic/Pid.cxx
   src/smear/Acceptance.cxx
   src/smear/Bremsstrahlung.cxx
   src/smear/Detector.cxx
   src/smear/Device.cxx
   src/smear/Distributor.cxx
   src/smear/EventDisFactory.cxx
   src/smear/EventSmear.cxx
   src/smear/FormulaString.cxx
   src/smear/ParticleID.cxx
   src/smear/ParticleMCS.cxx
   src/smear/PerfectID.cxx
   src/smear/PlanarTracker.cxx
   src/smear/RadialTracker.cxx
   src/smear/Smear.cxx
   src/smear/SmearTree.cxx
   src/smear/Tracker.cxx
)

# Generate various CINT dictionaries.
# Arguments are:
#  Input header files
#  LinkDef file
#  Output dictionary file
#  Include director(y/ies)

## UPDATED to new syntax

include_directories(${CMAKE_SOURCE_DIR}/cint)

# Main eRHIC dictionary.
ROOT_GENERATE_DICTIONARY( erhicDict functions.h ${ERHIC_HEADERS} LINKDEF LinkDef.h )
	# ${ERHIC_HEADERS}
	# ${ROOT_INCLUDE_DIR}
	# ${CMAKE_SOURCE_DIR}/include

# Smearing code dictionary.
ROOT_GENERATE_DICTIONARY( smearDict ${SMEAR_HEADERS} LINKDEF smearLinkDef.h )

# Default complete set of source file (including dictionary files).
# We will add to this if PYTHIA support was requested.
SET(
   ALL_SOURCES
   ${CMAKE_CURRENT_BINARY_DIR}/erhicDict.cxx
   ${CMAKE_CURRENT_BINARY_DIR}/smearDict.cxx
   ${ERHIC_SOURCE} # Includes smear source files
)

## PCM files for installation
list( APPEND
   ALL_DICTS
   ${CMAKE_CURRENT_BINARY_DIR}/erhicDict_rdict.pcm
   ${CMAKE_CURRENT_BINARY_DIR}/smearDict_rdict.pcm
)


# Build PYTHIA6-dependent code if a library directory was specified.
# This comprises PYTHIA-related stuff in the main eRHIC directory,
# plus the hadronic code.
if(PYTHIA6_LIBDIR)
   message ( "PYTHIA6_LIBDIR specified --> Adding PYTHIA sources and headers" )
   # PYTHIA-related header files.
   SET(
      PYTHIA_HEADERS
      eicsmear/erhic/Pythia6.h
      eicsmear/erhic/Pythia6EventBuilder.h
      eicsmear/erhic/Pythia6ParticleBuilder.h
   )

   # Headers related to hadronic code.
   # Because they are referenced by a different LinkDef
   # they are listed separately from PYTHIA_HEADERS.
   SET(
      HADRONIC_HEADERS
      eicsmear/hadronic/EventMC.h
      eicsmear/hadronic/EventPythia.h
      eicsmear/hadronic/ParticleMC.h
      eicsmear/hadronic/Pythia6EventFactory.h
      eicsmear/hadronic/EventSmear.h
   )
   INSTALL(DIRECTORY ${CMAKE_SOURCE_DIR}/include/eicsmear/hadronic DESTINATION include/eicsmear )

   # PYTHIA-related source files (including for hadronic code).
   SET(
      PYTHIA_SOURCE
      src/erhic/Pythia6.cxx
      src/erhic/Pythia6EventBuilder.cxx
      src/erhic/Pythia6ParticleBuilder.cxx
      src/hadronic/EventMC.cxx
      src/hadronic/EventPythia.cxx
      src/hadronic/EventSmear.cxx
      src/hadronic/ParticleMC.cxx
      src/hadronic/Pythia6EventFactory.cxx
   )

   # Create the pythia dictionary.
   ROOT_GENERATE_DICTIONARY( pythiaDict ${PYTHIA_HEADERS} LINKDEF LinkDefPythia.h )

   # Create the hadronic dictionary.
   ROOT_GENERATE_DICTIONARY( hadronicDict ${HADRONIC_HEADERS} LINKDEF hadronicLinkDef.h )

   # Update dictionaries to bne installed   
   list( APPEND
   	 ALL_DICTS
	 ${CMAKE_CURRENT_BINARY_DIR}/pythiaDict_rdict.pcm
   	 ${CMAKE_CURRENT_BINARY_DIR}/hadronicDict_rdict.pcm
	 )

   # Update sources with pythia-dependent files.
   SET(
      ALL_SOURCES
      ${ALL_SOURCES}
      ${CMAKE_CURRENT_BINARY_DIR}/pythiaDict.cxx
      ${CMAKE_CURRENT_BINARY_DIR}/hadronicDict.cxx
      ${PYTHIA_SOURCE} # Includes hadronic source files
   )
endif(PYTHIA6_LIBDIR)

include_directories(include)

add_library(
   eicsmear
   SHARED
   ${ALL_SOURCES}
)

# Add ROOT elements
if(ROOT_FOUND)
   include_directories(${ROOT_INCLUDE_DIR})
   if(PYTHIA6_LIBDIR)
   	target_link_libraries(eicsmear ${ROOT_LIBRARIES} -L${ROOT_LIBRARY_DIR} -lEG -lEGPythia6 -lEve )
   else (PYTHIA6_LIBDIR)
   	target_link_libraries(eicsmear ${ROOT_LIBRARIES} -L${ROOT_LIBRARY_DIR} -lEG )
   endif(PYTHIA6_LIBDIR)
endif(ROOT_FOUND)

INSTALL(TARGETS eicsmear LIBRARY DESTINATION lib)

# Dictionaries for ROOT
INSTALL(FILES ${ALL_DICTS} DESTINATION lib)
