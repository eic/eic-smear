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

cmake_minimum_required(VERSION 3.10)
project(eicsmear VERSION 1.1.0 LANGUAGES CXX )

# cmake needs a bit of help to find modules
set(
   CMAKE_MODULE_PATH
   ${CMAKE_MODULE_PATH}
   ${CMAKE_SOURCE_DIR}/cmake/Modules/
   ${CMAKE_SOURCE_DIR}/cmake/
   )

# Need ROOT first to define some macros
find_package(ROOT 6.10 REQUIRED)

# Is HepMC available?
find_package(HepMC3 QUIET HINTS ${HepMC3_DIR} ${HepMC3} )

# Make zlib required. It's almost everywhere, and this avoids ugly ifdefs and such
find_package(ZLIB REQUIRED)

##############################################################################################################

# Main target is the libeicsmear library
# staying away from wildcards and variables
add_library(
   eicsmear
   SHARED
   ${CMAKE_CURRENT_BINARY_DIR}/erhicDict.cxx
   ${CMAKE_CURRENT_BINARY_DIR}/smearDict.cxx
   src/gzstream.cc
   src/functions.cxx
   src/erhic/BeamParticles.cxx
   src/erhic/BuildTree.cxx
   src/erhic/EventDis.cxx
   src/erhic/EventDjangoh.cxx
   src/erhic/EventDpmjet.cxx
   src/erhic/EventFactory.cxx
   src/erhic/EventGmcTrans.cxx
   src/erhic/EventHepMC.cxx
   src/erhic/EventMC.cxx
   src/erhic/EventMilou.cxx
   src/erhic/EventPepsi.cxx
   src/erhic/EventPythia.cxx
   src/erhic/EventRapgap.cxx
   src/erhic/EventSimple.cxx
   src/erhic/EventSartre.cxx
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


# include directories
target_include_directories(eicsmear
  PUBLIC 
  $<INSTALL_INTERFACE:include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>

  PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/cint
  ${ROOT_INCLUDE_DIR}
  )


# Require at least c++11. Could require individual features instead too
target_compile_features(eicsmear PUBLIC cxx_std_11)

# Can add addtional options here
target_compile_options(eicsmear PRIVATE -Wall -g)

# Optional sources if HepMC3 is installed
if(HepMC3_FOUND)
  message ( "-- HepMC3 library found at ${HEPMC3_LIB} --> Adding HepMC sources and headers" )
  target_sources( eicsmear
    PRIVATE
    src/erhic/EventFactoryHepMC.cxx
    )
  target_include_directories(eicsmear
    PRIVATE
    ${HEPMC3_INCLUDE_DIR}
    )
else()
  message ( "-- HepMC3 library not found" )
endif()
message ( "---- You can specify a HepMC3 installation using \"-DHepMC3=[path/to/lib]\" " )


##############################################################################################################
# Generate various CINT dictionaries.

#----------------------------------------------------------------------------
# Arguments are:
# function custom_generate_root_dict( dictionary
#                                    header1 header2 ...
#                                    LINKDEF linkdef1 ...
#                                    OPTIONS opt1...
#                                    INCPATH includepath1...)
# I tried my best to pull the include directories out of the target
# but it chokes on $<INSTALL_INTERFACE>,   $<BUILD_INTERFACE>
# so we're passing them by hand
include( custom_generate_root_dict )


# Main eRHIC dictionary.
# To make hepmc sources optional, build up arguments
set(
  ERHIC_HEADERS
  eicsmear/functions.h
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
  eicsmear/erhic/EventSimple.h
  eicsmear/erhic/EventSartre.h
  eicsmear/erhic/EventHepMC.h
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
set( ERHIC_INCPATHS
  ${CMAKE_SOURCE_DIR}/include
  ${CMAKE_SOURCE_DIR}/cint
  ${ROOT_INCLUDE_DIR}
  )

if(HepMC3_FOUND)
  list( APPEND
    ERHIC_HEADERS
    eicsmear/erhic/EventFactoryHepMC.h
    )
  list( APPEND
    ERHIC_INCPATHS 
    ${HEPMC3_INCLUDE_DIR}
    )
endif()
   
custom_generate_root_dict( erhicDict
  ${ERHIC_HEADERS}
  LINKDEF LinkDef.h
  OPTIONS -rml libeicsmear -rmf libeicsmear.rootmap -s libeicsmear -multiDict
  INCPATH ${ERHIC_INCPATHS}
  )

# Smearing code dictionary.
custom_generate_root_dict( smearDict
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
  LINKDEF smearLinkDef.h
  OPTIONS -rml libeicsmear -rmf libeicsmear.rootmap -s libeicsmear -multiDict
  INCPATH ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/cint ${ROOT_INCLUDE_DIR}
  )

## Build PYTHIA6-dependent functionality if a library directory was specified.
## This comprises PYTHIA-related stuff in the main eRHIC directory,
## plus the hadronic code.
if(PYTHIA6_LIBDIR)
  message ( "-- PYTHIA6_LIBDIR specified --> Adding PYTHIA sources and headers" )

  # Create the pythia dictionary.
  custom_generate_root_dict( pythiaDict
    eicsmear/erhic/Pythia6.h
    eicsmear/erhic/Pythia6EventBuilder.h
    eicsmear/erhic/Pythia6ParticleBuilder.h
    LINKDEF LinkDefPythia.h
    OPTIONS -rml libeicsmear -rmf libeicsmear.rootmap -s libeicsmear -multiDict
    INCPATH ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/cint ${ROOT_INCLUDE_DIR}
  )
  
  # Create the hadronic dictionary.
  custom_generate_root_dict( hadronicDict
    eicsmear/hadronic/EventMC.h
    eicsmear/hadronic/EventPythia.h
    eicsmear/hadronic/ParticleMC.h
    eicsmear/hadronic/Pythia6EventFactory.h
    eicsmear/hadronic/EventSmear.h
    LINKDEF hadronicLinkDef.h
    OPTIONS -rml libeicsmear -rmf libeicsmear.rootmap -s libeicsmear -multiDict
    INCPATH ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/cint ${ROOT_INCLUDE_DIR}
    )
  
  # Additional PYTHIA-related and hadronic source files
  target_sources( eicsmear
    PRIVATE 
    ${CMAKE_CURRENT_BINARY_DIR}/pythiaDict.cxx
    ${CMAKE_CURRENT_BINARY_DIR}/hadronicDict.cxx
    src/erhic/Pythia6.cxx
    src/erhic/Pythia6EventBuilder.cxx
    src/erhic/Pythia6ParticleBuilder.cxx
    src/hadronic/EventMC.cxx
    src/hadronic/EventPythia.cxx
    src/hadronic/EventSmear.cxx
    src/hadronic/ParticleMC.cxx
    src/hadronic/Pythia6EventFactory.cxx
    )
endif(PYTHIA6_LIBDIR)


##############################################################################################################

## Build the library
## this can be improved with newer root versions

target_link_libraries(eicsmear -lz ${ROOT_LIBRARIES} -L${ROOT_LIBRARY_DIR} -lEG -lEGPythia6 -lEve )
if(PYTHIA6_LIBDIR)
  target_link_libraries(eicsmear ${ROOT_LIBRARIES} -L${ROOT_LIBRARY_DIR} -lEG -lEGPythia6 -lEve )
else (PYTHIA6_LIBDIR)
  target_link_libraries(eicsmear ${ROOT_LIBRARIES} -L${ROOT_LIBRARY_DIR} -lEG )
endif(PYTHIA6_LIBDIR)

if(HepMC3_FOUND)
  target_link_libraries(eicsmear ${HEPMC3_LIB} )
endif()

##############################################################################################################

## Build executables

# eic-smear interface
add_executable(eic-smear scripts/eic-smear.cxx)
target_link_libraries(eic-smear eicsmear )
target_include_directories(eic-smear
  PRIVATE
  ${ROOT_INCLUDE_DIR}
  )

# if(${CMAKE_SYSTEM_NAME} MATCHES Linux)
#   target_compile_options(eic-smear PUBLIC -fPIC)
# endif()
install(TARGETS eic-smear DESTINATION bin)


## pythia6:
if ( PYTHIA6_LIBDIR )
  add_executable(compiled_runpythia
    scripts/compiled_runpythia.cpp
    scripts/runpythia.cpp
    )

  target_include_directories(compiled_runpythia
    PRIVATE
    ${ROOT_INCLUDE_DIR}
    )
  
target_link_libraries(compiled_runpythia
    eicsmear ${ROOT_LIBRARIES} -L${ROOT_LIBRARY_DIR} -lEG -lEGPythia6 -lEve
    )
endif(PYTHIA6_LIBDIR)

##############################################################################################################

## Install library
# Could also use include(GNUInstallDirs)
# and then destinations of the form ${CMAKE_INSTALL_INCLUDEDIR}
install(TARGETS eicsmear
  EXPORT eicsmear-export
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  )

## Ideally, one would have a small list of headers
## matched to the target via set_target_properties(eicsmear PROPERTIES  PUBLIC_HEADER eicsmear/functions.h  )
## and then use PUBLIC_HEADER DESTINATION include above.
## However, ROOT wants to see all headers, so we just copy them into place by hand
## Note that this means the pythia6-related headers get installed regardless of
## whether that functionality is compiled in (because they live in the same directories).

## Install headers and pcm files
install (DIRECTORY ${CMAKE_SOURCE_DIR}/include/eicsmear/erhic
  DESTINATION  include/eicsmear
  )
INSTALL(DIRECTORY ${CMAKE_SOURCE_DIR}/include/eicsmear/smear
  DESTINATION DESTINATION  include/eicsmear
  )

INSTALL(FILES ${CMAKE_SOURCE_DIR}/include/eicsmear/functions.h
  DESTINATION include/eicsmear
  )

INSTALL(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/libeicsmear_erhicDict_rdict.pcm
  ${CMAKE_CURRENT_BINARY_DIR}/libeicsmear_smearDict_rdict.pcm
  ${CMAKE_CURRENT_BINARY_DIR}/libeicsmear.rootmap
  DESTINATION lib
  )

if(PYTHIA6_LIBDIR)
  INSTALL(DIRECTORY ${CMAKE_SOURCE_DIR}/include/eicsmear/hadronic
    DESTINATION include/eicsmear
    )
  INSTALL(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/libeicsmear_pythiaDict_rdict.pcm
    ${CMAKE_CURRENT_BINARY_DIR}/libeicsmear_hadronicDict_rdict.pcm
    DESTINATION lib
  )

endif()

## Generate configuration file
install(EXPORT eicsmear-export
  FILE
    eicsmearConfig.cmake
  NAMESPACE
    eicsmear::
  DESTINATION
  cmake
  )

## Final message
message( " !! After installation, you should also install detectors and tests from https://github.com/eic/eicsmeardetectors. ")

  
# Version check no longer needed b/c c++11 is now always required.
# But leaving it in for future needs
# string(SUBSTRING ${ROOT_VERSION} 0 1 ROOT_MAJOR)
# message ( " ROOT_VERSION = " ${ROOT_VERSION})
# message ( " ROOT_MAJOR = " ${ROOT_MAJOR})
# IF(${ROOT_MAJOR} GREATER 5)
#   if(NOT DEFINED CMAKE_CXX_STANDARD)
#     set(CMAKE_CXX_STANDARD 11)
#   endif()
#ENDIF()

