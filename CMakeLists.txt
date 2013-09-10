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

cmake_minimum_required(VERSION 2.7)

project(EICSMEAR)

# Add this project's module path to the default path.
# Needed to locate this project's FindROOT.cmake.
set(
   CMAKE_MODULE_PATH
   ${CMAKE_MODULE_PATH}
   ${CMAKE_SOURCE_DIR}/cmake/Modules/
)

# Find ROOT package sets the following.
# ROOT_FOUND          If the ROOT is found
# ROOT_INCLUDE_DIR    PATH to the include directory
# ROOT_LIBRARIES      Most common libraries
# ROOT_LIBRARY_DIR    PATH to the library directory 
# It also provides ROOT_GENERATE_DICTIONARY, which is used below.
find_package(ROOT)

# List all main eRHIC headers EXCEPT PYTHIA-related files.
# We don't need to prefix "${CMAKE_SOURCE_DIR}/include/" before all
# the header names for the main library (cmake figures it out) BUT
# we do for it to work nicely with the ROOT dictionary creation.
# Otherwise the dictionary-generation code doesn't seem to locate
# the header files properly.
SET(
   ERHIC_HEADERS
   ${CMAKE_SOURCE_DIR}/include/eicsmear/erhic/BeamParticles.h
   ${CMAKE_SOURCE_DIR}/include/eicsmear/erhic/DisKinematics.h
   ${CMAKE_SOURCE_DIR}/include/eicsmear/erhic/EventBase.h
   ${CMAKE_SOURCE_DIR}/include/eicsmear/erhic/EventDis.h
   ${CMAKE_SOURCE_DIR}/include/eicsmear/erhic/EventDjangoh.h
   ${CMAKE_SOURCE_DIR}/include/eicsmear/erhic/EventDpmjet.h
   ${CMAKE_SOURCE_DIR}/include/eicsmear/erhic/EventFactory.h
   ${CMAKE_SOURCE_DIR}/include/eicsmear/erhic/EventGmcTrans.h
   ${CMAKE_SOURCE_DIR}/include/eicsmear/erhic/EventMC.h
   ${CMAKE_SOURCE_DIR}/include/eicsmear/erhic/EventMCFilterABC.h
   ${CMAKE_SOURCE_DIR}/include/eicsmear/erhic/EventMilou.h
   ${CMAKE_SOURCE_DIR}/include/eicsmear/erhic/EventPepsi.h
   ${CMAKE_SOURCE_DIR}/include/eicsmear/erhic/EventPythia.h
   ${CMAKE_SOURCE_DIR}/include/eicsmear/erhic/EventRapgap.h
   ${CMAKE_SOURCE_DIR}/include/eicsmear/erhic/File.h
   ${CMAKE_SOURCE_DIR}/include/eicsmear/erhic/Forester.h
   ${CMAKE_SOURCE_DIR}/include/eicsmear/erhic/Kinematics.h
   ${CMAKE_SOURCE_DIR}/include/eicsmear/erhic/Particle.h
   ${CMAKE_SOURCE_DIR}/include/eicsmear/erhic/ParticleIdentifier.h
   ${CMAKE_SOURCE_DIR}/include/eicsmear/erhic/ParticleMC.h
   ${CMAKE_SOURCE_DIR}/include/eicsmear/erhic/Pid.h
   ${CMAKE_SOURCE_DIR}/include/eicsmear/erhic/VirtualEvent.h
   ${CMAKE_SOURCE_DIR}/include/eicsmear/erhic/VirtualParticle.h
)
INSTALL(FILES ${ERHIC_HEADERS} DESTINATION include/eicsmear/erhic)
INSTALL(FILES ${CMAKE_SOURCE_DIR}/include/eicsmear/functions.h
        DESTINATION include/eicsmear)

# Smearing-related headers.
# These are referenced by a different LinkDef to ERHIC_HEADERS so
# we list the separately.
SET(
   SMEAR_HEADERS
   ${CMAKE_SOURCE_DIR}/include/eicsmear/smear/Acceptance.h
   ${CMAKE_SOURCE_DIR}/include/eicsmear/smear/Bremsstrahlung.h
   ${CMAKE_SOURCE_DIR}/include/eicsmear/smear/Detector.h
   ${CMAKE_SOURCE_DIR}/include/eicsmear/smear/Device.h
   ${CMAKE_SOURCE_DIR}/include/eicsmear/smear/Distributor.h
   ${CMAKE_SOURCE_DIR}/include/eicsmear/smear/EventS.h
   ${CMAKE_SOURCE_DIR}/include/eicsmear/smear/EventSmear.h
   ${CMAKE_SOURCE_DIR}/include/eicsmear/smear/EventDisFactory.h
   ${CMAKE_SOURCE_DIR}/include/eicsmear/smear/FormulaString.h
   ${CMAKE_SOURCE_DIR}/include/eicsmear/smear/functions.h
   ${CMAKE_SOURCE_DIR}/include/eicsmear/smear/ParticleID.h
   ${CMAKE_SOURCE_DIR}/include/eicsmear/smear/ParticleMCS.h
   ${CMAKE_SOURCE_DIR}/include/eicsmear/smear/PerfectID.h
   ${CMAKE_SOURCE_DIR}/include/eicsmear/smear/PlanarTracker.h
   ${CMAKE_SOURCE_DIR}/include/eicsmear/smear/RadialTracker.h
   ${CMAKE_SOURCE_DIR}/include/eicsmear/smear/Smear.h
   ${CMAKE_SOURCE_DIR}/include/eicsmear/smear/Smearer.h
   ${CMAKE_SOURCE_DIR}/include/eicsmear/smear/Tracker.h
)
INSTALL(FILES ${SMEAR_HEADERS} DESTINATION include/eicsmear/smear)

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

# Main eRHIC dictionary.
ROOT_GENERATE_DICTIONARY(
   "${CMAKE_SOURCE_DIR}/include/eicsmear/functions.h;${ERHIC_HEADERS}"
   ${CMAKE_SOURCE_DIR}/cint/LinkDef.h
   ${CMAKE_CURRENT_BINARY_DIR}/erhicDict.cxx
   "${ROOT_INCLUDE_DIR};${CMAKE_SOURCE_DIR}/include"
)

# Smearing code dictionary.
ROOT_GENERATE_DICTIONARY(
   "${SMEAR_HEADERS}"
   ${CMAKE_SOURCE_DIR}/cint/smearLinkDef.h
   ${CMAKE_CURRENT_BINARY_DIR}/smearDict.cxx
   "${ROOT_INCLUDE_DIR};${CMAKE_SOURCE_DIR}/include"
)

# Default complete set of source file (including dictionary files).
# We will add to this if PYTHIA support was requested.
SET(
   ALL_SOURCES
   ${CMAKE_CURRENT_BINARY_DIR}/erhicDict.cxx
   ${CMAKE_CURRENT_BINARY_DIR}/smearDict.cxx
   ${ERHIC_SOURCE} # Includes smear source files
)

# Build PYTHIA6-dependent code if a library directory was specified.
# This comprises PYTHIA-related stuff in the main eRHIC directory,
# plus the hadronic code.
if(PYTHIA6_LIBDIR)
   # PYTHIA-related header files.
   SET(
      PYTHIA_HEADERS
      ${CMAKE_SOURCE_DIR}/include/eicsmear/erhic/Pythia6.h
      ${CMAKE_SOURCE_DIR}/include/eicsmear/erhic/Pythia6EventBuilder.h
      ${CMAKE_SOURCE_DIR}/include/eicsmear/erhic/Pythia6ParticleBuilder.h
   )
   INSTALL(FILES ${PYTHIA_HEADERS} DESTINATION include/eicsmear/erhic)

   # Headers related to hadronic code.
   # Because they are referenced by a different LinkDef
   # they are listed separately from PYTHIA_HEADERS.
   SET(
      HADRONIC_HEADERS
      ${CMAKE_SOURCE_DIR}/include/eicsmear/hadronic/EventMC.h
      ${CMAKE_SOURCE_DIR}/include/eicsmear/hadronic/EventPythia.h
      ${CMAKE_SOURCE_DIR}/include/eicsmear/hadronic/ParticleMC.h
      ${CMAKE_SOURCE_DIR}/include/eicsmear/hadronic/Pythia6EventFactory.h
      ${CMAKE_SOURCE_DIR}/include/eicsmear/hadronic/EventSmear.h
   )
   INSTALL(FILES ${HADRONIC_HEADERS} DESTINATION include/eicsmear/hadronic)

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
   ROOT_GENERATE_DICTIONARY(
      "${PYTHIA_HEADERS}"
      ${CMAKE_SOURCE_DIR}/cint/LinkDefPythia.h
      ${CMAKE_CURRENT_BINARY_DIR}/pythiaDict.cxx
      "${ROOT_INCLUDE_DIR};${CMAKE_SOURCE_DIR}/include"
   )

   # Create the hadronic dictionary.
   ROOT_GENERATE_DICTIONARY(
      "${HADRONIC_HEADERS}"
      ${CMAKE_SOURCE_DIR}/cint/hadronicLinkDef.h
      ${CMAKE_CURRENT_BINARY_DIR}/hadronicDict.cxx
      "${ROOT_INCLUDE_DIR};${CMAKE_SOURCE_DIR}/include"
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
   target_link_libraries(eicsmear ${ROOT_LIBRARIES} -lEGPythia6 -lEve -lEG)
endif(ROOT_FOUND)

INSTALL(TARGETS eicsmear LIBRARY DESTINATION lib)
