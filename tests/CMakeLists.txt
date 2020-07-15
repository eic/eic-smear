message(STATUS "Adding tests")

## detector script sources
SET(
  DET_HEAD
  ePHENIXDetector.h
  )
SET(
  DET_SRC  
  ${PROJECT_SOURCE_DIR}/scripts/smearHandBook.cxx
  ${PROJECT_SOURCE_DIR}/scripts/PerfectDetector.cxx
  ${PROJECT_SOURCE_DIR}/scripts/smearNaiveHandBook.cxx
  ${PROJECT_SOURCE_DIR}/scripts/smearBeAST.cxx
  ${PROJECT_SOURCE_DIR}/scripts/ZEUSDetector.cpp
  ${PROJECT_SOURCE_DIR}/scripts/ePHENIXDetector.cpp
  )

## Particle gun
add_executable(particlegun particlegun.cxx)
target_sources(particlegun PRIVATE ${DET_SRC} PUBLIC ${DET_HEAD})

## QA plots
add_executable(qaplots qaplots.cxx qaplots.hh )
target_sources(qaplots PRIVATE ${DET_SRC} PUBLIC ${DET_HEAD} )

## simple tree builder
add_executable(test_simple_buildtree test_simple_buildtree.cxx)

# Add libraries
include_directories(${ROOT_INCLUDE_DIR})
target_link_libraries(particlegun ${ROOT_LIBRARIES} -L${ROOT_LIBRARY_DIR} -lEG )
target_link_libraries(qaplots ${ROOT_LIBRARIES} -L${ROOT_LIBRARY_DIR} -lEG )
target_link_libraries(test_simple_buildtree ${ROOT_LIBRARIES} -L${ROOT_LIBRARY_DIR} -lEG )     

if(PYTHIA6_LIBDIR)
  target_link_libraries(particlegun -lEGPythia6 -lEve )
  target_link_libraries(qaplots -lEGPythia6 -lEve )
  target_link_libraries(test_simple_buildtree -lEGPythia6 -lEve )     
endif(PYTHIA6_LIBDIR)

# Add libraries
target_link_libraries(particlegun eicsmear)
target_link_libraries(qaplots eicsmear)
target_link_libraries(test_simple_buildtree eicsmear)

# copy input files
configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/ep_lowQ2.20x250.small.txt ${CMAKE_CURRENT_BINARY_DIR}/ep_lowQ2.20x250.small.txt COPYONLY)
configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/ep_hiQ2.20x250.small.txt ${CMAKE_CURRENT_BINARY_DIR}/ep_hiQ2.20x250.small.txt COPYONLY)
configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/epref_qaplotsHANDBOOK.pdf ${CMAKE_CURRENT_BINARY_DIR}/epref_qaplotsHANDBOOK.pdf COPYONLY)

configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/beagle_eD.txt ${CMAKE_CURRENT_BINARY_DIR}/beagle_eD.txt COPYONLY)
configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/simple_gen.txt ${CMAKE_CURRENT_BINARY_DIR}/simple_gen.txt COPYONLY)
configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/ep_sartre.txt ${CMAKE_CURRENT_BINARY_DIR}/ep_sartre.txt COPYONLY)

configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/DIS_hepmc.txt ${CMAKE_CURRENT_BINARY_DIR}/DIS_hepmc.txt COPYONLY)
configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/NC_hepmc.txt ${CMAKE_CURRENT_BINARY_DIR}/NC_hepmc.txt COPYONLY)

# execute_process(
#   COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/ep_noradcorr.20x250.quicktest.txt.gz ${CMAKE_CURRENT_BINARY_DIR}
#   COMMAND gunzip ${CMAKE_CURRENT_BINARY_DIR}/ep_noradcorr.20x250.quicktest.txt.gz
#   )
