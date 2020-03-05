message(STATUS "test_missing_particles starts")

# name to refer everywhere
set(EXEC_NAME test_missing_particles)

add_executable(${EXEC_NAME} "")

# (!) GLOB is not recommended by CMake. So it is better to use:
# target_sources_local(${PLUGIN_NAME} PRIVATE your_file.cc ... PUBLIC your_file.h ...)
# But... GLOB here makes this file just hot pluggable
file(GLOB SRC_FILES *.cc)
file(GLOB HEADER_FILES *.h)

## Add detectors
list( APPEND
   SRC_FILES
   ${PROJECT_SOURCE_DIR}/scripts/smearBeAST.cxx
   ${PROJECT_SOURCE_DIR}/scripts/ZEUSDetector.cpp
)

# Add sources to our target
target_sources(${EXEC_NAME} PRIVATE ${SRC_FILES} PUBLIC ${HEADER_FILES})

#target_include_directories(${EXEC_NAME} PUBLIC ${EIC_SMEAR_INCLUDE_DIR})

# Add ROOT elements
if(ROOT_FOUND)
    include_directories(${ROOT_INCLUDE_DIR})
    if(PYTHIA6_LIBDIR)
        target_link_libraries(${EXEC_NAME} ${ROOT_LIBRARIES} -L${ROOT_LIBRARY_DIR} -lEG -lEGPythia6 -lEve )
    else (PYTHIA6_LIBDIR)
        target_link_libraries(${EXEC_NAME} ${ROOT_LIBRARIES} -L${ROOT_LIBRARY_DIR} -lEG )
    endif(PYTHIA6_LIBDIR)
endif(ROOT_FOUND)


# Libraries used for this plugin
target_link_libraries(${EXEC_NAME} eicsmear)
