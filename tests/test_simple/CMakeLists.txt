message(STATUS "test_simple starts")

add_executable(test_simple_buildtree test_simple_buildtree.cxx)

# Add ROOT elements
if(ROOT_FOUND)
    include_directories(${ROOT_INCLUDE_DIR})
    if(PYTHIA6_LIBDIR)
        target_link_libraries(test_simple_buildtree ${ROOT_LIBRARIES} -L${ROOT_LIBRARY_DIR} -lEG -lEGPythia6 -lEve )
    else (PYTHIA6_LIBDIR)
        target_link_libraries(test_simple_buildtree ${ROOT_LIBRARIES} -L${ROOT_LIBRARY_DIR} -lEG )
    endif(PYTHIA6_LIBDIR)
endif(ROOT_FOUND)


# Libraries used for this plugin
target_link_libraries(test_simple_buildtree eicsmear)

# copy input file
configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/beagle_eD.txt ${CMAKE_CURRENT_BINARY_DIR}/beagle_eD.txt COPYONLY)

