##############################################################################################################
# Helper pulled out of FindROOT and customized.
#----------------------------------------------------------------------------
# Arguments are:
# function custom_generate_root_dict( dictionary
#                                    header1 header2 ...
#                                    LINKDEF linkdef1 ...
#                                    OPTIONS opt1...
#                                    INCPATH includepath1...)

function(custom_generate_root_dict dictionary)
  # CMAKE_PARSE_ARGUMENTS(ARG "" "" "LINKDEF;OPTIONS" "" ${ARGN})
  CMAKE_PARSE_ARGUMENTS(ARG "" "" "INCPATH;LINKDEF;OPTIONS" "" ${ARGN})
  
  #---Get the list of include directories------------------
  #get_directory_property(incdirs INCLUDE_DIRECTORIES)
  #get_target_property(incdirs eicsmear INCLUDE_DIRECTORIES)
  
  set(includedirs)
  # foreach( d ${incdirs})
  foreach( d ${ARG_INCPATH})
    set(includedirs ${includedirs} -I${d})
  endforeach()
  #---Get the list of header files-------------------------
  ## Changes Kolja Kauder:
  ## - Catch notFound - for some reason REQUIRED doesn't work
  ## - do NOT use the find_file result --> we want no hard-coded paths
  ## - keep track of absolute paths for dependency generation
  ##   (which magically transforms it back to a relative path - but it works)
  ## Note: Not fixed in the glob version, but I can't figure out how to invoke that anyway
  ## If fixable, insert this back into the loop
  #   if(${fp} MATCHES "[*?]") # Is this header a globbing expression?
  #   file(GLOB files ${fp})
  #   foreach(f ${files})
  #     if(NOT f MATCHES LinkDef) # skip LinkDefs from globbing result
  #       set(headerfiles ${headerfiles} ${f})
  #     endif()
  #   endforeach()
  # else()
      
  set(headerfiles)
  set(absheaderfiles)
  foreach(fp ${ARG_UNPARSED_ARGUMENTS})
    find_file(headerFile ${fp} HINTS ${ARG_INCPATH} NO_DEFAULT_PATH )
    if ( ${headerFile} STREQUAL "headerFile-NOTFOUND" )
      message (FATAL_ERROR "GENERATE-DICTIONARY ${fp} not found." )
    else()
      set(headerfiles ${headerfiles} ${fp})
      set(absheaderfiles ${absheaderfiles} ${headerFile})	
    endif()
    unset(headerFile CACHE)
  endforeach()

  #---Get LinkDef.h file------------------------------------
  set(linkdefs)
  foreach( f ${ARG_LINKDEF})
    find_file(linkFile ${f} HINTS ${ARG_INCPATH})
    set(linkdefs ${linkdefs} ${linkFile})
    unset(linkFile CACHE)
  endforeach()
  
  #---call rootcling------------------------------------------
  add_custom_command(OUTPUT ${dictionary}.cxx
    COMMAND ${ROOTCLING_EXECUTABLE}
    -f ${dictionary}.cxx
    ${ARG_OPTIONS}
    ${includedirs}
    ${headerfiles}
    ${linkdefs}
    DEPENDS
    ${absheaderfiles} ${linkdefs} VERBATIM
    )

  ## Not even trying to maintain root5 compatibility here
  # if ( ${ROOT_MAJOR} SMALLER 5 )
  #   add_custom_command(OUTPUT ${dictionary}.cxx
  #     COMMAND ${ROOTCLING_EXECUTABLE} -f ${dictionary}.cxx
  #     -c -DHAVE_CONFIG_H -inlineInputHeader 
  #     ${ARG_OPTIONS} ${includedirs} ${headerfiles} ${linkdefs}
  #     DEPENDS ${absheaderfiles} ${linkdefs} VERBATIM)

    
endfunction()
