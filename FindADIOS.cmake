# - Find ADIOS based on environment setup by adios module

# Based on FindHDF5.cmake to parse environment variable
# Looks up all the libraries in the environment variable for the ADIOS
# module. 



# Parse a compile line for definitions, includes, library paths, and libraries.
# This from CMake-2.8.0-rc6/Module/FindHDF5.cmake
macro( _ADIOS_parse_module_line 
        module_line 
        include_paths
        library_paths
        libraries )

    # Match the include paths
    string( REGEX MATCHALL "-I([^\" ]+)" include_path_flags 
        "${module_line}" 
        )
    foreach( IPATH ${include_path_flags} )
        string( REGEX REPLACE "^-I" "" IPATH ${IPATH} )
        string( REGEX REPLACE "//" "/" IPATH ${IPATH} )
        list( APPEND ${include_paths} ${IPATH} )
    endforeach()

    # Match the definitions
    string( REGEX MATCHALL "-D[^ ]*" definition_flags "${module_line}" )
    foreach( DEF ${definition_flags} )
        list( APPEND ${definitions} ${DEF} )
    endforeach()

    # Match the library paths
    string( REGEX MATCHALL "-L([^\" ]+|\"[^\"]+\")" library_path_flags
        "${module_line}"
        )

    foreach( LPATH ${library_path_flags} )
        string( REGEX REPLACE "^-L" "" LPATH ${LPATH} )
        string( REGEX REPLACE "//" "/" LPATH ${LPATH} )
        list( APPEND ${library_paths} ${LPATH} )
    endforeach()

    # now search for the library names specified in the compile line (match -l...)
    # match only -l's preceded by a space or comma
    # this is to exclude directory names like xxx-linux/
    string( REGEX MATCHALL "[, ]-l([^\", ]+)" library_name_flags
        "${module_line}" )
    # strip the -l from all of the library flags and add to the search list
    foreach( LIB ${library_name_flags} )
        string( REGEX REPLACE "^[, ]-l" "" LIB ${LIB} )
        list( APPEND ${libraries} ${LIB} )
    endforeach()
endmacro()

if( ADIOS_INCLUDE_DIR AND ADIOS_LIBRARIES )
    # Already setup
else()

    # Parse environment variables for a list of libraries and a 
    # list of places to look

    _ADIOS_parse_module_line( $ENV{ADIOS_FLIB} 
    ADIOS_INCLUDE_FLAGS
    ADIOS_LIBRARY_DIRS
    ADIOS_LIBRARY_NAMES
    )

_ADIOS_parse_module_line( $ENV{ADIOS_INC} 
ADIOS_INCLUDE_FLAGS
ADIOS_LIBRARY_DIRS
ADIOS_LIBRARY_NAMES
)

    _ADIOS_parse_module_line( $ENV{ADIOSREAD_FLIB} 
    ADIOS_INCLUDE_FLAGS
    ADIOS_LIBRARY_DIRS
    ADIOS_LIBRARY_NAMES
    )

#MESSAGE(" Parsing resulted in include:" ${ADIOS_INCLUDE_FLAGS} )
#MESSAGE(" Parsing resulted in library dirs:" ${ADIOS_LIBRARY_DIRS} )
#MESSAGE(" Parsing resulted in library names:" ${ADIOS_LIBRARY_NAMES} )

SET( ADIOS_FOUND TRUE )
# Look for each library identified
foreach ( LIB ${ADIOS_LIBRARY_NAMES} )
    if(  ${LIB} MATCHES "m$"  )
	list( APPEND ADIOS_LIBRARIES "-lm" )
    else()
	if(  ${LIB} MATCHES "z$"  )
	    list( APPEND ADIOS_LIBRARIES "-lz" )
	else()
	    # Look for shared library first
	    find_library( ADIOS_${LIB}_LIBRARY 
		NAMES lib${LIB}.a
		HINTS ${ADIOS_LIBRARY_DIRS}
		ENV ${ADIOS_DIR}
		PATH_SUFFIXES lib64 lib )
	    # If no shared library, look for static library
	    if( NOT ADIOS_${LIB}_LIBRARY )
		find_library( ADIOS_${LIB}_LIBRARY 
		    NAMES ${LIB}
		    HINTS ${ADIOS_LIBRARY_DIRS}
		    ENV ${ADIOS_DIR}
		    PATH_SUFFIXES lib64 lib )
	    endif( NOT ADIOS_${LIB}_LIBRARY )
	    if( ADIOS_${LIB}_LIBRARY )
		message(" Added ${ADIOS_${LIB}_LIBRARY} ")
		list( APPEND ADIOS_LIBRARIES ${ADIOS_${LIB}_LIBRARY} )
	    else()
		MESSAGE( "Unable to find library lib${LIB}.a" )
		SET( ADIOS_FOUND FALSE )
	    endif()

	endif()
    endif()
endforeach()

foreach (INC ${ADIOS_INCLUDE_FLAGS})
    list( APPEND ADIOS_INCLUDE_DIR ${INC} )
endforeach()

list( APPEND ${ADIOS_LIBRARIES}  )

endif()

# Give some feedback
IF (ADIOS_FOUND)
    IF (NOT ADIOS_FIND_QUIETLY)
        MESSAGE(STATUS "Found ADIOS: ${ADIOS_LIBRARIES}")
        MESSAGE(STATUS "Found ADIOS INCLUDE: ${ADIOS_INCLUDE_DIR}")
    ENDIF (NOT ADIOS_FIND_QUIETLY)
ELSE (ADIOS_FOUND)
    IF (ADIOS_FIND_REQUIRED)
        MESSAGE(FATAL_ERROR "Could not find ADIOS or a required component")
    ENDIF (ADIOS_FIND_REQUIRED)
ENDIF (ADIOS_FOUND)
