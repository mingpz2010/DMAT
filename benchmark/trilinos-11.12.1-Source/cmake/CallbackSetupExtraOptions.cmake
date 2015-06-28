
# We need to inject the Trilinos/cmake directory to find
# TrilinosCreateClientTemplateHeaders.cmake
SET(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH} "${Trilinos_SOURCE_DIR}/cmake")


MACRO(TRIBITS_REPOSITORY_SETUP_EXTRA_OPTIONS)


  #
  # Change compiler warning flags for just Trilinos-based projects!
  #

  MULTILINE_SET(Trilinos_COMMON_STRONG_COMPILE_WARNING_FLAGS
    " -pedantic" # Adds more static checking to remove non-ANSI GNU extensions
    " -Wall" # Enable a bunch of default warnings
    " -Wno-long-long" # Allow long long int since it is used by MPI, SWIG, etc.
  )

  MULTILINE_SET(Trilinos_C_STRONG_COMPILE_WARNING_FLAGS
    "${Trilinos_COMMON_STRONG_COMPILE_WARNING_FLAGS}"
    " -std=c99" # Check for C99
  )

  MULTILINE_SET(Trilinos_CXX_STRONG_COMPILE_WARNING_FLAGS
    "${Trilinos_COMMON_STRONG_COMPILE_WARNING_FLAGS}"
    " -std=c++98" # C++98 standard code
    " -Wwrite-strings" # Checks for non-const char * copy of string constants
  )

  # NOTE: By using "Trilinos_" and *NOT* "${PROJECT_NAME}_", we ensure that we
  # only change the flags for Trilinos-based builds.

  #MESSAGE("TRIBITS_REPOSITORY_SETUP_EXTRA_OPTIONS got called!")

  SET(TPL_ENABLE_MPI OFF CACHE BOOL "Enable MPI support.")

  ADVANCED_SET(Trilinos_DATA_DIR  NOTFOUND
    CACHE PATH
    "Path TrilinosData directory to find more tests and other stuff" )
    
  IF (NOT ${PROJECT_NAME}_ENABLE_Fortran)
    MESSAGE(
      "\n***"
      "\n*** Warning: Setting ${PROJECT_NAME}_ENABLE_ForTrilinos=OFF"
      " because ${PROJECT_NAME}_ENABLE_Fortran=OFF!"
      "\n***\n"
      )
    SET(${PROJECT_NAME}_ENABLE_ForTrilinos OFF)
  ENDIF()

  IF ("${${PROJECT_NAME}_ENABLE_PyTrilinos}" STREQUAL "" AND NOT BUILD_SHARED_LIBS)
    MESSAGE(
      "\n***"
      "\n*** Warning: Setting ${PROJECT_NAME}_ENABLE_PyTrilinos=OFF"
      " because BUILD_SHARED_LIBS=OFF!"
      "\n***\n"
      )
    SET(${PROJECT_NAME}_ENABLE_PyTrilinos OFF)
  ENDIF()

  IF (NOT EXISTS "${Trilinos_SOURCE_DIR}/packages/TriKota/Dakota")
    MESSAGE("-- " "  Setting ${PROJECT_NAME}_ENABLE_TriKota=OFF"
      " because '${Trilinos_SOURCE_DIR}/packages/TriKota/Dakota' does not exit!")
    SET(${PROJECT_NAME}_ENABLE_TriKota OFF)
  ENDIF()
    
  # ToDo: What is this and why is it needed?
  SET(TRILINOS_BUILD_SHARED_LIBS "@BUILD_SHARED_LIBS@")

ENDMACRO()
