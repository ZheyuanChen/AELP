#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "sdfc" for configuration "Release"
set_property(TARGET sdfc APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(sdfc PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/SDFC_14.4.7/libsdfc.a"
  )

list(APPEND _cmake_import_check_targets sdfc )
list(APPEND _cmake_import_check_files_for_sdfc "${_IMPORT_PREFIX}/lib/SDFC_14.4.7/libsdfc.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
