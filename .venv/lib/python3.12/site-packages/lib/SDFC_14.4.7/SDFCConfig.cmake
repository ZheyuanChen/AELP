set(SDFC_VERSION 14.4.7)
include("${CMAKE_CURRENT_LIST_DIR}/SDFCTargets.cmake")

####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was SDFCConfig.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

####################################################################################
set_and_check(SDFC_INCLUDE_DIR "${PACKAGE_PREFIX_DIR}/include/SDFC_14.4.7")
