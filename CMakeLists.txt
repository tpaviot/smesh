# SMESH cmake build file
project(SMESH CXX Fortran)
cmake_minimum_required(VERSION 2.8)

set(SMESH_VERSION_MAJOR 6)
set(SMESH_VERSION_MINOR 5)
set(SMESH_VERSION_PATCH 3)
set(SMESH_VERSION_TWEAK 1-dev)


#####################
# Set ABI SOVERSION #
#####################
set(SMESH_ABI_SOVERSION 1)
set(SMESH_ABI_VERSION ${SMESH_ABI_SOVERSION}.0.0)

set(CMAKE_BUILD_TYPE "Release") # By default, build in release mode
set(BUILD_SHARED_LIBS ON)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(inc)


#################
# Setup options #
#################
option(MONOLITHIC_BUILD "Build one monolithic smesh library." ON)


################################
# OCE include and lib location #
################################
set(OCE_REQUIRED_LIBS TKernel TKMath TKBRep TKG2d TKG3d TKMeshVS TKPrim TKGeomBase TKGeomAlgo TKTopAlgo TKShHealing TKAdvTools)
find_package(OCE
    COMPONENTS ${OCE_REQUIRED_LIBS}
    CONFIG)
if(OCE_FOUND)
    message(STATUS "OpenCASCADE Community Edition found.")
    include_directories(${OCE_INCLUDE_DIRS})
    set(OCE_LIBRARIES ${OCE_REQUIRED_LIBS})
else()
    set( OCE_INCLUDE_PATH /usr/local/include/oce CACHE
        PATH "OCE headers location" )
    set( OCE_LIB_PATH /usr/local/lib CACHE
        PATH "OCE libraries location" )
    set(OCE_LIBRARIES TKernel TKMath TKBRep TKG2d TKG3d TKMeshVS TKPrim TKGeomBase TKGeomAlgo TKTopAlgo TKShHealing TKAdvTools)
    # Tells CMake where are the required headers
    include_directories(${OCE_INCLUDE_PATH})
    LINK_DIRECTORIES(${OCE_LIB_PATH})
endif()


############################
# Default to release build #
############################
if(NOT DEFINED CMAKE_BUILD_TYPE )
  set(CMAKE_BUILD_TYPE "Release")
endif()
message(STATUS "Built type will be ${CMAKE_BUILD_TYPE}")


##############################
# Boost headers are required #
##############################
find_package(Boost)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
else()
    message(STATUS "Boost not found, using bundled headers.")
    include_directories(./boost-1.52.0)
endif()


#################
# Check Fortran #
#################
#if(UNIX)
#  include(CMakeDetermineFortranCompiler)
#endif(UNIX)


############################
# Check 32/64 bit platform #
############################
if(${CMAKE_SIZEOF_VOID_P} MATCHES "8") # It is 64bit, 32 bit systems returns 4
    add_definitions(-D_OCC64)
    message(STATUS "Build 64bit")
    set(BIT 64)
else()
    message(STATUS "Build 32bit")
    set(BIT 32)
endif()


######################################################
# Configure Build settings according to the platform #
######################################################
if(UNIX)
    # Same settings are used for both MacOSX and Unix/Linux
    add_definitions(-DHAVE_CONFIG_H -DLIN -DCSFDB -DOCC_CONVERT_SIGNALS)
else(UNIX)
  if(WIN32)
    ######### Windows ###########
    message(STATUS "Windows system detected")
    add_definitions(-DWNT -DWIN32 -D_WINDOWS -DCSFDB -DHAVE_CONFIG_H
                    -DSMDS_EXPORTS -DSMESHimpl_EXPORTS
                    -DSMESH_EXPORTS -DSMESHDS_EXPORTS -DMESHDRIVERUNV_EXPORTS
                    -DMESHDRIVERSTL_EXPORTS -DMESHDRIVERDAT_EXPORTS
                    -DSTDMESHERS_EXPORTS -DSMESHCONTROLS_EXPORTS
                    -DMESHDRIVER_EXPORTS -DMEFISTO2D_EXPORTS)
  else(WIN32)
    message("Unknown platform")
  endif(WIN32)
endif(UNIX)


#################
# SMESH library #
#################
add_subdirectory(src)


#######################
# SMESH Documentation #
#######################
add_subdirectory(doc/src)


########################
# Packaging with CPack #
########################
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "SMESH library")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENCE.lgpl.txt")
set(CPACK_SET_DESTDIR "ON")

# Define version
set(CPACK_PACKAGE_VERSION_MAJOR ${SMESH_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${SMESH_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${SMESH_VERSION_PATCH})
set(CPACK_PACKAGE_NAME ${PROJECT_NAME} )
set(CPACK_PACKAGE_CONTACT "Thomas Paviot <tpaviot@gmail.com>")
include(InstallRequiredSystemLibraries)

include(CPack)
