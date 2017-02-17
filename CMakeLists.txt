project(SMESH CXX C)

cmake_minimum_required(VERSION 2.8)

set(SMESH_VERSION_MAJOR 6)
set(SMESH_VERSION_MINOR 7)
set(SMESH_VERSION_PATCH 1)
set(SMESH_VERSION_TWEAK)


#####################
# Set ABI SOVERSION #
#####################
set(SMESH_ABI_SOVERSION 3)
set(SMESH_ABI_VERSION ${SMESH_ABI_SOVERSION}.0.0)

set(BUILD_SHARED_LIBS ON)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(inc)


#################
# Setup options #
#################
option(SMESH_TESTING "Enable unittest" ON)
option(SMESH_BUILD_NETGENPLUGIN "NETGEN support" OFF)

#############
# Unittest  #
#############
if(SMESH_TESTING)
    include(CTest)
    enable_testing()
    subdirs(test)
endif(SMESH_TESTING)

################################
# OCE include and lib location #
################################
set(OCE_REQUIRED_LIBS TKernel TKMath TKBRep TKG2d TKG3d TKMeshVS TKPrim TKGeomBase TKGeomAlgo TKTopAlgo TKShHealing)
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
    set(OCE_LIBRARIES TKernel TKMath TKBRep TKG2d TKG3d TKMeshVS TKPrim TKGeomBase TKGeomAlgo TKTopAlgo TKShHealing)
    # Tells CMake where are the required headers
    include_directories(${OCE_INCLUDE_PATH})
    LINK_DIRECTORIES(${OCE_LIB_PATH})
endif()

###################################
# NETGEN include and lib location #
###################################
if(SMESH_BUILD_NETGENPLUGIN)
    set(NETGEN_NGLIB_INCLUDE_PATH /usr/local/netgen CACHE PATH "Netgen headers location")
    set(NETGEN_LIBRARIES nglib)
    include_directories(${NETGEN_NGLIB_INCLUDE_PATH})
    include_directories(${NETGEN_NGLIB_INCLUDE_PATH}/occ)
    include_directories(${NETGEN_NGLIB_INCLUDE_PATH}/meshing)
    include_directories(${NETGEN_NGLIB_INCLUDE_PATH}/general)
    include_directories(${NETGEN_NGLIB_INCLUDE_PATH}/linalg)
endif(SMESH_BUILD_NETGENPLUGIN)

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
    IF(MSVC)
        ADD_DEFINITIONS(-DWNT -DWIN32 -D_WINDOWS -DCSFDB -DUSE_CLOCK -DMSDOS -DNO_ONEXIT -DNO_My_ctype -DNO_ISATTY -DNO_FPINIT /wd4290 /wd4251 /wd4018 /wd4800 /wd4996 /wd4244 /wd4805 /wd4806 /wd4275 /wd4005 /wd4099 /wd4101 /wd4146 /wd4267 /wd4390 /wd4503 /wd4436)
    ELSE(MSVC)
        ADD_DEFINITIONS(-DWNT -DWIN32 -D_WINDOWS -DCSFDB)
    ENDIF(MSVC)
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
