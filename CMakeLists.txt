project(SMESH CXX C)

cmake_minimum_required(VERSION 2.8)

set(SMESH_VERSION_MAJOR 6)
set(SMESH_VERSION_MINOR 7)
set(SMESH_VERSION_PATCH 6)
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
option(SMESH_EXTRA_WARNINGS "Enable Maximum Warnings Level" OFF)
option(SMESH_USE_BUNDLED_BOOST "Use boost-1.52 headers bundled with SMESH" OFF)
#option(SMESH_BUILD_NETGENPLUGIN "NETGEN support" OFF)

##################
# Extra warnings #
##################
# Enable Extra warnings for each compiler
if(SMESH_EXTRA_WARNINGS)
    if(MSVC)
        add_definitions("/W4")
    endif(MSVC)
    if(CMAKE_COMPILER_IS_GNUCXX)
        add_definitions("-Wall -Wextra")
    endif(CMAKE_COMPILER_IS_GNUCXX)
    if(CLANG)
        add_definitions("-Weverything")
    endif(CLANG)
endif(SMESH_EXTRA_WARNINGS)

# raise exception if MEFISTO2 mesh fails
if(RAISE_EXCEPTION_ON_FAILURE)
    add_definitions("-DRAISE_EXCEPTION_ON_FAILURE")
endif(RAISE_EXCEPTION_ON_FAILURE)

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
if(MSVC)
    message(STATUS "Don't look for installed OCE, use provided paths")
else(MSVC)
    message(STATUS "Looking for OCE ...")
    find_package(OCE
          COMPONENTS ${OCE_REQUIRED_LIBS}
          CONFIG)
endif(MSVC)
if(OCE_FOUND)
    message(STATUS "OCE found.")
    include_directories(${OCE_INCLUDE_DIRS})
    set(OCE_LIBRARIES ${OCE_REQUIRED_LIBS})
else()
    set(OCE_INCLUDE_PATH /usr/local/include/oce CACHE
        PATH "OCE headers location" )
    set(OCE_LIB_PATH /usr/local/lib CACHE
        PATH "OCE libraries location" )
    set(OCE_LIBRARIES TKernel TKMath TKBRep TKG2d TKG3d TKMeshVS TKPrim TKGeomBase TKGeomAlgo TKTopAlgo TKShHealing)
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
if(SMESH_USE_BUNDLED_BOOST)
    message(STATUS "Using bundled 1.52.0 headers.")
    include_directories(${CMAKE_SOURCE_DIR}/boost-1.52.0)
else(SMESH_USE_BUNDLED_BOOST)
    find_package(Boost)
    if(Boost_FOUND)
        include_directories(${Boost_INCLUDE_DIRS})
    else()
        message(STATUS "Boost not found, using bundled 1.52.0 headers.")
        include_directories(${CMAKE_SOURCE_DIR}/boost-1.52.0)
    endif(Boost_FOUND)
endif(SMESH_USE_BUNDLED_BOOST)

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

##############################
# Specific Appe OSX settings #
##############################
# RPATH
if(APPLE)
    SET(CMAKE_MACOSX_RPATH ON)
    # use, i.e. don't skip the full RPATH for the build tree
    SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
    # when building, don't use the install RPATH already
    # (but later on when installing)
    SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
    SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
    # add the automatically determined parts of the RPATH
    # which point to directories outside the build tree to the install RPATH
    SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
    # the RPATH to be used when installing, but only if it's not a system directory
    LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
    IF("${isSystemDir}" STREQUAL "-1")
       SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
    ENDIF("${isSystemDir}" STREQUAL "-1")
endif(APPLE)

# MACOSX_DEPLOYEMENT_TARGET
if(APPLE)
    if(DEFINED ENV{MACOSX_DEPLOYMENT_TARGET})
      message(STATUS "MACOSX_DEPLOYMENT_TARGET set to: " $ENV{MACOSX_DEPLOYMENT_TARGET})
      if($ENV{MACOSX_DEPLOYMENT_TARGET} MATCHES "10.9")
          message(STATUS "Link to libc++")
          add_definitions(-stdlib=libc++)
      endif($ENV{MACOSX_DEPLOYMENT_TARGET} MATCHES "10.9")
    endif(DEFINED ENV{MACOSX_DEPLOYMENT_TARGET})
endif(APPLE)

#################
# SMESH library #
#################
add_subdirectory(src)


#######################
# SMESH Documentation #
#######################
add_subdirectory(doc/src)
