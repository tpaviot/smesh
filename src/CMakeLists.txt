##########################
# SMESH lib source files #
##########################

set(SMESH_LIBS ${OCE_LIBRARIES})

##########
# Driver #
##########
FILE(GLOB Driver_source_files Driver/*.cpp inc/Driver_*.hxx)
ADD_LIBRARY(Driver SHARED ${Driver_source_files})
TARGET_LINK_LIBRARIES(Driver ${SMESH_LIBS})
set_target_properties(Driver PROPERTIES SOVERSION ${${PROJECT_NAME}_ABI_SOVERSION}
                                         VERSION ${${PROJECT_NAME}_ABI_VERSION})
if(WIN32)
    set_target_properties(Driver PROPERTIES COMPILE_FLAGS "-DMESHDRIVER_EXPORTS -DSMESHUtils_EXPORTS")
endif(WIN32)

#############
# DriverDAT #
#############
FILE(GLOB DriverDAT_source_files DriverDAT/*.cpp ../inc/DriverDAT*.h*)
ADD_LIBRARY(DriverDAT SHARED ${DriverDAT_source_files})
TARGET_LINK_LIBRARIES(DriverDAT ${SMESH_LIBS} SMDS Driver)
set_target_properties(DriverDAT PROPERTIES SOVERSION ${${PROJECT_NAME}_ABI_SOVERSION}
                                         VERSION ${${PROJECT_NAME}_ABI_VERSION})
if(WIN32)
    set_target_properties(DriverDAT PROPERTIES COMPILE_FLAGS "-DMESHDRIVERDAT_EXPORTS")
endif(WIN32)

#############
# DriverSTL #
#############
FILE(GLOB DriverSTL_source_files DriverSTL/*.cpp ../inc/DriverSTL*.h*)
ADD_LIBRARY(DriverSTL SHARED ${DriverSTL_source_files})
TARGET_LINK_LIBRARIES(DriverSTL ${SMESH_LIBS} Driver SMDS ${Boost_LIBRARIES})
set_target_properties(DriverSTL PROPERTIES SOVERSION ${${PROJECT_NAME}_ABI_SOVERSION}
                                         VERSION ${${PROJECT_NAME}_ABI_VERSION})
if(WIN32)
    set_target_properties(DriverSTL PROPERTIES COMPILE_FLAGS "-DMESHDRIVERSTL_EXPORTS -DBASICS_EXPORT -DSMESHUtils_EXPORTS -DBASICS_EXPORTS")
endif(WIN32)

#############
# DriverUNV #
#############
FILE(GLOB DriverUNV_source_files DriverUNV/*.cpp ../inc/DriverUNV*.h*)
ADD_LIBRARY(DriverUNV SHARED ${DriverUNV_source_files})
TARGET_LINK_LIBRARIES(DriverUNV ${SMESH_LIBS} SMESHDS SMDS Driver)
set_target_properties(DriverUNV PROPERTIES SOVERSION ${${PROJECT_NAME}_ABI_SOVERSION}
                                         VERSION ${${PROJECT_NAME}_ABI_VERSION})
if(WIN32)
    set_target_properties(DriverUNV PROPERTIES COMPILE_FLAGS "-DMESHDRIVERUNV_EXPORTS")
endif(WIN32)

#########
# SMDS #
#########
FILE(GLOB SMDS_source_files SMDS/SMDS*.cpp ../inc/SMDS*.hxx)
ADD_LIBRARY(SMDS SHARED ${SMDS_source_files})
TARGET_LINK_LIBRARIES(SMDS ${SMESH_LIBS})
set_target_properties(SMDS PROPERTIES SOVERSION ${${PROJECT_NAME}_ABI_SOVERSION}
                                         VERSION ${${PROJECT_NAME}_ABI_VERSION})
if(WIN32)
    set_target_properties(SMDS PROPERTIES COMPILE_FLAGS "-DSMDS_EXPORTS -DUTILS_EXPORTS")
endif(WIN32)

###########
# SMESHDS #
###########
FILE(GLOB SMESHDS_source_files SMESHDS/*.cpp ../inc/SMESHDS_*.h*)
ADD_LIBRARY(SMESHDS SHARED ${SMESHDS_source_files})
set_target_properties(SMESHDS PROPERTIES SOVERSION ${${PROJECT_NAME}_ABI_SOVERSION}
                                         VERSION ${${PROJECT_NAME}_ABI_VERSION})
TARGET_LINK_LIBRARIES(SMESHDS ${SMESH_LIBS} SMDS)

#########
# SMESH #
#########
FILE(GLOB SMESH_source_files SMESH/*.c SMESH/*.cpp Controls/*.cpp ../inc/SMESH_*.h*)
ADD_LIBRARY(SMESH SHARED ${SMESH_source_files})
TARGET_LINK_LIBRARIES(SMESH SMDS SMESHDS  Driver DriverSTL DriverDAT DriverUNV ${SMESH_LIBS})
set_target_properties(SMESH PROPERTIES SOVERSION ${${PROJECT_NAME}_ABI_SOVERSION}
                                         VERSION ${${PROJECT_NAME}_ABI_VERSION})

if(WIN32)
    set_target_properties(SMESH PROPERTIES COMPILE_FLAGS "-DSMESH_EXPORTS -DSMESHCONTROLS_EXPORTS -DSMESHimpl_EXPORTS -DMESHDriverGMF_EXPORTS -DMESHDRIVERMED_EXPORTS -DMEDWRAPPER_EXPORTS -DMEDWRAPPER_BASE_EXPORTS -DMEDWRAPPER_FACTORY_EXPORTS -DSMESHUtils_EXPORTS -DUTILS_EXPORTS -DMEDWRAPPER_V2_2_EXPORTS")
endif(WIN32)

#############
# MEFISTO2D #
#############
FILE(GLOB MEFISTO2_source_files MEFISTO2/aptrte.cpp MEFISTO2/trte.c)
ADD_LIBRARY(MEFISTO2 SHARED ${MEFISTO2_source_files})
TARGET_LINK_LIBRARIES(MEFISTO2 SMESH ${SMESH_LIBS})
set_target_properties(MEFISTO2 PROPERTIES SOVERSION ${${PROJECT_NAME}_ABI_SOVERSION}
                                         VERSION ${${PROJECT_NAME}_ABI_VERSION})
if(WIN32)
    set_target_properties(MEFISTO2 PROPERTIES COMPILE_FLAGS "-DMEFISTO2D_EXPORTS")
endif(WIN32)

##############
# StdMeshers #
##############
FILE(GLOB StdMeshers_source_files StdMeshers/*.cpp ../inc/StdMeshers_*.h*)
ADD_LIBRARY(StdMeshers SHARED ${StdMeshers_source_files})
TARGET_LINK_LIBRARIES(StdMeshers SMESH MEFISTO2 ${SMESH_LIBS})
set_target_properties(StdMeshers PROPERTIES SOVERSION ${${PROJECT_NAME}_ABI_SOVERSION}
                                         VERSION ${${PROJECT_NAME}_ABI_VERSION})
if(WIN32)
    set_target_properties(StdMeshers PROPERTIES COMPILE_FLAGS "-DSTDMESHERS_EXPORTS -DSMESHUtils_EXPORTS")
endif(WIN32)

################
# NETGENPlugin #
################
if(SMESH_BUILD_NETGENPLUGIN)
    FILE(GLOB NETGENPlugin_source_files NETGENPlugin/*.cpp ../inc/NETGENPlugin_*.h* )
    ADD_LIBRARY(NETGENPlugin SHARED ${NETGENPlugin_source_files})
    TARGET_LINK_LIBRARIES(NETGENPlugin SMDS SMESHDS SMESH StdMeshers ${SMESH_LIBS})
    set_target_properties(NETGENPlugin PROPERTIES SOVERSION ${${PROJECT_NAME}_ABI_SOVERSION}
                                         VERSION ${${PROJECT_NAME}_ABI_VERSION})
    if(WIN32)
        set_target_properties(NETGENPlugin PROPERTIES COMPILE_FLAGS "-DNETGENPLUGIN_EXPORTS -DNO_PARALLEL_THREADS -DOCCGEOMETRY -DNETGEN_VERSION=${NETGEN_VERSION}")
    else()
        set_target_properties(NETGENPlugin PROPERTIES COMPILE_FLAGS "${NETGEN_CXX_FLAGS}")
    endif(WIN32)
endif (SMESH_BUILD_NETGENPLUGIN)

####################
# Install location #
####################
if(WIN32 AND NOT DEFINED CMAKE_INSTALL_PREFIX)
    set(CMAKE_INSTALL_PREFIX "c:/${PROJECT_NAME}-${${PROJECT_NAME}_VERSION}")
endif()

###########
# INSTALL #
###########
include(GNUInstallDirs)
install(TARGETS Driver DriverUNV DriverDAT DriverSTL SMDS SMESHDS SMESH StdMeshers MEFISTO2
        RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
        ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}")
if(SMESH_BUILD_NETGENPLUGIN)
  install(TARGETS NETGENPlugin
          DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif(SMESH_BUILD_NETGENPLUGIN)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/inc/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/smesh)
# install boost bundled headers if needed
if(SMESH_USE_BUNDLED_BOOST OR (NOT SMESH_USE_BUNDLED_BOOST AND NOT Boost_FOUND))
    install(DIRECTORY ${CMAKE_SOURCE_DIR}/boost-1.52.0/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/smesh)
endif(SMESH_USE_BUNDLED_BOOST OR (NOT SMESH_USE_BUNDLED_BOOST AND NOT Boost_FOUND))

