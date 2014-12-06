##########################
# SMESH lib source files #
##########################

set(Driver_SRCS
    Driver/Driver_Document.cpp
    Driver/Driver_Mesh.cpp
    Driver/Driver_SMDS_Mesh.cpp
    Driver/Driver_SMESHDS_Mesh.cpp
)
set(DriverDAT_SRCS
    #DriverDAT/DAT_Test.cpp
    DriverDAT/DriverDAT_R_SMDS_Mesh.cpp
    DriverDAT/DriverDAT_R_SMESHDS_Document.cpp
    DriverDAT/DriverDAT_R_SMESHDS_Mesh.cpp
    DriverDAT/DriverDAT_W_SMDS_Mesh.cpp
    DriverDAT/DriverDAT_W_SMESHDS_Document.cpp
    DriverDAT/DriverDAT_W_SMESHDS_Mesh.cpp
)
set(DriverSTL_SRCS
    #DriverSTL/STL_Test.cpp
    DriverSTL/DriverSTL_R_SMDS_Mesh.cpp
    DriverSTL/DriverSTL_W_SMDS_Mesh.cpp
)
set(DriverUNV_SRCS
    DriverUNV/DriverUNV_R_SMDS_Mesh.cpp
    DriverUNV/DriverUNV_R_SMESHDS_Document.cpp
    DriverUNV/DriverUNV_R_SMESHDS_Mesh.cpp
    DriverUNV/DriverUNV_W_SMDS_Mesh.cpp
    DriverUNV/DriverUNV_W_SMESHDS_Document.cpp
    DriverUNV/DriverUNV_W_SMESHDS_Mesh.cpp
    DriverUNV/UNV2411_Structure.cpp
    DriverUNV/UNV2412_Structure.cpp
    DriverUNV/UNV2417_Structure.cpp
    #UNV_Test.cpp
    #UNV_Utilities.cpp
)
set(SMDS_SRCS
    #SMDS/Find_SMDS_MemoryLimit.cpp
    SMDS/SMDS_EdgePosition.cpp
    SMDS/SMDS_FaceOfEdges.cpp
    SMDS/SMDS_FaceOfNodes.cpp
    SMDS/SMDS_FacePosition.cpp
    SMDS/SMDS_IteratorOfElements.cpp
    #SMDS/SMDS_MemoryLimit.cpp
    SMDS/SMDS_Mesh.cpp
    SMDS/SMDS_Mesh0DElement.cpp
    SMDS/SMDS_MeshEdge.cpp
    SMDS/SMDS_MeshElement.cpp
    SMDS/SMDS_MeshElementIDFactory.cpp
    SMDS/SMDS_MeshFace.cpp
    SMDS/SMDS_MeshGroup.cpp
    SMDS/SMDS_MeshIDFactory.cpp
    SMDS/SMDS_MeshNode.cpp
    SMDS/SMDS_MeshObject.cpp
    SMDS/SMDS_MeshVolume.cpp
    SMDS/SMDS_PolygonalFaceOfNodes.cpp
    SMDS/SMDS_PolyhedralVolumeOfNodes.cpp
    SMDS/SMDS_Position.cpp
    SMDS/SMDS_QuadraticEdge.cpp
    SMDS/SMDS_QuadraticFaceOfNodes.cpp
    SMDS/SMDS_QuadraticVolumeOfNodes.cpp
    SMDS/SMDS_SpacePosition.cpp
    SMDS/SMDS_VertexPosition.cpp
    SMDS/SMDS_VolumeOfFaces.cpp
    SMDS/SMDS_VolumeOfNodes.cpp
    SMDS/SMDS_VolumeTool.cpp
)
set(SMESH_SRCS
    SMESH/SMESH_0D_Algo.cpp
    SMESH/SMESH_1D_Algo.cpp
    SMESH/SMESH_2D_Algo.cpp
    SMESH/SMESH_3D_Algo.cpp
    SMESH/SMESH_Algo.cpp
    SMESH/SMESH_Block.cpp
    SMESH/SMESH_Exception.cpp
    SMESH/SMESH_Gen.cpp
    SMESH/SMESH_Group.cpp
    SMESH/SMESH_HypoFilter.cpp
    SMESH/SMESH_Hypothesis.cpp
    SMESH/SMESH_Mesh.cpp
    SMESH/SMESH_MeshEditor.cpp
    SMESH/SMESH_MesherHelper.cpp
    SMESH/SMESH_MeshVSLink.cpp
    SMESH/SMESH_Octree.cpp
    SMESH/SMESH_OctreeNode.cpp
    SMESH/SMESH_Pattern.cpp
    SMESH/SMESH_subMesh.cpp
    Controls/SMESH_Controls.cpp
    Controls/SMESHControls.cpp
)
set(SMESHDS_SRCS
    SMESHDS/SMESHDS_Command.cpp
    SMESHDS/SMESHDS_Document.cpp
    SMESHDS/SMESHDS_Group.cpp
    SMESHDS/SMESHDS_GroupBase.cpp
    SMESHDS/SMESHDS_GroupOnGeom.cpp
    SMESHDS/SMESHDS_Hypothesis.cpp
    SMESHDS/SMESHDS_Mesh.cpp
    SMESHDS/SMESHDS_Script.cpp
    SMESHDS/SMESHDS_SubMesh.cpp
)
set(NETGENPlugin_SRCS
    NETGENPlugin/NETGENPlugin_Hypothesis.cpp
    NETGENPlugin/NETGENPlugin_Hypothesis_2D.cpp
    NETGENPlugin/NETGENPlugin_Mesher.cpp
    NETGENPlugin/NETGENPlugin_NETGEN_2D.cpp
    NETGENPlugin/NETGENPlugin_NETGEN_2D3D.cpp
    NETGENPlugin/NETGENPlugin_NETGEN_2D_ONLY.cpp
    NETGENPlugin/NETGENPlugin_NETGEN_3D.cpp
    NETGENPlugin/NETGENPlugin_SimpleHypothesis_2D.cpp
    NETGENPlugin/NETGENPlugin_SimpleHypothesis_3D.cpp
)
set(StdMeshers_SRCS
    #MEFISTO2/trte.for
    MEFISTO2/aptrte.cpp
    MEFISTO2/trte.f
    StdMeshers/StdMeshers_Arithmetic1D.cpp
    StdMeshers/StdMeshers_AutomaticLength.cpp
    StdMeshers/StdMeshers_CompositeHexa_3D.cpp
    StdMeshers/StdMeshers_CompositeSegment_1D.cpp
    StdMeshers/StdMeshers_Deflection1D.cpp
    StdMeshers/StdMeshers_Distribution.cpp
    StdMeshers/StdMeshers_FaceSide.cpp
    StdMeshers/StdMeshers_FixedPoints1D.cpp
    StdMeshers/StdMeshers_Hexa_3D.cpp
    StdMeshers/StdMeshers_LayerDistribution.cpp
    StdMeshers/StdMeshers_LayerDistribution2D.cpp
    StdMeshers/StdMeshers_LengthFromEdges.cpp
    StdMeshers/StdMeshers_LocalLength.cpp
    StdMeshers/StdMeshers_MaxElementArea.cpp
    StdMeshers/StdMeshers_MaxElementVolume.cpp
    StdMeshers/StdMeshers_MaxLength.cpp
    StdMeshers/StdMeshers_MEFISTO_2D.cpp
    StdMeshers/StdMeshers_NotConformAllowed.cpp
    StdMeshers/StdMeshers_NumberOfLayers.cpp
    StdMeshers/StdMeshers_NumberOfLayers2D.cpp
    StdMeshers/StdMeshers_NumberOfSegments.cpp
    StdMeshers/StdMeshers_Penta_3D.cpp
    StdMeshers/StdMeshers_Prism_3D.cpp
    StdMeshers/StdMeshers_Projection_1D.cpp
    StdMeshers/StdMeshers_Projection_2D.cpp
    StdMeshers/StdMeshers_Projection_3D.cpp
    StdMeshers/StdMeshers_ProjectionSource1D.cpp
    StdMeshers/StdMeshers_ProjectionSource2D.cpp
    StdMeshers/StdMeshers_ProjectionSource3D.cpp
    StdMeshers/StdMeshers_ProjectionUtils.cpp
    StdMeshers/StdMeshers_Propagation.cpp
    StdMeshers/StdMeshers_Quadrangle_2D.cpp
    StdMeshers/StdMeshers_QuadrangleParams.cpp
    StdMeshers/StdMeshers_QuadranglePreference.cpp
    StdMeshers/StdMeshers_QuadraticMesh.cpp
    StdMeshers/StdMeshers_QuadToTriaAdaptor.cpp
    StdMeshers/StdMeshers_RadialPrism_3D.cpp
    StdMeshers/StdMeshers_RadialQuadrangle_1D2D.cpp
    StdMeshers/StdMeshers_Regular_1D.cpp
    StdMeshers/StdMeshers_SegmentAroundVertex_0D.cpp
    StdMeshers/StdMeshers_SegmentLengthAroundVertex.cpp
    StdMeshers/StdMeshers_StartEndLength.cpp
    StdMeshers/StdMeshers_TrianglePreference.cpp
    StdMeshers/StdMeshers_UseExisting_1D2D.cpp
)

if(MONOLITHIC_BUILD)
    set(SMESH_LIBRARIES SMESH)
    add_library(SMESH
        ${Driver_SRCS}
        ${DriverDAT_SRCS}
        ${DriverSTL_SRCS}
        ${DriverUNV_SRCS}
        ${SMDS_SRCS}
        ${SMESH_SRCS}
        ${SMESHDS_SRCS}
        ${StdMeshers_SRCS})
    target_link_libraries(SMESH ${OCE_LIBRARIES})
    set_target_properties(SMESH PROPERTIES
        SOVERSION ${${PROJECT_NAME}_ABI_SOVERSION}
        VERSION ${${PROJECT_NAME}_ABI_VERSION})
else()
    set(SMESH_LIBRARIES Driver;DriverDAT;DriverSTL;DriverUNV;SMDS;SMESH;SMESHDS;StdMeshers)

    foreach(_lib ${SMESH_LIBRARIES})
        add_library(${_lib} ${${_lib}_SRCS})
        target_link_libraries(${_lib} ${OCE_LIBRARIES})
        set_target_properties(${_lib} PROPERTIES
            SOVERSION ${${PROJECT_NAME}_ABI_SOVERSION}
            VERSION ${${PROJECT_NAME}_ABI_VERSION})
    endforeach()
    target_link_libraries(DriverDAT SMDS)
    target_link_libraries(DriverSTL SMDS Driver)
    target_link_libraries(DriverUNV SMDS Driver)
    target_link_libraries(SMESH Driver DriverDAT DriverUNV DriverSTL SMESHDS)
    target_link_libraries(SMESHDS SMDS)
    target_link_libraries(StdMeshers SMESH)
endif()

if(UNIX)
    if(${CMAKE_Fortran_COMPILER} MATCHES "f77")
        target_link_libraries(SMESH f2c)
    endif()
endif()

################################################
# Define output path for generated libraries:  #
# platform/compiler-build_type-bits            #
# for instance:                                #
#     ./win32/bin/vc7.1-release-64             #
#     ./win32/bin/vc9-debug-32                 #
#     ./Unix/i386-debug-64                     #
################################################
if(WIN32)
  if(MSVC) # The compiler used is MSVC
    MESSAGE(STATUS "Found MSVC compiler: ${MSVC} ${MSVC_VERSION}")
    set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin/)
  elseif(BORLAND) # The compiler used is BORLAND
    set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/win${BIT}/bin/${CMAKE_BUILD_TYPE})
  else()
    set(LIBRARY_OUTPUT_PATH win${BIT}/bin/${CMAKE_BUILD_TYPE})
  endif()
  set(EXECUTABLE_OUTPUT_PATH ${LIBRARY_OUTPUT_PATH})
else(WIN32)
  if(APPLE)
  set(PLATFORM Darwin)
    else(APPLE)
  set(PLATFORM Unix)
    endif(APPLE)
  set(LIBRARY_OUTPUT_PATH ${PLATFORM}-${CMAKE_SYSTEM_PROCESSOR}-${CMAKE_BUILD_TYPE}-${BIT})
endif(WIN32)
message(STATUS "Library output path: ${LIBRARY_OUTPUT_PATH}")

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
install(TARGETS ${SMESH_LIBRARIES} DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY ${CMAKE_SOURCE_DIR}/inc/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/smesh)
if(NOT Boost_FOUND)
  install(DIRECTORY ${CMAKE_SOURCE_DIR}/boost-1.52.0/
      DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/smesh)
endif()
