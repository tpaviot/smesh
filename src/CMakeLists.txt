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

if(NETGEN_PLUGIN)
    set(NETGENPlugin_SRCS
        NETGENPlugin/NETGENPlugin_Hypothesis.cpp
        NETGENPlugin/NETGENPlugin_Hypothesis_2D.cpp
        NETGENPlugin/NETGENPlugin_Mesher.cpp
        NETGENPlugin/NETGENPlugin_NETGEN_2D.cpp
        NETGENPlugin/NETGENPlugin_NETGEN_2D3D.cpp
        NETGENPlugin/NETGENPlugin_NETGEN_2D_ONLY.cpp
        NETGENPlugin/NETGENPlugin_NETGEN_3D.cpp
        NETGENPlugin/NETGENPlugin_SimpleHypothesis_2D.cpp
        NETGENPlugin/NETGENPlugin_SimpleHypothesis_3D.cpp)
else()
    set(NETGENPlugin_SRCS "")
endif()

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

if(NETGEN_PLUGIN AND NETGEN_INTERNAL)
    set(NETGEN_SRCS
      ../netgen-5.3.1/libsrc/csg/algprim.cpp
      ../netgen-5.3.1/libsrc/csg/brick.cpp
      ../netgen-5.3.1/libsrc/csg/bspline2d.cpp
      ../netgen-5.3.1/libsrc/csg/csgeom.cpp
      ../netgen-5.3.1/libsrc/csg/csgparser.cpp
      #../netgen-5.3.1/libsrc/csg/csgpkg.cpp
      ../netgen-5.3.1/libsrc/csg/curve2d.cpp
      ../netgen-5.3.1/libsrc/csg/edgeflw.cpp
      ../netgen-5.3.1/libsrc/csg/explicitcurve2d.cpp
      ../netgen-5.3.1/libsrc/csg/extrusion.cpp
      ../netgen-5.3.1/libsrc/csg/gencyl.cpp
      ../netgen-5.3.1/libsrc/csg/genmesh.cpp
      ../netgen-5.3.1/libsrc/csg/identify.cpp
      ../netgen-5.3.1/libsrc/csg/manifold.cpp
      ../netgen-5.3.1/libsrc/csg/meshsurf.cpp
      ../netgen-5.3.1/libsrc/csg/polyhedra.cpp
      ../netgen-5.3.1/libsrc/csg/revolution.cpp
      ../netgen-5.3.1/libsrc/csg/singularref.cpp
      ../netgen-5.3.1/libsrc/csg/solid.cpp
      ../netgen-5.3.1/libsrc/csg/specpoin.cpp
      ../netgen-5.3.1/libsrc/csg/spline3d.cpp
      ../netgen-5.3.1/libsrc/csg/surface.cpp
      ../netgen-5.3.1/libsrc/csg/triapprox.cpp
      #../netgen-5.3.1/libsrc/csg/vscsg.cpp
      ../netgen-5.3.1/libsrc/geom2d/genmesh2d.cpp
      ../netgen-5.3.1/libsrc/geom2d/geom2dmesh.cpp
      #../netgen-5.3.1/libsrc/geom2d/geom2dpkg.cpp
      ../netgen-5.3.1/libsrc/geom2d/geometry2d.cpp
      #../netgen-5.3.1/libsrc/geom2d/vsgeom2d.cpp
      ../netgen-5.3.1/libsrc/general/array.cpp
      ../netgen-5.3.1/libsrc/general/bitarray.cpp
      ../netgen-5.3.1/libsrc/general/dynamicmem.cpp
      ../netgen-5.3.1/libsrc/general/flags.cpp
      ../netgen-5.3.1/libsrc/general/gzstream.cpp
      ../netgen-5.3.1/libsrc/general/hashtabl.cpp
      ../netgen-5.3.1/libsrc/general/mpi_interface.cpp
      ../netgen-5.3.1/libsrc/general/mystring.cpp
      ../netgen-5.3.1/libsrc/general/ngexception.cpp
      ../netgen-5.3.1/libsrc/general/optmem.cpp
      ../netgen-5.3.1/libsrc/general/parthreads.cpp
      ../netgen-5.3.1/libsrc/general/profiler.cpp
      ../netgen-5.3.1/libsrc/general/seti.cpp
      ../netgen-5.3.1/libsrc/general/sort.cpp
      ../netgen-5.3.1/libsrc/general/spbita2d.cpp
      ../netgen-5.3.1/libsrc/general/symbolta.cpp
      ../netgen-5.3.1/libsrc/general/table.cpp
      ../netgen-5.3.1/libsrc/gprim/adtree.cpp
      ../netgen-5.3.1/libsrc/gprim/geom2d.cpp
      ../netgen-5.3.1/libsrc/gprim/geom3d.cpp
      ../netgen-5.3.1/libsrc/gprim/geomfuncs.cpp
      ../netgen-5.3.1/libsrc/gprim/geomtest3d.cpp
      ../netgen-5.3.1/libsrc/gprim/spline.cpp
      ../netgen-5.3.1/libsrc/gprim/splinegeometry.cpp
      ../netgen-5.3.1/libsrc/gprim/transform3d.cpp
      ../netgen-5.3.1/libsrc/linalg/bfgs.cpp
      ../netgen-5.3.1/libsrc/linalg/densemat.cpp
      ../netgen-5.3.1/libsrc/linalg/linopt.cpp
      ../netgen-5.3.1/libsrc/linalg/linsearch.cpp
      ../netgen-5.3.1/libsrc/linalg/polynomial.cpp
      ../netgen-5.3.1/libsrc/interface/nginterface.cpp
      ../netgen-5.3.1/libsrc/interface/nginterface_v2.cpp
      ../netgen-5.3.1/libsrc/interface/read_fnf_mesh.cpp
      ../netgen-5.3.1/libsrc/interface/readtetmesh.cpp
      ../netgen-5.3.1/libsrc/interface/readuser.cpp
      ../netgen-5.3.1/libsrc/interface/writeabaqus.cpp
      ../netgen-5.3.1/libsrc/interface/writediffpack.cpp
      ../netgen-5.3.1/libsrc/interface/writedolfin.cpp
      ../netgen-5.3.1/libsrc/interface/writeelmer.cpp
      ../netgen-5.3.1/libsrc/interface/writefeap.cpp
      ../netgen-5.3.1/libsrc/interface/writefluent.cpp
      ../netgen-5.3.1/libsrc/interface/writegmsh.cpp
      ../netgen-5.3.1/libsrc/interface/writegmsh2.cpp
      ../netgen-5.3.1/libsrc/interface/writejcm.cpp
      ../netgen-5.3.1/libsrc/interface/writeOpenFOAM15x.cpp
      ../netgen-5.3.1/libsrc/interface/writepermas.cpp
      ../netgen-5.3.1/libsrc/interface/writetecplot.cpp
      ../netgen-5.3.1/libsrc/interface/writetet.cpp
      ../netgen-5.3.1/libsrc/interface/writetochnog.cpp
      ../netgen-5.3.1/libsrc/interface/writeuser.cpp
      ../netgen-5.3.1/libsrc/interface/wuchemnitz.cpp
      ../netgen-5.3.1/libsrc/meshing/adfront2.cpp
      ../netgen-5.3.1/libsrc/meshing/adfront3.cpp
      ../netgen-5.3.1/libsrc/meshing/basegeom.cpp
      ../netgen-5.3.1/libsrc/meshing/bcfunctions.cpp
      ../netgen-5.3.1/libsrc/meshing/bisect.cpp
      ../netgen-5.3.1/libsrc/meshing/boundarylayer.cpp
      ../netgen-5.3.1/libsrc/meshing/clusters.cpp
      ../netgen-5.3.1/libsrc/meshing/curvedelems.cpp
      ../netgen-5.3.1/libsrc/meshing/delaunay.cpp
      ../netgen-5.3.1/libsrc/meshing/delaunay2d.cpp
      ../netgen-5.3.1/libsrc/meshing/geomsearch.cpp
      ../netgen-5.3.1/libsrc/meshing/global.cpp
      ../netgen-5.3.1/libsrc/meshing/hprefinement.cpp
      ../netgen-5.3.1/libsrc/meshing/improve2.cpp
      ../netgen-5.3.1/libsrc/meshing/improve2gen.cpp
      ../netgen-5.3.1/libsrc/meshing/improve3.cpp
      ../netgen-5.3.1/libsrc/meshing/localh.cpp
      ../netgen-5.3.1/libsrc/meshing/meshclass.cpp
      ../netgen-5.3.1/libsrc/meshing/meshfunc.cpp
      ../netgen-5.3.1/libsrc/meshing/meshfunc2d.cpp
      ../netgen-5.3.1/libsrc/meshing/meshing2.cpp
      ../netgen-5.3.1/libsrc/meshing/meshing3.cpp
      ../netgen-5.3.1/libsrc/meshing/meshtool.cpp
      ../netgen-5.3.1/libsrc/meshing/meshtype.cpp
      ../netgen-5.3.1/libsrc/meshing/msghandler.cpp
      ../netgen-5.3.1/libsrc/meshing/netrule2.cpp
      ../netgen-5.3.1/libsrc/meshing/netrule3.cpp
      ../netgen-5.3.1/libsrc/meshing/parallelmesh.cpp
      ../netgen-5.3.1/libsrc/meshing/paralleltop.cpp
      ../netgen-5.3.1/libsrc/meshing/parser2.cpp
      ../netgen-5.3.1/libsrc/meshing/parser3.cpp
      ../netgen-5.3.1/libsrc/meshing/prism2rls.cpp
      ../netgen-5.3.1/libsrc/meshing/pyramid2rls.cpp
      ../netgen-5.3.1/libsrc/meshing/pyramidrls.cpp
      ../netgen-5.3.1/libsrc/meshing/quadrls.cpp
      ../netgen-5.3.1/libsrc/meshing/refine.cpp
      ../netgen-5.3.1/libsrc/meshing/ruler2.cpp
      ../netgen-5.3.1/libsrc/meshing/ruler3.cpp
      ../netgen-5.3.1/libsrc/meshing/secondorder.cpp
      ../netgen-5.3.1/libsrc/meshing/smoothing2.5.cpp
      ../netgen-5.3.1/libsrc/meshing/smoothing2.cpp
      ../netgen-5.3.1/libsrc/meshing/smoothing3.cpp
      ../netgen-5.3.1/libsrc/meshing/specials.cpp
      ../netgen-5.3.1/libsrc/meshing/tetrarls.cpp
      ../netgen-5.3.1/libsrc/meshing/topology.cpp
      ../netgen-5.3.1/libsrc/meshing/triarls.cpp
      ../netgen-5.3.1/libsrc/meshing/validate.cpp
      ../netgen-5.3.1/libsrc/meshing/zrefine.cpp
      ../netgen-5.3.1/libsrc/occ/occconstruction.cpp
      ../netgen-5.3.1/libsrc/occ/occgenmesh.cpp
      ../netgen-5.3.1/libsrc/occ/occgeom.cpp
      ../netgen-5.3.1/libsrc/occ/occmeshsurf.cpp
      #../netgen-5.3.1/libsrc/occ/occpkg.cpp
      #../netgen-5.3.1/libsrc/occ/vsocc.cpp
      ../netgen-5.3.1/libsrc/occ/Partition_Inter2d.cxx
      ../netgen-5.3.1/libsrc/occ/Partition_Inter3d.cxx
      ../netgen-5.3.1/libsrc/occ/Partition_Loop.cxx
      ../netgen-5.3.1/libsrc/occ/Partition_Loop2d.cxx
      ../netgen-5.3.1/libsrc/occ/Partition_Loop3d.cxx
      ../netgen-5.3.1/libsrc/occ/Partition_Spliter.cxx
      ../netgen-5.3.1/libsrc/stlgeom/meshstlsurface.cpp
      ../netgen-5.3.1/libsrc/stlgeom/stlgeom.cpp
      ../netgen-5.3.1/libsrc/stlgeom/stlgeomchart.cpp
      ../netgen-5.3.1/libsrc/stlgeom/stlgeommesh.cpp
      ../netgen-5.3.1/libsrc/stlgeom/stlline.cpp
      #../netgen-5.3.1/libsrc/stlgeom/stlpkg.cpp
      ../netgen-5.3.1/libsrc/stlgeom/stltool.cpp
      ../netgen-5.3.1/libsrc/stlgeom/stltopology.cpp
      #../netgen-5.3.1/libsrc/stlgeom/vsstl.cpp
      ## NGLIB
      ../netgen-5.3.1/nglib/nglib.cpp)

    add_definitions(-DOCCGEOMETRY -DPACKAGE_VERSION="5.3.1")
else()
    set(NETGEN_SRCS "")
endif()

if(MONOLITHIC_BUILD)
    set(SMESH_LIBRARIES SMESH)
    add_library(SMESH
        ${NETGEN_SRCS}
        ${NETGENPlugin_SRCS}
        ${Driver_SRCS}
        ${DriverDAT_SRCS}
        ${DriverSTL_SRCS}
        ${DriverUNV_SRCS}
        ${SMDS_SRCS}
        ${SMESH_SRCS}
        ${SMESHDS_SRCS}
        ${StdMeshers_SRCS})
    target_link_libraries(SMESH ${OCE_LIBRARIES} TKSTL TKXCAF TKXDESTEP TKXDEIGES ${ZLIB_LIBRARIES})
    if(NETGEN_PLUGIN AND NETGEN_INTERNAL)
        target_link_libraries(SMESH ${ZLIB_LIBRARIES})
    endif()
    set_target_properties(SMESH PROPERTIES
        SOVERSION ${${PROJECT_NAME}_ABI_SOVERSION}
        VERSION ${${PROJECT_NAME}_ABI_VERSION})
else()
    # Build separate libraries
    set(SMESH_LIBRARIES Driver;DriverDAT;DriverSTL;DriverUNV;SMDS;SMESH;SMESHDS;StdMeshers)
    if(NETGEN_PLUGIN)
        list(APPEND SMESH_LIBRARIES NETGENPlugin)
        if(NETGEN_INTERNAL)
            list(APPEND ${NETGENPlugin_SRCS} ${NETGEN_SRCS})
        endif()
    endif()
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
    if(NETGEN_PLUGIN AND NOT NETGEN_INTERNAL)
        target_link_libraries(NETGENPlugin ${NETGEN_LIBRARY})
   endif()
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
