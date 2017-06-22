// Copyright (C) 2007-2016  CEA/DEN, EDF R&D, OPEN CASCADE
//
// Copyright (C) 2003-2007  OPEN CASCADE, EADS/CCR, LIP6, CEA/DEN,
// CEDRAT, EDF R&D, LEG, PRINCIPIA R&D, BUREAU VERITAS
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
// See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//
// File      : DriverCGNS_Write.cxx
// Created   : Fri Aug  5 17:43:54 2011
// Author    : Edward AGAPOV (eap)

#include "DriverCGNS_Write.hxx"

#include "SMDS_MeshNode.hxx"
#include "SMDS_VolumeTool.hxx"
#include "SMESHDS_GroupBase.hxx"
#include "SMESHDS_Mesh.hxx"
#include "SMESH_Comment.hxx"

#include <limits>
#include <cgnslib.h>

#if CGNS_VERSION < 3100
# define cgsize_t int
#endif

using namespace std;

namespace
{
  //================================================================================
  /*!
   * \brief Return interlace and type of CGNS element for the given SMDSAbs_EntityType
   */
  //================================================================================

  const int* getInterlaceAndType( const SMDSAbs_EntityType      smType,
                                  CGNS_ENUMT( ElementType_t ) & cgType )
  {
    static vector< const int* >                 interlaces;
    static vector< CGNS_ENUMT( ElementType_t )> cgTypes; 
    if ( interlaces.empty() )
    {
      interlaces.resize( SMDSEntity_Last, 0 );
      cgTypes.resize( SMDSEntity_Last, CGNS_ENUMV( ElementTypeNull ));
      {
        static int ids[] = {0};
        interlaces[SMDSEntity_0D] = ids;
        cgTypes   [SMDSEntity_0D] = CGNS_ENUMV( NODE );
      }
      {
        static int ids[] = { 0, 1 };
        interlaces[SMDSEntity_Edge] = ids;
        cgTypes   [SMDSEntity_Edge] = CGNS_ENUMV( BAR_2 );
      }
      {
        static int ids[] = { 0, 1, 2 };
        interlaces[SMDSEntity_Quad_Edge] = ids;
        cgTypes   [SMDSEntity_Quad_Edge] = CGNS_ENUMV( BAR_3 );
      }
      {
        static int ids[] = { 0, 2, 1 };
        interlaces[SMDSEntity_Triangle] = ids;
        cgTypes   [SMDSEntity_Triangle] = CGNS_ENUMV( TRI_3 );
      }
      {
        static int ids[] = { 0, 2, 1, 5, 4, 3 };
        interlaces[SMDSEntity_Quad_Triangle] = ids;
        cgTypes   [SMDSEntity_Quad_Triangle] = CGNS_ENUMV( TRI_6 );
        interlaces[SMDSEntity_BiQuad_Triangle] = ids;
        cgTypes   [SMDSEntity_BiQuad_Triangle] = CGNS_ENUMV( TRI_6 );
      }
      {
        static int ids[] = { 0, 3, 2, 1 };
        interlaces[SMDSEntity_Quadrangle] = ids;
        cgTypes   [SMDSEntity_Quadrangle] = CGNS_ENUMV( QUAD_4 );
      }
      {
        static int ids[] = { 0,3,2,1,7,6,5,4 };
        interlaces[SMDSEntity_Quad_Quadrangle] = ids;
        cgTypes   [SMDSEntity_Quad_Quadrangle] = CGNS_ENUMV( QUAD_8 );
      }
      {
        static int ids[] = { 0,3,2,1,7,6,5,4,8 };
        interlaces[SMDSEntity_BiQuad_Quadrangle] = ids;
        cgTypes   [SMDSEntity_BiQuad_Quadrangle] = CGNS_ENUMV( QUAD_9 );
      }
      {
        static int ids[] = { 0, 2, 1, 3 };
        interlaces[SMDSEntity_Tetra] = ids;
        cgTypes   [SMDSEntity_Tetra] = CGNS_ENUMV( TETRA_4 );
      }
      {
        static int ids[] = { 0,2,1,3,6,5,4,7,9,8 };
        interlaces[SMDSEntity_Quad_Tetra] = ids;
        cgTypes   [SMDSEntity_Quad_Tetra] = CGNS_ENUMV( TETRA_10 );
      }
      {
        static int ids[] = { 0,3,2,1,4 };
        interlaces[SMDSEntity_Pyramid] = ids;
        cgTypes   [SMDSEntity_Pyramid] = CGNS_ENUMV( PYRA_5 );
      }
      {
        static int ids[] = { 0,3,2,1,4,8,7,6,5,9,12,11,10 };
        interlaces[SMDSEntity_Quad_Pyramid] = ids;
        cgTypes   [SMDSEntity_Quad_Pyramid] = CGNS_ENUMV( PYRA_13 );
      }
      {
        static int ids[] = { 0,2,1,3,5,4 };
        interlaces[SMDSEntity_Penta] = ids;
        cgTypes   [SMDSEntity_Penta] = CGNS_ENUMV( PENTA_6 );
      }
      {
        static int ids[] = { 0,2,1,3,5,4,8,7,6,9,11,10,14,13,12 };
        interlaces[SMDSEntity_Quad_Penta] = ids;
        cgTypes   [SMDSEntity_Quad_Penta] = CGNS_ENUMV( PENTA_15 );
      }
      {
        static int ids[] = { 0,3,2,1,4,7,6,5 };
        interlaces[SMDSEntity_Hexa] = ids;
        cgTypes   [SMDSEntity_Hexa] = CGNS_ENUMV( HEXA_8 );
      }
      {
        static int ids[] = { 0,3,2,1,4,7,6,5,11,10,9,8,12,15,14,13,19,18,17,16 };
        interlaces[SMDSEntity_Quad_Hexa] = ids;
        cgTypes   [SMDSEntity_Quad_Hexa] = CGNS_ENUMV( HEXA_20 );
      }
      {
        static int ids[] = { 0,3,2,1,4,7,6,5,11,10,9,8,12,15,14,13,19,18,17,16,
                             20, 24,23,22,21, 25, 26};
        interlaces[SMDSEntity_TriQuad_Hexa] = ids;
        cgTypes   [SMDSEntity_TriQuad_Hexa] = CGNS_ENUMV( HEXA_27 );
      }
      {
        cgTypes[SMDSEntity_Polygon]         = CGNS_ENUMV( NGON_n );
        cgTypes[SMDSEntity_Quad_Polygon]    = CGNS_ENUMV( NGON_n );
        cgTypes[SMDSEntity_Polyhedra]       = CGNS_ENUMV( NFACE_n );
        cgTypes[SMDSEntity_Hexagonal_Prism] = CGNS_ENUMV( NFACE_n );
      }
    }
    cgType  = cgTypes[ smType ];
    return interlaces[ smType ];
  }

  //================================================================================
  /*!
   * \brief Cut off type of boundary condition from the group name
   */
  //================================================================================

  CGNS_ENUMT( BCType_t ) getBCType( string& groupName )
  {
    CGNS_ENUMT( BCType_t ) bcType = CGNS_ENUMV( BCGeneral ); // default type

    // boundary condition type starts from "BC"
    size_t bcBeg = groupName.find("BC");
    if ( bcBeg != string::npos )
    {
      for ( int t = 0; t < NofValidBCTypes; ++t )
      {
        CGNS_ENUMT( BCType_t ) type = CGNS_ENUMT( BCType_t )( t );
        string typeName = cg_BCTypeName( type );
        if ( typeName == &groupName[0] + bcBeg )
        {
          bcType = type;
          while ( bcBeg > 0 && isspace( bcBeg-1 ))
            --bcBeg;
          if ( bcBeg == 0 )
            groupName = "Group";
          else
            groupName = groupName.substr( 0, bcBeg-1 );
        }
      }
    }
    return bcType;
  }

  //================================================================================
  /*!
   * \brief Sortable face of a polyhedron
   */
  struct TPolyhedFace
  {
    int _id; // id of NGON_n
    vector< int > _nodes; // lowest node IDs used for sorting

    TPolyhedFace( const SMDS_MeshNode** nodes, const int nbNodes, int ID):_id(ID)
    {
      set< int > ids;
      for ( int i = 0; i < nbNodes; ++i )
        ids.insert( nodes[i]->GetID() );

      _nodes.resize( 3 ); // std::min( nbNodes, 4 )); hope 3 nodes is enough
      set< int >::iterator idIt = ids.begin();
      for ( size_t j = 0; j < _nodes.size(); ++j, ++idIt )
        _nodes[j] = *idIt;
    }
    bool operator< (const TPolyhedFace& o ) const
    {
      return _nodes < o._nodes;
    }
  };
  //================================================================================
  /*!
   * \brief Return CGNS id of an element
   */
  //================================================================================

  cgsize_t cgnsID( const SMDS_MeshElement*                         elem,
                   const map< const SMDS_MeshElement*, cgsize_t >& elem2cgID )
  {
    map< const SMDS_MeshElement*, cgsize_t >::const_iterator e2id = elem2cgID.find( elem );
    return ( e2id == elem2cgID.end() ? elem->GetID() : e2id->second );
  }

} // namespace

//================================================================================
/*!
 * \brief Write the mesh into the CGNS file
 */
//================================================================================

Driver_Mesh::Status DriverCGNS_Write::Perform()
{
  myErrorMessages.clear();

  if ( !myMesh || myMesh->GetMeshInfo().NbElements() < 1 )
    return addMessage( !myMesh ? "NULL mesh" : "Empty mesh (no elements)", /*fatal = */true );

  // open the file
  if ( cg_open(myFile.c_str(), CG_MODE_MODIFY, &_fn) != CG_OK &&
       cg_open(myFile.c_str(), CG_MODE_WRITE,  &_fn) != CG_OK )
    return addMessage( cg_get_error(), /*fatal = */true );

  // create a Base
  // --------------

  const int spaceDim = 3;
  int        meshDim = 1;
  if ( myMesh->NbFaces()   > 0 ) meshDim = 2;
  if ( myMesh->NbVolumes() > 0 ) meshDim = 3;

  if ( myMeshName.empty() )
  {
    int nbases = 0;
    if ( cg_nbases( _fn, &nbases) == CG_OK )
      myMeshName = ( SMESH_Comment("Base_") << nbases+1 );
    else
      myMeshName = "Base_0";
  }
  int iBase;
  if ( cg_base_write( _fn, myMeshName.c_str(), meshDim, spaceDim, &iBase ))
    return addMessage( cg_get_error(), /*fatal = */true );

  // create a Zone
  // --------------

  int nbCells = myMesh->NbEdges();
  if ( meshDim == 3 )
    nbCells = myMesh->NbVolumes();
  else if ( meshDim == 2 )
    nbCells = myMesh->NbFaces();

  cgsize_t size[9] = { myMesh->NbNodes(), nbCells, /*NBoundVertex=*/0, 0,0,0,0,0,0 };
  int iZone;
  if ( cg_zone_write( _fn, iBase, "SMESH_Mesh", size,
                      CGNS_ENUMV( Unstructured ), &iZone) != CG_OK )
    return addMessage( cg_get_error(), /*fatal = */true );

  // Map to store only elements whose an SMDS ID differs from a CGNS one
  typedef map< const SMDS_MeshElement*, cgsize_t > TElem2cgIDMap;
  vector< TElem2cgIDMap > elem2cgIDByEntity( SMDSEntity_Last );
  TElem2cgIDMap::iterator elem2cgIDIter;

  TElem2cgIDMap & n2cgID = elem2cgIDByEntity[ SMDSEntity_Node ];

  // Write nodes
  // ------------
  {
    vector< double > coords( myMesh->NbNodes() );
    int iC;
    // X
    SMDS_NodeIteratorPtr nIt = myMesh->nodesIterator( /*idInceasingOrder=*/true );
    for ( int i = 0; nIt->more(); ++i ) coords[i] = nIt->next()->X();
    if ( cg_coord_write( _fn, iBase, iZone, CGNS_ENUMV(RealDouble),
                          "CoordinateX", &coords[0], &iC) != CG_OK )
      return addMessage( cg_get_error(), /*fatal = */true );
    // Y
    nIt = myMesh->nodesIterator( /*idInceasingOrder=*/true );
    for ( int i = 0; nIt->more(); ++i ) coords[i] = nIt->next()->Y();
    if ( cg_coord_write( _fn, iBase, iZone, CGNS_ENUMV(RealDouble),
                          "CoordinateY", &coords[0], &iC) != CG_OK )
      return addMessage( cg_get_error(), /*fatal = */true );
    // Z
    nIt = myMesh->nodesIterator( /*idInceasingOrder=*/true );
    for ( int i = 0; nIt->more(); ++i ) coords[i] = nIt->next()->Z();
    if ( cg_coord_write( _fn, iBase, iZone, CGNS_ENUMV(RealDouble),
                          "CoordinateZ", &coords[0], &iC) != CG_OK )
      return addMessage( cg_get_error(), /*fatal = */true );

    // store CGNS ids of nodes
    nIt = myMesh->nodesIterator( /*idInceasingOrder=*/true );
    for ( int i = 0; nIt->more(); ++i )
    {
      const SMDS_MeshElement* n = nIt->next();
      if ( n->GetID() != i+1 )
        n2cgID.insert( n2cgID.end(), make_pair( n, i+1 ));
    }
  }
  // Write elements
  // ---------------
  
  cgsize_t cgID = 1, startID;

  // write into a section all successive elements of one geom type
  int iSec;
  vector< cgsize_t > elemData;
  SMDS_ElemIteratorPtr  elemIt = myMesh->elementsIterator();
  const SMDS_MeshElement* elem = elemIt->next();
  while ( elem )
  {
    const SMDSAbs_EntityType elemType = elem->GetEntityType();
    CGNS_ENUMT( ElementType_t ) cgType;
    const int* interlace = getInterlaceAndType( elemType, cgType );

    TElem2cgIDMap & elem2cgID = elem2cgIDByEntity[ elemType ];

    elemData.clear();
    startID = cgID;

    if ( interlace ) // STANDARD elements
    {
      int cgnsNbNodes; // get nb nodes by element type, that can be less that elem->NbNodes()
      cg_npe( cgType, &cgnsNbNodes );
      do
      {
        for ( int i = 0; i < cgnsNbNodes; ++i )
          elemData.push_back( cgnsID( elem->GetNode( interlace[i] ), n2cgID ));
        if ( elem->GetID() != cgID )
          elem2cgID.insert( elem2cgID.end(), make_pair( elem, cgID ));
        ++cgID;
        elem = elemIt->more() ? elemIt->next() : 0;
      }
      while ( elem && elem->GetEntityType() == elemType );
    }
    else if ( elemType == SMDSEntity_Polygon ) // POLYGONS
      do
      {
        elemData.push_back( elem->NbNodes() );
        for ( int i = 0, nb = elem->NbNodes(); i < nb; ++i )
          elemData.push_back( cgnsID( elem->GetNode(i), n2cgID ));
        if ( elem->GetID() != cgID )
          elem2cgID.insert( elem2cgID.end(), make_pair( elem, cgID ));
        ++cgID;
        elem = elemIt->more() ? elemIt->next() : 0;
      }
      while ( elem && elem->GetEntityType() == elemType );

    else if ( elemType == SMDSEntity_Quad_Polygon ) // QUADRATIC POLYGONS
      do // write as linear NGON_n
      {
        elemData.push_back( elem->NbNodes() );
        interlace = & SMDS_MeshCell::interlacedSmdsOrder( SMDSEntity_Quad_Polygon,
                                                          elem->NbNodes() )[0];
        for ( int i = 0, nb = elem->NbNodes(); i < nb; ++i )
          elemData.push_back( cgnsID( elem->GetNode( interlace[i] ), n2cgID ));
        if ( elem->GetID() != cgID )
          elem2cgID.insert( elem2cgID.end(), make_pair( elem, cgID ));
        ++cgID;
        elem = elemIt->more() ? elemIt->next() : 0;
      }
      while ( elem && elem->GetEntityType() == elemType );

    else if ( elemType == SMDSEntity_Polyhedra ||
              elemType == SMDSEntity_Hexagonal_Prism) // POLYHEDRA
    {
      // to save polyhedrons after all
      const SMDS_MeshInfo& meshInfo = myMesh->GetMeshInfo();
      if ( meshInfo.NbPolyhedrons() == meshInfo.NbElements() - cgID + 1 )
        break; // only polyhedrons remain
      while ( elem && elem->GetEntityType() == elemType )
        elem = elemIt->more() ? elemIt->next() : 0;
      continue;
    }

    SMESH_Comment sectionName( cg_ElementTypeName( cgType ));
    sectionName << " " << startID << " - " << cgID-1;

    if ( cg_section_write(_fn, iBase, iZone, sectionName.c_str(), cgType, startID,
                          cgID-1, /*nbndry=*/0, &elemData[0], &iSec) != CG_OK )
      return addMessage( cg_get_error(), /*fatal = */true );
  }
  // Write polyhedral volumes
  // -------------------------

  if ( myMesh->GetMeshInfo().NbElements() > cgID-1 ) // polyhedra or hexagonal prisms remain
  {
    // the polyhedron (NFACE_n) is described as a set of signed face IDs,
    // so first we are to write all polygones (NGON_n) bounding polyhedrons

    vector< cgsize_t > faceData;
    set< TPolyhedFace > faces;
    set< TPolyhedFace >::iterator faceInSet;
    vector<const SMDS_MeshNode *> faceNodesVec;
    int nbPolygones = 0, faceID;

    SMDS_VolumeTool vol;

    elemData.clear();

    int nbPolyhTreated = 0;

    TElem2cgIDMap * elem2cgID = 0;
    TElem2cgIDMap & n2cgID    = elem2cgIDByEntity[ SMDSEntity_Node ];

    SMDS_ElemIteratorPtr elemIt = myMesh->elementsIterator();
    while ( elemIt->more() )
    {
      elem = elemIt->next();
      SMDSAbs_EntityType type = elem->GetEntityType();
      if ( type == SMDSEntity_Polyhedra ||
           type == SMDSEntity_Hexagonal_Prism )
      {
        ++nbPolyhTreated;
        vol.Set( elem );
        vol.SetExternalNormal();
        const int nbFaces = vol.NbFaces();
        elemData.push_back( nbFaces );
        for ( int iF = 0; iF < nbFaces; ++iF )
        {
          const int nbNodes = vol.NbFaceNodes( iF );
          const SMDS_MeshNode** faceNodes = vol.GetFaceNodes( iF );
          faceNodesVec.assign( faceNodes, faceNodes + nbNodes );
          if (( elem = myMesh->FindElement( faceNodesVec, SMDSAbs_Face, /*noMedium=*/false)))
          {
            // a face of the polyhedron is present in the mesh
            faceID = cgnsID( elem, elem2cgIDByEntity[ elem->GetEntityType() ]);
          }
          else if ( vol.IsFreeFace( iF ))
          {
            // the face is not shared by volumes
            faceID = cgID++;
            ++nbPolygones;
            faceData.push_back( nbNodes );
            for ( int i = 0; i < nbNodes; ++i )
              faceData.push_back( cgnsID( faceNodes[i], n2cgID ));
          }
          else
          {
            TPolyhedFace face( faceNodes, nbNodes, cgID );
            faceInSet = faces.insert( faces.end(), face );
            if ( faceInSet->_id == cgID ) // the face encounters for the 1st time
            {
              faceID = cgID++;
              ++nbPolygones;
              faceData.push_back( nbNodes );
              for ( int i = 0; i < nbNodes; ++i )
                faceData.push_back( cgnsID( faceNodes[i], n2cgID ));
            }
            else
            {
              // the face encounters for the 2nd time; we hope it won't encounter once more,
              // for that we can erase it from the set of faces
              faceID = -faceInSet->_id;
              faces.erase( faceInSet );
            }
          }
          elemData.push_back( faceID );
        }
      }
    }

    if ( nbPolygones > 0 )
    {
      if ( cg_section_write(_fn, iBase, iZone, "Faces of Polyhedrons",
                             CGNS_ENUMV( NGON_n ), cgID - nbPolygones, cgID-1,
                             /*nbndry=*/0, &faceData[0], &iSec) != CG_OK )
        return addMessage( cg_get_error(), /*fatal = */true );
    }
    
    if ( cg_section_write(_fn, iBase, iZone, "Polyhedrons", 
                             CGNS_ENUMV( NFACE_n ), cgID, cgID+nbPolyhTreated-1,
                           /*nbndry=*/0, &elemData[0], &iSec) != CG_OK )
      return addMessage( cg_get_error(), /*fatal = */true );

    if ( !myMesh->GetGroups().empty() )
    {
      // store CGNS ids of polyhedrons
      elem2cgID = &elem2cgIDByEntity[ SMDSEntity_Polyhedra ];
      elemIt = myMesh->elementsIterator();
      while ( elemIt->more() )
      {
        elem = elemIt->next();
        if ( elem->GetEntityType() == SMDSEntity_Polyhedra )
        {
          if ( elem->GetID() != cgID )
            elem2cgID->insert( elem2cgID->end(), make_pair( elem, cgID ));
          ++cgID;
        }
      }
    }
  } // write polyhedral volumes


  // Write groups as boundary conditions
  // ------------------------------------

  const set<SMESHDS_GroupBase*>& groups = myMesh->GetGroups();
  set<SMESHDS_GroupBase*>::const_iterator grpIt = groups.begin();
  set< string > groupNames; groupNames.insert(""); // to avoid duplicated and empty names
  for ( ; grpIt != groups.end(); ++grpIt )
  {
    const SMESHDS_GroupBase* group = *grpIt;

    // write BC location (default is Vertex)
    CGNS_ENUMT( GridLocation_t ) location = CGNS_ENUMV( Vertex );
    if ( group->GetType() != SMDSAbs_Node )
    {
      switch ( meshDim ) {
      case 3:
        switch ( group->GetType() ) {
        case SMDSAbs_Volume: location = CGNS_ENUMV( FaceCenter ); break; // !!!
        case SMDSAbs_Face:   location = CGNS_ENUMV( FaceCenter ); break; // OK
        case SMDSAbs_Edge:   location = CGNS_ENUMV( EdgeCenter ); break; // OK
        default:;
        }
        break;
      case 2:
        switch ( group->GetType() ) {
        case SMDSAbs_Face: location = CGNS_ENUMV( FaceCenter ); break; // ???
        case SMDSAbs_Edge: location = CGNS_ENUMV( EdgeCenter ); break; // OK
        default:;
        }
        break;
      case 1:
        location = CGNS_ENUMV( EdgeCenter ); break; // ???
        break;
      }
    }

    // try to extract type of boundary condition from the group name
    string name = group->GetStoreName();
    CGNS_ENUMT( BCType_t ) bcType = getBCType( name );
    while ( !groupNames.insert( name ).second )
      name = (SMESH_Comment( "Group_") << groupNames.size());

    // write IDs of elements
    vector< cgsize_t > pnts;
    pnts.reserve( group->Extent() );
    SMDS_ElemIteratorPtr elemIt = group->GetElements();
    while ( elemIt->more() )
    {
      const SMDS_MeshElement* elem = elemIt->next();
      pnts.push_back( cgnsID( elem, elem2cgIDByEntity[ elem->GetEntityType() ]));
    }
    int iBC;
    if ( cg_boco_write( _fn, iBase, iZone, name.c_str(), bcType,
                        CGNS_ENUMV( PointList ), pnts.size(), &pnts[0], &iBC) != CG_OK )
      return addMessage( cg_get_error(), /*fatal = */true);

    // write BC location
    if ( location != CGNS_ENUMV( Vertex ))
    {
      if ( cg_boco_gridlocation_write( _fn, iBase, iZone, iBC, location) != CG_OK )
        return addMessage( cg_get_error(), /*fatal = */false);
    }
  }
  return DRS_OK;
}

//================================================================================
/*!
 * \brief Constructor
 */
//================================================================================

DriverCGNS_Write::DriverCGNS_Write(): _fn(0)
{
}

//================================================================================
/*!
 * \brief Close the cgns file at destruction
 */
//================================================================================

DriverCGNS_Write::~DriverCGNS_Write()
{
  if ( _fn > 0 )
    cg_close( _fn );
}
