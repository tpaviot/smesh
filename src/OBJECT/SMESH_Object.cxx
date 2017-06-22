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

//  SMESH OBJECT : interactive object for SMESH visualization
//  File   : SMESH_Grid.cxx
//  Author : Nicolas REJNERI
//  Module : SMESH
//
#include "SMESH_ObjectDef.h"
#include "SMESH_ActorUtils.h"

#include "SMDS_BallElement.hxx"
#include "SMDS_Mesh.hxx"
#include "SMDS_MeshCell.hxx"
#include "SMDS_PolyhedralVolumeOfNodes.hxx"
#include "SMESH_Actor.h"
#include "SMESH_ControlsDef.hxx"

#include <SalomeApp_Application.h>
#include <VTKViewer_ExtractUnstructuredGrid.h>
#include <VTKViewer_CellLocationsArray.h>

#include CORBA_SERVER_HEADER(SMESH_Gen)
#include CORBA_SERVER_HEADER(SALOME_Exception)

#include <vtkCell.h>
#include <vtkIdList.h>
#include <vtkCellArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkCellData.h>
#include <vtkUnstructuredGrid.h>

#include <memory>
#include <sstream>      
#include <stdexcept>
#include <set>

#include "utilities.h"

using namespace std;

#ifndef EXCEPTION
#define EXCEPTION(TYPE, MSG) {\
  std::ostringstream aStream;\
  aStream<<__FILE__<<"["<<__LINE__<<"]::"<<MSG;\
  throw TYPE(aStream.str());\
}
#endif

#ifdef _DEBUG_
static int MYDEBUG = 0;
static int MYDEBUGWITHFILES = 0;//1;
#else
static int MYDEBUG = 0;
static int MYDEBUGWITHFILES = 0;
#endif


/*
  Class       : SMESH_VisualObjDef
  Description : Base class for all mesh objects to be visuilised
*/

//=================================================================================
// function : getCellType
// purpose  : Get type of VTK cell
//=================================================================================
// static inline vtkIdType getCellType( const SMDSAbs_ElementType theType,
//                                      const bool                thePoly,
//                                      const int                 theNbNodes )
// {
//   switch( theType )
//   {
//     case SMDSAbs_0DElement:         return VTK_VERTEX;

//     case SMDSAbs_Ball:              return VTK_POLY_VERTEX;

//     case SMDSAbs_Edge: 
//       if( theNbNodes == 2 )         return VTK_LINE;
//       else if ( theNbNodes == 3 )   return VTK_QUADRATIC_EDGE;
//       else return VTK_EMPTY_CELL;

//     case SMDSAbs_Face  :
//       if (thePoly && theNbNodes>2 ) return VTK_POLYGON;
//       else if ( theNbNodes == 3 )   return VTK_TRIANGLE;
//       else if ( theNbNodes == 4 )   return VTK_QUAD;
//       else if ( theNbNodes == 6 )   return VTK_QUADRATIC_TRIANGLE;
//       else if ( theNbNodes == 8 )   return VTK_QUADRATIC_QUAD;
//       else if ( theNbNodes == 9 )   return VTK_BIQUADRATIC_QUAD;
//       else if ( theNbNodes == 7 )   return VTK_BIQUADRATIC_TRIANGLE;
//       else return VTK_EMPTY_CELL;
      
//     case SMDSAbs_Volume:
//       if (thePoly && theNbNodes>3 ) return VTK_POLYHEDRON; //VTK_CONVEX_POINT_SET;
//       else if ( theNbNodes == 4 )   return VTK_TETRA;
//       else if ( theNbNodes == 5 )   return VTK_PYRAMID;
//       else if ( theNbNodes == 6 )   return VTK_WEDGE;
//       else if ( theNbNodes == 8 )   return VTK_HEXAHEDRON;
//       else if ( theNbNodes == 12 )  return VTK_HEXAGONAL_PRISM;
//       else if ( theNbNodes == 10 )  return VTK_QUADRATIC_TETRA;
//       else if ( theNbNodes == 20 )  return VTK_QUADRATIC_HEXAHEDRON;
//       else if ( theNbNodes == 27 )  return VTK_TRIQUADRATIC_HEXAHEDRON;
//       else if ( theNbNodes == 15 )  return VTK_QUADRATIC_WEDGE;
//       else if ( theNbNodes == 13 )  return VTK_QUADRATIC_PYRAMID; //VTK_CONVEX_POINT_SET;
//       else return VTK_EMPTY_CELL;

//     default: return VTK_EMPTY_CELL;
//   }
// }

//=================================================================================
// functions : SMESH_VisualObjDef
// purpose   : Constructor
//=================================================================================
SMESH_VisualObjDef::SMESH_VisualObjDef()
{
  if ( MYDEBUG ) MESSAGE("-------------------------------SMESH_VisualObjDef::SMESH_VisualObjDef");
  myGrid = vtkUnstructuredGrid::New();
  myLocalGrid = false;
  ClearEntitiesFlags();
  SMESH::GetEntitiesFromObject(NULL);
}
SMESH_VisualObjDef::~SMESH_VisualObjDef()
{
  if ( MYDEBUG ) MESSAGE("--------------------------------SMESH_VisualObjDef::~SMESH_VisualObjDef");
  if ( MYDEBUG ) MESSAGE( "myGrid->GetReferenceCount() = " << myGrid->GetReferenceCount() );
  myGrid->Delete();
}

//=================================================================================
// functions : GetNodeObjId, GetNodeVTKId, GetElemObjId, GetElemVTKId
// purpose   : Methods for retrieving VTK IDs by SMDS IDs and  vice versa
//=================================================================================
vtkIdType SMESH_VisualObjDef::GetNodeObjId( int theVTKID )
{
  if (myLocalGrid)
  {
    TMapOfIds::const_iterator i = myVTK2SMDSNodes.find(theVTKID);
    return i == myVTK2SMDSNodes.end() ? -1 : i->second;
  }
  const SMDS_MeshNode* aNode = 0;
  if( this->GetMesh() )
    aNode = this->GetMesh()->FindNodeVtk( theVTKID );

  return aNode ? aNode->GetID() : -1;
}

vtkIdType SMESH_VisualObjDef::GetNodeVTKId( int theObjID )
{
  if (myLocalGrid)
  {
    TMapOfIds::const_iterator i = mySMDS2VTKNodes.find(theObjID);
    return i == mySMDS2VTKNodes.end() ? -1 : i->second;
  }

  const SMDS_MeshNode* aNode = 0;
  if( this->GetMesh() ) {
    aNode = this->GetMesh()->FindNode(theObjID);
  }
  return aNode ? aNode->getVtkId() : -1;
}

vtkIdType SMESH_VisualObjDef::GetElemObjId( int theVTKID )
{
  if (myLocalGrid)
  {
    TMapOfIds::const_iterator i = myVTK2SMDSElems.find(theVTKID);
    return i == myVTK2SMDSElems.end() ? -1 : i->second;
  }
  return this->GetMesh()->fromVtkToSmds(theVTKID);
}

vtkIdType SMESH_VisualObjDef::GetElemVTKId( int theObjID )
{
  if (myLocalGrid)
  {
    TMapOfIds::const_iterator i = mySMDS2VTKElems.find(theObjID);
    return i == mySMDS2VTKElems.end() ? -1 : i->second;
  }

  const SMDS_MeshElement* e = 0;
  if ( this->GetMesh() )
    e = this->GetMesh()->FindElement(theObjID);

  return e ? e->getVtkId() : -1;
}

//=================================================================================
// function : SMESH_VisualObjDef::createPoints
// purpose  : Create points from nodes
//=================================================================================
/*! fills a vtkPoints structure for a submesh.
 *  fills a std::list of SMDS_MeshElements*, then extract the points.
 *  fills also conversion id maps between SMDS and VTK.
 */
void SMESH_VisualObjDef::createPoints( vtkPoints* thePoints )
{
  if ( thePoints == 0 )
    return;

  TEntityList aNodes;
  vtkIdType nbNodes = GetEntities( SMDSAbs_Node, aNodes );
  thePoints->SetNumberOfPoints( nbNodes );

  int nbPoints = 0;

  TEntityList::const_iterator anIter;
  for ( anIter = aNodes.begin(); anIter != aNodes.end(); ++anIter )
  {
    const SMDS_MeshNode* aNode = ( const SMDS_MeshNode* )(*anIter);
    if ( aNode != 0 )
    {
      thePoints->SetPoint( nbPoints, aNode->X(), aNode->Y(), aNode->Z() );
      int anId = aNode->GetID();
      mySMDS2VTKNodes.insert( mySMDS2VTKNodes.end(), std::make_pair( anId, nbPoints ));
      myVTK2SMDSNodes.insert( myVTK2SMDSNodes.end(), std::make_pair( nbPoints, anId ));
      nbPoints++;
    }
  }

  if ( nbPoints != nbNodes )
    thePoints->SetNumberOfPoints( nbPoints );
}

//=================================================================================
// function : buildPrs
// purpose  : create VTK cells( fill unstructured grid )
//=================================================================================
void SMESH_VisualObjDef::buildPrs(bool buildGrid)
{
  if ( MYDEBUG ) MESSAGE("---------------------------SMESH_VisualObjDef::buildPrs " << buildGrid);
  if (buildGrid)
  {
    myLocalGrid = true;
    try
    {
      mySMDS2VTKNodes.clear();
      myVTK2SMDSNodes.clear();
      mySMDS2VTKElems.clear();
      myVTK2SMDSElems.clear();

      if ( IsNodePrs() )
        buildNodePrs();
      else
        buildElemPrs();
    }
    catch(...)
    {
      mySMDS2VTKNodes.clear();
      myVTK2SMDSNodes.clear();
      mySMDS2VTKElems.clear();
      myVTK2SMDSElems.clear();

      myGrid->SetPoints( 0 );
      myGrid->SetCells( 0, 0, 0, 0, 0 );
      throw;
    }
  }
  else
  {
    myLocalGrid = false;
    if (!GetMesh()->isCompacted())
    {
      NulData(); // detach from the SMDS grid to allow immediate memory de-allocation in compactMesh()
      if ( MYDEBUG ) MESSAGE("*** buildPrs ==> compactMesh!");
      GetMesh()->compactMesh();
    }
    vtkUnstructuredGrid *theGrid = GetMesh()->getGrid();
    updateEntitiesFlags();
    myGrid->ShallowCopy(theGrid);
    //MESSAGE(myGrid->GetReferenceCount());
    //MESSAGE( "Update - myGrid->GetNumberOfCells() = "<<myGrid->GetNumberOfCells() );
    //MESSAGE( "Update - myGrid->GetNumberOfPoints() = "<<myGrid->GetNumberOfPoints() );
    if( MYDEBUGWITHFILES ) {
      SMESH::WriteUnstructuredGrid( myGrid,"myPrs.vtu" );
    }
  }
}

//=================================================================================
// function : buildNodePrs
// purpose  : create VTK cells for nodes
//=================================================================================

void SMESH_VisualObjDef::buildNodePrs()
{
  // PAL16631: without swap, bad_alloc is not thrown but hung up and crash instead,
  // so check remaining memory size for safety
  SMDS_Mesh::CheckMemory(); // PAL16631
  vtkPoints* aPoints = vtkPoints::New();
  createPoints( aPoints );
  SMDS_Mesh::CheckMemory();
  myGrid->SetPoints( aPoints );
  aPoints->Delete();

  myGrid->SetCells( 0, 0, 0, 0, 0 );
}

//=================================================================================
// function : buildElemPrs
// purpose  : Create VTK cells for elements
//=================================================================================

namespace{
  typedef std::vector<const SMDS_MeshElement*> TConnect;

  int GetConnect(const SMDS_ElemIteratorPtr& theNodesIter, 
                 TConnect& theConnect)
  {
    theConnect.clear();
    for(; theNodesIter->more();)
      theConnect.push_back(theNodesIter->next());
    return theConnect.size();
  }
  
  inline 
  void SetId(vtkIdList *theIdList, 
             const SMESH_VisualObjDef::TMapOfIds& theSMDS2VTKNodes, 
             const TConnect& theConnect, 
             int thePosition,
             int theId)
  {
    theIdList->SetId(thePosition,theSMDS2VTKNodes.find(theConnect[theId]->GetID())->second);
  }

}


void SMESH_VisualObjDef::buildElemPrs()
{
  // Create points

  vtkPoints* aPoints = vtkPoints::New();
  createPoints( aPoints );
  myGrid->SetPoints( aPoints );
  aPoints->Delete();

  if ( MYDEBUG )
    MESSAGE("Update - myGrid->GetNumberOfPoints() = "<<myGrid->GetNumberOfPoints());

  // Calculate cells size

  const int nbTypes = 5;
  static SMDSAbs_ElementType aTypes[ nbTypes ] =
    { SMDSAbs_Edge, SMDSAbs_Face, SMDSAbs_Volume, SMDSAbs_Ball, SMDSAbs_0DElement };

  // get entity data
  map<SMDSAbs_ElementType,int> nbEnts;
  map<SMDSAbs_ElementType,TEntityList> anEnts;

  vtkIdType aNbCells = 0;

  for ( int i = 0; i < nbTypes; i++ )
  {
    nbEnts[ aTypes[ i ] ] = GetEntities( aTypes[ i ], anEnts[ aTypes[ i ] ] );
    aNbCells += nbEnts[ aTypes [ i ]];
  }
  // PAL16631: without swap, bad_alloc is not thrown but hung up and crash instead,
  // so check remaining memory size for safety
  SMDS_Mesh::CheckMemory(); // PAL16631

  vtkIdType aCellsSize =  2 * nbEnts[ SMDSAbs_0DElement ] + 3 * nbEnts[ SMDSAbs_Edge ];
  aCellsSize += 2 * nbEnts[ SMDSAbs_Ball ];

  for ( int i = 1; i <= 2; i++ ) // iterate through faces and volumes
  {
    if ( nbEnts[ aTypes[ i ] ] )
    {
      const TEntityList& aList = anEnts[ aTypes[ i ] ];
      TEntityList::const_iterator anIter;
      for ( anIter = aList.begin(); anIter != aList.end(); ++anIter ) {
        if((*anIter)->GetEntityType() != SMDSEntity_Polyhedra &&
           (*anIter)->GetEntityType() != SMDSEntity_Quad_Polyhedra) {
          aCellsSize += (*anIter)->NbNodes() + 1;
        } 
        // Special case for the VTK_POLYHEDRON:
        // itsinput cellArray is of special format.
        //  [nCellFaces, nFace0Pts, i, j, k, nFace1Pts, i, j, k, ...]   
        else {
          if( const SMDS_VtkVolume* ph = dynamic_cast<const SMDS_VtkVolume*>(*anIter) ) {
            int nbFaces = ph->NbFaces();
            aCellsSize += (1 + ph->NbFaces());
            for( int i = 1; i <= nbFaces; i++ ) {
              aCellsSize += ph->NbFaceNodes(i);
            }
          }
        }
      }
    }
  }
  if ( MYDEBUG )
    MESSAGE( "Update - aNbCells = "<<aNbCells<<"; aCellsSize = "<<aCellsSize );

  // Create cells

  vtkCellArray* aConnectivity = vtkCellArray::New();
  aConnectivity->Allocate( aCellsSize, 0 );

  SMDS_Mesh::CheckMemory(); // PAL16631

  vtkUnsignedCharArray* aCellTypesArray = vtkUnsignedCharArray::New();
  aCellTypesArray->SetNumberOfComponents( 1 );
  aCellTypesArray->Allocate( aNbCells * aCellTypesArray->GetNumberOfComponents() );

  SMDS_Mesh::CheckMemory(); // PAL16631

  vtkIdList *anIdList = vtkIdList::New();
  vtkIdType iElem = 0;

  TConnect aConnect;
  aConnect.reserve(VTK_CELL_SIZE);

  SMDS_Mesh::CheckMemory(); // PAL16631
  bool hasBalls = nbEnts[ SMDSAbs_Ball ] > 0;
  vtkDataArray* aScalars = 0;
  if(hasBalls) {
    aScalars = vtkDataArray::CreateDataArray(VTK_DOUBLE);
    aScalars->SetNumberOfComponents(1);
    aScalars->SetNumberOfTuples(aNbCells);
  }
  for ( int i = 0; i < nbTypes; i++ ) // iterate through all types of elements
  {
    if ( nbEnts[ aTypes[ i ] ] > 0 ) {

      const SMDSAbs_ElementType& aType = aTypes[ i ];
      const TEntityList& aList = anEnts[ aType ];
      TEntityList::const_iterator anIter;
      for ( anIter = aList.begin(); anIter != aList.end(); ++anIter )
      {
        const SMDS_MeshElement* anElem = *anIter;

        vtkIdType aNbNodes = anElem->NbNodes();
        anIdList->SetNumberOfIds( aNbNodes );
        const vtkIdType vtkElemType = SMDS_MeshCell::toVtkType( anElem->GetEntityType() );

        int anId = anElem->GetID();

        mySMDS2VTKElems.insert( mySMDS2VTKElems.end(), std::make_pair( anId, iElem ));
        myVTK2SMDSElems.insert( myVTK2SMDSElems.end(), std::make_pair( iElem, anId ));

        SMDS_ElemIteratorPtr aNodesIter = anElem->nodesIterator();
        {
          // Convertions connectivities from SMDS to VTK

          if (aType == SMDSAbs_Volume && anElem->IsPoly() && aNbNodes > 3) { // POLYEDRE
            anIdList->Reset();
            if ( const SMDS_VtkVolume* ph = dynamic_cast<const SMDS_VtkVolume*>(anElem) ) {
              int nbFaces = ph->NbFaces();
              anIdList->InsertNextId(nbFaces);
              for( int i = 1; i <= nbFaces; i++ ) {
                anIdList->InsertNextId(ph->NbFaceNodes(i));
                for(int j = 1; j <= ph->NbFaceNodes(i); j++) {
                  const SMDS_MeshNode* n = ph->GetFaceNode(i,j);
                  if(n) {
                    anIdList->InsertNextId(mySMDS2VTKNodes[n->GetID()]);
                  }
                }
              }
            }
          }
          else {
            const std::vector<int>& aConnectivities =
              SMDS_MeshCell::toVtkOrder( VTKCellType( vtkElemType ));
            if (aConnectivities.size() > 0) {
              aConnect.clear();
              GetConnect(aNodesIter,aConnect);
              for (vtkIdType aNodeId = 0; aNodeId < aNbNodes; aNodeId++)
                SetId(anIdList,mySMDS2VTKNodes,aConnect,aNodeId,aConnectivities[aNodeId]);
            }
            else {
              for( vtkIdType aNodeId = 0; aNodesIter->more(); aNodeId++ ){
                const SMDS_MeshElement* aNode = aNodesIter->next();
                anIdList->SetId( aNodeId, mySMDS2VTKNodes[aNode->GetID()] );
              }
            }
          }
        }
        vtkIdType aCurId = aConnectivity->InsertNextCell( anIdList );
        aCellTypesArray->InsertNextValue( vtkElemType );
        
        //Store diameters of the balls
        if(aScalars) {
          double aDiam = 0;
          if(aType == SMDSAbs_Ball) {
            if (const SMDS_BallElement* ball = dynamic_cast<const SMDS_BallElement*>(anElem) ) {
              aDiam = ball->GetDiameter();
            }
          }
          aScalars->SetTuple(aCurId,&aDiam);
        }

        iElem++;
      }
    }
    SMDS_Mesh::CheckMemory(); // PAL16631
  }

  // Insert cells in grid

  VTKViewer_CellLocationsArray* aCellLocationsArray = VTKViewer_CellLocationsArray::New();
  aCellLocationsArray->SetNumberOfComponents( 1 );
  aCellLocationsArray->SetNumberOfTuples( aNbCells );

  SMDS_Mesh::CheckMemory(); // PAL16631

  aConnectivity->InitTraversal();
  for( vtkIdType idType = 0, *pts, npts; aConnectivity->GetNextCell( npts, pts ); idType++ )
    aCellLocationsArray->SetValue( idType, aConnectivity->GetTraversalLocation( npts ) );

  myGrid->SetCells( aCellTypesArray, aCellLocationsArray,aConnectivity );
  myGrid->GetCellData()->SetScalars(aScalars);

  aCellLocationsArray->Delete();
  aCellTypesArray->Delete();
  aConnectivity->Delete();
  anIdList->Delete();

  SMDS_Mesh::CheckMemory(); // PAL16631
}

//=================================================================================
// function : GetEdgeNodes
// purpose  : Retrieve ids of nodes from edge of elements ( edge is numbered from 0 )
//=================================================================================
bool SMESH_VisualObjDef::GetEdgeNodes( const int theElemId,
                                       const int theEdgeNum,
                                       int&      theNodeId1,
                                       int&      theNodeId2 ) const
{
  const SMDS_Mesh* aMesh = GetMesh();
  if ( aMesh == 0 )
    return false;
    
  const SMDS_MeshElement* anElem = aMesh->FindElement( theElemId );
  if ( anElem == 0 )
    return false;
    
  int nbNodes = anElem->NbCornerNodes();

  if (( theEdgeNum < 0 || theEdgeNum > 3 ) ||
      ( nbNodes != 3 && nbNodes != 4 ) ||
      ( theEdgeNum >= nbNodes ))
    return false;

  theNodeId1 = anElem->GetNode(  theEdgeNum                 )->GetID();
  theNodeId2 = anElem->GetNode(( theEdgeNum + 1 ) % nbNodes )->GetID();

  return true;
}

vtkUnstructuredGrid* SMESH_VisualObjDef::GetUnstructuredGrid()
{
  if ( !myLocalGrid && !GetMesh()->isCompacted() )
  {
    NulData(); // detach from the SMDS grid to allow immediate memory de-allocation in compactMesh()
    GetMesh()->compactMesh();
    updateEntitiesFlags();
    vtkUnstructuredGrid *theGrid = GetMesh()->getGrid();
    myGrid->ShallowCopy(theGrid);
  }
  return myGrid;
}


//=================================================================================
// function : IsValid
// purpose  : Return true if there are some entities
//=================================================================================
bool SMESH_VisualObjDef::IsValid() const
{
  return ( GetNbEntities(SMDSAbs_0DElement) > 0 ||
           GetNbEntities(SMDSAbs_Ball     ) > 0 ||
           GetNbEntities(SMDSAbs_Edge     ) > 0 ||
           GetNbEntities(SMDSAbs_Face     ) > 0 ||
           GetNbEntities(SMDSAbs_Volume   ) > 0 ||
           GetNbEntities(SMDSAbs_Node     ) > 0 );
}

//=================================================================================
// function : updateEntitiesFlags
// purpose  : Update entities flags
//=================================================================================
void SMESH_VisualObjDef::updateEntitiesFlags()
{
  unsigned int tmp = myEntitiesState;
  ClearEntitiesFlags();

  map<SMDSAbs_ElementType,int> entities = SMESH::GetEntitiesFromObject(this);


  if( myEntitiesCache[SMDSAbs_0DElement] != 0 ||
      myEntitiesCache[SMDSAbs_0DElement] >= entities[SMDSAbs_0DElement] )
    myEntitiesState &= ~SMESH_Actor::e0DElements;

  if( myEntitiesCache[SMDSAbs_Ball] != 0 ||
      myEntitiesCache[SMDSAbs_Ball] >= entities[SMDSAbs_Ball] )
    myEntitiesState &= ~SMESH_Actor::eBallElem;

  if( myEntitiesCache[SMDSAbs_Edge] != 0 ||
      myEntitiesCache[SMDSAbs_Edge] >= entities[SMDSAbs_Edge] )
    myEntitiesState &= ~SMESH_Actor::eEdges;

  if( myEntitiesCache[SMDSAbs_Face] != 0 ||
      myEntitiesCache[SMDSAbs_Face] >= entities[SMDSAbs_Face] )
    myEntitiesState &= ~SMESH_Actor::eFaces;

  if( myEntitiesCache[SMDSAbs_Volume] != 0 ||
      myEntitiesCache[SMDSAbs_Volume] >= entities[SMDSAbs_Volume] )
    myEntitiesState &= ~SMESH_Actor::eVolumes;

  if( tmp != myEntitiesState ) {
    myEntitiesFlag = true;
  }

  myEntitiesCache = entities;
}

//=================================================================================
// function : ClearEntitiesFlags
// purpose  : Clear the entities flags
//=================================================================================
void SMESH_VisualObjDef::ClearEntitiesFlags()
{
  myEntitiesState = SMESH_Actor::eAllEntity;
  myEntitiesFlag = false;
}

//=================================================================================
// function : GetEntitiesFlag
// purpose  : Return the entities flag
//=================================================================================
bool SMESH_VisualObjDef::GetEntitiesFlag()
{
  return myEntitiesFlag;
}

//=================================================================================
// function : GetEntitiesState
// purpose  : Return the entities state
//=================================================================================
unsigned int SMESH_VisualObjDef::GetEntitiesState()
{
  return myEntitiesState;
}

/*
  Class       : SMESH_MeshObj
  Description : Class for visualisation of mesh
*/

//=================================================================================
// function : SMESH_MeshObj
// purpose  : Constructor
//=================================================================================
SMESH_MeshObj::SMESH_MeshObj(SMESH::SMESH_Mesh_ptr theMesh):
  myClient(SalomeApp_Application::orb(),theMesh)
{
        myEmptyGrid = 0;
  if ( MYDEBUG ) 
    MESSAGE("SMESH_MeshObj - this = "<<this<<"; theMesh->_is_nil() = "<<theMesh->_is_nil());
}

//=================================================================================
// function : ~SMESH_MeshObj
// purpose  : Destructor
//=================================================================================
SMESH_MeshObj::~SMESH_MeshObj()
{
  if ( MYDEBUG ) 
    MESSAGE("SMESH_MeshObj - this = "<<this<<"\n");
  if ( myEmptyGrid )
    myEmptyGrid->Delete();
}

//=================================================================================
// function : Update
// purpose  : Update mesh and fill grid with new values if necessary 
//=================================================================================
bool SMESH_MeshObj::Update( int theIsClear )
{
  // Update SMDS_Mesh on client part
  if ( MYDEBUG ) MESSAGE("SMESH_MeshObj::Update " << this);
  if ( myClient.Update(theIsClear) || GetUnstructuredGrid()->GetNumberOfPoints()==0) {
    if ( MYDEBUG ) MESSAGE("buildPrs");
    buildPrs();  // Fill unstructured grid
    return true;
  }
  return false;
}

bool SMESH_MeshObj::NulData()
{
  if ( MYDEBUG ) MESSAGE ("SMESH_MeshObj::NulData() =============================================");
  if (!myEmptyGrid)
  {
    myEmptyGrid = SMDS_UnstructuredGrid::New();
    myEmptyGrid->Initialize();
    myEmptyGrid->Allocate();
    vtkPoints* points = vtkPoints::New();
    points->SetNumberOfPoints(0);
    myEmptyGrid->SetPoints( points );
    points->Delete();
    //myEmptyGrid->BuildLinks();
  }
  myGrid->ShallowCopy(myEmptyGrid);
  return true;
}
//=================================================================================
// function : GetElemDimension
// purpose  : Get dimension of element
//=================================================================================
int SMESH_MeshObj::GetElemDimension( const int theObjId )
{
  const SMDS_MeshElement* anElem = myClient->FindElement( theObjId );
  if ( anElem == 0 )
    return 0;

  int aType = anElem->GetType();
  switch ( aType )
  {
    case SMDSAbs_0DElement : return 0;
    case SMDSAbs_Ball : return 0;
    case SMDSAbs_Edge  : return 1;
    case SMDSAbs_Face  : return 2;
    case SMDSAbs_Volume: return 3;
    default            : return 0;
  }
}

//=================================================================================
// function : GetEntities
// purpose  : Get entities of specified type. Return number of entities
//=================================================================================
int SMESH_MeshObj::GetNbEntities( const SMDSAbs_ElementType theType) const
{
  switch ( theType )
  {
    case SMDSAbs_Node:
    {
      return myClient->NbNodes();
    }
    break;
    case SMDSAbs_0DElement:
    {
      return myClient->Nb0DElements();
    }
    case SMDSAbs_Ball:
    {
      return myClient->NbBalls();
    }
    break;
    case SMDSAbs_Edge:
    {
      return myClient->NbEdges();
    }
    break;
    case SMDSAbs_Face:
    {
      return myClient->NbFaces();
    }
    break;
    case SMDSAbs_Volume:
    {
      return myClient->NbVolumes();
    }
    break;
    default:
      return 0;
    break;
  }
}

int SMESH_MeshObj::GetEntities( const SMDSAbs_ElementType theType, TEntityList& theObjs ) const
{
  theObjs.clear();

  switch ( theType )
  {
    case SMDSAbs_Node:
    {
      SMDS_NodeIteratorPtr anIter = myClient->nodesIterator();
      while ( anIter->more() ) theObjs.push_back( anIter->next() );
    }
    break;
    case SMDSAbs_0DElement:
    {
      SMDS_ElemIteratorPtr anIter = myClient->elementsIterator(SMDSAbs_0DElement);
      while ( anIter->more() ) theObjs.push_back( anIter->next() );
    }
    break;
    case SMDSAbs_Ball:
    {
      SMDS_ElemIteratorPtr anIter = myClient->elementGeomIterator(SMDSGeom_BALL);
      while ( anIter->more() ) theObjs.push_back( anIter->next() );
    }
    break;
    case SMDSAbs_Edge:
    {
      SMDS_EdgeIteratorPtr anIter = myClient->edgesIterator();
      while ( anIter->more() ) theObjs.push_back( anIter->next() );
    }
    break;
    case SMDSAbs_Face:
    {
      SMDS_FaceIteratorPtr anIter = myClient->facesIterator();
      while ( anIter->more() ) theObjs.push_back( anIter->next() );
    }
    break;
    case SMDSAbs_Volume:
    {
      SMDS_VolumeIteratorPtr anIter = myClient->volumesIterator();
      while ( anIter->more() ) theObjs.push_back( anIter->next() );
    }
    break;
    default:
    break;
  }

  return theObjs.size();
}

//=================================================================================
// function : UpdateFunctor
// purpose  : Update functor in accordance with current mesh
//=================================================================================
void SMESH_MeshObj::UpdateFunctor( const SMESH::Controls::FunctorPtr& theFunctor )
{
  theFunctor->SetMesh( GetMesh() );
}

//=================================================================================
// function : IsNodePrs
// purpose  : Return true if node presentation is used
//=================================================================================
bool SMESH_MeshObj::IsNodePrs() const
{
  return myClient->Nb0DElements() + myClient->NbEdges() + myClient->NbFaces() + myClient->NbVolumes() + myClient->NbBalls() == 0 ;
}


/*
  Class       : SMESH_SubMeshObj
  Description : Base class for visualisation of submeshes and groups
*/

//=================================================================================
// function : SMESH_SubMeshObj
// purpose  : Constructor
//=================================================================================
SMESH_SubMeshObj::SMESH_SubMeshObj( SMESH_MeshObj* theMeshObj )
{
  if ( MYDEBUG ) MESSAGE( "SMESH_SubMeshObj - theMeshObj = " << theMeshObj );
  
  myMeshObj = theMeshObj;
}

SMESH_SubMeshObj::~SMESH_SubMeshObj()
{
}

//=================================================================================
// function : GetElemDimension
// purpose  : Get dimension of element
//=================================================================================
int SMESH_SubMeshObj::GetElemDimension( const int theObjId )
{
  return myMeshObj == 0 ? 0 : myMeshObj->GetElemDimension( theObjId );
}

//=================================================================================
// function : UpdateFunctor
// purpose  : Update functor in accordance with current mesh
//=================================================================================

void SMESH_SubMeshObj::UpdateFunctor( const SMESH::Controls::FunctorPtr& theFunctor )
{
  theFunctor->SetMesh( myMeshObj->GetMesh() );
}

//=================================================================================
// function : Update
// purpose  : Update mesh object and fill grid with new values 
//=================================================================================
bool SMESH_SubMeshObj::Update( int theIsClear )
{
  if ( MYDEBUG ) MESSAGE("SMESH_SubMeshObj::Update " << this)
  bool changed = myMeshObj->Update( theIsClear );
  buildPrs(true);
  return changed;
}


/*
  Class       : SMESH_GroupObj
  Description : Class for visualisation of groups
*/

//=================================================================================
// function : SMESH_GroupObj
// purpose  : Constructor
//=================================================================================
SMESH_GroupObj::SMESH_GroupObj( SMESH::SMESH_GroupBase_ptr theGroup, 
                                SMESH_MeshObj*             theMeshObj )
: SMESH_SubMeshObj( theMeshObj ),
  myGroupServer( SMESH::SMESH_GroupBase::_duplicate(theGroup) )
{
  if ( MYDEBUG ) MESSAGE("SMESH_GroupObj - theGroup->_is_nil() = "<<theGroup->_is_nil());
  myGroupServer->Register();
}

SMESH_GroupObj::~SMESH_GroupObj()
{
  if ( MYDEBUG ) MESSAGE("~SMESH_GroupObj");
  myGroupServer->UnRegister();
}

//=================================================================================
// function : IsNodePrs
// purpose  : Return true if node presentation is used
//=================================================================================
bool SMESH_GroupObj::IsNodePrs() const
{
  return myGroupServer->GetType() == SMESH::NODE;
}

//=================================================================================
// function : GetElementType
// purpose  : Return type of elements of group
//=================================================================================
SMDSAbs_ElementType SMESH_GroupObj::GetElementType() const
{
  return SMDSAbs_ElementType(myGroupServer->GetType());
}

//=================================================================================
// function : getNodesFromElems
// purpose  : Retrieve nodes from elements
//=================================================================================
static int getNodesFromElems( SMESH::long_array_var&              theElemIds,
                              const SMDS_Mesh*                    theMesh,
                              std::list<const SMDS_MeshElement*>& theResList )
{
  set<const SMDS_MeshElement*> aNodeSet;

  for ( CORBA::Long i = 0, n = theElemIds->length(); i < n; i++ )
  {
    const SMDS_MeshElement* anElem = theMesh->FindElement( theElemIds[ i ] );
    if ( anElem != 0 )
    {
      SMDS_ElemIteratorPtr anIter = anElem->nodesIterator();
      while ( anIter->more() )
      {
        const SMDS_MeshElement* aNode = anIter->next();
        if ( aNode != 0 )
          aNodeSet.insert( aNode );
      }
    }
  }

  set<const SMDS_MeshElement*>::const_iterator anIter;
  for ( anIter = aNodeSet.begin(); anIter != aNodeSet.end(); ++anIter )
    theResList.push_back( *anIter );

  return theResList.size();    
}

//=================================================================================
// function : getPointers
// purpose  : Get std::list<const SMDS_MeshElement*> from list of IDs
//=================================================================================
static int getPointers( const SMDSAbs_ElementType           theRequestType,
                        SMESH::long_array_var&              theElemIds,
                        const SMDS_Mesh*                    theMesh,
                        std::list<const SMDS_MeshElement*>& theResList )
{
  for ( CORBA::Long i = 0, n = theElemIds->length(); i < n; i++ )
  {
    const SMDS_MeshElement* anElem = theRequestType == SMDSAbs_Node
      ? theMesh->FindNode( theElemIds[ i ] ) : theMesh->FindElement( theElemIds[ i ] );

    if ( anElem != 0 )
      theResList.push_back( anElem );
  }

  return theResList.size();
}


//=================================================================================
// function : GetEntities
// purpose  : Get entities of specified type. Return number of entities
//=================================================================================
int SMESH_GroupObj::GetNbEntities( const SMDSAbs_ElementType theType) const
{
  if(SMDSAbs_ElementType(myGroupServer->GetType()) == theType) {
    return myGroupServer->Size();
  }
  if ( theType == SMDSAbs_Node ) {
    return myGroupServer->GetNumberOfNodes();
  }
  return 0;
}

int SMESH_GroupObj::GetEntities( const SMDSAbs_ElementType theType, TEntityList& theResList ) const
{
  theResList.clear();
  SMDS_Mesh* aMesh = myMeshObj->GetMesh();
  
  if ( aMesh == 0 )
    return 0;

  SMDSAbs_ElementType aGrpType = SMDSAbs_ElementType(myGroupServer->GetType());
  if ( aGrpType != theType && theType != SMDSAbs_Node )
    return 0;

  SMESH::long_array_var anIds = myGroupServer->GetListOfID();
  if ( anIds->length() == 0 )
    return 0;

  if ( aGrpType == theType )
    return getPointers( theType, anIds, aMesh, theResList );
  else if ( theType == SMDSAbs_Node )
    return getNodesFromElems( anIds, aMesh, theResList );
  else
    return 0;
}



/*
  Class       : SMESH_subMeshObj
  Description : Class for visualisation of submeshes
*/

//=================================================================================
// function : SMESH_subMeshObj
// purpose  : Constructor
//=================================================================================
SMESH_subMeshObj::SMESH_subMeshObj( SMESH::SMESH_subMesh_ptr theSubMesh,
                                    SMESH_MeshObj*           theMeshObj )
: SMESH_SubMeshObj( theMeshObj ),
  mySubMeshServer( SMESH::SMESH_subMesh::_duplicate( theSubMesh ) )
{
  if ( MYDEBUG ) MESSAGE( "SMESH_subMeshObj - theSubMesh->_is_nil() = " << theSubMesh->_is_nil() );
  
  mySubMeshServer->Register();
}

SMESH_subMeshObj::~SMESH_subMeshObj()
{
  if ( MYDEBUG ) MESSAGE( "~SMESH_subMeshObj" );
  mySubMeshServer->UnRegister();
}

//=================================================================================
// function : GetEntities
// purpose  : Get entities of specified type. Return number of entities
//=================================================================================
int SMESH_subMeshObj::GetNbEntities( const SMDSAbs_ElementType theType) const
{
  switch ( theType )
  {
    case SMDSAbs_Node:
    {
      return mySubMeshServer->GetNumberOfNodes( /*all=*/true );
    }
    break;
    case SMDSAbs_Ball:
    case SMDSAbs_0DElement:
    case SMDSAbs_Edge:
    case SMDSAbs_Face:
    case SMDSAbs_Volume:
    {
      SMESH::long_array_var anIds = 
        mySubMeshServer->GetElementsByType( SMESH::ElementType(theType) );
      return anIds->length();
    }
    default:
      return 0;
    break;
  }
}

int SMESH_subMeshObj::GetEntities( const SMDSAbs_ElementType theType, TEntityList& theResList ) const
{
  theResList.clear();

  SMDS_Mesh* aMesh = myMeshObj->GetMesh();
  if ( aMesh == 0 )
    return 0;

  bool isNodal = IsNodePrs();

  if ( isNodal )
  {
    if ( theType == SMDSAbs_Node )
    {
      SMESH::long_array_var anIds = mySubMeshServer->GetNodesId();
      return getPointers( SMDSAbs_Node, anIds, aMesh, theResList );
    }
  }
  else
  {
    if ( theType == SMDSAbs_Node )
    {
      SMESH::long_array_var anIds = mySubMeshServer->GetElementsId();
      return getNodesFromElems( anIds, aMesh, theResList );
    }
    else
    {
      SMESH::long_array_var anIds = 
        mySubMeshServer->GetElementsByType( SMESH::ElementType(theType) );
      return getPointers( theType, anIds, aMesh, theResList );
    }
  }

  return 0;
}

//=================================================================================
// function : IsNodePrs
// purpose  : Return true if node presentation is used
//=================================================================================
bool SMESH_subMeshObj::IsNodePrs() const
{
  return mySubMeshServer->GetNumberOfElements() == 0;
}
