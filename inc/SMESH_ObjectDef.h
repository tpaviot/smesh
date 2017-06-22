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
//  File   : SMESH_Object.h
//  Author : Nicolas REJNERI
//  Module : SMESH
//
#ifndef SMESH_OBJECTDEF_H
#define SMESH_OBJECTDEF_H

#include "SMESH_Controls.hxx"
#include "SMESH_Object.h"
#include "SMESH_Client.hxx"
#include "SMESH_Actor.h"

// IDL Headers
#include <SALOMEconfig.h>
#include CORBA_SERVER_HEADER(SMESH_Mesh)
#include CORBA_SERVER_HEADER(SMESH_Group)

#include <map>
#include <list>

class vtkPoints;
class SALOME_ExtractUnstructuredGrid;

class SMDS_MeshNode;
class SMDS_MeshElement;

/*
  Class       : SMESH_VisualObj
  Description : Base class for all mesh objects to be visuilised
*/
class SMESHOBJECT_EXPORT SMESH_VisualObjDef: public SMESH_VisualObj
{
public:
  typedef std::list<const SMDS_MeshElement*>   TEntityList;
  typedef std::map<vtkIdType,vtkIdType>  TMapOfIds;
  
                            SMESH_VisualObjDef();
  virtual                   ~SMESH_VisualObjDef();
  
  virtual bool              Update( int theIsClear = true ) = 0;
  virtual bool              NulData() {return 0; };
  virtual void              UpdateFunctor( const SMESH::Controls::FunctorPtr& theFunctor ) = 0;
  virtual int               GetElemDimension( const int theObjId ) = 0;

  virtual int               GetNbEntities( const SMDSAbs_ElementType theType) const = 0;
  virtual int               GetEntities( const SMDSAbs_ElementType, TEntityList& ) const = 0;
  virtual bool              IsNodePrs() const = 0;
  virtual SMDS_Mesh*        GetMesh() const = 0;
  virtual SMESH::SMESH_Mesh_ptr GetMeshServer() = 0;

  virtual bool              IsValid() const;

  virtual bool              GetEdgeNodes( const int theElemId,
                                          const int theEdgeNum,
                                          int&      theNodeId1,
                                          int&      theNodeId2 ) const;

  virtual vtkUnstructuredGrid* GetUnstructuredGrid();
  
  virtual vtkIdType         GetNodeObjId( int theVTKID );
  virtual vtkIdType         GetNodeVTKId( int theObjID );
  virtual vtkIdType         GetElemObjId( int theVTKID );
  virtual vtkIdType         GetElemVTKId( int theObjID );
  
  virtual void              ClearEntitiesFlags();
  virtual bool              GetEntitiesFlag();
  virtual unsigned int      GetEntitiesState();
  
protected:

  void                      createPoints( vtkPoints* );
  void                      buildPrs(bool buildGrid = false);
  void                      buildNodePrs();
  void                      buildElemPrs();
  void                      updateEntitiesFlags();
//private:

  TMapOfIds                 mySMDS2VTKNodes;
  TMapOfIds                 myVTK2SMDSNodes;
  TMapOfIds                 mySMDS2VTKElems;
  TMapOfIds                 myVTK2SMDSElems;
  bool                      myLocalGrid;

  bool                      myEntitiesFlag;
  unsigned int              myEntitiesState;

  vtkUnstructuredGrid*      myGrid;
  std::map<SMDSAbs_ElementType,int> myEntitiesCache;
};


/*
  Class       : SMESH_MeshObj
  Description : Class for visualisation of mesh
*/

class SMESHOBJECT_EXPORT SMESH_MeshObj: public SMESH_VisualObjDef
{
public:

                            SMESH_MeshObj( SMESH::SMESH_Mesh_ptr );
  virtual                   ~SMESH_MeshObj();
  
  virtual bool              Update( int theIsClear = true );
  virtual bool              NulData();
  
  virtual int               GetNbEntities( const SMDSAbs_ElementType) const;
  virtual int               GetEntities( const SMDSAbs_ElementType, TEntityList& ) const;
  virtual bool              IsNodePrs() const;

  virtual int               GetElemDimension( const int theObjId );

  virtual void              UpdateFunctor( const SMESH::Controls::FunctorPtr& theFunctor );
  
  virtual SMESH::SMESH_Mesh_ptr GetMeshServer() { return myClient.GetMeshServer(); }
  virtual SMDS_Mesh*        GetMesh() const { return myClient.GetMesh(); }

protected:

  SMESH_Client              myClient;
  vtkUnstructuredGrid*      myEmptyGrid;
};


/*
  Class       : SMESH_SubMeshObj
  Description : Base class for visualisation of submeshes and groups
*/

class SMESHOBJECT_EXPORT SMESH_SubMeshObj: public SMESH_VisualObjDef
{
public:

                            SMESH_SubMeshObj(SMESH_MeshObj* theMeshObj);
  virtual                   ~SMESH_SubMeshObj();

  virtual bool              Update( int theIsClear = true );
  
  virtual void              UpdateFunctor( const SMESH::Controls::FunctorPtr& theFunctor );
  virtual int               GetElemDimension( const int theObjId );
  virtual SMDS_Mesh*        GetMesh() const { return myMeshObj->GetMesh(); }
  virtual SMESH::SMESH_Mesh_ptr GetMeshServer() { return myMeshObj->GetMeshServer(); }
  
protected:

  SMESH_MeshObj*            myMeshObj;
};


/*
  Class       : SMESH_GroupObj
  Description : Class for visualisation of groups
*/

class SMESHOBJECT_EXPORT SMESH_GroupObj: public SMESH_SubMeshObj
{
public:
                            SMESH_GroupObj( SMESH::SMESH_GroupBase_ptr, SMESH_MeshObj* );
  virtual                   ~SMESH_GroupObj();

  virtual int               GetNbEntities( const SMDSAbs_ElementType) const;
  virtual int               GetEntities( const SMDSAbs_ElementType, TEntityList& ) const;
  virtual bool              IsNodePrs() const;

  virtual SMDSAbs_ElementType GetElementType() const;

private:

  SMESH::SMESH_GroupBase_var    myGroupServer;
};


/*
  Class       : SMESH_subMeshObj
  Description : Class for visualisation of submeshes
*/

class SMESHOBJECT_EXPORT SMESH_subMeshObj : public SMESH_SubMeshObj
{
public:

                            SMESH_subMeshObj( SMESH::SMESH_subMesh_ptr, 
                                              SMESH_MeshObj* );
  virtual                   ~SMESH_subMeshObj();

  virtual int               GetNbEntities( const SMDSAbs_ElementType) const;
  virtual int               GetEntities( const SMDSAbs_ElementType, TEntityList& ) const;
  virtual bool              IsNodePrs() const;    
  
protected:

  SMESH::SMESH_subMesh_var  mySubMeshServer;
};


#endif
