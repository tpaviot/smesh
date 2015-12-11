//  SMESH  SMESH_MeshVSLink2 : Connection of SMESH with MeshVS from OCC 
//
//  Copyright (C) 2003  OPEN CASCADE, EADS/CCR, LIP6, CEA/DEN,
//  CEDRAT, EDF R&D, LEG, PRINCIPIA R&D, BUREAU VERITAS
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
// See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//
// File      : SMESH_MeshVSLink2.cxx
// Created   : 7/1/2012
// Author    : Mark Blome
// Module    : SMESH
// Based on the implementation in SMESH_MeshVSLink.cxx/hxx
// Extended functionality: 
//    - allows individual submeshes to be displayed independantly
//    - supports slicing volume/surface meshes by a cut-plane
//    - considerable efficiency improvement for volume meshes in shaded display mode:
//      only boundary faces of volume grids are displayed taking a defined cut plane into 
//      account. In Wireframe/shrink display modes all elements are shown as it is done in 
//      SMESH_MeshVSLink class

#ifndef _SMESH_MeshVSLink2_HeaderFile
#define _SMESH_MeshVSLink2_HeaderFile

#include <Standard.hxx>
#include <Handle_SMESH_MeshVSLink2.hxx>
#include <TColStd_PackedMapOfInteger.hxx>
#include <Handle_TColStd_HArray2OfInteger.hxx>
#include <Handle_TColStd_HArray2OfReal.hxx>
#include <MeshVS_DataSource3D.hxx>
#include <Standard_Boolean.hxx>
#include <Standard_Integer.hxx>
#include <MeshVS_EntityType.hxx>
#include <Standard_Address.hxx>
#include <Handle_TColStd_HArray1OfInteger.hxx>
#include <TColStd_HPackedMapOfInteger.hxx>
#include <Standard_Real.hxx>
#include <SMESH_Mesh.hxx>
#include <SMESH_subMesh.hxx>

#include <boost/functional/hash.hpp>

#include <map>
#include <vector>
#include <set>
#include <algorithm>

#include "SMDSAbs_ElementType.hxx"
#include <TColgp_Array1OfXYZ.hxx>
#include <TColStd_HPackedMapOfInteger.hxx>
#include <MeshVS_HArray1OfSequenceOfInteger.hxx>
#include <boost/functional/hash.hpp>

enum VisTypeofElem {SMESHElem, FaceofSMESHElem};

struct XYZ_ {
  double x;
  double y;
  double z;
  XYZ_()                               { x = 0; y = 0; z = 0; }
  XYZ_( double X, double Y, double Z ) { x = X; y = Y; z = Z; }
  XYZ_( const SMDS_MeshNode* n )       { x = n->X(); y = n->Y(); z = n->Z(); }
};

class VisMeshElementList
{
public:
  VisMeshElementList(const SMESH_Mesh *aMesh, int displaymode, SMESH_subMesh *submesh, bool supportsmoothshading=true);
  void SetDisplayMode(int displaymode); // 0:MeshVS_DMF_Shading, 1:MeshVS_DMF_WireFrame, 2:MeshVS_DMF_Shrink
  
  const TColStd_PackedMapOfInteger& GetAllElements();
  bool GetElement(int elemnumber, Standard_Integer& NbNodes, Handle(MeshVS_HArray1OfSequenceOfInteger)& Data); 			    
  bool GetElementNodes(int elemnumber, TColStd_Array1OfReal& Coords, Standard_Integer& NbNodes); 
  bool GetElementType(int elemnumber, MeshVS_EntityType& elementtype);
  void GetHiddenElements(int planedir, double planevalue, TColStd_PackedMapOfInteger& hiddenelements);  
  bool GetNodesByElement(int elemnumber, TColStd_Array1OfInteger& NodeIDs, int& NbNodes );
  bool GetNodeNormal(int nodenumber, int elemnumber, double &nx, double &ny, double &nz);
  bool GetNormal(int elemnumber, double& nx, double& ny, double& nz);
  
private:
  bool GetFaceofSMESHElement(int elemnumber, Standard_Integer& NbNodes, Handle(MeshVS_HArray1OfSequenceOfInteger)& Data);
  bool GetSMESHElement(int smeshid, Standard_Integer& NbNodes, Handle(MeshVS_HArray1OfSequenceOfInteger)& Data);
  void ComputeMeshElements();
  bool ComputeNormal(int smeshelem, double& nx, double& ny,double& nz);
  void addedges();
  void addfaces();
  void addvolumes(); 
  void addreversemapping(int faceid, int volsmeshid);
  int addSMESHElement(const SMDS_MeshElement * elem);                        // returns element number
  int addFaceofSMESHElement(std::vector<int>& smeshnodeids, int volsmeshid); // returns element number 
  void addFacesofSMESHElement(const SMDS_MeshElement * elem);
  XYZ_ ComputeCentroid(const SMDS_MeshElement* elem);
  void ComputeNodalFaceNormals();
  double ElemArea(const SMDS_MeshElement * elem);
  std::size_t TriHashKey(int nodeid1, int nodeid2, int nodeid3);
  std::size_t QuadHashKey(int nodeid1, int nodeid2, int nodeid3, int nodeid4);  
  bool GetCentroid(int elemnumber, XYZ_& centroid);
  bool IsValidElementNumber(int elemnumber) const;  
  bool IsHidden(int elemnumber, int planedir, double planevalue);
  bool IsElemHidden(int smeshelemid, int planedir, double planevalue);
  void clear();
  std::map<int, std::set<int> > reversemap; // maps volume face id to a set of smesh element ids (the volume elements 
                                            // the volume face is a neighbor to in the same submesh (if applicable))
  
  std::map<int, XYZ_> ElemCentroids; // maps smesh element ids to element centroids  
                                     // no centroids are computed for volume faces because their hidden state is being decided
                                     // based on the hidden state of the volume elements they are connected to

  std::map<int, XYZ_> ElemNormals;   // maps smesh element ids to element normal vectors
  std::map<int, XYZ_> NodalNormals;  // maps smesh node ids to (area weighted mean) nodal normal vectors, only comp. if _supportsmoothshading 

  std::map<int, std::pair<int, VisTypeofElem> > elemnum2idmap;  // maps element numbers (as returned by GetAllElements() and to be used 
                                                                // with the public member functions above) to internally used ids and element types
                                                                // for element types == SMESHElem ids are smesh element ids (as used in mesh/submeshds)
                                                                // for element types == FaceofSMESHElem ids are volume face ids (index to volumefaces)

  std::map<std::size_t, int> facehashkeys;
  TColStd_PackedMapOfInteger elements; 
  std::vector<std::vector<int> > volumefaces; // smesh node ids for each face element, idx is volume face id

  const SMESH_Mesh *mesh;
  SMESH_subMesh *submesh;
  SMESHDS_SubMesh* submeshds;   
  int _displaymode;
  bool _supportsmoothshading;

  std::vector<int> int3, int4;
  boost::hash<std::string> string_hash;
};

class SMESH_MeshVSLink2 : public MeshVS_DataSource3D 
{
  
public:
  
  Standard_EXPORT SMESH_MeshVSLink2(const SMESH_Mesh *aMesh, SMESH_subMesh *aSubMesh, int displaymode);

  Standard_EXPORT void SetDisplayMode(int displaymode);  // 0:MeshVS_DMF_Shading, 1:MeshVS_DMF_WireFrame, 2:MeshVS_DMF_Shrink  
  
  Standard_EXPORT void GetHiddenNodes(int planedir, double value, TColStd_PackedMapOfInteger& hiddennodes);
  Standard_EXPORT void GetHiddenElements(int planedir, double value, TColStd_PackedMapOfInteger& hiddenelems); 

  Standard_EXPORT   Standard_Boolean GetGeom(const Standard_Integer ID,const Standard_Boolean IsElement,TColStd_Array1OfReal& Coords,Standard_Integer& NbNodes,MeshVS_EntityType& Type) const;
  Standard_EXPORT   Standard_Boolean Get3DGeom(const Standard_Integer ID,Standard_Integer& NbNodes,Handle(MeshVS_HArray1OfSequenceOfInteger)& Data) const;
  Standard_EXPORT   Standard_Boolean GetGeomType(const Standard_Integer ID,const Standard_Boolean IsElement,MeshVS_EntityType& Type) const; 
  Standard_EXPORT   Standard_Address GetAddr(const Standard_Integer ID,const Standard_Boolean IsElement) const;
  Standard_EXPORT   Standard_Boolean GetNodesByElement(const Standard_Integer ID,TColStd_Array1OfInteger& NodeIDs,Standard_Integer& NbNodes) const; 
  Standard_EXPORT Standard_Boolean GetNormal(const Standard_Integer Id,const Standard_Integer Max,Standard_Real& nx,Standard_Real& ny,Standard_Real& nz) const;
  Standard_EXPORT Standard_Boolean GetNodeNormal(const Standard_Integer ranknode, const Standard_Integer Id, Standard_Real &nx, Standard_Real &ny, Standard_Real &nz) const;

  Standard_EXPORT void GetAllGroups(TColStd_PackedMapOfInteger& Ids) const;
  
  Standard_EXPORT const Handle(Standard_Type)& DynamicType() const;

  Standard_EXPORT  const TColStd_PackedMapOfInteger& GetAllNodes() const;
  Standard_EXPORT  const TColStd_PackedMapOfInteger& GetAllElements() const;
  
private:
 
  const SMESH_Mesh *mesh;
  SMESH_subMesh *submesh;
  SMESHDS_SubMesh* submeshds; 
  int _displaymode;
  mutable VisMeshElementList myElements; 
  TColStd_PackedMapOfInteger myNodes;
  TColStd_PackedMapOfInteger myGroups;
};


#endif
