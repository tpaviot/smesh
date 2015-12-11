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

//local headers
#include <SMESH_MeshVSLink2.ixx>
#include <SMESHDS_Group.hxx>
#include <SMESHDS_SubMesh.hxx>
#include <SMESH_subMesh.hxx>


//occ headers
#include <TColgp_Array1OfXYZ.hxx>
#include <MeshVS_HArray1OfSequenceOfInteger.hxx>
#include <SMDS_VolumeTool.hxx>
#include <SMDS_VolumeOfFaces.hxx>
#include <MeshVS_DisplayModeFlags.hxx>

#include <gp_XYZ.hxx>

#include <TColgp_Array1OfXYZ.hxx>
#include <TColgp_Array1OfVec.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TColStd_HPackedMapOfInteger.hxx>
#include <TColStd_MapIteratorOfPackedMapOfInteger.hxx>
#include <utilities.h>

#define MAX_SORT_NODE_COUNT 12
#define RETURN_BAD_RESULT(msg) { MESSAGE(")-: Error: " << msg); return false; }

typedef std::map<double, int> T_Double_NodeID_Map;


bool sortNodes2(const SMDS_MeshElement* theTool, const int* idNodes, int theNodeCount, int *myResult)
{
  if (theNodeCount < 3) RETURN_BAD_RESULT("Number of nodes < 3");
  //INITIAL VARS
  TColgp_Array1OfXYZ myNodeList(1, theNodeCount);
  TColgp_Array1OfVec myVectList(1, theNodeCount);
  TColStd_Array1OfReal myAngleList(1, theNodeCount);
  gp_XYZ BaryCenter(0.,0.,0.);
  //int myResult[MAX_SORT_NODE_COUNT];
  //INITIALIZE THE POINTS
  for (int i = 1; i <= theNodeCount; i++ ) {
	const SMDS_MeshNode *n = theTool->GetNode( idNodes[i-1] );
	gp_XYZ aPoint(n->X(), n->Y(), n->Z());
	myNodeList.SetValue(i, aPoint);
  }
  //CALCULATE THE BARYCENTER
  for (int i = 1; i <= theNodeCount; i++ )
	BaryCenter += myNodeList.Value(i);
  BaryCenter /= theNodeCount;
  //CREATE THE VECTORS
  for (int i = 1; i <= theNodeCount; i++ ) {
	gp_Vec aVector(BaryCenter, myNodeList.Value(i));
	myVectList.SetValue(i, aVector);
  }
  //CALCULATE THE NORMAL USING FIRST THREE POINTS
  gp_XYZ q1 = myNodeList.Value(2)-myNodeList.Value(1);
  gp_XYZ q2 = myNodeList.Value(3)-myNodeList.Value(1);
  gp_XYZ normal  = q1 ^ q2;
  double modul = normal.Modulus();
  if ( modul > 0 )
	normal /= modul;
  //COUNT THE ANGLE OF THE FIRST WITH EACH
  for (int i = 1; i <= theNodeCount; i++ )
	myAngleList.SetValue(i, myVectList.Value(1).AngleWithRef(myVectList.Value(i), normal));
  //CREATE THE RESULT MAP (WILL SORT THE VERTICES)
  T_Double_NodeID_Map myMap;
  for (int i = 1; i <= theNodeCount; i++ )
	myMap.insert( make_pair(myAngleList.Value(i), idNodes[i-1]));
  int resID = 0;
  T_Double_NodeID_Map::iterator it;
  for(it = myMap.begin(); it!= myMap.end(); ++it)
	myResult[resID++] = it->second;

  return true;
}

VisMeshElementList::VisMeshElementList(const SMESH_Mesh *aMesh, int displaymode, SMESH_subMesh *aSubMesh, bool supportsmoothshading)
  : mesh(aMesh), submesh(aSubMesh), submeshds(NULL), _displaymode(displaymode), int3(3,0), int4(4,0), _supportsmoothshading(supportsmoothshading)
{
  if (submesh) submeshds = submesh->CreateSubMeshDS();
  ComputeMeshElements();    
}

int VisMeshElementList::addSMESHElement(const SMDS_MeshElement * elem)      
{
  // returns element number
  if (!elem) return -1;
  int newid = elemnum2idmap.size(); elements.Add(newid);  
  elemnum2idmap[newid]=std::pair<int, VisTypeofElem>(elem->GetID(), SMESHElem);
  return newid;
} 

void VisMeshElementList::addreversemapping(int faceid, int volsmeshid)
{
  if (reversemap.find(faceid)==reversemap.end())
    {
      std::set<int> smeshids;  smeshids.insert(volsmeshid);
      reversemap[faceid]=smeshids; 
    }
  else reversemap[faceid].insert(volsmeshid); 
}

int VisMeshElementList::addFaceofSMESHElement(std::vector<int>& smeshnodeids, int volsmeshid)
{
  // returns element number if face was added, otherwise (if it already exists) -1
  int nbnodes=smeshnodeids.size(); std::size_t hashkey=0;
  if (nbnodes==3) hashkey = TriHashKey(smeshnodeids[0], smeshnodeids[1], smeshnodeids[2]);
  else if (nbnodes==4) hashkey = QuadHashKey(smeshnodeids[0], smeshnodeids[1], smeshnodeids[2], smeshnodeids[3]);
  else
    { 
      MESSAGE("addFaceofSMESHElement: Unknown volume mesh element face type (#nodes: " << nbnodes << ")"); 
      return -1;
    }
  if (facehashkeys.find(hashkey)==facehashkeys.end()) // face was not added yet
    {
      int newid = elemnum2idmap.size(); volumefaces.push_back(smeshnodeids); elements.Add(newid); 
      int volumefaceid = volumefaces.size()-1;
      elemnum2idmap[newid]=std::pair<int, VisTypeofElem>(volumefaceid, FaceofSMESHElem);
      addreversemapping(volumefaceid, volsmeshid); facehashkeys[hashkey]=volumefaceid;      
      //MESSAGE("New face: ID: " << newid << "; #nodes: "<< nbnodes << "; reverse mapping: " << volumefaceid << " -> " << volsmeshid << "; hashkey: " << hashkey);
      return newid; // newid = element number 
    }
  else
    {
      // face element already added, need to add entry to reverse map nevertheless
      addreversemapping(facehashkeys[hashkey], volsmeshid);
      return -1;
    }
}

void VisMeshElementList::addFacesofSMESHElement(const SMDS_MeshElement * elem)
{
  SMDS_VolumeTool aTool; aTool.Set(elem);   
  int NbNodes = aTool.NbNodes(); int NbFaces = aTool.NbFaces();
  for (unsigned int facenum=0;  facenum < NbFaces; facenum++) 
    {
      int NbFaceNodes = aTool.NbFaceNodes(facenum);
      std::vector<int> facenodeids(NbFaceNodes, 0);
      const int *FaceIndices = aTool.GetFaceNodesIndices(facenum);
      for (unsigned int n=0; n<NbFaceNodes; n++) facenodeids[n]=elem->GetNode(FaceIndices[n])->GetID();
      addFaceofSMESHElement(facenodeids, elem->GetID());
    }  
}

void VisMeshElementList::addedges()
{
  int nbedges=0; 
  //add the edges
  SMDS_EdgeIteratorPtr 	anEdgeIter = mesh->GetMeshDS()->edgesIterator();
  for(;anEdgeIter->more();) {
    const SMDS_MeshEdge* anElem = anEdgeIter->next();
    if (submeshds!=NULL) 
      {
	if (submeshds->Contains(dynamic_cast<const SMDS_MeshElement *>(anElem)))	    
	  {
	    addSMESHElement(dynamic_cast<const SMDS_MeshElement *>(anElem)); nbedges++;
	  }
      }
    else 
      {
	addSMESHElement(dynamic_cast<const SMDS_MeshElement *>(anElem)); nbedges++;
      }
  }
}

double VisMeshElementList::ElemArea(const SMDS_MeshElement * elem)
{
  if (!elem) RETURN_BAD_RESULT("ElemArea: elem is null.");
  if (elem->GetType() != SMDSAbs_Face) RETURN_BAD_RESULT("ElemArea: face expected.");
  int NbNodes = elem->NbNodes(); 
  if ((NbNodes!=3) && (NbNodes!=4)) RETURN_BAD_RESULT("ElemArea: face with 3 or 4 nodes expected.");
  
  gp_XYZ p1 ( elem->GetNode(0)->X(), elem->GetNode(0)->Y(), elem->GetNode(0)->Z() );
  gp_XYZ p2 ( elem->GetNode(1)->X(), elem->GetNode(1)->Y(), elem->GetNode(1)->Z() );
  gp_XYZ p3 ( elem->GetNode(2)->X(), elem->GetNode(2)->Y(), elem->GetNode(2)->Z() );
  gp_XYZ aVec12( p2 - p1 );
  gp_XYZ aVec13( p3 - p1 );
  gp_XYZ aVectmp = aVec12.Crossed( aVec13 );
  double area = sqrt (aVectmp.X() * aVectmp.X() + aVectmp.Y() * aVectmp.Y() + aVectmp.Z() * aVectmp.Z())*0.5;
  if ( NbNodes == 4 ) {
    gp_XYZ p4 ( elem->GetNode(3)->X(), elem->GetNode(3)->Y(), elem->GetNode(3)->Z() );
    gp_XYZ aVec14( p4 - p1 );
    gp_XYZ aVectmp2 = aVec14.Crossed( aVec13 );
    area += sqrt (aVectmp2.X() * aVectmp2.X() + aVectmp2.Y() * aVectmp2.Y() + aVectmp2.Z() * aVectmp2.Z())*0.5;
  }
  return area; 
}

void VisMeshElementList::ComputeNodalFaceNormals()
{
  SMDS_FaceIteratorPtr 	aFaceIter = mesh->GetMeshDS()->facesIterator();
  std::map<int, std::vector<XYZ_> > nodefacenormals;
  std::map<int, std::vector<double> > nodefaceweights;
  for(;aFaceIter->more();) 
    {
      const SMDS_MeshFace* elem = aFaceIter->next();
      bool includeelem=true;
      if (submeshds!=NULL) includeelem=(submeshds->Contains(dynamic_cast<const SMDS_MeshElement *>(elem)));
      if (includeelem)
	{
	  SMDS_ElemIteratorPtr nodeIt = elem->nodesIterator();
	  while ( nodeIt->more() ) {
	    const SMDS_MeshNode* node = static_cast<const SMDS_MeshNode*>( nodeIt->next() );
	    if (node)
	      {
		int nodeid = node->GetID();
		if (nodefacenormals.find(nodeid)==nodefacenormals.end())
		  { 
		    std::vector<XYZ_> normals; std::vector<double> weights;
		    nodefacenormals.insert(std::make_pair(nodeid, normals)); 
		    nodefaceweights.insert(std::make_pair(nodeid, weights));
		  }		   
		nodefacenormals[nodeid].push_back(ElemNormals[elem->GetID()]);
		nodefaceweights[nodeid].push_back(ElemArea(elem));
	      }
	  }
	}
    }
  std::map<int, std::vector<XYZ_> >::const_iterator itr;
  for (itr=nodefacenormals.begin(); itr!=nodefacenormals.end(); ++itr)
    {
      const std::vector<XYZ_> & facenormals=itr->second;
      int nodeid=itr->first; int nbfacenormals=facenormals.size();
      const std::vector<double> & faceweights=nodefaceweights[nodeid];      
      XYZ_ N(0,0,0);
      for (int i=0; i<nbfacenormals; i++)
	{
	  N.x+=faceweights[i]*facenormals[i].x;
	  N.y+=faceweights[i]*facenormals[i].y;
	  N.z+=faceweights[i]*facenormals[i].z;
	}
      gp_XYZ N_(N.x, N.y, N.z); N_.Normalize();
      NodalNormals[nodeid]=XYZ_(N_.X(), N_.Y(), N_.Z());
    }
}

void VisMeshElementList::addfaces()
{
  //add the faces
  int nbfaces=0; 
  SMDS_FaceIteratorPtr 	aFaceIter = mesh->GetMeshDS()->facesIterator();
  for(;aFaceIter->more();) 
    {
      const SMDS_MeshFace* anElem = aFaceIter->next();
      bool includeelem=true;
      if (submeshds!=NULL) includeelem=(submeshds->Contains(dynamic_cast<const SMDS_MeshElement *>(anElem)));
      if (includeelem)
	{	
	  addSMESHElement(dynamic_cast<const SMDS_MeshElement *>(anElem)); nbfaces++;
	  ElemCentroids[anElem->GetID()]=ComputeCentroid(dynamic_cast<const SMDS_MeshElement *>(anElem));
	  double nx=0.0, ny=0.0, nz=0.0; ComputeNormal(anElem->GetID(), nx, ny, nz);
	  ElemNormals[anElem->GetID()] = XYZ_(nx, ny, nz);
	}
    }
  if (_supportsmoothshading) ComputeNodalFaceNormals();	    
}

void VisMeshElementList::addvolumes()
{
  //add the volumes 
  int nbvolumes = 0;
  SMDS_VolumeIteratorPtr aVolumeIter = mesh->GetMeshDS()->volumesIterator();
  for(;aVolumeIter->more();) 
    {
      const SMDS_MeshVolume* anElem = aVolumeIter->next();
      if (submeshds!=NULL) 
	{
	  if (submeshds->Contains(dynamic_cast<const SMDS_MeshElement *>(anElem)))	    
	    {
	      if (_displaymode==MeshVS_DMF_Shading)
		{
		  addFacesofSMESHElement(dynamic_cast<const SMDS_MeshElement *>(anElem)); ; nbvolumes++;
		  ElemCentroids[anElem->GetID()]=ComputeCentroid(dynamic_cast<const SMDS_MeshElement *>(anElem));		  
		}
	      else // if _displaymode==MeshVS_DMF_Shading
		{
		  addSMESHElement(dynamic_cast<const SMDS_MeshElement *>(anElem)); nbvolumes++;
		  ElemCentroids[anElem->GetID()]=ComputeCentroid(dynamic_cast<const SMDS_MeshElement *>(anElem));		  
		}
	    }
	} 	// if (submeshds!=NULL)
      else 
	{
	  addSMESHElement(dynamic_cast<const SMDS_MeshElement *>(anElem)); nbvolumes++;
	  ElemCentroids[anElem->GetID()]=ComputeCentroid(dynamic_cast<const SMDS_MeshElement *>(anElem));		  	
	}
    }
}


void VisMeshElementList::clear()
{
  reversemap.clear(); ElemCentroids.clear(); elemnum2idmap.clear();
  facehashkeys.clear(); elements.Clear(); volumefaces.clear(); 
  ElemNormals.clear(); NodalNormals.clear();
}

void VisMeshElementList::ComputeMeshElements()
{
  clear();
  addedges();
  addfaces();
  addvolumes(); 
}

void VisMeshElementList::SetDisplayMode(int displaymode)
{
  if (_displaymode != displaymode)
    { 
      _displaymode = displaymode; // 1:MeshVS_DMF_WireFrame, 2:MeshVS_DMF_Shading, 3:MeshVS_DMF_Shrink
      ComputeMeshElements(); //TODO: more fine-grained check if this is necessary
    }
}

std::size_t VisMeshElementList::TriHashKey(int nodeid1, int nodeid2, int nodeid3)
{
  int3[0]=nodeid1; int3[1]=nodeid2; int3[2]=nodeid3;
  std::sort(int3.begin(), int3.end()); 
  std::stringstream is; is << int3[0] << '/' << int3[1] << '/' << int3[2];  
  return string_hash(is.str());
}

std::size_t VisMeshElementList::QuadHashKey(int nodeid1, int nodeid2, int nodeid3, int nodeid4)
{
  int4[0]=nodeid1; int4[1]=nodeid2; int4[2]=nodeid3; int4[3]=nodeid4;
  std::sort(int4.begin(), int4.end()); 
  std::stringstream is; is << int4[0] << '/' << int4[1] << '/' << int4[2] << '/' << int4[3];
  return string_hash(is.str());
}

bool VisMeshElementList::IsValidElementNumber(int elemnumber) const
{
  if (elemnum2idmap.find(elemnumber)==elemnum2idmap.end())
    {
      MESSAGE("Element number "<< elemnumber <<" does not exist.");
      RETURN_BAD_RESULT("Invalid element number!");
    }
  return true;
}

const TColStd_PackedMapOfInteger& VisMeshElementList::GetAllElements()
{
  // returns a list of all element numbers to be displayed
  // the numbers are ids to be used with methods GetElement/GetElementNodes/GetElementType
  // for obtaining element geometry information 
  return elements;
}

bool VisMeshElementList::GetNodesByElement(int elemnumber, TColStd_Array1OfInteger& NodeIDs, int& NbNodes ) 
{
  if (!IsValidElementNumber(elemnumber)) RETURN_BAD_RESULT("GetNodesByElement: Invalid element number");
  std::pair<int, int> id_type = elemnum2idmap[elemnumber];
  if (id_type.second==FaceofSMESHElem)
    {
      std::vector<int>& nodeids = volumefaces[id_type.first]; NbNodes=nodeids.size();
      for (int i = 0; i < NbNodes; i++ ) NodeIDs.SetValue(i+1, nodeids[i]);	
      return true;
    }
  else
    {
      const SMDS_MeshElement* elem = mesh->GetMeshDS()->FindElement(id_type.first);
      if (!elem) RETURN_BAD_RESULT("GetNodesByElement: elem is NULL");
      NbNodes = elem->NbNodes();
      for (int i = 0; i < NbNodes; i++ ) 
	{
	  const SMDS_MeshNode* node = elem->GetNode(i);
	  if (!node) RETURN_BAD_RESULT("GetNodesByElement: node is NULL");
	  NodeIDs.SetValue(i+1, node->GetID());
	}
      return true;
    }
}

bool VisMeshElementList::GetNodeNormal(int nodenumber, int elemnumber, double &nx, double &ny, double &nz) 
{
  if (!IsValidElementNumber(elemnumber)) RETURN_BAD_RESULT("GetNodeNormal: Invalid element number");
  std::pair<int, int> id_type = elemnum2idmap[elemnumber];
  if (id_type.second==FaceofSMESHElem) return false; // smooth shading not supported for volume grids
  else
    {
      const SMDS_MeshElement* elem = mesh->GetMeshDS()->FindElement(id_type.first);
      if(!elem) RETURN_BAD_RESULT("GetNodeNormal: elem is NULL");
      const SMDS_MeshNode* node = elem->GetNode(nodenumber);
      if (!node) RETURN_BAD_RESULT("GetNodeNormal: node is NULL");
      if (NodalNormals.find(node->GetID())==NodalNormals.end()) RETURN_BAD_RESULT("GetNodeNormal: Nodal normal not computed.");
      XYZ_ N=NodalNormals[node->GetID()]; nx=N.x; ny=N.y; nz=N.z;
      return true;
    }
}

bool VisMeshElementList::ComputeNormal(int smeshelem, double& nx, double& ny,double& nz)
{
  const SMDS_MeshElement* elem = mesh->GetMeshDS()->FindElement(smeshelem);
  if(!elem) RETURN_BAD_RESULT("ComputeNormal: elem is NULL");
  std::vector<const SMDS_MeshNode*> nodes; 
  int NbNodes = elem->NbNodes();
  if ((NbNodes==3) || (NbNodes==4))
    {
      // note: first 3 nodes for quad faces sufficient
      for (int i = 0; i < 3; i++ ) nodes.push_back(elem->GetNode(i)); 	  
    }
  else RETURN_BAD_RESULT("GetNormal: Number of Nodes must be 3 or 4");
  gp_XYZ normal; gp_XYZ points[3];
  for (int i = 0; i < 3; i++) 
    {
      if (!nodes[i]) RETURN_BAD_RESULT("GetNormal: Node is NULL");
      points[i] = gp_XYZ(nodes[i]->X(), nodes[i]->Y(), nodes[i]->Z());
    }
  normal = (points[1]-points[0]) ^ (points[2]-points[0]);
  if ( normal.Modulus() > 0 ) normal /= normal.Modulus();
  nx = normal.X();
  ny = normal.Y();
  nz = normal.Z();
  return true;
}

bool VisMeshElementList::GetNormal(int elemnumber, double& nx, double& ny,double& nz )
{		
  if (!IsValidElementNumber(elemnumber)) RETURN_BAD_RESULT("GetNormal: Invalid element number");
  std::pair<int, int> id_type = elemnum2idmap[elemnumber];
  std::vector<const SMDS_MeshNode*> nodes; 
  if (id_type.second==FaceofSMESHElem)
    {
      std::vector<int>& nodeids = volumefaces[id_type.first];
      int NbNodes=nodeids.size();
      if ((NbNodes==3) || (NbNodes==4))
	{
	  // note: first 3 nodes for quad faces sufficient
	  for (int i = 0; i < 3; i++ ) nodes.push_back(mesh->GetMeshDS()->FindNode(nodeids[i])); 	  
	}
      else RETURN_BAD_RESULT("GetNormal: Number of nodes must be 3 or 4");
      gp_XYZ normal; gp_XYZ points[3];
      for (int i = 0; i < 3; i++) 
	{
	  if (!nodes[i]) RETURN_BAD_RESULT("GetNormal: Node is NULL");
	  points[i] = gp_XYZ(nodes[i]->X(), nodes[i]->Y(), nodes[i]->Z());
	}
      normal = (points[1]-points[0]) ^ (points[2]-points[0]);
      if ( normal.Modulus() > 0 ) normal /= normal.Modulus();
      nx = normal.X();
      ny = normal.Y();
      nz = normal.Z();
      return true;
    }
  else 
    {
      if (ElemNormals.find(id_type.first)==ElemNormals.end()) return false;
      XYZ_ N=ElemNormals[id_type.first]; nx=N.x; ny=N.y; nz=N.z;
      return true;
    }
}

bool VisMeshElementList::GetCentroid(int elemnumber, XYZ_& centroid)
{
  if (!IsValidElementNumber(elemnumber)) RETURN_BAD_RESULT("GetCentroid: Invalid element number");
  std::pair<int, int> id_type = elemnum2idmap[elemnumber];
  // no centroids are being computed for these elements because
  // their hidden state is being decided based on connected volume cells 
  if (id_type.second==FaceofSMESHElem) RETURN_BAD_RESULT("GetCentroid: Face of SMESH element!"); 
  if (ElemCentroids.find(id_type.first)==ElemCentroids.end()) RETURN_BAD_RESULT("GetCentroid: Centroid was not computed");
  centroid=ElemCentroids[id_type.first]; return true;
}

XYZ_ VisMeshElementList::ComputeCentroid(const SMDS_MeshElement* elem)
{  
  XYZ_ centroid(0,0,0); int count=0;
  SMDS_ElemIteratorPtr nodeIt = elem->nodesIterator();
  while ( nodeIt->more() ) {
    const SMDS_MeshNode* node = static_cast<const SMDS_MeshNode*>( nodeIt->next() );
    if (node)
      {
	centroid.x += node->X(); centroid.y += node->Y(); centroid.z += node->Z(); count++;
      }
    else MESSAGE("ComputeCentroid: node not found.");
  }
  centroid.x /= count; centroid.y /= count; centroid.z /= count;
  return centroid;
}

bool VisMeshElementList::GetFaceofSMESHElement(int volumefaceid, Standard_Integer& NbNodes, Handle(MeshVS_HArray1OfSequenceOfInteger)& Data)
{
  if ((volumefaceid<0) || (volumefaceid>=volumefaces.size())) RETURN_BAD_RESULT("Invalid volume face id");
  int NbFaces = 1; NbNodes = volumefaces[volumefaceid].size();
  if (Data.IsNull())
    Data = new MeshVS_HArray1OfSequenceOfInteger(1, NbFaces);
  else if (Data->Length() != NbFaces) 
    {
      Data.Nullify();
      Data = new MeshVS_HArray1OfSequenceOfInteger(1, NbFaces);
    }
  TColStd_SequenceOfInteger aSeq;
  for (unsigned int v=0; v<volumefaces[volumefaceid].size(); v++)
    aSeq.Append(volumefaces[volumefaceid][v]); 
  Data->SetValue(1, aSeq);
  return true;
}

bool VisMeshElementList::GetSMESHElement(int smeshid, Standard_Integer& NbNodes, Handle(MeshVS_HArray1OfSequenceOfInteger)& Data)
{
  const SMDS_MeshElement* elem = mesh->GetMeshDS()->FindElement(smeshid);
  if (!elem) RETURN_BAD_RESULT("GetSMESHElement: elem is NULL");
  if (elem->GetType() != SMDSAbs_Volume) return false;
  //initialize VolumeTool
  SMDS_VolumeTool aTool;
  aTool.Set(elem);
  //set the nodes number
  NbNodes = aTool.NbNodes();
  //check validity or create Data
  int NbFaces = aTool.NbFaces();
  if (Data.IsNull())
    Data = new MeshVS_HArray1OfSequenceOfInteger(1, NbFaces);
  else if (Data->Length() != NbFaces) {
    Data.Nullify();
    Data = new MeshVS_HArray1OfSequenceOfInteger(1, NbFaces);
  }
  //iterate the faces and their nodes and add them to Data
  for (int itr=0;itr < NbFaces;itr++) {
    int NbThisFaceNodeCount = aTool.NbFaceNodes(itr);
    const int *FaceIndices = aTool.GetFaceNodesIndices(itr);
    int sortedFaceIndices[MAX_SORT_NODE_COUNT];
    TColStd_SequenceOfInteger aSeq;
    if (sortNodes2(elem, FaceIndices, NbThisFaceNodeCount, sortedFaceIndices)) {
      for (int itrX=0;itrX < NbThisFaceNodeCount;itrX++)
	aSeq.Append(sortedFaceIndices[itrX]);
    } else {
      for (int itrX=0;itrX < NbThisFaceNodeCount;itrX++)
	aSeq.Append(FaceIndices[itrX]);
    }
    Data->SetValue(itr+1, aSeq);
  }
  return true;
}

bool VisMeshElementList::GetElement(int elemnumber, Standard_Integer& NbNodes, Handle(MeshVS_HArray1OfSequenceOfInteger)& Data)
{
  if (!IsValidElementNumber(elemnumber)) RETURN_BAD_RESULT("GetElement: Invalid element number");
  std::pair<int, int> id_type = elemnum2idmap[elemnumber];
  if (id_type.second==FaceofSMESHElem)
    return GetFaceofSMESHElement(id_type.first, NbNodes, Data);
  else
    return GetSMESHElement(id_type.first, NbNodes, Data);  
}

bool VisMeshElementList::GetElementType(int elemnumber, MeshVS_EntityType& elementtype) 
{
  if (!IsValidElementNumber(elemnumber)) RETURN_BAD_RESULT("GetElementType: Invalid element number");
  std::pair<int, int> id_type = elemnum2idmap[elemnumber];
  if (id_type.second==FaceofSMESHElem) {elementtype=MeshVS_ET_Face; return true;}
  else
  {
    const SMDS_MeshElement* elem = mesh->GetMeshDS()->FindElement(id_type.first);
    if (!elem) RETURN_BAD_RESULT("GetElementType: Elem is NULL");
    if (elem->GetType() == SMDSAbs_Edge)
      elementtype = MeshVS_ET_Link;
    else if (elem->GetType() == SMDSAbs_Face)
      elementtype = MeshVS_ET_Face;
    else if (elem->GetType() == SMDSAbs_Volume)
      elementtype = MeshVS_ET_Volume;
    else
      elementtype = MeshVS_ET_Element;
    return true;
  }
}

bool VisMeshElementList::GetElementNodes(int elemnumber, TColStd_Array1OfReal& Coords, Standard_Integer& NbNodes) 
{
  if (!IsValidElementNumber(elemnumber)) RETURN_BAD_RESULT("GetElementNodes: Invalid element number");
  std::pair<int, int> id_type = elemnum2idmap[elemnumber];
  if (id_type.second==FaceofSMESHElem)
    {
      std::vector<int>& nodeids = volumefaces[id_type.first];
      int nbCoord = 1;  NbNodes = nodeids.size();
      for(int i=0; i<NbNodes; i++)
	{
	  const SMDS_MeshNode * node = mesh->GetMeshDS()->FindNode(nodeids[i]);
	  if (!node) RETURN_BAD_RESULT("GetElementNodes: Node is NULL");
	  Coords(nbCoord++) = node->X();
	  Coords(nbCoord++) = node->Y();
	  Coords(nbCoord++) = node->Z();	  
	}
      return true;
    }
  else
    {
      const SMDS_MeshElement* elem = mesh->GetMeshDS()->FindElement(id_type.first);
      if (!elem) RETURN_BAD_RESULT("GetElementNodes: Elem is NULL");
      NbNodes = elem->NbNodes();
      int nbCoord = 1; 
      for(int i = 0; i < NbNodes; i++ ) 
	{
	  Coords(nbCoord++) = elem->GetNode(i)->X();
	  Coords(nbCoord++) = elem->GetNode(i)->Y();
	  Coords(nbCoord++) = elem->GetNode(i)->Z();
	}
      return true;
    }
}

void VisMeshElementList::GetHiddenElements(int planedir, double planevalue, TColStd_PackedMapOfInteger& hiddenelements) 
{
  hiddenelements.Clear();
  TColStd_MapIteratorOfPackedMapOfInteger anIter (elements);
  for ( ; anIter.More(); anIter.Next() )
    {
      Standard_Integer elemnumber = anIter.Key();
      if (IsHidden(elemnumber, planedir, planevalue)) hiddenelements.Add(elemnumber);
    }
}

bool VisMeshElementList::IsHidden(int elemnumber, int planedir, double planevalue)
{
  if (!IsValidElementNumber(elemnumber)) RETURN_BAD_RESULT("IsHidden: Invalid element number");
  std::pair<int, int> id_type = elemnum2idmap[elemnumber];
  if (id_type.second==FaceofSMESHElem)
    {
      std::set<int> & volelems = reversemap[id_type.first]; 
      int nbdisplayedvolelems = 0;
      int nbconncectedvolelems = volelems.size();
      for (std::set<int>::iterator i = volelems.begin(); i!=volelems.end(); ++i)
	if (!IsElemHidden(*i, planedir, planevalue)) nbdisplayedvolelems++;    
      // 2 possible cases for which this face must be displayed:
      // 1) face is on boundary: nbconncectedvolelems=1 and that volume element is displayed
      // 2) face is on boundary introduced by cutplane:
      //    nbconncectedvolelems=1 or 2, but among them only one is displayed
      // --> criterium for 1) and 2): nbdisplayedvolelems==1
      return (nbdisplayedvolelems!=1); 
    }
  else return IsElemHidden(id_type.first, planedir, planevalue);
}

bool VisMeshElementList::IsElemHidden(int smeshelemid, int planedir, double planevalue)
{
  if (ElemCentroids.find(smeshelemid)==ElemCentroids.end()) RETURN_BAD_RESULT("IsElemHidden: Centroid for element not found!");
  XYZ_ centroid=ElemCentroids[smeshelemid];
  if (planedir==1)
    {
      if (centroid.x < planevalue) return true;
    }
  else if (planedir==2)
    {
      if (centroid.y < planevalue) return true;
    }
  else if (planedir==3)
    {
      if (centroid.z < planevalue) return true;
    }
  else if (planedir==4)
    {
      if (centroid.x > planevalue) return true;
    }
  else if (planedir==5)
    {
      if (centroid.y > planevalue) return true;
    }
  else if (planedir==6)
    {
      if (centroid.z > planevalue) return true;
    }
  else return false; // cutplane deactivated (planedir==0)
  return false;     
}


SMESH_MeshVSLink2::SMESH_MeshVSLink2(const SMESH_Mesh *aMesh, SMESH_subMesh *aSubMesh, int displaymode)
  : myElements(aMesh, displaymode, aSubMesh), mesh(aMesh), submesh(aSubMesh), submeshds(NULL), _displaymode(displaymode)
{
  // displaymode: 0:MeshVS_DMF_Shading, 1:MeshVS_DMF_WireFrame, 2:MeshVS_DMF_Shrink  
  if (submesh) submeshds = submesh->CreateSubMeshDS();   

  //add the nodes
  int nbnodes=0;
  SMDS_NodeIteratorPtr aNodeIter = mesh->GetMeshDS()->nodesIterator();
  for(;aNodeIter->more();) {
    const SMDS_MeshNode* aNode = aNodeIter->next();
    if (submeshds!=NULL) 
      {
	if (submeshds->Contains(dynamic_cast<const SMDS_MeshElement *>(aNode)))
	  {
	    myNodes.Add( aNode->GetID() ); nbnodes++;
	  }
      }
    else myNodes.Add( aNode->GetID() );
  }
  
  //add the groups
  const std::set<SMESHDS_GroupBase*>& groups = mesh->GetMeshDS()->GetGroups();
  if (!groups.empty()) {
    std::set<SMESHDS_GroupBase*>::const_iterator GrIt = groups.begin();
    for (; GrIt != groups.end(); GrIt++) {
      SMESHDS_Group* grp = dynamic_cast<SMESHDS_Group*>(*GrIt);
      if (!grp || grp->IsEmpty()) continue;
      myGroups.Add(grp->GetID());
    }
  }
}

void SMESH_MeshVSLink2::SetDisplayMode(int displaymode)
{
  _displaymode = displaymode; // 1:MeshVS_DMF_WireFrame, 2:MeshVS_DMF_Shading, 3:MeshVS_DMF_Shrink
  myElements.SetDisplayMode(displaymode);
}


void SMESH_MeshVSLink2::GetHiddenNodes(int planedir, double value, TColStd_PackedMapOfInteger& hiddennodes) 
{  
  TColStd_MapIteratorOfPackedMapOfInteger anIter (myNodes);
  for ( ; anIter.More(); anIter.Next() )
    {
      Standard_Integer aKey = anIter.Key();
      // get the node 
      const SMDS_MeshNode* myNode = mesh->GetMeshDS()->FindNode(aKey);
      if (myNode)
	{
	  if (planedir==1)
	    {
	      if (myNode->X() < value) hiddennodes.Add(aKey);
	    }
	  else if (planedir==2)
	    {
	      if (myNode->Y() < value) hiddennodes.Add(aKey);
	    }
	  else if (planedir==3)
	    {
	      if (myNode->Z() < value) hiddennodes.Add(aKey);
	    }
	  else if (planedir==4)
	    {
	      if (myNode->X() > value) hiddennodes.Add(aKey);
	    }
	  else if (planedir==5)
	    {
	      if (myNode->Y() > value) hiddennodes.Add(aKey);
	    }
	  else if (planedir==6)
	    {
	      if (myNode->Z() > value) hiddennodes.Add(aKey);
	    }
	}
    }
}

void SMESH_MeshVSLink2::GetHiddenElements(int planedir, double planevalue, TColStd_PackedMapOfInteger& hiddenelements) 
{
  myElements.GetHiddenElements(planedir, planevalue, hiddenelements);
}



Standard_Boolean SMESH_MeshVSLink2::GetGeom(const Standard_Integer ID, const Standard_Boolean IsElement,   
					    TColStd_Array1OfReal& Coords, Standard_Integer& NbNodes,
					    MeshVS_EntityType& Type) const
{
  if( IsElement )
    {
      if (!myElements.GetElementType(ID, Type) ) return false;
      return myElements.GetElementNodes(ID, Coords, NbNodes);
    }
  else 
    {
      const SMDS_MeshNode* myNode = mesh->GetMeshDS()->FindNode(ID);
      if (!myNode) return Standard_False;
      if (myNode->GetType() == SMDSAbs_Node)
	Type = MeshVS_ET_Node;
      else
	Type = MeshVS_ET_0D;
      NbNodes = 1;
      Coords(1) = myNode->X();
      Coords(2) = myNode->Y();
      Coords(3) = myNode->Z();
      return Standard_True;
    }
}

Standard_Boolean  SMESH_MeshVSLink2::Get3DGeom(const Standard_Integer ID, Standard_Integer& NbNodes,
					       Handle(MeshVS_HArray1OfSequenceOfInteger)& Data) const
{  
  return myElements.GetElement(ID, NbNodes, Data);
}



Standard_Boolean SMESH_MeshVSLink2::GetGeomType(const Standard_Integer ID, const Standard_Boolean IsElement,						
						MeshVS_EntityType& Type) const
{
  if( IsElement ) return myElements.GetElementType(ID, Type);
  else 
    {
      const SMDS_MeshNode* myNode = mesh->GetMeshDS()->FindNode(ID);
      if (!myNode) return Standard_False;
      if (myNode->GetType() == SMDSAbs_Node)
	Type = MeshVS_ET_Node;
      else
	Type = MeshVS_ET_0D;
      return Standard_True;
    }
}


Standard_Boolean SMESH_MeshVSLink2::GetNodesByElement( const Standard_Integer ID,TColStd_Array1OfInteger& NodeIDs,Standard_Integer& NbNodes ) const
{
  return myElements.GetNodesByElement(ID, NodeIDs, NbNodes);
}

const TColStd_PackedMapOfInteger& SMESH_MeshVSLink2::GetAllNodes() const
{
  return myNodes;
}

const TColStd_PackedMapOfInteger& SMESH_MeshVSLink2::GetAllElements() const
{
  return myElements.GetAllElements();
}

void SMESH_MeshVSLink2::GetAllGroups(TColStd_PackedMapOfInteger& Ids) const
{
  Ids = myGroups;
}

Standard_Address SMESH_MeshVSLink2::GetAddr( const Standard_Integer, const Standard_Boolean ) const
{
  return NULL;
}

Standard_Boolean SMESH_MeshVSLink2::GetNodeNormal( const Standard_Integer ranknode, const Standard_Integer Id, Standard_Real &nx, Standard_Real &ny, Standard_Real &nz) const
{
  return myElements.GetNodeNormal(ranknode-1, Id, nx, ny, nz);
}

Standard_Boolean SMESH_MeshVSLink2::GetNormal( const Standard_Integer Id, const Standard_Integer Max,
					       Standard_Real& nx, Standard_Real& ny,Standard_Real& nz ) const
{
  if(Max<3) return Standard_False;
  return myElements.GetNormal(Id, nx, ny, nz);
}
