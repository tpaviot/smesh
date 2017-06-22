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
//  File   : SMESH_DeviceActor.cxx
//  Author : 
//  Module : SMESH
//
#include "SMESH_DeviceActor.h"
#include "SMESH_ScalarBarActor.h"
#include "SMESH_ExtractGeometry.h"
#include "SMESH_ControlsDef.hxx"
#include "SMESH_ActorUtils.h"
#include "SMESH_FaceOrientationFilter.h"
#include "VTKViewer_CellLocationsArray.h"
#include "VTKViewer_PolyDataMapper.h"

#include <VTKViewer_Transform.h>
#include <VTKViewer_TransformFilter.h>
#include <VTKViewer_ExtractUnstructuredGrid.h>

// VTK Includes
#include <vtkObjectFactory.h>
#include <vtkShrinkFilter.h>
#include <vtkShrinkPolyData.h>

#include <vtkProperty.h>
#include <vtkPolyData.h>
#include <vtkMergeFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkUnstructuredGrid.h>

#include <vtkLookupTable.h>
#include <vtkDoubleArray.h>
#include <vtkCellData.h>

#include <vtkCell.h>
#include <vtkIdList.h>
#include <vtkCellArray.h>
#include <vtkUnsignedCharArray.h>

#include <vtkImplicitBoolean.h>
#include <vtkPassThroughFilter.h>

#include <vtkRenderer.h>

#include <vtkPlaneCollection.h>

#include "utilities.h"

#ifdef _DEBUG_
static int MYDEBUG = 0;
#else
static int MYDEBUG = 0;
#endif

using namespace std;


vtkStandardNewMacro(SMESH_DeviceActor);


SMESH_DeviceActor
::SMESH_DeviceActor()
{
  if(MYDEBUG) MESSAGE("SMESH_DeviceActor - "<<this);

  myIsShrinkable = false;
  myIsShrunk = false;
  myIsHighlited = false;

  myRepresentation = SMESH_DeviceActor::EReperesent(-1);

  myProperty = vtkProperty::New();
  myMapper = VTKViewer_PolyDataMapper::New();
  myPlaneCollection = vtkPlaneCollection::New();

  vtkMapper::GetResolveCoincidentTopologyPolygonOffsetParameters(myPolygonOffsetFactor,
                                                                 myPolygonOffsetUnits);

  myMapper->UseLookupTableScalarRangeOn();
  myMapper->SetColorModeToMapScalars();

  myShrinkFilter = vtkShrinkFilter::New();

  myStoreClippingMapping = false;

  myExtractGeometry = SMESH_ExtractGeometry::New();
  myExtractGeometry->SetReleaseDataFlag(true);
  myIsImplicitFunctionUsed = false;

  myExtractUnstructuredGrid = VTKViewer_ExtractUnstructuredGrid::New();
    
  myMergeFilter = vtkMergeFilter::New();

  myGeomFilter = VTKViewer_GeometryFilter::New();

  myTransformFilter = VTKViewer_TransformFilter::New();

  for(int i = 0; i < 6; i++)
    myPassFilter.push_back(vtkPassThroughFilter::New());

  // Orientation of faces
  myIsFacesOriented = false;

  double anRGB[3] = { 1, 1, 1 };
  SMESH::GetColor( "SMESH", "orientation_color", anRGB[0], anRGB[1], anRGB[2], QColor( 255, 255, 255 ) );

  myFaceOrientationFilter = SMESH_FaceOrientationFilter::New();

  myFaceOrientationDataMapper = vtkPolyDataMapper::New();
  myFaceOrientationDataMapper->SetInputConnection(myFaceOrientationFilter->GetOutputPort());

  myFaceOrientation = vtkActor::New();
  myFaceOrientation->SetMapper(myFaceOrientationDataMapper);
  myFaceOrientation->GetProperty()->SetColor(anRGB[0], anRGB[1], anRGB[2]);
}


SMESH_DeviceActor
::~SMESH_DeviceActor()
{
  if(MYDEBUG) MESSAGE("~SMESH_DeviceActor - "<<this);

  myMapper->Delete();
  // myPlaneCollection->Delete(); -- it is vtkSmartPointer
  myProperty->Delete();

  myExtractGeometry->Delete();

  myMergeFilter->Delete();
  myExtractUnstructuredGrid->Delete();

  // Orientation of faces
  myFaceOrientationFilter->Delete();
  myFaceOrientationDataMapper->RemoveAllInputs();
  myFaceOrientationDataMapper->Delete();
  myFaceOrientation->Delete();

  myGeomFilter->Delete();

  myTransformFilter->Delete();

  for(int i = 0, iEnd = myPassFilter.size(); i < iEnd; i++)
    myPassFilter[i]->Delete();

  myShrinkFilter->Delete();
}


void
SMESH_DeviceActor
::SetStoreGemetryMapping(bool theStoreMapping)
{
  myGeomFilter->SetStoreMapping(theStoreMapping);
  // for optimization, switch the mapping explicitly in each filter/algorithm
  //SetStoreClippingMapping(theStoreMapping);
}


void
SMESH_DeviceActor
::SetStoreClippingMapping(bool theStoreMapping)
{
  myStoreClippingMapping = theStoreMapping;
  myExtractGeometry->SetStoreMapping(theStoreMapping && myIsImplicitFunctionUsed);
  // EAP, 23315
  // Mapping in myExtractUnstructuredGrid and myGeomFilter is ON in the pickable DeviceActor only.
  // To show labels, the mapping is computed explicitly via myExtractUnstructuredGrid->BuildOut2InMap();
  //SetStoreIDMapping(theStoreMapping);
}


void
SMESH_DeviceActor
::SetStoreIDMapping(bool theStoreMapping)
{
  myExtractUnstructuredGrid->SetStoreMapping(theStoreMapping);
}


void 
SMESH_DeviceActor
::Init(TVisualObjPtr theVisualObj, 
       vtkImplicitBoolean* theImplicitBoolean)
{
  myVisualObj = theVisualObj;
  myExtractGeometry->SetImplicitFunction(theImplicitBoolean);
  SetUnstructuredGrid(myVisualObj->GetUnstructuredGrid());
}

void
SMESH_DeviceActor
::SetImplicitFunctionUsed(bool theIsImplicitFunctionUsed)
{
  int anId = 0;
  if(theIsImplicitFunctionUsed)
    myPassFilter[ anId ]->SetInputConnection( myExtractGeometry->GetOutputPort() );
  else
    myPassFilter[ anId ]->SetInputConnection( myMergeFilter->GetOutputPort() );
    
  myIsImplicitFunctionUsed = theIsImplicitFunctionUsed;
  SetStoreClippingMapping(myStoreClippingMapping);
}


void
SMESH_DeviceActor
::SetUnstructuredGrid(vtkUnstructuredGrid* theGrid)
{
  myExtractUnstructuredGrid->SetInputData(theGrid);

  if ( theGrid )
  {
    myIsShrinkable = true;

    myMergeFilter->SetGeometryConnection(myExtractUnstructuredGrid->GetOutputPort());

    //Pass diameters of the balls
    if(myMapper->GetBallEnabled()) {
      myMergeFilter->SetScalarsConnection(myExtractUnstructuredGrid->GetOutputPort());
    }

    myExtractGeometry->SetInputConnection(myMergeFilter->GetOutputPort());

    int anId = 0;
    SetImplicitFunctionUsed(myIsImplicitFunctionUsed);
    myPassFilter[ anId + 1]->SetInputConnection( myPassFilter[ anId ]->GetOutputPort() );

    anId++; // 1
    myTransformFilter->SetInputConnection( myPassFilter[ anId ]->GetOutputPort() );

    anId++; // 2
    myPassFilter[ anId ]->SetInputConnection( myTransformFilter->GetOutputPort() );
    myPassFilter[ anId + 1 ]->SetInputConnection( myPassFilter[ anId ]->GetOutputPort() );

    anId++; // 3
    myGeomFilter->SetInputConnection( myPassFilter[ anId ]->GetOutputPort() );

    anId++; // 4
    myPassFilter[ anId ]->SetInputConnection( myGeomFilter->GetOutputPort() );
    myPassFilter[ anId + 1 ]->SetInputConnection( myPassFilter[ anId ]->GetOutputPort() );

    anId++; // 5
    myMapper->SetInputConnection( myPassFilter[ anId ]->GetOutputPort() );
    if( myPlaneCollection->GetNumberOfItems() )
      myMapper->SetClippingPlanes( myPlaneCollection );

    vtkLODActor::SetMapper( myMapper );
  }
  Modified();
}

void
SMESH_DeviceActor
::SetPlaneCollection( vtkPlaneCollection* theCollection )
{
  myPlaneCollection = theCollection;
}

VTKViewer_ExtractUnstructuredGrid* 
SMESH_DeviceActor
::GetExtractUnstructuredGrid()
{
  return myExtractUnstructuredGrid;
}

#include "SMDS_Mesh.hxx"

vtkUnstructuredGrid* 
SMESH_DeviceActor
::GetUnstructuredGrid()
{
  myExtractUnstructuredGrid->Update();
  return myExtractUnstructuredGrid->GetOutput();
}


void
SMESH_DeviceActor
::SetControlMode(SMESH::Controls::FunctorPtr theFunctor,
                 SMESH_ScalarBarActor* theScalarBarActor,
                 vtkLookupTable* theLookupTable)
{
  bool anIsInitialized = theFunctor != NULL;
  if(anIsInitialized){
    vtkUnstructuredGrid* aDataSet = vtkUnstructuredGrid::New();

    // SetStoreIDMapping(true);
    myExtractUnstructuredGrid->Update();
    vtkUnstructuredGrid* aGrid = myExtractUnstructuredGrid->GetOutput();

    aDataSet->ShallowCopy(aGrid);
    
    vtkDoubleArray *aScalars = vtkDoubleArray::New();
    vtkIdType aNbCells = aGrid->GetNumberOfCells();
    aScalars->SetNumberOfComponents(1);
    aScalars->SetNumberOfTuples(aNbCells);
    double* range = 0;// = aScalars->GetRange();
    
    myVisualObj->UpdateFunctor(theFunctor);

    using namespace SMESH::Controls;
    if(NumericalFunctor* aNumericalFunctor = dynamic_cast<NumericalFunctor*>(theFunctor.get()))
    {
      myExtractUnstructuredGrid->BuildOut2InMap();
      for(vtkIdType i = 0; i < aNbCells; i++)
      {
        vtkIdType anId = myExtractUnstructuredGrid->GetInputId(i);
        vtkIdType anObjId = myVisualObj->GetElemObjId(anId);
        double aValue = aNumericalFunctor->GetValue(anObjId);
        aScalars->SetValue(i,aValue);
      }
      range = aScalars->GetRange();
      if ( range[1] - range[0] < ( qMax(qAbs(range[0]),qAbs(range[1])) + 1e-100 ) * 1e-6 )
      {
        range[1] = range[0];
        for(vtkIdType i = 0; i < aNbCells; i++)
          aScalars->SetValue(i,range[0]);
      }
    }
    else if(Predicate* aPredicate = dynamic_cast<Predicate*>(theFunctor.get()))
    {
      myExtractUnstructuredGrid->BuildOut2InMap();
      for(vtkIdType i = 0; i < aNbCells; i++)
      {
        vtkIdType anId = myExtractUnstructuredGrid->GetInputId(i);
        vtkIdType anObjId = myVisualObj->GetElemObjId(anId);
        bool aValue = aPredicate->IsSatisfy(anObjId);
        aScalars->SetValue(i,aValue);
      }
      range = aScalars->GetRange();
    }

    aDataSet->GetCellData()->SetScalars(aScalars);
    aScalars->Delete();

    theLookupTable->SetRange( range );
    theLookupTable->SetNumberOfTableValues(theScalarBarActor->GetMaximumNumberOfColors());
    theLookupTable->Build();
    
    myMergeFilter->SetScalarsData(aDataSet);
    aDataSet->Delete();
  }
  GetMapper()->SetScalarVisibility(anIsInitialized);
  theScalarBarActor->SetVisibility(anIsInitialized);
}

void
SMESH_DeviceActor
::SetExtControlMode(SMESH::Controls::FunctorPtr theFunctor,
                    SMESH_ScalarBarActor* theScalarBarActor,
                    vtkLookupTable* theLookupTable)
{
  bool anIsInitialized = theFunctor != NULL;
  myExtractUnstructuredGrid->ClearRegisteredCells();
  myExtractUnstructuredGrid->ClearRegisteredCellsWithType();
  myExtractUnstructuredGrid->SetModeOfChanging(VTKViewer_ExtractUnstructuredGrid::ePassAll);
  myVisualObj->UpdateFunctor(theFunctor);

  using namespace SMESH::Controls;
  if (anIsInitialized){
    if (Length2D* aLength2D = dynamic_cast<Length2D*>(theFunctor.get())){
      SMESH::Controls::Length2D::TValues aValues;

      aLength2D->GetValues(aValues);
      vtkUnstructuredGrid* aDataSet = vtkUnstructuredGrid::New();
      vtkUnstructuredGrid* aGrid = myVisualObj->GetUnstructuredGrid();

      aDataSet->SetPoints(aGrid->GetPoints());
      
      vtkIdType aNbCells = aValues.size();
      
      vtkDoubleArray *aScalars = vtkDoubleArray::New();
      aScalars->SetNumberOfComponents(1);
      aScalars->SetNumberOfTuples(aNbCells);

      vtkIdType aCellsSize = 3*aNbCells;
      vtkCellArray* aConnectivity = vtkCellArray::New();
      aConnectivity->Allocate( aCellsSize, 0 );
      
      vtkUnsignedCharArray* aCellTypesArray = vtkUnsignedCharArray::New();
      aCellTypesArray->SetNumberOfComponents( 1 );
      aCellTypesArray->Allocate( aNbCells * aCellTypesArray->GetNumberOfComponents() );

      vtkIdList *anIdList = vtkIdList::New();
      anIdList->SetNumberOfIds(2);

      Length2D::TValues::const_iterator anIter = aValues.begin();
      aNbCells = 0;
      for(; anIter != aValues.end(); anIter++){
        const Length2D::Value& aValue = *anIter;
        int aNode[2] = {
          myVisualObj->GetNodeVTKId(aValue.myPntId[0]),
          myVisualObj->GetNodeVTKId(aValue.myPntId[1])
        };
        if(aNode[0] >= 0 && aNode[1] >= 0){
          anIdList->SetId( 0, aNode[0] );
          anIdList->SetId( 1, aNode[1] );
          aConnectivity->InsertNextCell( anIdList );
          aCellTypesArray->InsertNextValue( VTK_LINE );
          aScalars->SetValue(aNbCells,aValue.myLength);
          aNbCells++;
        }
      }
      aCellTypesArray->SetNumberOfTuples( aNbCells );
      aScalars->SetNumberOfTuples( aNbCells );

      VTKViewer_CellLocationsArray* aCellLocationsArray = VTKViewer_CellLocationsArray::New();
      aCellLocationsArray->SetNumberOfComponents( 1 );
      aCellLocationsArray->SetNumberOfTuples( aNbCells );

      aConnectivity->InitTraversal();
      for( vtkIdType idType = 0, *pts, npts; aConnectivity->GetNextCell( npts, pts ); idType++ )
        aCellLocationsArray->SetValue( idType, aConnectivity->GetTraversalLocation( npts ) );

      aDataSet->SetCells( aCellTypesArray, aCellLocationsArray, aConnectivity );
      SetUnstructuredGrid(aDataSet);

      aDataSet->GetCellData()->SetScalars(aScalars);
      aScalars->Delete();

      theLookupTable->SetRange(aScalars->GetRange());
      theLookupTable->Build();

      myMergeFilter->SetScalarsData(aDataSet);
      aDataSet->Delete();
    }
    else if (MultiConnection2D* aMultiConnection2D = dynamic_cast<MultiConnection2D*>(theFunctor.get())){
      SMESH::Controls::MultiConnection2D::MValues aValues;

      aMultiConnection2D->GetValues(aValues);
      vtkUnstructuredGrid* aDataSet = vtkUnstructuredGrid::New();
      vtkUnstructuredGrid* aGrid = myVisualObj->GetUnstructuredGrid();
      aDataSet->SetPoints(aGrid->GetPoints());
      
      vtkIdType aNbCells = aValues.size();
      vtkDoubleArray *aScalars = vtkDoubleArray::New();
      aScalars->SetNumberOfComponents(1);
      aScalars->SetNumberOfTuples(aNbCells);

      vtkIdType aCellsSize = 3*aNbCells;
      vtkCellArray* aConnectivity = vtkCellArray::New();
      aConnectivity->Allocate( aCellsSize, 0 );

      vtkUnsignedCharArray* aCellTypesArray = vtkUnsignedCharArray::New();
      aCellTypesArray->SetNumberOfComponents( 1 );
      aCellTypesArray->Allocate( aNbCells * aCellTypesArray->GetNumberOfComponents() );

      vtkIdList *anIdList = vtkIdList::New();
      anIdList->SetNumberOfIds(2);

      MultiConnection2D::MValues::const_iterator anIter = aValues.begin();
      aNbCells = 0;
      for(; anIter != aValues.end(); anIter++){
        const MultiConnection2D::Value& aValue = (*anIter).first;
        int aNode[2] = {
          myVisualObj->GetNodeVTKId(aValue.myPntId[0]),
          myVisualObj->GetNodeVTKId(aValue.myPntId[1])
        };
        if(aNode[0] >= 0 && aNode[1] >= 0){
          anIdList->SetId( 0, aNode[0] );
          anIdList->SetId( 1, aNode[1] );
          aConnectivity->InsertNextCell( anIdList );
          aCellTypesArray->InsertNextValue( VTK_LINE );
          aScalars->SetValue( aNbCells,(*anIter).second);
          aNbCells++;
        }
      }
      aCellTypesArray->SetNumberOfTuples( aNbCells );
      aScalars->SetNumberOfTuples( aNbCells );

      VTKViewer_CellLocationsArray* aCellLocationsArray = VTKViewer_CellLocationsArray::New();
      aCellLocationsArray->SetNumberOfComponents( 1 );
      aCellLocationsArray->SetNumberOfTuples( aNbCells );

      aConnectivity->InitTraversal();
      for( vtkIdType idType = 0, *pts, npts; aConnectivity->GetNextCell( npts, pts ); idType++ )
        aCellLocationsArray->SetValue( idType, aConnectivity->GetTraversalLocation( npts ) );

      aDataSet->SetCells( aCellTypesArray, aCellLocationsArray,aConnectivity );
      SetUnstructuredGrid(aDataSet);

      aDataSet->GetCellData()->SetScalars(aScalars);
      aScalars->Delete();

      theLookupTable->SetRange(aScalars->GetRange());
      theLookupTable->Build();

      myMergeFilter->SetScalarsData(aDataSet);
      aDataSet->Delete();
    }
  }
  GetMapper()->SetScalarVisibility(anIsInitialized);
  theScalarBarActor->SetVisibility(anIsInitialized);
}

void
SMESH_DeviceActor
::SetExtControlMode(SMESH::Controls::FunctorPtr theFunctor)
{
  myExtractUnstructuredGrid->ClearRegisteredCells();
  myExtractUnstructuredGrid->ClearRegisteredCellsWithType();
  myExtractUnstructuredGrid->SetModeOfChanging(VTKViewer_ExtractUnstructuredGrid::ePassAll);
  myVisualObj->UpdateFunctor(theFunctor);

  using namespace SMESH::Controls;
  Predicate* aPredicate = 0;
  if (( aPredicate =  dynamic_cast<FreeBorders          *>(theFunctor.get())) ||
      ( aPredicate =  dynamic_cast<FreeFaces            *>(theFunctor.get())) ||
      ( aPredicate =  dynamic_cast<BareBorderVolume     *>(theFunctor.get())) ||
      ( aPredicate =  dynamic_cast<BareBorderFace       *>(theFunctor.get())) ||
      ( aPredicate =  dynamic_cast<OverConstrainedVolume*>(theFunctor.get())) ||
      ( aPredicate =  dynamic_cast<CoincidentElements1D *>(theFunctor.get())) ||
      ( aPredicate =  dynamic_cast<CoincidentElements2D *>(theFunctor.get())) ||
      ( aPredicate =  dynamic_cast<CoincidentElements3D *>(theFunctor.get())) ||
      ( aPredicate =  dynamic_cast<OverConstrainedFace  *>(theFunctor.get())))
  {
    myExtractUnstructuredGrid->SetModeOfChanging(VTKViewer_ExtractUnstructuredGrid::eAdding);
    vtkUnstructuredGrid* aGrid = myVisualObj->GetUnstructuredGrid();
    vtkIdType aNbCells = aGrid->GetNumberOfCells();
    for( vtkIdType i = 0; i < aNbCells; i++ ){
      vtkIdType anObjId = myVisualObj->GetElemObjId(i);
      if(aPredicate->IsSatisfy(anObjId))
        myExtractUnstructuredGrid->RegisterCell(i);
    }
    if(!myExtractUnstructuredGrid->IsCellsRegistered())
      myExtractUnstructuredGrid->RegisterCell(-1);
    SetUnstructuredGrid(myVisualObj->GetUnstructuredGrid());
  }
  else if(FreeEdges* aFreeEdges = dynamic_cast<FreeEdges*>(theFunctor.get()))
  {
    SMESH::Controls::FreeEdges::TBorders aBorders;
    aFreeEdges->GetBoreders(aBorders);
    vtkUnstructuredGrid* aDataSet = vtkUnstructuredGrid::New();
    vtkUnstructuredGrid* aGrid = myVisualObj->GetUnstructuredGrid();
    aDataSet->SetPoints(aGrid->GetPoints());

    vtkIdType aNbCells = aBorders.size();
    vtkIdType aCellsSize = 3*aNbCells;
    vtkCellArray* aConnectivity = vtkCellArray::New();
    aConnectivity->Allocate( aCellsSize, 0 );
    
    vtkUnsignedCharArray* aCellTypesArray = vtkUnsignedCharArray::New();
    aCellTypesArray->SetNumberOfComponents( 1 );
    aCellTypesArray->Allocate( aNbCells * aCellTypesArray->GetNumberOfComponents() );
    
    vtkIdList *anIdList = vtkIdList::New();
    anIdList->SetNumberOfIds(2);
    
    FreeEdges::TBorders::const_iterator anIter = aBorders.begin();
    for(; anIter != aBorders.end(); anIter++){
      const FreeEdges::Border& aBorder = *anIter;
      int aNode[2] = {
        myVisualObj->GetNodeVTKId(aBorder.myPntId[0]),
        myVisualObj->GetNodeVTKId(aBorder.myPntId[1])
      };
      //cout<<"aNode = "<<aBorder.myPntId[0]<<"; "<<aBorder.myPntId[1]<<endl;
      if(aNode[0] >= 0 && aNode[1] >= 0){
        anIdList->SetId( 0, aNode[0] );
        anIdList->SetId( 1, aNode[1] );
        aConnectivity->InsertNextCell( anIdList );
        aCellTypesArray->InsertNextValue( VTK_LINE );
      }
    }
    
    VTKViewer_CellLocationsArray* aCellLocationsArray = VTKViewer_CellLocationsArray::New();
    aCellLocationsArray->SetNumberOfComponents( 1 );
    aCellLocationsArray->SetNumberOfTuples( aNbCells );
    
    aConnectivity->InitTraversal();
    for( vtkIdType idType = 0, *pts, npts; aConnectivity->GetNextCell( npts, pts ); idType++ )
      aCellLocationsArray->SetValue( idType, aConnectivity->GetTraversalLocation( npts ) );
    
    aDataSet->SetCells( aCellTypesArray, aCellLocationsArray,aConnectivity );

    SetUnstructuredGrid(aDataSet);
    aDataSet->Delete();
  }
  else if (( aPredicate = dynamic_cast<FreeNodes      *>(theFunctor.get())) ||
           ( aPredicate = dynamic_cast<CoincidentNodes*>(theFunctor.get())))
  {
    myExtractUnstructuredGrid->SetModeOfChanging(VTKViewer_ExtractUnstructuredGrid::eAdding);
    vtkIdType aNbNodes = myVisualObj->GetNbEntities(SMDSAbs_Node);
    for( vtkIdType i = 0; i < aNbNodes; i++ ){
      vtkIdType anObjId = myVisualObj->GetNodeObjId(i);
      if(aPredicate->IsSatisfy(anObjId))
        myExtractUnstructuredGrid->RegisterCell(i);
    }
    if(!myExtractUnstructuredGrid->IsCellsRegistered())
      myExtractUnstructuredGrid->RegisterCell(-1);
    SetUnstructuredGrid(myVisualObj->GetUnstructuredGrid());
  }
}




vtkMTimeType
SMESH_DeviceActor
::GetMTime()
{
  // cout << "DA " << this
  //      << " GF " << myGeomFilter;
  // if ( this->Property )
  //   cout << " P " << this->Property->GetMTime();
  // if ( this->BackfaceProperty != NULL )
  //   cout << " BP " << BackfaceProperty->GetMTime();
  // if ( this->Texture != NULL )
  //   cout << " T " << this->Texture->GetMTime();
  // cout << " U " << this->GetUserTransformMatrixMTime()
  //      << " M " << this->MTime.GetMTime() << endl;

  // cout << "DA " << this
  //      << " GF " << myGeomFilter
  //      << " " << this->Superclass::GetMTime()
  //      << " " << myExtractGeometry->GetMTime()
  //      << " " << myExtractUnstructuredGrid->GetMTime()
  //      << " " << myMergeFilter->GetMTime()
  //      << " " << myGeomFilter->GetMTime()
  //      << " " << myTransformFilter->GetMTime()
  //      << " " << myFaceOrientationFilter->GetMTime() << endl;

  vtkMTimeType mTime = this->Superclass::GetMTime();
  mTime = max(mTime,myExtractGeometry->GetMTime());
  mTime = max(mTime,myExtractUnstructuredGrid->GetMTime());
  mTime = max(mTime,myMergeFilter->GetMTime());
  mTime = max(mTime,myGeomFilter->GetMTime());
  mTime = max(mTime,myTransformFilter->GetMTime());
  mTime = max(mTime,myFaceOrientationFilter->GetMTime());
  return mTime;
}


void
SMESH_DeviceActor
::SetTransform(VTKViewer_Transform* theTransform)
{
  myTransformFilter->SetTransform(theTransform);
}


void
SMESH_DeviceActor
::SetShrink() 
{
  if ( !myIsShrinkable ) return;
  if ( vtkAlgorithmOutput* aDataSet = myPassFilter[ 0 ]->GetOutputPort() )
  {
    myShrinkFilter->SetInputConnection( aDataSet );
    myPassFilter[ 1 ]->SetInputConnection( myShrinkFilter->GetOutputPort() );
    myIsShrunk = true;
  }
}

void
SMESH_DeviceActor
::UnShrink() 
{
  if ( !myIsShrunk ) return;
  if ( vtkAlgorithmOutput* aDataSet = myPassFilter[ 0 ]->GetOutputPort() )
  {    
    myPassFilter[ 1 ]->SetInputConnection( aDataSet );
    myPassFilter[ 1 ]->Modified();
    myIsShrunk = false;
    Modified();
  }
}


void
SMESH_DeviceActor
::SetFacesOriented(bool theIsFacesOriented) 
{
  if ( vtkAlgorithmOutput* aDataSet = myTransformFilter->GetOutputPort() )
  {
    myIsFacesOriented = theIsFacesOriented;
    if( theIsFacesOriented )
      myFaceOrientationFilter->SetInputConnection( aDataSet );
    UpdateFaceOrientation();
  }
}

void
SMESH_DeviceActor
::SetFacesOrientationColor(double r,double g,double b)
{
  myFaceOrientation->GetProperty()->SetColor( r, g, b );
}

void
SMESH_DeviceActor
::GetFacesOrientationColor(double& r,double& g,double& b)
{
  myFaceOrientation->GetProperty()->GetColor( r, g, b );
}

void
SMESH_DeviceActor
::SetFacesOrientationScale(double theScale)
{
  myFaceOrientationFilter->SetOrientationScale( theScale );
}

double
SMESH_DeviceActor
::GetFacesOrientationScale()
{
  return myFaceOrientationFilter->GetOrientationScale();
}

void
SMESH_DeviceActor
::SetFacesOrientation3DVectors(bool theState)
{
  myFaceOrientationFilter->Set3dVectors( theState );
}

bool
SMESH_DeviceActor
::GetFacesOrientation3DVectors()
{
  return myFaceOrientationFilter->Get3dVectors();
}

void
SMESH_DeviceActor
::UpdateFaceOrientation()
{
  bool aShowFaceOrientation = myIsFacesOriented;
  aShowFaceOrientation &= vtkLODActor::GetVisibility(); //GetVisibility(); -- avoid calling GetUnstructuredGrid()  
  aShowFaceOrientation &= myRepresentation == eSurface;
  myFaceOrientation->SetVisibility(aShowFaceOrientation);
}


void
SMESH_DeviceActor
::SetRepresentation(EReperesent theMode)
{
  if ( myRepresentation == theMode )
    return;
  switch(theMode){
  case ePoint:
    myGeomFilter->SetInside(true);
    myGeomFilter->SetWireframeMode(false);
    GetProperty()->SetRepresentation(0);
    break;
  case eWireframe:
    myGeomFilter->SetInside(false);
    myGeomFilter->SetWireframeMode(true);
    GetProperty()->SetRepresentation(theMode);
    break;
  case eInsideframe:
    myGeomFilter->SetInside(true);
    myGeomFilter->SetWireframeMode(true);
    GetProperty()->SetRepresentation(1);
    break;
  case eSurface:
    myGeomFilter->SetInside(false);
    myGeomFilter->SetWireframeMode(false);
    GetProperty()->SetRepresentation(theMode);
  }
  SetMarkerEnabled(theMode == ePoint);
  myRepresentation = theMode;
  UpdateFaceOrientation();
  GetProperty()->Modified();
  myMapper->Modified();
  Modified();
}


void
SMESH_DeviceActor
::SetVisibility(int theMode)
{
  if(( theMode ) &&
     ( !myExtractUnstructuredGrid->GetInput() || 
       GetUnstructuredGrid()->GetNumberOfCells()))
  {
    vtkLODActor::SetVisibility(theMode);
  }else{
    vtkLODActor::SetVisibility(false);
  }
  UpdateFaceOrientation();
}


int
SMESH_DeviceActor
::GetVisibility()
{
  int visibi = vtkLODActor::GetVisibility();
  if(visibi && !GetUnstructuredGrid()->GetNumberOfCells()){
    vtkLODActor::SetVisibility(false);
    visibi = 0;
  }
  return visibi;
}


void
SMESH_DeviceActor
::AddToRender(vtkRenderer* theRenderer)
{
  theRenderer->AddActor(this);
  theRenderer->AddActor(myFaceOrientation);
}

void
SMESH_DeviceActor
::RemoveFromRender(vtkRenderer* theRenderer)
{
  theRenderer->RemoveActor(this);
  theRenderer->RemoveActor(myFaceOrientation);
}


int
SMESH_DeviceActor
::GetNodeObjId(int theVtkID)
{
  vtkIdType anID = theVtkID;

  if(IsImplicitFunctionUsed())
    anID = myExtractGeometry->GetNodeObjId(theVtkID);

  vtkIdType aRetID = myVisualObj->GetNodeObjId(anID);
  if(MYDEBUG) MESSAGE("GetNodeObjId - theVtkID = "<<theVtkID<<"; anID = "<<anID<<"; aRetID = "<<aRetID);
  return aRetID;
}

double* 
SMESH_DeviceActor
::GetNodeCoord(int theObjID)
{
  vtkDataSet* aDataSet = myMergeFilter->GetOutput();
  vtkIdType anID = myVisualObj->GetNodeVTKId(theObjID);
  double* aCoord = (anID >=0 && anID < aDataSet->GetNumberOfPoints()) ? aDataSet->GetPoint(anID) : NULL;
  if(MYDEBUG) MESSAGE("GetNodeCoord - theObjID = "<<theObjID<<"; anID = "<<anID);
  return aCoord;
}


int
SMESH_DeviceActor
::GetElemObjId(int theVtkID)
{
  vtkIdType anId = myGeomFilter->GetElemObjId(theVtkID);
  if(anId < 0) 
    return -1;

  vtkIdType anId2 = anId;
  if(IsImplicitFunctionUsed())
    anId2 = myExtractGeometry->GetElemObjId(anId);
  if(anId2 < 0) 
    return -1;

  vtkIdType anId3 = myExtractUnstructuredGrid->GetInputId(anId2);
  if(anId3 < 0) 
    return -1;

  vtkIdType aRetID = myVisualObj->GetElemObjId(anId3);
  if(MYDEBUG) 
     MESSAGE("GetElemObjId - theVtkID = "<<theVtkID<<"; anId2 = "<<anId2<<"; anId3 = "<<anId3<<"; aRetID = "<<aRetID);
  return aRetID;
}

vtkCell* 
SMESH_DeviceActor
::GetElemCell(int theObjID)
{
  vtkDataSet* aDataSet = myVisualObj->GetUnstructuredGrid();
  vtkIdType aGridID = myVisualObj->GetElemVTKId(theObjID);
  vtkCell* aCell = (aGridID >= 0 ) ? aDataSet->GetCell(aGridID) : NULL;
  if(MYDEBUG) 
    MESSAGE("GetElemCell - theObjID = "<<theObjID<<"; aGridID = "<<aGridID);
  return aCell;
}


double 
SMESH_DeviceActor
::GetShrinkFactor()
{
  return myShrinkFilter->GetShrinkFactor();
}

void
SMESH_DeviceActor
::SetShrinkFactor(double theValue)
{
  theValue = theValue > 0.1? theValue: 0.8;
  myShrinkFilter->SetShrinkFactor(theValue);
  Modified();
}


void
SMESH_DeviceActor
::SetHighlited(bool theIsHighlited)
{
  if ( myIsHighlited == theIsHighlited )
    return;
  myIsHighlited = theIsHighlited;
  Modified();
}

void
SMESH_DeviceActor
::Render(vtkRenderer *ren, vtkMapper* m)
{
  int aResolveCoincidentTopology = vtkMapper::GetResolveCoincidentTopology();
  double aStoredFactor, aStoredUnit; 
  vtkMapper::GetResolveCoincidentTopologyPolygonOffsetParameters(aStoredFactor,aStoredUnit);

  vtkMapper::SetResolveCoincidentTopologyToPolygonOffset();
  double aFactor = myPolygonOffsetFactor, aUnits = myPolygonOffsetUnits;
  if(myIsHighlited){
    static double EPS = .01;
    aUnits *= (1.0-EPS);
  }
  vtkMapper::SetResolveCoincidentTopologyPolygonOffsetParameters(aFactor,aUnits);
  vtkLODActor::Render(ren,m);

  vtkMapper::SetResolveCoincidentTopologyPolygonOffsetParameters(aStoredFactor,aStoredUnit);
  vtkMapper::SetResolveCoincidentTopology(aResolveCoincidentTopology);
}


void
SMESH_DeviceActor
::SetPolygonOffsetParameters(double factor, 
                             double units)
{
  myPolygonOffsetFactor = factor;
  myPolygonOffsetUnits = units;
}

/*!
 * On/Off representation 2D quadratic element as arked polygon
 */
void SMESH_DeviceActor::SetQuadraticArcMode(bool theFlag){
  myGeomFilter->SetQuadraticArcMode(theFlag);
}

/*!
 * Return true if 2D quadratic element displayed as arked polygon
 */
bool SMESH_DeviceActor::GetQuadraticArcMode(){
  return myGeomFilter->GetQuadraticArcMode();
}
/*!
 * Set Max angle for representation 2D quadratic element as arked polygon
 */
void SMESH_DeviceActor::SetQuadraticArcAngle(double theMaxAngle){
  myGeomFilter->SetQuadraticArcAngle(theMaxAngle);
}

/*!
 * Return Max angle of the representation 2D quadratic element as arked polygon
 */
double SMESH_DeviceActor::GetQuadraticArcAngle(){
  return myGeomFilter->GetQuadraticArcAngle();
}

/*!
 * Set point marker enabled
 * \param theMarkerEnabled flag to enable/disable point marker
 */
void SMESH_DeviceActor::SetMarkerEnabled( bool theMarkerEnabled )
{
  myMapper->SetMarkerEnabled( theMarkerEnabled );
}

/*!
 * Set point marker enabled
 * \param theBallEnabled flag to enable/disable ball drawing
 */
void SMESH_DeviceActor::SetBallEnabled( bool theBallEnabled ) {
  myMapper->SetBallEnabled( theBallEnabled );
}

/*!
 * Set point marker scale factor
 * \param theBallScale double value which specifies a scale factor of ball element
 */
void SMESH_DeviceActor::SetBallScale( double theBallScale )
{
  myMapper->SetBallScale( theBallScale );
  myMapper->Modified();
}

/*!
 * Set standard point marker
 * \param theMarkerType type of the marker
 */
void SMESH_DeviceActor::SetMarkerStd( VTK::MarkerType theMarkerType, VTK::MarkerScale theMarkerScale )
{
  myMapper->SetMarkerStd( theMarkerType, theMarkerScale );
}

/*!
 * Set custom point marker
 * \param theMarkerId id of the marker texture
 * \param theMarkerTexture marker texture
 */
void SMESH_DeviceActor::SetMarkerTexture( int theMarkerId, VTK::MarkerTexture theMarkerTexture )
{
  myMapper->SetMarkerTexture( theMarkerId, theMarkerTexture );
}

/*!
 * Get type of the point marker
 * \return type of the point marker
 */
VTK::MarkerType SMESH_DeviceActor::GetMarkerType()
{
  return myMapper->GetMarkerType();
}

/*!
  Get scale of the point marker
  \return scale of the point marker
*/
VTK::MarkerScale SMESH_DeviceActor::GetMarkerScale()
{
  return myMapper->GetMarkerScale();
}

/*!
 * Get texture identifier of the point marker
 * \return texture identifier of the point marker
 */
int SMESH_DeviceActor::GetMarkerTexture()
{
  return myMapper->GetMarkerTexture();
}

/*!
 * Get scale factor of ball element
 * \return scale factor of ball element
 */
double SMESH_DeviceActor::GetBallScale()
{
  return myMapper->GetBallScale();
}

void SMESH_DeviceActor::SetCoincident3DAllowed(bool theFlag) {
  myGeomFilter->SetAppendCoincident3D(theFlag);
}

bool SMESH_DeviceActor::IsCoincident3DAllowed() const {
  return myGeomFilter->GetAppendCoincident3D();
}
