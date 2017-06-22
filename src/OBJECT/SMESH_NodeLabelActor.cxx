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
//  File   : SMESH_NodeLabelActor.cxx
//  Author : Roman NIKOLAEV
//  Module : SMESH
//
#include "SMESH_NodeLabelActor.h"

#include <VTKViewer_TransformFilter.h>

#include <vtkObjectFactory.h>
#include <vtkCallbackCommand.h>
#include <vtkMaskPoints.h>
#include <vtkSelectVisiblePoints.h>
#include <vtkLabeledDataMapper.h>
#include <vtkActor2D.h>
#include <vtkTextProperty.h>
#include <vtkPointData.h>
#include <vtkProperty2D.h>
#include <vtkRenderer.h>
#include <vtkUnstructuredGrid.h>

vtkStandardNewMacro(SMESH_NodeLabelActor);

/*!
  Constructor.
*/
SMESH_NodeLabelActor::SMESH_NodeLabelActor()
{
  //Definition of points numbering pipeline
  //---------------------------------------
  myPointsNumDataSet = vtkUnstructuredGrid::New();

  myPtsMaskPoints = vtkMaskPoints::New();
  myPtsMaskPoints->SetInputData(myPointsNumDataSet);
  myPtsMaskPoints->SetOnRatio(1);

  myPtsSelectVisiblePoints = vtkSelectVisiblePoints::New();
  myPtsSelectVisiblePoints->SetInputConnection(myPtsMaskPoints->GetOutputPort());
  myPtsSelectVisiblePoints->SelectInvisibleOff();
  myPtsSelectVisiblePoints->SetTolerance(0.1);

  myPtsLabeledDataMapper = vtkLabeledDataMapper::New();
  myPtsLabeledDataMapper->SetInputConnection(myPtsSelectVisiblePoints->GetOutputPort());
  myPtsLabeledDataMapper->SetLabelFormat("%d");
  myPtsLabeledDataMapper->SetLabelModeToLabelScalars();

  myPtsTextProp = vtkTextProperty::New();
  myPtsTextProp->SetFontFamilyToTimes();
  myPtsTextProp->SetFontSize(10);
  myPtsTextProp->SetBold(1);
  myPtsTextProp->SetItalic(0);
  myPtsTextProp->SetShadow(0);
  myPtsTextProp->SetColor( 1, 1, 1 );
  myPtsLabeledDataMapper->SetLabelTextProperty(myPtsTextProp);

  myIsPointsLabeled = false;

  myPointLabels = vtkActor2D::New();
  myPointLabels->SetMapper(myPtsLabeledDataMapper);
  myPointLabels->SetVisibility(myIsPointsLabeled);

  vtkCallbackCommand* callBackCommand = vtkCallbackCommand::New();
  callBackCommand->SetClientData(this);
  callBackCommand->SetCallback(SMESH_NodeLabelActor::ProcessEvents);

  myTransformFilter->AddObserver("VTKViewer_TransformFilter::TransformationFinished",
                                 callBackCommand);
  callBackCommand->Delete();
}

/*!
  Destructor
*/
SMESH_NodeLabelActor::~SMESH_NodeLabelActor()
{
  //Deleting of points numbering pipeline
  //---------------------------------------
  myPointsNumDataSet->Delete();

  // commented: porting to vtk 5.0
  //  myPtsLabeledDataMapper->RemoveAllInputs();
  myPtsLabeledDataMapper->Delete();

  // commented: porting to vtk 5.0
  //  myPtsSelectVisiblePoints->UnRegisterAllOutputs();
  myPtsSelectVisiblePoints->Delete();

  // commented: porting to vtk 5.0
  //  myPtsMaskPoints->UnRegisterAllOutputs();
  myPtsMaskPoints->Delete();
  myPointLabels->Delete();
  myPtsTextProp->Delete();
}

void SMESH_NodeLabelActor::SetFontProperties( SMESH::LabelFont family, int size, 
                                              bool bold, bool italic, bool shadow,
                                              double r, double g, double b )
{
  switch ( family ) {
  case SMESH::FntArial:
    myPtsTextProp->SetFontFamilyToArial(); break;
  case SMESH::FntCourier:
    myPtsTextProp->SetFontFamilyToCourier(); break;
  case SMESH::FntTimes:
  default:
    myPtsTextProp->SetFontFamilyToTimes(); break;
  }
  myPtsTextProp->SetFontSize( size );
  myPtsTextProp->SetBold( bold );
  myPtsTextProp->SetItalic( italic );
  myPtsTextProp->SetShadow( shadow ); 
  myPtsTextProp->SetColor( r, g, b ); 
}

void SMESH_NodeLabelActor::SetPointsLabeled(bool theIsPointsLabeled)
{
  myIsPointsLabeled = theIsPointsLabeled;

  myPointLabels->SetVisibility( false );

  myTransformFilter->Update();
  vtkDataSet* aGrid = vtkUnstructuredGrid::SafeDownCast(myTransformFilter->GetOutput());

  if ( myIsPointsLabeled && aGrid )
  {
    myPointsNumDataSet->ShallowCopy(aGrid);
    vtkUnstructuredGrid *aDataSet = myPointsNumDataSet;
    
    int aNbElem = aDataSet->GetNumberOfPoints();
    
    vtkIntArray *anArray = vtkIntArray::New();
    anArray->SetNumberOfValues( aNbElem );
    
    for ( vtkIdType anId = 0; anId < aNbElem; anId++ )
    {
      int aSMDSId = myVisualObj->GetNodeObjId( anId );
      anArray->SetValue( anId, aSMDSId );
    }
    
    aDataSet->GetPointData()->SetScalars( anArray );
    myPtsMaskPoints->SetInputData( aDataSet );
    myPointLabels->SetVisibility( GetVisibility() );
    anArray->Delete();
  }
}


void SMESH_NodeLabelActor::SetVisibility(int theMode)
{
  SMESH_DeviceActor::SetVisibility(theMode);
  myPointLabels->VisibilityOff();
  if(myIsPointsLabeled && theMode)
    myPointLabels->VisibilityOn();
}


void SMESH_NodeLabelActor::AddToRender(vtkRenderer* theRenderer)
{
  SMESH_DeviceActor::AddToRender(theRenderer);
  myPtsSelectVisiblePoints->SetRenderer(theRenderer);
  theRenderer->AddActor2D(myPointLabels);
}

void SMESH_NodeLabelActor::RemoveFromRender(vtkRenderer* theRenderer)
{
  theRenderer->RemoveActor(myPointLabels);
  SMESH_DeviceActor::RemoveFromRender(theRenderer);
}

void SMESH_NodeLabelActor::UpdateLabels()
{
  if(myIsPointsLabeled)
    SetPointsLabeled(myIsPointsLabeled);
}


void SMESH_NodeLabelActor::ProcessEvents(vtkObject* vtkNotUsed(theObject),
                                         unsigned long theEvent,
                                         void* theClientData,
                                         void* vtkNotUsed(theCallData))
{
  SMESH_NodeLabelActor* self = reinterpret_cast<SMESH_NodeLabelActor*>(theClientData);
  if(self)
    self->UpdateLabels();
}
