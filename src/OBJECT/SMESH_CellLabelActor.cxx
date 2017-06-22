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
//  File   : SMESH_CellLabelActor.cxx
//  Author : Roman NIKOLAEV
//  Module : SMESH
//
#include "SMESH_CellLabelActor.h"

#include "SMESH_ExtractGeometry.h"

#include <VTKViewer_TransformFilter.h>
#include <VTKViewer_CellCenters.h>
#include <VTKViewer_ExtractUnstructuredGrid.h>

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
#include <vtkCellData.h>

vtkStandardNewMacro(SMESH_CellLabelActor);

/*!
  Constructor.
*/
SMESH_CellLabelActor::SMESH_CellLabelActor()
{
  //Definition of cells numbering pipeline
  //---------------------------------------
  myCellsNumDataSet = vtkUnstructuredGrid::New();

  myCellCenters = VTKViewer_CellCenters::New();
  myCellCenters->SetInputData(myCellsNumDataSet);

  myClsMaskPoints = vtkMaskPoints::New();
  myClsMaskPoints->SetInputConnection(myCellCenters->GetOutputPort());
  myClsMaskPoints->SetOnRatio(1);

  myClsSelectVisiblePoints = vtkSelectVisiblePoints::New();
  myClsSelectVisiblePoints->SetInputConnection(myClsMaskPoints->GetOutputPort());
  myClsSelectVisiblePoints->SelectInvisibleOff();
  myClsSelectVisiblePoints->SetTolerance(0.1);

  myClsLabeledDataMapper = vtkLabeledDataMapper::New();
  myClsLabeledDataMapper->SetInputConnection(myClsSelectVisiblePoints->GetOutputPort());

  myClsLabeledDataMapper->SetLabelFormat("%d");
  myClsLabeledDataMapper->SetLabelModeToLabelScalars();

  myClsTextProp = vtkTextProperty::New();
  myClsTextProp->SetFontFamilyToTimes();
  myClsTextProp->SetFontSize(12);
  myClsTextProp->SetBold(1);
  myClsTextProp->SetItalic(0);
  myClsTextProp->SetShadow(0);
  myClsTextProp->SetColor( 0, 1, 0 );
  myClsLabeledDataMapper->SetLabelTextProperty(myClsTextProp);
    
  myIsCellsLabeled = false;

  myCellsLabels = vtkActor2D::New();
  myCellsLabels->SetMapper(myClsLabeledDataMapper);
  myCellsLabels->SetVisibility(myIsCellsLabeled);

  vtkCallbackCommand* callBackCommand = vtkCallbackCommand::New();
  callBackCommand->SetClientData(this);
  callBackCommand->SetCallback(SMESH_CellLabelActor::ProcessEvents);

  myTransformFilter->AddObserver("VTKViewer_TransformFilter::TransformationFinished",
                                 callBackCommand);
  callBackCommand->Delete();
}


/*!
  Destructor.
*/
SMESH_CellLabelActor::~SMESH_CellLabelActor()
{
  //Deleting of cells numbering pipeline
  //---------------------------------------
  myCellsNumDataSet->Delete();
  myCellsLabels->Delete();
  // commented: porting to vtk 5.0
  //  myClsMaskPoints->UnRegisterAllOutputs();
  myClsMaskPoints->Delete();
  // commented: porting to vtk 5.0
  //  myCellCenters->UnRegisterAllOutputs();
  myCellCenters->Delete();

  myClsLabeledDataMapper->RemoveAllInputs();
  myClsLabeledDataMapper->Delete();
  // commented: porting to vtk 5.0
  //  myClsSelectVisiblePoints->UnRegisterAllOutputs();
  myClsSelectVisiblePoints->Delete();
  myClsTextProp->Delete();
}


void SMESH_CellLabelActor::SetFontProperties( SMESH::LabelFont family, int size,
                                              bool bold, bool italic, bool shadow,
                                              double r, double g, double b  )
{
  switch ( family ) {
  case SMESH::FntArial:
    myClsTextProp->SetFontFamilyToArial(); break;
  case SMESH::FntCourier:
    myClsTextProp->SetFontFamilyToCourier(); break;
  case SMESH::FntTimes:
  default:
    myClsTextProp->SetFontFamilyToTimes(); break;
  }    
  myClsTextProp->SetFontSize( size );
  myClsTextProp->SetBold( bold );
  myClsTextProp->SetItalic( italic );
  myClsTextProp->SetShadow( shadow );
  myClsTextProp->SetColor( r, g, b );
}

void SMESH_CellLabelActor::SetCellsLabeled(bool theIsCellsLabeled)
{
  myIsCellsLabeled = theIsCellsLabeled;

  myCellsLabels->SetVisibility(false);

  myTransformFilter->Update();
  vtkUnstructuredGrid* aGrid = vtkUnstructuredGrid::SafeDownCast(myTransformFilter->GetOutput());

  if ( myIsCellsLabeled && aGrid )
  {
    myCellsNumDataSet->ShallowCopy(aGrid);
    vtkUnstructuredGrid *aDataSet = myCellsNumDataSet;
    int aNbElem = aDataSet->GetNumberOfCells();
    vtkIntArray *anArray = vtkIntArray::New();
    anArray->SetNumberOfValues(aNbElem);
    myExtractUnstructuredGrid->BuildOut2InMap();
    for(int anId = 0; anId < aNbElem; anId++)
    {
      vtkIdType id = anId;
      if(IsImplicitFunctionUsed())
        id = myExtractGeometry->GetElemObjId(id);
      id = myExtractUnstructuredGrid->GetInputId(id);
      id = (id >=0) ? id : anId;
      int aSMDSId = myVisualObj->GetElemObjId(id);
      anArray->SetValue(anId,aSMDSId);
    }
    aDataSet->GetCellData()->SetScalars(anArray);
    myCellCenters->SetInputData(aDataSet);
    myCellsLabels->SetVisibility(GetVisibility());
  }
}

void SMESH_CellLabelActor::SetVisibility(int theMode)
{
  SMESH_DeviceActor::SetVisibility(theMode);
  myCellsLabels->VisibilityOff();
  if(myIsCellsLabeled && theMode)
    myCellsLabels->VisibilityOn();
}

void SMESH_CellLabelActor::AddToRender(vtkRenderer* theRenderer)
{
  SMESH_DeviceActor::AddToRender(theRenderer);
  myClsSelectVisiblePoints->SetRenderer(theRenderer);
  theRenderer->AddActor2D(myCellsLabels);
}

void SMESH_CellLabelActor::RemoveFromRender(vtkRenderer* theRenderer)
{
  theRenderer->RemoveActor(myCellsLabels);
  SMESH_DeviceActor::RemoveFromRender(theRenderer);
}

void SMESH_CellLabelActor::UpdateLabels()
{
  if(myIsCellsLabeled)
    SetCellsLabeled(myIsCellsLabeled);
}


void SMESH_CellLabelActor::ProcessEvents(vtkObject* vtkNotUsed(theObject),
                                         unsigned long theEvent,
                                         void* theClientData,
                                         void* vtkNotUsed(theCallData))
{
  SMESH_CellLabelActor* self = reinterpret_cast<SMESH_CellLabelActor*>(theClientData);
  if(self)
    self->UpdateLabels();
}
