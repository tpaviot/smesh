// Copyright (C) 2007-2016  CEA/DEN, EDF R&D, OPEN CASCADE
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

#include "SMESH_FaceOrientationFilter.h"
#include "SMESH_ActorUtils.h"

#include "SUIT_Session.h"
#include "SUIT_ResourceMgr.h"

#include <VTKViewer_CellCenters.h>

#include <vtkCellData.h>
#include <vtkDataSet.h>
#include <vtkPolyData.h>
#include <vtkObjectFactory.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>

#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkMaskPoints.h>
#include <vtkGlyph3D.h>
#include <vtkGlyphSource2D.h>

#include <QColor>

#define PI   3.14159265359

vtkStandardNewMacro(SMESH_FaceOrientationFilter);

/*!
 * \class SMESH_FaceOrientationFilter
 * Passive filter take a polydata as input and create a dataset as output.
 */

SMESH_FaceOrientationFilter::SMESH_FaceOrientationFilter()
{
  SUIT_ResourceMgr* mgr = SUIT_Session::session()->resourceMgr();
  myOrientationScale = mgr->doubleValue( "SMESH", "orientation_scale", 0.1 );
  my3dVectors = mgr->booleanValue( "SMESH", "orientation_3d_vectors", false );

  myArrowPolyData = CreateArrowPolyData();

  myFacePolyData = vtkPolyData::New();

  myFaceCenters = VTKViewer_CellCenters::New();
  myFaceCenters->SetInputData(myFacePolyData);

  myFaceMaskPoints = vtkMaskPoints::New();
  myFaceMaskPoints->SetInputConnection(myFaceCenters->GetOutputPort());
  myFaceMaskPoints->SetOnRatio(1);

  myGlyphSource = vtkGlyphSource2D::New();
  myGlyphSource->SetGlyphTypeToThickArrow();
  myGlyphSource->SetFilled(0);
  myGlyphSource->SetCenter(0.5, 0.0, 0.0);

  myBaseGlyph = vtkGlyph3D::New();
  myBaseGlyph->SetInputConnection(myFaceMaskPoints->GetOutputPort());
  myBaseGlyph->SetVectorModeToUseVector();
  myBaseGlyph->SetScaleModeToDataScalingOff();
  myBaseGlyph->SetColorModeToColorByScalar();
  if( my3dVectors )
    myBaseGlyph->SetSourceData(myArrowPolyData);
  else
    myBaseGlyph->SetSourceConnection(myGlyphSource->GetOutputPort());
}

SMESH_FaceOrientationFilter::~SMESH_FaceOrientationFilter()
{
  myArrowPolyData->Delete();
  myFacePolyData->Delete();
  myFaceCenters->Delete();
  myFaceMaskPoints->Delete();
  myGlyphSource->Delete();
  myBaseGlyph->Delete();
}

void SMESH_FaceOrientationFilter::SetOrientationScale( double theScale )
{
  myOrientationScale = theScale;
  Modified();
}

void SMESH_FaceOrientationFilter::Set3dVectors( bool theState )
{
  my3dVectors = theState;
  if( my3dVectors )
    myBaseGlyph->SetSourceData(myArrowPolyData);
  else
    myBaseGlyph->SetSourceConnection(myGlyphSource->GetOutputPort());
  Modified();
}

vtkPolyData* SMESH_FaceOrientationFilter::CreateArrowPolyData()
{
  vtkPoints* points = vtkPoints::New();
  vtkCellArray* polys = vtkCellArray::New();

  float l1 = 0.8;
  float l2 = 1.0;
  int n = 16;
  float r1 = 0.04;
  float r2 = 0.08;
  float angle = 2. * PI / n;
  float p[3];
  vtkIdType c3[3];
  vtkIdType c4[4];

  float p0[3] = { 0.0, 0.0, 0.0 };
  float p1[3] = {  l1, 0.0, 0.0 };
  float p2[3] = {  l2, 0.0, 0.0 };

  points->InsertPoint( 0, p0 );
  points->InsertPoint( 1, p1 );
  points->InsertPoint( 2, p2 );

  // shaft
  for( int i = 0; i < n; i++ )
  {
    p[0] = 0;
    p[1] = r1 * sin( i * angle );
    p[2] = r1 * cos( i * angle );
    points->InsertPoint( i + 3, p );

    p[0] = l1;
    points->InsertPoint( i + 3 + n, p );
  }

  // insert the last cells outside a loop
  {
    c3[0] = 0;
    c3[1] = 3;
    c3[2] = 3 + n - 1;
    polys->InsertNextCell( 3, c3 );

    c4[0] = 3;
    c4[1] = 3 + n - 1;
    c4[2] = 3 + 2 * n - 1;
    c4[3] = 3 + n;
    polys->InsertNextCell( 4, c4 );
  }
  for( int i = 0; i < n - 1; i++ )
  {
    c3[0] = 0;
    c3[1] = i + 3;
    c3[2] = i + 4;
    polys->InsertNextCell( 3, c3 );

    c4[0] = i + 3;
    c4[1] = i + 4;
    c4[2] = i + 4 + n;
    c4[3] = i + 3 + n;
    polys->InsertNextCell( 4, c4 );
  }

  // cone
  for( int i = 0; i < n; i++ )
  {
    p[0] = l1;
    p[1] = r2 * sin( i * angle );
    p[2] = r2 * cos( i * angle );
    points->InsertPoint( i + 3 + 2 * n, p );
  }

  // insert the last cells outside a loop
  {
    c3[0] = 1;
    c3[1] = 3 + 2 * n;
    c3[2] = 3 + 2 * n + n - 1;
    polys->InsertNextCell( 3, c3 );

    c3[0] = 2;
    polys->InsertNextCell( 3, c3 );
  }
  for( int i = 0; i < n - 1; i++ )
  {
    c3[0] = 1;
    c3[1] = 3 + i + 2 * n;
    c3[2] = 3 + i + 2 * n + 1;
    polys->InsertNextCell( 3, c3 );

    c3[0] = 2;
    polys->InsertNextCell( 3, c3 );
  }

  vtkPolyData* aPolyData = vtkPolyData::New();

  aPolyData->SetPoints(points);
  points->Delete();

  aPolyData->SetPolys(polys);
  polys->Delete();

  return aPolyData;
}

void GetFaceParams( vtkCell* theFace, double theNormal[3], double& theSize ) 
{
  vtkPoints* aPoints = theFace->GetPoints();

  // here we get first 3 points from the face and calculate the normal as a cross-product of vectors
  double x0 = aPoints->GetPoint(0)[0], y0 = aPoints->GetPoint(0)[1], z0 = aPoints->GetPoint(0)[2];
  double x1 = aPoints->GetPoint(1)[0], y1 = aPoints->GetPoint(1)[1], z1 = aPoints->GetPoint(1)[2];
  double x2 = aPoints->GetPoint(2)[0], y2 = aPoints->GetPoint(2)[1], z2 = aPoints->GetPoint(2)[2];

  theNormal[0] = ( y1 - y0 ) * ( z2 - z0 ) - ( z1 - z0 ) * ( y2 - y0 );
  theNormal[1] = ( z1 - z0 ) * ( x2 - x0 ) - ( x1 - x0 ) * ( z2 - z0 );
  theNormal[2] = ( x1 - x0 ) * ( y2 - y0 ) - ( y1 - y0 ) * ( x2 - x0 );

  double* aBounds = theFace->GetBounds();
  theSize = pow( pow( aBounds[1] - aBounds[0], 2 ) +
                 pow( aBounds[3] - aBounds[2], 2 ) +
                 pow( aBounds[5] - aBounds[4], 2 ), 0.5 );
}

/*!
 * Execute method. Output calculation.
 */
int SMESH_FaceOrientationFilter::RequestData(
  vtkInformation *request,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the input and ouptut
  vtkDataSet *input = vtkDataSet::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData *output = vtkPolyData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  myFacePolyData->Initialize();
  myFacePolyData->ShallowCopy(input);

  vtkCellArray* aFaces = vtkCellArray::New();

  vtkFloatArray* aVectors = vtkFloatArray::New();
  aVectors->SetNumberOfComponents(3);

  int anAllFaces = 0;
  double anAverageSize = 0;

  vtkIdList* aNeighborIds = vtkIdList::New();

  for(int aCellId = 0, aNbCells = input->GetNumberOfCells(); aCellId < aNbCells; aCellId++)
  {
    vtkCell* aCell = input->GetCell(aCellId);

    if( aCell->GetNumberOfFaces() == 0 && aCell->GetNumberOfPoints() > 2 ) // cell is a face
    {
      double aSize, aNormal[3];
      GetFaceParams( aCell, aNormal, aSize );

      aFaces->InsertNextCell(aCell);
      aVectors->InsertNextTuple(aNormal);

      anAllFaces++;
      anAverageSize += aSize;

      continue;
    }

    for(int aFaceId = 0, aNbFaces = aCell->GetNumberOfFaces(); aFaceId < aNbFaces; aFaceId++)
    {
      vtkCell* aFace = aCell->GetFace(aFaceId);

      input->GetCellNeighbors( aCellId, aFace->PointIds, aNeighborIds );
      if( aNeighborIds->GetNumberOfIds() > 0 )
        continue;

      double aSize, aNormal[3];
      GetFaceParams( aFace, aNormal, aSize );

      aFaces->InsertNextCell(aFace->GetPointIds());
      aVectors->InsertNextTuple(aNormal);

      anAllFaces++;
      anAverageSize += aSize;
    }
  }
  aNeighborIds->Delete();

  myFacePolyData->SetPolys(aFaces);
  aFaces->Delete();

  myFacePolyData->GetCellData()->SetScalars(0);
  myFacePolyData->GetCellData()->SetVectors(aVectors);
  aVectors->Delete();

  if( anAllFaces == 0 )
    return 0;

  anAverageSize /= anAllFaces;
  anAverageSize *= myOrientationScale;

  myBaseGlyph->SetScaleFactor( anAverageSize );
  myBaseGlyph->Update();

  output->ShallowCopy( myBaseGlyph->GetOutput() );

  return 1;
}

int SMESH_FaceOrientationFilter::FillInputPortInformation(int, vtkInformation *info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
  return 1;
}
