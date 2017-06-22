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

#ifndef SMESH_FACEORIENTATIONFILTER_H
#define SMESH_FACEORIENTATIONFILTER_H

#include "SMESH_Object.h"

#include <vtkPolyDataAlgorithm.h>

class vtkGlyph3D;
class vtkGlyphSource2D;
class vtkMaskPoints;

class VTKViewer_CellCenters;

class SMESHOBJECT_EXPORT SMESH_FaceOrientationFilter : public vtkPolyDataAlgorithm
{
public:
  vtkTypeMacro( SMESH_FaceOrientationFilter, vtkPolyDataAlgorithm );

  /*!Create a new SMESH_FaceOrientationFilter.*/
  static SMESH_FaceOrientationFilter *New();

  void SetOrientationScale( double );
  double GetOrientationScale() const { return myOrientationScale; }

  void Set3dVectors( bool );
  bool Get3dVectors() const { return my3dVectors; }

protected:
  SMESH_FaceOrientationFilter();
  virtual ~SMESH_FaceOrientationFilter();

  virtual int RequestData(vtkInformation *, vtkInformationVector **,
                          vtkInformationVector *); //generate output data

  virtual int FillInputPortInformation(int port, vtkInformation *info);

  vtkPolyData* CreateArrowPolyData();

private:
  SMESH_FaceOrientationFilter( const SMESH_FaceOrientationFilter& );  //!< Not implemented.
  void operator=( const SMESH_FaceOrientationFilter& );               //!< Not implemented.

private:
  bool my3dVectors;
  double myOrientationScale;
  vtkPolyData* myArrowPolyData;
  vtkPolyData* myFacePolyData;
  VTKViewer_CellCenters* myFaceCenters;
  vtkMaskPoints* myFaceMaskPoints;
  vtkGlyphSource2D* myGlyphSource;
  vtkGlyph3D* myBaseGlyph;
};

#endif
