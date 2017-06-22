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

#ifndef SALOME_ExtractGeometry_HeaderFile
#define SALOME_ExtractGeometry_HeaderFile

#include "SMESH_Object.h"

#include <vtkExtractGeometry.h>
#include <vector>

#include "VTKViewer.h"

class SMESHOBJECT_EXPORT SMESH_ExtractGeometry : public vtkExtractGeometry{
public:
  vtkTypeMacro(SMESH_ExtractGeometry,vtkExtractGeometry);

  static SMESH_ExtractGeometry *New();

  void SetStoreMapping(bool theStoreMapping){
    myStoreMapping = theStoreMapping;
    Modified();
  }
  bool GetStoreMapping(){ return myStoreMapping;}

  virtual vtkIdType GetNodeObjId(int theVtkID);
  virtual vtkIdType GetElemObjId(int theVtkID);

protected:
  SMESH_ExtractGeometry();
  ~SMESH_ExtractGeometry();

  // Usual data generation method
  virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
  bool myStoreMapping;   
  typedef std::vector<vtkIdType> TVectorId;
  TVectorId myElemVTK2ObjIds;
  TVectorId myNodeVTK2ObjIds;

  SMESH_ExtractGeometry(const SMESH_ExtractGeometry&);  // Not implemented.
  void operator=(const SMESH_ExtractGeometry&);  // Not implemented.
};


#endif


