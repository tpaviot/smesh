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
//  File   : SMESH_SVTKActor.h
//  Author : Roman NIKOLAEV
//  Module : SMESH
//
#ifndef SMESH_SVTKACTOR_H
#define SMESH_SVTKACTOR_H

#include "SMESH_Object.h"
#include <SVTK_Actor.h>


class SVTK_Actor;
class vtkUnstructureGrid;
class vtkDataSetMapper;

class SMESHOBJECT_EXPORT SMESH_SVTKActor : public SVTK_Actor {

public:
  static SMESH_SVTKActor* New();

  vtkTypeMacro(SMESH_SVTKActor, SVTK_Actor);

  void SetBallScale(double theSize);
  void SetBallSize(float theSize);
  void Set0DSize(float theSize);

  //! To publish the actor an all its internal devices
  virtual
  void
  AddToRender(vtkRenderer* theRendere); 

  virtual void SetVisibility( int theVisibility );

  //! Initialiaze the instance completely
  virtual void
  Initialize();

  //! Allow to recostruct selected cells from source SALOME_Actor and map of subindexes
  virtual void
  MapCells(SALOME_Actor* theMapActor,
           const TColStd_IndexedMapOfInteger& theMapIndex);


  //! To remove the actor an all its internal devices
  virtual
  void
  RemoveFromRender(vtkRenderer* theRendere);

  void SetVisualObject(TVisualObjPtr theVisualObj);

 protected:
  SVTK_DeviceActor* my0DActor;
  SVTK_DeviceActor* myBallActor;

  vtkUnstructuredGrid* my0DGrid;
  vtkUnstructuredGrid* myBallGrid;
  TVisualObjPtr myVisualObj;
  
  SMESH_SVTKActor();
  virtual ~SMESH_SVTKActor();
};

#endif
