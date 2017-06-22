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
//  File   : SMESH_CellLabelActor.h
//  Author : Roman NIKOLAEV
//  Module : SMESH
//
#ifndef SMESH_CELL_LABEL_ACTOR_H
#define SMESH_CELL_LABEL_ACTOR_H

#include "SMESH_DeviceActor.h"
#include "SMESH_ActorUtils.h"

class vtkSelectVisiblePoints;
class vtkLabeledDataMapper;
class vtkActor2D;
class vtkMaskPoints;
class vtkUnstructuredGrid;
class vtkTextProperty;

class VTKViewer_CellCenters;


class SMESHOBJECT_EXPORT SMESH_CellLabelActor : public SMESH_DeviceActor {
public:
  static SMESH_CellLabelActor* New();

  static void ProcessEvents(vtkObject* theObject,
                            unsigned long theEvent,
                            void* theClientData,
                            void* theCallData);


  vtkTypeMacro(SMESH_CellLabelActor, SMESH_DeviceActor);


  virtual void SetCellsLabeled(bool theIsCellsLabeled);
  virtual bool GetCellsLabeled(){ return myIsCellsLabeled;}

  virtual void SetVisibility(int theMode);

  virtual void AddToRender(vtkRenderer* theRenderer);
  virtual void RemoveFromRender(vtkRenderer* theRenderer);
  
  virtual void SetFontProperties( SMESH::LabelFont family, int size,
                                  bool bold, bool italic, bool shadow,
                                  double r, double g, double b );

  void UpdateLabels();
  
protected:
  SMESH_CellLabelActor();
  ~SMESH_CellLabelActor();

  bool myIsCellsLabeled;
  vtkUnstructuredGrid* myCellsNumDataSet;
  vtkActor2D *myCellsLabels;
  vtkMaskPoints* myClsMaskPoints;
  VTKViewer_CellCenters* myCellCenters;
  vtkLabeledDataMapper* myClsLabeledDataMapper;
  vtkSelectVisiblePoints* myClsSelectVisiblePoints;  
  SMESH_DeviceActor* myBaseActor; //Pointer to the base actor
  vtkTextProperty* myClsTextProp;

protected:
  // Not implemented.
  SMESH_CellLabelActor(const SMESH_CellLabelActor&);
  void operator=(const SMESH_CellLabelActor&);
};

#endif //SMESH_NODE_LABEL_ACTOR_H
