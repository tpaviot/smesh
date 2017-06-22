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

// ---
// File    : GeomSelectionTools.h
// Authors : Nicolas GEIMER (OCC)
// ---

#ifndef _GEOMSELECTIONTOOLS_H_
#define _GEOMSELECTIONTOOLS_H_

#include "SMESH_PluginUtils.h"

#include <SALOMEDSClient.hxx>
#include <SALOME_InteractiveObject.hxx>
#include <SALOME_ListIO.hxx>
#include <SalomeApp_Application.h>

#include <TopoDS_Shape.hxx>
#include <GeomAbs_SurfaceType.hxx>

class LightApp_SelectionMgr;


/*!
 * The GeomSelectionTools class gives high level tools to select Geom (and other objects)
 * A specific attention has been given to analyze selected GEOM objects.
 *
 * @param myStudy This class is specific to the study !
 *
 */

class PLUGINUTILS_EXPORT GeomSelectionTools
{

private:

  _PTR(Study) myStudy;

public:

  GeomSelectionTools(_PTR(Study));
  static SalomeApp_Application*  GetSalomeApplication();
  static LightApp_SelectionMgr* selectionMgr();
  SALOME_ListIO* getSelectedSalomeObjects();
  Handle(SALOME_InteractiveObject) getFirstSelectedSalomeObject();
  std::string getFirstSelectedEntry();
  std::string getEntryOfObject(Handle(SALOME_InteractiveObject));
  std::string getNameFromEntry(std::string);
  std::string getFirstSelectedComponentDataType();
  TopAbs_ShapeEnum getFirstSelectedShapeType();
  TopAbs_ShapeEnum entryToShapeType(std::string );
  GeomAbs_SurfaceType getFaceInformation(TopoDS_Shape);
  _PTR(Study) getMyStudy();
};

//////////////////////////////////////////
// Utility functions
//////////////////////////////////////////

namespace PluginUtils
{
  PLUGINUTILS_EXPORT QString PrintDoubleValue( double, int = 16 );
};

#endif // _GEOMSELECTIONTOOLS_H_

