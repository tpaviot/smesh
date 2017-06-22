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

#ifndef SMESH_ACTORUTILS_H
#define SMESH_ACTORUTILS_H

#include "SMESH_Object.h"
#include <map>

#include <QColor>

class vtkUnstructuredGrid;
class SMESH_Actor;

namespace SMESH
{
  //! Label font family
  enum LabelFont { 
    FntArial,    //!< arial font family
    FntCourier,  //!< courier font family
    FntTimes     //!< times font family
  };

SMESHOBJECT_EXPORT  
  double 
  GetFloat( const QString& theValue, 
            double theDefault = 0 );

SMESHOBJECT_EXPORT
  double 
  GetFloat( const QString& theName, 
            const QString& theSection, 
            double theDefault = 0 );

SMESHOBJECT_EXPORT
  QColor 
  GetColor( const QString& theSect, 
            const QString& theName, 
            const QColor& = QColor() );

SMESHOBJECT_EXPORT
  void
  GetColor( const QString& theSect, 
            const QString& theName, 
            int&, 
            int&, 
            int&, 
            const QColor& = QColor() );

SMESHOBJECT_EXPORT
  void
  GetColor( const QString& theSect, 
            const QString& theName, 
            double&, 
            double&, 
            double&, 
            const QColor& = QColor() );

 SMESHOBJECT_EXPORT
   void
   GetColor(  const QString& theSect, 
              const QString& theName, 
              QColor& color,
              int& delta,
              QString def);
   
 SMESHOBJECT_EXPORT
   std::map<SMDSAbs_ElementType,int>
   GetEntitiesFromObject(SMESH_VisualObj *theObject);
   
SMESHOBJECT_EXPORT
  void 
  WriteUnstructuredGrid(vtkUnstructuredGrid* theGrid, 
                        const char* theFileName);


#ifndef DISABLE_PLOT2DVIEWER
 
 typedef enum {UpdateIn2dViewer = 0, RemoveFrom2dViewer } Viewer2dActionType;
 
 SMESHOBJECT_EXPORT
   void ProcessIn2DViewers( SMESH_Actor* theActor, Viewer2dActionType = UpdateIn2dViewer );
 
#endif


}

#endif
