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

#include "SMESH_ActorUtils.h"
#include "SMESH_Actor.h"

#include "SUIT_Tools.h"
#include "SUIT_Session.h"
#include "SUIT_ResourceMgr.h"
#include <SALOMEconfig.h> // To fix some redefinition
#include "SalomeApp_Application.h"

#ifndef DISABLE_PLOT2DVIEWER
#include <SPlot2d_ViewModel.h>
#include <SPlot2d_Histogram.h>
#include <Plot2d_ViewManager.h>
#endif

#include <Qtx.h>


#include "utilities.h"

#include <vtkUnstructuredGrid.h>
#include <vtkCellType.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkUnsignedCharArray.h>

//#ifdef _DEBUG_
//static int MYDEBUG = 1;
//#else
//static int MYDEBUG = 0;
//#endif

namespace SMESH
{

  double
  GetFloat( const QString& theValue, 
            double theDefault )
  {
    int pos = theValue.indexOf( ":" );
    double val = theDefault;
    if( pos>=0 ) 
    {
      QString name = theValue.right( theValue.length()-pos-1 ),
              sect = theValue.left( pos );
      if( !name.isEmpty() && !sect.isEmpty() )
        val = GetFloat( name, sect, theDefault );
    }
    return val;
  }

  double
  GetFloat( const QString& theValue, 
            const QString& theSection, 
            double theDefault )
  {
    double val = theDefault;
    SUIT_ResourceMgr* mgr = SUIT_Session::session()->resourceMgr();
    if( mgr )
      val = (double) mgr->doubleValue( theSection, theValue, theDefault );

    return val;
  }

  void
  WriteUnstructuredGrid(vtkUnstructuredGrid* theGrid, 
                        const char* theFileName)
  {
    vtkXMLUnstructuredGridWriter* aWriter = vtkXMLUnstructuredGridWriter::New();
    aWriter->SetFileName(theFileName);
    aWriter->SetInputData(theGrid);
    aWriter->SetDataModeToAscii();
    if(theGrid->GetNumberOfCells()){
      aWriter->Write();
    }
    aWriter->Delete();
  }

  QColor
  GetColor( const QString& theSect, 
            const QString& theName, 
            const QColor& def )
  {
    QColor c = def;
    SUIT_ResourceMgr* mgr = SUIT_Session::session()->resourceMgr();
    if ( mgr )
      c = mgr->colorValue( theSect, theName, def );
    return c;
  }

  void
  GetColor( const QString& theSect, 
            const QString& theName, 
            int& r, 
            int& g, 
            int& b, 
            const QColor& def )
  {
    QColor c = def;
    SUIT_ResourceMgr* mgr = SUIT_Session::session()->resourceMgr();
    if ( mgr )
      c = mgr->colorValue( theSect, theName, def );

    SUIT_Tools::rgbSet( SUIT_Tools::rgbSet( c ), r, g, b );
  }

  void
  GetColor( const QString& theSect, 
            const QString& theName, 
            double& r, 
            double& g, 
            double& b, 
            const QColor& def )
  {
    int ir( 0 ), ig( 0 ), ib( 0 );
    GetColor( theSect, theName, ir, ig, ib, def );
    r = ir / 255.;
    g = ig / 255.;
    b = ib / 255.;
  }


  void
  GetColor(  const QString& theSect, 
             const QString& theName, 
             QColor& color,
             int& delta,
             QString def) 
  {
    
    SUIT_ResourceMgr* mgr = SUIT_Session::session()->resourceMgr();
    if ( mgr ) {
      QString str = mgr->stringValue( theSect, theName, def );
      Qtx::stringToBiColor(str,color,delta);
    }
  }

  std::map<SMDSAbs_ElementType,int> GetEntitiesFromObject(SMESH_VisualObj *theObject) {
    std::map<SMDSAbs_ElementType,int> entities;
    entities.insert(std::pair<SMDSAbs_ElementType,int>(SMDSAbs_0DElement,
                theObject ? theObject->GetNbEntities(SMDSAbs_0DElement) : 0));
    entities.insert(std::pair<SMDSAbs_ElementType,int>(SMDSAbs_Ball,
                theObject ? theObject->GetNbEntities(SMDSAbs_Ball) : 0));
    entities.insert(std::pair<SMDSAbs_ElementType,int>(SMDSAbs_Edge,
                theObject ? theObject->GetNbEntities(SMDSAbs_Edge) : 0));
    entities.insert(std::pair<SMDSAbs_ElementType,int>(SMDSAbs_Face,
                theObject ? theObject->GetNbEntities(SMDSAbs_Face) : 0));
    entities.insert(std::pair<SMDSAbs_ElementType,int>(SMDSAbs_Volume,
                theObject ? theObject->GetNbEntities(SMDSAbs_Volume) : 0));
    return entities;
  }
  


#ifndef DISABLE_PLOT2DVIEWER
  //=======================================================================
  /**
     Get histogram from the input actor
     Repaint/Remove the histogram in/from each opened Plot2D Viewer 
  */
  //=======================================================================
  void ProcessIn2DViewers( SMESH_Actor *theActor, Viewer2dActionType aType ) {
    SalomeApp_Application* anApp = dynamic_cast<SalomeApp_Application*>(SUIT_Session::session()->activeApplication());
    
    if(!anApp || !theActor)
      return;
    
    SPlot2d_Histogram* aHistogram = 0;
    
    if(theActor->GetPlot2Histogram())
      if(aType == UpdateIn2dViewer)
        aHistogram = theActor->UpdatePlot2Histogram();
      else
        aHistogram = theActor->GetPlot2Histogram();
    else 
      return;
    
    ViewManagerList aViewManagerList;
    anApp->viewManagers(SPlot2d_Viewer::Type(), aViewManagerList);
    
    aType = aHistogram->getPointList().empty() ? RemoveFrom2dViewer : aType;
    
    SUIT_ViewManager* aViewManager;
    foreach( aViewManager, aViewManagerList ) {
      if (Plot2d_ViewManager* aManager = dynamic_cast<Plot2d_ViewManager*>(aViewManager)) {
        if (SPlot2d_Viewer* aViewer = dynamic_cast<SPlot2d_Viewer*>(aManager->getViewModel())) {
          if (Plot2d_ViewFrame* aViewFrame = aViewer->getActiveViewFrame()) {
            if(aType == UpdateIn2dViewer )
              aViewFrame->displayObject(aHistogram, true);
            else if (aType == RemoveFrom2dViewer)
              aViewFrame->eraseObject(aHistogram, true);
          }
        }
      }
    }
  }
#endif //DISABLE_PLOT2DVIEWER
  
}
