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

//  SMESH SCALAR BAR : 2D Actor for the visualization scalar bar with the distribution diagram
//                     it is customized vtkScalarBarActor.
//  File   : SMESH_ScalarBarActor.h
//  Author : Roman NIKOLAEV
//  Module : SMESH


// .NAME vtkScalarBarActor - Create a scalar bar with labels 
// .SECTION Description
// vtkScalarBarActor creates a scalar bar with annotation text. A scalar
// bar is a legend that indicates to the viewer the correspondence between
// color value and data value. The legend consists of a rectangular bar 
// made of rectangular pieces each colored a constant value. Since 
// vtkScalarBarActor is a subclass of vtkActor2D, it is drawn in the image 
// plane (i.e., in the renderer's viewport) on top of the 3D graphics window.
//
// To use vtkScalarBarActor you must associate a vtkScalarsToColors (or
// subclass) with it. The lookup table defines the colors and the
// range of scalar values used to map scalar data.  Typically, the
// number of colors shown in the scalar bar is not equal to the number
// of colors in the lookup table, in which case sampling of
// the lookup table is performed. 
//
// Other optional capabilities include specifying the fraction of the
// viewport size (both x and y directions) which will control the size
// of the scalar bar and the number of annotation labels. The actual position
// of the scalar bar on the screen is controlled by using the
// vtkActor2D::SetPosition() method (by default the scalar bar is
// centered in the viewport).  Other features include the ability to
// orient the scalar bar horizontally of vertically and controlling
// the format (printf style) with which to print the labels on the
// scalar bar. Also, the vtkScalarBarActor's property is applied to
// the scalar bar and annotation (including layer, and
// compositing operator).
//
// Set the text property/attributes of the title and the labels through the 
// vtkTextProperty objects associated to this actor.
//
// .SECTION Caveats
// If a vtkLogLookupTable is specified as the lookup table to use, then the
// labels are created using a logarithmic scale.
//
// .SECTION See Also
// vtkActor2D vtkTextProperty vtkTextMapper vtkPolyDataMapper2D

#ifndef SMESH_SCALAR_BAR_ACTOR_H
#define SMESH_SCALAR_BAR_ACTOR_H

#include <vtkActor2D.h>

#include "SMESH_Object.h"

#include <vector>

class vtkPolyData;
class vtkPolyDataMapper2D;
class vtkScalarsToColors;
class vtkTextMapper;
class vtkTextProperty;

#define VTK_ORIENT_HORIZONTAL 0
#define VTK_ORIENT_VERTICAL 1

#define SMESH_MONOCOLOR_TYPE 0
#define SMESH_MULTICOLOR_TYPE 1


class SMESHOBJECT_EXPORT SMESH_ScalarBarActor: public vtkActor2D {
 public:
  void PrintSelf(ostream& os, vtkIndent indent);

  vtkTypeMacro(SMESH_ScalarBarActor,vtkActor2D);

  // Description:
  // Instantiate object with 64 maximum colors; 5 labels; %%-#6.3g label
  // format, no title, and vertical orientation. The initial scalar bar
  // size is (0.05 x 0.8) of the viewport size.
  static SMESH_ScalarBarActor *New();

  // Description:
  // Draw the scalar bar and annotation text to the screen.
  int RenderOpaqueGeometry(vtkViewport* viewport);
  int RenderTranslucentGeometry(vtkViewport*) { return 0; };
  int RenderOverlay(vtkViewport* viewport);

  // Description:
  // Release any graphics resources that are being consumed by this actor.
  // The parameter window could be used to determine which graphic
  // resources to release.
  virtual void ReleaseGraphicsResources(vtkWindow *);

  // Description:
  // Set/Get the vtkLookupTable to use. The lookup table specifies the number
  // of colors to use in the table (if not overridden), as well as the scalar
  // range.
  virtual void SetLookupTable(vtkScalarsToColors*);
  vtkGetObjectMacro(LookupTable,vtkScalarsToColors);

  // Description:
  // Set/Get the maximum number of scalar bar segments to show. This may
  // differ from the number of colors in the lookup table, in which case
  // the colors are samples from the lookup table.
  vtkSetClampMacro(MaximumNumberOfColors, int, 2, VTK_INT_MAX);
  vtkGetMacro(MaximumNumberOfColors, int);
  
  // Description:
  // Set/Get the number of annotation labels to show.
  vtkSetClampMacro(NumberOfLabels, int, 0, 64);
  vtkGetMacro(NumberOfLabels, int);
  
  // Description:
  // Control the orientation of the scalar bar.
  vtkSetClampMacro(Orientation,int,VTK_ORIENT_HORIZONTAL, VTK_ORIENT_VERTICAL);
  vtkGetMacro(Orientation, int);
  void SetOrientationToHorizontal()
       {this->SetOrientation(VTK_ORIENT_HORIZONTAL);};
  void SetOrientationToVertical() {this->SetOrientation(VTK_ORIENT_VERTICAL);};

  // Description:
  // Set/Get the title text property.
  virtual void SetTitleTextProperty(vtkTextProperty *p);
  vtkGetObjectMacro(TitleTextProperty,vtkTextProperty);
  
  // Description:
  // Set/Get the labels text property.
  virtual void SetLabelTextProperty(vtkTextProperty *p);
  vtkGetObjectMacro(LabelTextProperty,vtkTextProperty);
    
  // Description:
  // Set/Get the format with which to print the labels on the scalar
  // bar.
  vtkSetStringMacro(LabelFormat);
  vtkGetStringMacro(LabelFormat);

  // Description:
  // Set/Get the title of the scalar bar actor,
  vtkSetStringMacro(Title);
  vtkGetStringMacro(Title);

  // Description:
  // Shallow copy of a scalar bar actor. Overloads the virtual vtkProp method.
  void ShallowCopy(vtkProp *prop);

  // Description:
  // Set visibility of the distribution histogram
  // rnv: Customization of the vtkScalarBarActor to show distribution histogram:
  virtual void SetDistributionVisibility(int flag);

  // Description:
  // Set visibility of the distribution histogram
  // rnv: Customization of the vtkScalarBarActor to show distribution histogram:
  virtual int GetDistributionVisibility();
  // Description:
  // Set distribution
  virtual void SetDistribution(std::vector<int> theNbValues);
  
  // Description: 
  // Set distribution coloring type (SMESH_MONOCOLOR_TYPE or SMESH_MULTICOLOR_TYPE)
  void SetDistributionColoringType(int theDistributionColoringType) {myDistributionColoringType = theDistributionColoringType;Modified();}

  // Description: 
  // Get distribution coloring type ((SMESH_MONOCOLOR_TYPE or SMESH_MULTICOLOR_TYPE))  
  int GetDistributionColoringType() {return myDistributionColoringType;}

  // Description:
  // Set Distribution Color
  void SetDistributionColor (double rgb[3]);

  // Description:
  // Get Distribution Color
  void GetDistributionColor (double rgb[3]);
  
  // Description:
  // Set visibility status of scalar map 
  void SetTitleOnlyVisibility( bool );

  // Description:
  // Get visibility status of scalar map 
  bool GetTitleOnlyVisibility();

 protected:
  SMESH_ScalarBarActor();
  ~SMESH_ScalarBarActor();
  
  vtkScalarsToColors *LookupTable;
  vtkTextProperty *TitleTextProperty;
  vtkTextProperty *LabelTextProperty;
  
  int   MaximumNumberOfColors;
  int   NumberOfLabels;
  int   NumberOfLabelsBuilt;
  int   Orientation;
  char  *Title;
  char  *LabelFormat;

  vtkTextMapper **TextMappers;
  virtual void AllocateAndSizeLabels(int *labelSize, int *size,
                                     vtkViewport *viewport, double *range);

  
  
 private:
  vtkTextMapper *TitleMapper;
  vtkActor2D    *TitleActor;

  vtkActor2D    **TextActors;

  vtkPolyData         *ScalarBar;
  vtkPolyDataMapper2D *ScalarBarMapper;
  vtkActor2D          *ScalarBarActor;

  vtkTimeStamp  BuildTime;
  int LastSize[2];
  int LastOrigin[2];

  void SizeTitle(int *titleSize, int *size, vtkViewport *viewport);
  
  // rnv: Customization of the vtkScalarBarActor to show distribution histogram:
  vtkPolyData*           myDistribution;             //Distribution polygonal data
  vtkActor2D*            myDistributionActor;        //Distribution actor
  vtkPolyDataMapper2D*   myDistributionMapper;       //Distribution mapper
  std::vector<int>       myNbValues;                 //Nb values for the range
  int                    myDistributionColoringType; //Distribution color type (monocolor or multicolor)
  bool                   myTitleOnlyVisibility;      //Show scalar map or not
  
 private:
  SMESH_ScalarBarActor(const SMESH_ScalarBarActor&);  // Not implemented.
  void operator=(const SMESH_ScalarBarActor&);  // Not implemented.
};

#endif //SMESH_SCALAR_BAR_ACTOR_H
