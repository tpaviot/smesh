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
//  File   : SMESH_ScalarBarActor.cxx
//  Author : Roman NIKOLAEV
//  Module : SMESH
//

#include "SMESH_ScalarBarActor.h"

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkObjectFactory.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper2D.h>
#include <vtkScalarsToColors.h>
#include <vtkTextMapper.h>
#include <vtkTextProperty.h>
#include <vtkViewport.h>
#include <vtkWindow.h>
#include <vtkLookupTable.h>
#include <vtkProperty2D.h>

#define SHRINK_COEF 0.08;

vtkStandardNewMacro(SMESH_ScalarBarActor);

vtkCxxSetObjectMacro(SMESH_ScalarBarActor,LookupTable,vtkScalarsToColors);
vtkCxxSetObjectMacro(SMESH_ScalarBarActor,LabelTextProperty,vtkTextProperty);
vtkCxxSetObjectMacro(SMESH_ScalarBarActor,TitleTextProperty,vtkTextProperty);

//----------------------------------------------------------------------------
// Instantiate object with 64 maximum colors; 5 labels; %%-#6.3g label
// format, no title, and vertical orientation. The initial scalar bar
// size is (0.05 x 0.8) of the viewport size.
SMESH_ScalarBarActor::SMESH_ScalarBarActor() {
  this->LookupTable = NULL;
  this->Position2Coordinate->SetValue(0.17, 0.8);
  
  this->PositionCoordinate->SetCoordinateSystemToNormalizedViewport();
  this->PositionCoordinate->SetValue(0.82,0.1);
  
  this->MaximumNumberOfColors = 64;
  this->NumberOfLabels = 5;
  this->NumberOfLabelsBuilt = 0;
  this->Orientation = VTK_ORIENT_VERTICAL;
  this->Title = NULL;

  this->LabelTextProperty = vtkTextProperty::New();
  this->LabelTextProperty->SetFontSize(12);
  this->LabelTextProperty->SetBold(1);
  this->LabelTextProperty->SetItalic(1);
  this->LabelTextProperty->SetShadow(1);
  this->LabelTextProperty->SetFontFamilyToArial();

  this->TitleTextProperty = vtkTextProperty::New();
  this->TitleTextProperty->ShallowCopy(this->LabelTextProperty);

  this->LabelFormat = new char[8]; 
  sprintf(this->LabelFormat,"%s","%-#6.3g");

  this->TitleMapper = vtkTextMapper::New();
  this->TitleActor = vtkActor2D::New();
  this->TitleActor->SetMapper(this->TitleMapper);
  this->TitleActor->GetPositionCoordinate()->
    SetReferenceCoordinate(this->PositionCoordinate);
  
  this->TextMappers = NULL;
  this->TextActors = NULL;

  this->ScalarBar = vtkPolyData::New();
  this->ScalarBarMapper = vtkPolyDataMapper2D::New();
  this->ScalarBarMapper->SetInputData(this->ScalarBar);
  this->ScalarBarActor = vtkActor2D::New();
  this->ScalarBarActor->SetMapper(this->ScalarBarMapper);
  this->ScalarBarActor->GetPositionCoordinate()->
    SetReferenceCoordinate(this->PositionCoordinate);
  this->LastOrigin[0] = 0;
  this->LastOrigin[1] = 0;
  this->LastSize[0] = 0;
  this->LastSize[1] = 0;


  // rnv begin
  // Customization of the vtkScalarBarActor to show distribution histogram.
  myDistribution = vtkPolyData::New();
  myDistributionMapper = vtkPolyDataMapper2D::New();
  myDistributionMapper->SetInputData(this->myDistribution);
  
  myDistributionActor = vtkActor2D::New();
  myDistributionActor->SetMapper(this->myDistributionMapper);
  myDistributionActor->GetPositionCoordinate()->
    SetReferenceCoordinate(this->PositionCoordinate);

  // By default distribution histogram is invisible
  myDistributionActor->SetVisibility(0);

  // By default monocolor
  myDistributionColoringType = SMESH_MONOCOLOR_TYPE;

  // By default scalar map is shown
  myTitleOnlyVisibility = false;
  // rnv end
}

//----------------------------------------------------------------------------
// Release any graphics resources that are being consumed by this actor.
// The parameter window could be used to determine which graphic
// resources to release.
void SMESH_ScalarBarActor::ReleaseGraphicsResources(vtkWindow *win)
{
  this->TitleActor->ReleaseGraphicsResources(win);
  if (this->TextMappers != NULL )
    {
    for (int i=0; i < this->NumberOfLabelsBuilt; i++)
      {
      this->TextActors[i]->ReleaseGraphicsResources(win);
      }
    }
  this->ScalarBarActor->ReleaseGraphicsResources(win);
  // rnv begin
  // Customization of the vtkScalarBarActor to show distribution histogram.
  myDistributionActor->ReleaseGraphicsResources(win);
}


/*--------------------------------------------------------------------------*/
SMESH_ScalarBarActor::~SMESH_ScalarBarActor() {
  if (this->LabelFormat) 
    {
    delete [] this->LabelFormat;
    this->LabelFormat = NULL;
    }

  this->TitleMapper->Delete();
  this->TitleActor->Delete();

  if (this->TextMappers != NULL )
    {
    for (int i=0; i < this->NumberOfLabelsBuilt; i++)
      {
      this->TextMappers[i]->Delete();
      this->TextActors[i]->Delete();
      }
    delete [] this->TextMappers;
    delete [] this->TextActors;
    }

  this->ScalarBar->Delete();
  this->ScalarBarMapper->Delete();
  this->ScalarBarActor->Delete();

  if (this->Title)
    {
    delete [] this->Title;
    this->Title = NULL;
    }
  
  this->SetLookupTable(NULL);
  this->SetLabelTextProperty(NULL);
  this->SetTitleTextProperty(NULL);
  
  // rnv begin
  // Customization of the vtkScalarBarActor to show distribution histogram:
  myDistribution->Delete();
  myDistributionMapper->Delete();
  myDistributionActor->Delete();
  // rnv end
}

//----------------------------------------------------------------------------
int SMESH_ScalarBarActor::RenderOverlay(vtkViewport *viewport)
{
  int renderedSomething = 0;
  int i;
  
  // Everything is built, just have to render
  if (this->Title != NULL)
    {
    renderedSomething += this->TitleActor->RenderOverlay(viewport);
    }
  if (!myTitleOnlyVisibility) {
    this->ScalarBarActor->RenderOverlay(viewport);
    this->myDistributionActor->RenderOverlay(viewport);
    if( this->TextActors == NULL)
      {
       vtkWarningMacro(<<"Need a mapper to render a scalar bar");
       return renderedSomething;
      }
  
    for (i=0; i<this->NumberOfLabels; i++)
      {
      renderedSomething += this->TextActors[i]->RenderOverlay(viewport);
      }
  }  
  renderedSomething = (renderedSomething > 0)?(1):(0);

  return renderedSomething;
}


//----------------------------------------------------------------------------
int SMESH_ScalarBarActor::RenderOpaqueGeometry(vtkViewport *viewport)
{
  int renderedSomething = 0;
  int i;
  int size[2];
  
  if (!this->LookupTable)
    {
    vtkWarningMacro(<<"Need a mapper to render a scalar bar");
    return 0;
    }

  if (!this->TitleTextProperty)
    {
    vtkErrorMacro(<<"Need title text property to render a scalar bar");
    return 0;
    }

  if (!this->LabelTextProperty)
    {
    vtkErrorMacro(<<"Need label text property to render a scalar bar");
    return 0;
    }

  // Check to see whether we have to rebuild everything
  int positionsHaveChanged = 0;
  if (viewport->GetMTime() > this->BuildTime || 
      (viewport->GetVTKWindow() && 
       viewport->GetVTKWindow()->GetMTime() > this->BuildTime))
    {
    // if the viewport has changed we may - or may not need
    // to rebuild, it depends on if the projected coords chage
    int *barOrigin;
    barOrigin = this->PositionCoordinate->GetComputedViewportValue(viewport);
    size[0] = 
      this->Position2Coordinate->GetComputedViewportValue(viewport)[0] -
      barOrigin[0];
    size[1] = 
      this->Position2Coordinate->GetComputedViewportValue(viewport)[1] -
      barOrigin[1];
    if (this->LastSize[0] != size[0] || 
        this->LastSize[1] != size[1] ||
        this->LastOrigin[0] != barOrigin[0] || 
        this->LastOrigin[1] != barOrigin[1])
      {
      positionsHaveChanged = 1;
      }
    }
  
  // Check to see whether we have to rebuild everything
  if (positionsHaveChanged ||
      this->GetMTime() > this->BuildTime || 
      this->LookupTable->GetMTime() > this->BuildTime ||
      this->LabelTextProperty->GetMTime() > this->BuildTime ||
      this->TitleTextProperty->GetMTime() > this->BuildTime)
    {
    vtkDebugMacro(<<"Rebuilding subobjects");

    // Delete previously constructed objects
    //
    if (this->TextMappers != NULL )
      {
      for (i=0; i < this->NumberOfLabelsBuilt; i++)
        {
        this->TextMappers[i]->Delete();
        this->TextActors[i]->Delete();
        }
      delete [] this->TextMappers;
      delete [] this->TextActors;
      }

    // Build scalar bar object; determine its type
    //
    // is this a vtkLookupTable or a subclass of vtkLookupTable 
    // with its scale set to log
    // NOTE: it's possible we could to without the 'lut' variable
    // later in the code, but if the vtkLookupTableSafeDownCast operation
    // fails for some reason, this code will break in new ways. So, the 'LUT'
    // variable is used for this operation only
    vtkLookupTable *LUT = vtkLookupTable::SafeDownCast( this->LookupTable );
    int isLogTable = 0;
    if ( LUT )
      {
      if ( LUT->GetScale() == VTK_SCALE_LOG10 )
        {
        isLogTable = 1; 
        }
      }
    
    // we hard code how many steps to display
    vtkScalarsToColors *lut = this->LookupTable;
    int numColors = this->MaximumNumberOfColors;
    double *range = lut->GetRange();

    int numPts = 2*(numColors + 1);
    vtkPoints *pts = vtkPoints::New();
    pts->SetNumberOfPoints(numPts);
    vtkCellArray *polys = vtkCellArray::New();
    polys->Allocate(polys->EstimateSize(numColors,4));
    vtkUnsignedCharArray *colors = vtkUnsignedCharArray::New();
    colors->SetNumberOfComponents(3);
    colors->SetNumberOfTuples(numColors);


    // rnv begin
    // Customization of the vtkScalarBarActor to show distribution histogram.
    bool distrVisibility =  (numColors == (int)this->myNbValues.size());
    vtkPoints *distrPts = 0;
    vtkCellArray *distrPolys = 0;
    vtkUnsignedCharArray *distColors = 0;
    int numDistrPts = 0, numPositiveVal=0, maxValue=0;
    if(!distrVisibility)
      vtkDebugMacro(<<" Distribution invisible, because numColors == this->myNbValues.size()");

    if ( distrVisibility && GetDistributionVisibility() ) {
      for ( i = 0 ; i < (int)myNbValues.size(); i++ ) {
        if ( myNbValues[i] ) {
          numPositiveVal++;
          maxValue = std::max(maxValue,myNbValues[i]);
        }
      }
      numDistrPts = 4*(numPositiveVal);
      distrPts = vtkPoints::New();
      distrPolys = vtkCellArray::New();
      distrPts->SetNumberOfPoints(numDistrPts);
      distrPolys->Allocate(distrPolys->EstimateSize(numPositiveVal,4));
      this->myDistribution->Initialize();
      this->myDistribution->SetPoints(distrPts);
      this->myDistribution->SetPolys(distrPolys);
      distrPts->Delete();
      distrPolys->Delete();
      if ( myDistributionColoringType == SMESH_MULTICOLOR_TYPE ) {
        distColors = vtkUnsignedCharArray::New();
        distColors->SetNumberOfComponents(3);
        distColors->SetNumberOfTuples(numPositiveVal);
        this->myDistribution->GetCellData()->SetScalars(distColors);
        distColors->Delete();
      } else if( myDistributionColoringType == SMESH_MONOCOLOR_TYPE ){
        this->myDistribution->GetCellData()->SetScalars(NULL);
      }
    } else {
      myDistribution->Reset();
    }
    // rnv end

    this->ScalarBarActor->SetProperty(this->GetProperty());
    this->ScalarBar->Initialize();
    this->ScalarBar->SetPoints(pts);
    this->ScalarBar->SetPolys(polys);
    this->ScalarBar->GetCellData()->SetScalars(colors);
    pts->Delete(); polys->Delete(); colors->Delete();

    // get the viewport size in display coordinates
    int *barOrigin, barWidth, barHeight, distrHeight;
    barOrigin = this->PositionCoordinate->GetComputedViewportValue(viewport);
    size[0] = 
      this->Position2Coordinate->GetComputedViewportValue(viewport)[0] -
      barOrigin[0];
    size[1] = 
      this->Position2Coordinate->GetComputedViewportValue(viewport)[1] -
      barOrigin[1];
    this->LastOrigin[0] = barOrigin[0];
    this->LastOrigin[1] = barOrigin[1];
    this->LastSize[0] = size[0];
    this->LastSize[1] = size[1];
    
    // Update all the composing objects
    this->TitleActor->SetProperty(this->GetProperty());
    this->TitleMapper->SetInput(this->Title);
    if (this->TitleTextProperty->GetMTime() > this->BuildTime)
      {
      // Shallow copy here so that the size of the title prop is not affected
      // by the automatic adjustment of its text mapper's size (i.e. its
      // mapper's text property is identical except for the font size
      // which will be modified later). This allows text actors to
      // share the same text property, and in that case specifically allows
      // the title and label text prop to be the same.
      this->TitleMapper->GetTextProperty()->ShallowCopy(this->TitleTextProperty);
      this->TitleMapper->GetTextProperty()->SetJustificationToCentered();
      }
    
    // find the best size for the title font
    int titleSize[2];
    this->SizeTitle(titleSize, size, viewport);
    
    // find the best size for the ticks
    int labelSize[2];
    this->AllocateAndSizeLabels(labelSize, size, viewport,range);
    this->NumberOfLabelsBuilt = this->NumberOfLabels;
    
    // generate points
    double x[3]; x[2] = 0.0;
    double delta, itemH, shrink;
    if ( this->Orientation == VTK_ORIENT_VERTICAL ) {
      // rnv begin
      // Customization of the vtkScalarBarActor to show distribution histogram.
      double delimeter=0.0;
      if(GetDistributionVisibility() && distrVisibility) {
        delimeter=0.01*size[0]; //1 % from horizontal size of the full presentation size.
        barWidth = size[0] - 4 - labelSize[0];
        distrHeight = barWidth/2;
      } else {
        barWidth = size[0] - 4 - labelSize[0];
        distrHeight = 0;
      }

      barHeight = (int)(0.86*size[1]);
      delta=(double)barHeight/numColors;
      
      for ( i=0; i<numPts/2; i++ ) {
        x[0] = distrHeight+delimeter/2.0;
        x[1] = i*delta;
        pts->SetPoint(2*i,x);
        x[0] = barWidth;
        pts->SetPoint(2*i+1,x);
      }

      if(GetDistributionVisibility() && distrVisibility) {
        // Distribution points 
        shrink = delta*SHRINK_COEF;
        vtkIdType distPtsId=0;
        vtkIdType distPtsIds[4];
        for(i=0; i<numColors; i++) {
          if(myNbValues[i]) {
            itemH = distrHeight*((double)myNbValues[i]/maxValue);
            
            if(distrHeight == itemH) 
              itemH = itemH - delimeter/2;

            x[1] = i*delta+shrink;

            // first point of polygon (quadrangle)
            x[0] = 0; 
            distPtsIds[0] = distPtsId;
            distrPts->SetPoint(distPtsId++,x);

            // second point of polygon (quadrangle)
            x[0] = itemH;
            distPtsIds[1] = distPtsId;
            distrPts->SetPoint(distPtsId++,x);

            x[1] = i*delta+delta-shrink;

            // third point of polygon (quadrangle)
            x[0] = 0; 
            distPtsIds[3] = distPtsId;
            distrPts->SetPoint(distPtsId++,x);

            // fourth point of polygon (quadrangle)
            x[0] = itemH;
            distPtsIds[2] = distPtsId;
            distrPts->SetPoint(distPtsId++,x);

            //Inser Quadrangle
            distrPolys->InsertNextCell(4,distPtsIds);
          }
        }
      }    
    }
    // rnv end
    else {
      barWidth = size[0];
      
      // rnv begin
      // Customization of the vtkScalarBarActor to show distribution histogram.
      double coef1, delimeter=0.0;
      if(GetDistributionVisibility() && distrVisibility) {
        coef1=0.62;
        distrHeight = (int)((coef1/2)*size[1]);
        //delimeter between distribution diagram and scalar bar 
        delimeter=0.02*size[1];
      }
      else {
        coef1=0.4;
        barHeight = (int)(coef1*size[1]);
        distrHeight = 0;
      }
      
      barHeight = (int)(coef1*size[1]);
      
      delta=(double)barWidth/numColors;
      for (i=0; i<numPts/2; i++) {
        x[0] = i*delta;
        x[1] = barHeight;
        pts->SetPoint(2*i,x);                        
        x[1] = distrHeight + delimeter;
        pts->SetPoint(2*i+1,x);
      }
      
      if(GetDistributionVisibility() && distrVisibility) {
        // Distribution points 
        shrink = delta*SHRINK_COEF;
        vtkIdType distPtsId=0;
        vtkIdType distPtsIds[4];
        for(i=0; i<numColors; i++) {
          if(myNbValues[i]) {
            itemH = distrHeight*((double)myNbValues[i]/maxValue);
            
            // first point of polygon (quadrangle)
            x[0] = i*delta+shrink; 
            x[1] = 0;
            distPtsIds[0] = distPtsId;
            distrPts->SetPoint(distPtsId++,x);
            
            // second point of polygon (quadrangle)
            x[0] = i*delta+shrink; 
            x[1] = itemH;
            distPtsIds[3] = distPtsId;
            distrPts->SetPoint(distPtsId++,x);
            
            // third point of polygon (quadrangle)
            x[0] = i*delta+delta-shrink; 
            x[1] = 0;
            distPtsIds[1] = distPtsId;
            distrPts->SetPoint(distPtsId++,x);
            
            // fourth point of polygon (quadrangle)
            x[0] = i*delta+delta-shrink; 
            x[1] = itemH;
            distPtsIds[2] = distPtsId;
            distrPts->SetPoint(distPtsId++,x);
            
            // Add polygon into poly data
            distrPolys->InsertNextCell(4,distPtsIds);
          }
        } 
      }
      // rnv end
    }
    
    //polygons & cell colors
    unsigned char *rgba, *rgb;
    vtkIdType ptIds[4], dcCount=0;
    for (i=0; i<numColors; i++)
      {
      ptIds[0] = 2*i;
      ptIds[1] = ptIds[0] + 1;
      ptIds[2] = ptIds[1] + 2;
      ptIds[3] = ptIds[0] + 2;
      polys->InsertNextCell(4,ptIds);

      if ( isLogTable )
        {
        double rgbval = log10(range[0]) + 
          i*(log10(range[1])-log10(range[0]))/(numColors -1);
        rgba = lut->MapValue(pow(10.0,rgbval));
        }
      else
        {
        rgba = lut->MapValue(range[0] + (range[1] - range[0])*
                             ((double)i /(numColors-1.0)));
        }

      rgb = colors->GetPointer(3*i); //write into array directly
      rgb[0] = rgba[0];
      rgb[1] = rgba[1];
      rgb[2] = rgba[2];
      
      // rnv begin
      // Customization of the vtkScalarBarActor to show distribution histogram.
      if(myDistributionColoringType == SMESH_MULTICOLOR_TYPE && GetDistributionVisibility() && distrVisibility)
        {
          rgb = distColors->GetPointer(3*dcCount); //write into array directly
          rgb[0] = rgba[0];
          rgb[1] = rgba[1];
          rgb[2] = rgba[2];
          dcCount++;
        }
      }

    // Now position everything properly
    //
    double val;
    if (this->Orientation == VTK_ORIENT_VERTICAL)
      {
      int sizeTextData[2];
      
      // center the title
      this->TitleActor->SetPosition(size[0]/2, 0.9*size[1]);
      
      for (i=0; i < this->NumberOfLabels; i++)
        {
        if (this->NumberOfLabels > 1)
          {
          val = (double)i/(this->NumberOfLabels-1) *barHeight;
          }
        else 
          {
          val = 0.5*barHeight;
          }
        this->TextMappers[i]->GetSize(viewport,sizeTextData);
        this->TextMappers[i]->GetTextProperty()->SetJustificationToLeft();
        this->TextActors[i]->SetPosition(barWidth+3,
                                         val - sizeTextData[1]/2);
        }
      }
    else
      {
      this->TitleActor->SetPosition(size[0]/2, 
                                    barHeight + labelSize[1] + 0.1*size[1]);
      for (i=0; i < this->NumberOfLabels; i++)
        {
        this->TextMappers[i]->GetTextProperty()->SetJustificationToCentered();
        if (this->NumberOfLabels > 1)
          {
          val = (double)i/(this->NumberOfLabels-1) * barWidth;
          }
        else
          {
          val = 0.5*barWidth;
          }
        this->TextActors[i]->SetPosition(val, barHeight + 0.05*size[1]);
        }
      }

    this->BuildTime.Modified();
    }

  // Everything is built, just have to render
  if (this->Title != NULL)
    {
    renderedSomething += this->TitleActor->RenderOpaqueGeometry(viewport);
    }
  this->ScalarBarActor->RenderOpaqueGeometry(viewport);
  this->myDistributionActor->RenderOpaqueGeometry(viewport);
  for (i=0; i<this->NumberOfLabels; i++)
    {
    renderedSomething += this->TextActors[i]->RenderOpaqueGeometry(viewport);
    }

  renderedSomething = (renderedSomething > 0)?(1):(0);

  return renderedSomething;
}

//----------------------------------------------------------------------------
void SMESH_ScalarBarActor::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  if ( this->LookupTable )
    {
    os << indent << "Lookup Table:\n";
    this->LookupTable->PrintSelf(os,indent.GetNextIndent());
    }
  else
    {
    os << indent << "Lookup Table: (none)\n";
    }

  if (this->TitleTextProperty)
    {
    os << indent << "Title Text Property:\n";
    this->TitleTextProperty->PrintSelf(os,indent.GetNextIndent());
    }
  else
    {
    os << indent << "Title Text Property: (none)\n";
    }

  if (this->LabelTextProperty)
    {
    os << indent << "Label Text Property:\n";
    this->LabelTextProperty->PrintSelf(os,indent.GetNextIndent());
    }
  else
    {
    os << indent << "Label Text Property: (none)\n";
    }

  os << indent << "Title: " << (this->Title ? this->Title : "(none)") << "\n";
  os << indent << "Maximum Number Of Colors: " 
     << this->MaximumNumberOfColors << "\n";
  os << indent << "Number Of Labels: " << this->NumberOfLabels << "\n";
  os << indent << "Number Of Labels Built: " << this->NumberOfLabelsBuilt << "\n";

  os << indent << "Orientation: ";
  if ( this->Orientation == VTK_ORIENT_HORIZONTAL )
    {
    os << "Horizontal\n";
    }
  else
    {
    os << "Vertical\n";
    }

  os << indent << "Label Format: " << this->LabelFormat << "\n";
}

//----------------------------------------------------------------------------
void SMESH_ScalarBarActor::ShallowCopy(vtkProp *prop)
{
  SMESH_ScalarBarActor *a = SMESH_ScalarBarActor::SafeDownCast(prop);
  if ( a != NULL )
    {
    this->SetPosition2(a->GetPosition2());
    this->SetLookupTable(a->GetLookupTable());
    this->SetMaximumNumberOfColors(a->GetMaximumNumberOfColors());
    this->SetOrientation(a->GetOrientation());
    this->SetLabelTextProperty(a->GetLabelTextProperty());
    this->SetTitleTextProperty(a->GetTitleTextProperty());
    this->SetLabelFormat(a->GetLabelFormat());
    this->SetTitle(a->GetTitle());
    this->GetPositionCoordinate()->SetCoordinateSystem(
      a->GetPositionCoordinate()->GetCoordinateSystem());    
    this->GetPositionCoordinate()->SetValue(
      a->GetPositionCoordinate()->GetValue());
    this->GetPosition2Coordinate()->SetCoordinateSystem(
      a->GetPosition2Coordinate()->GetCoordinateSystem());    
    this->GetPosition2Coordinate()->SetValue(
      a->GetPosition2Coordinate()->GetValue());
    }

  // Now do superclass
  this->vtkActor2D::ShallowCopy(prop);
}

//----------------------------------------------------------------------------
void SMESH_ScalarBarActor::AllocateAndSizeLabels(int *labelSize, 
                                              int *size,
                                              vtkViewport *viewport,
                                              double *range)
{
  labelSize[0] = labelSize[1] = 0;

  this->TextMappers = new vtkTextMapper * [this->NumberOfLabels];
  this->TextActors = new vtkActor2D * [this->NumberOfLabels];

  char string[512];

  double val;
  int i;
  
  // TODO: this should be optimized, maybe by keeping a list of
  // allocated mappers, in order to avoid creation/destruction of
  // their underlying text properties (i.e. each time a mapper is
  // created, text properties are created and shallow-assigned a font size
  // which value might be "far" from the target font size).

  // is this a vtkLookupTable or a subclass of vtkLookupTable 
  // with its scale set to log
  vtkLookupTable *LUT = vtkLookupTable::SafeDownCast( this->LookupTable );
  int isLogTable = 0;
  if ( LUT )
    {
    if ( LUT->GetScale() == VTK_SCALE_LOG10 )
      {
      isLogTable = 1; 
      }
    }

  for (i=0; i < this->NumberOfLabels; i++)
    {
    this->TextMappers[i] = vtkTextMapper::New();

    if ( isLogTable )
      {
      double lval;
      if (this->NumberOfLabels > 1)
        {
        lval = log10(range[0]) + (double)i/(this->NumberOfLabels-1) *
          (log10(range[1])-log10(range[0]));
        }
      else
        {
        lval = log10(range[0]) + 0.5*(log10(range[1])-log10(range[0]));
        }
      val = pow(10.0,lval);
      }
    else
      {
      if (this->NumberOfLabels > 1)
        {
        val = range[0] + 
          (double)i/(this->NumberOfLabels-1) * (range[1]-range[0]);
        }
      else
        {
        val = range[0] + 0.5*(range[1]-range[0]);
        }
      }

    sprintf(string, this->LabelFormat, val);
    this->TextMappers[i]->SetInput(string);

    // Shallow copy here so that the size of the label prop is not affected
    // by the automatic adjustment of its text mapper's size (i.e. its
    // mapper's text property is identical except for the font size
    // which will be modified later). This allows text actors to
    // share the same text property, and in that case specifically allows
    // the title and label text prop to be the same.
    this->TextMappers[i]->GetTextProperty()->ShallowCopy(
      this->LabelTextProperty);

    this->TextActors[i] = vtkActor2D::New();
    this->TextActors[i]->SetMapper(this->TextMappers[i]);
    this->TextActors[i]->SetProperty(this->GetProperty());
    this->TextActors[i]->GetPositionCoordinate()->
      SetReferenceCoordinate(this->PositionCoordinate);
    }

  if (this->NumberOfLabels)
    {
    int targetWidth, targetHeight;
    // rnv begin
    // Customization of the vtkScalarBarActor to show distribution histogram.
    bool distrVisibility = ( this->MaximumNumberOfColors == (int) this->myNbValues.size() );
    double coef;
    if( GetDistributionVisibility() && distrVisibility )
      if(this->Orientation == VTK_ORIENT_VERTICAL)
        coef = 0.4;
      else 
        coef = 0.18;
    else 
      if(this->Orientation == VTK_ORIENT_VERTICAL)
        coef = 0.6;
      else 
        coef=0.25;


    if ( this->Orientation == VTK_ORIENT_VERTICAL )
      {
      targetWidth = (int)(coef*size[0]);
      targetHeight = (int)(0.86*size[1]/this->NumberOfLabels);
      }
    else
      {
      targetWidth = (int)(size[0]*0.8/this->NumberOfLabels);
      targetHeight = (int)(coef*size[1]);
      }
    // rnv end
    
    vtkTextMapper::SetMultipleConstrainedFontSize(viewport, 
                                                  targetWidth, 
                                                  targetHeight,
                                                  this->TextMappers,
                                                  this->NumberOfLabels,
                                                  labelSize);
    }
}

//----------------------------------------------------------------------------
void SMESH_ScalarBarActor::SizeTitle(int *titleSize,
                                     int *size,
                                     vtkViewport *viewport)
{
  titleSize[0] = titleSize[1] = 0;

  if (this->Title == NULL || !strlen(this->Title))
  {
    return;
  }

  int targetWidth, targetHeight;

  targetWidth = size[0];
  // rnv begin
  // Customization of the vtkScalarBarActor to show distribution histogram.
  bool distrVisibility = ( this->MaximumNumberOfColors == (int) this->myNbValues.size() );
  double coef;
  if ( GetDistributionVisibility() && distrVisibility )
    coef=0.18;
  else
    coef=0.25;

  if ( this->Orientation == VTK_ORIENT_VERTICAL )
  {
    targetHeight = (int)(0.1*size[1]);
  }
  else
  {
    targetHeight = (int)(coef*size[1]);
  }

  this->TitleMapper->SetConstrainedFontSize(viewport, targetWidth, targetHeight);

  this->TitleMapper->GetSize(viewport, titleSize);
}


/*--------------------------------------------------------------------------*/
void SMESH_ScalarBarActor::SetDistributionVisibility(int flag) {
  myDistributionActor->SetVisibility(flag);
  Modified();
}


/*--------------------------------------------------------------------------*/
int SMESH_ScalarBarActor::GetDistributionVisibility() {
  return myDistributionActor->GetVisibility();
}


void SMESH_ScalarBarActor::SetDistribution(std::vector<int> theNbValues) {
  myNbValues = theNbValues;
} 


void SMESH_ScalarBarActor::SetDistributionColor (double rgb[3]) {
  myDistributionActor->GetProperty()->SetColor(rgb);
  Modified();
}

void SMESH_ScalarBarActor::GetDistributionColor (double rgb[3]) {
  myDistributionActor->GetProperty()->GetColor(rgb);
}

void SMESH_ScalarBarActor::SetTitleOnlyVisibility( bool theTitleOnlyVisibility) {
  myTitleOnlyVisibility = theTitleOnlyVisibility;
}

bool SMESH_ScalarBarActor::GetTitleOnlyVisibility() {
  return myTitleOnlyVisibility;
}
