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

#include "SMESH_PreviewActorsCollection.h"

#include <utilities.h>

// OCC includes
#include <TopoDS.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>

// VTK includes
#include <vtkUnstructuredGrid.h>
#include <vtkPlane.h>
#include <vtkRenderer.h>
#include <vtkProperty.h>

#include <QSet>

// GEOM includes
#include <GEOM_Actor.h>

// GUI includes
#include <VTKViewer_Actor.h>
#include <SVTK_DeviceActor.h>
#include <SALOME_InteractiveObject.hxx>
#include <SUIT_ResourceMgr.h>
#include <SUIT_Session.h>

#ifdef _DEBUG_
static int MYDEBUG = 0;
#else
static int MYDEBUG = 0;
#endif


SMESH_PreviewActorsCollection::SMESH_PreviewActorsCollection() :
  mySelector( 0 ), myRenderer( 0 ), myCurrentChunk( 0 ), myChunkSize( 0 ), myIsShown( true )
{
  if(MYDEBUG) MESSAGE("SMESH_PreviewActorsCollection - "<<this);
}


SMESH_PreviewActorsCollection::~SMESH_PreviewActorsCollection()
{
  if(MYDEBUG) MESSAGE("~SMESH_PreviewActorsCollection - "<<this);
  clearActors();
}

bool SMESH_PreviewActorsCollection::Init( const TopoDS_Shape& theShape,
                                          const TopoDS_Shape& theMainShape,
                                          TopAbs_ShapeEnum    theType,
                                          const QString&      theEntry )
{
  SUIT_ResourceMgr* mgr = SUIT_Session::session()->resourceMgr();

  myType =  theType;
  myEntry = theEntry;
  myMainShape = theShape;
  myMapOfActors.clear();
  myMapOfShapes.Clear();
  myIndices.clear();
  myCurrentChunk = 0;
  myChunkSize = qMax(1, mgr->integerValue( "SMESH", "preview_actor_chunk_size", 100 ) );

  if ( theShape.IsNull() )
    return false;

  // Handle( SALOME_InteractiveObject ) anIO = new SALOME_InteractiveObject();
  // anIO->setEntry( theEntry.toLatin1().constData() );
  
  // get indexes of seleted elements
  TopExp::MapShapes( theMainShape, myMapOfShapes );
  TopExp_Explorer exp( theShape, theType );
  QSet<int> indices;
  for ( ; exp.More(); exp.Next() )
    indices << myMapOfShapes.FindIndex( exp.Current() );
  myIndices = indices.toList();
  //qSort(myIndices);

  // show current chunk
  showCurrentChunk();

  return count() > 0;
}

GEOM_Actor* SMESH_PreviewActorsCollection::createActor(const TopoDS_Shape& shape)
{
  GEOM_Actor* actor = GEOM_Actor::New();
  actor->SetShape( shape, 0, 0 );
  return actor;
}

GEOM_Actor* SMESH_PreviewActorsCollection::GetActorByIndex(int index)
{
  return myMapOfActors.value( index );
}

bool SMESH_PreviewActorsCollection::IsValidIndex( int index )
{
  return 0 < index && index <= myMapOfShapes.Extent();
}

int SMESH_PreviewActorsCollection::GetIndexByShape( const TopoDS_Shape& theShape )
{
  return myMapOfShapes.FindIndex( theShape );
}

TopoDS_Shape SMESH_PreviewActorsCollection::GetShapeByIndex( int index )
{
  return IsValidIndex( index ) ? myMapOfShapes.FindKey( index ) : TopoDS_Shape();
}

int SMESH_PreviewActorsCollection::NbShapesOfType( TopAbs_ShapeEnum type )
{
  if ( type == TopAbs_SHAPE ) return myMapOfShapes.Extent();

  int nb = 0;
  for ( int i = 1; i <= myMapOfShapes.Extent(); ++i )
    nb += ( myMapOfShapes(i).ShapeType() == type );

  return nb;
}

void SMESH_PreviewActorsCollection::SetIndices( const QList<int>& indices)
{
  if ( myIndices != indices )
  {
    myIndices = indices;
    showCurrentChunk();
  }
}

void SMESH_PreviewActorsCollection::AddToRender(vtkRenderer* theRenderer)
{
  myRenderer = theRenderer;

  QMap<int, GEOM_Actor*>::iterator iter = myMapOfActors.begin();
  for ( ; iter != myMapOfActors.end(); ++iter ) {
    iter.value()->SetVisibility( myIsShown );
    iter.value()->AddToRender( theRenderer );
  }
}

void SMESH_PreviewActorsCollection::RemoveFromRender(vtkRenderer* theRenderer)
{
  QMap<int, GEOM_Actor*>::iterator iter = myMapOfActors.begin();
  for ( ; iter != myMapOfActors.end(); ++iter )
    iter.value()->RemoveFromRender( theRenderer );
}

void SMESH_PreviewActorsCollection::SetSelector(SVTK_Selector* theSelector)
{
  mySelector = theSelector;
}

void SMESH_PreviewActorsCollection::HighlightAll( bool theHighlight )
{
  QMap<int, GEOM_Actor*>::iterator iter = myMapOfActors.begin();
  for ( ; iter != myMapOfActors.end(); ++iter )
    iter.value()->Highlight( theHighlight );
}

void SMESH_PreviewActorsCollection::HighlightID( int index )
{
  GEOM_Actor* anActor = GetActorByIndex( index );
  if ( anActor && !anActor->isHighlighted() )
    anActor->Highlight( true );
}

void SMESH_PreviewActorsCollection::SetShown( bool shown )
{
  myIsShown = shown;
  QMap<int, GEOM_Actor*>::iterator iter = myMapOfActors.begin();
  for ( ; iter != myMapOfActors.end(); ++iter )
    iter.value()->SetVisibility( shown );
}

int SMESH_PreviewActorsCollection::count() const
{
  return myIndices.count();
}

int SMESH_PreviewActorsCollection::chunkSize() const
{
  return myChunkSize;
}

int SMESH_PreviewActorsCollection::currentChunk() const
{
  return myCurrentChunk;
}

bool SMESH_PreviewActorsCollection::hasPrevious() const
{
  return chunkSize() > 0 && currentChunk() > 0;
}

bool SMESH_PreviewActorsCollection::hasNext() const
{
  return chunkSize() > 0 && (currentChunk()+1)*chunkSize() < count();
}

void SMESH_PreviewActorsCollection::previous()
{
  if ( !hasPrevious() ) return;
  myCurrentChunk--;
  showCurrentChunk();
}

void SMESH_PreviewActorsCollection::next()
{
  if ( !hasNext() ) return;
  myCurrentChunk++;
  showCurrentChunk();
}

void SMESH_PreviewActorsCollection::showCurrentChunk()
{
  clearActors();
  int imin = currentChunk() * chunkSize();
  int imax = std::min( (currentChunk()+1) * chunkSize(), count() );
  for ( int i = imin; i < imax; i++ ) {
    int index = myIndices[i];
    if ( !index || myMapOfActors.contains( index ) ) continue;
    TopoDS_Shape s = myMapOfShapes.FindKey( index );
    if ( s.IsNull() ) continue;
    // create actor if the index is present
    if ( GEOM_Actor* anActor = createActor( s.Oriented(TopAbs_FORWARD))) {
      // Create new entry for actor
      QString entry = QString( "%1_%2" ).arg( myEntry ).arg( index );
      // Create interactive object
      Handle( SALOME_InteractiveObject ) anIO = new SALOME_InteractiveObject();
      anIO->setEntry( entry.toLatin1().constData() );
      // Init Actor
      anActor->SetVectorMode( myType==TopAbs_EDGE );
      anActor->setIO( anIO );
      anActor->SetSelector( mySelector );
      anActor->SetPickable( true );
      anActor->SetResolveCoincidentTopology( true );

      // Add Actor to the Actors Map
      myMapOfActors.insert(index, anActor);
    }
  }
  if ( mySelector )
    mySelector->ClearIObjects();
  if ( myRenderer )
    AddToRender( myRenderer );
}

void SMESH_PreviewActorsCollection::clearActors()
{
  if (myRenderer)
    RemoveFromRender(myRenderer);

  QMap<int, GEOM_Actor*>::iterator iter = myMapOfActors.begin();
  for ( ; iter != myMapOfActors.end(); ++iter )
    if ( GEOM_Actor* anActor = iter.value() )
      anActor->Delete();
  myMapOfActors.clear();
}
