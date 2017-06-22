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
// File      : DriverCGNS_Read.cxx
// Created   : Thu Jun 30 10:33:31 2011
// Author    : Edward AGAPOV (eap)

#include "DriverCGNS_Read.hxx"

#include "SMDS_MeshNode.hxx"
#include "SMESHDS_Group.hxx"
#include "SMESHDS_Mesh.hxx"
#include "SMESH_Comment.hxx"
#include "SMESH_TypeDefs.hxx"

#include <gp_XYZ.hxx>

#include <cgnslib.h>

#include <map>

#if CGNS_VERSION < 3100
# define cgsize_t int
#endif

#define NB_ZONE_SIZE_VAL 9
#define CGNS_NAME_SIZE 33
#define CGNS_STRUCT_RANGE_SZ 6

using namespace std;

namespace
{
  //================================================================================
  /*!
   * \brief Data of a zone
   */
  struct TZoneData
  {
    int                    _id;
    int                    _nodeIdShift; // nb nodes in previously read zones
    int                    _elemIdShift; // nb faces in previously read zones
    int                    _nbNodes, _nbElems;
    int                    _meshDim;
    int                    _sizeX, _sizeY, _sizeZ, _nbCells; // structured
    cgsize_t               _sizes[NB_ZONE_SIZE_VAL];
    CGNS_ENUMT(ZoneType_t) _type;
    map< int, int >        _nodeReplacementMap;/* key:   id of node to replace (in this zone),
                                                  value: id of node to replace by (in another zone)
                                                  id values include _nodeIdShift of the zones */
    void SetSizeAndDim( cgsize_t* sizes, int meshDim )
    {
      _meshDim = meshDim;
      memcpy( _sizes, sizes, NB_ZONE_SIZE_VAL*sizeof(cgsize_t));
      _sizeX = _sizes[0];
      _sizeY = _meshDim > 1 ? _sizes[1] : 0;
      _sizeZ = _meshDim > 2 ? _sizes[2] : 0;
      _nbCells = (_sizeX - 1) * ( _meshDim > 1 ? _sizeY : 1 ) * ( _meshDim > 2 ? _sizeZ : 1 );
    }
    bool IsStructured() const { return ( _type == CGNS_ENUMV( Structured )); }
    int IndexSize() const { return IsStructured() ? _meshDim : 1; }
    string ReadZonesConnection(int file, int base,
                               const map< string, TZoneData >& zonesByName,
                               SMESHDS_Mesh*                   mesh);
    void ReplaceNodes( cgsize_t* ids, int nbIds, int idShift = 0 ) const;

    // Methods for a structured zone

    int NodeID( int i, int j, int k = 1 ) const
    {
      return _nodeIdShift + (k-1)*_sizeX*_sizeY + (j-1)*_sizeX + i;
    }
    int NodeID( const gp_XYZ& ijk ) const
    {
      return NodeID( int(ijk.X()), int(ijk.Y()), int(ijk.Z()));
    }
    void CellNodes( int i, int j, int k, cgsize_t* ids ) const
    {
      ids[0] = NodeID( i  , j  , k  );
      ids[1] = NodeID( i  , j+1, k  );
      ids[2] = NodeID( i+1, j+1, k  );
      ids[3] = NodeID( i+1, j  , k  );
      ids[4] = NodeID( i  , j  , k+1);
      ids[5] = NodeID( i  , j+1, k+1);
      ids[6] = NodeID( i+1, j+1, k+1);
      ids[7] = NodeID( i+1, j  , k+1);
    }
    void CellNodes( int i, int j, cgsize_t* ids ) const
    {
      ids[0] = NodeID( i  , j   );
      ids[1] = NodeID( i  , j+1 );
      ids[2] = NodeID( i+1, j+1 );
      ids[3] = NodeID( i+1, j   );
    }
    void IFaceNodes( int i, int j, int k, cgsize_t* ids ) const // face perpendiculaire to X (3D)
    {
      ids[0] = NodeID( i, j, k );
      ids[1] = ids[0] + _sizeX*( i==_sizeX ? 1 : _sizeY );
      ids[2] = ids[0] + _sizeX*( _sizeY + 1 );
      ids[3] = ids[0] + _sizeX*( i==_sizeX ? _sizeY : 1 );
    }
    void JFaceNodes( int i, int j, int k, cgsize_t* ids ) const
    {
      ids[0] = NodeID( i, j, k );
      ids[1] = ids[0] + ( j==_sizeY ? _sizeX*_sizeY : 1);
      ids[2] = ids[0] + _sizeX*_sizeY + 1;
      ids[3] = ids[0] + ( j==_sizeY ? 1 : _sizeX*_sizeY);
    }
    void KFaceNodes( int i, int j, int k, cgsize_t* ids ) const
    {
      ids[0] = NodeID( i, j, k );
      ids[1] = ids[0] + ( k==_sizeZ ? 1 : _sizeX);
      ids[2] = ids[0] + _sizeX + 1;
      ids[3] = ids[0] + ( k==_sizeZ ? _sizeX : 1);
    }
    void IEdgeNodes( int i, int j, int k, cgsize_t* ids ) const // edge perpendiculaire to X (2D)
    {
      ids[0] = NodeID( i, j, 0 );
      ids[1] = ids[0] + _sizeX;
    }
    void JEdgeNodes( int i, int j, int k, cgsize_t* ids ) const
    {
      ids[0] = NodeID( i, j, 0 );
      ids[1] = ids[0] + 1;
    }
#define gpXYZ2IJK(METHOD)                                     \
    void METHOD( const gp_XYZ& ijk, cgsize_t* ids ) const {        \
      METHOD( int(ijk.X()), int(ijk.Y()), int(ijk.Z()), ids); \
    }
    gpXYZ2IJK( IFaceNodes )
    gpXYZ2IJK( JFaceNodes )
    gpXYZ2IJK( KFaceNodes )
    gpXYZ2IJK( IEdgeNodes )
    gpXYZ2IJK( JEdgeNodes )
  };

  //================================================================================
  /*!
   * \brief Iterator over nodes of the structired grid using FORTRAN multidimensional
   * array ordering.
   */
  class TPointRangeIterator
  {
    int _beg[3], _end[3], _cur[3], _dir[3], _dim;
    bool _more;
  public:
    TPointRangeIterator( const cgsize_t* range, int dim ):_dim(dim)
    {
      _more = false;
      for ( int i = 0; i < dim; ++i )
      {
        _beg[i] = range[i];
        _end[i] = range[i+dim];
        _dir[i] = _end[i] < _beg[i] ? -1 : 1;
        _end[i] += _dir[i];
        _cur[i] = _beg[i];
        if ( _end[i] - _beg[i] )
          _more = true;
      }
//       for ( int i = dim; i < 3; ++i )
//         _cur[i] = _beg[i] = _end[i] = _dir[i] = 0;
    }
    bool More() const
    {
      return _more;
    }
    gp_XYZ Next()
    {
      gp_XYZ res( _cur[0], _cur[1], _cur[2] );
      for ( int i = 0; i < _dim; ++i )
      {
        _cur[i] += _dir[i];
        if ( _cur[i]*_dir[i] < _end[i]*_dir[i] )
          break;
        if ( i+1 < _dim )
          _cur[i] = _beg[i];
        else
          _more = false;
      }
      return res;
    }
    size_t Size()  const
    {
      size_t size = 1;
      for ( int i = 0; i < _dim; ++i )
        size *= _dir[i]*(_end[i]-_beg[i]);
      return size;
    }
    gp_XYZ Begin() const { return gp_XYZ( _beg[0], _beg[1], _beg[2] ); }
    //gp_XYZ End() const { return gp_XYZ( _end[0]-1, _end[1]-1, _end[2]-1 ); }
  };

  //================================================================================
  /*!
   * \brief Checks if the two arrays of node IDs describe nodes with equal coordinates
   */
  //================================================================================

  bool isEqualNodes( const int* nIds1, const int* nIds2, int nbNodes, SMESHDS_Mesh* mesh )
  {
    if ( nbNodes > 0 )
    {
      SMESH_TNodeXYZ nn1[2], nn2[2];
      nn1[0] = mesh->FindNode( nIds1[0] );
      nn2[0] = mesh->FindNode( nIds2[0] );
      if ( !nn1[0]._node || !nn2[0]._node )
        return false;
      double dist1 = ( nn1[0] - nn2[0] ).Modulus();
      double dist2 = 0, tol = 1e-7;
      if ( nbNodes > 1 )
      {
        nn1[1] = mesh->FindNode( nIds1[1] );
        nn2[1] = mesh->FindNode( nIds2[1] );
        if ( !nn1[1]._node || !nn2[1]._node )
          return false;
        dist2 = ( nn1[1] - nn2[1] ).Modulus();
        tol   = 1e-5 * ( nn1[0] - nn1[1] ).Modulus();
      }
      return ( dist1 < tol && dist2 < tol );
    }
    return false;
  }

  //================================================================================
  /*!
   * \brief Reads zone interface connectivity
   *  \param file - file to read
   *  \param base - base to read
   *  \param zone - zone to replace nodes in
   *  \param zonesByName - TZoneData by name
   *  \retval string - warning message
   *
   * see // http://www.grc.nasa.gov/WWW/cgns/CGNS_docs_current/sids/cnct.html
   */
  //================================================================================

  string TZoneData::ReadZonesConnection( int                             file,
                                         int                             base,
                                         const map< string, TZoneData >& zonesByName,
                                         SMESHDS_Mesh*                   mesh)
  {
    string error;

    char connectName[ CGNS_NAME_SIZE ], donorName [ CGNS_NAME_SIZE ];

    // ----------------------------
    // read zone 1 to 1 interfaces
    // ----------------------------
    if ( IsStructured() )
    {
      int nb1to1 = 0;
      if ( cg_n1to1 ( file, base, _id, &nb1to1) == CG_OK )
      {
        cgsize_t range[CGNS_STRUCT_RANGE_SZ], donorRange[CGNS_STRUCT_RANGE_SZ];
        int transform[3] = {0,0,0};

        for ( int I = 1; I <= nb1to1; ++I )
        {
          if ( cg_1to1_read(file, base, _id, I, connectName,
                            donorName, range, donorRange, transform) == CG_OK )
          {
            map< string, TZoneData>::const_iterator n_z = zonesByName.find( donorName );
            if ( n_z == zonesByName.end() )
              continue; // donor zone not yet read
            const TZoneData& zone2 = n_z->second;

            // set up matrix to transform ijk of the zone to ijk of the zone2
            gp_Mat T;
            for ( int i = 0; i < _meshDim; ++i )
              if ( transform[i] )
              {
                int row = Abs(transform[i]);
                int col = i+1;
                int val = transform[i] > 0 ? +1 : -1;
                T( row, col ) = val;
              }

            // fill nodeReplacementMap
            TPointRangeIterator rangeIt1( range, _meshDim );
            TPointRangeIterator rangeIt2( donorRange, _meshDim );
            gp_XYZ begin1 = rangeIt1.Begin(), begin2 = rangeIt2.Begin(), index1, index2;
            if ( &zone2 == this )
            {
              // not to read twice the same interface with self
              TPointRangeIterator rangeIt1bis( range, _meshDim );
              if ( rangeIt1bis.More() )
              {
                index1 = rangeIt1bis.Next();
                index2 = T * ( index1 - begin1 ) + begin2;
                int node1 = NodeID( index1 );
                int node2 = zone2.NodeID( index2 );
                if ( _nodeReplacementMap.count( node2 ) &&
                     _nodeReplacementMap[ node2 ] == node1 )
                  continue; // this interface already read
              }
            }
            // check if range and donorRange describe the same nodes
            {
              int ids1[2], ids2[2], nbN = 0;
              TPointRangeIterator rangeIt1bis( range, _meshDim );
              index1 = rangeIt1bis.Next();
              index2 = T * ( index1 - begin1 ) + begin2;
              ids1[0] = NodeID( index1 );
              ids2[0] = zone2.NodeID( index2 );
              ++nbN;
              if ( rangeIt1bis.More() )
              {
                index1 = rangeIt1bis.Next();
                index2 = T * ( index1 - begin1 ) + begin2;
                ids1[1] = NodeID( index1 );
                ids2[1] = zone2.NodeID( index2 );
                ++nbN;
              }
              if ( !isEqualNodes( &ids1[0], &ids2[0], nbN, mesh ))
                continue;
            }
            while ( rangeIt1.More() )
            {
              index1 = rangeIt1.Next();
              index2 = T * ( index1 - begin1 ) + begin2;
              int node1 = NodeID( index1 );
              int node2 = zone2.NodeID( index2 );
              _nodeReplacementMap.insert( make_pair( node1, node2 ));
            }
          }
          else
          {
            error = cg_get_error();
          }
        }
      }
      else
      {
        error = cg_get_error();
      }
    }

    // ---------------------------------
    // read general zone connectivities
    // ---------------------------------
    int nbConn = 0;
    if ( cg_nconns( file, base, _id, &nbConn) == CG_OK )
    {
      cgsize_t nb, donorNb;
      CGNS_ENUMT(GridLocation_t) location;
      CGNS_ENUMT(GridConnectivityType_t) connectType;
      CGNS_ENUMT(PointSetType_t) ptype, donorPtype;
      CGNS_ENUMT(ZoneType_t) donorZonetype;
      CGNS_ENUMT(DataType_t) donorDatatype;

      for ( int I = 1; I <= nbConn; ++I )
      {
        if ( cg_conn_info(file, base, _id, I, connectName, &location, &connectType,
                          &ptype, &nb, donorName, &donorZonetype, &donorPtype,
                          &donorDatatype, &donorNb ) == CG_OK )
        {
          if ( location != CGNS_ENUMV( Vertex ))
            continue; // we do not support cell-to-cell connectivity
          if ( ptype != CGNS_ENUMV( PointList ) &&
               ptype != CGNS_ENUMV( PointRange ))
            continue;
          if ( donorPtype != CGNS_ENUMV( PointList ) &&
               donorPtype != CGNS_ENUMV( PointRange ))
            continue;
          
          map< string, TZoneData>::const_iterator n_z = zonesByName.find( donorName );
          if ( n_z == zonesByName.end() )
            continue; // donor zone not yet read
          const TZoneData& zone2 = n_z->second;

          vector< cgsize_t > ids( nb * IndexSize() );
          vector< cgsize_t > donorIds( donorNb * zone2.IndexSize() );
          if (cg_conn_read ( file, base, _id, I,
                             &ids[0], CGNS_ENUMV(Integer), &donorIds[0]) == CG_OK )
          {
            for ( int isThisZone = 0; isThisZone < 2; ++isThisZone )
            {
              const TZoneData&           zone = isThisZone ? *this : zone2;
              CGNS_ENUMT(PointSetType_t) type = isThisZone ? ptype : donorPtype;
              vector< cgsize_t >&      points = isThisZone ? ids : donorIds;
              if ( type == CGNS_ENUMV( PointRange ))
              {
                TPointRangeIterator rangeIt( &points[0], zone._meshDim );
                points.clear();
                while ( rangeIt.More() )
                  points.push_back ( NodeID( rangeIt.Next() ));
              }
              else if ( zone.IsStructured() )
              {
                vector< cgsize_t > resIDs; resIDs.reserve( points.size() / IndexSize() );
                for ( size_t i = 0; i < points.size(); i += IndexSize() )
                  resIDs.push_back( zone.NodeID( points[i+0], points[i+1], points[i+2] ));
                resIDs.swap( points );
              }
              else if ( zone._nodeIdShift > 0 )
              {
                for ( size_t i = 0; i < points.size(); ++i )
                  points[i] += zone._nodeIdShift;
              }
            }
            size_t nbN = std::min( ids.size(), donorIds.size());
            if ( isEqualNodes( &ids[0], &donorIds[0], nbN, mesh ))
              for ( size_t i = 0; i < nbN; ++i )
                _nodeReplacementMap.insert( make_pair( ids[i], donorIds[i] ));
          }
          else
          {
            error = cg_get_error();
          }
        }
        else
        {
          error = cg_get_error();
        }
      }
    }
    else
    {
      error = cg_get_error();
    }
    return error;
  }

  //================================================================================
  /*!
   * \brief Replaces node ids according to nodeReplacementMap to take into account
   *        connection of zones
   */
  //================================================================================

  void TZoneData::ReplaceNodes( cgsize_t* ids, int nbIds, int idShift/* = 0*/ ) const
  {
    if ( !_nodeReplacementMap.empty() )
    {
      map< int, int >::const_iterator it, end = _nodeReplacementMap.end();
      for ( int i = 0; i < nbIds; ++i )
        if (( it = _nodeReplacementMap.find( ids[i] + idShift)) != end )
          ids[i] = it->second;
        else
          ids[i] += idShift;
    }
    else if ( idShift )
    {
      for ( int i = 0; i < nbIds; ++i )
        ids[i] += idShift;
    }
  }
  //================================================================================
  /*!
   * \brief functions adding an element of a particular type
   */
  SMDS_MeshElement* add_0D(cgsize_t* ids, SMESHDS_Mesh* mesh, int ID)
  {
    return mesh->Add0DElementWithID( ids[0], ID );
  }
  SMDS_MeshElement* add_BAR_2(cgsize_t* ids, SMESHDS_Mesh* mesh, int ID)
  {
    return mesh->AddEdgeWithID( ids[0], ids[1], ID );
  }
  SMDS_MeshElement* add_BAR_3(cgsize_t* ids, SMESHDS_Mesh* mesh, int ID)
  {
    return mesh->AddEdgeWithID( ids[0], ids[1], ids[2], ID );
  }
  SMDS_MeshElement* add_TRI_3(cgsize_t* ids, SMESHDS_Mesh* mesh, int ID)
  {
    return mesh->AddFaceWithID( ids[0], ids[2], ids[1], ID );
  }
  SMDS_MeshElement* add_TRI_6(cgsize_t* ids, SMESHDS_Mesh* mesh, int ID)
  {
    return mesh->AddFaceWithID( ids[0], ids[2], ids[1], ids[5], ids[4], ids[3], ID );
  }
  SMDS_MeshElement* add_QUAD_4(cgsize_t* ids, SMESHDS_Mesh* mesh, int ID)
  {
    return mesh->AddFaceWithID( ids[0], ids[3], ids[2], ids[1], ID );
  }
  SMDS_MeshElement* add_QUAD_8(cgsize_t* ids, SMESHDS_Mesh* mesh, int ID)
  {
    return mesh->AddFaceWithID( ids[0],ids[3],ids[2],ids[1],ids[7],ids[6],ids[5],ids[4], ID );
  }
  SMDS_MeshElement* add_QUAD_9(cgsize_t* ids, SMESHDS_Mesh* mesh, int ID)
  {
    return mesh->AddFaceWithID( ids[0],ids[3],ids[2],ids[1],ids[7],ids[6],ids[5],ids[4],ids[8], ID);
  }
  SMDS_MeshElement* add_TETRA_4(cgsize_t* ids, SMESHDS_Mesh* mesh, int ID)
  {
    return mesh->AddVolumeWithID( ids[0], ids[2], ids[1], ids[3], ID );
  }
  SMDS_MeshElement* add_TETRA_10(cgsize_t* ids, SMESHDS_Mesh* mesh, int ID)
  {
    return mesh->AddVolumeWithID( ids[0],ids[2],ids[1],ids[3],ids[6],
                                  ids[5],ids[4],ids[7],ids[9],ids[8], ID );
  }
  SMDS_MeshElement* add_PYRA_5(cgsize_t* ids, SMESHDS_Mesh* mesh, int ID)
  {
    return mesh->AddVolumeWithID( ids[0],ids[3],ids[2],ids[1],ids[4],ID );
  }
  SMDS_MeshElement* add_PYRA_13(cgsize_t* ids, SMESHDS_Mesh* mesh, int ID)
  {
    return mesh->AddVolumeWithID( ids[0],ids[3],ids[2],ids[1],ids[4],ids[8],ids[7],
                                  ids[6],ids[5],ids[9],ids[12],ids[11],ids[10], ID );
  }
  SMDS_MeshElement* add_PENTA_6(cgsize_t* ids, SMESHDS_Mesh* mesh, int ID)
  {
    return mesh->AddVolumeWithID( ids[0],ids[2],ids[1],ids[3],ids[5],ids[4], ID );
  }
  SMDS_MeshElement* add_PENTA_15(cgsize_t* ids, SMESHDS_Mesh* mesh, int ID)
  {
    return mesh->AddVolumeWithID( ids[0],ids[2],ids[1],ids[3],ids[5],ids[4],ids[8],ids[7],
                                  ids[6],ids[9],ids[11],ids[10],ids[14],ids[13],ids[12], ID );
  }
  SMDS_MeshElement* add_HEXA_8(cgsize_t* ids, SMESHDS_Mesh* mesh, int ID)
  {
    return mesh->AddVolumeWithID( ids[0],ids[3],ids[2],ids[1],ids[4],ids[7],ids[6],ids[5], ID );
  }
  SMDS_MeshElement* add_HEXA_20(cgsize_t* ids, SMESHDS_Mesh* mesh, int ID)
  {
    return mesh->AddVolumeWithID( ids[0],ids[3],ids[2],ids[1],ids[4],ids[7],ids[6],
                                  ids[5],ids[11],ids[10],ids[9],ids[8],ids[12],ids[15],
                                  ids[14],ids[13],ids[19],ids[18],ids[17],ids[16], ID );
  }
  SMDS_MeshElement* add_HEXA_27(cgsize_t* ids, SMESHDS_Mesh* mesh, int ID)
  {
    return mesh->AddVolumeWithID( ids[0],ids[3],ids[2],ids[1],ids[4],ids[7],ids[6],
                                  ids[5],ids[11],ids[10],ids[9],ids[8],ids[12],ids[15],
                                  ids[14],ids[13],ids[19],ids[18],ids[17],ids[16],
                                  ids[20],ids[24],ids[23],ids[22],ids[21],ids[25],ids[26], ID );
  }
  SMDS_MeshElement* add_NGON(cgsize_t* ids, SMESHDS_Mesh* mesh, int ID)
  {
    vector<int> idVec( ids[0] );
    for ( int i = 0; i < ids[0]; ++i )
      idVec[ i ] = (int) ids[ i + 1];
    return mesh->AddPolygonalFaceWithID( idVec, ID );
  }

  typedef SMDS_MeshElement* (* PAddElemFun) (cgsize_t* ids, SMESHDS_Mesh* mesh, int ID);
  
  //================================================================================
  /*!
   * \brief Return an array of functions each adding an element of a particular type
   */
  //================================================================================

  PAddElemFun* getAddElemFunTable()
  {
    static vector< PAddElemFun > funVec;
    if ( funVec.empty() )
    {
      funVec.resize( NofValidElementTypes, (PAddElemFun)0 );
      funVec[ CGNS_ENUMV( NODE     )] = add_0D      ;
      funVec[ CGNS_ENUMV( BAR_2    )] = add_BAR_2   ;
      funVec[ CGNS_ENUMV( BAR_3    )] = add_BAR_3   ;
      funVec[ CGNS_ENUMV( TRI_3    )] = add_TRI_3   ;
      funVec[ CGNS_ENUMV( TRI_6    )] = add_TRI_6   ;
      funVec[ CGNS_ENUMV( QUAD_4   )] = add_QUAD_4  ;
      funVec[ CGNS_ENUMV( QUAD_8   )] = add_QUAD_8  ;
      funVec[ CGNS_ENUMV( QUAD_9   )] = add_QUAD_9  ;
      funVec[ CGNS_ENUMV( TETRA_4  )] = add_TETRA_4 ;
      funVec[ CGNS_ENUMV( TETRA_10 )] = add_TETRA_10;
      funVec[ CGNS_ENUMV( PYRA_5   )] = add_PYRA_5  ;
      funVec[ CGNS_ENUMV( PYRA_13  )] = add_PYRA_13 ;
      funVec[ CGNS_ENUMV( PYRA_14  )] = add_PYRA_13 ;
      funVec[ CGNS_ENUMV( PENTA_6  )] = add_PENTA_6 ;
      funVec[ CGNS_ENUMV( PENTA_15 )] = add_PENTA_15;
      funVec[ CGNS_ENUMV( PENTA_18 )] = add_PENTA_15;
      funVec[ CGNS_ENUMV( HEXA_8   )] = add_HEXA_8  ;
      funVec[ CGNS_ENUMV( HEXA_20  )] = add_HEXA_20 ;
      funVec[ CGNS_ENUMV( HEXA_27  )] = add_HEXA_27 ;
      funVec[ CGNS_ENUMV( NGON_n   )] = add_NGON    ;
    }
    return &funVec[0];
  }

  //================================================================================
  /*!
   * \brief Finds an existing boundary element
   */
  //================================================================================

  const SMDS_MeshElement* findElement(const cgsize_t*     nodeIDs,
                                      const int           nbNodes,
                                      const SMESHDS_Mesh* mesh)
  {
    const SMDS_MeshNode* nn[4]; // look for quad4 or seg2
    if (( nn[0] = mesh->FindNode( nodeIDs[0] )))
    {
      SMDSAbs_ElementType eType = nbNodes==4 ? SMDSAbs_Face : SMDSAbs_Edge;
      SMDS_ElemIteratorPtr eIt = nn[0]->GetInverseElementIterator( eType );
      if ( eIt->more() )
        for ( int i = 1; i < nbNodes; ++i )
          nn[i] = mesh->FindNode( nodeIDs[i] );
      while ( eIt->more() )
      {
        const SMDS_MeshElement* e = eIt->next();
        if ( e->NbNodes() == nbNodes )
        {
          bool elemOK = true;
          for ( int i = 1; i < nbNodes && elemOK; ++i )
            elemOK = ( e->GetNodeIndex( nn[i] ) >= 0 );
          if ( elemOK )
            return e;
        }
      } 
    }
    return 0;
  }

} // namespace

//================================================================================
/*!
 * \brief Perform reading a myMeshId-th mesh
 */
//================================================================================

Driver_Mesh::Status DriverCGNS_Read::Perform()
{
  myErrorMessages.clear();

  Status aResult;
  if (( aResult = open() ) != DRS_OK )
    return aResult;

  // read nb of meshes (CGNSBase_t)
  if ( myMeshId < 0 || myMeshId >= GetNbMeshes(aResult))
    return addMessage( SMESH_Comment("Invalid mesh index :") << myMeshId );

  // read a name and a dimension of the mesh
  const int cgnsBase = myMeshId + 1;
  char meshName[CGNS_NAME_SIZE];
  int meshDim, spaceDim;
  if ( cg_base_read( _fn, cgnsBase, meshName, &meshDim, &spaceDim) != CG_OK )
    return addMessage( cg_get_error() );

  if ( spaceDim < 1 || spaceDim > 3 )
    return addMessage( SMESH_Comment("Invalid space dimension: ") << spaceDim
                       << " in mesh '" << meshName << "'");

  myMeshName = meshName;

  // read nb of domains (Zone_t) in the mesh
  int nbZones = 0;
  if ( cg_nzones (_fn, cgnsBase, &nbZones) != CG_OK )
    return addMessage( cg_get_error() );

  if ( nbZones < 1 )
    return addMessage( SMESH_Comment("Empty mesh: '") << meshName << "'");

  // read the domains (zones)
  // ------------------------
  map< string, TZoneData > zonesByName;
  char name[CGNS_NAME_SIZE];
  cgsize_t sizes[NB_ZONE_SIZE_VAL];
  memset(sizes, 0, NB_ZONE_SIZE_VAL * sizeof(cgsize_t));

  const SMDS_MeshInfo& meshInfo = myMesh->GetMeshInfo();
  int groupID = myMesh->GetGroups().size();

  for ( int iZone = 1; iZone <= nbZones; ++iZone )
  {
    // size and name of a zone
    if ( cg_zone_read( _fn, cgnsBase, iZone, name, sizes) != CG_OK) {
      addMessage( cg_get_error() );
      continue;
    }
    TZoneData& zone = zonesByName[ name ];
    zone._id          = iZone;
    zone._nodeIdShift = meshInfo.NbNodes();
    zone._elemIdShift = meshInfo.NbElements();
    zone.SetSizeAndDim( sizes, meshDim );

    // mesh type of the zone
    if ( cg_zone_type ( _fn, cgnsBase, iZone, &zone._type) != CG_OK) {
      addMessage( cg_get_error() );
      continue;
    }

    switch ( zone._type )
    {
    case CGNS_ENUMV( Unstructured ):
    case CGNS_ENUMV( Structured ):
      break;
    case CGNS_ENUMV( ZoneTypeNull ):
      addMessage( "Meshes with ZoneTypeNull are not supported");
      continue;
    case CGNS_ENUMV( ZoneTypeUserDefined ):
      addMessage( "Meshes with ZoneTypeUserDefined are not supported");
      continue;
    default:
      addMessage( "Unknown ZoneType_t");
      continue;
    }

    // -----------
    // Read nodes
    // -----------

    if ( cg_ncoords( _fn, cgnsBase, iZone, &spaceDim) != CG_OK ) {
      addMessage( cg_get_error() );
      continue;
    }
    if ( spaceDim < 1 ) {
      addMessage( SMESH_Comment("No coordinates defined in zone ")
                  << iZone << " of Mesh " << myMeshId );
      continue;
    }
    // read coordinates

    cgsize_t rmin[3] = {1,1,1}; // range of nodes to read
    cgsize_t rmax[3] = {1,1,1};
    int nbNodes = rmax[0] = zone._sizes[0];
    if ( zone.IsStructured())
      for ( int i = 1; i < meshDim; ++i )
        nbNodes *= rmax[i] = zone._sizes[i];

    vector<double> coords[3];
    for ( int c = 1; c <= spaceDim; ++c)
    {
      coords[c-1].resize( nbNodes );

      CGNS_ENUMV( DataType_t ) type;
      if ( cg_coord_info( _fn, cgnsBase, iZone, c, &type, name) != CG_OK ||
           cg_coord_read( _fn, cgnsBase, iZone, name, CGNS_ENUMV(RealDouble),
                          rmin, rmax, (void*)&(coords[c-1][0])) != CG_OK)
      {
        addMessage( cg_get_error() );
        coords[c-1].clear();
        break;
      }
    }
    if ( coords[ spaceDim-1 ].empty() )
      continue; // there was an error while reading coordinates 

    // fill coords with zero if spaceDim < 3
    for ( int c = 2; c <= 3; ++c)
      if ( coords[ c-1 ].empty() )
        coords[ c-1 ].resize( nbNodes, 0.0 );

    // create nodes
    try {
      for ( int i = 0; i < nbNodes; ++i )
        myMesh->AddNodeWithID( coords[0][i], coords[1][i], coords[2][i], i+1+zone._nodeIdShift );
    }
    catch ( std::exception& exc ) // expect std::bad_alloc
    {
      addMessage( exc.what() );
      break;
    }

    // Read connectivity between zones. Nodes of the zone interface will be
    // replaced withing the zones read later
    string err = zone.ReadZonesConnection( _fn, cgnsBase, zonesByName, myMesh );
    if ( !err.empty() )
      addMessage( err );

    // --------------
    // Read elements
    // --------------
    if ( zone.IsStructured())
    {
      int nbI = zone._sizeX - 1, nbJ = zone._sizeY - 1, nbK = zone._sizeZ - 1;
      cgsize_t nID[8];
      if ( meshDim > 2 && nbK > 0 )
      {
        for ( int k = 1; k <= nbK; ++k )
          for ( int j = 1; j <= nbJ; ++j )
            for ( int i = 1; i <= nbI; ++i )
            {
              zone.CellNodes( i, j, k, nID );
              zone.ReplaceNodes( nID, 8 );
              myMesh->AddVolumeWithID(nID[0],nID[1],nID[2],nID[3],nID[4],nID[5],nID[6],nID[7],
                                      meshInfo.NbElements()+1);
            }
      }
      else if ( meshDim > 1 && nbJ > 0 )
      {
        for ( int j = 1; j <= nbJ; ++j )
          for ( int i = 1; i <= nbI; ++i )
          {
            zone.CellNodes( i, j, nID );
            zone.ReplaceNodes( nID, 4 );
            myMesh->AddFaceWithID(nID[0],nID[1],nID[2],nID[3], meshInfo.NbElements()+1);
          }
      }
      else if ( meshDim > 0 && nbI > 0 )
      {
        nID[0] = zone.NodeID( 1, 0, 0 );
        for ( int i = 1; i <= nbI; ++i, ++nID[0] )
        {
          nID[1] = nID[0]+1;
          zone.ReplaceNodes( nID, 2 );
          myMesh->AddEdgeWithID(nID[0],nID[1], meshInfo.NbElements()+1);
        }
      }
    }
    else
    {
      // elements can be stored in different sections each dedicated to one element type
      int nbSections = 0;
      if ( cg_nsections( _fn, cgnsBase, iZone, &nbSections) != CG_OK)
      {
        addMessage( cg_get_error() );
        continue;
      }
      PAddElemFun* addElemFuns = getAddElemFunTable(), curAddElemFun = 0;
      int nbNotSuppElem = 0; // nb elements of not supported types
      bool polyhedError = false; // error at polyhedron creation

      // read element data

      CGNS_ENUMT( ElementType_t ) elemType;
      cgsize_t start, end; // range of ids of elements of a zone
      cgsize_t eDataSize = 0;
      int nbBnd, parent_flag;
      for ( int iSec = 1; iSec <= nbSections; ++iSec )
      {
        if ( cg_section_read( _fn, cgnsBase, iZone, iSec, name, &elemType,
                              &start, &end, &nbBnd, &parent_flag) != CG_OK ||
             cg_ElementDataSize( _fn, cgnsBase, iZone, iSec, &eDataSize ) != CG_OK )
        {
          addMessage( cg_get_error() );
          continue;
        }
        vector< cgsize_t > elemData( eDataSize );
        if ( cg_elements_read( _fn, cgnsBase, iZone, iSec, &elemData[0], NULL ) != CG_OK )
        {
          addMessage( cg_get_error() );
          continue;
        }
        // store elements

        int pos = 0, cgnsNbNodes = 0, elemID = start + zone._elemIdShift;
        cg_npe( elemType, &cgnsNbNodes ); // get nb nodes by element type
        curAddElemFun = addElemFuns[ elemType ];
        SMDS_MeshElement* newElem = 0;
        const SMDS_MeshElement* face;

        while ( pos < eDataSize )
        {
          CGNS_ENUMT( ElementType_t ) currentType = elemType;
          if ( currentType == CGNS_ENUMV( MIXED )) {
            //ElementConnectivity = Etype1, Node11, Node21, ... NodeN1,
            //                      Etype2, Node12, Node22, ... NodeN2,
            //                      ...
            //                      EtypeM, Node1M, Node2M, ... NodeNM
            currentType = (CGNS_ENUMT(ElementType_t)) elemData[ pos++ ];
            cg_npe( currentType, &cgnsNbNodes );
            curAddElemFun = addElemFuns[ currentType ];
          }
          if ( cgnsNbNodes < 1 ) // poly elements
          {
            if ( currentType == CGNS_ENUMV( NFACE_n )) // polyhedron
            {
              //ElementConnectivity = Nfaces1, Face11, Face21, ... FaceN1,
              //                      Nfaces2, Face12, Face22, ... FaceN2,
              //                      ...
              //                      NfacesM, Face1M, Face2M, ... FaceNM
              const int nbFaces = elemData[ pos++ ];
              vector<int> quantities( nbFaces );
              vector<const SMDS_MeshNode*> nodes, faceNodes;
              nodes.reserve( nbFaces * 4 );
              for ( int iF = 0; iF < nbFaces; ++iF )
              {
                const int faceID = std::abs( elemData[ pos++ ]) + zone._elemIdShift; 
                if (( face = myMesh->FindElement( faceID )) && face->GetType() == SMDSAbs_Face )
                {
                  const bool reverse = ( elemData[ pos-1 ] < 0 );
                  const int    iQuad = face->IsQuadratic() ? 1 : 0;
                  SMDS_ElemIteratorPtr nIter = face->interlacedNodesElemIterator();
                  faceNodes.assign( SMDS_MeshElement::iterator( nIter ),
                                    SMDS_MeshElement::iterator());
                  if ( iQuad && reverse )
                    nodes.push_back( faceNodes[0] );
                  if ( reverse )
                    nodes.insert( nodes.end(), faceNodes.rbegin(), faceNodes.rend() - iQuad );
                  else
                    nodes.insert( nodes.end(), faceNodes.begin(), faceNodes.end() );

                  quantities[ iF ] = face->NbNodes();
                }
                else {
                  polyhedError = true;
                  break;
                }
              }
              if ( quantities.back() )
              {
                myMesh->AddPolyhedralVolumeWithID( nodes, quantities, elemID );
              }
            }
            else if ( currentType == CGNS_ENUMV( NGON_n )) // polygon
            {
              // ElementConnectivity = Nnodes1, Node11, Node21, ... NodeN1,
              //                       Nnodes2, Node12, Node22, ... NodeN2,
              //                       ...
              //                       NnodesM, Node1M, Node2M, ... NodeNM
              const int nbNodes = elemData[ pos ];
              zone.ReplaceNodes( &elemData[pos+1], nbNodes, zone._nodeIdShift );
              newElem = add_NGON( &elemData[pos  ], myMesh, elemID );
              pos += nbNodes + 1;
            }
          }
          else // standard elements
          {
            zone.ReplaceNodes( &elemData[pos], cgnsNbNodes, zone._nodeIdShift );
            newElem = curAddElemFun( &elemData[pos], myMesh, elemID );
            pos += cgnsNbNodes;
            nbNotSuppElem += int( newElem && newElem->NbNodes() != cgnsNbNodes );
          }
          elemID++;

        } // loop on elemData
      } // loop on cgns sections

      if ( nbNotSuppElem > 0 )
        addMessage( SMESH_Comment(nbNotSuppElem) << " elements of not supported types"
                    << " have beem converted to close types");
      if ( polyhedError )
        addMessage( "Some polyhedral elements have been skipped due to internal(?) errors" );

    } // reading unstructured elements

    zone._nbNodes = meshInfo.NbNodes() - zone._nodeIdShift;
    zone._nbElems = meshInfo.NbElements() - zone._elemIdShift;

    // -------------------------------------------
    // Read Boundary Conditions into SMESH groups
    // -------------------------------------------
    int nbBC = 0;
    if ( cg_nbocos( _fn, cgnsBase, iZone, &nbBC) == CG_OK )
    {
      CGNS_ENUMT( BCType_t ) bcType;
      CGNS_ENUMT( PointSetType_t ) psType;
      CGNS_ENUMT( DataType_t ) normDataType;
      cgsize_t nbPnt, normFlag;
      int normIndex[3], nbDS;
      for ( int iBC = 1; iBC <= nbBC; ++iBC )
      {
        if ( cg_boco_info( _fn, cgnsBase, iZone, iBC, name, &bcType, &psType,
                           &nbPnt, normIndex, &normFlag, &normDataType, &nbDS ) != CG_OK )
        {
          addMessage( cg_get_error() );
          continue;
        }
        vector< cgsize_t > ids( nbPnt * zone.IndexSize() );
        CGNS_ENUMT( GridLocation_t ) location;
        if ( cg_boco_read( _fn, cgnsBase, iZone, iBC, &ids[0], NULL ) != CG_OK ||
             cg_boco_gridlocation_read( _fn, cgnsBase, iZone, iBC, &location) != CG_OK )
        {
          addMessage( cg_get_error() );
          continue;
        }
        SMDSAbs_ElementType elemType = SMDSAbs_All;
        switch ( location ) {
        case CGNS_ENUMV( Vertex      ): elemType = SMDSAbs_Node; break;
        case CGNS_ENUMV( FaceCenter  ): elemType = SMDSAbs_Face; break;
        case CGNS_ENUMV( IFaceCenter ): elemType = SMDSAbs_Face; break;
        case CGNS_ENUMV( JFaceCenter ): elemType = SMDSAbs_Face; break;
        case CGNS_ENUMV( KFaceCenter ): elemType = SMDSAbs_Face; break;
        case CGNS_ENUMV( EdgeCenter  ): elemType = SMDSAbs_Edge; break;
        default:;
        }
        SMESHDS_Group* group = new SMESHDS_Group ( groupID++, myMesh, elemType );
        myMesh->AddGroup( group );
        SMESH_Comment groupName( name ); groupName << " " << cg_BCTypeName( bcType );
        group->SetStoreName( groupName.c_str() );
        SMDS_MeshGroup& groupDS = group->SMDSGroup();

        if ( elemType == SMDSAbs_Node )
        {
          if ( zone.IsStructured() )
          {
            vector< cgsize_t > nodeIds;
            if ( psType == CGNS_ENUMV( PointRange ))
            {
              // nodes are given as (ijkMin, ijkMax)
              TPointRangeIterator idIt( & ids[0], meshDim );
              nodeIds.reserve( idIt.Size() );
              while ( idIt.More() )
                nodeIds.push_back( zone.NodeID( idIt.Next() ));
            }
            else
            {
              // nodes are given as (ijk1, ijk2, ..., ijkN)
              nodeIds.reserve( ids.size() / meshDim );
              for ( size_t i = 0; i < ids.size(); i += meshDim )
                nodeIds.push_back( zone.NodeID( ids[i], ids[i+1], ids[i+2] ));
            }
            ids.swap( nodeIds );
          }
          else if ( zone._nodeIdShift )
          {
            for ( size_t i = 0; i < ids.size(); ++i )
              ids[i] += zone._nodeIdShift;
          }
          zone.ReplaceNodes( &ids[0], ids.size() );

          for ( size_t i = 0; i < ids.size(); ++i )
            if ( const SMDS_MeshNode* n = myMesh->FindNode( ids[i] ))
              groupDS.Add( n );
        }
        else // BC applied to elements
        {
          if ( zone.IsStructured() )
          {
            int axis = 0; // axis perpendiculaire to which boundary elements are oriented
            if ( (int) ids.size() >= meshDim * 2 )
            {
              for ( ; axis < meshDim; ++axis )
                if ( ids[axis] - ids[axis+meshDim] == 0 )
                  break;
            }
            else
            {
              for ( ; axis < meshDim; ++axis )
                if ( normIndex[axis] != 0 )
                  break;
            }
            if ( axis == meshDim )
            {
              addMessage( SMESH_Comment("Invalid NormalIndex in BC ") << name );
              continue;
            }
            const int nbElemNodesByDim[] = { 1, 2, 4, 8 };
            const int nbElemNodes = nbElemNodesByDim[ meshDim ];

            if ( psType == CGNS_ENUMV( PointRange ) ||
                 psType == CGNS_ENUMV( ElementRange ))
            {
              // elements are given as (ijkMin, ijkMax)
              typedef void (TZoneData::*PGetNodesFun)( const gp_XYZ& ijk, cgsize_t* ids ) const;
              PGetNodesFun getNodesFun = 0;
              if ( elemType == SMDSAbs_Face  && meshDim == 3 )
                switch ( axis ) {
                case 0: getNodesFun = & TZoneData::IFaceNodes;
                case 1: getNodesFun = & TZoneData::JFaceNodes;
                case 2: getNodesFun = & TZoneData::KFaceNodes;
                }
              else if ( elemType == SMDSAbs_Edge && meshDim == 2 )
                switch ( axis ) {
                case 0: getNodesFun = & TZoneData::IEdgeNodes;
                case 1: getNodesFun = & TZoneData::JEdgeNodes;
                }
              if ( !getNodesFun )
              {
                addMessage( SMESH_Comment("Unsupported BC location in BC ") << name
                            << " " << cg_GridLocationName( location )
                            << " in " << meshDim << " mesh");
                continue;
              }
              TPointRangeIterator rangeIt( & ids[0], meshDim );
              vector< cgsize_t > elemNodeIds( rangeIt.Size() * nbElemNodes );
              for ( int i = 0; rangeIt.More(); i+= nbElemNodes )
                (zone.*getNodesFun)( rangeIt.Next(), &elemNodeIds[i] );

              ids.swap( elemNodeIds );
            }
            else
            {
              // elements are given as (ijk1, ijk2, ..., ijkN)
              typedef void (TZoneData::*PGetNodesFun)( int i, int j, int k, cgsize_t* ids ) const;
              PGetNodesFun getNodesFun = 0;
              if ( elemType == SMDSAbs_Face )
                switch ( axis ) {
                case 0: getNodesFun = & TZoneData::IFaceNodes;
                case 1: getNodesFun = & TZoneData::JFaceNodes;
                case 2: getNodesFun = & TZoneData::KFaceNodes;
                }
              else if ( elemType == SMDSAbs_Edge && meshDim == 2 )
                switch ( axis ) {
                case 0: getNodesFun = & TZoneData::IEdgeNodes;
                case 1: getNodesFun = & TZoneData::JEdgeNodes;
                }
              if ( !getNodesFun )
              {
                addMessage( SMESH_Comment("Unsupported BC location in BC ") << name
                            << " " << cg_GridLocationName( location )
                            << " in " << meshDim << " mesh");
                continue;
              }
              vector< cgsize_t > elemNodeIds( ids.size()/meshDim * nbElemNodes );
              for ( size_t i = 0, j = 0; i < ids.size(); i += meshDim, j += nbElemNodes )
                (zone.*getNodesFun)( ids[i], ids[i+1], ids[i+2], &elemNodeIds[j] );

              ids.swap( elemNodeIds );
            }
            zone.ReplaceNodes( &ids[0], ids.size() );

            PAddElemFun addElemFun = 0;
            switch ( meshDim ) {
            case 1: addElemFun = & add_BAR_2;
            case 2: addElemFun = & add_QUAD_4;
            case 3: addElemFun = & add_HEXA_8;
            }
            int elemID = meshInfo.NbElements();
            const SMDS_MeshElement* elem = 0;
            for ( size_t i = 0; i < ids.size(); i += nbElemNodes )
            {
              if ( iZone == 1 || !( elem = findElement( &ids[i], nbElemNodes, myMesh )))
                elem = addElemFun( &ids[i], myMesh, ++elemID );
              groupDS.Add( elem );
            }
          }
          else // unstructured zone
          {
            if ( zone._elemIdShift )
              for ( size_t i = 0; i < ids.size(); ++i )
                ids[i] += zone._elemIdShift;

            if ( psType == CGNS_ENUMV( PointRange ) && ids.size() == 2 )
            {
              for ( cgsize_t i = ids[0]; i <= ids[1]; ++i )
                if ( const SMDS_MeshElement* e = myMesh->FindElement( i ))
                  groupDS.Add( e );
            }
            else
            {
              for ( size_t i = 0; i < ids.size(); ++i )
                if ( const SMDS_MeshElement* e = myMesh->FindElement( ids[i] ))
                  groupDS.Add( e );
            }
          }
        } // end "BC applied to elements"

        // to have group type according to a real elem type
        group->SetType( groupDS.GetType() );

      } // loop on BCs of the zone
    }
    else
    {
      addMessage( cg_get_error() );
    }
  } // loop on the zones of a mesh


  // ------------------------------------------------------------------------
  // Make groups for multiple zones and remove free nodes at zone interfaces
  // ------------------------------------------------------------------------
  map< string, TZoneData >::iterator nameZoneIt = zonesByName.begin();
  for ( ; nameZoneIt != zonesByName.end(); ++nameZoneIt )
  {
    TZoneData& zone = nameZoneIt->second;
    if ( zone._nbElems == 0 ) continue;
    if ( zone._nbElems == meshInfo.NbElements() ) break; // there is only one non-empty zone

    // make a group
    SMDSAbs_ElementType elemType = myMesh->GetElementType( zone._elemIdShift + 1,
                                                           /*iselem=*/true );
    SMESHDS_Group* group = new SMESHDS_Group ( groupID++, myMesh, elemType );
    myMesh->AddGroup( group );
    group->SetStoreName( nameZoneIt->first.c_str() );
    SMDS_MeshGroup& groupDS = group->SMDSGroup();

    for ( int i = 1; i <= zone._nbElems; ++i )
      if ( const SMDS_MeshElement* e = myMesh->FindElement( i + zone._elemIdShift ))
        groupDS.Add( e );

    // remove free nodes
    map< int, int >::iterator nnRmKeepIt = zone._nodeReplacementMap.begin();
    for ( ; nnRmKeepIt != zone._nodeReplacementMap.end(); ++nnRmKeepIt )
      if ( const SMDS_MeshNode* n = myMesh->FindNode( nnRmKeepIt->first ))
        if ( n->NbInverseElements() == 0 )
          myMesh->RemoveFreeNode( n, (SMESHDS_SubMesh *)0, /*fromGroups=*/false );
  }

  aResult = myErrorMessages.empty() ? DRS_OK : DRS_WARN_SKIP_ELEM;

  return aResult;
}

//================================================================================
/*!
 * \brief Constructor
 */
//================================================================================

DriverCGNS_Read::DriverCGNS_Read()
{
  _fn = -1;
}
//================================================================================
/*!
 * \brief Close the cgns file at destruction
 */
//================================================================================

DriverCGNS_Read::~DriverCGNS_Read()
{
  if ( _fn > 0 )
    cg_close( _fn );
}

//================================================================================
/*!
 * \brief Opens myFile
 */
//================================================================================

Driver_Mesh::Status DriverCGNS_Read::open()
{
  if ( _fn < 0 )
  {
    
#ifdef CG_MODE_READ
    int res = cg_open(myFile.c_str(), CG_MODE_READ, &_fn);
#else
    int res = cg_open(myFile.c_str(), MODE_READ, &_fn);
#endif
    if ( res != CG_OK)
    {
      addMessage( cg_get_error(), /*fatal = */true );
    }
  }
  return _fn >= 0 ? DRS_OK : DRS_FAIL;
}

//================================================================================
/*!
 * \brief Reads nb of meshes in myFile
 */
//================================================================================

int DriverCGNS_Read::GetNbMeshes(Status& theStatus)
{
  if (( theStatus = open()) != DRS_OK )
    return 0;

  int nbases = 0;
  if(cg_nbases( _fn, &nbases) != CG_OK)
    theStatus = addMessage( cg_get_error(), /*fatal = */true );

  return nbases;
}
