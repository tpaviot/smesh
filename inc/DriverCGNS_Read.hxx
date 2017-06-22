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
// File      : DriverCGNS_Read.hxx
// Created   : Thu Jun 30 10:25:09 2011
// Author    : Edward AGAPOV (eap)

#ifndef __DriverCGNS_Read_HXX__
#define __DriverCGNS_Read_HXX__

#include "SMESH_DriverCGNS.hxx"

#include "Driver_SMESHDS_Mesh.h"

#include <vector>
#include <string>

/*!
 * \brief Driver reading a mesh from the CGNS file. The mesh to read is selected by 
 *  an index (counted from 0) set via SetMeshId()
 */
class MESHDriverCGNS_EXPORT DriverCGNS_Read : public Driver_SMESHDS_Mesh
{
public:

  DriverCGNS_Read();
  ~DriverCGNS_Read();

  virtual Status Perform();

  int GetNbMeshes(Status& theStatus);


private:

  Status open();

  int _fn; //!< file index
};

#endif
