/*===========================================================================*\
 *                                                                           *
 *                               OpenMesh                                    *
 *      Copyright (C) 2001-2013 by Computer Graphics Group, RWTH Aachen      *
 *                           www.openmesh.org                                *
 *                                                                           *
 *---------------------------------------------------------------------------* 
 *  This file is part of OpenMesh.                                           *
 *                                                                           *
 *  OpenMesh is free software: you can redistribute it and/or modify         * 
 *  it under the terms of the GNU Lesser General Public License as           *
 *  published by the Free Software Foundation, either version 3 of           *
 *  the License, or (at your option) any later version with the              *
 *  following exceptions:                                                    *
 *                                                                           *
 *  If other files instantiate templates or use macros                       *
 *  or inline functions from this file, or you compile this file and         *
 *  link it with other files to produce an executable, this file does        *
 *  not by itself cause the resulting executable to be covered by the        *
 *  GNU Lesser General Public License. This exception does not however       *
 *  invalidate any other reasons why the executable file might be            *
 *  covered by the GNU Lesser General Public License.                        *
 *                                                                           *
 *  OpenMesh is distributed in the hope that it will be useful,              *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 *  GNU Lesser General Public License for more details.                      *
 *                                                                           *
 *  You should have received a copy of the GNU LesserGeneral Public          *
 *  License along with OpenMesh.  If not,                                    *
 *  see <http://www.gnu.org/licenses/>.                                      *
 *                                                                           *
\*===========================================================================*/ 

/*===========================================================================*\
 *                                                                           *             
 *   $Revision: 922 $                                                         *
 *   $Date: 2013-08-11 12:26:11 +0200 (Ne, 11 aug 2013) $                   *
 *                                                                           *
\*===========================================================================*/

#ifndef OPENMESH_APPS_VDPMSTREAMING_CLIENT_MYMESH_HH
#define OPENMESH_APPS_VDPMSTREAMING_CLIENT_MYMESH_HH

#include <qthread.h>
#include <qmutex.h>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Tools/VDPM/MeshTraits.hh>


using OpenMesh::VDPM::MeshTraits;


//== CLASS DEFINITION =========================================================


typedef OpenMesh::TriMesh_ArrayKernelT<MeshTraits>  MyMesh;

static QMutex  mutex_;


//== CLASS DEFINITION =========================================================

#endif //OPENMESH_APPS_VDPMSTREAMING_CLIENT_MYMESH_HH defined
