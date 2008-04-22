// -*- C++ -*-

//
// C++ Interface: types
//
// Description: 
//
//
// Author: Manuel Baptista <mbaptist@gmail.com>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//

/*
Copyright (C) 2008 Manuel Baptista <mbaptist@gmail.com>

This file is part of eke.

Eke is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Eke is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with eke.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef TYPES_HPP
#define TYPES_HPP

typedef double Real;


//// BLITZ ////
#define USE_BLITZ
//#define USE_CAT

#ifdef USE_BLITZ

#include <blitz/blitz.h>
#include <blitz/array.h>
#include <blitz/tiny.h>
#include <blitz/tinyvec.h>

using namespace blitz;

BZ_DECLARE_FUNCTION_RET(norm,Real)
BZ_DECLARE_FUNCTION2_RET(dot,Real)
//BZ_DECLARE_FUNCTION_RET(log,Real)

//Types
//Integer vector
typedef blitz::TinyVector<int,3> IV;
//Integer scalar field
typedef blitz::Array<int,3> ISF;
//Real vector
typedef blitz::TinyVector<Real,3> RV;
//Real scalar field
typedef blitz::Array<Real,3> RSF;
//Real vector field
typedef blitz::Array<RV,3> RVF;



#endif
////    ////


//// CAT ////
#ifdef USE_CAT



#endif
////    ////



#endif
