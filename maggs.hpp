// -*- C++ -*-

//
// C++ Interface: maggs
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

#ifndef MAGGS_HPP
#define MAGGS_HPP

#include "grid.hpp"

#include <blitz/blitz.h>
#include <blitz/array.h>
#include <blitz/tiny.h>
#include <blitz/tinyvec.h>

void charge_move(blitz::Array<blitz::TinyVector<double,3>,3> & electric_field,
                 blitz::Array<double,3> & charges,
                 const Grid & grid,
                 const blitz::TinyVector<int,3> & node, const int & dir);

void sequential_sweep_charge_moves(blitz::Array<blitz::TinyVector<double,3>,3> & electric_field, 
                                   blitz::Array<double,3> & charges,
                                   const Grid & grid);

void loop_move(blitz::Array<blitz::TinyVector<double,3>,3> & electric_field,
               const Grid & grid,
               const Loop & loop);

void sequential_sweep_loop_moves(blitz::Array<blitz::TinyVector<double,3>,3> & electric_field, 
                                 const Grid & grid);
void random_sweep_loop_moves(blitz::Array<blitz::TinyVector<double,3>,3> & electric_field, 
                             const Grid & grid);


#if 0
void loop_move(VField & electric_field,const Grid & grid,const Loop & loop);

void sequential_sweep_loop_moves(Array<TinyVector<double,3>,3> & electric_field, 
                                 const Grid & grid);
void random_sweep_loop_moves(Array<TinyVector<double,3>,3> & electric_field, 
                             const Grid & grid);
#endif


#endif