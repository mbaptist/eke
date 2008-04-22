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

#include "types.hpp"

//// BASIC MOVES ////

void charge_move(RVF & electric_field,
                 RSF & charges,
                 const Grid & grid,
                 const IV & node, const int & dir);

void loop_move(RVF & electric_field,
               const Grid & grid,
               const Loop & loop);

//// LATTICE SWEEPS ////

void sequential_sweep_concentration_moves(RVF & electric_field, 
                                   RSF & charges,
                                   const Grid & grid);

void sequential_sweep_loop_moves(RVF & electric_field, 
                                 const Grid & grid);
void random_sweep_loop_moves(RVF & electric_field, 
                             const Grid & grid);

//// INITIALISATIONS ////

void initialise_electric_field(RVF & electric_field,
                               const RSF & charges,
                               const Grid & grid);


//// FUNCTIONALS ////

//Functional for electrostitcs
Real functional(const RVF & electric_field);
//Functional for Poisson-Boltzmann
Real functional(const RSF & concentration,const RVF & electric_field);


#endif
