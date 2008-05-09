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

#include <vector>

//// BASIC MOVES ////

Real density_move(RVF & electric_field,
		  RSF & density,
		  const int & ion_valence,
		  const Grid & grid,
		  const IV & node, const int & dir);

Real loop_move(RVF & electric_field,
               const Grid & grid,
               const Loop & loop);

//// LATTICE SWEEPS ////

Real sequential_sweep_concentration_moves(RVF & electric_field, 
					  RSF & density,
					  const int & ion_valence,
					  const Grid & grid);

Real sequential_sweep_loop_moves(RVF & electric_field, 
                                 const Grid & grid);


//// INITIALISATIONS ////

void initialise_electric_field(RVF & electric_field,
                               const RSF & charges,
                               const Grid & grid);

//// FUNCTIONAL ////

Real functional(RVF & electric_field,std::vector<RSF> ion_concentration,std::vector<int> ion_valence);

Real d_deltafunc_d_deltac(const Real & deltac,const Real & c1, const Real & c2, const Real & e,const Real & deltas);

#endif
