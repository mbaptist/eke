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

#include "coord.hpp"

//// INITIALISATIONS ////
//Distribute ionic species
void distribute_ionic_species(std::vector<RSF> & ion_density,
                              const std::vector<double> & ion_number,
                              const string & runsname,
                              const Grid & grid);

//Initialise the electric field
void initialise_electric_field(RVF & electric_field,
                               const RSF & fixed_charge_density,
                               const std::vector<RSF> & ion_density,
                               const std::vector<int> & ion_valence,
                               const Real & eps,
                               const string & runsname,
                               const Grid & grid);


//// FUNCTIONAL ////
Real functional(const RVF & electric_field,const std::vector<RSF> ion_concentration);


//// LOOP MOVES ////
//Single loop move
Real loop_move(RVF & electric_field,const Loop & loop);
//Lattice sweep
Real sequential_sweep_loop_moves(RVF & electric_field,const Grid & grid);


//// ION MOVES ////
// ION MOVE ALONG SOME PATH //
Real ion_move(RVF & electric_field,
              RSF & concentration,
              const int & ion_valence,
              const IV & node1,
              const int & dir,
              const Real & eps,
              const Grid & grid);
//Variation of the functional for ion moves
Real deltafunc(const Real & deltac,const Real & deltae,const Real & c1, const Real & c2, const Real & e);
//Derivative with respect to the charge variation
//of the variation of the functional for ion moves
Real d_deltafunc_d_deltac(const Real & deltac,const Real & c1, const Real & c2, const Real & e, const Real & d_deltae_d_deltac);
//Lattice sweep
Real sequential_sweep_ion_moves(RVF & electric_field,
                                RSF & concentration,
                                const int & ion_valence,
                                const Real & eps,
                                const Grid & grid);
//// MINIMISE ////

void minimise(RVF & electric_field,
              std::vector<RSF> & ion_density,
              const std::vector<int> & ion_valence,
              const int & savingstep, 
              const Real & eps,
              const string & runsname,
              const Grid & grid);

//// DIFFERENTIAL OPERATORS ////
//Divergence
RSF divergence(const RVF & field,const Grid & grid);

#endif
