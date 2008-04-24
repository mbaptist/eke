// -*- C++ -*-

//
// C++ Interface: eke
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

#ifndef EKE_HPP
#define EKE_HPP



#include "types.hpp"
#include "grid.hpp"

#include <string>

void poisson_boltzmann(std::string run_name);

//Distribute charges on the colloidal sphere
void distribute_colloidal_charge(RSF & density_colloid,Grid & grid,const Real & colloid_valence);
//Distribute ions of each ionic species	
void distribute_ions(RSF & density,Grid & grid,const double & ion_number);

//Evaluate total charge density
RSF total_charge_density(const RSF & density_colloid,
			 const std::vector<RSF> & ion_density,
			 const std::vector<int> & ion_valence);

#endif
