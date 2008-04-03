// -*- C++ -*-

//
// C++ Implementation: grid
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


#include "grid.hpp"

#include "types.hpp"


//Ctor
Grid::Grid(const int & _nx_, const int & _ny_,const int & _nz_,
           const double & _lx_, const double & _ly_,const double & _lz_,
           const double & _rs_):
    nx_(_nx_),ny_(_ny_),nz_(_nz_),
    lx_(_lx_),ly_(_ly_),lz_(_lz_),
    rs_(_rs_),
    coordinates_(nx_,ny_,nz_),
    point_type_(nx_,ny_,nz_)
{
  for (int i=0;i<nx_;++i)
    for (int j=0;j<ny_;++j)
      for (int k=0;k<nz_;++k)
        {
          coordinates_(i,j,k)=TinyVector<double,3>(-.5*lx_+i*lx_/nx_,
                              -.5*ly_+j*ly_/ny_,
                              -.5*lz_+k*lz_/nz_);
          if (norm(coordinates_(i,j,k))>=2*rs_)
            point_type_(i,j,k)=1;
          else if ((norm(coordinates_(i,j,k))<2*rs_)&&(norm(coordinates_(i,j,k))>=rs_))
            point_type_(i,j,k)=3;
          else
            point_type_(i,j,k)=5;
        }
}

//Dtor
Grid::~Grid()
{
}

