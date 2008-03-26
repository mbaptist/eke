// -*- C++ -*-

//
// C++ Implementation: eke
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


#include "maggs.hpp"

#include <blitz/blitz.h>
#include <blitz/array.h>
#include <blitz/tiny.h>
#include <blitz/tinyvec.h>

using namespace blitz;

#include <blitz/random.h>
#include <random/uniform.h>
#include <random/discrete-uniform.h>

using namespace ranlib;



void charge_move(Array<TinyVector<double,3>,3> & electric_field,
                 Array<double,3> & charges,
                 const Grid & grid,
                 const TinyVector<int,3> & node, const int & dir)
{
  if (grid.point_type()(node[0],node[1],node[2])==1)
    {
      TinyVector<int,3> node2(node);
      node2[dir]+=1;

      double & q = charges(node[0],node[1],node[2]);
      double & q2 = charges(node2[0],node2[1],node2[2]);
      double & e =   electric_field(node[0],node[1],node[2])[dir];
      double ep=0;
      //Evaluate the change in the functional
      double delta_func=ep*ep-e*e;
      if (delta_func<0)
        {
          e=ep;
          q-=e;
          q2+=e;
        }
    }
}

void sequential_sweep_charge_moves(Array<TinyVector<double,3>,3> & electric_field,
                                   Array<double,3> & charges,
                                   const Grid & grid)
{
  double delta_func;
  for (int i=0;i<grid.nx();++i)
    for (int j=0;j<grid.ny();++j)
      for (int k=0;k<grid.nz();++k)
        for (int dir=0;dir<3;++dir)
          {
//Perform the charge move
            charge_move(electric_field,charges,grid,TinyVector<int,3>(i,j,k),dir);
          }
}

void loop_move(Array<TinyVector<double,3>,3> & electric_field,
               const Grid & grid,
               const Loop & loop)
{
//Evaluate the trial field
  double & e1 = electric_field(loop.node1()[0],loop.node1()[1],loop.node1()[2])[loop.dir1()];
  double & e2 = electric_field(loop.node2()[0],loop.node2()[1],loop.node2()[2])[loop.dir2()];
  double & e3 = electric_field(loop.node3()[0],loop.node3()[1],loop.node3()[2])[loop.dir1()];
  double & e4 = electric_field(loop.node4()[0],loop.node4()[1],loop.node4()[2])[loop.dir2()];
  double delta_field=-.25*(e1+e2-e3-e4);
  double e1p=e1+delta_field;
  double e2p=e2+delta_field;
  double e3p=e3-delta_field;
  double e4p=e4-delta_field;
  //Evaluate the change in the functional
  double delta_func=(e1p*e1p+e2p*e2p+e3p*e3p+e4p*e4p
                     -e1*e1-e2*e2-e3*e3-e4*e4);
  //cout << "Dfunc=" << delta_func << endl;
  //Accept/reject move
  if (delta_func<0)
    {
      e1=e1p;
      e2=e2p;
      e3=e3p;
      e4=e4p;
    }
  //cout << "Circulation: " << e1+e2-e3-e4 << endl;
};

void sequential_sweep_loop_moves(Array<TinyVector<double,3>,3> & electric_field,
                                 const Grid & grid)
{
  double delta_func;
  for (int i=0;i<grid.nx();++i)
    for (int j=0;j<grid.ny();++j)
      for (int k=0;k<grid.nz();++k)
        for (int loop_number=0;loop_number<3;++loop_number)
          {
//Create the loop
            Loop loop(blitz::TinyVector<int,3>(i,j,k),loop_number,grid);
//Perform the loop move
            loop_move(electric_field,grid,loop);
          }
}

void random_sweep_loop_moves(Array<TinyVector<double,3>,3> & electric_field,
                             const Grid & grid)
{
  //Choose 3x the number of nodes
  for (int nn=0;nn<3*electric_field.size();++nn)
    {
      //Choose a node
      DiscreteUniform<int> rgx(grid.nx());
      DiscreteUniform<int> rgy(grid.ny());
      DiscreteUniform<int> rgz(grid.nz());
      int ni = rgx.random();
      int nj = rgy.random();
      int nk = rgz.random();
      //cout << "Node: " << ni << " " << nj << " " << nk << endl;
      //Choose a loop
      DiscreteUniform<int> rgp(3);
      int loop_number=rgp.random();
      //cout << "Loop number: " << loop_number << endl;
      //Create the loop
      Loop loop(blitz::TinyVector<int,3>(ni,nj,nk),loop_number,grid);
      //Perform the loop move
      loop_move(electric_field,grid,loop);
    }
}
