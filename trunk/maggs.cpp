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
#include "types.hpp"
#include "random.hpp"

using namespace ranlib;



//// BASIC MOVES ////


Real concentration_move(RVF & electric_field,
                        RSF & concentration,
                        const int & ion_valence,
                        const Grid & grid,
                        const IV & node1, const int & dir)
{
  //Define the second node
  //from the first and the direction
  IV node2(node1);
  node2[dir]+=1;
  //Apply periodic boundary
  if(node2[0]==grid.nx())
    node2[0]=0;
  if(node2[1]==grid.ny())
    node2[1]=0;
  if(node2[2]==grid.nz())
    node2[2]=0;
  //Do move only if both nodes
  //belong to region 1
  if ((grid.point_type()(node1[0],node1[1],node1[2])!=1)
      ||(grid.point_type()(node2[0],node2[1],node2[2])!=1))
    return 0;
  else
  {
    Real & c1 = concentration(node1[0],node1[1],node1[2]);
    Real & c2 = concentration(node2[0],node2[1],node2[2]);
    Real & e = electric_field(node1[0],node1[1],node1[2])[dir];
    //Check if concentration is negative
    if((c1<0)||(c2<0))
    {
      cout << "Concentrations cannot be negative. Aborting..." << endl;
      exit(1);
    }
    Real deltas=grid.deltas(dir);
#if 1
    //Optimal concentration change
    //solve for the optimal concentration change with successive bisections
    Real deltac;
    Real a=-c2;
    Real b=c1;
    Real fa=d_deltafunc_d_deltac(a,c1,c2,e,deltas);
    Real fb=d_deltafunc_d_deltac(b,c1,c2,e,deltas);
    while(1)
    {
      Real c=.5*(a+b);
      Real fc=d_deltafunc_d_deltac(c,c1,c2,e,deltas);
      if(fa*fc<0)
        b=c;
      else
        a=c;
      if(b-a<1e-16)
      {
        deltac=c;
        //cout << deltac << endl;
        break;
      }
    }
#endif
#if 0
    //Optimal concentration change
    //solve for the optimal concentration change with fixed point iterations
    //Real deltac;
    Real deltac0=deltac;
    while(1)
    {
      deltac=deltac0-d_deltafunc_d_deltac(deltac0,c1,c2,e,deltas);
      if(fabs(deltac-deltac0)<.5e-16)
        break;
      deltac0=deltac;
    }
#endif
#if 0
    //Change the concentration by a random amount
    //in the interval (-c2,c1)
    UniformOpen<Real> rg;
    Real deltac=-c2+rg.random()*(c2+c1);
#endif
      //Evaluate the change in the functional
    Real delta_func=deltafunc(deltac,c1,c2,e,deltas);
      //cout << delta_func << endl;
//Accept the move if if minimises the functional
    if (delta_func<0)
    {
      e-=deltac/deltas;
      c1-=deltac;
      c2+=deltac;
      return delta_func;
    }
    else
      return 0;
  }
}


Real loop_move(RVF & electric_field,
               const Grid & grid,
               const Loop & loop)
{
  //Evaluate the trial field
  Real & e1 = electric_field(loop.node1()[0],loop.node1()[1],loop.node1()[2])[loop.dir1()];
  Real & e2 = electric_field(loop.node2()[0],loop.node2()[1],loop.node2()[2])[loop.dir2()];
  Real & e3 = electric_field(loop.node3()[0],loop.node3()[1],loop.node3()[2])[loop.dir1()];
  Real & e4 = electric_field(loop.node4()[0],loop.node4()[1],loop.node4()[2])[loop.dir2()];
  //Optimal move
  //This is valid both for electrostatics and Poisson-Boltzmann
  //It is characteristic of loop moves
  Real delta_field=-.25*(e1+e2-e3-e4);
  Real e1p=e1+delta_field;
  Real e2p=e2+delta_field;
  Real e3p=e3-delta_field;
  Real e4p=e4-delta_field;
  //Evaluate the change in the functional
  //This is valid both for electrostatics and Poisson-Boltzmann
  //It is characteristic of loop moves
  Real delta_func=(e1p*e1p+e2p*e2p+e3p*e3p+e4p*e4p
                   -e1*e1-e2*e2-e3*e3-e4*e4);
  //cout << "Dfunc=" << delta_func << endl;
  //Accept/reject move
  if (delta_func<0)
  {
    e1=e1p;
    e2=e2p;
    e3=e3p;
    e4=e4p;
    return delta_func;
      //cout << "Circulation: " << e1+e2-e3-e4 << endl;
  }
  else 
    return 0;
};



//// LATTICE SWEEPS ////


Real sequential_sweep_concentration_moves(RVF & electric_field,
                                          RSF & concentration,
                                          const int & ion_valence,
                                          const Grid & grid)
{
  Real delta_func=0;
  for (int i=0;i<grid.nx();++i)
    for (int j=0;j<grid.ny();++j)
      for (int k=0;k<grid.nz();++k)
        for (int dir=0;dir<3;++dir)
        {
	    //Perform the concentration move
          delta_func+=concentration_move(electric_field,concentration,ion_valence,grid,IV(i,j,k),dir);
        }
  return delta_func;
}


Real sequential_sweep_loop_moves(RVF & electric_field,
                                 const Grid & grid)
{
  Real delta_func=0;
  for (int i=0;i<grid.nx();++i)
    for (int j=0;j<grid.ny();++j)
      for (int k=0;k<grid.nz();++k)
        for (int loop_number=0;loop_number<3;++loop_number)
        {
	    //Create the loop
          Loop loop(blitz::TinyVector<int,3>(i,j,k),loop_number,grid);
	    //Perform the loop move
          delta_func+=loop_move(electric_field,grid,loop);
        }
  return delta_func;
}



//// INITIALISATIONS ////

void initialise_electric_field(RVF & electric_field,
                               const RSF & total_density,
                               const Grid & grid)
{
  int nx(grid.nx());
  int ny(grid.ny());
  int nz(grid.nz());
  electric_field=RV(0,0,0);
  RSF density(total_density.copy());
  Real mean_density;
  double dv=grid.deltav();
  //Ex
  Real ds=grid.deltasx();
  mean_density=mean(density(0,Range::all(),Range::all()));
  density(0,Range::all(),Range::all())-=mean_density;
  for(int j=0;j<ny;++j)
    for(int k=0;k<nz;++k)
      electric_field[0](0,j,k)=mean_density*dv/ds;
  for(int i=1;i<nx-1;++i)
  {
    mean_density=mean(density(i,Range::all(),Range::all()));
    density(i,Range::all(),Range::all())-=mean_density;
    for(int j=0;j<ny;++j)
      for(int k=0;k<nz;++k)
        electric_field[0](i,j,k)=electric_field[0](i-1,j,k)+mean_density*dv/ds;
  }
  mean_density=mean(density(nx-1,Range::all(),Range::all()));
  density(nx-1,Range::all(),Range::all())-=mean_density;
  for(int j=0;j<ny;++j)
    for(int k=0;k<nz;++k)
      electric_field[0](nx-1,j,k)=0.;
  //Ey
  ds=grid.deltasy();
  for(int i=0;i<nx;++i)
  {
    mean_density=mean(density(i,0,Range::all()));
    density(i,0,Range::all())-=mean_density;
    for(int k=0;k<nz;++k)
      electric_field[1](i,0,k)=mean_density*dv/ds;
    for(int j=1;j<ny-1;++j)
    {
      mean_density=mean(density(i,j,Range::all()));
      density(i,j,Range::all())-=mean_density;
      for(int k=0;k<nz;++k)
        electric_field[1](i,j,k)=electric_field[1](i,j-1,k)+mean_density*dv/ds;;
    }
    mean_density=mean(density(i,ny-1,Range::all()));
    density(i,ny-1,Range::all())-=mean_density;
    for(int k=0;k<nz;++k)
      electric_field[1](i,ny-1,k)=0.;
  }
  //Ez
  ds=grid.deltasz();
  for(int i=0;i<nx;++i)
    for(int j=0;j<ny;++j)
    {
      electric_field[2](i,j,0)=density(i,j,0)*dv/ds;
      for(int k=1;k<nz-1;++k)
        electric_field[2](i,j,k)=electric_field[2](i,j,k-1)+density(i,j,k)*dv/ds;
      electric_field[2](i,j,nz-1)=0;
    }
  
#if 1
  //Test for div E = \rho
  for(int i=0;i<nx;++i)
    for(int j=0;j<ny;++j)
      for(int k=0;k<nz;++k)
      {
        int n1,n2;
        n1=i-1;
        if (n1==-1)
          n1=nx-1;
        n2=i;
        Real divx=electric_field[0](n2,j,k)-electric_field[0](n1,j,k);
        divx/=grid.deltax();
        n1=j-1;
        if (n1==-1)
          n1=ny-1;
        n2=j;
        Real divy=electric_field[1](i,n2,k)-electric_field[1](i,n1,k);
        divy/=grid.deltay();
        n1=k-1;
        if (n1==-1)
          n1=nz-1;
        n2=k;
        Real divz=electric_field[2](i,j,n2)-electric_field[2](i,j,n1);
        divz/=grid.deltaz();
        Real div=divx+divy+divz;
        if (fabs(div-total_density(i,j,k))>.5e-10)
        {
          cout << i << " " << j << " " << k << " " 
            << divx << " " << divy <<" " << divz << "\n " 
            << div << " " << total_density(i,j,k) << endl;
          exit(1);
        }
      }
#endif 
  
  
}


//// FUNCTIONAL ////

Real functional(RVF & electric_field,std::vector<RSF> ion_concentration,std::vector<int> ion_valence)
{
  Real func=.5*sum(dot(electric_field,electric_field));
  for (int n=0;n<ion_valence.size();++n)
    func-=sum(where(ion_concentration[n]>0,ion_concentration[n]*log(ion_concentration[n]),0));
  return func;
}


Real deltafunc(const Real & deltac,const Real & c1, const Real & c2, const Real & e,const Real & deltas)
{
  return (-e*deltac/deltas+.5*deltac*deltac)/deltas+(c1*log(1-deltac/c1)+c2*log(1+deltac/c2)
    -deltac*log(c1-deltac)+deltac*log(c2+deltac));
}

Real d_deltafunc_d_deltac(const Real & deltac,const Real & c1, const Real & c2, const Real & e,const Real & deltas)
{
  Real earg=(-e+deltac/deltas)/deltas;
  Real exp_earg=exp(-fabs(earg));
  //cout << earg << " " << exp_earg << endl;
  if (earg<0)
    return deltac-(c1-c2*exp_earg)/(1+exp_earg);
  else if(earg>0)
    return deltac-(c1*exp_earg-c2)/(exp_earg+1);
  else
    return deltac-(c1-c2)/2.;
}


