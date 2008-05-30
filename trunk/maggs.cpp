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
    Real deltac;
    Real deltae;
    //Optimal concentration change
    //solve for the optimal concentration change with successive bisections
    Real a=-c2;
    Real fa=d_deltafunc_d_deltac(a,c1,c2,e,-ion_valence*a/deltas);
    Real b=c1;
    Real fb=d_deltafunc_d_deltac(b,c1,c2,e,-ion_valence*b/deltas);
    //cout << fa << " " << fb << endl;
    if(fa*fb<0)
    {
      while(1)
      {
        Real c=.5*(a+b);
        Real fc=d_deltafunc_d_deltac(c,c1,c2,e,-ion_valence*c/deltas);
        if((fc==0)||(b-a<1e-16))
        {
          deltac=c;
          deltae=-ion_valence*deltac/deltas;
        //cout << deltac << endl;
          break;
        }
        else if(fa*fc<0)
        {
          b=c;
          fb=fc;
        }
        else if(fb*fc<0)
        {
          a=c;
          fa=fc;
        }
        else
        {
        //cout << "Successive bisection method failed" << endl;
          deltac=0;
          break;
        }    
      }   
    }
    else
      deltac=0;
  //Try a random change if the optimal change
    //could not be determined or is equal to 0
    if(deltac==0)
    {
       UniformOpen<Real> rg;
      deltac=-c2+rg.random()*(c1+c2);
    }
      //Evaluate the change in the functional
    Real delta_func=deltafunc(deltac,c1,c2,e,deltae);
    //cout << "delta_func: " << delta_func << endl;
//Accept the move if if minimises the functional
    if (delta_func<0)
    {
      //cout << "delta_func: " << delta_func << endl;
      e+=deltae;
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
  Real & e3 = electric_field(loop.node4()[0],loop.node4()[1],loop.node4()[2])[loop.dir1()];
  Real & e4 = electric_field(loop.node1()[0],loop.node1()[1],loop.node1()[2])[loop.dir2()];
  //Optimal move
  Real delta_field=-.25*(e1+e2-e3-e4);
  //Evaluate the change in the functional
  Real delta_func=delta_field*(e1+e2-e3-e4+2.*delta_field);
  //cout << "Dfunc=" << delta_func << endl;
  //Accept/reject move
  if (delta_func<0)
  {
    e1+=delta_field;
    e2+=delta_field;
    e3-=delta_field;
    e4-=delta_field;
    return delta_func;
  }
  else 
    return 0;
}



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
                               const RSF & total_charge_density,
                               const Grid & grid)
{
  electric_field=RV(0,0,0);
  RSF density(total_charge_density.copy());
  int nx(grid.nx());
  int ny(grid.ny());
  int nz(grid.nz());
  Real dv=grid.deltav();
  //Ex
  Real mean_density=mean(density(0,Range::all(),Range::all()));
  density(0,Range::all(),Range::all())-=mean_density;
  for(int j=0;j<ny;++j)
    for(int k=0;k<nz;++k)
      electric_field[0](0,j,k)=0;
  for(int i=1;i<nx;++i)
  {
    mean_density=mean(density(i,Range::all(),Range::all()));
    density(i,Range::all(),Range::all())-=mean_density;
    for(int j=0;j<ny;++j)
      for(int k=0;k<nz;++k)
        electric_field[0](i,j,k)=electric_field[0](i-1,j,k)+mean_density;
  }
  //geometric factors for Ex
  Real ds=grid.deltasx();
  electric_field[0]*=dv/ds;
  //Ey
  for(int i=0;i<nx;++i)
  {
    mean_density=mean(density(i,0,Range::all()));
    density(i,0,Range::all())-=mean_density;
    for(int k=0;k<nz;++k)
      electric_field[1](i,0,k)=0;
    for(int j=1;j<ny;++j)
    {
      mean_density=mean(density(i,j,Range::all()));
      density(i,j,Range::all())-=mean_density;
      for(int k=0;k<nz;++k)
        electric_field[1](i,j,k)=electric_field[1](i,j-1,k)+mean_density;
    }
  }
    //geometric factors for Ey
  ds=grid.deltasy();
  electric_field[1]*=dv/ds;
  //Ez
  for(int i=0;i<nx;++i)
    for(int j=0;j<ny;++j)
    {
      electric_field[2](i,j,0)=0;
      for(int k=1;k<nz;++k)
        electric_field[2](i,j,k)=electric_field[2](i,j,k-1)+density(i,j,k);
    }  
  //geometric factors for Ez
  ds=grid.deltasz();
  electric_field[2]*=dv/ds;
  
  //cout << sum(divergence(electric_field,grid))*dv << endl;
  //cout << sum(total_charge_density)*dv << endl;
      
  //Test for div E = \rho
  if(count(fabs(divergence(electric_field,grid)-total_charge_density)>1e-10)>0)
  {
    cout << count(fabs(divergence(electric_field,grid)-total_charge_density)>1e-10) << endl;
    cout << sum(divergence(electric_field,grid)) << endl;
    cout << sum(total_charge_density) << endl;
    cout << "The electric field doesn't obey the Gauss law. \n"
      << "Aborting..." << endl;
    exit(1);
  }
  
  cout << count(fabs(divergence(electric_field,grid)-total_charge_density)>1e-15) << endl;
  
  //exit(0);
  
}


//// FUNCTIONAL ////

Real functional(RVF & electric_field,std::vector<RSF> ion_concentration)
{
  Real func=.5*sum(dot(electric_field,electric_field));
  for (int n=0;n<ion_concentration.size();++n)
    func+=sum(where(ion_concentration[n]>0,ion_concentration[n]*log(ion_concentration[n]),0));
  return func;
}


Real deltafunc(const Real & deltac,const Real & c1, const Real & c2, const Real & e,const Real & deltae)
{
  return (e+.5*deltae)*deltae+(c1*log(1-deltac/c1)+c2*log(1+deltac/c2)
                                 -deltac*log(c1-deltac)+deltac*log(c2+deltac));
}

Real d_deltafunc_d_deltac(const Real & deltac,const Real & c1, const Real & c2, const Real & e,const Real & deltae)
{
  Real earg=(e+deltae)*deltae/deltac;
  Real exp_earg=exp(-fabs(earg));
  //cout << earg << " " << exp_earg << endl;
  if (earg<0)
    return deltac-(c1-c2*exp_earg)/(1+exp_earg);
  else if(earg>0)
    return deltac-(c1*exp_earg-c2)/(exp_earg+1);
  else
    return deltac-(c1-c2)/2.;
}

//// DIFFERENTIAL OPERATORS ////

RSF divergence(const RVF & field,const Grid & grid)
{
  RSF div(field.shape());
  div=0;
  int nx=grid.nx();
  int ny=grid.ny();
  int nz=grid.nz();
  for(int i=0;i<nx;++i)
    for(int j=0;j<ny;++j)
      for(int k=0;k<nz;++k)
      {
        Real divx,divy,divz;
        if (i==0)
          divx=field[0](i,j,k)-field[0](nx-1,j,k);
        else
          divx=field[0](i,j,k)-field[0](i-1,j,k);
        divx/=grid.deltax();
        if (j==0)
          divy=field[1](i,j,k)-field[1](i,ny-1,k);
        else
          divy=field[1](i,j,k)-field[1](i,j-1,k);
        divy/=grid.deltay();
        if (k==0)
          divz=field[2](i,j,k)-field[2](i,j,nz-1);
        else
          divz=field[2](i,j,k)-field[2](i,j,k-1);
        divz/=grid.deltaz();
        div(i,j,k)=divx+divy+divz;
      }
  return div;
}


