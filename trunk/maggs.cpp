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
#include "io.hpp"

using namespace ranlib;


//// INITIALISATIONS ////

//Distribute ions of each ionic species	
void distribute_ions(RSF & density,Grid & grid,const double & ion_number)
{
  cout << "Entering distribute_ions..." << endl;
  int nx(grid.nx());
  int ny(grid.ny());
  int nz(grid.nz());
  int np=count(grid.point_type()==1);
  density=0;
  Real delta_number=ion_number/np;
  cout << delta_number*np << endl;
  for (int i=0;i<grid.nx();++i)
    for (int j=0;j<grid.ny();++j)
      for (int k=0;k<grid.nz();++k)
        if(grid.point_type()(i,j,k)==1)
          density(i,j,k)=delta_number;
  density/=grid.deltav();  
  cout << "Ion number: " << ion_number << endl;
  cout << "Number of points of type 1: " << np << endl;
  cout << "delta_number: " << delta_number << endl;
  cout << delta_number*np << endl;
  cout << sum(density)*grid.deltav() << endl;
  cout << "...leaving distribute_ions." << endl;
}

//Initialise the elctric field
void initialise_electric_field(RVF & electric_field,
                               const RSF & fixed_charge_density,
                               const std::vector<RSF> & ion_density,
                               const std::vector<int> & ion_valence,
                               const Grid & grid)
{
  //Initilise the field with zeros
  electric_field=RV(0,0,0);
  //Evaluate the total charge density (tcd)
  RSF tcd(fixed_charge_density.copy());
  for (int n=0;n<ion_density.size();++n)
    tcd+=RSF(ion_valence[n]*ion_density[n]);
  RSF density(tcd.copy());
  //Check for charge neutrality
  double total_charge=sum(tcd)*grid.deltav();
  cout << "Total charge: " << total_charge << endl;
  if (fabs(total_charge)>1e-10)
  {
    cout << "The total charge must be zero, instead of "
      << total_charge << ".\n"
      << "Aborting..." << endl;
    exit(1);
  }
  //Get some necessary grid parameters
  int nx(grid.nx());
  int ny(grid.ny());
  int nz(grid.nz());
  Real dv=grid.deltav();
  //Initialise Ex
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
  //Initialise Ey
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
  //Initialise Ez
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
  if(count(fabs(divergence(electric_field,grid)-tcd)>1e-10)>0)
  {
    cout << count(fabs(divergence(electric_field,grid)-tcd)>1e-10) << endl;
    cout << sum(divergence(electric_field,grid)) << endl;
    cout << sum(tcd) << endl;
    cout << "The electric field doesn't obey the Gauss law. \n"
      << "Aborting..." << endl;
    exit(1);
  }
  
  cout << count(fabs(divergence(electric_field,grid)-tcd)>1e-15) << endl;
  
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

///////////////////////


//// LOOP MOVES ////

//Single loop move
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

//Lattice sweep
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

//////////////////////////////


//// ION MOVES ////

// ION MOVE ALONG SOME PATH //
Real ion_move(RVF & electric_field,
              RSF & concentration,
              const int & ion_valence,
              const Path & path,
              const Grid & grid)
{
//Define nodes
  IV node1(path(0));
  //  IV node2(path(path.size()-1));
  IV node2(path(1));
  IV dif=node2;
  dif-=node1;
  int dir;
  for(int ii=0;ii<3;++ii)
    if(dif[ii]==1)
      dir=ii;
  
//Do move only if both nodes
  //belong to region 1
  if ((grid.point_type()(node1[0],node1[1],node1[2])!=1)
      ||(grid.point_type()(node2[0],node2[1],node2[2])!=1))
    return 0;
  else
  {
    Real & c1 = concentration(node1[0],node1[1],node1[2]);
    Real & c2 = concentration(node2[0],node2[1],node2[2]);
    
    
    Real & e = electric_field(node1[0],node1[1],node1[2])[path.direction_axis(0)];
    

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
        //cout << a << " " << b << " " << b-a <<" " << c << " " << fc << endl;
        if((fc==0)||(b-a<1e-15))
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
#if 1
    if(deltac==0)
    {
      UniformOpen<Real> rg;
      deltac=-c2+rg.random()*(c1+c2);
      deltae=-ion_valence*deltac/deltas;
    }
#endif
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

//Variation of the functional for ion moves
Real deltafunc(const Real & deltac,const Real & c1, const Real & c2, const Real & e,const Real & deltae)
{
  if(c1==0&&c2==0)
    return (e+.5*deltae)*deltae-deltac*log(-deltac)+deltac*log(deltac);
  else if(c1==0)
    return (e+.5*deltae)*deltae+c2*log(1+deltac/c2)-deltac*log(-deltac)+deltac*log(c2+deltac);
  else if(c2==0)
    return (e+.5*deltae)*deltae+c1*log(1-deltac/c1)-deltac*log(c1-deltac)+deltac*log(deltac);
  
  else
    return (e+.5*deltae)*deltae+c1*log(1-deltac/c1)+c2*log(1+deltac/c2)-deltac*log(c1-deltac)+deltac*log(c2+deltac);
}

//Derivative with respect to the charge variation
//of the variation of the functional for ion moves
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

//Lattice sweep
Real sequential_sweep_ion_moves(RVF & electric_field,
                                RSF & concentration,
                                const int & ion_valence,
                                const Grid & grid)
{
  Real delta_func=0;
  for (int i=0;i<grid.nx();++i)
    for (int j=0;j<grid.ny();++j)
      for (int k=0;k<grid.nz();++k) 
      {
        DiscreteUniform<int> rg(3);
        int dir=rg.random();
        for (int m=0;m<3;++m)
        {
          Path path(IV(i,j,k),dir%3,grid);
	    //Perform the concentration move
          delta_func+=ion_move(electric_field,concentration,ion_valence,path,grid);
          ++dir;
        }
      }
  return delta_func;
}

////////////////////////////////////////////


//// MINIMISE ////

void minimise(RVF & electric_field, std::vector<RSF> & ion_density,const std::vector<int> & ion_valence,const Grid & grid, const int & savingstep, const Real & tolerance)
{
  int num_steps=0;
  Real func0=functional(electric_field,ion_density);
  cout << "Initial value of the functional: " << func0 << endl;
  Real func=func0;
  while(1)
  {
    ++num_steps;
    Real delta_func=0;
      //Field moves
    delta_func+=sequential_sweep_loop_moves(electric_field,grid);
      //Ion moves
    for (int n=0;n<ion_valence.size();++n)
      delta_func+=sequential_sweep_ion_moves(electric_field,ion_density[n],ion_valence[n],grid);
      //Update the functional
    func+=delta_func;
      //Print iteration infos
    cout << "Minimisation step: " << num_steps << "\t"
      << "Variation in functional: " << delta_func << "\t"
      << "Functional: " << func << endl;
  //Save fields
    if(num_steps%savingstep==0)
    {
      for (int n=0;n<ion_density.size();++n)
      {
        std::stringstream ss;
        ss << "ion_density_" << n << "_intermediate";
        vtkSave(ss.str(),ion_density[0],"density",grid);
      }
      vtkSave("electric_field_intermediate",electric_field,"electric_field",grid);
    }
      //Check for stop criterium
    if(fabs(delta_func)<tolerance)
      break;
  }
}


//// DIFFERENTIAL OPERATORS ////

//Divergence
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

///////////////////////////////////////////////
