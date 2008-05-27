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


#include<iostream>
#include<fstream>
#include<string>
#include<vector>

using namespace std;

#include "types.hpp"
#include "random.hpp"
#include "eke.hpp"
#include "maggs.hpp"
#include "io.hpp"


///////////////////////////////////////////

int main()
{	
  poisson_boltzmann("default.run");
}


void poisson_boltzmann(std::string run_name)
{
  
  cout << fixed << setprecision(16);
  
  //Input
  //Read input file
  PyInputParser input(run_name.c_str());
  //Physical parameters
  //Geometry parameters
  double lx=input.parse_double("lx");
  double ly=input.parse_double("ly");
  double lz=input.parse_double("lz");
  double rs=input.parse_double("rs");
  //Colloid
  double lb=input.parse_double("lb");
  double colloid_valence=input.parse_double("colloid_valence");
  //Ions
  vector<int> ion_valence=input.parse_vector_int("ion_valence");
  vector<double> ion_number=input.parse_vector_double("ion_number");
  //Numerical parameters
  //Grid parameters
  int nx=input.parse_int("nx");
  int ny=input.parse_int("ny");
  int nz=input.parse_int("nz");
  //Tolerance
  double eps=input.parse_double("eps");
  
  //Rescaling to non-dimensional units
  double k=sqrt(4.*M_PI*lb*fabs(colloid_valence)/(lx*ly*lz));
  colloid_valence*=lb/rs;
  for (int n=0;n<ion_number.size();++n)
    ion_number[n]*=lb/rs;
  lx*=k;
  ly*=k;
  lz*=k;
  rs*=k;
  
  //Grid
  Grid grid(nx,ny,nz,lx,ly,lz,rs);
  
  //Colloidal charges	
  RSF colloid_charge_density(nx,ny,nz);
  colloid_charge_density=0;
  distribute_colloidal_charge(colloid_charge_density,grid,colloid_valence);
  
  //Distribute ionic species	
  std::vector<RSF> ion_density;
  for (int n=0;n<ion_number.size();++n)
    {
      ion_density.push_back(RSF(nx,ny,nz));
      ion_density[n]=0;
      distribute_ions(ion_density[n],grid,ion_number[n]);
      std::stringstream ss;
      ss << "ion_density_" << n << "_initial";
      vtkSave(ss.str(),ion_density[0],"density",grid);
    }

  //Verify charge neutrality
  double total_charge=sum(total_charge_density(colloid_charge_density,ion_density,ion_valence))*grid.deltav();
  if (total_charge>1e-10)
    {
      cout << "The total charge must be zero, instead of "
	   << total_charge << ".\n"
	   << "Aborting..." << endl;
      exit(1);
    }
  
  //electric field
  RVF electric_field(nx,ny,nz);
  initialise_electric_field(electric_field,total_charge_density(colloid_charge_density,ion_density,ion_valence),grid);
  vtkSave("electric_field_initial",electric_field,"electric_field",grid);

  //Minimise
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
	delta_func+=sequential_sweep_concentration_moves(electric_field,ion_density[n],ion_valence[n],grid);
      //Update the functional
      func+=delta_func;
      //Print iteration infos
      cout << "Minimisation step: " << num_steps << "\t"
	   << "Variation in functional: " << delta_func << "\n"
	   << func << " " << func-func0 << endl;
      if(num_steps%100==0)
	{
	  for (int n=0;n<ion_number.size();++n)
	    {
	      std::stringstream ss;
	      ss << "ion_density_" << n << "_intermediate";
	      vtkSave(ss.str(),ion_density[0],"density",grid);
	    }
	  vtkSave("electric_field_intermediate",electric_field,"electric_field",grid);
	}
      //Check for stop criterium
      if(fabs(delta_func)<1e-16)
	break;
    }


  for (int n=0;n<ion_number.size();++n)
    {
      std::stringstream ss;
      ss << "ion_density_" << n;
      vtkSave(ss.str(),ion_density[0],"density",grid);
    }
  vtkSave("electric_field",electric_field,"electric_field",grid);
  
  cout << sum(electric_field[0]) << endl;
  cout << sum(electric_field[1]) << endl;
  cout << sum(electric_field[2]) << endl;
}

////////////////////////////////////////////////////////////////////////////////////////////

//Distribute charges on the colloidal sphere
void distribute_colloidal_charge(RSF & colloid_charge_density,Grid & grid,const Real & colloid_valence)
{
#if 1
  int nx(grid.nx());
  int ny(grid.ny());
  int nz(grid.nz());
  Real lx(grid.lx());
  Real ly(grid.ly());
  Real lz(grid.lz());
  Real rs(grid.rs());
  colloid_charge_density=0;
  UniformClosed<Real> rg;
  //number of points to represent the sphere surface
  //100 per grid unit area
  int np=100*static_cast<int>(4*M_PI*pow(rs,2)/grid.deltasx());
  //cout << np << endl;
  Real delta_charge=colloid_valence/np;
  //cout << colloid_valence << endl;
  //cout << np*delta_density << endl;
  for(int n=0;n<np;++n)
    {
      Real theta=rg.random()*M_PI;
      Real phi=rg.random()*2.*M_PI;
      RV coord_ps(rs*cos(phi)*sin(theta),
		  rs*sin(phi)*sin(theta),
		  rs*cos(theta));
      IV ind_np=grid.nearest_point_index(coord_ps);	
      colloid_charge_density(ind_np[0],ind_np[1],ind_np[2])+=delta_charge/grid.deltav();
      grid.point_type()(ind_np[0],ind_np[1],ind_np[2])=2;
    }
#endif
  
  //Distribute the charges on a double plane
  //For testing
  //Redifines the point types so that the ionic species are
  //distributed inside
#if 0
  int nx(grid.nx());
  int ny(grid.ny());
  int nz(grid.nz());
  Real lx(grid.lx());
  Real ly(grid.ly());
  Real lz(grid.lz());
  grid.point_type()=0;
  int izero=nx/2;
  int ileft=izero-3;
  int iright=izero+3;
  //cout << izero << endl;
  //cout << ileft << endl;
  //cout << iright << endl;
  for (int i=0;i<grid.nx();++i)
    for (int j=0;j<grid.ny();++j)
      for (int k=0;k<grid.nz();++k)
        if((i>ileft)&&(i<iright))
          grid.point_type()(i,j,k)=1;
	else if ((i==ileft)||(i==iright))
	  grid.point_type()(i,j,k)=2;
	else
	  grid.point_type()(i,j,k)=3;
  int np=2*ny*nz;
  colloid_charge_density=0;
  Real delta_charge=colloid_valence/np;
  for (int i=0;i<grid.nx();++i)
    for (int j=0;j<grid.ny();++j)
      for (int k=0;k<grid.nz();++k)
        if(grid.point_type()(i,j,k)==2)
          colloid_charge_density(i,j,k)+=delta_charge/grid.deltav();	
#endif
}
//Distribute ions of each ionic species	
void distribute_ions(RSF & density,Grid & grid,const double & ion_number)
{
  int nx(grid.nx());
  int ny(grid.ny());
  int nz(grid.nz());
  Real lx(grid.lx());
  Real ly(grid.ly());
  Real lz(grid.lz());
  int np=count(grid.point_type()==1);
  //cout << np << endl;
  density=0;
  Real delta_number=ion_number/np;
  for (int i=0;i<grid.nx();++i)
    for (int j=0;j<grid.ny();++j)
      for (int k=0;k<grid.nz();++k)
        if(grid.point_type()(i,j,k)==1)
          density(i,j,k)+=delta_number/grid.deltav();
  //cout << sum(density) << endl;
}

////////////////////////////////////////////////////////////////////////////////

//Evaluate total charge density
RSF total_charge_density(const RSF & colloid_charge_density,
                         const std::vector<RSF> & ion_density,
                         const std::vector<int> & ion_valence)
{
  cout << sum(ion_density[0]) << " " << ion_valence[0] << " " << sum(colloid_charge_density) << endl;
  RSF tcd(colloid_charge_density.copy());
  for (int n=0;n<ion_density.size();++n)
    tcd+=RSF(ion_valence[n]*ion_density[n]);
  return tcd;
}
