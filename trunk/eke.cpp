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
  colloid_valence*=lb/rs;
  for (int n=0;n<ion_number.size();++n)
  	ion_number[n]*=lb/rs;
  double k=sqrt(4.*M_PI*lb*colloid_valence/(lx*ly*lz));
  lx*=k;
  ly*=k;
  lz*=k;
  rs*=k;
	
  //Grid
  Grid grid(nx,ny,nz,lx,ly,lz,rs);
	
  //Colloidal charges	
  RSF density_colloid(nx,ny,nz);
  density_colloid=0;
  distribute_colloidal_charge(density_colloid,grid,colloid_valence);
	
  //Distribute ionic species	
  std::vector<RSF> ion_density;
  for (int n=0;n<ion_valence.size();++n)
    {
      ion_density.push_back(RSF(nx,ny,nz));
      ion_density[n]=0;
      distribute_ions(ion_density[n],grid,ion_number[n]);
    }
		
#if 0	
  //Verify charge neutrality
  double total_charge=sum(total_charge_density(density_colloid,ion_density,ion_valence))*grid.deltav();
  if (total_charge>1e-10)
    {
      cout << "The total charge must be zero, instead of "
	   << total_charge << ".\n"
	   << "Aborting..." << endl;
      exit(1);
    }
#endif	
	
  //electric field
  RVF electric_field(nx,ny,nz);
  initialise_electric_field(electric_field,total_charge_density(density_colloid,ion_density,ion_valence),grid);
	
  //Minimise
  int num_steps=0;
  while(1)
    {
      ++num_steps;
      cout << "Minimisation step: " << num_steps << "\t";
      //Ion moves
      Real delta_func=0;
      for (int n=0;n<ion_valence.size();++n)
	delta_func+=sequential_sweep_concentration_moves(electric_field,ion_density[n],ion_valence[n],grid);
      //Field moves
      delta_func+=sequential_sweep_loop_moves(electric_field,grid);
		
      cout << "Variation in functional: " << delta_func << endl;
      if((fabs(delta_func)<eps)&&(num_steps>100))
	break;
    }
}

////////////////////////////////////////////////////////////////////////////////////////////

//Distribute charges on the colloidal sphere
void distribute_colloidal_charge(RSF & density_colloid,Grid & grid,const Real & colloid_valence)
{
#if 1
  int nx(grid.nx());
  int ny(grid.ny());
  int nz(grid.nz());
  Real lx(grid.lx());
  Real ly(grid.ly());
  Real lz(grid.lz());
  Real rs(grid.rs());
  density_colloid=0;
  UniformClosed<Real> rg;
  //number of points to represent the sphere surface
  //100 per grid unit area
  int np=100*static_cast<int>(4*M_PI*pow(rs,2)/grid.deltasx());
  //cout << np << endl;
  Real delta_density=colloid_valence/np;
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
      density_colloid(ind_np[0],ind_np[1],ind_np[2])+=delta_density/grid.deltav();
      grid.point_type()(ind_np[0],ind_np[1],ind_np[2])=2;
    }
#endif
	
  //Distrubute the charges on a double plane
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
  int ileft=izero-nx/4;
  int iright=izero+nx/4;
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
  density_colloid=0;
  Real delta_density=colloid_valence/np;
  for (int i=0;i<grid.nx();++i)
    for (int j=0;j<grid.ny();++j)
      for (int k=0;k<grid.nz();++k)
	if(grid.point_type()(i,j,k)==2)
	  density_colloid(i,j,k)+=delta_density/grid.deltav();	
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
  int np=0;
  for (int i=0;i<grid.nx();++i)
    for (int j=0;j<grid.ny();++j)
      for (int k=0;k<grid.nz();++k)
	if(grid.point_type()(i,j,k)==1)
	  ++np;
  //cout << np << endl;
  density=0;
  Real delta_density=static_cast<double>(ion_number)/np;
  for (int i=0;i<grid.nx();++i)
    for (int j=0;j<grid.ny();++j)
      for (int k=0;k<grid.nz();++k)
	if(grid.point_type()(i,j,k)==1)
	  density(i,j,k)+=delta_density/grid.deltav();
  //cout << sum(density) << endl;
}

////////////////////////////////////////////////////////////////////////////////

//Evaluate total charge density
RSF total_charge_density(const RSF & density_colloid,
			 const std::vector<RSF> & ion_density,
			 const std::vector<int> & ion_valence)
{
  RSF tcd(density_colloid.copy());
  for (int n=0;n<ion_valence.size();++n)
    tcd+=RSF(ion_valence[n]*ion_density[n]);
  return tcd;
}
