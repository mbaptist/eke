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
	//Physical parameters
	int number_of_ionic_species=input.parse_int("number_of_ionic_species");
	double lb=input.parse_double("lb");
	double colloid_valence=input.parse_double("colloid_valence");
	//Numerical parameters
	//Grid parameters
	int nx=input.parse_int("nx");
	int ny=input.parse_int("ny");
	int nz=input.parse_int("nz");
	

	//Rescaling to non-dimensional units
	double k=sqrt(4.*M_PI*lb*colloid_valence/(lx*ly*lz));
	colloid_valence*=lb/rs;
	double colloid_charge_density=colloid_valence/(k*rs);
	lx*=k;
	ly*=k;
	lz*=k;
	rs*=k;
	double colloid_charge=colloid_charge_density*4.*M_PI*pow(rs,2);
	
   //Grid
	Grid grid(nx,ny,nz,lx,ly,lz,rs);
	
  //Colloidal charges	
	RSF density_colloid(nx,ny,nz);
	density_colloid=0;
	distribute_colloidal_charge(density_colloid,grid,colloid_charge);
	
	
  //Distribute ionic species
	RSF density_counterions(nx,ny,nz);
	density_counterions=0;
	distribute_ionic_species(density_counterions,grid,colloid_charge);
//	const N=number_of_ionic_species
	//Concentratrion
//	blitz::TinyVector<RSF,N>(RSF(nx,ny,nz));
	
#if 0	
	//Verify charge neutrality
	if (sum(density_colloid+density_counterions)>1e-10)
	{
		cout << "The total charge must be zero, instead of "
			<< sum(density_colloid+density_counterions) << ".\n"
			<< "Aborting..." << endl;
		exit(1);
	}
#endif	
	
  //electric field
	RVF electric_field(nx,ny,nz);
	initialise_electric_field(electric_field,RSF(density_colloid+density_counterions),grid);
	
	
  //Minimise
	int num_steps=0;
	Real func0=functional(density_counterions,electric_field);
	while(1)
	{
		++num_steps;
		cout << "Minimisation step " << num_steps << endl;
		sequential_sweep_loop_moves(electric_field,grid);
		sequential_sweep_concentration_moves(electric_field,density_counterions,grid);
		Real func=functional(density_counterions,electric_field);      
		Real delta_func=func-func0;
		cout << "Functional before // Functional after // Variation" << endl;
		cout <<  func0 << " // " << func << " // " << delta_func << endl << endl;
		if((delta_func<=0)&&(-delta_func<1e-16))
			break;
		func0=func;
	}
	
	
}

////////////////////////////////////////////////////////////////////////////////////////////

//Distribute charges on the colloidal sphere
void distribute_colloidal_charge(RSF & density_colloid,
                                 Grid & grid,
                                 const Real & total_charge)
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
	int np=10*static_cast<int>(4*M_PI*pow(rs,2)/grid.deltasx());
	//cout << np << endl;
	Real delta_charge=total_charge/np;
	cout << total_charge << endl;
	//cout << np*delta_charge << endl;
	for(int n=0;n<np;++n)
	{
		Real theta=rg.random()*M_PI;
		Real phi=rg.random()*2.*M_PI;
		RV coord_ps(rs*cos(phi)*sin(theta),
		            rs*sin(phi)*sin(theta),
		            rs*cos(theta));
		IV ind_np=grid.nearest_point_index(coord_ps);	
		density_colloid(ind_np[0],ind_np[1],ind_np[2])+=delta_charge/grid.deltav();
		grid.point_type()(ind_np[0],ind_np[1],ind_np[2])=2;
	//cout << sum(charges) << endl;
	}
#endif
	
#if 0
	int nx(grid.nx());
	int ny(grid.ny());
	int nz(grid.nz());
	Real lx(grid.lx());
	Real ly(grid.ly());
	Real lz(grid.lz());
	Real rs(grid.rs());
	
	grid.point_type()=0;
	
	int izero=nx/2;
	int ileft=izero-nx/4;
	int iright=izero+nx/4;
	
	cout << izero << endl;
	cout << ileft << endl;
	cout << iright << endl;
	
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
	Real delta_charge=total_charge/np;
	
	ISF density_count(nx,ny,nx);
	
	for (int i=0;i<grid.nx();++i)
		for (int j=0;j<grid.ny();++j)
			for (int k=0;k<grid.nz();++k)
				if(grid.point_type()(i,j,k)==2)
					density_count(i,j,k)+=1;
	
	density_colloid=density_count*delta_charge/grid.deltav();
	
	cout << sum(density_colloid) << endl;
	
#endif
	
	
	
}
//Distribute ionic species	
void distribute_ionic_species(RSF & density,
                              Grid & grid,
                              const Real & total_charge)
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
	Real delta_charge=-total_charge/np;
	for (int i=0;i<grid.nx();++i)
		for (int j=0;j<grid.ny();++j)
			for (int k=0;k<grid.nz();++k)
				if(grid.point_type()(i,j,k)==1)
					density(i,j,k)+=delta_charge/grid.deltav();
	cout << sum(density) << endl;
}

