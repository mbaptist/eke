// -*- C++ -*-

//
// C++ Implementation: p_spherical_capacitor
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
#include "maggs.hpp"
#include "io.hpp"


///////////////////////////////////////////

int main(int argc, char * argv[])
{
  string run_name;
  if(argc>1)
    run_name=argv[1];
  else
    run_name="p_spherical_capacitor.cfg";
  
  cout << fixed << setprecision(16);
  
  //Input
  //Open input file
  PyInputParser input(run_name.c_str());
  //Runsname
  string runsname=input.parse_string("runsname");
  //Physical parameters
  //Geometry parameters
  double lx=input.parse_double("lx");
  double ly=input.parse_double("ly");
  double lz=input.parse_double("lz");
  double ris=input.parse_double("ris");
  double ros=input.parse_double("ros");
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
  //Saving step
  int savingstep=input.parse_int("savingstep");
  //Close input file
  input.close();
  
  //Rescaling to non-dimensional units
  double k=sqrt(4.*M_PI*lb*fabs(colloid_valence)/(lx*ly*lz));
  colloid_valence*=4.*M_PI*lb*k;
  for (int n=0;n<ion_number.size();++n)
    ion_number[n]*=4.*M_PI*lb*k;
  lx*=k;
  ly*=k;
  lz*=k;
  ris*=k;
  ros*=k;

  //Grid
  Grid grid(nx,ny,nz,lx,ly,lz);
  
  //Classify Grid Points
  //There are no movable charges (ionic species)
  
  //Fixed charges	
  RSF fixed_charge_density(nx,ny,nz);
  fixed_charge_density=0;
  UniformClosedOpen<Real> rgco;
  UniformClosed<Real> rgc;
  //Inner sphere  
  //number of points to represent the sphere surface
  //1000 per grid unit area
  int np=1000*static_cast<int>(4*M_PI*pow(ris,2)/grid.deltasx());
  Real delta_density=(colloid_valence/grid.deltav())/np;
  for(int n=0;n<np;++n)
  {
    Real z=ris*(-1.+rgc.random()*2.);
    Real phi=rgco.random()*2.*M_PI;
    Real rsintheta=sqrt(ris*ris-z*z);
    RV coord_ps(rsintheta*cos(phi),
                rsintheta*sin(phi),
                z);
    IV ind_np=grid.nearest_point_index(coord_ps);	
    fixed_charge_density(ind_np)+=delta_density;
    grid.point_type()(ind_np)=2;
  } 
  //Outer sphere  
  //number of points to represent the sphere surface
  //1000 per grid unit area
  np=1000*static_cast<int>(4*M_PI*pow(ros,2)/grid.deltasx());
  delta_density=-(colloid_valence/grid.deltav())/np;
  for(int n=0;n<np;++n)
  {
    Real z=ros*(-1.+rgc.random()*2.);
    Real phi=rgco.random()*2.*M_PI;
    Real rsintheta=sqrt(ros*ros-z*z);
    RV coord_ps(rsintheta*cos(phi),
                rsintheta*sin(phi),
                z);
    IV ind_np=grid.nearest_point_index(coord_ps);	
    fixed_charge_density(ind_np)+=delta_density;
    grid.point_type()(ind_np)=2;
  } 

  
  //Save the fixed charge density
  std::stringstream ss;
  ss << runsname << "_fixed_charge_density";
  vtkSave(ss.str(),fixed_charge_density,"fixed_charge_density",grid);
  
  //Distribute ionic species	
  //There are no movable charges (ionic species)
  //Passing an empty vector to the routines,
//with the size zero defined in the config file
  std::vector<RSF> ion_density;
//  ion_density.push_back(RSF(nx,ny,nz));
//  ion_density[0]=0;
  
  //electric field
  RVF electric_field(nx,ny,nz);
  initialise_electric_field(electric_field,fixed_charge_density,ion_density,ion_valence,eps,runsname,grid);
  
//  RV avef(sum(electric_field[0]),sum(electric_field[1]),sum(electric_field[2]));
//  avef*=grid.deltav()/(lx*ly*lz);
//  electric_field-=avef;

  RVF electric_field_init(electric_field.copy());

  //Minimise
minimise(electric_field,ion_density,ion_valence,savingstep,eps,runsname,grid);
  
  
  //Test
  
  RVF electric_field_a(nx,ny,nz);
  electric_field_a=RV(0,0,0);;
  for (int i=0;i<grid.nx();++i)
    for (int j=0;j<grid.ny();++j)
      for (int k=0;k<grid.nz();++k)
      {
        RV coord=grid.coordinates(i,j,k);
        Real ncoord=norm(coord);
        RV unit_vector_r=coord/ncoord;
        //cout << coord << endl;
        if(ncoord>ris&&ncoord<ros)
        {
          electric_field_a(i,j,k)=colloid_valence/(4.*M_PI*ncoord*ncoord)*unit_vector_r;
          //cout << electric_field_a(i,j,k) << endl;
        }
        else
        {
          electric_field_a(i,j,k)=RV(0,0,0);
        }
      }
  std::stringstream ssss;
  ssss << runsname << "_electric_field_analytic";
  vtkSave(ssss.str(),electric_field_a,"electric_field_analytic",grid);
  

//   RV avef(sum(electric_field[0]),sum(electric_field[1]),sum(electric_field[2]));
//   avef*=grid.deltav()/(lx*ly*lz);
//   electric_field-=avef;


  cout << electric_field(Range::all(),ny/2,nz/2)[0] << endl;
  cout << endl;
  cout << electric_field_a(Range::all(),ny/2,nz/2)[0] << endl;
  
  
  ofstream efcut("efcut.dat");
  
  for(int n=0;n<nx;++n)
    efcut << grid.coordinates(n,ny/2,nz/2)[0] << " " 
	<< electric_field(n,ny/2,nz/2)[0] << " " 
  	<< electric_field_a(n,ny/2,nz/2)[0] << endl;
  
cout << sum(electric_field_init[0](Range::all(),ny/2,nz/2)) << endl;
cout << sum(electric_field[0](Range::all(),ny/2,nz/2)) << endl;

  
  
  return 0;
}
