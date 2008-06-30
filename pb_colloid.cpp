// -*- C++ -*-

//
// C++ Implementation: pb_colloid
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
  string cfg_file;
  if(argc>1)
    cfg_file=argv[1];
  else
    cfg_file="pb_colloid.cfg";
  
  cout << fixed << setprecision(16);
  
  //Input
  //Open input file
  PyInputParser input(cfg_file.c_str());
  //Runsname
  string runsname=input.parse_string("runsname");
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
  rs*=k; 
  
  //Grid
  Grid grid(nx,ny,nz,lx,ly,lz);
  
  //Classify Grid Points
  for (int i=0;i<grid.nx();++i)
    for (int j=0;j<grid.ny();++j)
      for (int k=0;k<grid.nz();++k)
        if(norm(grid.coordinates(i,j,k))>rs)
          grid.point_type()(i,j,k)=1;
  else
    grid.point_type()(i,j,k)=0;
  
  //Fixed charges (colloidal charges)
  RSF colloid_charge_density(nx,ny,nz);
  colloid_charge_density=0;
#if 1
  UniformClosedOpen<Real> rgco;
  UniformClosed<Real> rgc;
  //number of points to represent the sphere surface
  //100 per grid unit area
  int np=100*static_cast<int>(4*M_PI*pow(rs,2)/grid.deltasx());
  Real delta_density=(colloid_valence/grid.deltav())/np;
  for(int n=0;n<np;++n)
  {
    Real z=-1.+rgc.random()*2.;
    Real phi=rgco.random()*2.*M_PI;
    Real rsintheta=sqrt(rs*rs-z*z);
    RV coord_ps(rsintheta*cos(phi),
                rsintheta*sin(phi),
                z);
    IV ind_np=grid.nearest_point_index(coord_ps);	
    colloid_charge_density(ind_np[0],ind_np[1],ind_np[2])+=delta_density;
    grid.point_type()(ind_np[0],ind_np[1],ind_np[2])=2;
  } 
#endif
# if 0
int NN=1000;

Real delta_density=(colloid_valence/grid.deltav());
int np=0;
 for(int nn=0;nn<NN;++nn)
{
Real theta=nn*M_PI/(NN-1);
int MM=static_cast<int>(2.*NN*sin(theta));
for(int mm=0;mm<MM;++mm)
{
Real phi=mm*2.*M_PI/MM;
	RV coord_ps(rs*cos(phi)*sin(theta),
                rs*sin(phi)*sin(theta),
                rs*cos(theta));
    IV ind_np=grid.nearest_point_index(coord_ps);
    colloid_charge_density(ind_np[0],ind_np[1],ind_np[2])+=delta_density;
    grid.point_type()(ind_np[0],ind_np[1],ind_np[2])=2;
++np;
  }
}
colloid_charge_density/=np;
#endif

  //Save the fixed charge density
  std::stringstream ss;
  ss << runsname << "_fixed_charge_density";
  vtkSave(ss.str(),colloid_charge_density,"fixed_charge_density",grid);
  
//Distribute ionic species	
  std::vector<RSF> ion_density;
  distribute_ionic_species(ion_density,ion_number,runsname,grid);
  
  //electric field
  RVF electric_field(nx,ny,nz);
  initialise_electric_field(electric_field,colloid_charge_density,ion_density,ion_valence,eps,runsname,grid);
  
  //Minimise
  minimise(electric_field,ion_density,ion_valence,savingstep,eps,runsname,grid);
  
  return 0;
}
