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
  UniformClosed<Real> rg;
  //number of points to represent the sphere surface
  //100 per grid unit area
  int np=100*static_cast<int>(4*M_PI*pow(ros,2)/grid.deltasx());
  //cout << np << endl;
  Real delta_density=(colloid_valence/grid.deltav())/np;
  for(int n=0;n<np;++n)
  {
    Real theta=rg.random()*M_PI;
    Real phi=rg.random()*2.*M_PI;
    //Inner sphere
    RV coord_ps(ris*cos(phi)*sin(theta),
                ris*sin(phi)*sin(theta),
                ris*cos(theta));
    IV ind_np=grid.nearest_point_index(coord_ps);	
    fixed_charge_density(ind_np[0],ind_np[1],ind_np[2])+=delta_density;
    //Outer sphere
    coord_ps=RV(ros*cos(phi)*sin(theta),
                ros*sin(phi)*sin(theta),
                ros*cos(theta));
    ind_np=grid.nearest_point_index(coord_ps);	
    fixed_charge_density(ind_np[0],ind_np[1],ind_np[2])-=delta_density;
  }
  //Save the fixed charge density
  std::stringstream ss;
  ss << runsname << "_fixed_charge_density";
  vtkSave(ss.str(),fixed_charge_density,"fixed_charge_density",grid);
  
  //Distribute ionic species	
  //There are no movable charges (ionic species)
  //Setting ionic density to zero in order
  //to be able to use the electric field
  //initialisation function
  std::vector<RSF> ion_density;
  ion_density.push_back(RSF(nx,ny,nz));
  ion_density[0]=0;
  
  //electric field
  RVF electric_field(nx,ny,nz);
  initialise_electric_field(electric_field,fixed_charge_density,ion_density,ion_valence,eps,runsname,grid);
  
  //Minimise
  //minimise(electric_field,ion_density,ion_valence,100,eps,grid);
  //Minimise only for the electric field
  cout << "Minimising... " << endl;
  int num_steps=0;
  Real func0=.5*sum(dot(electric_field,electric_field));
  cout << " Initial value of the functional: " << func0 << endl;
  Real func=func0;
  while(1)
  {
    ++num_steps;
    Real delta_func=0;
      //Field moves
    delta_func+=sequential_sweep_loop_moves(electric_field,grid);
      //Update the functional
    func+=delta_func;
      //Print iteration infos
    cout << " Minimisation step: " << num_steps << "\t"
      << "Variation in functional: " << delta_func << "\t"
      << "Functional: " << func << endl;
  //Save field
    if(num_steps%savingstep==0)
    {
      std::stringstream ss;
      ss << runsname << "_electric_field_" << num_steps;
      vtkSave(ss.str(),electric_field,ss.str(),grid);
    }
      //Check for stop criterium
    if(fabs(delta_func)<eps)
      break;
    //Save final values of the field and concentrations
    std::stringstream ss;
    ss << runsname << "_electric_field_" << num_steps;
    vtkSave(ss.str(),electric_field,ss.str(),grid);
  }
  cout << "...done." << endl;
    
  return 0;
}
