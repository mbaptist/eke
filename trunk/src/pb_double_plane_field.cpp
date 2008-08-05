// -*- C++ -*-

//
// C++ Implementation: pb_double_plane.cpp
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
    run_name="pb_double_plane.cfg";
  
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
  double xleft=input.parse_double("xleft");
  double xright=input.parse_double("xright");
  //Double plane
  double lb=input.parse_double("lb");
  double numberoffixedcharges=input.parse_double("numberoffixedcharges");
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
  double k=sqrt(4.*M_PI*lb*fabs(numberoffixedcharges)/(lx*ly*lz));
  numberoffixedcharges*=4.*M_PI*lb*k;
  for (int n=0;n<ion_number.size();++n)
    ion_number[n]*=4.*M_PI*lb*k;
  lx*=k;
  ly*=k;
  lz*=k;
  xleft*=k;
  xright*=k;
  
  //Grid
  Grid grid(nx,ny,nz,lx,ly,lz);
  
  //Classify Grid Points
  for (int i=0;i<grid.nx();++i)
    for (int j=0;j<grid.ny();++j)
      for (int k=0;k<grid.nz();++k)
      {
        double x=grid.coordinates(i,j,k)[0];
        if(x<xleft||x>xright)
          grid.point_type()(i,j,k)=0;
        else if(x>xleft&&x<xright)
          grid.point_type()(i,j,k)=1;
        else
          grid.point_type()(i,j,k)=2;
      }
  
  //Fixed charges on the planes	
  RSF fixed_charge_density(nx,ny,nz);
  fixed_charge_density=0;
  //Determine the x indices of the planes
  int ileft=grid.nearest_point_index(RV(xleft,0,0))[0];
  int iright=grid.nearest_point_index(RV(xright,0,0))[0];
  //Number of points in both planes
  int np=grid.ny()*grid.nz();
  Real delta_density=(numberoffixedcharges/grid.deltav())/np;
  for (int j=0;j<grid.ny();++j)
    for (int k=0;k<grid.nz();++k)
    {
      fixed_charge_density(ileft,j,k)=delta_density;
      fixed_charge_density(iright,j,k)=-delta_density;
    }
  //Save the fixed charge density
  std::stringstream ss;
  ss << runsname << "_fixed_charge_density";
  vtkSave(ss.str(),fixed_charge_density,"fixed_charge_density",grid);
  
  //Distribute ionic species	
  std::vector<RSF> ion_density;
  //distribute_ionic_species(ion_density,ion_number,runsname,grid);
  ion_density.push_back(RSF(nx,ny,nz));
  ion_density[0]=0;
  
  //electric field
  RVF electric_field(nx,ny,nz);
  initialise_electric_field(electric_field,fixed_charge_density,ion_density,ion_valence,eps,runsname,grid);
 
#if 1
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
    //Field moves
    Real delta_func=sequential_sweep_loop_moves(electric_field,grid);
      //Update the functional
    func+=delta_func;
      //Print iteration infos
    cout << " Minimisation step: " << num_steps << "\t"
      << "Variation in functional: " << delta_func << "\t"
      << "Functional: " << func << endl;
    // cout << sum(divergence(electric_field,grid)) << endl;
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
  }
    //Save final values of the field and concentrations
  std::stringstream sss;
  sss << runsname << "_electric_field_" << num_steps;
  vtkSave(sss.str(),electric_field,ss.str(),grid);
  cout << "...done." << endl;
#endif
  
    //Test
  
  RVF electric_field_a(nx,ny,nz);
  electric_field_a=RV(0,0,0);;
  for (int i=0;i<grid.nx();++i)
    for (int j=0;j<grid.ny();++j)
      for (int k=0;k<grid.nz();++k)
      {
        Real x=grid.coordinates(i,j,k)[0];
        if(x>xleft&&x<xright)
        {
          electric_field_a(i,j,k)=numberoffixedcharges/(lx*ly)*unit_vector(0);
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
  
  cout << electric_field(Range::all(),ny/2,nz/2)[0] << endl;
  cout << endl;
  cout << electric_field_a(Range::all(),ny/2,nz/2)[0] << endl;
  
  ofstream efcut("efcut.dat");
  
  for(int n=0;n<electric_field(Range::all(),ny/2,nz/2)[0].size();++n)
    efcut << electric_field(Range::all(),ny/2,nz/2)[0](n) << " " 
    << electric_field_a(Range::all(),ny/2,nz/2)[0](n) << endl;
  
  RVF electric_field_comp(electric_field_a.copy());
  electric_field_comp-=electric_field;
  cout << sum(dot(electric_field_comp,electric_field_comp))*grid.deltav() << endl;
  
  return 0;
}

