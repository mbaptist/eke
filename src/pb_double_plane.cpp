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
  int np=2*grid.ny()*grid.nz();
  Real delta_density=(numberoffixedcharges/grid.deltav())/np;
  for (int j=0;j<grid.ny();++j)
    for (int k=0;k<grid.nz();++k)
    {
      fixed_charge_density(ileft,j,k)=delta_density;
      fixed_charge_density(iright,j,k)=delta_density;
    }
  //Save the fixed charge density
  std::stringstream ss;
  ss << runsname << "_fixed_charge_density";
  vtkSave(ss.str(),fixed_charge_density,"fixed_charge_density",grid);
  
  //Distribute ionic species	
  std::vector<RSF> ion_density;
  distribute_ionic_species(ion_density,ion_number,runsname,grid);
  
  //electric field
  RVF electric_field(nx,ny,nz);
  initialise_electric_field(electric_field,fixed_charge_density,ion_density,ion_valence,eps,runsname,grid);
  
  //Minimise
  minimise(electric_field,ion_density,ion_valence,savingstep,eps,runsname,grid);
  
  blitz::Array<Real,1> c_eke(nx);
  c_eke=ion_density[0](Range::all(),ny/2,nz/2);
  blitz::Array<Real,1> c_a(nx);
  blitz::Array<Real,1> c_a_av(nx);
  
  Real sigma=numberoffixedcharges/2./(ly*lz);
  cout << "Sigma: " << sigma << endl;
  Real s;
  Real s0=1.;
  while(1)
  {
    s=atan(-.5*sigma*ion_valence[0]*xright/s0);
    //cout << s << endl;
    if(fabs(s-s0)<eps)
      break;
    s0=s;
  }
  //cout << s << endl;

  std::stringstream sss;
  sss << runsname << "_analytic.dat";
  ofstream ofs(sss.str().c_str());
  
  //Maximum relative error
  Real mre;

  for(int i=0;i<nx;++i)
  {
    Real x=-.5*lx+grid.deltax()*i;
    if(x>xleft&&x<xright)
    {
      c_a(i)=2.*pow(s/(ion_valence[0]*xright),2)/pow(cos(s*x/xright),2);
	Real re=fabs(c_eke(i)-c_a(i))/c_a(i);
	if(re>mre)
		mre=re;
      ofs << x << " " << c_a(i) << " " << c_eke(i) << " " << re << endl;
    }
    else
    {
      c_a(i)=0;
	ofs << x << " " << c_a(i) << " " << c_eke(i) << " " << 0 << endl;
    }
  }

cout << "Maximum relative error: " << mre << endl;  

  return 0;
}

