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

using namespace std;

#include <blitz/blitz.h>
#include <blitz/array.h>
#include <blitz/tiny.h>
#include <blitz/tinyvec.h>

using namespace blitz;

BZ_DECLARE_FUNCTION_RET(norm,double)
BZ_DECLARE_FUNCTION2_RET(dot,double)


#include <blitz/random.h>
#include <random/uniform.h>
#include <random/discrete-uniform.h>

using namespace ranlib;

#include "eke.hpp"
#include "maggs.hpp"


///////////////////////////////////////////

int main()
{

  cout << fixed << setprecision(20);

  //Input (in future via python)
  int nx,ny,nz;
  double lx,ly,lz;
  double rs;
  double total_charge;
  nx=32;
  ny=nx;
  nz=ny;
  lx=32.;
  ly=32.;
  lz=32.;
  rs=6.;
  total_charge=1000000000000.;


   //Grid
  Grid grid(nx,ny,nz,lx,ly,lz,rs);

  //Charges	
  Array<double,3> charges(nx,ny,nz);
  //Distribute the charges
  distribute_charges(charges,grid,total_charge);
  
  //Save charges
  ofstream ofs("charges.dat");
  //ofs << charge_density << endl;
  for (int n=0;n<charges.size();++n)
    ofs << charges.data()[n] << endl;
   //cout << charge_density << endl;

  //electric field
  Array<TinyVector<double,3>,3> electric_field(nx,ny,nz);
  initialise_electric_field(electric_field,charges,grid);

  //Save initial electric field
  ofstream ofs2("efi.dat");
  //ofs2 << electric_field << endl;
  for (int n=0;n<electric_field.size();++n)
    {
      for (int m=0;m<3;++m)
	ofs2 << electric_field.data()[n][m] <<" ";
      ofs2 << endl;
    }

  //minimise
  int num_steps=0;
  double func0=functional(electric_field);
  while(1)
    {
      ++num_steps;
      cout << "Minimisation step " << num_steps << endl;
      sequential_sweep_loop_moves(electric_field,grid);
      sequential_sweep_charge_moves(electric_field,charges,grid);
      double func=functional(electric_field);      
      double delta_func=func-func0;
      cout << "Functional before // Functional after // Variation" << endl;
      cout <<  func0 << " // " << func << " // " << delta_func << endl << endl;
      if((delta_func<=0)&&(-delta_func<1e-16))
        break;
      func0=func;
    }
  
  //write electric field after minimisation  
  ofstream ofs3("ef.dat");
  ofs3 << fixed << setprecision(20);
  //ofs3 << electric_field << endl;
  for (int n=0;n<electric_field.size();++n)
    {
      for (int m=0;m<3;++m)
	ofs3 << electric_field.data()[n][m] <<" ";
      ofs3 << endl;
    }



 
  return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////

void distribute_charges(Array<double,3> & charges,
                        Grid & grid,
                        const double & total_charge)
{
  int nx(grid.nx());
  int ny(grid.ny());
  int nz(grid.nz());
  double lx(grid.lx());
  double ly(grid.ly());
  double lz(grid.lz());
  double rs(grid.rs());
  const blitz::Array<blitz::TinyVector<double,3>,3> & coordinates(grid.coordinates());
  blitz::Array<int,3> & point_type(grid.point_type());
  charges=0;
  int np=100000;
  double delta_charge=total_charge/np;
  UniformClosed<double> rg;
  //Distribute charges on the central sphere
  for(int n=0;n<np;++n)
  {
    double theta=rg.random()*M_PI;
    double phi=rg.random()*2.*M_PI;
    TinyVector<double,3> coord_ps(rs*cos(phi)*sin(theta),
                                  rs*sin(phi)*sin(theta),
                                  rs*cos(theta));
    int ind_lip_x=static_cast<int>((coord_ps[0]+.5*lx)/(lx/nx));
    int ind_lip_y=static_cast<int>((coord_ps[1]+.5*ly)/(ly/ny));
    int ind_lip_z=static_cast<int>((coord_ps[2]+.5*lz)/(lz/nz));
    TinyVector<double,3> dist_vec=coord_ps;
    dist_vec-=coordinates(ind_lip_x,ind_lip_y,ind_lip_z);
    double distance0=norm(dist_vec);
    int ind_np_x=ind_lip_x;
    int ind_np_y=ind_lip_y;
    int ind_np_z=ind_lip_z;
    for(int i=0;i<2;++i)
      for(int j=0;j<2;++j)
        for(int k=0;k<2;++k)
        {
          dist_vec=coord_ps;
          dist_vec-=coordinates(ind_lip_x+i,ind_lip_y+j,ind_lip_z+k);
          double distance=norm(dist_vec); 
          if(distance<distance0)
          {
            ind_np_x=ind_lip_x+i;
            ind_np_y=ind_lip_y+j;
            ind_np_z=ind_lip_z+k;
            distance0=distance;
          }
        }
    charges(ind_np_x,ind_np_y,ind_np_z)+=delta_charge;
    point_type(ind_np_x,ind_np_y,ind_np_z)=4;     
  }
  
  
for(int n=0;n<np;++n)
{
  DiscreteUniform<int> rgx(grid.nx());
  DiscreteUniform<int> rgy(grid.ny());
  DiscreteUniform<int> rgz(grid.nz());
  while(1)
  {
    int ni = rgx.random();
    int nj = rgy.random();
    int nk = rgz.random();
    if(grid.point_type()(ni,nj,nk)==1)
       {
       charges(ni,nj,nk)+=-delta_charge;
         break;
       }
  }
}
  
#if 0  
  
  //Put negative charges at the outer sphere
  for(int n=0;n<np;++n)
    {
      double theta=rg.random()*M_PI;
      double phi=rg.random()*2.*M_PI;
      TinyVector<double,3> coord_ps(2*rs*cos(phi)*sin(theta),
				    2*rs*sin(phi)*sin(theta),
				    2*rs*cos(theta));
      int ind_lip_x=static_cast<int>((coord_ps[0]+.5*lx)/(lx/nx));
      int ind_lip_y=static_cast<int>((coord_ps[1]+.5*ly)/(ly/ny));
      int ind_lip_z=static_cast<int>((coord_ps[2]+.5*lz)/(lz/nz));
      TinyVector<double,3> dist_vec=coord_ps;
      dist_vec-=coordinates(ind_lip_x,ind_lip_y,ind_lip_z);
      double distance0=norm(dist_vec);
      int ind_np_x=ind_lip_x;
      int ind_np_y=ind_lip_y;
      int ind_np_z=ind_lip_z;
      for(int i=0;i<2;++i)
	for(int j=0;j<2;++j)
	  for(int k=0;k<2;++k)
	    {
	      dist_vec=coord_ps;
	      dist_vec-=coordinates(ind_lip_x+i,ind_lip_y+j,ind_lip_z+k);
	      double distance=norm(dist_vec); 
	      if(distance<distance0)
		{
		  ind_np_x=ind_lip_x+i;
		  ind_np_y=ind_lip_y+j;
		  ind_np_z=ind_lip_z+k;
		  distance0=distance;
		}
	    }

      charges(ind_np_x,ind_np_y,ind_np_z)+=-delta_charge;
      point_type(ind_np_x,ind_np_y,ind_np_z)=2;     
    }     
#endif

}







