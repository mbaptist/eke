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

  
#if 0
  //Geometry
  
  Geometry geometry();
  
  Cube cube();
  Sphere sphere();
  
  geometry.add_domain(cube);
  geometry.add_outer_boundary(cube);
  
  
  //Grid
  Grid (geometry,shape);
  
  
#endif
  
  
#if 1
  
  
   //Grid
  Grid grid(nx,ny,nz,lx,ly,lz,rs);

  //Charges
  Array<double,3> charges(nx,ny,nz);
  //Distribute the charges over the surface of the colloid
  initialise_colloid(charges,grid,total_charge);
  
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

  //Test
  Array<TinyVector<double,3>,3> ef_test(nx,ny,nz);
  for(int i=0;i<nx;++i)
    for(int j=0;j<ny;++j)
      for(int k=0;k<nz;++k)
        if(grid.point_type()(i,j,k)==3)
	  {
	    //cout << "EF: " << electric_field(i,j,k) << endl;
	    //cout << "TT" << endl;
          ef_test(i,j,k)=grid.coordinates()(i,j,k);
          ef_test(i,j,k)*=.25/M_PI*total_charge/pow(norm(grid.coordinates()(i,j,k)),3);
	    ef_test(i,j,k)-=electric_field(i,j,k);
          //ef_test(i,j,k)/=norm(electric_field(i,j,k));
	    //cout << "EFT: " << ef_test(i,j,k) << endl;
	  }
	else
	  ef_test(i,j,k)=0;

  cout << "TEST: " << sum(dot(ef_test,ef_test))/(nx*ny*nz) << endl;
  cout << max(dot(ef_test,ef_test)) << endl;
  //cout << "TEST: " << sum(dot(electric_field,electric_field))/(nx*ny*nz) << endl;

#endif
  
  
  return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////

void initialise_colloid(Array<double,3> & charges,
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

void initialise_electric_field(Array<TinyVector<double,3>,3> & electric_field,
                               const Array<double,3> & charges,
                               const Grid & grid)
{
  int nx(grid.nx());
  int ny(grid.ny());
  int nz(grid.nz());
  electric_field=TinyVector<double,3>(0,0,0);
  Array<double,3> charge(charges.copy());
  double mean_charge;
  //Ex
  mean_charge=mean(charge(0,Range::all(),Range::all()));
  charge(0,Range::all(),Range::all())-=mean_charge;
  for(int j=0;j<ny;++j)
    for(int k=0;k<nz;++k)
      electric_field[0](0,j,k)=mean_charge;
  for(int i=1;i<nx-1;++i)
    {
      mean_charge=mean(charge(i,Range::all(),Range::all()));
      charge(i,Range::all(),Range::all())-=mean_charge;
      double efield=electric_field[0](i-1,0,0)+mean_charge;
       for(int j=0;j<ny;++j)
	for(int k=0;k<nz;++k)
	  electric_field[0](i,j,k)=efield;
    }
  mean_charge=mean(charge(nx-1,Range::all(),Range::all()));
  charge(nx-1,Range::all(),Range::all())-=mean_charge;
  for(int j=0;j<ny;++j)
    for(int k=0;k<nz;++k)
      electric_field[0](nx-1,j,k)=0.;
  //Ey
  for(int i=0;i<nx;++i)
    {
      mean_charge=mean(charge(i,0,Range::all()));
      charge(i,0,Range::all())-=mean_charge;
      for(int k=0;k<nz;++k)
	electric_field[1](i,0,k)=mean_charge;
      for(int j=1;j<ny-1;++j)
	{
	  mean_charge=mean(charge(i,j,Range::all()));
	  charge(i,j,Range::all())-=mean_charge;
	  double efield=electric_field[1](i,j-1,0)+mean_charge;
	  for(int k=0;k<nz;++k)
	    electric_field[1](i,j,k)=efield;
	}
      mean_charge=mean(charge(i,ny-1,Range::all()));
      charge(i,ny-1,Range::all())-=mean_charge;
      for(int k=0;k<nz;++k)
	electric_field[1](i,ny-1,k)=0.;
    }
  //Ez
  for(int i=0;i<nx;++i)
    for(int j=0;j<ny;++j)
      {
	electric_field[2](i,j,0)=charge(i,j,0);
	for(int k=1;k<nz-1;++k)
	  electric_field[2](i,j,k)=electric_field[2](i,j,k-1)+charge(i,j,k);
	electric_field[2](i,j,nz-1)=0;
      }

  cout << "Mean: " << mean(electric_field[0]) << endl;
  cout << "Mean: " << mean(electric_field[1]) << endl;
  cout << "Mean: " << mean(electric_field[2]) << endl;
  electric_field[0]-=mean(electric_field[0]);
  electric_field[1]-=mean(electric_field[1]);
  electric_field[2]-=mean(electric_field[2]);
  cout << "Mean: " << mean(electric_field[0]) << endl;
  cout << "Mean: " << mean(electric_field[1]) << endl;
  cout << "Mean: " << mean(electric_field[2]) << endl;
  for(int i=0;i<nx;++i)
    for(int j=0;j<ny;++j)
      for(int k=0;k<nz;++k)
        for(int mm=0;mm<3;++mm)
          if(electric_field(i,j,k)[mm]>1e-2)
            cout << electric_field(i,j,k)[mm] << endl;

  

#if 1
  //Test for div E = \rho
  for(int i=0;i<nx;++i)
    for(int j=0;j<ny;++j)
      for(int k=0;k<nz;++k)
	{
	  int n1,n2;
	  n1=i-1;
	  if (n1==-1)
	    n1=nx-1;
	  n2=i;
	  double divx=electric_field[0](n2,j,k)-electric_field[0](n1,j,k);
	  n1=j-1;
	  if (n1==-1)
	    n1=ny-1;
	  n2=j;
	  double divy=electric_field[1](i,n2,k)-electric_field[1](i,n1,k);
	  n1=k-1;
	  if (n1==-1)
	    n1=nz-1;
	  n2=k;
	  double divz=electric_field[2](i,j,n2)-electric_field[2](i,j,n1);
	  double div=divx+divy+divz;
	  if (fabs(div-charges(i,j,k))>.5e-10)
	    {
	      cout << i << " " << j << " " << k << " " 
		   << divx << " " << divy <<" " << divz << " " 
		   << div << " " << charges(i,j,k) << endl;
	      exit(1);
	    }
	}
#endif 
  
  
}

double functional(Array<TinyVector<double,3>,3> & electric_field)
{
  return double(.5*sum(dot(electric_field,electric_field)));
}



