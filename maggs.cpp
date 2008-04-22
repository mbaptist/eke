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

#include "maggs.hpp"
#include "types.hpp"
#include "random.hpp"

using namespace ranlib;



//// BASIC MOVES ////


void concentration_move(RVF & electric_field,
                        RSF & concentration,
                        const Grid & grid,
                        const IV & node1, const int & dir)
{
	IV node2(node1);
	node2[dir]+=1;
	if ((grid.point_type()(node1[0],node1[1],node1[2])==1)
	    &&(grid.point_type()(node2[0],node2[1],node2[2])==1))
	{
		Real & c1 = concentration(node1[0],node1[1],node1[2]);
		Real & c2 = concentration(node2[0],node2[1],node2[2]);
		Real & e = electric_field(node1[0],node1[1],node1[2])[dir];
		//Optimal concentration change
		//solve for the optimal concentration change with fixed point iterations
		Real deltac;
		Real deltac0=c1;
		Real area=grid.deltas(dir);
		while(1)
		{
			Real earg=(-e+deltac0/area)/area;
			//cout << earg << endl;
			if (earg<0)
				deltac=(c1-c2*exp(earg))/(1+exp(earg));
			else
				deltac=(c1*exp(-earg)-c2)/(exp(-earg)+1);
			if(fabs(deltac-deltac0)<.5e-16)
				break;
			deltac0=deltac;
		}
		//Evaluate the change in the functional
		Real delta_func=-e*deltac/area+.5*pow(deltac/area,2)
			+c1*log(1-deltac/c1)+c2*log(1+deltac/c2)
			-deltac*log(c1-deltac)+deltac*log(c2+deltac);
		//Accept the move if if minimises the functional
		if (delta_func<0)
		{
			e-=deltac/area;
			c1-=deltac;
			c2+=deltac;
		}
	}
}


void loop_move(RVF & electric_field,
               const Grid & grid,
               const Loop & loop)
{
	//Evaluate the trial field
	Real & e1 = electric_field(loop.node1()[0],loop.node1()[1],loop.node1()[2])[loop.dir1()];
	Real & e2 = electric_field(loop.node2()[0],loop.node2()[1],loop.node2()[2])[loop.dir2()];
	Real & e3 = electric_field(loop.node3()[0],loop.node3()[1],loop.node3()[2])[loop.dir1()];
	Real & e4 = electric_field(loop.node4()[0],loop.node4()[1],loop.node4()[2])[loop.dir2()];
	//Optimal move
	//This is valid both for electrostatics and Poisson-Boltzmann
	//It is characteristic of loop moves
	Real delta_field=-.25*(e1+e2-e3-e4);
	Real e1p=e1+delta_field;
	Real e2p=e2+delta_field;
	Real e3p=e3-delta_field;
	Real e4p=e4-delta_field;
	//Evaluate the change in the functional
	//This is valid both for electrostatics and Poisson-Boltzmann
	//It is characteristic of loop moves
	Real delta_func=(e1p*e1p+e2p*e2p+e3p*e3p+e4p*e4p
	                   -e1*e1-e2*e2-e3*e3-e4*e4);
	//cout << "Dfunc=" << delta_func << endl;
	//Accept/reject move
	if (delta_func<0)
	{
		e1=e1p;
		e2=e2p;
		e3=e3p;
		e4=e4p;
	}
	//cout << "Circulation: " << e1+e2-e3-e4 << endl;
};



//// LATTICE SWEEPS ////


void sequential_sweep_concentration_moves(RVF & electric_field,
                                          RSF & concentrations,
                                          const Grid & grid)
{
	Real delta_func;
	for (int i=0;i<grid.nx();++i)
		for (int j=0;j<grid.ny();++j)
			for (int k=0;k<grid.nz();++k)
				for (int dir=0;dir<3;++dir)
				{
//Perform the concentration move
					concentration_move(electric_field,concentrations,grid,TinyVector<int,3>(i,j,k),dir);
				}
}


void sequential_sweep_loop_moves(RVF & electric_field,
                                 const Grid & grid)
{
	Real delta_func;
	for (int i=0;i<grid.nx();++i)
		for (int j=0;j<grid.ny();++j)
			for (int k=0;k<grid.nz();++k)
				for (int loop_number=0;loop_number<3;++loop_number)
				{
//Create the loop
					Loop loop(blitz::TinyVector<int,3>(i,j,k),loop_number,grid);
//Perform the loop move
					loop_move(electric_field,grid,loop);
				}
}

void random_sweep_loop_moves(RVF & electric_field,
                             const Grid & grid)
{
  //Choose 3x the number of nodes
	for (int nn=0;nn<3*electric_field.size();++nn)
	{
      //Choose a node
		DiscreteUniform<int> rgx(grid.nx());
		DiscreteUniform<int> rgy(grid.ny());
		DiscreteUniform<int> rgz(grid.nz());
		int ni = rgx.random();
		int nj = rgy.random();
		int nk = rgz.random();
      //cout << "Node: " << ni << " " << nj << " " << nk << endl;
      //Choose a loop
		DiscreteUniform<int> rgp(3);
		int loop_number=rgp.random();
      //cout << "Loop number: " << loop_number << endl;
      //Create the loop
		Loop loop(blitz::TinyVector<int,3>(ni,nj,nk),loop_number,grid);
      //Perform the loop move
		loop_move(electric_field,grid,loop);
	}
}


//// INITIALISATIONS ////

void initialise_electric_field(RVF & electric_field,
                               const RSF & total_density,
                               const Grid & grid)
{
	int nx(grid.nx());
	int ny(grid.ny());
	int nz(grid.nz());
	electric_field=RV(0,0,0);
	RSF density(total_density.copy());
	Real mean_density;
	double dv=grid.deltav();
//Ex
	Real ds=grid.deltasx();
	mean_density=mean(density(0,Range::all(),Range::all()));
	density(0,Range::all(),Range::all())-=mean_density;
	for(int j=0;j<ny;++j)
		for(int k=0;k<nz;++k)
			electric_field[0](0,j,k)=mean_density*dv/ds;
	for(int i=1;i<nx-1;++i)
	{
		mean_density=mean(density(i,Range::all(),Range::all()));
		density(i,Range::all(),Range::all())-=mean_density;
		for(int j=0;j<ny;++j)
			for(int k=0;k<nz;++k)
				electric_field[0](i,j,k)=electric_field[0](i-1,j,k)+mean_density*dv/ds;;
	}
	mean_density=mean(density(nx-1,Range::all(),Range::all()));
	density(nx-1,Range::all(),Range::all())-=mean_density;
	for(int j=0;j<ny;++j)
		for(int k=0;k<nz;++k)
			electric_field[0](nx-1,j,k)=0.;
  //Ey
	ds=grid.deltasy();
	for(int i=0;i<nx;++i)
	{
		mean_density=mean(density(i,0,Range::all()));
		density(i,0,Range::all())-=mean_density;
		for(int k=0;k<nz;++k)
			electric_field[1](i,0,k)=mean_density*dv/ds;
		for(int j=1;j<ny-1;++j)
		{
			mean_density=mean(density(i,j,Range::all()));
			density(i,j,Range::all())-=mean_density;
			for(int k=0;k<nz;++k)
				electric_field[1](i,j,k)=electric_field[1](i,j-1,k)+mean_density*dv/ds;;
		}
		mean_density=mean(density(i,ny-1,Range::all()));
		density(i,ny-1,Range::all())-=mean_density;
		for(int k=0;k<nz;++k)
			electric_field[1](i,ny-1,k)=0.;
	}
  //Ez
	ds=grid.deltasz();
	for(int i=0;i<nx;++i)
		for(int j=0;j<ny;++j)
		{
			electric_field[2](i,j,0)=density(i,j,0)*dv/ds;
			for(int k=1;k<nz-1;++k)
				electric_field[2](i,j,k)=electric_field[2](i,j,k-1)+density(i,j,k)*dv/ds;
			electric_field[2](i,j,nz-1)=0;
		}
	
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
				Real divx=electric_field[0](n2,j,k)-electric_field[0](n1,j,k);
				divx/=grid.deltax();
				n1=j-1;
				if (n1==-1)
					n1=ny-1;
				n2=j;
				Real divy=electric_field[1](i,n2,k)-electric_field[1](i,n1,k);
				divy/=grid.deltay();
				n1=k-1;
				if (n1==-1)
					n1=nz-1;
				n2=k;
				Real divz=electric_field[2](i,j,n2)-electric_field[2](i,j,n1);
				divz/=grid.deltaz();
				Real div=divx+divy+divz;
				if (fabs(div-total_density(i,j,k))>.5e-10)
				{
					cout << i << " " << j << " " << k << " " 
						<< divx << " " << divy <<" " << divz << "\n " 
						<< div << " " << total_density(i,j,k) << endl;
					exit(1);
				}
			}
#endif 
	
	
}


//// FUNCTIONALS ////

//Functional for electrostitcs
Real functional(const RVF & electric_field)
{
	return Real(.5*sum(dot(electric_field,electric_field)));
}

//Functional for Poisson-Boltzmann
Real functional(const RSF & concentration,
                  const RVF & electric_field)
{
	Real result=.5*sum(dot(electric_field,electric_field));
	for(int n=0;n<concentration.size();++n)
	{
		Real c=concentration.data()[n];
		if (c!=0)
			result+=-c*log(-c);
	}
	//return Real(.5*sum(dot(electric_field,electric_field))+sum(concentration*log(concentration)));
}

