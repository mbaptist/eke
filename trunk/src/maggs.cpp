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
#include "io.hpp"

using namespace ranlib;


//// INITIALISATIONS ////

//Distribute ionic species
void distribute_ionic_species(std::vector<RSF> & ion_density,
                              const std::vector<double> & ion_number,
                              const string & runsname,
                              const Grid & grid)
{
  cout << "Distributing ions species..." << endl;
  int nx(grid.nx());
  int ny(grid.ny());
  int nz(grid.nz());
  int npt1=count(grid.point_type()==1);
  cout << " number of point of type 1 (where charges can move): " <<  npt1 << endl;
  for (int n=0;n<ion_number.size();++n)
  {
    cout << " distributing ions of species " << n << ":" << endl;
    ion_density.push_back(RSF(nx,ny,nz));
    ion_density[n]=0;
    Real delta_density=(ion_number[n]/grid.deltav())/npt1;
    cout << "  delta_density: " << delta_density << endl;
    cout << "  running nodes..." << endl;
    for (int i=0;i<grid.nx();++i)
      for (int j=0;j<grid.ny();++j)
        for (int k=0;k<grid.nz();++k)
          if(grid.point_type()(i,j,k)==1)
            ion_density[n](i,j,k)=delta_density;  
    cout << "  saving density..." << endl;
    std::stringstream ss;
    ss << runsname << "_ion_" << n << "_density_0";
    vtkSave(ss.str(),ion_density[n],ss.str(),grid);
  }
  cout << "... done." << endl;
}


//Initialise the elctric field
void initialise_electric_field(RVF & electric_field,
                               const RSF & fixed_charge_density,
                               const std::vector<RSF> & ion_density,
                               const std::vector<int> & ion_valence,
                               const Real & eps,
                               const string & runsname,
                               const Grid & grid)
{
  cout << "Initialising the electric field..." << endl;
  
  //Initilise the field with zeros
  electric_field=RV(0,0,0);
  //Evaluate the total charge density (tcd)
  RSF tcd(fixed_charge_density.copy());
  for (int n=0;n<ion_density.size();++n)
    tcd+=RSF(ion_valence[n]*ion_density[n]);
  
  RSF density(tcd.copy());
  //Check charge neutrality
  cout << " Checking charge neutrality..." << endl;
  Real total_charge=sum(tcd)*grid.deltav();
  cout << "  total charge: " << total_charge << endl;
#if 1
  if (fabs(total_charge)>eps)
  {
    cout << "  The total charge must be zero, instead of "
      << total_charge << ".\n"
      << "  Aborting..." << endl;
    exit(1);
  }
#endif


  //Get some necessary grid parameters
  int nx(grid.nx());
  int ny(grid.ny());
  int nz(grid.nz());
#if 0
//Initialize Ez
Real mean_density=mean(density(Range::all(),Range::all(),0));
density(Range::all(),Range::all(),0)-=mean_density;
for(int i=0;i<nx;++i)
    for(int j=0;j<ny;++j)
      electric_field[2](i,j,0)=0;
  for(int k=1;k<nz;++k)
  {
    mean_density=mean(density(Range::all(),Range::all(),k));
    density(Range::all(),Range::all(),k)-=mean_density;
for(int i=0;i<nx;++i)
    for(int j=0;j<ny;++j)
        electric_field[2](i,j,k)=electric_field[2](i,j,k-1)+mean_density;
  }
  //geometric factors for Ez
  electric_field[2]*=grid.deltal(2);
  //Initialise Ey
  for(int k=0;k<nz;++k)
  {
    mean_density=mean(density(Range::all(),0,k));
    density(Range::all(),0,k)-=mean_density;
    for(int i=0;i<nx;++i)
      electric_field[1](i,0,k)=0;
    for(int j=1;j<ny;++j)
    {
      mean_density=mean(density(Range::all(),j,k));
      density(Range::all(),j,k)-=mean_density;
      for(int i=0;i<nx;++i)
        electric_field[1](i,j,k)=electric_field[1](i,j-1,k)+mean_density;
    }
  }
    //geometric factors for Ey
  electric_field[1]*=grid.deltal(1);
  //Initialise Ex
for(int j=0;j<ny;++j)
  for(int k=0;k<nz;++k)
    {
      electric_field[0](0,j,k)=0;
      for(int i=1;i<nx;++i)
        electric_field[0](i,j,k)=electric_field[0](i-1,j,k)+density(i,j,k);
    }  
  //geometric factors for Ex
  electric_field[0]*=grid.deltal(0);
#endif

#if 1
  //Initialise Ex
  Real mean_density=mean(density(0,Range::all(),Range::all()));
  density(0,Range::all(),Range::all())-=mean_density;
  
  for(int j=0;j<ny;++j)
    for(int k=0;k<nz;++k)
      electric_field[0](0,j,k)=0;
  for(int i=1;i<nx;++i)
  {
    mean_density=mean(density(i,Range::all(),Range::all()));
    density(i,Range::all(),Range::all())-=mean_density;
    for(int j=0;j<ny;++j)
      for(int k=0;k<nz;++k)
        electric_field[0](i,j,k)=electric_field[0](i-1,j,k)+mean_density;
  }
  //geometric factors for Ex
  electric_field[0]*=grid.deltal(0);
  //Initialise Ey
  for(int i=0;i<nx;++i)
  {
    mean_density=mean(density(i,0,Range::all()));
    density(i,0,Range::all())-=mean_density;
    for(int k=0;k<nz;++k)
      electric_field[1](i,0,k)=0;
    for(int j=1;j<ny;++j)
    {
      mean_density=mean(density(i,j,Range::all()));
      density(i,j,Range::all())-=mean_density;
      for(int k=0;k<nz;++k)
        electric_field[1](i,j,k)=electric_field[1](i,j-1,k)+mean_density;
    }
  }
    //geometric factors for Ey
  electric_field[1]*=grid.deltal(1);
  //Initialise Ez
  for(int i=0;i<nx;++i)
    for(int j=0;j<ny;++j)
    {
      electric_field[2](i,j,0)=0;
      for(int k=1;k<nz;++k)
        electric_field[2](i,j,k)=electric_field[2](i,j,k-1)+density(i,j,k);
    }  
  //geometric factors for Ez
  electric_field[2]*=grid.deltal(2);
#endif
	//Remove the average
	RV avef(sum(electric_field[0]),sum(electric_field[1]),sum(electric_field[2]));
  	avef*=grid.deltav()/(grid.lx()*grid.ly()*grid.lz());
  	electric_field-=avef;


  //Test for div E = \rho
  cout << " Checking Gauss law..." << endl;
  if(count(fabs(divergence(electric_field,grid)-tcd)>eps)>0)
  {
    cout << count(fabs(divergence(electric_field,grid)-tcd)>eps) << endl;
    cout << sum(divergence(electric_field,grid)) << endl;
    cout << sum(tcd) << endl;
    cout << "  The electric field doesn't obey the Gauss law. \n"
      << "  Aborting..." << endl;
    exit(1);
  }
  cout << " Saving..." << endl;
  std::stringstream ss;
  ss << runsname << "_electric_field_0";
  vtkSave(ss.str(),electric_field,ss.str(),grid);
  cout << "...done." << endl;
}


//// FUNCTIONAL ////

Real functional(const RVF & electric_field,const std::vector<RSF> ion_concentration,const Grid & grid)
{
  Real func=.5*sum(dot(electric_field,electric_field));
  for (int n=0;n<ion_concentration.size();++n)
    func+=sum(where(ion_concentration[n]>0,ion_concentration[n]*log(ion_concentration[n]),0));
  func*=grid.deltav();
  return func;
}

///////////////////////


//// LOOP MOVES ////

//Single loop move
inline Real loop_move(RVF & electric_field,const Loop & loop)
{
  //References to the fields
  Real & e1 = electric_field(loop.node1())[loop.dir1()];
  Real & e2 = electric_field(loop.node2())[loop.dir2()];
  Real & e3 = electric_field(loop.node4())[loop.dir1()];
  Real & e4 = electric_field(loop.node1())[loop.dir2()];
  //Optimal change
  Real delta_field=-.25*(e1+e2-e3-e4);
  //Evaluate the change in the functional
  Real delta_func=delta_field*(e1+e2-e3-e4+2.*delta_field);

  //Return the change in the functional
if(delta_func<0)
{
  //Update the field
  e1+=delta_field;
  e2+=delta_field;
  e3-=delta_field;
  e4-=delta_field;


  return delta_func;
}
else
return 0;
}

//Lattice sweep
Real sequential_sweep_loop_moves(RVF & electric_field,const Grid & grid)
{
  Real delta_func=0;
  for (int i=0;i<grid.nx();++i)
    for (int j=0;j<grid.ny();++j)
      for (int k=0;k<grid.nz();++k)
        for (int loop_number=0;loop_number<3;++loop_number)
        {
          //Create the loop
          Loop loop(IV(i,j,k),loop_number,grid);
          //Perform the loop move
          delta_func+=loop_move(electric_field,loop);
        }
  return delta_func;
}

//////////////////////////////


//// ION MOVES ////

// ION MOVE ALONG SOME PATH //
inline Real ion_move(RVF & electric_field,
              RSF & concentration,
              const int & ion_valence,
              const IV & node1,
              const int & dir,
              const Real & eps,
              const Grid & grid)
{
  //Note that node1 must be a valid node
  //Define node2
  IV node2(grid.intobox(node1+unit_vector(dir)));
  //Check if node 2 belongs to region 1
  if (grid.point_type()(node2)!=1)
    return 0;
  //Define references to the field and concentrations
  Real & c1 = concentration(node1);
  Real & c2 = concentration(node2);
  Real & e = electric_field(node1)[dir];
//Check if concentration is negative
  if((c1<0)||(c2<0))
  {
    cout << "Concentrations cannot be negative. Aborting..." << endl;
    exit(1);
  }
  //Define concentration change
  Real deltac;
    //Determine the optimal concentration change
#if 0
    //solve for the optimal concentration change with successive bisections
  Real a=-c2;
  Real fa=d_deltafunc_d_deltac(a,c1,c2,e,-ion_valence*grid.deltal(dir));
  Real b=c1;
  Real fb=d_deltafunc_d_deltac(b,c1,c2,e,-ion_valence*grid.deltal(dir));
    //cout << fa << " " << fb << endl;
  if(fa*fb<0)
  {
    while(1)
    {
      Real c=.5*(a+b);
      Real fc=d_deltafunc_d_deltac(c,c1,c2,e,-ion_valence*grid.deltal(dir));
        //cout << a << " " << b << " " << b-a <<" " << c << " " << fc << endl;
      if((fc==0)||(b-a<eps))
      {
        deltac=c;
        break;
      }
      else if(fa*fc<0)
      {
        b=c;
        fb=fc;
      }
      else
      {
        a=c;
        fa=fc;
      }
    }   
  }
  else
  {
    cout << "Successive bisection method failed" << endl;
    exit(1);
  }
#endif
#if 0
  //Secant method
  Real deltac0=-c2;
  Real f0=d_deltafunc_d_deltac(deltac0,c1,c2,e,-ion_valence*grid.deltal(dir));
  Real deltac1=c1;
  while(1)
  {
    Real f1=d_deltafunc_d_deltac(deltac1,c1,c2,e,-ion_valence*grid.deltal(dir));
    deltac=deltac1-(deltac1-deltac0)*f1/(f1-f0);
    if(fabs(deltac-deltac1)<eps*1e-2)
      break;
    //cout << deltac << endl;
    deltac0=deltac1;
    f0=f1;
    deltac1=deltac;
  } 
#endif
#if 1
	//Newton's Method
	Real ivdl=ion_valence*grid.deltal(dir);
	Real cp=.5*(c1+c2);
	Real cm=.5*(c1-c2);	
	Real a=2./(pow(ivdl,2)*cp);
	Real b=e/(ivdl*cp)-cm/cp;
	Real csi0=-.5*ivdl*(e-ivdl*cm);
	Real csi;
	while(1)
	{
		Real th=tanh(csi0);
		csi=csi0-(th+a*csi0+b)/(1.-pow(th,2)+a);
		if(fabs(csi-csi0)<2./pow(ivdl,2)*eps*1e-2)
      		break;	
    //cout << deltac << endl;
    csi0=csi;
	}
	deltac=(e+2*csi/ivdl)/ivdl;
	
#endif
	
  //Evaluate the change in the electric field
  Real deltae=-ion_valence*deltac*grid.deltal(dir);
  //Evaluate the change in the functional
  Real delta_func=deltafunc(deltac,deltae,c1,c2,e);
  //cout << "delta_func: " << delta_func << endl;

  //Return the change in the functional
if (delta_func<0)
{
  //Update the concentrations
  c1-=deltac;
  c2+=deltac;
  //Update the field
  e+=deltae;

  return delta_func;

}
else 
	return 0;
}

//Variation of the functional for ion moves
inline Real deltafunc(const Real & deltac,const Real & deltae,const Real & c1, const Real & c2, const Real & e)
{
  return (e+.5*deltae)*deltae+c1*log(1.-deltac/c1)+c2*log(1.+deltac/c2)-deltac*log((c1-deltac)/(c2+deltac));
}

//Derivative with respect to the charge variation
//of the variation of the functional for ion moves
inline Real d_deltafunc_d_deltac(const Real & deltac,const Real & c1, const Real & c2, const Real & e, const Real & d_deltae_d_deltac)
{
  Real earg=(e+deltac*d_deltae_d_deltac)*d_deltae_d_deltac;
  Real exp_earg=exp(-fabs(earg));
  //cout << earg << " " << exp_earg << endl;
  if (earg<0)
    return deltac-(c1-c2*exp_earg)/(1.+exp_earg);
else

//  else if(earg>0)
    return deltac-(c1*exp_earg-c2)/(exp_earg+1.);
//  else
//    return deltac-.5*(c1-c2);
}

//Lattice sweep
Real sequential_sweep_ion_moves(RVF & electric_field,
                                RSF & concentration,
                                const int & ion_valence,
                                const Real & eps,
                                const Grid & grid)
{
  Real delta_func=0;
for (int i=0;i<grid.nx();++i)  
    for (int j=0;j<grid.ny();++j)
	 for (int k=0;k<grid.nz();++k)
        for (int m=2;m>-1;--m)
        {
          //If the node is moveable
          //Perform the concentration move
          if(grid.point_type()(i,j,k)==1)
            delta_func+=ion_move(electric_field,concentration,ion_valence,IV(i,j,k),m,eps,grid);
        }
  return delta_func;
}

////////////////////////////////////////////


//// MINIMISE ////

void minimise(RVF & electric_field,
              std::vector<RSF> & ion_density,
              const std::vector<int> & ion_valence,
              const int & savingstep, 
              const Real & eps,
              const string & runsname,
              const Grid & grid)
{
  cout << "Minimising... " << endl;
  int num_steps=0;
  Real func0=functional(electric_field,ion_density,grid);
  cout << " Initial value of the functional: " << func0 << endl;
  Real func=func0;
  while(1)
  {
    ++num_steps;
    Real delta_func=0;
      //Field moves
    delta_func+=sequential_sweep_loop_moves(electric_field,grid);
     //Ion moves
   for (int n=0;n<ion_valence.size();++n)
      delta_func+=sequential_sweep_ion_moves(electric_field,ion_density[n],ion_valence[n],eps,grid);
	//Multiply by volume element
	delta_func*=grid.deltav();
      //Update the functional
    func+=delta_func;
      //Print iteration infos
       cout << " Minimisation step: " << num_steps << "\t"
           << "Variation in functional: " << delta_func << "\t"
           << "Functional: " << func << endl;
  //Save fields
    if(savingstep>0&&num_steps%savingstep==0)
    {
      for (int n=0;n<ion_density.size();++n)
      {
        std::stringstream ss;
        ss << runsname << "_ion_" << n << "_density_" << num_steps;
        vtkSave(ss.str(),ion_density[n],ss.str(),grid);
      }
      std::stringstream ss;
      ss << runsname << "_electric_field_" << num_steps;
      vtkSave(ss.str(),electric_field,ss.str(),grid);
    }
      //Check for stop criterium
    if(fabs(delta_func)<eps)
      break;
  }
  cout << " Saving final values..." << endl;
  for (int n=0;n<ion_density.size();++n)
  {
    std::stringstream ss;
    ss << runsname << "_ion_" << n << "_density_" << num_steps;
    vtkSave(ss.str(),ion_density[n],ss.str(),grid);
  }
  std::stringstream ss;
  ss << runsname << "_electric_field_" << num_steps;
  vtkSave(ss.str(),electric_field,ss.str(),grid);
  cout << "...done." << endl;




}


//// DIFFERENTIAL OPERATORS ////

//Divergence
RSF divergence(const RVF & field,const Grid & grid)
{
  RSF div(field.shape());
  div=0;
  int nx=grid.nx();
  int ny=grid.ny();
  int nz=grid.nz();
  for(int i=0;i<nx;++i)
    for(int j=0;j<ny;++j)
      for(int k=0;k<nz;++k)
      {
        Real divx,divy,divz;
        if (i==0)
          divx=field[0](i,j,k)-field[0](nx-1,j,k);
        else
          divx=field[0](i,j,k)-field[0](i-1,j,k);
        divx/=grid.deltax();
        if (j==0)
          divy=field[1](i,j,k)-field[1](i,ny-1,k);
        else
          divy=field[1](i,j,k)-field[1](i,j-1,k);
        divy/=grid.deltay();
        if (k==0)
          divz=field[2](i,j,k)-field[2](i,j,nz-1);
        else
          divz=field[2](i,j,k)-field[2](i,j,k-1);
        divz/=grid.deltaz();
        div(i,j,k)=divx+divy+divz;
      }
  return div;
}

///////////////////////////////////////////////

