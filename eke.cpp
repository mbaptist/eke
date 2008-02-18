
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

int main()
{

#if 0
  //cout << fixed << setprecision(20);

  //Input (in future via python)
  int nx,ny,nz;
  double lx,ly,lz;
  double rs;
  double charge;
  nx=32;
  ny=nx;
  nz=ny;
  lx=10.;
  ly=10.;
  lz=10.;
  rs=1.;
  charge=1.;

  
  GeometryObject<block>();
  

  Geometry geometry();


  Grid grid(n1,n2,n3,geometry);

  
#endif



//cout << fixed << setprecision(20);

  //Input (in future via python)
  int nx,ny,nz;
  double lx,ly,lz;
  double rs;
  double charge;
  nx=32;
  ny=nx;
  nz=ny;
  lx=10.;
  ly=10.;
  lz=10.;
  rs=1.;
  charge=1.;

  //Grid
  //Coordinates of grid points
  Array<TinyVector<double,3>,3> coordinates(nx,ny,nz);
  //Classification of grid points
  Array<int,3> point_type(nx,ny,nz);
  //Initialise the arrays above
  initialise_grid(coordinates,point_type,nx,ny,nz,lx,ly,lz,rs);

  //Charge density
  Array<double,3> charge_density(nx,ny,nz);
  //Distribute the charges over the surface of the colloid
  initialise_colloid(charge_density,point_type,coordinates,
		     nx,ny,nz,lx,ly,lz,rs,charge);
  

  ofstream ofs("chd.dat");
  //ofs << charge_density << endl;
  for (int n=0;n<charge_density.size();++n)
    ofs << charge_density.data()[n] << endl;
 
  //cout << charge_density << endl;



  //electric field
  Array<TinyVector<double,3>,3> electric_field(nx,ny,nz);

#if 0
  for(int i=0;i<nx;++i)
    for(int j=0;j<ny;++j)
      for(int k=0;k<nz;++k)
	if(point_type(i,j,k)==3||point_type(i,j,k)==2||point_type(i,j,k)==3)
	  {
	    electric_field(i,j,k)=coordinates(i,j,k);
	    electric_field(i,j,k)*=.25/M_PI*charge/pow(norm(coordinates(i,j,k)),3);
	  }
	else
	  electric_field(i,j,k)=TinyVector<double,3>(0,0,0);
#endif


  initialise_electric_field(electric_field,charge_density,point_type,nx,ny,nz,lx,ly,lz);



  ofstream ofs2("efi.dat");
  //ofs2 << electric_field << endl;
  for (int n=0;n<electric_field.size();++n)
    {
      for (int m=0;m<3;++m)
	ofs2 << electric_field.data()[n][m] <<" ";
      ofs2 << endl;
    }

  //minimise
  int npm=0;
  double delta_functional=1;
  double func0=.5*sum(dot(electric_field,electric_field));
  while(1)
    {
      ++npm;
      //cout << ++npm << endl;
#if 0
      for (int nn=0;nn<300*(nx*ny*nz);++nn)
	{
	  delta_functional=0;
	  delta_functional+=loop_move(electric_field,nx,ny,nz,lx,ly,lz);
	  //delta_functional+=charge_move(electric_field,charge_density);
	}
#endif
      delta_functional=0;
      delta_functional+=sequential_loop_move(electric_field,nx,ny,nz,lx,ly,lz);
      

      double func=.5*sum(dot(electric_field,electric_field));
      cout << "F=" << func << endl;
      delta_functional=func-func0;
      func0=func;
      cout << "DF=" << delta_functional << endl;
      if(delta_functional<0)
	cout << npm << endl;
      if((delta_functional<=0)&&(-delta_functional<1e-20))
	break;
    }
  
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
	if(point_type(i,j,k)==3)
	  {
	    //cout << "EF: " << electric_field(i,j,k) << endl;
	    //cout << "TT" << endl;
	    ef_test(i,j,k)=coordinates(i,j,k);
	    ef_test(i,j,k)*=.25/M_PI*charge/pow(norm(coordinates(i,j,k)),3);
	    ef_test(i,j,k)-=electric_field(i,j,k);
	    ef_test(i,j,k)/=norm(electric_field(i,j,k));
	    //cout << "EFT: " << ef_test(i,j,k) << endl;
	  }
	else
	  ef_test(i,j,k)=0;

  cout << "TEST: " << sum(dot(ef_test,ef_test))/(nx*ny*nz) << endl;
  cout << max(dot(ef_test,ef_test)) << endl;
  //cout << "TEST: " << sum(dot(electric_field,electric_field))/(nx*ny*nz) << endl;

  return 0;
}



void initialise_grid(Array<TinyVector<double,3>,3> & coordinates,
		     Array<int,3> & point_type,
		     const int & nx, const int & ny,const int & nz,
		     const double & lx, const double & ly,const double & lz,
		     const double & rs)
{
  for(int i=0;i<nx;++i)
    for(int j=0;j<ny;++j)
      for(int k=0;k<nz;++k)
	{
	  coordinates(i,j,k)=TinyVector<double,3>(-.5*lx+i*lx/nx,
						  -.5*ly+j*ly/ny,-
						  .5*lz+k*lz/nz);
	  if (norm(coordinates(i,j,k))>=2*rs)
	    point_type(i,j,k)=1;
	  else if ((norm(coordinates(i,j,k))<2*rs)&&(norm(coordinates(i,j,k))>=rs))
	    point_type(i,j,k)=3;
	  else
	    point_type(i,j,k)=5;
	}
}
 

void initialise_colloid(Array<double,3> & charge_density,
			Array<int,3> & point_type,
			const Array<TinyVector<double,3>,3> & coordinates,
			const int & nx, const int & ny,const int & nz,
			const double & lx, const double & ly,const double & lz,
			const double & rs,const double & charge)
{
  charge_density=0;
  int np=10000;
  double delta_charge=charge/np;
  UniformClosed<double> rg;
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
      charge_density(ind_np_x,ind_np_y,ind_np_z)+=delta_charge;
      point_type(ind_np_x,ind_np_y,ind_np_z)=4;     
    }
  
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

      charge_density(ind_np_x,ind_np_y,ind_np_z)+=-delta_charge;
      point_type(ind_np_x,ind_np_y,ind_np_z)=2;     
    }     
      

}


void initialise_electric_field(Array<TinyVector<double,3>,3> & electric_field,
			       const Array<double,3> & charge_density,
			       const Array<int,3> & point_type,
			       const int & nx, const int & ny,const int & nz,
			       const double & lx, const double & ly,const double & lz)
{


#if 1

  electric_field=TinyVector<double,3>(0,0,0);

  Array<double,3> charge(charge_density.copy());

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

#endif

  cout << "Mean: " << mean(electric_field[0]) << endl;
  cout << "Mean: " << mean(electric_field[1]) << endl;
  cout << "Mean: " << mean(electric_field[2]) << endl;
  electric_field[0]-=mean(electric_field[0]);
  electric_field[1]-=mean(electric_field[1]);
  electric_field[2]-=mean(electric_field[2]);
  cout << "Mean: " << mean(electric_field[0]) << endl;
  cout << "Mean: " << mean(electric_field[1]) << endl;
  cout << "Mean: " << mean(electric_field[2]) << endl;

  //exit(0);

  

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
	  if (fabs(div-charge_density(i,j,k))>.5e-16)
	    {
	      cout << i << " " << j << " " << k << " " 
		   << divx << " " << divy <<" " << divz << " " 
		   << div << " " << charge_density(i,j,k) << endl;
	      exit(1);
	    }
	}
#endif 

}

double random_loop_move(Array<TinyVector<double,3>,3> & electric_field,
		 const int & nx,const int & ny,const int & nz,
		 const double & lx,const double & ly,const double & lz)
{
  //Choose a node
  DiscreteUniform<int> rgx(nx);
  DiscreteUniform<int> rgy(ny);
  DiscreteUniform<int> rgz(nz);
  int ni = rgx.random();
  int nj = rgy.random();
  int nk = rgz.random();
  //cout << "Node: " << ni << " " << nj << " " << nk << endl;
  //Choose a plaquete
  DiscreteUniform<int> rgp(3);
  int plaquete=rgp.random();
  //cout << "Plaquete: " << plaquete << endl;
  return try_move(electric_field,nx,ny,nz,ni,nj,nk,plaquete);
}


double try_move(Array<TinyVector<double,3>,3> & electric_field,
		const int & nx,const int & ny,const int & nz,
		const int & ni,const int & nj,const int & nk,
		const int & plaquete)
{
  int pn1i=ni,pn1j=nj,pn1k=nk;
  int pn2i=ni,pn2j=nj,pn2k=nk;
  int pn3i=ni,pn3j=nj,pn3k=nk;
  int pn4i=ni,pn4j=nj,pn4k=nk;
  int dir1,dir2;
  if (plaquete==0)
    {
      if(nj<ny-1)
	pn2j=nj+1;
      else
	pn2j=0;
      if(ni<nx-1)
	pn3i=ni+1;
      else
	pn3i=0;
      dir1=1;
      dir2=0;
    }
  else if (plaquete==1)
    {
       if(nk<nz-1)
	 pn2k=nk+1;
       else
	 pn2k=0;
       if(nj<ny-1)
	 pn3j=nj+1;
       else
	 pn3j=0;
      dir1=2;
      dir2=1;
    }
  else if (plaquete==2)
    {
       if(nk<nz-1)
	 pn2k=nk+1;
       else
	 pn2k=0;
       if(ni<nx-1)
	 pn3i=ni+1;
       else
	 pn3i=0;
      dir1=2;
      dir2=0;
    }	
  //Change in the field
  //UniformClosed<double> rgdf;
  //double mdf=.1;
  //double delta_field=mdf*(-1.+2.*rgdf.random());
  //cout << "Dfield=" << delta_field << endl;
  //Evaluate the trial field
  double & e1 = electric_field(pn1i,pn1j,pn1k)[dir1];
  double & e2 = electric_field(pn2i,pn2j,pn2k)[dir2];
  double & e3 = electric_field(pn3i,pn3j,pn3k)[dir1];
  double & e4 = electric_field(pn4i,pn4j,pn4k)[dir2];
  double delta_field=-.25*(e1+e2-e3-e4);
  double e1p=e1+delta_field;
  double e2p=e2+delta_field;
  double e3p=e3-delta_field;
  double e4p=e4-delta_field;
  //Evaluate the change in the functional
  double delta_func=(e1p*e1p+e2p*e2p+e3p*e3p+e4p*e4p
		     -e1*e1-e2*e2-e3*e3-e4*e4);
  //cout << "Dfunc=" << delta_func << endl;
  //Accept/reject move
  if(delta_func<0)
    {
      e1=e1p;
      e2=e2p;
      e3=e3p;
      e4=e4p;
      //cout << delta_func << endl;
    }
  //Return the change in the functional
  return(delta_func);
};




double sequential_loop_move(Array<TinyVector<double,3>,3> & electric_field,
		 const int & nx,const int & ny,const int & nz,
		 const double & lx,const double & ly,const double & lz)
{
  double delta_func;
  for(int i=0;i<nx;++i)
    for(int j=0;j<ny;++j)
      for(int k=0;k<nz;++k)
	for (int plaquete=0;plaquete<3;++plaquete)
	  delta_func=try_move(electric_field,nx,ny,nz,i,j,k,plaquete);
  return delta_func;
}
