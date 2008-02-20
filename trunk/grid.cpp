
#include "grid.hpp"

#include <blitz/blitz.h>
#include <blitz/array.h>
#include <blitz/tiny.h>
#include <blitz/tinyvec.h>
using namespace blitz;


//Ctor
Grid::Grid(const int & _nx_, const int & _ny_,const int & _nz_,
           const double & _lx_, const double & _ly_,const double & _lz_,
           const double & _rs_):
nx_(_nx_),ny_(_ny_),nz_(_nz_),
lx_(_lx_),ly_(_ly_),lz_(_lz_),
rs_(_rs_),
coordinates_(nx_,ny_,nz_),
point_type_(nx_,ny_,nz_)
{
  for(int i=0;i<nx_;++i)
    for(int j=0;j<ny_;++j)
      for(int k=0;k<nz_;++k)
      {
        coordinates_(i,j,k)=TinyVector<double,3>(-.5*lx_+i*lx_/nx_,
                                                -.5*ly_+j*ly_/ny_,
                                                -.5*lz_+k*lz_/nz_);
        if (norm(coordinates_(i,j,k))>=2*rs_)
          point_type_(i,j,k)=1;
        else if ((norm(coordinates_(i,j,k))<2*rs_)&&(norm(coordinates_(i,j,k))>=rs_))
          point_type_(i,j,k)=3;
        else
          point_type_(i,j,k)=5;
      }
}

//Dtor
Grid::~Grid()
{
}

