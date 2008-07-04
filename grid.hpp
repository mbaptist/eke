// -*- C++ -*-

//
// C++ Interface: grid
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

#ifndef GRID_HPP
#define GRID_HPP

#include "types.hpp"
#include "coord.hpp"
#include "random.hpp"

#include <vector>

class Grid
{
  
private:
  //Members
  int nx_,ny_,nz_;
  Real lx_,ly_,lz_;
  ISF point_type_;
  Real deltax_,deltay_,deltaz_;
  Real deltasx_,deltasy_,deltasz_;
  Real deltav_;
public:
  //Accessors
  int & nx(){return nx_;};
  int & ny(){return ny_;};
  int & nz(){return nz_;};
  Real & lx(){return lx_;};
  Real & ly(){return ly_;};
  Real & lz(){return lz_;};
  Real & deltax(){return deltax_;};
  Real & deltay(){return deltay_;};
  Real & deltaz(){return deltaz_;};
  Real & deltasx(){return deltasx_;};
  Real & deltasy(){return deltasy_;};
  Real & deltasz(){return deltasz_;};
  Real & deltav(){return deltav_;};
  Real & deltas(const int & normal)
  {
    if (normal==0)
      return deltasx_; 
    else if (normal==1)
      return deltasy_;
    else
      return deltasz_;
  };
  Real & deltal(const int & normal)
  {
    if (normal==0)
      return deltax_; 
    else if (normal==1)
      return deltay_;
    else
      return deltaz_;
  };
  const RV coordinates(const int & i,const int & j,const int & k)
  {
    return RV(-.5*lx_+i*deltax_,-.5*ly_+j*deltay_,-.5*lz_+k*deltaz_);
  };
  ISF & point_type(){return point_type_;};
  
  //Constant Accessors
  const int & nx() const {return nx_;};
  const int & ny() const {return ny_;};
  const int & nz() const {return nz_;};
  const Real & lx() const {return lx_;};
  const Real & ly() const {return ly_;};
  const Real & lz() const {return lz_;};
  const Real & deltax() const {return deltax_;};
  const Real & deltay() const {return deltay_;};
  const Real & deltaz() const {return deltaz_;};
  const Real & deltasx() const {return deltasx_;};
  const Real & deltasy() const {return deltasy_;};
  const Real & deltasz() const {return deltasz_;};
  const Real & deltav() const {return deltav_;};	
  const Real & deltas(const int & normal) const
  {
    if (normal==0)
      return deltasx_; 
    else if (normal==1)
      return deltasy_;
    else
      return deltasz_;
  };
  const Real & deltal(const int & normal) const
  {
    if (normal==0)
      return deltax_; 
    else if (normal==1)
      return deltay_;
    else
      return deltaz_;
  };
  const RV coordinates(const int & i,const int & j,const int & k) const 
  {
    return RV(-.5*lx_+i*deltax_,-.5*ly_+j*deltay_,-.5*lz_+k*deltaz_);
  };
  const ISF & point_type() const {return point_type_;};
  
  
  //Methods
  void intobox(IV & node) const
  {
    node[0]%=nx_;
    node[1]%=ny_;
    node[2]%=nz_;
  };
  
  const IV nearest_point_index(const RV & coord)
  {
    int ind_lip_x=static_cast<int>((coord[0]+.5*lx_)/deltax_);
    int ind_lip_y=static_cast<int>((coord[1]+.5*ly_)/deltay_);
    int ind_lip_z=static_cast<int>((coord[2]+.5*lz_)/deltaz_);
    RV dist_vec=coord;
    dist_vec-=coordinates(ind_lip_x,ind_lip_y,ind_lip_z);
    Real distance0=norm(dist_vec);
    int ind_np_x=ind_lip_x;
    int ind_np_y=ind_lip_y;
    int ind_np_z=ind_lip_z;
    for(int i=0;i<2;++i)
      for(int j=0;j<2;++j)
        for(int k=0;k<2;++k)
        {
          dist_vec=coord;
          dist_vec-=coordinates(ind_lip_x+i,ind_lip_y+j,ind_lip_z+k);
          Real distance=norm(dist_vec); 
          if(distance<distance0)
          {
            ind_np_x=ind_lip_x+i;
            ind_np_y=ind_lip_y+j;
            ind_np_z=ind_lip_z+k;
            distance0=distance;
          }
        }
    return IV(ind_np_x,ind_np_y,ind_np_z);
  };
  
public:
  //Ctors
  Grid(const int & _nx_, const int & _ny_,const int & _nz_,
       const Real & _lx_, const Real & _ly_,const Real & _lz_):
    nx_(_nx_),ny_(_ny_),nz_(_nz_),
    lx_(_lx_),ly_(_ly_),lz_(_lz_),
    point_type_(nx_,ny_,nz_),
    deltax_(lx_/nx_),
    deltay_(ly_/ny_),
    deltaz_(lz_/nz_),
    deltasx_(deltay_*deltaz_),
    deltasy_(deltax_*deltaz_),
    deltasz_(deltay_*deltaz_),
    deltav_(deltax_*deltay_*deltaz_)
  {};
  //Dtor
  ~Grid(){};
  
private:
  //Forbidden Ctors
  Grid();
  
};



class Loop
{
private:
  //Members
  IV node1_;
  IV node2_;
  IV node3_;
  IV node4_;
  int loop_number_;
  int dir1_,dir2_;
public:
  //Accessors  
  IV & node1(){return node1_;};
  IV & node2(){return node2_;};
  IV & node3(){return node3_;};
  IV & node4(){return node4_;};
  int & loop_number(){return loop_number_;};
  int & dir1(){return dir1_;};
  int & dir2(){return dir2_;};
  //Const Accessors
  const IV & node1() const {return node1_;};
  const IV & node2() const {return node2_;};
  const IV & node3() const {return node3_;};
  const IV & node4() const {return node4_;};
  const int & loop_number() const {return loop_number_;};
  const int & dir1() const {return dir1_;};
  const int & dir2() const {return dir2_;};
public:
  //Ctor
  Loop(const IV & _node_,const int & _loop_number_,const Grid & _grid_):
    node1_(_node_),
    node2_(_node_),
    node3_(_node_),
    node4_(_node_),
    loop_number_(_loop_number_)
  {
    if (loop_number_==0)
    {
      dir1_=1;
      dir2_=2;
    }
    else if (loop_number_==1)
    {
      dir1_=0;
      dir2_=2;
    }
    else if (loop_number_==2)
    {
      dir1_=0;
      dir2_=1;
    }
    node2_+=unit_vector(dir1_);
    node3_+=unit_vector(dir1_)+unit_vector(dir2_);
    node4_+=unit_vector(dir2_);
    _grid_.intobox(node2_);
    _grid_.intobox(node3_);
    _grid_.intobox(node4_);
  };
  
  //Dtor
  ~Loop(){};
private:
  //Forbidden Ctors  
  Loop();
  
};



#endif
