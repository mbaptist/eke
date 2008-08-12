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
	
IV intobox(const IV & node) const
{
	return IV(node[0]%nx_,node[1]%ny_,node[2]%nz_);
}

  const IV nearest_point_index(const RV & coord)
  { 
#if 0
	IV np(static_cast<int>(floor((coord[0]+.5*lx_)/deltax_)),
               static_cast<int>(floor((coord[1]+.5*ly_)/deltay_)),
               static_cast<int>(floor((coord[2]+.5*lz_)/deltaz_)));
	RV dist=coord-coordinates(np[0],np[1],np[2]);
	for(int n=0;n<3;++n)
		if(dist[n]>=.5*deltal(n))
			np[n]+=1;
	return np;
#endif
   return IV(static_cast<int>(round((coord[0]+.5*lx_)/deltax_))%nx_,
               static_cast<int>(round((coord[1]+.5*ly_)/deltay_))%ny_,
              static_cast<int>(round((coord[2]+.5*lz_)/deltaz_))%nz_);
  };
  
public:
  //Ctors
  explicit Grid(const int & _nx_, const int & _ny_,const int & _nz_,
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
  explicit Loop(const IV & _node_,const int & _loop_number_,const Grid & _grid_):
    node1_(_node_),
    loop_number_(_loop_number_),
    dir1_((loop_number_+1)%3),
    dir2_((loop_number_+2)%3)
  {
    node2_=_grid_.intobox(node1_+unit_vector(dir1_));
    node3_=_grid_.intobox(node2_+unit_vector(dir2_));
    node4_=_grid_.intobox(node1_+unit_vector(dir2_));
  };
  
  //Dtor
  ~Loop(){};
private:
  //Forbidden Ctors  
  Loop();
  
};



#endif
