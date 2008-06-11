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
#include <vector>
#include "random.hpp"

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
  const RV coordinates(const int & i,const int & j,const int & k){return RV(-.5*lx_+i*deltax_,-.5*ly_+j*deltay_,-.5*lz_+k*deltaz_);};
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
  const RV coordinates(const int & i,const int & j,const int & k) const {return RV(-.5*lx_+i*deltax_,-.5*ly_+j*deltay_,-.5*lz_+k*deltaz_);};
  const ISF & point_type() const {return point_type_;};
  

  //Methods
  IV intobox(IV & node) const
  {
    return IV(node[0]%nx_,node[1]%ny_,node[2]%nz_);
  };
  
  const IV nearest_point_index(const RV & coord)
  {
    int ind_lip_x=static_cast<int>((coord[0]+.5*lx_)/(lx_/nx_));
    int ind_lip_y=static_cast<int>((coord[1]+.5*ly_)/(ly_/ny_));
    int ind_lip_z=static_cast<int>((coord[2]+.5*lz_)/(lz_/nz_));

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
    int nx(_grid_.nx());
    int ny(_grid_.ny());
    int nz(_grid_.nz());
    if (loop_number_==0)
    {
      node2_+=blitz::TinyVector<int,3>(0,1,0);
      node3_+=blitz::TinyVector<int,3>(1,1,0);
      node4_+=blitz::TinyVector<int,3>(1,0,0);
      if(node2_[1]==ny)
        node2_[1]=0;
      if(node3_[0]==nx)
        node3_[0]=0;
      if(node3_[1]==ny)
        node3_[1]=0;
      if(node4_[0]==nx)
        node4_[0]=0;
      dir1_=1;
      dir2_=0;
    }
    else if (loop_number_==1)
    {
      node2_+=blitz::TinyVector<int,3>(0,0,1);
      node3_+=blitz::TinyVector<int,3>(1,0,1);
      node4_+=blitz::TinyVector<int,3>(1,0,0);
      if(node2_[2]==nz)
        node2_[2]=0;
      if(node3_[0]==nx)
        node3_[0]=0;
      if(node3_[2]==nz)
        node3_[2]=0;
      if(node4_[0]==nx)
        node4_[0]=0;
      dir1_=2;
      dir2_=0;
    }
    else if (loop_number_==2)
    {
      node2_+=blitz::TinyVector<int,3>(0,0,1);
      node3_+=blitz::TinyVector<int,3>(0,1,1);
      node4_+=blitz::TinyVector<int,3>(0,1,0);
      if(node2_[2]==nz)
        node2_[2]=0;
      if(node3_[1]==ny)
        node3_[1]=0;
      if(node3_[2]==nz)
        node3_[2]=0;
      if(node4_[1]==ny)
        node4_[1]=0;
      dir1_=2;
      dir2_=1;
    }	 
  };
  
  //Dtor
  ~Loop(){};
private:
  //Forbidden Ctors  
  Loop();
  
};



class Path
{
  //Members
private:
  std::vector<IV> path_;
  int size_;
  std::vector<IV> direction_;
  std::vector<int> direction_axis_;
  
  //Accessors
  
public:
  const std::vector<IV> & operator()(){return path_;};
  const IV & operator()(const int & n){return path_[n];}
  const int & size(){return size_;};
  const std::vector<IV> & direction(){return direction_;};
  const IV & direction(const int & n){return direction_[n];};
  const int & direction_axis(const int & n){direction_axis_[n];};
    
  const std::vector<IV> & operator()() const {return path_;};
  const IV & operator()(const int & n) const {return path_[n];}
  const int & size() const {return size_;};
  const std::vector<IV> & direction() const {return direction_;};
  const IV & direction(const int & n) const {return direction_[n];};
  const int & direction_axis(const int & n) const {direction_axis_[n];};
  
  //Ctors
public:
#if 0
  Path(const IV & _node1_,const IV & _node2_)
  {
    path_.push_back(_node1_);
    for(int i=_node1_[0]+1;i<=_node2_[0];++i)
    {
      path_.push_back(IV(_node1_[0]+i,_node1_[1],_node1_[2]));
      direction_axis_.push_back(0);
      direction_.push_back(IV(1,0,0));
    }
    for(int j=_node1_[1]+1;j<=_node2_[1];++j)
    {
      path_.push_back(IV(_node1_[0],_node1_[1]+j,_node1_[2]));
      direction_axis_.push_back(1);
      direction_.push_back(IV(0,1,0));
    }
    for(int k=_node1_[2]+1;k<_node2_[2];++k)
    {
      path_.push_back(IV(_node1_[0],_node1_[1],_node1_[2]+k)); 
      direction_axis_.push_back(2);
      direction_.push_back(IV(0,0,1));
    }
    path_.push_back(IV(_node1_[0],_node1_[1],_node1_[2]+k));
    
  };
#endif
  Path(const IV & _node1_,const int & _dir_,const Grid & grid)
  {
    path_.push_back(_node1_);
    direction_axis_.push_back(_dir_);
    direction_.push_back(IV(0,0,0)); 
    direction_[0][direction_axis_[0]]=1;
    IV node2(_node1_);
    node2+=direction_[0];
    //Apply periodic boundary
    node2=grid.intobox(node2);
    path_.push_back(node2);
    size_=2;
  };
  //Dtor
     ~Path(){};

private:
  //Forbidden Ctors
  Path();
  
  
};



#endif
