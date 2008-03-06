
#ifndef GRID_HPP
#define GRID_HPP


#include <blitz/blitz.h>
#include <blitz/array.h>
#include <blitz/tiny.h>
#include <blitz/tinyvec.h>

class Grid
{
  
private:
  //Members
  int nx_,ny_,nz_;
  double lx_,ly_,lz_,rs_;
  blitz::Array<blitz::TinyVector<double,3>,3> coordinates_;
  blitz::Array<int,3> point_type_;

public:
  //Accessors
  int & nx(){return nx_;};
  int & ny(){return ny_;};
  int & nz(){return nz_;};
  double & lx(){return lx_;};
  double & ly(){return ly_;};
  double & lz(){return lz_;};
  double & rs(){return rs_;};
  blitz::Array<blitz::TinyVector<double,3>,3> & coordinates(){return coordinates_;};
  blitz::Array<int,3> & point_type(){return point_type_;};
  //Constant Accessors
  const int & nx() const {return nx_;};
  const int & ny() const {return ny_;};
  const int & nz() const {return nz_;};
  const double & lx() const {return lx_;};
  const double & ly() const {return ly_;};
  const double & lz() const {return lz_;};
  const double & rs() const {return rs_;};
  const blitz::Array<blitz::TinyVector<double,3>,3> & coordinates() const {return coordinates_;};
  const blitz::Array<int,3> & point_type() const {return point_type_;};
  
public:
  //Ctors
  Grid(const int & _nx_, const int & _ny_,const int & _nz_,
       const double & _lx_, const double & _ly_,const double & _lz_,
       const double & _rs_);
  //Dtor
  ~Grid();
  
private:
  //Forbidden Ctors
  Grid();

};



class Loop
{
private:
  //Members
  blitz::TinyVector<int,3> node1_;
  blitz::TinyVector<int,3> node2_;
  blitz::TinyVector<int,3> node3_;
  blitz::TinyVector<int,3> node4_;
  int loop_number_;
  int dir1_,dir2_;
public:
  //Accessors  
  blitz::TinyVector<int,3> & node1(){return node1_;};
  blitz::TinyVector<int,3> & node2(){return node2_;};
  blitz::TinyVector<int,3> & node3(){return node3_;};
  blitz::TinyVector<int,3> & node4(){return node4_;};
  int & loop_number(){return loop_number_;};
  int & dir1(){return dir1_;};
  int & dir2(){return dir2_;};
  //Const Accessors
  const blitz::TinyVector<int,3> & node1() const {return node1_;};
  const blitz::TinyVector<int,3> & node2() const {return node2_;};
  const blitz::TinyVector<int,3> & node3() const {return node3_;};
  const blitz::TinyVector<int,3> & node4() const {return node4_;};
  const int & loop_number() const {return loop_number_;};
  const int & dir1() const {return dir1_;};
  const int & dir2() const {return dir2_;};
  

public:
  //Ctor
  Loop(const blitz::TinyVector<int,3> & _node_,const int & _loop_number_,const Grid & _grid_):
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
	node3_+=blitz::TinyVector<int,3>(1,0,0);
	if(node2_[1]==ny)
	  node2_[1]=0;
	if(node3_[0]==nx)
	  node3_[0]=0;
	dir1_=1;
	dir2_=0;
      }
    else if (loop_number_==1)
      {
	node2_+=blitz::TinyVector<int,3>(0,0,1);
	node3_+=blitz::TinyVector<int,3>(0,1,0);
	if(node2_[2]==nz)
	  node2_[2]=0;
	if(node3_[1]==ny)
	  node3_[1]=0;
	dir1_=2;
	dir2_=1;
      }
    else if (loop_number_==2)
      {
	node2_+=blitz::TinyVector<int,3>(0,0,1);
	node3_+=blitz::TinyVector<int,3>(1,0,0);
	if(node2_[2]==nz)
	  node2_[2]=0;
	if(node3_[0]==nx)
	  node3_[0]=0;
	dir1_=2;
	dir2_=0;
      }	 
  };
  
  //Dtor
  ~Loop(){};
private:
  //Forbidden Ctors  
  Loop();
  
};


#endif
