
class Grid
{
  
private:
  //Members
  int nx_,ny_,nz_;
  double lx_,ly_,lz_,rs_;
  Array<TinyVector<double,3>,3> coordinates_;
  Array<int,3> point_type_;

public:
  //Accessors
  int & nx(){return & nx_;};
  int & ny(){return & ny_;};
  int & nz(){return & nz_;};
  int & lx(){return & lx_;};
  int & ly(){return & ly_;};
  int & lz(){return & lz_;};
  int & rs(){return & rs_;};
  Array<TinyVector<double,3>,3> & coordinates(){return & coordinates_;};
  Array<TinyVector<double,3>,3> & point_type(){return & point_type_;};
  
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










class Plaquete
{
private:
  blitz::TinyVector<int,3> node1_;
  blitz::TinyVector<int,3> node2_;
  blitz::TinyVector<int,3> node3_;
  blitz::TinyVector<int,3> node4_;
  int plaquete_;
  int dir1_,dir2_;

  

public:
  Plaquete(const blitz::TinyVector<int,3> & _node_,const int & _plaquete_):
    node1_(_node_),
    node2_(_node_),
    node3_(_node_),
    node4_(_node_),
    plaquete_(_plaquete_)
  {
    if (plaquete_==0)
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
    else if (plaquete_==1)
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
    else if (plaquete_==2)
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


private:
  Plaquete(){};
  ~Plaquete(){};
public:


};
