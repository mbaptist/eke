

#include <blitz/blitz.h>
#include <blitz/array.h>
#include <blitz/tiny.h>
#include <blitz/tinyvec.h>

using namespace blitz;


void initialise_grid(Array<TinyVector<double,3>,3> & coordinates,
		     Array<int,3> & point_type,
		     const int & nx, const int & ny,const int & nz,
		     const double & lx, const double & ly,const double & lz,
		     const double & rs);

void initialise_colloid(Array<double,3> & charge_density,
			Array<int,3> & point_type,
			const Array<TinyVector<double,3>,3> & coordinates,
			const int & nx, const int & ny,const int & nz,
			const double & lx, const double & ly,const double & lz,
			const double & rs, const double & charge);


void initialise_electric_field(Array<TinyVector<double,3>,3> & electric_field,
			       const Array<double,3> & charge_density,
			       const Array<int,3> & point_type,
			       const int & nx, const int & ny,const int & nz,
			       const double & lx, const double & ly,const double & lz);  


double random_loop_move(Array<TinyVector<double,3>,3> & electric_field,
		 const int & nx,const int & ny,const int & nz,
		 const double & lx, const double & ly,const double & lz);

double sequential_loop_move(Array<TinyVector<double,3>,3> & electric_field,
		 const int & nx,const int & ny,const int & nz,
			    const double & lx,const double & ly,const double & lz);


double try_move(Array<TinyVector<double,3>,3> & electric_field,
		const int & nx,const int & ny,const int & nz,
		const int & ni,const int & nj,const int & nk,
		const int & plaquete);




class Grid()
{
  

}










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
    if (plaquete==0)
      {
	node2_+=blitz::TinyVector<int,3>(0,1,0);
	node3_+=blitz::TinyVector<int,3>(1,0,0);
	if(node2_[1]==ny)
	  node2_[1]=0;
	if(node3_[0]==nx)
	  node3_[0]=0;
	dir1=1;
	dir2=0;
      }
    else if (plaquete==1)
      {
	node2_+=blitz::TinyVector<int,3>(0,0,1);
	node3_+=blitz::TinyVector<int,3>(0,1,0);
	if(node2_[2]==nz)
	  node2_[2]=0;
	if(node3_[1]==ny)
	  node3_[1]=0;
	dir1=2;
	dir2=1;
      }
    else if (plaquete==2)
      {

	node2_+=blitz::TinyVector<int,3>(0,0,1);
	node3_+=blitz::TinyVector<int,3>(1,0,0);
	if(node2_[2]==nz)
	  node2_[2]=0;
	if(node3_[0]==nx)
	  node3_[0]=0;
	dir1=2;
	dir2=0;
      }	 
  };


private:
  Plaquete(){};
  ~Plaquete(){};
public:


};
