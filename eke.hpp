

#include <blitz/blitz.h>
#include <blitz/array.h>
#include <blitz/tiny.h>
#include <blitz/tinyvec.h>

using namespace blitz;

#include "grid.hpp"

void initialise_colloid(Array<double,3> & charge_density, const Grid & grid, const double & charge);


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




