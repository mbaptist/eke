

#include <blitz/blitz.h>
#include <blitz/array.h>
#include <blitz/tiny.h>
#include <blitz/tinyvec.h>

using namespace blitz;

#include "grid.hpp"

void initialise_colloid(Array<double,3> & charges, 
                        Grid & grid,
                        const double & total_charge);

void initialise_electric_field(Array<TinyVector<double,3>,3> & electric_field,
                               const Array<double,3> & charges,
                               const Grid & grid);  

double functional(Array<TinyVector<double,3>,3> & electric_field);

void loop_move(Array<TinyVector<double,3>,3> & electric_field,
               const Grid & grid,
               const Loop & loop);

void sequential_sweep_loop_moves(Array<TinyVector<double,3>,3> & electric_field, 
                            const Grid & grid);
void random_sweep_loop_moves(Array<TinyVector<double,3>,3> & electric_field, 
                            const Grid & grid);

