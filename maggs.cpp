

#include "maggs.hpp"

#include <blitz/blitz.h>
#include <blitz/array.h>
#include <blitz/tiny.h>
#include <blitz/tinyvec.h>

using namespace blitz;

#include <blitz/random.h>
#include <random/uniform.h>
#include <random/discrete-uniform.h>

using namespace ranlib;



void loop_move(Array<TinyVector<double,3>,3> & electric_field,
               const Grid & grid,
               const Loop & loop)
{
//Evaluate the trial field
  double & e1 = electric_field(loop.node1()[0],loop.node1()[1],loop.node1()[2])[loop.dir1()];
  double & e2 = electric_field(loop.node2()[0],loop.node2()[1],loop.node2()[2])[loop.dir2()];
  double & e3 = electric_field(loop.node3()[0],loop.node3()[1],loop.node3()[2])[loop.dir1()];
  double & e4 = electric_field(loop.node4()[0],loop.node4()[1],loop.node4()[2])[loop.dir2()];
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
  }
  //cout << "Circulation: " << e1+e2-e3-e4 << endl;
};

void sequential_sweep_loop_moves(Array<TinyVector<double,3>,3> & electric_field, 
                                 const Grid & grid)
{
  double delta_func;
  for(int i=0;i<grid.nx();++i)
    for(int j=0;j<grid.ny();++j)
      for(int k=0;k<grid.nz();++k)
        for (int loop_number=0;loop_number<3;++loop_number)
        { 
//Create the loop
          Loop loop(blitz::TinyVector<int,3>(i,j,k),loop_number,grid);
//Perform the loop move
          loop_move(electric_field,grid,loop);
        }
}

void random_sweep_loop_moves(Array<TinyVector<double,3>,3> & electric_field,
                             const Grid & grid)
{
  //Choose 3x the number of nodes
  for(int nn=0;nn<3*electric_field.size();++nn)
  {
  //Choose a node
    DiscreteUniform<int> rgx(grid.nx());
    DiscreteUniform<int> rgy(grid.ny());
    DiscreteUniform<int> rgz(grid.nz());
    int ni = rgx.random();
    int nj = rgy.random();
    int nk = rgz.random();
  //cout << "Node: " << ni << " " << nj << " " << nk << endl;
  //Choose a loop
    DiscreteUniform<int> rgp(3);
    int loop_number=rgp.random();
  //cout << "Loop number: " << loop_number << endl;
  //Create the loop
    Loop loop(blitz::TinyVector<int,3>(ni,nj,nk),loop_number,grid);
  //Perform the loop move
    loop_move(electric_field,grid,loop);
  }
}
