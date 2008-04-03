// -*- C++ -*-
//
// //
// // C++ Interface: input
// //
// // Description:
// //
// //
// // Author: Manuel Baptista <mbaptist@gmail.com>, (C) 2008
// //
// // Copyright: See COPYING file that comes with this distribution
// //
// //
//
// /*
// Copyright (C) 2008 Manuel Baptista <mbaptist@gmail.com>
//
// This file is part of eke.
//
// Eke is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Eke is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with eke.  If not, see <http://www.gnu.org/licenses/>.
// */
//


//input.hpp

//Class Input declares the Input parameters and is inherited
//by the problem class ns3d
//It is responsible for reading the configuration file and
//seting the Input parameters

#ifndef INPUT_HPP
#define INPUT_HPP

#include <Python.h>

#include <iostream>
#include <string>

using namespace std;

//Forward declaration of class py_input_parser
class py_input_parser;


///////////////////////////////////////////
// Declaration of class Input
///////////////////////////////////////////
class Input
  {
  public:
    //constructor from filename
    Input(std::string cfg_fname);
    //constructor from python object
    Input(PyObject * module);
    //destructor
    ~Input();
  public:
    string runsname;
    //Input parameters
    int n1,n2,n3;//Spatial grid
    double l1,l2,l3,ell1,ell2;//Spatial extent
    //physical Parameters
    double visc,omegaz,compress,g,diff,deltat,econd,tcond;

    //refine and resume
    bool refine,resume;
    std::string lr_runsname;
    int lr_n1,lr_n2,lr_n3;

    //basic fields
    std::string basic_mode;
    //Random Mode
    std::string br_spectrum;
    int br_seed,br_ki,br_kf;
    double 	br_alpha,br_rms_norm;
    bool	br_kind,br_sym;


    //large scale stability

    double ls_eps;
    double qq;
    double qq_adj;
    int kk;
    double small;
    double small_adj;

    //short scale stability
    int sym_sub,mp,nseq,sss_seed,sss_basic_restore_seed;
    double ep,thr,sc;
    std::string sss_ifname,sss_int_ofbname;

  private:
    //load from config file
    void py_load_Input(py_input_parser & parse);
  };


///////////////////////////////////////////
// Declaration of class py_input_parser
///////////////////////////////////////////
class py_input_parser
  {
  private:
    //Filename to parse
    std::string filename;
    //Python objects
    PyObject * module;
    PyObject * dict;
    PyObject * pval(const string & item);
    PyThreadState * inter;
    bool python_initialised;
  public:
    //Constructor from filename
    py_input_parser(string filename_);
    //Constructor from dict
    py_input_parser(PyObject * module_);
    //Destructor
    ~py_input_parser();
    //Python to C++ conversion function
    void operator()(int & val,const std::string & item);
    void operator()(bool & val,const std::string & item);
    void operator()(double & val,const std::string & item);
    void operator()(std::string & val,const std::string & item);
  private:
    //Forbidden constructors
    py_input_parser();//default
    py_input_parser(const py_input_parser &);//copy
  };




#endif
