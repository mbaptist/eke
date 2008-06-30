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

#ifndef IO_HPP
#define IO_HPP

#include <Python.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "grid.hpp"

///////////////////////////////////////////
// Declaration of class PyInputParser
///////////////////////////////////////////
class PyInputParser
  {
  private:
    //Filename to parse
    std::string filename;
    //Python objects
    PyObject * module;
    PyObject * dict;
    PyObject * pval(const std::string & item);
    PyThreadState * inter;
    bool python_initialised;
    bool state;
  public:
    //Constructor from filename
    PyInputParser(std::string filename_);
    //Constructor from dict
    PyInputParser(PyObject * module_);
    //Destructor
    ~PyInputParser();
    //Methods
    //void open();
    void close();
    //Python to C++ conversion function
    int parse_int(const std::string & item);
    bool parse_bool(const std::string & item);
    double parse_double(const std::string & item);
    std::string parse_string(const std::string & item);
    std::vector<int> parse_vector_int(const std::string & item);
    std::vector<double> parse_vector_double(const std::string & item);
  private:
    //Forbidden constructors
    PyInputParser();//default
    PyInputParser(const PyInputParser &);//copy
  };



void vtkSave(const std::string & filename,const RSF & field,const std::string & fieldname,const Grid & grid);

void vtkSave(const std::string & filename,const RVF & field,const std::string & fieldname,const Grid & grid);


#endif
