// -*- C++ -*-
//
// //
// // C++ Implementation: input
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



#include <Python.h>

#include "input.hpp"

#include <iostream>
#include <sstream>
#include <string>
using namespace std;


///////////////////////////////////////////
// Implementation of class PyInputParser
///////////////////////////////////////////

PyInputParser::PyInputParser(string filename_):
    filename(filename_),
    python_initialised(Py_IsInitialized())
{
  Py_Initialize();
  stringstream ss;
  ss << "import sys\n"
  << "sys.path.insert(0,'')\n"
  << "import imp\n"
  << "imp.load_source('input','"
  << filename
  << "')\n";
  PyRun_SimpleString(ss.str().c_str());
  module=PyImport_ImportModule("input");
  dict=PyModule_GetDict(module);
}

PyInputParser::PyInputParser(PyObject * dict_):
    python_initialised(Py_IsInitialized())
{
  dict=dict_;
}

PyInputParser::~PyInputParser()
{
  if(!python_initialised)
    {
      Py_DECREF(dict);
      Py_DECREF(module);
      PyRun_SimpleString("if(imp.lock_held()):imp.release_lock()\n");
      Py_Finalize();
    }
}

int PyInputParser::parse_int(const string & item)
{
  //val=PyInt_AsLong(PyDict_GetItem(dict,PyString_FromString(item.c_str())));
  return PyInt_AsLong(pval(item));
}
bool PyInputParser::parse_bool(const string & item)
{
  //val=PyInt_AsLong(PyDict_GetItem(dict,PyString_FromString(item.c_str())));
  return PyInt_AsLong(pval(item));
}
double PyInputParser::parse_double(const string & item)
{
  // val=PyFloat_AsDouble(PyDict_GetItem(dict,PyString_FromString(item.c_str())));
  return PyFloat_AsDouble(pval(item));
}

std::string PyInputParser::parse_string(const string & item)
{
  //  val=PyString_AsString(PyDict_GetItem(dict,PyString_FromString(item.c_str())));
  return PyString_AsString(pval(item));
}

PyObject * PyInputParser::pval(const string & item)
{
  PyObject * pyvalue=PyDict_GetItem(dict,PyString_FromString(item.c_str()));
  if (pyvalue!=0)
    {
      return(pyvalue);
    }
  else
    {
      cout << item
      << " not defined in "
      << filename
      << " . Aborting... "
      << endl;
      exit(1);
    }
}
