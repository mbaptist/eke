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

#include "io.hpp"

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;

#include "grid.hpp"


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

int PyInputParser::parse_int(const std::string & item)
{
  //val=PyInt_AsLong(PyDict_GetItem(dict,PyString_FromString(item.c_str())));
  return PyInt_AsLong(pval(item));
}
bool PyInputParser::parse_bool(const std::string & item)
{
  //val=PyInt_AsLong(PyDict_GetItem(dict,PyString_FromString(item.c_str())));
  return PyInt_AsLong(pval(item));
}
double PyInputParser::parse_double(const std::string & item)
{
  // val=PyFloat_AsDouble(PyDict_GetItem(dict,PyString_FromString(item.c_str())));
  return PyFloat_AsDouble(pval(item));
}

std::string PyInputParser::parse_string(const std::string & item)
{
  //  val=PyString_AsString(PyDict_GetItem(dict,PyString_FromString(item.c_str())));
  return PyString_AsString(pval(item));
}

std::vector<int> PyInputParser::parse_vector_int(const std::string & item)
{
	vector<int> cppvector;
	PyObject * pyvector;
	pyvector=pval(item);
	for(int n=0;n<PyList_Size(pyvector);++n)
		cppvector.push_back(PyInt_AsLong(PyList_GetItem(pyvector,n)));
	return cppvector;		
}

std::vector<double> PyInputParser::parse_vector_double(const std::string & item)
{
	vector<double> cppvector;
	PyObject * pyvector;
	pyvector=pval(item);
	for(int n=0;n<PyList_Size(pyvector);++n)
		cppvector.push_back(PyFloat_AsDouble(PyList_GetItem(pyvector,n)));
	return cppvector;		
}

PyObject * PyInputParser::pval(const std::string & item)
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




#if 1

void vtkSave(const std::string & filename,const RSF & field,const std::string & fieldname,const Grid & grid)
{
  ofstream ofs((filename+".vtk").c_str());
  ofs << "# vtk DataFile Version 2.0 \n"
    << ". \n"
    << "ASCII \n"
    << "DATASET STRUCTURED_POINTS \n"
    << "DIMENSIONS " << field.shape()[0] << " " << field.shape()[1] << " " << field.shape()[2] << "\n"
    << "ORIGIN 0 0 0 \n "
    << "SPACING " << grid.deltax() << " " << grid.deltay() << " " << grid.deltaz() << "\n"
    << "POINT_DATA " << field.size() << "\n"
    << "SCALARS " << fieldname.c_str() << " float\n"
    << "LOOKUP_TABLE default" << endl;
  for (int n=0;n<field.size();++n)
    ofs << field.data()[n] << endl;
}

void vtkSave(const std::string & filename,const RVF & field,const std::string & fieldname,const Grid & grid)
{
  ofstream ofs((filename+".vtk").c_str());
  ofs << "# vtk DataFile Version 2.0 \n"
    << ". \n"
    << "ASCII \n"
    << "DATASET STRUCTURED_POINTS \n"
    << "DIMENSIONS " << field.shape()[0] << " " << field.shape()[1] << " " << field.shape()[2] << "\n"
    << "ORIGIN 0 0 0 \n "
    << "SPACING " << grid.deltax() << " " << grid.deltay() << " " << grid.deltaz() << "\n"
    << "POINT_DATA " << field.size() << "\n"
    << "VECTORS " << fieldname.c_str() << " float" << endl;
  for (int n=0;n<field.size();++n)
    {
      for (int m=0;m<3;++m)
        ofs << field.data()[n][m] <<" ";
      ofs << endl;
    }
}

#endif
