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
#include <string>
#include <vector>
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



#if 0
 //Save charges
ofstream ofs("charges.dat");
  //ofs << charge_density << endl;
for (int n=0;n<charges.size();++n)
	ofs << charges.data()[n] << endl;
   //cout << charge_density << endl;


  //Save initial electric field
	ofstream ofs2("efi.dat");
  //ofs2 << electric_field << endl;
	for (int n=0;n<electric_field.size();++n)
	{
		for (int m=0;m<3;++m)
			ofs2 << electric_field.data()[n][m] <<" ";
		ofs2 << endl;
	}



  //write electric field after minimisation  
	ofstream ofs3("ef.dat");
	ofs3 << fixed << setprecision(16);
  //ofs3 << electric_field << endl;
	for (int n=0;n<electric_field.size();++n)
	{
		for (int m=0;m<3;++m)
			ofs3 << electric_field.data()[n][m] <<" ";
		ofs3 << endl;
	}



#endif
