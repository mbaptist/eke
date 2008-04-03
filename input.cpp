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

#include "input.h"

#include <string>
#include <sstream>
#include <fstream>

using namespace std;

///////////////////////////////////////////
// Definition of class Input
///////////////////////////////////////////


//constructor from filename
Input::Input(std::string cfg_fname)
{
  py_input_parser parse(cfg_fname);
  py_load_Input(parse);
}
//constructor from python object
Input::Input(PyObject * module)
{
  py_input_parser parse(module);
  py_load_Input(parse);
}
//destructor
Input::~Input()
{}
;


//load from config file
void Input::py_load_Input(py_input_parser & parse)
{

  parse(runsname,"runsname");
  cout << runsname << endl;

  parse(n1,"n1");
  parse(n2,"n2");
  parse(n3,"n3");

  parse(l1,"l1");
  parse(l2,"l2");
  parse(l3,"l3");

  parse(ell1,"ell1");
  parse(ell2,"ell2");

  parse(visc,"visc");
  parse(omegaz,"omegaz");
  parse(compress,"compresss");

  cout << visc << " " << compress << endl;

  parse(g,"g");
  parse(diff,"diff");
  parse(deltat,"deltat");
  parse(econd,"econd");
  parse(tcond,"tcond");

  //refine and resume
  parse(refine,"refine");
  parse(resume,"resume");
  parse(lr_runsname,"lr_runsname");
  parse(lr_n1,"lr_n1");
  parse(lr_n2,"lr_n2");
  parse(lr_n3,"lr_n3");

  //Basic Fields
  parse(basic_mode,"basic_mode");
  //Random
  parse(br_spectrum,"br_spectrum");
  parse(br_seed,"br_seed");
  parse(br_ki,"br_ki");
  parse(br_kf,"br_kf");
  parse(br_alpha,"br_alpha");
  parse(br_rms_norm,"br_rms_norm");
  parse(br_kind,"br_kind");
  parse(br_sym,"br_sym");

  parse(ls_eps,"ls_eps");
  parse(qq,"qq");
  parse(qq_adj,"qq_adj");
  parse(kk,"kk");
  parse(small,"small");
  parse(small_adj,"small_adj");

  parse(ep,"ep");
  parse(thr,"thr");
  parse(mp,"mp");
  parse(sc,"sc");
  parse(nseq,"nseq");

  parse(sss_ifname,"sss_ifname");
  parse(sss_int_ofbname,"sss_int_ofbname");
  parse(sym_sub,"sym_sub");
  parse(sss_seed,"sss_seed");

}



///////////////////////////////////////////
// Definition of class py_input_parser
///////////////////////////////////////////

py_input_parser::py_input_parser(string filename_):
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

py_input_parser::py_input_parser(PyObject * dict_):
    python_initialised(Py_IsInitialized())
{
  dict=dict_;
}

py_input_parser::~py_input_parser()
{
  if(!python_initialised)
    {
      Py_DECREF(dict);
      Py_DECREF(module);
      PyRun_SimpleString("if(imp.lock_held()):imp.release_lock()\n");
      Py_Finalize();
    }
}

void py_input_parser::operator()(int & val, const string & item)
{
  //val=PyInt_AsLong(PyDict_GetItem(dict,PyString_FromString(item.c_str())));
  val=PyInt_AsLong(pval(item));
}
void py_input_parser::operator()(bool & val, const string & item)
{
  //val=PyInt_AsLong(PyDict_GetItem(dict,PyString_FromString(item.c_str())));
  val=PyInt_AsLong(pval(item));
}
void py_input_parser::operator()(double & val, const string & item)
{
  // val=PyFloat_AsDouble(PyDict_GetItem(dict,PyString_FromString(item.c_str())));
  val=PyFloat_AsDouble(pval(item));
}

void py_input_parser::operator()(string & val, const string & item)
{
  //  val=PyString_AsString(PyDict_GetItem(dict,PyString_FromString(item.c_str())));
  val=PyString_AsString(pval(item));
}

PyObject * py_input_parser::pval(const string & item)
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
