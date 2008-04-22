// -*- C++ -*-
//
// //
// // C++ Implementation: eke_python
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

#include "eke.hpp"

#include <string>

#include <iostream>

#include <complex>

#include <cat.h>


#include <goops.h>

using namespace cat;

using namespace std;

//Launch the lss code
PyObject *
mhdc3dl_python_lss_run(PyObject *self, PyObject *args)
{
  Real theta_min,theta_max;
  std::complex<Real> lambda_min,lambda_max;
  char * input_module_name;
  PyObject * dict;
  if (!PyArg_ParseTuple(args, "O", &dict ))
    return NULL;
  Input input_obj(dict);
  lss lss_obj(input_obj);
  lss_obj.run(theta_min,lambda_min,theta_max,lambda_max);
  /*  Py_complex plmin,plmax;
    plmin.real=lambda_min.real();
    plmin.imag=lambda_min.imag();
    plmax.real=lambda_max.real();
    plmax.imag=lambda_max.imag();*/
  return Py_BuildValue("dddddd",theta_min,lambda_min.real(),lambda_min.imag(),theta_max,lambda_max.real(),lambda_max.imag());
};

//Launch the sss code
PyObject *
mhdc3dl_python_sss_run(PyObject *self, PyObject *args)
{
  Real xp,eim;
  char * input_module_name;
  PyObject * dict;
  if (!PyArg_ParseTuple(args, "O", &dict ))
    return NULL;
  Input input_obj(dict);
  sss sss_obj(input_obj);
  sss_obj.run(xp,eim);
  return Py_BuildValue("dd",xp,eim);
};

//Add Methods to function
PyMethodDef mhdc3dl_python_Methods[] = {
                                         {"lss_run", mhdc3dl_python_lss_run , METH_VARARGS,
                                          "Run large scale stability for a magnetoconvective system"},
                                         {"sss_run", mhdc3dl_python_sss_run , METH_VARARGS,
                                          "Run short scale stability for a magnetoconvective system"},
                                         {NULL, NULL, 0, NULL}        /* Sentinel */
                                       };

PyMODINIT_FUNC
initmhdc3dl_python(void)
{
  (void) Py_InitModule("mhdc3dl_python", mhdc3dl_python_Methods);
};

