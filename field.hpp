// -*- C++ -*-

//
// C++ Interface: field
//
// Description: 
//
//
// Author: Manuel Baptista <mbaptist@gmail.com>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//

/*
Copyright (C) 2008 Manuel Baptista <mbaptist@gmail.com>

This file is part of eke.

Eke is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Eke is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with eke.  If not, see <http://www.gnu.org/licenses/>.
*/

#IFNDEF FIELD_HPP
#DEFINE FIELD_HPP

template <class GRID,class VALUE>
Field
{
private:
	GRID & grid_;
	VALUE & value_;
public:
	//Accessors
	GRID & grid(){return grid_;};
	VALUE & value(){return value_;};
	const GRID & grid() const {return grid_;};
	const VALUE & value() const {return value_;};
public:
	//Ctor
	Field(GRID & _grid_,VALUE & _value_):
		grid_(_grid_),
		value_(_value_){};
	//Dtor
	~Field(){};
private:
	//Forbidden Ctors
	Field();
};


#ENDIF