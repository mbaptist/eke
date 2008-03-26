#Makefile

#Copyright (C) 2008 Manuel Baptista <mbaptist@gmail.com>
#
#This file is part of eke.
#
#Eke is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#Eke is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with eke.  If not, see <http://www.gnu.org/licenses/>.


CXX = g++

DEBUG =
#-g

###############################################################################

IFLAGS = $(DEBUG) $(INCLUDE)
FLAGS = $(DEBUG) $(LIB)

OBJECTS = eke.o grid.o maggs.o

###############################################################################

all: eke

clean:
	@rm -rfv eke *.o *.so

distclean: clean
	@rm -rfv *~

%.o: %.cpp *.hpp
	$(CXX) $(IFLAGS) -c $<

eke: $(OBJECTS) 
	$(CXX) -o eke $(OBJECTS)

