//
// C++ Implementation: na
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

#include "na.hpp"
#include "types.hpp"

Real sucbis(Real f,const Real & aa, const Real & bb,const Real & eps)
{
 Real a=aa;
 Real fa=f(a);
 Real b=bb;
 Real fb=f(b);
 //cout << fa << " " << fb << endl;
 if(fa*fb<0)
    {
      while(1)
      {
        Real c=.5*(a+b);
        Real fc=f(c);
        //cout << a << " " << b << " " << b-a <<" " << c << " " << fc << endl;
        if((fc==0)||(b-a<eps))
        {
	  return c;
        }
        else if(fa*fc<0)
        {
          b=c;
          fb=fc;
        }
        else if(fb*fc<0)
        {
          a=c;
          fa=fc;
        }
        else
        {
        cout << "Successive bisection method failed" << endl;
          exit(1);
        }    
      }  
}
else
 {
        cout << "Successive bisection method failed" << endl;
          exit(1);
        }    

}