//-*- mode:c++; c-basic-offset:2;   -*-
/*  ------------------------------------------------------------------------------------------------

   Copyright (C) 2014  Taritree Wongjirad

   This file is part of qosc.

   Qosc is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   
   Qosc is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with Qosc.  If not, see <http://www.gnu.org/licenses/>.
   ---------------------------------------------------------------------------------------------- */

#include "SeparateString.hh"
#include <iostream>

void qosc::SeparateString( std::string strlist, std::vector< std::string >& list, std::string separators ) {
  size_t pos = 0;
  
  while ( pos<std::string::npos ) {
    size_t nextpos = strlist.find_first_of( separators, pos);
    if ( nextpos!=std::string::npos || pos<std::string::npos ) {
      std::string sysname;

      if ( nextpos!=std::string::npos ) {
	sysname = strlist.substr(pos,nextpos-pos); // found separating string
	pos = nextpos+1;
      }
      else {
	sysname = strlist.substr( pos ); // read off last sys term
	pos = nextpos; // npos
      }

      if (sysname!="") {      
	list.push_back( sysname );
      }
    }
    else {
      pos = nextpos;
    }
    // std::cout << "pos=" << pos << " nextpos=" << nextpos << " npos=" << std::string::npos << std::endl; // for debug
  }
  
}
