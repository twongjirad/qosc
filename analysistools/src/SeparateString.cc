//-*- mode:c++; c-basic-offset:2;   -*-
#include "SeparateString.hh"
#include <iostream>

void SeparateString( std::string strlist, std::vector< std::string >& list, std::string separators ) {
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
