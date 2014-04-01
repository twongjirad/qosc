/* ------------------------------------------------------------------------------------------------

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

#include "RootVariable.hh"
#include <sstream>

using namespace qosc;

int RootVariable::m_var_id = 0;

RootVariable::RootVariable() {
  m_instance_id = m_var_id;
  std::stringstream ss;
  ss.str("");
  ss << "var_" << m_instance_id;
  SetVariableName( ss.str() );
  SetChain( NULL );
  m_var_id++;
}

RootVariable::RootVariable( std::string varname, TChain* source ) {
  m_instance_id = m_var_id;
  std::stringstream ss;
  ss.str("");
  ss << "var_" << m_instance_id << "_" << varname;
  SetVariableName( ss.str() );
  SetChain( source );
  m_var_id++;
}
 
RootVariable::~RootVariable() {}


double RootVariable::Value( int ndims, ... ) {
  va_list args;
  double var_value = 0;
  va_start( args, ndims );
  var_value = Value( ndims, args ); // Call to concrete class's value function
  va_end( args );

  return var_value;
}

long RootVariable::GetEntry( unsigned long entry ) {
  // returns the bytes read. if zero, then no such entry number exists
  return GetChain()->GetEntry(entry);
}

bool RootVariable::HasTreeIDChanged() {
  int current_id = GetChain()->GetTreeNumber();
  if ( current_id != m_tree_id ) {
    m_tree_id = current_id;
    return true;
  }
  return false;
}
