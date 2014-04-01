#include "RootVariable.hh"
#include <sstream>

int RootVariable::m_var_id = 0;

RootVariable::RootVariable() {
  std::stringstream ss;
  ss.str("");
  ss << "var_" << m_var_id;
  SetVariableName( ss.str() );
  SetChain( NULL );
  m_var_id++;
}

RootVariable::RootVariable( std::string varname, TChain* source ) {
  std::stringstream ss;
  ss.str("");
  ss << "var_" << m_var_id << "_" << varname;
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
