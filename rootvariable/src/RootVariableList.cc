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
#include "RootVariableList.hh"
#include "TChain.h"
#include "TLeaf.h"
#include <iostream>
#include <cstdarg>
#include <assert.h>
#include <cmath>
#include "RootVariableFormula.hh"
#include "RootVariableFunction.hh"
#include "TTreeFormulaManager.h"
#include "RootVariableManager.hh"

using namespace qosc;

RootVariableList::RootVariableList() {
  // empty constructor
  m_source_tree = NULL;
}

RootVariableList::RootVariableList( std::string variables, TChain* source_tree ) {
  m_source_tree = source_tree;
  m_var_names.clear();
  ParseVariableList( variables, m_var_names );
  CreateRootVariables( m_var_names );
}

RootVariableList::~RootVariableList() {
  // Doesn't own any of the instances, so nothing to do.
}

void RootVariableList::SetChain( TChain* source_tree ) {
  if ( m_source_tree && source_tree!=m_source_tree ) {
    std::cout << "Trying to reset the chain! Previous Chain=" << m_source_tree << " new chain=" << source_tree << std::endl;
    assert(false);
  }
  m_source_tree = source_tree;
}

void RootVariableList::ParseVariableList( std::string variablenames, std::vector< std::string >& var_names ) {
  // separates list of variables divided by ;:,
  size_t pos = 0;
  while ( pos!=std::string::npos ) {
    size_t nextpos = variablenames.find_first_of(";:,",pos);
    if ( nextpos!=std::string::npos ) {
      // found variable name
      std::string varname = variablenames.substr( pos, nextpos-pos );
      if ( varname!="" ) var_names.push_back( varname );
      pos = nextpos+1;
    }
    else {
      if ( variablenames.size()>0 && pos!=std::string::npos ) {
	// last variable name
	std::string varname = variablenames.substr( pos, nextpos);
	if (varname!="") var_names.push_back( varname );
	pos = std::string::npos;
      }      
      pos = nextpos;
    }
  }

}

// ---------------------------------------------------------------------------
// Append new variables

int RootVariableList::Add( double (*VariableFunction)(RootVariableList&), std::string function_variables, TChain* source_tree ) {
  // Add function pointer
  RootVariableFunction* varfunc = RootVariableManager::GetTheRootVarManager()->MakeRootVariableFunction( VariableFunction, function_variables, source_tree );
  m_var_idmap[ varfunc->GetID() ] = varfunc;
  m_var_map[ varfunc->GetVariableName() ] = varfunc;
  return m_var_map.size();
}

int RootVariableList::Add( std::string varname, double (*VariableFunction)(RootVariableList&), std::string function_variables, TChain* source_tree ) {
  // add function pointer
  RootVariableFunction* varfunc =  RootVariableManager::GetTheRootVarManager()->MakeRootVariableFunction( varname, VariableFunction, function_variables, source_tree );
  m_var_idmap[ varfunc->GetID() ] = varfunc;
  //m_var_map[ varfunc->GetVariableName() ] = varfunc;
  m_var_map[ varname ] = varfunc; // default name
  return m_var_map.size();
}

int RootVariableList::Add( std::string varnames ) {
  // Add formula
  if ( m_source_tree==NULL) {
    std::cout << "Attempted to add a root variable to a variable list without specifying an initial source ROOT tree." << std::endl;
    assert(false);
  }

  std::vector< std::string > varname_list;
  ParseVariableList( varnames, varname_list );
  CreateRootVariables( varname_list );

  return m_var_map.size();
}

int RootVariableList::Add( RootVariable* variable ) {
  std::string variable_name = variable->GetVariableName();
  m_var_idmap[ variable->GetID() ] = variable;
  m_var_map[variable_name] = variable; 
  return m_var_map.size(); 
}

// ---------------------------------------------------------------------------
// Utility

void RootVariableList::CreateRootVariables( std::vector< std::string >& var_names ) {
  
  if ( GetChain()==NULL ) {
    std::cout << "Creating root variables before specifying chain for list!" << std::endl;
    assert(false);
  }
  for (unsigned int i=0; i<var_names.size(); i++) {
    std::string varname = var_names.at(i);
    RootVariable* varbase = (RootVariable*)RootVariableManager::GetTheRootVarManager()->MakeRootVariableFormula( varname, GetChain() );
    if (varbase==NULL) {
      std::cout << "Root variable creation produced NULL pointer!" << std::endl;
      std::cout << "  variable name: " << varname << std::endl;
      std::cout << "  source_tree: " << m_source_tree << std::endl;
      assert(false);
    }
    // store the variable in the root list
    m_var_idmap[varbase->GetID()] = varbase;
    m_var_map[varname] = varbase;
    //std::cout << "storing " << varname << ": " << varbase << ", nelemnts=" << m_var_map.size() << std::endl;
  }
  
}

// ---------------------------------------------------------------------------
// Access Functions

RootVariable* RootVariableList::GetVariable( std::string varname ) {
  //std::cout << "Retrieving " << variable;
  //std::cout << ": " << m_var_map[variable] << std::endl;
  std::map< std::string, RootVariable*>::iterator it_rootvar = m_var_map.find(varname);
  if (it_rootvar==m_var_map.end() ) {
    std::cout << "Asked for a rootvariable, " << varname << " not found in list" << std::endl;
    std::cout << "Variables in this list include:" << std::endl;
    for (it_rootvar=m_var_map.begin(); it_rootvar!=m_var_map.end(); it_rootvar++)
      std::cout << " " << (*it_rootvar).first << std::endl;
    assert(false);
  }
  return  (*it_rootvar).second;
}

RootVariable* RootVariableList::GetVariable( int variableid ) {
  //std::cout << "Retrieving " << variable;
  //std::cout << ": " << m_var_map[variable] << std::endl;
  std::map< int, RootVariable*>::iterator it_rootvar = m_var_idmap.find(variableid);
  if (it_rootvar==m_var_idmap.end() ) {
    std::cout << "Asked for a rootvariable with ID=" << variableid << " not found in list" << std::endl;
    std::cout << "Variables in this list include:" << std::endl;
    for (it_rootvar=m_var_idmap.begin(); it_rootvar!=m_var_idmap.end(); it_rootvar++)
      std::cout << " " << (*it_rootvar).first << " ID=" << (*it_rootvar).second->GetID() << std::endl;
    assert(false);
  }
  return (*it_rootvar).second;
}

RootVariable* RootVariableList::GetVariable( double (*VariableFunction)(RootVariableList&), std::string function_variables ) {
  /// I don't like this one bit.
  std::string varfuncname = RootVariableManager::GetTheRootVarManager()->GetCanonicalFunctionName( VariableFunction, function_variables, m_source_tree );
  return GetVariable( varfuncname );
}

void RootVariableList::GetNameList( std::vector< std::string >& varlist ) {
  std::map< std::string, RootVariable* >::iterator iter;
  for (iter=m_var_map.begin(); iter!=m_var_map.end(); iter++) {
    varlist.push_back( (*iter).first );
  }
}


double RootVariableList::GetVariableValue( RootVariable* varbase, int numdims, ... ) {

  double value = 0;
  if ( varbase==NULL ) {
    std::cout << "Somehow the variable address is NULL" << std::endl;
    assert(false);
  }
  std::string vartype = varbase->GetType();

  // unneccessary string check: can optimize by removing.
  if (vartype=="Char_t") {
    std::cout << "Variable is of type Char_t.  Call GetStrVariableValue instead." << std::endl;
    assert(false);
  }
  
  va_list args;
  va_start(  args, numdims );
  value =  varbase->Value(numdims, args);
  va_end( args );
  
  return value;

}

double RootVariableList::GetVariableValue( std::string varname, int numdims, ... ) {

  RootVariable* variable = GetVariable( varname );

  va_list args;
  va_start( args, numdims );
  double value =  GetVariableValue( variable, numdims, args);
  va_end( args );
  return value;
}

double RootVariableList::GetVariableValue( int varid, int numdims, ... ) {

  RootVariable* variable = GetVariable( varid );

  va_list args;
  va_start( args, numdims );
  double value =  GetVariableValue( variable, numdims, args);
  va_end( args );
  return value;
}

long RootVariableList::GetVariableValueInt( RootVariable* variable, int numdims, ... ) {

  // TTreeFormula always gives back a double. But I am worried about rounding
  va_list args;
  va_start( args, numdims );
  double value = GetVariableValue( variable, numdims );
  va_end( args );
  long lvalue = lrint( value );
  return lvalue;
}

long RootVariableList::GetVariableValueInt( std::string varname, int numdims, ... ) {
  RootVariable* var = GetVariable( varname );
  va_list args;
  va_start( args, numdims );
  long value = GetVariableValueInt( var, numdims, args );
  va_end( args );
  return value;
}

long RootVariableList::GetVariableValueInt( int varid, int numdims, ... ) {
  RootVariable* var = GetVariable( varid );
  va_list args;
  va_start( args, numdims );
  long value = GetVariableValueInt( var, numdims, args );
  va_end( args );
  return value;
}

void RootVariableList::Clear() {
  m_var_idmap.clear();
  m_var_map.clear();
}
