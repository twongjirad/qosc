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

#include <iostream>
#include <sstream>
#include <assert.h>
#include <cstdlib>
#include <cstdarg>
#include <sstream>

#include "TChain.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TTreeFormula.h"

#include "RootVariableFunction.hh"
#include "RootVariable.hh"
#include "RootVariableList.hh"

using namespace qosc;

RootVariableFunction::RootVariableFunction( std::string variable_name, double (*VariableFunction)(RootVariableList&), std::string root_vars, TChain* source_tree ) 
  : RootVariable( variable_name, source_tree ) {
  
  // initialize internal members
  RootVariable::SetVariableName( variable_name );
  RootVariable::SetArrayFlag( false ); //this is hacky -- I am going up one class level
  // save the address of the source tree
  RootVariable::SetChain( source_tree );

  // set function
  m_function_variables = root_vars;
  RootVariableFunction::SetVariableFunction( VariableFunction );
  // Make Variable List
  m_function_rootvars.SetChain( RootVariable::GetChain() );
  m_function_rootvars.Add( m_function_variables );

  // Set Type: Need to translate C++ type to ROOT Type (Note: this might not be machine portable)
  RootVariable::SetType( "Double_t" );

  lastCalculatedEntry = -1;
  useCachedValue = false;
  cachedValue = 0.;
}

RootVariableFunction::RootVariableFunction( std::string copyname, RootVariableFunction* orig ) 
  : RootVariable( copyname, NULL )
{
  // initialize internal members
  RootVariable::SetArrayFlag( false ); //this is hacky -- I am going up one class level
  // save the address of the source tree
  RootVariable::SetChain( orig->GetChain() );
  // copy type
  RootVariable::SetType( orig->GetType() );
  // Copy function pointer
  orig->CopyVariableFunctionPtr( (void*)m_variable_function );
  //RootVariableFunction::SetVariableFunction( orig->GetVariableFunction() );
  // copy function variables
  RootVariableFunction::SetVariableFunctionArguments( orig->GetRootVariablesString() );

}

RootVariableFunction::~RootVariableFunction() {};

void RootVariableFunction::SetVariableFunction( double (*VariableFunction)(RootVariableList&) ) 
{ m_variable_function = VariableFunction; };

double RootVariableFunction::Value( int ndims, va_list arg_list ) {
  
  // For single value data memebers
  if (RootVariable::IsVariableAnArray()==false) {

    if ( useCachedValue==true ) {
      if ( GetChain()->GetReadEntry()!=lastCalculatedEntry ) {
	lastCalculatedEntry = GetChain()->GetReadEntry();
	cachedValue = (*m_variable_function)( m_function_rootvars );
      }
      return cachedValue;	
    }
    else {
      return (*m_variable_function)( m_function_rootvars );
    }
  }
  
  std::cout << "Calling for the value of " << RootVariable::GetVariableName() << " with the wrong dimensions. Only non-array setup for now." << std::endl;
  assert(false);
 
}

std::string RootVariableFunction::PrintValue() {
  std::stringstream value;
  value << Value(1,0);
  return value.str();
}
