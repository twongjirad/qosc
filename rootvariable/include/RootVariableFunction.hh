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
/**
 * ------------------------------------------------------------------------------------------------
 * \class RootVariableFormula
 * \ingroup rootvariable
 * \brief This class stores a value that is calculated using data from a ROOT Tree and a user function pointer.
 *
 * Instances of this template class are suppose to be a representation of numbers in a ROOT tree.
 * However, unlike the RootVariableFunction class which attaches a block of memory to the ROOT tree's branch,
 * this class uses the TTreeFormula Class to evaluate a number for each entry.
 * Pass function variable of this form:
 * typdef (*VariableFunction) ( RootVariableList* ); // The prototypical function
 *
 * \note This class could really use a good cleaning. (Used to be a template class and has cruft that shows that fact.)
 *
 * -----------------------------------------------------------------------------------------------
 */

#ifndef __RootVariableFunction__
#define __RootVariableFunction__

#include <string>
#include <cstdarg>
#include "RootVariable.hh"
#include "RootVariableList.hh"

class TChain;

namespace qosc {

  class RootVariableFunction : public RootVariable {

    //typedef double (*VariableFunction) ( RootVariableList* );

  public:
  
    RootVariableFunction( std::string variable_name, double (*VariableFunction)(RootVariableList&), std::string root_vars, TChain* source_tree );
    RootVariableFunction( std::string copyname, RootVariableFunction* orig );
    virtual ~RootVariableFunction();
    virtual double Value( int ndims, va_list args ); // gets the value of the variable or formula.  number of arguments depends on dimensions of data member
    void SetVariableFunction( double (*VariableFunction)(RootVariableList&) );
    void SetVariableFunctionArguments( std::string funcvars ) { m_function_variables = funcvars; };
    void* GetVariableFunction() { return (void*)m_variable_function; };
    void CopyVariableFunctionPtr( void* funcptr ) { funcptr = (void*)m_variable_function; };
    std::string GetRootVariablesString() { return m_function_variables; };
    virtual std::string PrintValue(); // prints the value as a string. note, this is the way to get text items from tree

  protected:
  
  
  private:

    std::string m_function_variables;
    RootVariableList m_function_rootvars; // Root Variables for our function to work
    double (*m_variable_function)(RootVariableList&); // Root Function Pointer

    int lastCalculatedEntry;
    double cachedValue;
    bool useCachedValue;
  
  };

}
#endif
