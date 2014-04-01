/* 
   ------------------------------------------------------------------------------------------------
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
 * \class RootVariable
 * \ingroup RootVariable
 * \brief Abstract base class which defines the interface to ROOT variable types.
 *
 * It provides methods to help access ROOT data without the tedious and messy boiler plate.
 * Child classes that provide concrete instances
 *  (1) RootVariableFormula: based on TTreeFormula
 *  (2) RootVariableFunction: uses RootVariableFormulas + user-defined function pointer
 * Note that in both cases the variable represented by this object can be a function.
 *
 * Typically, one should create instances of this class only through
 * the RootVariableList object.  This class acts as container/factory and -- importantly --  registers 
 * the created instances to the RootVarManager, which is a singleton class that acts a cache 
 * that serves to avoid recalculating the variable value.
 *
 * One can refuse to use the cache by calling RootVarManager::ForbidCache()
 * 
 * -----------------------------------------------------------------------------------------------
 */

#ifndef __RootVariable__
#define __RootVariable__

#include "TChain.h"
#include <cstdarg>
#include <string>

namespace qosc {

  class RootVariable {

  public:
  
    RootVariable();
    RootVariable( std::string varname, TChain* source );
    virtual ~RootVariable();

    std::string GetVariableName() { return m_variable_name; };
  
    // Tree Management Calls
    // common both to formula and function, should do here
    long GetEntry(unsigned long entry);  // Loads the Tree entry number [entry]

    // Variable Value Calls
    virtual double Value() { return Value(1); }; ///< Shortcut for single l-value variable
    virtual double Value( int ndims, ... );  ///< Call to get l-value variables, wraps next function
    virtual double Value( int ndims, va_list args ) = 0; // gets the value of the variable or formula.  number of arguments depends on dimensions of data member (for container call)
    //virtual double String( int ndims, ... );  ///< Direct call to get either a string variable or a print out of the variable value for l-value variables
    //virtual double Value( int ndims, va_list args ) = 0; // gets the value of the variable or formula.  number of arguments depends on dimensions of data member (for container call)
    virtual std::string PrintValue() = 0; // prints the value as a string. note, this is the way to get text items from tree

    // Variable Info Calls
    bool IsVariableAnArray() { return kDataIsArray; };
    std::string GetType() { return m_variable_type; };

    static int GetTotalVariableInstances() { return RootVariable::m_var_id; };

  protected:

    // Internal setup calls
    void SetVariableName( std::string varname ) { m_variable_name = varname; };
    void SetArrayFlag( bool isarray ) { kDataIsArray = isarray; };
    bool HasTreeIDChanged();
    void SetType( std::string type ) { m_variable_type = type; };  
    TChain* GetChain() { return m_source_tree; };
    void SetChain( TChain* chain ) { m_source_tree = chain; if (m_source_tree) m_tree_id = m_source_tree->GetTreeNumber(); };

  private:
  
    bool kDataIsArray; // is the data member an array?
    std::string m_variable_type; // The variable's type
    std::string m_variable_name; // the variable's name
    int m_tree_id; // The tree number in our chain
    TChain* m_source_tree; // The ROOT Tree Chain that is the source of the data

    static int m_var_id;

  };
  
}

#endif
