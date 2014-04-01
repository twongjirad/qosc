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
 * ------------------------------------------------------------------------------------
 * \class RootVariableList
 * \ingroup rootvariable
 * \brief This is class is the container and object factory for RootVariables.  
 *
 * Its centered around a stl map object of type <string, RootVariable*>
 * Gives quick routines to declare and destroy variables.
 * The variable instances are all stored in the singleton RootVariableManager however.
 * They live there for optimization reasons.
 *
 * ------------------------------------------------------------------------------------
 */

#ifndef __RootVariableList__
#define __RootVariableList__

#include <vector>
#include <string>
#include <map>

class TChain;

namespace qosc {
  
  class RootVariable;

  typedef std::map< std::string, RootVariable* >::iterator RootVariableListIter;

  class RootVariableList {

  public:

    RootVariableList();
    RootVariableList( std::string variables, TChain* source_tree ); // provide a comma or colon sepearated list of formulas to load
    ~RootVariableList();

    void SetChain( TChain* source_tree ); // Need to make sure the formula leaves' chain gets updated properly.  this command is not enough.
    int GetEntries() { return m_var_map.size(); };
    int Add( std::string variables ); // add via formula
    int Add( RootVariable* variable );
    int Add( double (*VariableFunction)(RootVariableList&), std::string function_variables, TChain* source_tree );
    int Add( std::string variable_name, double (*VariableFunction)(RootVariableList&), std::string function_variables, TChain* source_tree );
    void CreateRootVariables(std::vector< std::string >& var_names ); // create a list of root variables by formula

    // Provide two methods for accessing variable:

    // By name for humans
    double GetVariableValue( std::string varname, int ndims=1, ... );
    long GetVariableValueInt( std::string varname, int ndims=1, ... );
    //double GetStrVariableValue( std::string varname, int ndims=1, ... );

    // By ID for machines
    double GetVariableValue( int varid, int ndims=1, ... );
    long GetVariableValueInt( int varid, int ndims=1, ... );
    //double GetStrVariableValue( std::string varname, int ndims=1, ... );

    RootVariable* GetVariable( int varid); // Get a variable object in the list by name
    RootVariable* GetVariable( std::string varname); // Get a variable object in the list by name
    RootVariable* GetVariable( double (*VariableFunction)(RootVariableList&), std::string function_variables );
    void GetNameList( std::vector< std::string >& varlist );
    void Clear();
    RootVariableListIter Begin() { return m_var_map.begin(); };
    RootVariableListIter End() { return m_var_map.end(); };
    
  protected:

    double GetVariableValue( RootVariable* variable, int ndims=1, ... );
    long GetVariableValueInt( RootVariable* variable, int ndims=1, ... );
    TChain* GetChain() { return m_source_tree; };
  
  private:

    void ParseVariableList( std::string variablenames, std::vector< std::string >& var_names );

    std::vector< std::string > m_var_names; ///< list of variable names  
    std::map< std::string, RootVariable* > m_var_map;
    std::map< int, RootVariable* > m_var_idmap;
    TChain* m_source_tree; // For this container, all variables created will hook to this chain
  
  };

}

#endif
