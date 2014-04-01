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
 * \class RootVariableManager
 * \ingroup rootvariable
 *
 * \brief This class is a singleton and creates and stores all the copies of the RootVariables created
 *
 * This class exists for optimization reasons.  We don't want multiple copies of a variable
 * that accesses the root tree. We only want to have one per formula.
 * 
 *
 * -----------------------------------------------------------------------------------------------
 */

#ifndef __RootVariableManager__
#define __RootVariableManager__

#include <set>
#include <map>
#include <string>
#include <sstream>

class TChain;

namespace qosc {

  class RootVariable;
  class RootVariableFormula;
  class RootVariableFunction;
  class RootVariableList;

  class RootVariableManager {
  
  private:
  
    RootVariableManager();
    virtual ~RootVariableManager();
  
    static RootVariableManager* m_theManager;
    std::set< TChain* > m_chain_list;
    std::map< TChain*, std::string > m_chain_name;
    std::map< TChain*, std::set< std::string >* > m_branch_lists;
  
    // cache formulas, so I don't repeat evaluations
    typedef std::map< std::string, RootVariable* > VariableDict; ///< dict of RootVariable instances we have created
    std::map< TChain*, std::set< std::string >* > m_formula_cache; ///< a list of the names of variables in the cache, separated into what chain they are associated to
    std::map< TChain*, std::set< int >* > m_id_cache; ///< a map of variable ID in the cache, separated into what chain they are associated to
    std::map< TChain*, VariableDict* > m_variable_dict;
    std::map< RootVariable*, int > m_var_uses; ///< counts the number of uses of this variable instance

  public:
  
    static RootVariableManager* GetTheRootVarManager() { return m_theManager; };
    void RegisterChain( TChain* chain );
    void RegisterBranch( TChain* chain, std::string branchname );
    void DeactivateUnusedBranches();
    void ResetManager();
    RootVariableFormula* MakeRootVariableFormula( std::string formula, TChain* source_tree ); ///< Ask the RootVariableManager() for a RootVariableFormula.  Will give you a repeat if one already exists.
    RootVariableFunction* MakeRootVariableFunction( double (*VariableFunction)(RootVariableList&), std::string root_var_arguments, TChain* source_tree ); ///< Ask the manager() for a RootVariableFunction.
    RootVariableFunction* MakeRootVariableFunction( std::string name, double (*VariableFunction)(RootVariableList&), std::string root_var_arguments, TChain* source_tree ); ///< Ask the manager for a RootVariableFunction.
    void DestroyRootVariable( RootVariable* rootvar );
    void DumpTheVariableCache();
    std::string GetCanonicalFunctionName( double (*VariableFunction)(RootVariableList&), std::string root_var_arguments, TChain* source_tree );

    // cache formulas, so I don't repeat evaluations
    //RootVariable* SearchForFormula( TChain* chain, std::string formula );

  protected:
    void SetupCaches( TChain* source_tree );
  
  private:

    //int m_entryNumber;

  };

}
#endif 
