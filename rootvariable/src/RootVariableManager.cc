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
#include "RootVariableManager.hh"
#include <assert.h>
#include <iostream>
#include <sstream>

#include "TChain.h"
#include "TFriendElement.h"
#include "TList.h"
#include "TBranch.h"
#include "RootVariableFormula.hh"
#include "RootVariableFunction.hh"

using namespace qosc;

RootVariableManager* RootVariableManager::m_theManager = new RootVariableManager();

RootVariableManager::RootVariableManager() {

}

RootVariableManager::~RootVariableManager() {

}

void RootVariableManager::RegisterChain( TChain* chain ) {
  
  std::set< TChain* >::iterator  it = m_chain_list.find( chain );
  // check if already registered
  if (it!=m_chain_list.end() ) return;

  m_chain_list.insert( chain );
  m_chain_name[chain] = chain->GetName();
  
  // Check if we've created a branch list for this chain
  std::map< TChain*, std::set< std::string >* >::iterator it_list = m_branch_lists.find( chain );
  if ( it_list==m_branch_lists.end() ) {
    m_branch_lists[ chain ] = new std::set< std::string >;
  }

  // Check if we've created a formula cache for this chain
  std::map< TChain*, std::set< std::string >* >::iterator it_cache = m_formula_cache.find( chain );
  if ( it_cache==m_formula_cache.end() ) {
    m_formula_cache[ chain ] = new std::set< std::string >;
  }

  // Check if we've created a variable dict for this chain
  std::map< TChain*, VariableDict* >::iterator it_vdict = m_variable_dict.find( chain );
  if ( it_vdict==m_variable_dict.end() ) {
    m_variable_dict[ chain ] = new VariableDict;
  }

}

void RootVariableManager::RegisterBranch( TChain* chain, std::string branchname ) {

  // Register the chain (and its friends)
  RegisterChain( chain );

  std::set< std::string >* branchlist = m_branch_lists[ chain ];
  if ( branchlist->find( branchname ) == branchlist->end() ) {
    branchlist->insert( branchname );
    //std::cout << "Registering " << branchname << " with chain=" << chain->GetName() << " (" << chain << ")" << std::endl;
  }
//   else {
//     std::cout << " already registered " << branchname << " with chain=" << chain->GetName() << " (" << chain << ")" << std::endl;
//   }
//   std::cin.get();
}



void RootVariableManager::DeactivateUnusedBranches() {
  
  // loop through all trees, loop 
  std::set< TChain* >::iterator it;
  for (it=m_chain_list.begin(); it!=m_chain_list.end(); it++) {
    TChain* chain = (*it);
    (*it)->SetBranchStatus("*",false);
    std::set< std::string >* branchlist = m_branch_lists[chain];
    std::set< std::string >::iterator it_b;
    for ( it_b=branchlist->begin(); it_b!=branchlist->end(); it_b++) {
      std::string branchname = *it_b;
      chain->SetBranchStatus(branchname.c_str(),true);
      //std::cout << "RootVariableManager::DeactiveUnusedBranches: Activating " << branchname << " (chain: " << chain->GetName() << " " << chain << ")" << std::endl;
    }
  }
  
}

void RootVariableManager::ResetManager() {
  // activate all branches
  std::set< TChain* >::iterator it;
  for (it=m_chain_list.begin(); it!=m_chain_list.end(); it++) {
    (*it)->SetBranchStatus("*",true);
  }  
}


RootVariableFormula* RootVariableManager::MakeRootVariableFormula( std::string formula, TChain* source_tree ) {
  // the (formula,tree) pair is used as a key to an instance of a RootVariableFormula object
  VariableDict* vardict = NULL;
  std::set< std::string >* cache = NULL;
  std::set< int >* idcache = NULL;

  SetupCaches( source_tree );

  // get variable dict and cache
  vardict = m_variable_dict[source_tree];
  cache = m_formula_cache[source_tree];
  idcache = m_id_cache[source_tree];

  // look to see if formula has been registered
  std::set< std::string >::iterator it_formula = cache->find( formula );

  if (it_formula!=cache->end()) {
    // formula found. return pointer to cached RootVariableFormula instance
    RootVariableFormula* cached_var = (RootVariableFormula*)(*vardict)[formula];
    m_var_uses[cached_var] += 1;
    return cached_var;
  }
  else {
    // no formula variable found.  create one, register it, then return it
    RootVariableFormula* var = new RootVariableFormula( formula, source_tree );
    (*vardict)[formula] = (RootVariable*)var;
    cache->insert( formula );
    idcache->insert( var->GetID() );
    m_var_uses[var] = 1;
    return var;
  }


  return NULL;

}

std::string RootVariableManager::GetCanonicalFunctionName( double (*VariableFunction)(RootVariableList&), std::string root_var_arguments, TChain* source_tree ) {
  // generate the tag we would use for this function
  std::stringstream functag;
  functag.str("");
  char funcptr[100];
  sprintf(funcptr, "%p",&(*VariableFunction));
    functag << "func_" << source_tree->GetName() << "_" << funcptr << "_" << root_var_arguments; // func_[chainname]_[function address]_[argumentlist]
    return functag.str();
}

RootVariableFunction* RootVariableManager::MakeRootVariableFunction( double(*VariableFunction)(RootVariableList&), std::string root_var_arguments, TChain* source_tree ) {

  // the (formula,tree) pair is used as a key to an instance of a RootVariableFormula object
  VariableDict* vardict = NULL;
  std::set< std::string >* cache = NULL;
  std::set< int >* idcache = NULL;

  SetupCaches( source_tree );

  // get variable dict and cache
  vardict = m_variable_dict[source_tree];
  cache = m_formula_cache[source_tree];
  idcache = m_id_cache[source_tree];

  // generate the tag we would use for this function
  std::string functag = GetCanonicalFunctionName( VariableFunction, root_var_arguments, source_tree );

  // look to see if function has been registered
  std::set< std::string >::iterator it_function = cache->find( functag );

  if (it_function!=cache->end()) {
    // function found. return pointer to cached RootVariableFunction instance
    RootVariableFunction* cached_var = (RootVariableFunction*)(*vardict)[functag];
    m_var_uses[(RootVariable*)cached_var] += 1;
    return cached_var;
  }
  else {
    // no function variable found.  create one, register it, then return it
    RootVariableFunction* var = new RootVariableFunction( functag, VariableFunction, root_var_arguments, source_tree );
    (*vardict)[functag] = (RootVariable*)var;
    cache->insert( functag );
    idcache->insert( var->GetID() );
    m_var_uses[(RootVariable*)var] = 1;
    return var;
  }
    
  return NULL;

}

RootVariableFunction* RootVariableManager::MakeRootVariableFunction( std::string varname, double(*VariableFunction)(RootVariableList&), std::string root_var_arguments, TChain* source_tree ) {

  // the (formula,tree) pair is used as a key to an instance of a RootVariableFormula object
  VariableDict* vardict = NULL;
  std::set< std::string >* cache = NULL;
  std::set< int >* idcache = NULL;

  SetupCaches( source_tree );

  // get variable dict and cache
  vardict = m_variable_dict[source_tree];
  cache = m_formula_cache[source_tree];
  idcache = m_id_cache[source_tree];

  // the user provided a name for this variable
  std::string functag = varname;

  // look to see if function has been registered
  std::set< std::string >::iterator it_function = cache->find( functag );

  if (it_function!=cache->end()) {
    // function found. return pointer to cached RootVariableFunction instance
    RootVariableFunction* cached_var = (RootVariableFunction*)(*vardict)[functag];
    m_var_uses[(RootVariable*)cached_var] += 1;
    return cached_var;
  }
  else {
    // no function variable found.  create one, register it, then return it
    RootVariableFunction* var = new RootVariableFunction( functag, VariableFunction, root_var_arguments, source_tree );
    (*vardict)[functag] = (RootVariable*)var;
    cache->insert( functag );
    idcache->insert( var->GetID() );
    m_var_uses[(RootVariable*)var] = 1;
    return var;
  }
    
  return NULL;

}

void RootVariableManager::SetupCaches( TChain* source_tree ) {
  std::map< TChain*, VariableDict* >::iterator it = m_variable_dict.find( source_tree );
  if ( it==m_variable_dict.end() ) {
    // a variable dictionary was not found for this chain.  Add one.
    m_variable_dict[source_tree] = new VariableDict;
    m_formula_cache[source_tree] = new std::set< std::string >;
    m_id_cache[source_tree] = new std::set< int >;
  }
}

void RootVariableManager::DestroyRootVariable( RootVariable* rootvar ) {
  if (rootvar==NULL) return;

  std::map< RootVariable*, int >::iterator it = m_var_uses.find( rootvar );

  if ( it==m_var_uses.end() ) return;
  else {
    int varuses = (*it).second;
    if (varuses==1) {
      /// ideally should destroy this variable.  but i am going to keep it around.
      m_var_uses[rootvar] = 0;
    }
    else {
      // decrement the variable
      m_var_uses[rootvar] -= 1;
    }
  }
}

void RootVariableManager::DumpTheVariableCache() {
  
  std::cout << "==============================================================================================================" << std::endl;
  std::cout << " RootVariableManager RootVariable Cache" << std::endl;
  std::cout << " ---------------------------------------------------------------------" << std::endl;
  std::cout << " ID   :     Name      :     Instances in use       :      Address     " << std::endl;
  std::map< TChain*, VariableDict* >::iterator it;
  for (it=m_variable_dict.begin(); it!=m_variable_dict.end(); it++) {
    TChain* chain = (*it).first;
    
    std::cout << "-- RootVariable Cache for chain '" << chain->GetName() << "' (" << chain << ") -----------------" << std::endl;
    std::set< std::string >::iterator it_var_id;
    for ( it_var_id = m_formula_cache[chain]->begin(); it_var_id!=m_formula_cache[chain]->end(); it_var_id++) {
      RootVariable* rootvar = (*m_variable_dict[chain])[(*it_var_id)];
      std::cout << " " << (*it_var_id) << " : " 
		<< "  " << rootvar->GetID() << "  :  "
		<< m_var_uses[ rootvar ] << " : " 
		<< (*m_variable_dict[chain])[(*it_var_id)] 
		<<  std::endl;
    }
  }
  std::cout << "==============================================================================================================" << std::endl;
}
