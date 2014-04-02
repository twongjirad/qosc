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

#include "WeightRootVariables.hh"
#include <iostream>
#include <assert.h>
#include "TChain.h"
#include "RootVariable.hh"

using namespace qosc;

WeightRootVariables::WeightRootVariables( TChain* datachain, std::string rootvars ) {

  m_weight_vars.SetChain( datachain );
  m_weight_vars.Add( rootvars );

  // check if multiplicity of all weights is 1: no arrays!
  std::vector< std::string > weight_names;
  m_weight_vars.GetNameList( weight_names );
  for (int i=0; i<weight_names.size(); i++) {
    RootVariable* weight_var = m_weight_vars.GetVariable( weight_names.at(i) );
    if ( weight_var->IsVariableAnArray() ) {
      std::cout << "The weight variable, " << weight_names.at(i) << ", is an array.  This is not allowed." << std::endl;
      assert(false);
    }
  }

}

WeightRootVariables::~WeightRootVariables() {

}


double WeightRootVariables::CalculateWeight() {

  double weight = 1.0;
  std::vector< std::string > weightnames;
  m_weight_vars.GetNameList( weightnames );
  for (int i=0; i<m_weight_vars.GetEntries(); i++) {
    std::string name = weightnames.at(i);
    weight *= m_weight_vars.GetVariableValue( name );
  }

  return weight;
}
