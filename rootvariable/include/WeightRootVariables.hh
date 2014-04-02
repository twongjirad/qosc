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
 * -----------------------------------------------------------------------
 * \class WeightRootVariables
 * \ingroup rootvariable
 * \brief This class calculates the weight for PDFRootVariable objects.
 * 
 * -----------------------------------------------------------------------
 */

#ifndef __WeightRootVariables__
#define __WeightRootVariables__

#include "Weight.hh"
#include <string>
#include "RootVariableList.hh"
class TChain;

namespace qosc {

  class WeightRootVariables : public Weight {

  public:

    WeightRootVariables( TChain* datachain, std::string rootvars );
    virtual ~WeightRootVariables();

    virtual double CalculateWeight();
    virtual double CalculateWeight( double variable ) { return CalculateWeight(); };

  private:
  
    TChain* m_source_chain;
    RootVariableList m_weight_vars;

  };

}

#endif
