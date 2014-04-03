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
/** --------------------------------------------------------------------------------------------
 * \class WeightOscEvent
 * \ingroup QoscFit
 * \brief Abstract base class for vent-by-event oscillations
 * 
 * ---------------------------------------------------------------------------------------------- */

#ifndef __WeightOscEvent__
#define __WeightOscEvent__

#include <string>

#include "Weight.hh" // base class which has function definitions needed to fill a HistRootVariable from ROOT tree
//#include "RootVariableList.hh" 
//#include "SinThetaForms.hh"

class TChain;
//class BargerPropagator;

namespace qosc {

  class WeightOscEvent : public Weight {

  public:

    WeightOscEvent( TChain* source_chain );
    virtual ~WeightOscEvent();
    virtual double CalculateWeight() = 0; ///< return oscillation weight
    virtual double CalculateWeight( double var ) = 0; ///< return oscillation weight given some variable
    void SetMapMode( bool mapmode ) { fMapMode = mapmode; }; ///< in map mode, all oscillation weights return 1.0. For building templates.

    bool HaveParsBeenSet() { return fParamsSet; };

  protected:

    // Flags
    bool fParamsSet; ///< indicate that parameters have been set
    void SetParamsSetFlag( bool state ) { fParamsSet = state; };

    bool fMapMode; ///< run oscillation wighting in map-mode

    TChain* m_source_chain;
    TChain* GetSourceChain() { return m_source_chain; };
    void SetSourceChain( TChain* source ) { m_source_chain = source; };
  };

}

#endif
