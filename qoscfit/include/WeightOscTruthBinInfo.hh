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
 * \class WeightOscTruthBinInfo
 * \ingroup QoscFit
 * \brief Abtract base class that implements oscillate effect via TruthBinInfo class
 *
 * ------------------------------------------------------------------------------------------- */

#ifndef __WeightOscTruthBinInfo__
#define __WeightOscTruthBinInfo__

#include <string>

class TH1D;
class TAxis;

namespace qosc {

  class TruthBinInfo;

  class WeightOscTruthBinInfo {

  public:

    WeightOscTruthBinInfo( TruthBinInfo* truth_bin_def );
    virtual ~WeightOscTruthBinInfo();

    virtual double CalcEventsOsc( TruthBinInfo* truth_info, bool modifyhist=false ) = 0; ///< Takes an instance of TruthBinInfo and uses it to calculate the number of events after oscillation
    void SetMapMode( bool mapmode ) { fMapMode = mapmode; }; ///< oscillation effect is 1.0. Used to build templates. [needed?]
    bool GetMapMode() { return fMapMode; };
    bool AreParamsSet() { return fParamsSet; };
    void Configure( TruthBinInfo* truth_bin_def ) { __configure__(truth_bin_def); };

  protected:
    // Flags
    bool fParamsSet; /// flag that turns on when osc parameters set for the first time
    bool fMapMode; /// in map mode, oscillation weight is always one

    virtual void __configure__( TruthBinInfo* truth_bin_def ) = 0;
    
  protected:
    // Verbosity flag and its get/set function
    int fVerbose;
  public:
    void SetVerbose( int verbose ) { fVerbose = verbose; };
    int GetVerbose() { return fVerbose; };

  };

}

#endif
