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

#include "WeightOscTruthBinInfo.hh"
#include <assert.h>

#include "TruthBinInfo.hh"

namespace qosc {

  WeightOscTruthBinInfo::WeightOscTruthBinInfo( TruthBinInfo* truth_bin_def ) { 
    /** 
	This takes a source chain file (and its friend trees) along with the leaf names 
	for the neutrino energy in GeV (nuE_GeV_var), the neutrino flux flavor (nuflux_var)
	and the neutrino cross section flavor (nuxsec_var).
	It then utilizes a RootVariableList to create a hook into the ROOT chain in order 
	to access the value of these variables event by event.
    */

    fParamsSet = false; // Parameters have not been set yet
    fMapMode = false;
    SetVerbose(0);

    Configure( truth_bin_def );

  }


  WeightOscTruthBinInfo::~WeightOscTruthBinInfo() {}

}
