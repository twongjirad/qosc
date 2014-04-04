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
/* --------------------------------------------------------------------------------------------
 * \class SampleSetupPars
 * \ingroup QoscFit
 * \brief Class containing information needed to define a set of a Sample instance.
 *
 * The class SampleSetupParser is responsible for extract the information stored in this class
 * -------------------------------------------------------------------------------------------*/

#ifndef __SampleSetupPars__
#define __SampleSetupPars__

#include <string>
#include <set>
#include <vector>

class TChain;

namespace qosc {

  class SampleBinInfo;

  class SampleSetupPars {
  public:

    SampleSetupPars( std::string sample_name );
    virtual ~SampleSetupPars();
    void Print();


    std::string name;
    std::string formula;
    std::string formulaY;
    std::string selection;
    std::string weight;
  
    std::string chain_name;
    TChain* chain;
  
    int ndims;
    int nbinsX;
    double* binedgesX;
    int nbinsY;
    double* binedgesY;
    SampleBinInfo* bininfo; /// This is the observble binning used for fits
  
    // These containers are used to store user-defined histogrammed quantities to associate with each observable bin.
    std::set< std::string > m_bininfo_names;
    std::vector< SampleBinInfo* > m_bininfo_list;
    std::string bininfo_name_for_osc; // special bin info object for oscillations. Assumes 1D histogram of true neutrino energy in GeV
    std::set< std::string > bininfo_names_for_splines; // names of the bin info objects used for spline response functions. Assumes 1D histogram of true neutrino energy in GeV

  };
}

#endif
