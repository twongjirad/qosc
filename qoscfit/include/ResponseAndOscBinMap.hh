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
 * --------------------------------------------------------------------------------------------
 * \class ResponseAndOscBinMap
 * \ingroup QoscFit
 * \brief Stores dictionary between oscillation bin and response bin
 *
 * This class is defined to be able to answer the question:
 * Given a true-reco energy oscillation bin, what spline truth bin should I use?
 * Assumption is that spline binning is coarser than oscillation binning -- good enough for now.
 * Reads the truth binning information from a SampleOscfitROOTPDF object.
 * -------------------------------------------------------------------------------------------*/

#ifndef __ResponseAndOscBinMap__
#define __ResponseAndOscBinMap__

#include <string>
#include <vector>
#include <map>

namespace qosc {

  class SampleOscfitROOTPDF;
  class TruthBinInfo;

  class ResponseAndOscBinMap {

  public:
    ResponseAndOscBinMap( SampleOscfitROOTPDF* sample );
    virtual ~ResponseAndOscBinMap();

    int GetSplineGlobalTruthBinFromOscGlobalTruthBin( int osc_truth_bin );
    std::string GetSplineTruthInfoNameFromOscGlobalTruthBin( int osc_truth_bin );
    //void GetOscGlobalTruthBinFromSplineGlobalTruthBin( std::vector< int >& oscbins );

    bool CheckSplineCoverage(); // If true -- all oscillation truth info bins are covered
    void PrintMap( bool stop_at_change=false );

    // below should be protected. but whatever.
    std::string samplename;

    std::map< std::string, int > spline_bininfo_id;
    int num_rec_bins;
    int num_osc_truthbins;
    int num_bininfo_objects;
    int num_largest_truthbins;
    int* osc_to_spline_map; // [osc_truth_bin][2]  : first entry is bin info object index, second is bin number
    TruthBinInfo** spline_bin_infos; // [num_bininfo_objects]
    TruthBinInfo* osc_bin_info;

  protected:
    void BuildMap( SampleOscfitROOTPDF* sample );

  };

}

#endif
