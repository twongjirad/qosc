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

#include "ResponseAndOscBinMap.hh"

#include <assert.h>
#include <iomanip>

#include "SampleOscfitROOTPDF.hh"
#include "TruthBinInfo.hh"
#include "SampleSetupPars.hh"

using namespace qosc;

ResponseAndOscBinMap::ResponseAndOscBinMap( SampleOscfitROOTPDF* sample ) {
  osc_bin_info = NULL;
  BuildMap( sample );  
}

ResponseAndOscBinMap::~ResponseAndOscBinMap() {
  delete [] osc_to_spline_map;
}

void ResponseAndOscBinMap::BuildMap( SampleOscfitROOTPDF* sample ) {
  
  // Gather all the information we need
  osc_bin_info = dynamic_cast<TruthBinInfo*>( sample->GetMasterListBinInfoInstance( sample->GetOscBinInfoName() ) );
  if ( !osc_bin_info ) 
    assert(false);
  std::vector< std::string > bininfo_names_for_splines;
  sample->GetResponseBinInfoNames( bininfo_names_for_splines );
  num_bininfo_objects = bininfo_names_for_splines.size();
  spline_bin_infos = new TruthBinInfo*[num_bininfo_objects];

  int ispline_info = 0;
  int largest_truths = 0;
  for ( std::vector< std::string >::iterator it=bininfo_names_for_splines.begin(); it!=bininfo_names_for_splines.end(); it++ ) {
    spline_bininfo_id[ *it ] = ispline_info;
    spline_bin_infos[ispline_info] = dynamic_cast<TruthBinInfo*>( sample->GetMasterListBinInfoInstance( *it ));
    if ( !spline_bin_infos[ispline_info] ){
      std::cout << "Did not find spline_info=" << *it << " in the Sample Master List of User Bin Info." << std::endl;
      assert(false);
    }
    if ( spline_bin_infos[ispline_info]->GetTotalTruthBins()>largest_truths )
      largest_truths = spline_bin_infos[ispline_info]->GetTotalTruthBins();
    ispline_info++;
  }

  num_rec_bins = sample->GetNumberOfHistogramBins();
  num_osc_truthbins = osc_bin_info->GetTotalTruthBins();
  num_largest_truthbins = largest_truths;
  
  osc_to_spline_map = new int[ num_osc_truthbins*2 ];
  memset( osc_to_spline_map, 0, sizeof(int)*num_osc_truthbins*2 );
  
  // OK Now we do the mapping
  for (int ioscbin=0; ioscbin<num_osc_truthbins; ioscbin++) {
    // Set default values
    *( osc_to_spline_map + ioscbin*(2) )    = -1;
    *( osc_to_spline_map + ioscbin*(2) + 1) = -1;
    // Get the cut name and definition for this bin
    std::string cutname = osc_bin_info->GetCutNameFromTruthBin( ioscbin );
    std::string cutdef  = osc_bin_info->GetCutDefinitionFromTruthBin( ioscbin );
    // Now we look at the spline bins and see if we can find a match in the cut
    for (int ispline_info=0; ispline_info<num_bininfo_objects; ispline_info++) {
      TruthBinInfo* bin_info = spline_bin_infos[ispline_info]; // spline bin info
      if ( !bin_info )
	assert(false);
      //bin_info->Print();
      //std::cin.get();
      for ( int icut=0; icut<bin_info->m_ncutlists; icut++ ) {
	// loop through the different cut definitions in this setup.
	//if ( bin_info->m_cut_list[icut].formula==cutdef ) {
	if ( bin_info->m_cut_list[icut].name==cutname ) {
	  double bin_energy = osc_bin_info->GetTruthHist( icut )->GetBinCenter( ioscbin%osc_bin_info->GetNumberOfTruthBinsPerHist()+1 );
	  int spline_bin = bin_info->m_cut_list[icut].hist->FindBin( bin_energy ) - 1;
	  *( osc_to_spline_map + ioscbin*(2) )     = ispline_info;
	  *( osc_to_spline_map + ioscbin*(2) + 1 ) = bin_info->m_cut_list[icut].start_bin + spline_bin;
// 	  if ( cutname=="numu_ccpi" ) {
// 	    std::cout << ioscbin << " " << bin_energy << " " << spline_bin << " " << std::endl;
// 	    std::cin.get();
// 	  }
	  break;
	}
      }
    }
  }


}


bool ResponseAndOscBinMap::CheckSplineCoverage() {
  for (int ioscbin=0; ioscbin<num_osc_truthbins; ioscbin++) {
    if ( *(osc_to_spline_map + 2*ioscbin )==-1 || *(osc_to_spline_map+2*ioscbin+1)==-1 )
      return false;
  }
  return true;
}

int ResponseAndOscBinMap::GetSplineGlobalTruthBinFromOscGlobalTruthBin( int osc_truth_bin ) {
  return *(osc_to_spline_map + 2*osc_truth_bin + 1 );
}

std::string ResponseAndOscBinMap::GetSplineTruthInfoNameFromOscGlobalTruthBin( int osc_truth_bin ) {
  return spline_bin_infos[  *(osc_to_spline_map + 2*osc_truth_bin + 0 ) ]->m_truth_name;
}

void ResponseAndOscBinMap::PrintMap( bool stop_at_change ) {
  std::cout << std::setw(30) << "OSCNAME" << std::setw(10) << "OSC BIN" << std::setw(40) << "SPL.NAME" << std::setw(20) << "SPL CUT" << std::setw(10) << "SPL. BIN" << std::endl;
  std::string lastcut = "";
  std::string lastcut_spl = "";
  for (int ioscbin=0; ioscbin<num_osc_truthbins; ioscbin++) {
    int isplinebin = GetSplineGlobalTruthBinFromOscGlobalTruthBin( ioscbin );
    int isplinebinfo = *(osc_to_spline_map + 2*ioscbin + 0 );    
    if ( stop_at_change || (lastcut!= osc_bin_info->GetCutNameFromTruthBin( ioscbin ) || lastcut_spl!=spline_bin_infos[ isplinebinfo ]->GetCutNameFromTruthBin( isplinebin ) ) ) {
      std::cout << std::setw(30) <<  osc_bin_info->GetCutNameFromTruthBin( ioscbin )
		<< std::setw(10) << ioscbin
		<< std::setw(40) << GetSplineTruthInfoNameFromOscGlobalTruthBin( ioscbin )
		<< std::setw(20) << spline_bin_infos[ isplinebinfo ]->GetCutNameFromTruthBin( GetSplineGlobalTruthBinFromOscGlobalTruthBin( ioscbin ) )
		<< std::setw(10) << GetSplineGlobalTruthBinFromOscGlobalTruthBin( ioscbin )
		<< std::endl;
    }
    if ( stop_at_change && (lastcut!= osc_bin_info->GetCutNameFromTruthBin( ioscbin ) || lastcut_spl!=spline_bin_infos[ isplinebinfo ]->GetCutNameFromTruthBin( isplinebin ) ) )
      std::cin.get();
    lastcut =  osc_bin_info->GetCutNameFromTruthBin( ioscbin );
    lastcut_spl = spline_bin_infos[ isplinebinfo ]->GetCutNameFromTruthBin( isplinebin );
  }
}

