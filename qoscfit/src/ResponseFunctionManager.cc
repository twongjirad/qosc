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

#include "ResponseFunctionManager.hh"

#include <assert.h>
#include <iomanip>

#include "TSpline.h"

#include "SampleManager.hh"
#include "ParameterManager.hh"
#include "TruthBinInfo.hh"
#include "Sample.hh"
#include "SampleOscfitROOTPDF.hh"
#include "ResponseFunction.hh"
#include "ResponseFunctionSpline.hh"
#include "ResponseFunctionPoly.hh"
#include "ResponseAndOscBinMap.hh"
//#include "SampleSetupParser.hh"


//#define DEBUG_ResponseFunctionManager

using namespace qosc;

ResponseFunctionManager::ResponseFunctionManager() {
  m_sampleman = NULL;
  m_numsamples  = 0;
  m_sample_func_list = NULL;
}

ResponseFunctionManager::ResponseFunctionManager( SampleManager* sampleman ) {
  Initialize( sampleman );
}

ResponseFunctionManager::~ResponseFunctionManager() {

  // Clear stored info
  if ( m_sampleman ) {
    for ( SampleManager::SampleListIter it=m_sampleman->SampleListBegin(); it!=m_sampleman->SampleListEnd(); it++ ) {
      Sample* asample = m_sampleman->GetSample( *it );
      int sampleid = asample->GetID();
      DestroyResponseFunctions( sampleid );
    }
  }
  
  delete [] m_sample_func_list;

}

void ResponseFunctionManager::Initialize( SampleManager* sampleman ) {
  // we load all the details of the sample's truth binning for spline functions
  std::cout << __PRETTY_FUNCTION__ << "[" << this << "]" << std::endl;

  m_sampleman = sampleman;
  m_numsamples = sampleman->GetNumberOfSamples();

  m_sample_func_list = new sample_response_ptrs[ m_numsamples  ];

  for ( SampleManager::SampleListIter it=sampleman->SampleListBegin(); it!=sampleman->SampleListEnd(); it++ ) {
    SampleOscfitROOTPDF* asample = dynamic_cast<SampleOscfitROOTPDF*>( sampleman->GetSample( *it ) );
    if ( !asample ) assert(false);
    Hist* pdf = asample->GetSampleHistogram();
    std::cout << "  Load Response Array for " << asample->GetName() << std::endl;
    m_sample_lookup[ asample->GetName() ] = asample->GetID();
    m_sample_func_list[ asample->GetID() ].sampleid = asample->GetID();
    m_sample_func_list[ asample->GetID() ].samplename = *it;
    m_sample_func_list[ asample->GetID() ].n_rec_bins = asample->GetNumberOfBins();
    m_sample_func_list[ asample->GetID() ].n_truth_bins = 0;
    m_sample_func_list[ asample->GetID() ].n_bin_infos = asample->GetNumOfResponseBinInfoNames();
    m_sample_func_list[ asample->GetID() ].resp_ptrs = NULL;
    TruthBinInfo* truthinfo = dynamic_cast< TruthBinInfo* >( pdf->GetMasterListBinInfo( asample->GetOscBinInfoName() ) ); // this makes this item special
    if ( !truthinfo ) {
      std::cout << "Did not find 'truthinfo' instance with name=" << m_truth_bin_info_name << " in sample=" << *it << std::endl;
      assert(false);
    }
    // Set up response array
    int ntruthbins = truthinfo->GetTotalTruthBins();
    m_sample_func_list[ asample->GetID() ].oscInfoTemplate = truthinfo;
    m_sample_func_list[ asample->GetID() ].n_truth_bins = ntruthbins;
    m_sample_func_list[ asample->GetID() ].oscToSplineMap = new ResponseAndOscBinMap( asample );

#ifdef DEBUG_ResponseFunctionManager
    // This is kind of a waste of memory.  This Map should only be created once per sample, but instead, we're createing one per parameter...
    // Option is to put it in Sample. The other is to use a singleton that stores all such maps.
    m_sample_func_list[ asample->GetID() ].oscToSplineMap->PrintMap();
    std::cin.get();
#endif

    std::vector< std::string > spline_info_objs;
    m_sample_func_list[ asample->GetID() ].binInfoTemplates = new TruthBinInfo*[ asample->GetNumOfResponseBinInfoNames() ];
    asample->GetResponseBinInfoNames( spline_info_objs );
    m_sample_func_list[ asample->GetID() ].n_max_truthinfo = 0;
    for ( std::vector< std::string >::iterator it=spline_info_objs.begin(); it!=spline_info_objs.end(); it++ ) {
      TruthBinInfo* info = dynamic_cast<TruthBinInfo*>( asample->GetMasterListBinInfoInstance( *it ) );
      if ( info ) {	
	m_sample_func_list[ asample->GetID() ].binInfoTemplates[ asample->GetResponseBinInfoIndex( *it ) ] = info;
	m_sample_func_list[ asample->GetID() ].spline_info_lookup[ *it ] = asample->GetResponseBinInfoIndex( *it );
      }
      else
	assert(false);
    }
    m_sample_func_list[ asample->GetID() ].n_max_truthinfo =  m_sample_func_list[ asample->GetID() ].oscToSplineMap->num_largest_truthbins;
    int total_ptrs = asample->GetNumberOfBins()*asample->GetNumOfResponseBinInfoNames()*m_sample_func_list[ asample->GetID() ].oscToSplineMap->num_largest_truthbins;
    m_sample_func_list[ asample->GetID() ].n_total_funcs = total_ptrs;
    m_sample_func_list[ asample->GetID() ].resp_ptrs = new ResponseFunction*[ total_ptrs ]; // create an array of pointers
    std::cout << "Allocating " << total_ptrs << " pointer spots" << std::endl;
    memset( m_sample_func_list[ asample->GetID() ].resp_ptrs, 0, sizeof( ResponseFunction* )*total_ptrs );
  }
}

void ResponseFunctionManager::DestroyResponseFunctions( int sampleid ) {
  if ( !m_sample_func_list ) return;
  if ( m_sample_func_list[ sampleid ].resp_ptrs ) {
    // Clean out response parameters
    for (int irec=0; irec<m_sample_func_list[sampleid].n_rec_bins; irec++) {
      for (int ibinfo=0; ibinfo<m_sample_func_list[sampleid].n_bin_infos; ibinfo++) {
	for (int itru=0; itru<m_sample_func_list[sampleid].n_max_truthinfo; itru++) {
	  ResponseFunction* response = GetResponseFunction( sampleid, irec, ibinfo, itru );
	  if ( response )
	    delete response;
	}
      }
    }
    delete [] m_sample_func_list[ sampleid ].resp_ptrs;
    delete [] m_sample_func_list[ sampleid ].binInfoTemplates;
    delete m_sample_func_list[ sampleid ].oscToSplineMap;
    m_sample_func_list[ sampleid ].resp_ptrs = NULL;
  }
}

bool ResponseFunctionManager::IsSampleDefinedInManager( std::string samplename ) {
  return m_sampleman->IsSampleDefined( samplename );
  if ( m_sample_lookup.find( samplename )!=m_sample_lookup.end() )
    return true;
  return false;
}

ResponseFunction* ResponseFunctionManager::GetResponseFunction( int sampleid, int recon_bin_index, int bin_info_index, int spline_truth_bin_index ) {
  int nbinfo = m_sample_func_list[sampleid].n_bin_infos;
  int ntruth = m_sample_func_list[sampleid].n_max_truthinfo;
  return *(m_sample_func_list[sampleid].resp_ptrs + recon_bin_index*(nbinfo*ntruth) + bin_info_index*(ntruth) + spline_truth_bin_index);  
}

ResponseFunction* ResponseFunctionManager::GetResponseFunction( int sampleid, int recon_bin_index, int truth_bin_index ) {

  int nbinfo = m_sample_func_list[sampleid].n_bin_infos;
  int ntruth = m_sample_func_list[sampleid].n_max_truthinfo;
  int spine_truth_bin = m_sample_func_list[sampleid].oscToSplineMap->GetSplineGlobalTruthBinFromOscGlobalTruthBin( truth_bin_index );
  std::string bin_info_name = m_sample_func_list[sampleid].oscToSplineMap->GetSplineTruthInfoNameFromOscGlobalTruthBin( truth_bin_index );
  int binfo_index =  m_sample_func_list[sampleid].spline_info_lookup[ bin_info_name ];
  //std::cout << "returning response for " << bin_info_name << " spline bin " << binfo_index << std::endl;
  return GetResponseFunction( sampleid, recon_bin_index, binfo_index, spine_truth_bin );

}

ResponseFunction* ResponseFunctionManager::GetResponseFunction( std::string samplename, int recon_bin_index, int truth_bin_index ) {
  int sampleid = m_sample_lookup[samplename];
  return GetResponseFunction( sampleid, recon_bin_index, truth_bin_index );
}

void ResponseFunctionManager::LoadResponseFunction( int sampleid, int recon_bin_index, std::string bin_info_name, int spline_truth_bin_index, ResponseFunction* response ) {

//   if ( recon_bin_index<0 || recon_bin_index>=nrec ) {
//     std::cout << "ResponseFunctionManager::LoadResponseFunction: Reconstructed bin out of range: " << recon_bin_index << std::endl;
//     assert(false);
//   }
//   if ( truth_bin_index<0 || truth_bin_index>=ntruth ) {
//     std::cout << "ResponseFunctionManager::LoadResponseFunction: Truth bin out of range: " << recon_bin_index << std::endl;
//     assert(false);
//   }

  int nbinfo = m_sample_func_list[sampleid].n_bin_infos;
  int ntruth = m_sample_func_list[sampleid].n_max_truthinfo;
  int binfo_index =  m_sample_func_list[sampleid].spline_info_lookup[ bin_info_name ];  
#ifdef DEBUG_ResponseFunctionManager
  std::cout << "Load response function: "
	    << " sampleid=" << sampleid 
	    << " recbin=" << recon_bin_index
	    << " bin-info=" << bin_info_name << " of " << nbinfo
	    << " truthbin=" << spline_truth_bin_index << " of " << ntruth
	    << std::endl;
  std::cout << "put " << response << " (Address " << &response << ") " 
	    << " into addres=" << &*(m_sample_func_list[sampleid].resp_ptrs + recon_bin_index*(nbinfo*ntruth) + binfo_index*(ntruth) + spline_truth_bin_index)
	    << " (current contents=" << *(m_sample_func_list[sampleid].resp_ptrs + recon_bin_index*(nbinfo*ntruth) + binfo_index*(ntruth) + spline_truth_bin_index) << ")"
	    << std::endl;
#endif
  *(m_sample_func_list[sampleid].resp_ptrs + recon_bin_index*(nbinfo*ntruth) + binfo_index*(ntruth) + spline_truth_bin_index) = response;
}

void ResponseFunctionManager::LoadResponseFunction( std::string samplename, int recon_bin_index, std::string bin_info_name, int spline_truth_bin_index, ResponseFunction* response ) {
  int sampleid = m_sample_lookup[samplename];
  LoadResponseFunction( sampleid, recon_bin_index, bin_info_name, spline_truth_bin_index, response );
}

int ResponseFunctionManager::GetNumberOfTruthBins( std::string samplename ) {
  if ( m_sample_lookup.find( samplename )==m_sample_lookup.end() )
    assert(false);

  int sample_index = m_sample_lookup[samplename];
  return m_sample_func_list[sample_index].n_truth_bins;
}

int ResponseFunctionManager::GetNumberOfReconBins( std::string samplename ) {
  if ( m_sample_lookup.find( samplename )==m_sample_lookup.end() )
    assert(false);

  int sample_index = m_sample_lookup[samplename];
  return m_sample_func_list[sample_index].n_rec_bins;
}

void ResponseFunctionManager::PrintShortReport() {
  size_t total_mem = 0;
  for (int isample=0; isample<m_numsamples; isample++) {
    std::cout << "Sample [" << isample << "] " << m_sample_func_list[isample].samplename << std::endl;
    sample_response_ptrs* resp_ptrs = m_sample_func_list + isample;

    int nloaded = 0;
    int nempty = 0;

    for (int irec=0; irec<resp_ptrs->n_rec_bins; irec++) {
      for (int itruth=0; itruth<resp_ptrs->n_truth_bins; itruth++) {
	ResponseFunction* respfunc = GetResponseFunction( isample, irec, itruth );
	if ( respfunc==NULL ) nempty++;
	else { 
	  total_mem += sizeof(*respfunc);
	  nloaded++;
	}
      }
    }
    
    std::cout << "  Num Pointers: " << resp_ptrs->n_total_funcs << " ( size of ptr array = " <<  resp_ptrs->n_total_funcs*sizeof(ResponseFunction*)*1.0e-6 << " MB )" << std::endl;
    std::cout << "  Num Empty Pointers: " << nempty << "  size=" << nempty*sizeof(ResponseFunction*)*1.0e-6 << " MB" << std::endl;
    std::cout << "  Num Loaded Pointers:  " << nloaded  << std::endl;
    std::cout << "  Total size of objects: " << total_mem << std::endl;
    
  }
}

int ResponseFunctionManager::GetSampleID( std::string samplename ) {
  if ( IsSampleDefinedInManager( samplename ) )
    return m_sample_lookup[samplename];
  assert(false);
}

void ResponseFunctionManager::Optimize( double spline_threshold ) {
  assert(false); // I've broken this
  /// This routine performs some spline optimizations
  /// (1) We go through all the truth bins and remove spline objects where the expected number of event is zero.
  /// (2) We also remove splines which have no effect below the spline_threshold
  /// This is teo help reduce the memory allocated by the analyses which is reaching several to tens of GBs depending on the number of bins and parameters
  int nresponses = 0;
  int ndeleted = 0;
  for ( SampleManager::SampleDictIter it=m_sampleman->SampleDictBegin(); it!=m_sampleman->SampleDictEnd(); it++ ) {
    Sample* asample = (*it).second;
    for (int reconbin=0; reconbin<asample->GetNumberOfHistogramBins(); reconbin++) {
      TruthBinInfo* truth_bins = dynamic_cast<TruthBinInfo*>( asample->GetSampleHistogram()->GetBinInfo( m_truth_bin_info_name, reconbin ) );
      int ntruebins = 0;
      for (int j=0; j<truth_bins->GetNumberOfTruthHists(); j++) {
	TH1D* intruhist = (TH1D*)truth_bins->GetTruthHist( j );
	for (int itrue=0; itrue<intruhist->GetNbinsX(); itrue++) {
	  double ntrueevents = intruhist->GetBinContent( itrue+1 );
	  ResponseFunction* respfunc = GetResponseFunction( (*it).first, reconbin, ntruebins );
	  nresponses++;
	  if ( ntrueevents==0 ) {
	    if ( respfunc ) {
	      delete respfunc;
	      LoadResponseFunction( (*it).first, reconbin, ntruebins, NULL );
	      ///int ntruth = m_sample_func_list[sampleid].n_truth_bins;
	      ///return *(m_sample_func_list[sampleid].resp_ptrs + recon_bin_index*ntruth + truth_bin_index);
	      ndeleted++;
	    }
	  }
	  ntruebins++;
	}// truth bins
      }// truth histograms
    }//end of recon bin
  }//end of sample loop
  std::cout << __PRETTY_FUNCTION__ << ": nresponses=" << nresponses << " ndeleted=" << ndeleted << " (" << float(ndeleted)/float(nresponses)*100.0 << "%)" << std::endl;
}

void ResponseFunctionManager::GetListOfSplineBinInfoNames( std::string sample, std::set< std::string >& bininfo ) {
  int sampleid = m_sample_lookup[sample];
  for ( std::map< std::string, int >::iterator it=m_sample_func_list[ sampleid ].spline_info_lookup.begin();
	it!=m_sample_func_list[ sampleid ].spline_info_lookup.end(); it++ ) {
    bininfo.insert( (*it).first );
  }
}
