#include "SampleManager.hh"
#include <assert.h>
#include <iostream>
#include <sstream>
#include "SeparateString.hh"
#include "ParameterManager.hh"

using namespace qosc;

SampleManager::SampleManager() {
  fSampleOrdered = false;
  fTotalNumberOfBins = 0;
  SetVerbose(0);
}

SampleManager::~SampleManager() {
}

// --------------------------------------------------------------------------------------------------------------------------------
// ------------------
// Sample Management
// ------------------

void SampleManager::RegisterSample( std::string name, Sample* asample ) {
  Sample* test = GetSample( name );
  int sampleid = asample->GetID();
  if ( test ) {
    std::cout << "SampleManager::RegisterSample -- WARNING! Re-registering sample with name '" << name << "'" << std::endl;
  }
  if ( !asample ) {
    std::cout << "SampleManager::RegisterSample -- ERROR! Sample instance given in NULL" << std::endl;
    assert(false);
  }
  m_sample_dict[name] = asample;
  if ( sampleid<0) {
    sampleid = m_sample_dict.size()-1;
    asample->SetID( sampleid );
  }
  m_sample_id[ name ] = sampleid;
  m_indexed_samples[ sampleid ] = name;
  
}

Sample* SampleManager::GetSample( std::string name ) {
  SampleDictIter it = m_sample_dict.find( name );
  if ( it!=SampleDictEnd() )
    return (*it).second;
  else 
    return NULL;
}

// --------------------------------------------------------------------------------------------------------------------------------
// ----------------------
// Sample Bin Management
// ----------------------

int SampleManager::GetNumberOfBins() {
  // The following is calculated in Initialize
  return fTotalNumberOfBins;
}

int SampleManager::GetNumberOfSampleBins( std::string samplename ) {
  return GetSample( samplename )->GetNumberOfBins();
}

void SampleManager::SetSampleOrder( std::string samplelist ) {
  m_sample_order.clear();
  fTotalNumberOfBins = 0;
  m_bin_index_dict.clear();

  int binindex = 0;
  std::vector< std::string > samples;
  SeparateString( samplelist, samples );
  std::vector< std::string >::iterator it;
  for (it=samples.begin(); it!=samples.end(); it++) {
    if ( !IsSampleDefined( (*it) ) ) {
      std::cout << "SampleManager::SetSampleOrder -- ERROR! Sample '" << (*it) << "' not defined." << std::endl;
      assert(false);
    }
    Sample* asample = GetSample( *it );
    if ( asample->IsActive() ) {
      asample->IndexBins();
      m_sample_order.push_back( *it );
      fTotalNumberOfBins += asample->GetNumberOfBins();
      for (int n=0; n<asample->GetNumberOfBins(); n++) {
	m_bin_index_dict[binindex] = *it;
	asample->PairGlobalBinToLocalBin( binindex, n );
	binindex++;
      }
    }
    else {
      std::cout << "SampleManager::SetSampleOrder -- WARNING! Sample '" << (*it) << "' not active and will be skipped." << std::endl;
    }
  }
  fSampleOrdered = true;
}

void SampleManager::SetupBins() {
  if ( !fSampleOrdered ) {
    // auto-order set. also sets up the m_bin_index_dict
    //std::cout << __PRETTY_FUNCTION__ << std::endl;
    std::string samplelist = "";
    for (SampleDictIter it=m_sample_dict.begin(); it!=m_sample_dict.end(); it++) {
      Sample* asample = (*it).second;
      if ( asample->IsActive() ) {
	samplelist += (*it).first + ";";
      }
    }
    SetSampleOrder( samplelist );
  }
}

double SampleManager::GetExpectedBinContent( int bin ) {
  return GetSampleFromBinIndex( bin )->GetExpectedContentByGlobalIndex( bin );
}

double SampleManager::GetDataBinContent( int bin ) {
  return GetSampleFromBinIndex( bin )->GetDataContentByGlobalIndex( bin );
}

void SampleManager::SetDataBinContent( int bin, double data ) {
  GetSampleFromBinIndex( bin )->SetDataContentByGlobalIndex( bin, data );
}

// --------------------------------------------------------------------------------------------------------------------------------
// Fill Bins
void SampleManager::FillSampleBins( ParameterManager* systerms, bool fillnominal ) {
  if ( GetVerbose()>=2 ) {
    std::cout << "--------------------------------" << std::endl;
    std::cout << "SampleManager::FillSampleBins" << std::endl;
  }
  // Default Fill Bins
  for ( SampleDictIter it=SampleDictBegin();it!=SampleDictEnd(); it++ ) {
    // Fill Sample Bins
    (*it).second->SetVerbose( GetVerbose() );
    (*it).second->FillBins( systerms, fillnominal );
  }
  if ( GetVerbose()>=2 ) {
    std::cout << "End of SampleManager::FillSampleBins" << std::endl;
    std::cout << "--------------------------------" << std::endl;
  }
}

void SampleManager::MakeFakeDataSet( ParameterManager* systerms, bool includeStatVars, bool includePullVars ) {
  std::cout << "========================================================================" << std::endl;
  std::cout << "[SampleManager::MakeFakeDataSet]" << std::endl;
  std::cout << "stat=" << includeStatVars << " pull=" << includePullVars << std::endl;

  std::map< std::string, double > parvalues;
  if ( includePullVars ) {
    systerms->GenerateRandomValues( parvalues );
    std::cout << "Pull Values Set to:" << std::endl;
    std::cout << "--------------------"  << std::endl;
    std::cout << "PARAMETER : VALUE : DISTANCE FROM CENTER (in sigma)" << std::endl;
    for ( std::map< std::string, double >::iterator it=parvalues.begin(); it!=parvalues.end(); it++ )
      std::cout << (*it).first 
		<< " : " << (*it).second 
		<< " : " << ((*it).second-systerms->GetParameter( (*it).first )->GetCentralValue() )/systerms->GetParameter( (*it).first )->GetSigma() 
		<< std::endl;
  }
  
  for ( SampleListIter it=SampleListBegin();it!=SampleListEnd(); it++ )
    GetSample( (*it) )->MakeFakeDataSet( systerms, includeStatVars, parvalues );

  std::cout << "========================================================================" << std::endl;
  
}

void SampleManager::SetSeed( int seed ) {
  for ( SampleListIter it=SampleListBegin();it!=SampleListEnd(); it++ ) {
    GetSample( (*it) )->SetSeed( seed );
  }
}

// --------------------------------------------------------------------------------------------------------------------------------
// Interface to Tree Output

void SampleManager::AddConfigurationTreeToFile( TFile* afile ) {

}

void SampleManager::AddExpectedBinsToTree( TTree* atree ) {

  // assign bin addresses to the tree
  // we do this instead of creating a mirror of values and trying to sync them
  // of course, now we've opened the way to pointer-address-magedon
  for ( int gbin=0; gbin<GetNumberOfBins(); gbin++) {
    Sample* asample = GetSampleFromBinIndex( gbin );
    double* pNbin = asample->GetExpectedBinAddressFromGlobalIndex( gbin );
    atree->Branch( GetExpectedBinBranchName(gbin).c_str(), pNbin, std::string( GetExpectedBinBranchName(gbin)+"/D" ).c_str() );
  }
  
}

void SampleManager::AddObservedBinsToTree( TTree* atree ) {
  
  // assign bin addresses to the tree
  // we do this instead of creating a mirror of values and trying to sync them
  // of course, now we've opened the way to pointer-address-magedon
  for ( int gbin=0; gbin<GetNumberOfBins(); gbin++) {
    Sample* asample = GetSampleFromBinIndex( gbin );
    double* pNbin = asample->GetObservedBinAddressFromGlobalIndex( gbin );
    atree->Branch( GetObservedBinBranchName(gbin).c_str(), pNbin, std::string( GetObservedBinBranchName(gbin)+"/D" ).c_str() );
  }
  
}

std::string SampleManager::GetExpectedBinBranchName( int global_index ) {
  std::stringstream branch_name;
  branch_name << "Nexp_" << GetSampleFromBinIndex( global_index )->GetName() << "_" << global_index;
  return branch_name.str();
}

std::string SampleManager::GetObservedBinBranchName( int global_index ) {
  std::stringstream branch_name;
  branch_name << "Nobs_" << GetSampleFromBinIndex(global_index)->GetName() << "_" << global_index;
  return branch_name.str();
}

// --------------------------------------------------------------------------------------------------------------------------------
// Misc Tools

void SampleManager::WriteHistograms() {
  for ( SampleDictIter it=m_sample_dict.begin(); it!=m_sample_dict.end(); it++ ) {
    Sample* asample = (*it).second;
    if ( asample ) asample->WriteHistogram();
  }
}

bool SampleManager::LoadSamplesFromFile( TFile* file, double scalefactor) {
  bool samplesok = true;
  for ( SampleDictIter it=SampleDictBegin(); it!=SampleDictEnd(); it++ ) {
    (*it).second->SetVerbose( GetVerbose() );
    bool loadedok = (*it).second->LoadSampleFromFile( file, scalefactor );
    if (!loadedok) {
      std::cout << "SampleManager::LoadSamplesFromFile -- WARNING!! Failed to load the sample, '" << (*it).first << "', from file" << std::endl;
      samplesok = false;
    }
  }
  return samplesok;
}

void SampleManager::Print() {
  std::cout << "========================================================" << std::endl;
  std::cout << "SampleManager: Sample Info" << std::endl;
  if ( !fSampleOrdered ) {
    std::cout << "Samples Unordered" << std::endl;
    for ( SampleDictIter it=SampleDictBegin(); it!=SampleDictEnd(); it++ ) {
      (*it).second->PrintSampleInfo();
    }
  }
  else {
    std::cout << "Samples ordered." << std::endl;
    for ( SampleListIter it=SampleListBegin(); it!=SampleListEnd(); it++ ) {
      GetSample( (*it) )->PrintSampleInfo();
    }
  }
  
  std::cout << "========================================================" << std::endl;
}

void SampleManager::PrintBins() {
  std::cout << "========================================================" << std::endl;
  std::cout << "SampleManager: Print Bin Info" << std::endl;
  if ( !fSampleOrdered ) {
    std::cout << "Sample Bins not yet indexed" << std::endl;
  }
  else {
    for ( SampleListIter it=SampleListBegin(); it!=SampleListEnd(); it++ ) {
      GetSample( (*it) )->PrintSampleInfo();
    }
  }
  std::cout << "========================================================" << std::endl;
}
