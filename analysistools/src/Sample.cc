#include "Sample.hh"
#include <assert.h>
#include <iostream>

#include "ParameterManager.hh"
#include "Hist.hh"

using namespace qosc;

Sample::Sample( std::string name ) {
  histogram = NULL;
  datahist = NULL;
  fSampleActive = true;
  m_sample_name = name;
  m_seed = 0;
  m_rand_gen = new TRandom3( m_seed );
  m_SampleID = -1;
  SetVerbose(0);
}

Sample::~Sample() {
  delete histogram;
}

// --------------------------------------------------------------------------------
// Histogram Management

void Sample::SetSampleHistogram( Hist* hist ) {
  m_histbin_active.clear();
  // create active bin flags for histogram
  // we list bins zero indexed. 
  // we also make no assumption about under or overflow. user defines this by setting pdfbinned::UsseOverUnderFlow( use );
  for (int bin=0; bin<=hist->GetNumberOfBins(); bin++)
    m_histbin_active[bin] = true;

  // store expectation histogram
  histogram = hist;

  // make mirror histogram for the data
  std::string name = hist->GetName()+"_data";
  if ( datahist ) delete datahist;
  if ( hist->GetHistDimensions()==1 ) {
    TH1D* newhist = (TH1D*)hist->GetHistogram()->Clone( name.c_str() );
    datahist = new Hist( name );
    datahist->SetHistogram( newhist );
  }
  else if ( hist->GetHistDimensions()==2 ) {
    TH2D* newhist = (TH2D*)hist->GetHistogram()->Clone( name.c_str() );
    datahist = new Hist( name );
    datahist->SetHistogram( newhist );
  }
}

void Sample::CheckHistogram() {
  // sanity checks
  if ( !histogram ) {
    std::cout << "Sample::IndexBins -- ERROR! The Sample histogram has not been defined yet!" << std::endl;
    assert(false);
  }
  if ( histogram->GetHistogramStatus()==Hist::kUndefined ) {
    std::cout << "Sample::IndexBins -- ERROR! The Sample histogram has not been defined yet!" << std::endl;
    assert(false);    
  }
}

// --------------------------------------------------------------------------------
// -------------------
// Bin Management
// -------------------

void Sample::IndexBins() {
  // This is where the sample takes whatever the user has created and formats/indexes for use by the analysis class.
  CheckHistogram();

  fNumberOfSampleBins = 0;
  int binindex = 0;
  int histbins = GetNumberOfHistogramBins();
  for (int n=0; n<histbins; n++) {
    if ( IsHistBinActive( n ) ) {
      fNumberOfSampleBins++;
      m_local_to_hist[binindex] = n;
      binindex++;
    }
  }

  assert( (unsigned int)fNumberOfSampleBins==m_local_to_hist.size() ); // another sanity check
}

// --------------------------------------------------------------------------------
// ----------------
// Bin Retrieval
// ----------------

double Sample::GetContentByHistIndex( int histindex, Hist* hist ) {
  CheckHistogram();  
  // Note to self. The Hist is going to need get bin content functions to enforce choice of bin numbering whatever it is.
  int ndims = hist->GetHistDimensions();
  if ( ndims==1 )
    return hist->GetBinContent( histindex );
  else if (ndims==2) {
    int binx, biny;
    hist->GetBinXYFromIndex( histindex, binx, biny );
    return hist->GetBinContent( binx, biny );
  }
  else
    assert(false);
}

double Sample::GetContentByLocalIndex( int localindex, Hist* hist ) {
  BinIndexDictIter it = m_local_to_hist.find( localindex );
  if ( it==m_local_to_hist.end() ) {
    std::cout << "Sample::GetContentByLocalIndex -- ERROR. Did not find a map from local index (" << localindex << ") to hist index" << std::endl;
    assert(false);
  }
  return GetContentByHistIndex( (*it).second, hist );
}

double Sample::GetContentByGlobalIndex( int globalindex, Hist* hist ) {
  BinIndexDictIter it = m_global_to_local.find( globalindex );
  if ( it==m_global_to_local.end() ) {
    std::cout << "Sample::GetContentByGlobalIndex -- ERROR. Did not find a map from global index (" << globalindex << ") to local index" << std::endl;
    assert(false);
  }
  return GetContentByLocalIndex( (*it).second, hist );
}

// Set Functions

void Sample::SetContentByHistIndex( int histindex, double value, Hist* hist ) {
  CheckHistogram();  
  // Note to self. The Hist is going to need get bin content functions to enforce choice of bin numbering whatever it is.
  int ndims = hist->GetHistDimensions();
  if ( ndims==1 )
    return hist->SetBinContent( histindex, value );
  else if (ndims==2) {
    int binx, biny;
    hist->GetBinXYFromIndex( histindex, binx, biny );
    return hist->SetBinContent( binx, biny, value );
  }
  else
    assert(false);
}

void Sample::SetContentByLocalIndex( int localindex, double value, Hist* hist ) {
  BinIndexDictIter it = m_local_to_hist.find( localindex );
  if ( it==m_local_to_hist.end() ) {
    std::cout << "Sample::GetContentByLocalIndex -- ERROR. Did not find a map from local index (" << localindex << ") to hist index" << std::endl;
    assert(false);
  }
  return SetContentByHistIndex( (*it).second, value, hist );
}

void Sample::SetContentByGlobalIndex( int globalindex, double value, Hist* hist ) {
  BinIndexDictIter it = m_global_to_local.find( globalindex );
  if ( it==m_global_to_local.end() ) {
    std::cout << "Sample::GetContentByGlobalIndex -- ERROR. Did not find a map from global index (" << globalindex << ") to local index" << std::endl;
    assert(false);
  }
  return SetContentByLocalIndex( (*it).second, value, hist );
}


// --------------------------------------------------------------------------------
// -------------------
// Data Management
// -------------------
void Sample::SetSeed(int seed ) {
  m_seed = seed;
  m_rand_gen->SetSeed( m_seed );
}

void Sample::MakeFakeDataSet( ParameterManager* systerms, bool includeStatVariation, std::map< std::string, double >& parvalues ) {
  
  std::cout << "Sample[" << GetName() << "]::MakeFakeDataSet" 
	    << "(stat=" << includeStatVariation << ", pull=" << parvalues.size() << ")"
	    << std::endl;
  assert( systerms!=NULL );
  
  if ( parvalues.size()!=0 ) {
    // include pull variations
    for (std::map< std::string, double>::iterator it=parvalues.begin(); it!=parvalues.end(); it++) {
      if ( systerms->GetParameterTypeName( (*it).first )=="EventReweightParameter" 
	   || systerms->GetParameterTypeName( (*it).first )=="FijParameter" 
	   || systerms->GetParameterTypeName( (*it).first )=="ResponseCurveParameter" )
	systerms->SetParameterValue( (*it).first, (*it).second );
    }
  }

  // Update Sample Hist
  FillBins( systerms );

  // create fake data with stat. variations
  for (int b=0; b<GetNumberOfHistogramBins(); b++) {
    double binvalue = GetExpectedContentByHistIndex( b );
    if ( includeStatVariation ) binvalue = (double)m_rand_gen->Poisson( binvalue );
    datahist->SetBinContent( b, binvalue );
  }
}

void Sample::MakeFakeDataSet( ParameterManager* systerms, bool includeStatVariation, bool includePullVariation ) {
  // generate data for each histogram bin
  if ( GetVerbose()>=2 ) std::cout << "Sample[" << GetName() << "]::MakeFakeDataSet" << std::endl;
  assert( systerms!=NULL );

  std::map< std::string, double > parvalues;
  if ( includePullVariation )
    systerms->GenerateRandomValues( parvalues );
  MakeFakeDataSet( systerms, includeStatVariation, parvalues );
}

Double_t* Sample::GetExpectedBinAddressFromGlobalIndex( int global_index ) {
  return GetSampleHistogram()->GetBinAddress( GetLocalIndexFromGlobal( global_index ) );
}


Double_t* Sample::GetExpectedBinAddressFromLocalIndex( int local_index ) {
  return GetSampleHistogram()->GetBinAddress( local_index );
}

Double_t* Sample::GetObservedBinAddressFromGlobalIndex( int global_index ) {
  return GetDataHistogram()->GetBinAddress( GetLocalIndexFromGlobal( global_index ) );
}


Double_t* Sample::GetObservedBinAddressFromLocalIndex( int local_index ) {
  return GetDataHistogram()->GetBinAddress( local_index );
}

// --------------------------------------------------------------------------------
// ------------------------
// User Bin Info Functions
// ------------------------
void Sample::GetListOfUserBinInfoNames( std::vector< std::string >& bin_info_names ) {
  Hist* pdf = GetSampleHistogram();
  for ( UserBinInfoListIter it=pdf->MasterBinInfoListBegin(); it!=pdf->MasterBinInfoListEnd(); it++ ) {
    bin_info_names.push_back( (*it).first );
  }
}

UserBinInfo* Sample::GetMasterListBinInfoInstance( std::string bin_info_name ) {
  Hist* pdf = GetSampleHistogram();
  return pdf->GetMasterListBinInfo( bin_info_name );
}

// --------------------------------------------------------------------------------
// -------------------
// Misc Tools
// -------------------

void Sample::WriteHistogram() {
  if ( histogram ) {
    histogram->GetHistogram()->Write();  /// this is the main histogram
    if ( histogram->GetStoredHistogram() ) histogram->GetStoredHistogram()->Write(); /// this is the 'stored' histogram used as the reference histogram for reweighting by oscillatin prob.
    histogram->WriteBinInfo(); /// this will cause the histograms in the UserBinInfo class (here OscBinInfo) to write to file.
    histogram->WriteHistInfo();
  }
}

void Sample::PrintSampleInfo( bool printBinInfo ) {
  std::cout << "========================================================================================" << std::endl;
  std::cout << "-- " << GetName() << " (" << this << ")" << " --" << std::endl;
  std::cout << "Sample stauts: " << ( ( fSampleActive ) ? "Active" : "Inactive" ) << std::endl;
  std::cout << "Histogram status: ";
  if ( !histogram ) {
    std::cout << "NULL" << std::endl;
    return;
  }
  else std::cout << ( ( histogram->GetHistogramStatus() ) ? " Defined " : " Undefined " ) << " (" << histogram << ")" << std::endl;
  std::cout << "Number of bins: " << GetNumberOfHistogramBins() << ", Number of active/indexed bins: " << GetNumberOfBins() << std::endl;
  for (int bin=0; bin<histogram->GetNumberOfBins(); bin++) {
    if ( printBinInfo ) std::cout << "----------------------------------------------------------------------------------------" << std::endl;
    std::cout << " [" << bin << "] (" << GetName() << ") ";
    if ( histogram->GetHistDimensions()==1 ) {
      std::cout << "[" << ((TH1D*)histogram->GetHistogram())->GetBinLowEdge(bin+1) << ", " << ((TH1D*)histogram->GetHistogram())->GetBinLowEdge(bin+2)  << "]";
    }
    else if ( histogram->GetHistDimensions()==2 ) {
      int binx, biny;
      histogram->GetBinXYFromIndex( bin, binx, biny );
      std::cout << "["
		<< "x=" << ((TH2D*)histogram->GetHistogram())->GetXaxis()->GetBinLowEdge(binx+1) << "," << ((TH2D*)histogram->GetHistogram())->GetXaxis()->GetBinLowEdge(binx+2)  
		<< ";"
		<< "y=" << ((TH2D*)histogram->GetHistogram())->GetYaxis()->GetBinLowEdge(biny+1) << "," << ((TH2D*)histogram->GetHistogram())->GetYaxis()->GetBinLowEdge(biny+2)  
		<< "]";
    }
    else {
      assert(false);
    }
    std::cout << ": Observed=" << datahist->GetBinContent( bin ) << " Expected= " <<  histogram->GetBinContent( bin ) << std::endl;
    if ( printBinInfo ) histogram->GetBinInfoList( bin )->Print();
  }
  std::cout << "------------------------------------------------------------------------" << std::endl;
  std::cout << "Sample Integral: Observed=" << datahist->GetHistogram()->Integral() << " Expected=" << histogram->GetHistogram()->Integral() << std::endl;
  std::cout << "========================================================================================" << std::endl;
}

