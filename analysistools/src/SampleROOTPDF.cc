#include "SampleROOTPDF.hh"
#include <iostream>
#include <assert.h>
#include <vector>

using namespace qosc;

SampleROOTPDF::SampleROOTPDF( std::string sampleName, std::string chainName ) 
  : Sample( sampleName )
{
  m_chain_name = chainName;
}

SampleROOTPDF::~SampleROOTPDF() {
};

void SampleROOTPDF::FillBins( ParameterManager* systerms, bool fillnominal ) { 
  std::cout << "Sample::FillBins called. This is unimplemented." << std::endl;
  assert(false);
}

void SampleROOTPDF::PrintSampleInfo() {
  Sample::PrintSampleInfo();
  std::cout << "Source chain: " << GetChainName() << std::endl;
}
