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

#include "ResponseFunctionParameter.hh"
#include <sstream>

#include "ResponseFunction.hh"
#include "ResponseFunctionManager.hh"
#include "SampleManager.hh"

using namespace qosc;

std::string qosc::generate_standard_response_name( std::string samplename, int recon_bin_index, int truth_bin_index ) {
  return ResponseFunctionParameter::GenerateStandardResponseName( samplename, recon_bin_index, truth_bin_index );
}

ResponseFunctionParameter::ResponseFunctionParameter( std::string name, double mean, double sigma, FitterParType t ) 
  : BasicParameter( name, mean, sigma, t ) {
  m_funcMan = NULL;
}


ResponseFunctionParameter::~ResponseFunctionParameter() {
}

void ResponseFunctionParameter::DefineSampleAndResponseStructure( SampleManager* sampleman ) {
  std::cout << __PRETTY_FUNCTION__ << "[" << GetName() << "]" <<  std::endl;
  m_funcMan = new ResponseFunctionManager( sampleman ); // manager tracks functions for a parameter for a given bin info
}

bool ResponseFunctionParameter::IsResponseLoaded( std::string samplename, int rec_bin, int truth_bin ) {
  ResponseFunction* response = GetResponseFunction( samplename, rec_bin, truth_bin );
  if ( response ) return true;
  else return false;
}

ResponseFunction* ResponseFunctionParameter::GetResponseFunction( int sampleid, int recon_bin_index, int truth_bin_index ) {
  return m_funcMan->GetResponseFunction( sampleid, recon_bin_index, truth_bin_index );
}

ResponseFunction* ResponseFunctionParameter::GetResponseFunction( std::string samplename, int recon_bin_index, int truth_bin_index ) {
  return m_funcMan->GetResponseFunction( samplename, recon_bin_index, truth_bin_index );
}

double ResponseFunctionParameter::GetResponse( std::string samplename, int recon_bin_index, int truth_bin_index, double x ) {
  double pull[1] = { x };
  ResponseFunction* r = GetResponseFunction( samplename, recon_bin_index, truth_bin_index );
  if ( !r )
    return 1.0;
  return r->GetResponse( pull );
}

double ResponseFunctionParameter::GetResponse( int sampleid, int recon_bin_index, int truth_bin_index, double x ) {
  double pull[1] = { x };
  ResponseFunction* r = GetResponseFunction( sampleid, recon_bin_index, truth_bin_index );
  if ( !r ) {
    //std::cout << "No response! " << sampleid << " " << recon_bin_index << " " << truth_bin_index << " x=" << x << std::endl;
    //std::cin.get();
    return 1.0;
  }
  return r->GetResponse( pull );
}

int ResponseFunctionParameter::GetNumberOfTruthBins( std::string samplename ) {
  return GetResponseFunctionManager()->GetNumberOfTruthBins( samplename );
}

int ResponseFunctionParameter::GetNumberOfReconBins( std::string samplename ) {
  return GetResponseFunctionManager()->GetNumberOfReconBins( samplename );
}

std::string ResponseFunctionParameter::GenerateStandardResponseName( std::string samplename, int recon_bin_index, int truth_bin_index ) {
  std::stringstream ss;
  ss << samplename << "_" << recon_bin_index << "_" << truth_bin_index;
  return ss.str();
}

int ResponseFunctionParameter::GetSampleID( std::string sample ) {
  return m_funcMan->GetSampleID( sample );
}
