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

#include "BinShiftParameter.hh"
#include "TH1D.h"
#include "ModelParameter.hh"

using namespace qosc;

BinShiftParameter::BinShiftParameter( std::string name, std::string sample_to_apply, double mean, double sig ) 
  : BasicParameter( name, mean, sig, kMinuitTerm )
{
  output_hist = NULL;
  AddAffectedSample( sample_to_apply );
}

BinShiftParameter::~BinShiftParameter() {}

void BinShiftParameter::AddAffectedSample( std::string sample ) {
  m_sample_applies.insert( sample );
}

bool BinShiftParameter::DoesSampleApply( std::string sample ) { 
  if ( m_sample_applies.find( sample )!=m_sample_applies.end() )
    return true;
  return false;
}

void BinShiftParameter::TransformHistogram( TH1D* h ) {
  double pvalue = GetValue();
  if ( pvalue==0.0 ) {
    // no change, just return.
    return;
  }

  if ( output_hist==NULL ) {
    output_hist = (TH1D*)h->Clone( std::string(h->GetName()+ModelParameter::GetName() ).c_str() );
  }
  output_hist->Reset();


  for (int ibin=0; ibin<h->GetNbinsX(); ibin++) {
    double ledge = h->GetBinLowEdge( ibin+1 );
    double hedge = h->GetBinLowEdge( ibin+2 );
    double nevents = 0.;
    for (int obin=0; obin<h->GetNbinsX(); obin++) {
    
      double ledge_cor = h->GetBinLowEdge( obin+1 )*(1+pvalue);
      double hedge_cor = h->GetBinLowEdge( obin+2 )*(1+pvalue);
      double width_cor = hedge_cor - ledge_cor;
      double rfrac = 0.;
      if ( hedge<=ledge_cor ) break;
      else if ( ledge>=hedge_cor ) continue;
      else if ( ledge < ledge_cor && hedge <= hedge_cor ) {
	// new bin higher, partially overlapping
	rfrac = ( hedge-ledge_cor )/width_cor;
      }else if ( ledge >= ledge_cor && hedge>hedge_cor ) {
	// new bin lower, partially overlapping
	rfrac = ( hedge_cor-ledge )/width_cor;
      }
      else if ( ledge<ledge_cor && hedge > hedge_cor ) {
	// new bin is subset of original bin
	rfrac = 1.0;
      }
      else if ( ledge>=ledge_cor && hedge<=hedge_cor ) {
	// new bin is superset of original bin
	rfrac = (hedge-ledge)/width_cor;
      }

      nevents += rfrac*h->GetBinContent( obin+1 );
    }
    output_hist->SetBinContent( ibin, nevents );
  }
  
  for (int bin=0; bin<h->GetNbinsX(); bin++)
    h->SetBinContent( bin+1, output_hist->GetBinContent( bin ) );
  
}

