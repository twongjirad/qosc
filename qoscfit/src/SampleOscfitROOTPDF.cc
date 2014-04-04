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

#include "SampleOscfitROOTPDF.hh"
#include <iostream>
#include <assert.h>

#include "TFile.h"

#include "Hist.hh"
#include "HistRootVariable.hh"
#include "TruthBinInfo.hh"
#include "ResponseFunction.hh"
#include "ResponseFunctionParameter.hh"
#include "ResponseFunctionManager.hh"
#include "BinShiftParameter.hh"
#include "SampleSetupParser.hh"
#include "WeightOscEvent.hh"
#include "WeightBargerOsc.hh"
#include "WeightBargerTruthHist.hh"

#include <vector>

using namespace qosc;

SampleOscfitROOTPDF::SampleOscfitROOTPDF( std::string sampleName, std::string chainName ) 
  : SampleROOTPDF( sampleName, chainName )
{
  m_chain_name = chainName;
  original_histogram = NULL;
  predistortion_histogram = NULL;
  oscillator = NULL;
  event_oscillator = NULL;
  Fij = NULL;
  m_sampleid_table = NULL;
  SampleROOTPDF::SetScalingFactor( 1.0 );
  fCheckedForShiftPars = false;
  m_shiftpar_array = NULL;
  fSaveTruthModeHists = false;
  truthmode_histograms = NULL;
  fUseEventReweightOsc = false; // Only expect to use this during checks of approximation model and for T2K spline making
  m_osc_bin_info_name = ""; // name of userbininfo object that will store true energy information used to oscillation the bins
  m_sample_truth_to_splinebin = NULL; // dictioanry between osc binning and respone binnig
  fUseDifferentOscAndResponseBinning = false;  
  fTruthToResponseMapMade = false;
}

SampleOscfitROOTPDF::~SampleOscfitROOTPDF() {
  if ( original_histogram ) delete original_histogram;
  if ( predistortion_histogram ) delete predistortion_histogram;
  if ( oscillator ) delete oscillator;
  if ( event_oscillator ) delete event_oscillator;
  if ( Fij ) delete [] Fij;
  if ( m_sampleid_table ) delete [] m_sampleid_table;
  if ( m_shiftpar_array ) delete [] m_shiftpar_array;
  if ( m_sample_truth_to_splinebin ) delete [] m_sample_truth_to_splinebin;
}


void SampleOscfitROOTPDF::SetSampleHistogram( HistRootVariable* hist ) {
  SampleROOTPDF::SetSampleHistogram( hist );
  if ( m_osc_bin_info_name=="" ) {
    std::cout << "OSC BIN INFO NAME NOT SPECIFIED" << std::endl;
    std::cout << "Must first call SampleOscfitROOTPDF::SetOscBinInfoName(...)" << std::endl;
    assert(false);
  }
  UserBinInfo* bininfo = hist->GetMasterListBinInfo( m_osc_bin_info_name ); // we must have already set this string
  if ( bininfo ) {
    m_truthbins = (TruthBinInfo*)bininfo->Copy( m_osc_bin_info_name+"_"+GetName() );
    m_truthbins_fij = (TruthBinInfo*)bininfo->Copy( m_osc_bin_info_name+"_"+GetName()+"_fij" );
    oscillator = new WeightBargerTruthHist( m_truthbins ); 
    predistortion_histogram = new Hist( "predist_" + GetName() );
    if (GetSampleHistogram()->GetHistDimensions()==1)
      predistortion_histogram->SetHistogram( (TH1D*)hist->GetHistogram()->Clone( std::string("predist_" + GetName()).c_str() ) );
    else
      predistortion_histogram->SetHistogram( (TH2D*)hist->GetHistogram()->Clone( std::string("predist_" + GetName()).c_str() ) );
  }
  else {
    std::cout << "NO OSC BIN INFO FOUND FOR " << GetName() << std::endl;
    assert(false);
  }
}

void SampleOscfitROOTPDF::SetSinThetaForm( SinThetaForms sintheta ) {
  oscillator->SetSinThetaForm( sintheta );
}

void SampleOscfitROOTPDF::FillBins( ParameterManager* parameters, bool fillnominal ) {

  // Set the oscillation weights
  if ( GetVerbose()>=2 ) {
    std::cout << "SampleOscfitROOTPDF::FillBins - Updating osc calculator" << std::endl;
    oscillator->SetVerbose(1);
  }
  UpdateOscPars( parameters );
  
  // Allocate space to store truth information
  Hist* pre_hist = GetPredistortionHist();
  double subtotals_preosc[ m_truthbins->GetNumberOfTruthHists() ];
  double subtotals_postosc[ m_truthbins->GetNumberOfTruthHists() ];
  memset( subtotals_preosc, 0, sizeof(double)*m_truthbins->GetNumberOfTruthHists() );
  memset( subtotals_postosc, 0, sizeof(double)*m_truthbins->GetNumberOfTruthHists() );

  /// --------------------------------------------------------------------------------------------------------------
  /// RESPONSE CALCULATIONS FOR EACH TRUTH BIN
  /// --------------------------------------------------------------------------------------------------------------
  // Get the Sample ID for this sample in order to look up response parameters. 
  // This is done in order to avoid string comparisons and using map lookup table.
  // Instead the id is a position in an array which should be much faster to access
  if ( !m_sampleid_table )
    StoreResponseParSampleID( parameters );

  // Begin the calculation of the expectation in each reconstructed energy bin
  for ( int bin=0; bin<GetSampleHistogram()->GetNumberOfBins(); bin++) {

    // Get the user bin info that is storing the truth bin template
    TruthBinInfo* true_energy_hists = dynamic_cast<TruthBinInfo*>( GetSampleHistogram()->GetBinInfo( m_osc_bin_info_name, bin ) ); // truth bins from the fill

    int ntruebins = 0; // truth bin counter

    // begin loop over the different true energy bins [they are split into flux and modes]
    for (int j=0; j<true_energy_hists->GetNumberOfTruthHists(); j++) {

      TH1D* intruhist = (TH1D*)true_energy_hists->GetTruthHist( j );

      // reset workspace array 
      // (reason we use a histogram object instead of simple array is that it allows us to pass the information to the oscillator in easy step later. )
      // (hopefuly the overhead is not terrible)
      TH1D* outtruhist = (TH1D*)m_truthbins->GetTruthHist( j );
      outtruhist->Reset();

      // if no events in this component, just skip.
      if ( intruhist->Integral()==0 ) {
	ntruebins += intruhist->GetNbinsX();
	continue;
      }

      //std::cout << "in=" << intruhist->GetName() << "  out=" << outtruhist->GetName() << std::endl;

      for (int itrue=0; itrue<intruhist->GetNbinsX(); itrue++) {
	double ntruevents = intruhist->GetBinContent( itrue+1 );
	ntruebins++;
	if ( ntruevents<=0 ) continue;
	double response = 1; // use response from all parameters

	for (int ipar=0; ipar<parameters->NumberOfParameters(); ipar++) {
	  if ( parameters->GetParameterFromIndex( ipar )->IsA()=="ResponseFunctionParameter" ) {

	    // When asked to fill nominal, we ignore newton terms as their effect will be passed to the fitter by Fij terms
	    if ( fillnominal && parameters->GetParameterFromIndex( ipar )->GetFitterParType()==kNewtonTerm ) {
// 	      std::cout << "Skipping " << parameters->GetParameterFromIndex( ipar )->GetName() << " nominal=" << fillnominal 
// 			<< " partype=" << ((parameters->GetParameterFromIndex( ipar )->GetFitterParType()==kNewtonTerm) ? "Newton" : "Minuit" )  << std::endl;
// 	      std::cin.get();
	      continue;
	    }

	    int ntruthbin_response = GetResponseTruthBin( ntruebins-1 );  // we allow the ability to have the binning for the response function bins and oscillation bins to be different.
	    double parvalue[1] = { parameters->GetParameterFromIndex( ipar )->GetValue() };
	    //double parresponse = ((ResponseFunctionParameter*)parameters->GetParameterFromIndex( ipar ))->GetResponse( Sample::GetName(), bin, j*intruhist->GetNbinsX()+itrue, parvalue[0] );
	    //double parresponse = ((ResponseFunctionParameter*)parameters->GetParameterFromIndex( ipar ))->GetResponse( Sample::GetName(), bin, ntruebins-1, parvalue[0] );
	    double parresponse = ((ResponseFunctionParameter*)parameters->GetParameterFromIndex( ipar ))->GetResponse( m_sampleid_table[ipar], bin, ntruthbin_response, parvalue[0] );
	    if ( parresponse!=1.0 && parresponse==parresponse) {
// 	      std::cout << "  " << parameters->GetParameterFromIndex( ipar )->GetName() << " response: " << parresponse 
// 			<< " irec=" << bin << " itru=" << itrue << " (" << ntruebins-1 << ", " << true_energy_hists->GetCutNameFromTruthBin( ntruthbin_response ) << ")" << std::endl;
// 	      std::cin.get();
	      //response *= parresponse;  
	      response += (parresponse-1);
	    }

// 	    if ( itrue==36 ) {
// 	      std::cout << "dN[" << bin << ", " << j << ", " << itrue << "]/sigma = " << ntruevents*(parresponse-1.0)/parvalue[0] << std::endl;
// 	    }
   
	  }
	}
	if ( response>0 )
	  outtruhist->SetBinContent( itrue+1, ntruevents*response );
	else
	  outtruhist->SetBinContent( itrue+1, 0.0 );
      }
//       if ( bin<6 ) {
// 	std::cout << "pre-osc xtruthhist[" << j << "] "
// 		  << " in(" << intruhist->GetName() << ")=" << intruhist->Integral()*SampleROOTPDF::GetScalingFactor() 
// 		  << " out(" << outtruhist->GetName() << ")=" << outtruhist->Integral()*SampleROOTPDF::GetScalingFactor() << std::endl;
//       }
      subtotals_preosc[j] += outtruhist->Integral(); /// pre osc
      
    }//end of loop over true hists
    
    if ( !fUseEventReweightOsc ) {
      // now that truth histograms are now prepared, we can oscillate them and include the scaling factor.
      // but we only do this if we aren't oscillating the expectation through the event-by-event mechanism.
      double nevents = oscillator->CalcEventsOsc( m_truthbins, true );
      pre_hist->SetBinContent( bin, nevents*SampleROOTPDF::GetScalingFactor() );
    }
    else {
      //the truth histograms should have already been oscillated when filling event-by-event
      //std::cout << "oscillation skipped due to event-by-event oscillations." << std::endl;
      double nevents = 0;
      for ( int j=0; j<true_energy_hists->GetNumberOfTruthHists(); j++) {
	nevents +=  m_truthbins->GetTruthHist( j )->Integral();
      }
      pre_hist->SetBinContent( bin, nevents*SampleROOTPDF::GetScalingFactor() );
    }
    
    // now we save some info on the accumulated totals for later output.
    // we can get rid of this later once we begin to look for optimizations.
    for ( int j=0; j<true_energy_hists->GetNumberOfTruthHists(); j++) {
      subtotals_postosc[j] += m_truthbins->GetTruthHist( j )->Integral();
      if ( fSaveTruthModeHists ) {
	truthmode_histograms[j]->GetHistogram()->SetBinContent( bin+1, truthmode_histograms[j]->GetHistogram()->GetBinContent(bin+1) + m_truthbins->GetTruthHist( j )->Integral() );
      }
    }
    
  }//end of loop over rec bins
  
  // dump out some information: removable when looking for speed-up
  double preosc_total = 0.;
  double postosc_total = 0.;
  for (int i=0; i<m_truthbins->GetNumberOfTruthHists(); i++) {
//     std::cout << m_truthbins->GetTruthHist(i)->GetName() 
// 	      << ": pre-osc=" << subtotals_preosc[i]*SampleROOTPDF::GetScalingFactor() 
// 	      << "  post-osc=" << subtotals_postosc[i]*SampleROOTPDF::GetScalingFactor() 
// 	      << std::endl;
    preosc_total += subtotals_preosc[i];
    postosc_total += subtotals_postosc[i];
  }
  std::cout << __PRETTY_FUNCTION__ << " - Total:"
	    << " pre-osc=" << preosc_total*SampleROOTPDF::GetScalingFactor() 
	    << " post-osc=" << postosc_total*SampleROOTPDF::GetScalingFactor() << std::endl;

  // ---------------------------------------------------------------------------
  // Apply Shift Variables
  if ( !fCheckedForShiftPars ) {
    // The first time this is run, we must make a list of bin shift parameters that apply to this sample.
    m_nshiftpars = 0;
    std::vector< BinShiftParameter* > shiftparlist;
    for (int ipar=0; ipar<parameters->NumberOfParameters(); ipar++) {
      if ( parameters->GetParameterFromIndex( ipar )->IsA()=="BinShiftParameter" 
	   && ((BinShiftParameter*)parameters->GetParameterFromIndex( ipar ))->DoesSampleApply( GetName() ) ) {
	shiftparlist.push_back( (BinShiftParameter*)parameters->GetParameterFromIndex( ipar ) );
	m_nshiftpars++;
      }
    }
    m_shiftpar_array = new BinShiftParameter*[m_nshiftpars];
    for (int n=0; n<m_nshiftpars; n++) {
      m_shiftpar_array[n] = shiftparlist.at(n);
    }
    fCheckedForShiftPars = true;
  }
  for (int n=0; n<m_nshiftpars; n++) {
    m_shiftpar_array[n]->TransformHistogram( (TH1D*)pre_hist->GetHistogram() );
    std::cout << "  after binshift par[" << n << "]: " << pre_hist->GetHistogram()->Integral() << std::endl;
  }

  // ---------------------------------------------------------------------------
  // Finally, Transfer final answer to sample bins
  for ( int bin=0; bin<GetSampleHistogram()->GetNumberOfBins(); bin++) {
    GetSampleHistogram()->SetBinContent(bin, pre_hist->GetBinContent(bin) );
  }
  // ---------------------------------------------------------------------------

}


void SampleOscfitROOTPDF::CalculateFijValues( ParameterManager* parameters ) {
  
  // Set the oscillation weights
  //if ( GetVerbose()>=2 ) {
  std::cout << __PRETTY_FUNCTION__ << std::endl;
  //}

  UpdateOscPars( parameters );

  if ( !m_sampleid_table )
    StoreResponseParSampleID( parameters );

  
  if ( Fij==NULL ) {
    // Allocate Fij array
    m_fij_npars = parameters->NumberOfParameters();
    m_fij_nelems = GetSampleHistogram()->GetNumberOfBins()*parameters->NumberOfParameters();
    Fij = new double[ m_fij_nelems ];
    memset( Fij, 0, sizeof(double)*m_fij_nelems );
  }
    

  // There are two sets of truth binnings
  // There is a set of truth bin info objects in each reconstructed bin. Set A.
  // There is an additional set of truth hists in the m_truthbins_fij object. Set B.
  // Set A will be filled using stored histograms from file, or by event-by-event filling from a root ttree object.
  // Set B uses the info from Set A and stores the respone function modified values.
  //   Note that the data in Set B is transitory. It resets every time we calculate the expectation in a recon bin.

  for ( int bin=0; bin<GetSampleHistogram()->GetNumberOfBins(); bin++) {
    TruthBinInfo* true_energy_hists = dynamic_cast<TruthBinInfo*>( GetSampleHistogram()->GetBinInfo( m_osc_bin_info_name, bin ) ); // truth bins from the fill
    int ntruebins = 0;
    //std::cout << "EREC BIN: " << bin << std::endl;

    for (int ipar=0; ipar<parameters->NumberOfParameters(); ipar++) {
      if ( parameters->GetParameterFromIndex( ipar )->IsA()!="ResponseFunctionParameter" )
	continue;
      if ( parameters->GetParameterFromIndex( ipar )->GetFitterParType()!=kNewtonTerm ) 
	continue;
      
      // Get parameter sigma value
      double parvalue[1] = { parameters->GetParameterFromIndex( ipar )->GetSigma() };
      for (int j=0; j<true_energy_hists->GetNumberOfTruthHists(); j++) {

	// reset workspace array 
	// (reason we use a histogram object instead of simple array is that it allows us to pass the information to the oscillator )
	TH1D* outtruhist = (TH1D*)m_truthbins_fij->GetTruthHist( j );
	outtruhist->Reset();

	TH1D* intruhist = (TH1D*)true_energy_hists->GetTruthHist( j );

	// if no events in this component, just skip.
	if ( intruhist->Integral()==0 ) {
	  ntruebins += intruhist->GetNbinsX();
	  continue;
	}
	
	//std::cout << "in=" << intruhist->GetName() << "  out=" << outtruhist->GetName() << std::endl;

	for (int itrue=0; itrue<intruhist->GetNbinsX(); itrue++) {
	  double ntruevents = intruhist->GetBinContent( itrue+1 );
	  ntruebins++;
	  if ( ntruevents<=0 ) continue;

	  // CHECK THIS!
	  int ntruthbin_response = GetResponseTruthBin( ntruebins-1 );  // we allow the ability to have the binning for the response function bins and oscillation bins to be different.
	  //double parresponse = ((ResponseFunctionParameter*)parameters->GetParameterFromIndex( ipar ))->GetResponse( Sample::GetName(), bin, j*intruhist->GetNbinsX()+itrue, parvalue[0] );
	  double parresponse = ((ResponseFunctionParameter*)parameters->GetParameterFromIndex( ipar ))->GetResponse( m_sampleid_table[ipar], bin, j*intruhist->GetNbinsX()+itrue, parvalue[0] );
	  //double parresponse = ((ResponseFunctionParameter*)parameters->GetParameterFromIndex( ipar ))->GetResponse( m_sampleid_table[ipar], bin, ntruthbin_response, parvalue[0] );
	  
	  if ( parresponse==parresponse ) {
	    if ( parresponse>0 )
	      outtruhist->SetBinContent( itrue+1, ntruevents*(parresponse-1.0) );
	    else
	      outtruhist->SetBinContent( itrue+1, -ntruevents );
	  }
	  
// 	  if ( itrue==36 ) {
// 	    std::cout << "F[" << bin << ", " << j << ", " << itrue << "] = " << ntruevents*(parresponse-1.0)/parvalue[0] << std::endl;
// 	  }
	}//end of loop over true energy bins
	//       std::cout << "truthhist[" << j << "] " << intruhist->GetName() 
	// 		<< " in=" << intruhist->Integral()*6.3933e20 << " out=" << outtruhist->Integral()*6.3933e20 << std::endl;
      }//end of loop over true hists
    
      // now that truth histograms are now prepared, we can oscillate them and include the scaling factor
      double dNbin_per_sigma = oscillator->CalcEventsOsc( m_truthbins_fij, false )/parvalue[0]*SampleROOTPDF::GetScalingFactor(); // change of events per sigma
      SetFij( bin, ipar, dNbin_per_sigma );
    } //end of loop over parameters
    
  }//end of loop over rec bins
  
}

void SampleOscfitROOTPDF::SetFij( int bin, int ipar, double fijvalue ) {
  *(Fij + bin*m_fij_npars + ipar ) = fijvalue;
}

double SampleOscfitROOTPDF::GetFij( int bin, int ipar ) {
  return *(Fij + bin*m_fij_npars + ipar );
}

void SampleOscfitROOTPDF::GetBinFijs( int bin, double* fijvalues ) {
  memcpy( fijvalues, Fij+bin*m_fij_npars, sizeof(double)*m_fij_npars );
}

bool SampleOscfitROOTPDF::LoadSampleFromFile( TFile* file, double scalefactor ) {
  if ( GetVerbose()>=1 ) std::cout << "SampleOscfitROOTPDF[" << GetName() << "]::LoadSampleFromFile called." << std::endl;
  if ( GetSampleHistogram() ) {

    // Load sample histogream
    if ( GetVerbose()>=2 )
      GetSampleHistogram()->SetVerbose( 1 );

    bool loadedok = GetSampleHistogram()->LoadHistogramsFromFile( file, scalefactor );
        
    if ( !loadedok ) return false;
    return true;
  }
  return false;
}

void SampleOscfitROOTPDF::WriteToFile() {
  GetSampleHistogram()->WriteBinInfo();
}

int SampleOscfitROOTPDF::GetNumberOfTruthBinsPerTruthHist() {
  return m_truthbins->GetNumberOfTruthBinsPerHist();
}

int SampleOscfitROOTPDF::GetNumberOfTotalTruthBins() {
  return m_truthbins->GetTotalTruthBins();
}

int SampleOscfitROOTPDF::GetNumberOfTruthHists() {
  return m_truthbins->GetNumberOfTruthHists();
}

void SampleOscfitROOTPDF::GetListOfTruthHistNames( std::vector< std::string >& names ) {
  m_truthbins->GetListOfCutNames( names );
}

void SampleOscfitROOTPDF::GetListOfTruthHistDefinitions( std::vector< std::string >& defs ) {
  m_truthbins->GetListOfCutDefinitions( defs );
}

double SampleOscfitROOTPDF::GetOscTruthBinExpectation( int irec, int itruth_global ) {
  return GetTruthInfoBinExpectation( m_osc_bin_info_name, irec, itruth_global );
}

std::string SampleOscfitROOTPDF::GetTruthBinHistName( int irec, int itruth_global ) {
  TruthBinInfo* true_energy_hists = dynamic_cast<TruthBinInfo*>( GetSampleHistogram()->GetBinInfo( m_osc_bin_info_name, irec ) ); // truth bins from the fill
  return true_energy_hists->GetCutNameFromTruthBin( itruth_global );
}

// ------------------------------------------------------------------------------------------
// INTERFACE TO GENERIC TRUTH INFO BINS

int SampleOscfitROOTPDF::GetNumberOfTruthBinsPerTruthHistForBinInfo( std::string bin_info_name ) {
  TruthBinInfo* truthinfo = dynamic_cast<TruthBinInfo*>( GetSampleHistogram()->GetMasterListBinInfo( bin_info_name ) );
  return truthinfo->GetNumberOfTruthBinsPerHist();
}

int SampleOscfitROOTPDF::GetNumberOfTotalTruthInfoBins( std::string bin_info_name ) {
  TruthBinInfo* truthinfo = dynamic_cast<TruthBinInfo*>( GetSampleHistogram()->GetMasterListBinInfo( bin_info_name ) );
  return truthinfo->GetTotalTruthBins();
}

int SampleOscfitROOTPDF::GetNumberOfTruthInfoHists( std::string bin_info_name ) {
  TruthBinInfo* truthinfo = dynamic_cast<TruthBinInfo*>( GetSampleHistogram()->GetMasterListBinInfo( bin_info_name ) );
  return truthinfo->GetNumberOfTruthHists();
}

double SampleOscfitROOTPDF::GetTruthInfoBinExpectation( std::string bin_info_name, int irec, int itruth_global ) {
  TruthBinInfo* true_energy_hists = dynamic_cast<TruthBinInfo*>( GetSampleHistogram()->GetBinInfo( bin_info_name, irec ) ); // truth bins from the fill                                                                                   
  return true_energy_hists->GetTruthBinContentFromGlobalIndex( itruth_global );
}

void SampleOscfitROOTPDF::GetListOfTruthInfoHistNames( std::string bin_info_name, std::vector< std::string >& names ) {
  TruthBinInfo* true_energy_hists = dynamic_cast<TruthBinInfo*>( GetSampleHistogram()->GetMasterListBinInfo( bin_info_name ) ); // truth bins from the fill                                                                                   
  return true_energy_hists->GetListOfCutNames( names );
}

std::string SampleOscfitROOTPDF::GetTruthBinInfoHistName( std::string bin_info_name, int irec, int itruth_global ) {
  TruthBinInfo* true_energy_hists = dynamic_cast<TruthBinInfo*>( GetSampleHistogram()->GetBinInfo( bin_info_name, irec ) ); // truth bins from the fill
  return true_energy_hists->GetCutNameFromTruthBin( itruth_global );  
}

// ------------------------------------------------------------------------------------------

void SampleOscfitROOTPDF::UpdateOscPars( ParameterManager* parameters ) {

  if ( fUseEventReweightOsc ) {
    // Will change this interface later to be less dependent on type of oscillation calculator
    WeightBargerOsc* eventosc = (WeightBargerOsc*)event_oscillator;
    // EVENT-BY-EVENT OSC WEIGHTING TECHNIQUE (MOST PRECISE)
    if ( parameters->IsParameterDefined( "sin2_2theta23" ) )  {
      eventosc->SetSinThetaForm( kSin2Theta );
      eventosc->SetMNS( parameters->GetParameter("sin2_2theta12")->GetUnscaledValue(),
			parameters->GetParameter("sin2_2theta13")->GetUnscaledValue(),
			parameters->GetParameter("sin2_2theta23")->GetUnscaledValue(),
			parameters->GetParameter("dm2_12")->GetUnscaledValue(),
			parameters->GetParameter("dm2_32")->GetUnscaledValue(),
			parameters->GetParameter("delta_CP")->GetUnscaledValue() );
    }
    else if ( parameters->IsParameterDefined( "sin2_theta23" ) )  {
      eventosc->SetSinThetaForm( kSinTheta );
      eventosc->SetMNS( parameters->GetParameter("sin2_theta12")->GetUnscaledValue(),
			parameters->GetParameter("sin2_theta13")->GetUnscaledValue(),
			parameters->GetParameter("sin2_theta23")->GetUnscaledValue(),
			parameters->GetParameter("dm2_12")->GetUnscaledValue(),
			parameters->GetParameter("dm2_32")->GetUnscaledValue(),
			parameters->GetParameter("delta_CP")->GetUnscaledValue() );
    }
    eventosc->PrintMNS();
  }
  else {
    // HISTOGRAM OSC WEIGHTING TECHNIQUE (FASTER)
    if ( parameters->IsParameterDefined( "sin2_2theta23" ) )  {
      oscillator->SetSinThetaForm( kSin2Theta );
      oscillator->SetMNS( parameters->GetParameter("sin2_2theta12")->GetUnscaledValue(),
			  parameters->GetParameter("sin2_2theta13")->GetUnscaledValue(),
			  parameters->GetParameter("sin2_2theta23")->GetUnscaledValue(),
			parameters->GetParameter("dm2_12")->GetUnscaledValue(),
			  parameters->GetParameter("dm2_32")->GetUnscaledValue(),
			  parameters->GetParameter("delta_CP")->GetUnscaledValue() );
    }
    else if ( parameters->IsParameterDefined( "sin2_theta23" ) )  {
    oscillator->SetSinThetaForm( kSinTheta );
    oscillator->SetMNS( parameters->GetParameter("sin2_theta12")->GetUnscaledValue(),
			parameters->GetParameter("sin2_theta13")->GetUnscaledValue(),
			parameters->GetParameter("sin2_theta23")->GetUnscaledValue(),
			parameters->GetParameter("dm2_12")->GetUnscaledValue(),
			parameters->GetParameter("dm2_32")->GetUnscaledValue(),
			parameters->GetParameter("delta_CP")->GetUnscaledValue() );
    }
    oscillator->PrintMNS();
  }
}

void SampleOscfitROOTPDF::UseOscEventReweight(  SinThetaForms sinthetaform, std::string nuE_GeV_var, std::string nuflux_var, std::string nuxsec_var, std::string mode_var, std::string weight_var ) {
  if ( event_oscillator ) {
    std::cout << "Tried to redfine WeightBargerOsc object!" << std::endl;
    assert(false);
  }
  HistRootVariable* pdfrootvar = (HistRootVariable*)GetSampleHistogram();
  event_oscillator = new WeightBargerOsc( sinthetaform, pdfrootvar->GetSourceChain(), nuE_GeV_var, nuflux_var, nuxsec_var, mode_var, weight_var );
  pdfrootvar->LoadWeightClass( "bargerosc_"+GetName(), event_oscillator );
  fUseEventReweightOsc = true;
  std::cout << "Sample, " << GetName() << ", will now use event-by-event oscillations." << std::endl;
}

void SampleOscfitROOTPDF::SetTemplateMode( bool ison ) {
  HistRootVariable* pdfrootvar = (HistRootVariable*)GetSampleHistogram();
  WeightBargerOsc* weightgen = (WeightBargerOsc*)pdfrootvar->GetWeightClass( "bargerosc_"+GetName() );
  if ( weightgen )
    weightgen->SetMapMode(ison);
}

void SampleOscfitROOTPDF::StoreResponseParSampleID( ParameterManager* parameters ) {
  // This is for an optimzation trick.
  // This is to avoid a string comparison and map lookup
  m_sampleid_table = new int[parameters->NumberOfParameters()];
  for (int ipar=0; ipar<parameters->NumberOfParameters(); ipar++) {
    if ( parameters->GetParameterFromIndex( ipar )->IsA()=="ResponseFunctionParameter" ) {
      m_sampleid_table[ipar] = ((ResponseFunctionParameter*)parameters->GetParameterFromIndex( ipar ))->GetSampleID( GetName() );
    }
  }
}

void SampleOscfitROOTPDF::SetSaveTruthModeHistsFlag( bool savehists ) {
  fSaveTruthModeHists = savehists;
  if ( fSaveTruthModeHists && truthmode_histograms==NULL ) {
    // get list of cuts
    Hist* hist = GetSampleHistogram();
    truthmode_histograms = new Hist*[m_truthbins->m_ncutlists];
    for (int i=0; i<m_truthbins->m_ncutlists; i++) {
      std::string cutname = m_truthbins->m_cut_list[i].name;
      std::string histname = GetName()+"_modes_"+cutname;
      truthmode_histograms[i] = new Hist( histname );
      if (GetSampleHistogram()->GetHistDimensions()==1)
	truthmode_histograms[i]->SetHistogram( (TH1D*)hist->GetHistogram()->Clone( histname.c_str() ) );
      else
	truthmode_histograms[i]->SetHistogram( (TH2D*)hist->GetHistogram()->Clone( histname.c_str() ) );
    }
  }
}

void SampleOscfitROOTPDF::WriteTruthModeHists() {
  if ( !fSaveTruthModeHists ) return;
  m_truthbins->Print();
  for (int i=0; i<m_truthbins->m_ncutlists; i++) {
    truthmode_histograms[i]->GetHistogram()->Write();
  }
}

int SampleOscfitROOTPDF::GetResponseTruthBin( int osc_truth_bin ) {
  if (!fTruthToResponseMapMade ) return osc_truth_bin;
  return m_sample_truth_to_splinebin[osc_truth_bin];
}

// void SampleOscfitROOTPDF::MakeTruthBinToSplineBinMap() {

//   // check that a separate response binning was defined.
//   fTruthToResponseMapMade = false;
//   if ( GetResponseBinInfoName()=="" ) {
//     fUseDifferentOscAndResponseBinning = false;
//     return;
//   }
//   if ( GetResponseBinInfoName()==GetOscBinInfoName() ) {
//     fUseDifferentOscAndResponseBinning = false;
//     return;
//   }
//   int ndims = GetSampleHistogram()->GetHistDimensions();

//   TruthBinInfo* truth_bins = dynamic_cast< TruthBinInfo* >( GetSampleHistogram()->GetMasterListBinInfo( GetOscBinInfoName() ) );
//   TruthBinInfo* spline_bins = dynamic_cast< TruthBinInfo* >( GetSampleHistogram()->GetMasterListBinInfo( GetResponseBinInfoName() ) );
//   if ( spline_bins==NULL ) {
//     std::cout << __PRETTY_FUNCTION__ << ": Error. Did not find spline user bin info object in Master BinInfo list." << std::endl;
//     assert(false);
//   }
//   if ( truth_bins==NULL ) {
//     std::cout << __PRETTY_FUNCTION__ << ": Error. Did not find truth(osc.) user bin info object in Master BinInfo list." << std::endl;
//     assert(false);
//   }
  
//   if ( spline_bins->GetNumberOfTruthHists()==truth_bins->GetNumberOfTruthHists() ) {
//     std::cout << "Spline and Truth bins should have same number of categories/cuts." << std::endl;
//     std::cout << "Spline bins: " << spline_bins->GetNumberOfTruthHists() << std::endl;
//     std::cout << "Truth (osc.) bins: " << truth_bins->GetNumberOfTruthHists() << std::endl;
//     assert(false);
//   }
//
//   if ( !m_sample_truth_to_splinebin ) {
//     // allocate the dictionary array
//     m_sample_truth_to_splinebin = new int[truth_bins->GetTotalTruthBins()];
//     memset( m_sample_truth_to_splinebin, 0, sizeof(int)*truth_bins->GetTotalTruthBins() );
//   }

//   int itruthbin=0;
//   for (int ihist=0; ihist<truth_bins->GetNumberOfTruthHists(); ihist++) {
//     TH1* truthhist = truth_bins->GetTruthHist( ihist );
//     TH1* splinehist = spline_bins->GetTruthHist( ihist );

//     for (int ibin=0; ibin<truthhist->GetNbinsX(); ibin++) {
// 	itruthbin++;
// 	if ( ndims==1 ) {
// 	  m_sample_truth_to_splinebin[itruthbin-1] = ((TH1D*)splinehist)->GetXaxis()->FindBin( ((TH1D*)truthhist)->GetXaxis()->GetBinCenter( itruthbin ) );
// 	}
// 	else {
// 	  std::cout << "2D truth binning not yet supported" << std::endl;
// 	  assert(false);
// 	}
//     }//loop over truth hists
//   }
//   fTruthToResponseMapMade = true;
// }
