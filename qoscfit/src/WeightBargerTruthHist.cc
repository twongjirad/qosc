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

#include "WeightBargerTruthHist.hh"
#include <assert.h>

#include "TChain.h"
#include "TMath.h"
#include "TH1D.h"

#include "TruthBinInfo.hh"
#include "RootVariableList.hh"

#include "BargerPropagator.h"

using namespace qosc;

WeightBargerTruthHist::WeightBargerTruthHist( TruthBinInfo* truth_bin_def ) { 
  /** 
      This takes a source chain file (and its friend trees) along with the leaf names 
      for the neutrino energy in GeV (nuE_GeV_var), the neutrino flux flavor (nuflux_var)
      and the neutrino cross section flavor (nuxsec_var).
      It then utilizes a RootVariableList to create a hook into the ROOT chain in order 
      to access the value of these variables event by event.
  */

  m_BP = new BargerPropagator( false ); // This is Roger's oscillation calculator
  fParamsSet = false; // Parameters have not been set yet
  
  SetSinThetaForm( kSinTheta );
  fMapMode = false;
  fUseCache = true;
  fCacheMade = false;
  fCalcGradient = false;
  fGradCacheMade = false;
  SetVerbose(0);

  Configure( truth_bin_def );

}


WeightBargerTruthHist::~WeightBargerTruthHist() {
  delete m_BP;
  delete [] EnergyBins;
  delete [] component_index;
  for (int n=0; n<18; n++) {
    if ( OscProbs[n] ) delete [] OscProbs[n];
    OscProbs[n] = NULL;
  }
}

void WeightBargerTruthHist::Configure( TruthBinInfo* truth_bin_def ) {
  // Things that need to be done
  // (1) Store energy of bins, don't want to have look this up constantly
  // (2) For each component histogram in the truthbininfo, we need to store what flux and interaction flavor is tied to it
  // (3) Build the probability cache
  BuildEnergyArray( truth_bin_def );
  BuildComponentIndex( truth_bin_def );
  BuildProbCache();
}

void WeightBargerTruthHist::BuildEnergyArray( TruthBinInfo* truth_bin_def  ) {
  // get truth bins
  TH1D* htemplate = dynamic_cast<TH1D*>( truth_bin_def->GetTemplateHist() );
  if ( !htemplate ) {
    assert(false);
  }
  nEnergyBins = htemplate->GetNbinsX()+2; // 2 for overflow and underflow
  EnergyBins = new TAxis( *htemplate->GetXaxis() );
}

void WeightBargerTruthHist::BuildComponentIndex( TruthBinInfo* truth_bin_def  ) {
  // get truth bins
  //std::cout << "WeightBargerTruthHist::BuildComponentIndex(" << truth_bin_def << ") ncomps = " << truth_bin_def->GetNumberOfTruthHists() << std::endl;
  component_index = new component_def[truth_bin_def->GetNumberOfTruthHists()];
  for ( unsigned int n=0; n<truth_bin_def->GetNumberOfTruthHists(); n++ ) {
    std::string name = truth_bin_def->m_cut_list[n].name;

    // set flux
    bool foundbar = false;
    if ( name.find("bar")!=std::string::npos ) foundbar = true;
    if ( name.find("numu")!=std::string::npos ) {
      component_index[n].flux = ( foundbar ) ? -2 : 2;
      component_index[n].xsec = ( foundbar ) ? -2 : 2;
    }
    else if ( name.find("nue")!=std::string::npos )  {
      component_index[n].flux = ( foundbar ) ? -1 : 1;
      component_index[n].xsec = ( foundbar ) ? -1 : 1;
    }
    else if ( name.find("nutau")!=std::string::npos )  {
      component_index[n].flux = ( foundbar ) ? -3 : 3;
      component_index[n].xsec = ( foundbar ) ? -3 : 3;
    }
    else if ( name.find("signal")!=std::string::npos ) {
      component_index[n].flux = ( foundbar ) ? -2 : 2;
      component_index[n].xsec = ( foundbar ) ? -1 : 1;
    }
    else {
      std::cout << "comp name not parsed correctly: " << name << std::endl;
      assert(false);
    }

    // set mode (CC or NC)
    if ( name.find("cc")!=std::string::npos ) {
      component_index[n].mode = 1;
    }
    else {
      component_index[n].mode = 33;
    }
  }
}

void WeightBargerTruthHist::BuildProbCache() {
  for (int n=0; n<18; n++) {
    OscProbs[n] = new double[nEnergyBins];
    memset( OscProbs[n], 0, sizeof(double)*nEnergyBins );
    for (int i=0; i<6; i++) {
      GradOscProbs[n][i] = new double[nEnergyBins];
      memset( GradOscProbs[n][i], 0, sizeof(double)*nEnergyBins );
    }
  }
}


void WeightBargerTruthHist::SetMNS( double s12, double s13, double s23, double dm12, double dm23, double CP, SinThetaForms sin_term ) {
  bool newvalues = false;
  if ( m_s12!=s12 ) {
    m_s12 = s12;
    newvalues = true;
  }
  if ( m_s13!=s13 ) {
    m_s13=s13;
    newvalues = true;
  }
  if ( m_s23!=s23 ) {
    m_s23 = s23;
    newvalues = true;
  }
  if ( m_dm12!=dm12 ) {
    m_dm12 = dm12;
    newvalues = true;
  }
  if ( m_dm23!=dm23 ) {
    m_dm23 = dm23;
    newvalues = true;
  }
  if ( m_CP!=CP ) {
    m_CP = CP;
    newvalues = true;
  }
  if ( sin_term!=fSinThetaForm ) {
    SetSinThetaForm( sin_term );
    newvalues = true;
  }

  if ( newvalues ) {

    if ( m_s23>1.0 ) m_s23 = 1.0;
    else if ( m_s23<0.0 ) m_s23 = 0.;
    
    if ( m_s13>1.0 ) m_s13 = 1.0;
    else if ( m_s13<0.0 ) m_s13 = 0.;
    
    if ( m_s12>1.0 ) m_s12 = 1.0;
    else if ( m_s12<0.0 ) m_s12 = 0.;
        
    fCacheMade = false;
    fGradCacheMade = false;
  }
  fParamsSet = true;
}

void WeightBargerTruthHist::SetMNS( double s12, double s13, double s23, double dm12, double dm23, double CP ) {
    SetMNS( s12, s13, s23, dm12, dm23, CP, fSinThetaForm);
}

void WeightBargerTruthHist::SetParameter( std::string parname, double value ) {
  if ( parname=="sin2_2theta12" ) m_s12 = value;
  else if ( parname=="sin2_2theta32" ) m_s23 = value;
  else if ( parname=="sin2_2theta13" ) m_s13 = value;
  else if ( parname=="dm2_12" ) m_dm12 = value;
  else if ( parname=="dm2_32" ) m_dm23 = value;
  else if ( parname=="delta_CP" ) m_CP = value;
  else assert(false);
  fCacheMade = false;
  fGradCacheMade = false;

  if ( m_s23>1.0 ) m_s23 = 1.0;
  else if ( m_s23<0.0 ) m_s23 = 0.;

  if ( m_s13>1.0 ) m_s13 = 1.0;
  else if ( m_s13<0.0 ) m_s13 = 0.;

  if ( m_s12>1.0 ) m_s12 = 1.0;
  else if ( m_s12<0.0 ) m_s12 = 0.;

}

void WeightBargerTruthHist::PrintMNS() {
  std::cout << "------------------------" << std::endl;
  std::cout << "WeightBargerTruthHist PMNS Matrix: " << std::endl;
  if ( fSinThetaForm==kSin2Theta ) {
    std::cout << "  sin2(2theta12) =" << m_s12 << std::endl;
    std::cout << "  sin2(2theta23) =" << m_s23 << std::endl;
    std::cout << "  sin2(2theta13) =" << m_s13 << std::endl;
  }
  else {
    std::cout << "  sin2(theta12) =" << m_s12 << std::endl;
    std::cout << "  sin2(theta23) =" << m_s23 << std::endl;
    std::cout << "  sin2(theta13) =" << m_s13 << std::endl;
  }
  std::cout << "  (dm12)^2 = " << m_dm12 << " eV^2" << std::endl;
  std::cout << "  (dm23)^2 = " << m_dm23 << " eV^2" << std::endl;
  std::cout << "  delta_CP = " << m_CP*180/3.14159 << " deg" << std::endl;
  std::cout << " Normal Hierarchy" << std::endl;
  std::cout << " UseCache: " << fUseCache << std::endl;
  std::cout << " CacheMade: " << fCacheMade << std::endl;
  std::cout << " MapMode: " << fMapMode << std::endl;
  std::cout << " SinThetaTerm: " << fSinThetaForm  << std::endl;
  std::cout << "------------------------" << std::endl;
}

void WeightBargerTruthHist::SetT2K( int flavortype ) {
  double baseline = 295.0; // km
  double density = 2.6; // g/cm3
  m_BP->propagateLinear( flavortype, baseline, density ); // probably initiates the oscillation machinery with the parameters for dealing with oscillations through uniform, dense matter
}

void WeightBargerTruthHist::SetNEUTmode( int mode ) {
  m_imode = mode;
}

void WeightBargerTruthHist::SetIncomingNeutrino( std::string flavor ) {
  // obsolete now?
  // neutrino types
  // 1:e 2:mu 3:tau   -1: e_bar -2: mu_bar -3: tau_bar
  if (flavor=="numu") m_iflux = 2;
  else if (flavor=="nue") m_iflux = 1;
  else if (flavor=="nutau") m_iflux = 3;
  else if (flavor=="numubar" || flavor=="numu_bar") m_iflux=-2;
  else if (flavor=="nuebar" || flavor=="nue_bar") m_iflux=-1;
  else if (flavor=="nutaubar" || flavor=="nutau_bar") m_iflux=-3;
  else {
    std::cout << "Unrecognized flavor, " << flavor << ", given to WeightBargerTruthHist::SetIncomingNeutrino "<< std::endl;
    assert(false);
  }
}

void WeightBargerTruthHist::SetOutgoingNeutrino( std::string flavor ) {
  // obsolete now?
  // neutrino types
  // 1:e 2:mu 3:tau   -1: e_bar -2: mu_bar -3: tau_bar
  if (flavor=="numu") m_ixsec = 2;
  else if (flavor=="nue") m_ixsec = 1;
  else if (flavor=="nutau") m_ixsec = 3;
  else if (flavor=="numubar" || flavor=="numu_bar") m_ixsec=-2;
  else if (flavor=="nuebar" || flavor=="nue_bar") m_ixsec=-1;
  else if (flavor=="nutaubar" || flavor=="nutau_bar") m_ixsec=-3;
  else {
    std::cout << "Unrecognized flavor, " << flavor << ", given to WeightBargerTruthHist::SetOutgoingNeutrino "<< std::endl;
    assert(false);
  }
}


void WeightBargerTruthHist::SetOscillationComp( int comp ) {
  // obsolete now?
  // for CC, first set incoming and outgoing neutrino flavors
  switch (comp) {
  case 0:
    SetIncomingNeutrino( "numu" ); 
    SetOutgoingNeutrino( "numu" );
    SetNEUTmode( 1 );
    break;
  case 1:
    SetIncomingNeutrino( "nue" );
    SetOutgoingNeutrino( "nue" );
    SetNEUTmode( 1 );
    break;
  case 2:
    SetIncomingNeutrino( "numu_bar" ); 
    SetOutgoingNeutrino( "numu_bar" );
    SetNEUTmode( 1 );
    break;
  case 3:
    SetIncomingNeutrino( "numu" ); 
    SetOutgoingNeutrino( "nue" );
    SetNEUTmode( 1 );
    break;
  case 4:
    SetIncomingNeutrino( "numu" ); 
    SetOutgoingNeutrino( "numu" );
    SetNEUTmode( 40 );
    break;
  };
  
}


double WeightBargerTruthHist::CalculateWeight( double energy_GeV ) {

  if (fParamsSet==false) {
    std::cout << "The oscillation parameters for the BargerPropagator instance has not been set. Cannot calculate weight." << std::endl;
    assert(false);
  }

  // Get Mode (for NC)
  if ( abs(m_imode)>30 ) {
    // no flavor-changing NC
    if (m_iflux!=m_ixsec ) 
      //weight = 0.0; //for debug
      return 0.0;
    else
      return 1.0;
  }


  //return 1.0; // HACK
  
  if ( fUseCache ) {
    return CalculateWeightFromCache( energy_GeV );
  }
  else {
    return CalculateWeightFromBarger( energy_GeV );
  }

  assert(false);
}

double WeightBargerTruthHist::CalculateWeightFromCache( double energy_GeV ) {
  int index=0, bin=0;
  if ( !fCacheMade ) MakeCache();
  // Use pre-calculated oscillation probabilities.
  index = GetProbArrayIndex( m_iflux, m_ixsec );
  bin = EnergyBins->FindBin( energy_GeV );
  return OscProbs[index][bin];
}

double WeightBargerTruthHist::CalculateWeightFromBarger( double energy_GeV ) {
  
  if (fMapMode==true)
    return 1.0; // needed when making signal pdfs
    
  // Calculate Oscillation Probability on-demand.

  // SetMNS
  bool kSquared = false; // { true: sin^2( theta ), false: sin^2( 2 theta ) }
  if ( fSinThetaForm==kSin2Theta ) kSquared = false;
  else kSquared = true;

  m_BP->SetMNS( m_s12, m_s13, m_s23, m_dm12, m_dm23, m_CP, energy_GeV,  kSquared, m_iflux );
  //m_BP->SetMNS( m_s12, m_s13, m_s23, m_dm12, m_dm23, m_CP, energy_GeV,  kSquared );
  SetT2K(m_iflux);

  // neutrino types
  // 1:e 2:mu 3:tau   -1: e_bar -2: mu_bar -3: tau_bar
  double prob = m_BP->GetProb( m_iflux, m_ixsec );
  return prob;
}

int WeightBargerTruthHist::GetProbArrayIndex( int incoming, int outgoing ) {
  if ( incoming*outgoing<0 ) return -1; /// neutrino and anti-neutrino combination. Invalid.

  int sign = 1;
  if ( incoming<0 ) sign = -1;
  
  return 3*(abs(incoming)-1) + (abs(outgoing)-1) + 9*((1-sign)/2);
}

void WeightBargerTruthHist::MakeCache() {
  if ( nEnergyBins<=0 ) return;
  if ( !fParamsSet ) return;

  if ( GetVerbose() ) {
    std::cout << "Making Weight Cache (" << nEnergyBins << " bins)" << std::endl;
    PrintMNS();
  }

  

  SetNEUTmode( 1 );

  for (int in=1; in<=3; in++) {
    for (int out=1; out<=3; out++) {

      for (int neu=1; neu>=-1; neu+=-2) {
	int index = GetProbArrayIndex( neu*in, neu*out );
	if ( index<0 ) {
	  continue;
	}

	//SetT2K( neu*in ); // Set Flux
	m_iflux = neu*in; // set flux
	m_ixsec = neu*out; // set cross section
	
	for (int bin=0; bin<nEnergyBins; bin++) {
	  OscProbs[index][bin] = CalculateWeightFromBarger( EnergyBins->GetBinCenter(bin) );
	}
      }//en fo loop over neutrino or anti-neutrino
    }//end of loop over outgoing flavor
  }//end of loop over incoming flavor
  fCacheMade = true;
}

void WeightBargerTruthHist::MakeGradCache() {

  if ( nEnergyBins<=0 ) return;
  if ( !fParamsSet ) return;

  if ( GetVerbose() ) {
    std::cout << "Making Gradient Probability Cache (" << nEnergyBins << " bins)" << std::endl;
    PrintMNS();
  }

  SetNEUTmode(1);
  
  for (int in=1; in<=3; in++) {
    for (int out=1; out<=3; out++) {
      
      for (int neu=1; neu>=-1; neu+=-2) {
	int index = GetProbArrayIndex( neu*in, neu*out );
	
	SetT2K( neu*in ); // Set Flux
	m_iflux = neu*in; // set flux
	m_ixsec = neu*out; // set cross section
	
	for (int par=0; par<kNumPars; par++) {  
	  if ( !GradOscProbs[index][par] ) GradOscProbs[index][par] = new double[nEnergyBins];
	
	  for (int bin=0; bin<nEnergyBins; bin++) {
	    
	    double energy_GeV = EnergyBins->GetBinCenter(bin);
	    	    
	    GradOscProbs[index][par][bin] = CalculateDerivativeFromBarger( par, energy_GeV );
	  }
	  
	}//end of loop over par
      }//en fo loop over neutrino or anti-neutrino
    }//end of loop over outgoing flavor
  }//end of loop over incoming flavor
  
  fGradCacheMade = true;

}

double WeightBargerTruthHist::CalculateDerivative( int par, double energy_GeV ) {

  if (fParamsSet==false) {
    std::cout << "The oscillation parameters for the BargerPropagator instance has not been set. Cannot calculate weight." << std::endl;
    assert(false);
  }

  // Derivative of constant is zero
  if ( abs(m_imode)>30 ) {
    return 0;
  }

  double dir = 0;
  if ( fUseCache ) {
    dir= CalculateDerivativeFromCache( par, energy_GeV );
  }
  else {
    dir= CalculateDerivativeFromBarger( par, energy_GeV );
  }
  std::cout << "oscP derivative for " << par << ": " << dir <<  std::endl;
  return dir;
  assert(false);
}

double WeightBargerTruthHist::CalculateDerivativeFromBarger( int par, double energy_GeV ) {

  // by assumption, flux, xsec and mode are implicitly setup
  // this just accesses the barger oscillator class twice to estimate the derivative

  // create variation in parameter
  bool kSquared = false;
  double dparup[6] = { m_s12, m_s13, m_s23, m_dm12, m_CP };
  double dpardown[6] = { m_s12, m_s13, m_s23, m_dm12, m_CP };
  double epsilon = dparup[par]*1.0e-6;
  if ( epsilon==0 )
    epsilon = 1.0e-6;
  dparup[par] += fabs(epsilon);
  dpardown[par] -= fabs(epsilon);

  if ( dparup[ks12]>1 ) { return 0; }//dparup[ks12] = 1.0; epsilon*=0.5; }
  if ( dparup[ks13]>1 ) { return 0; }//dparup[ks13] = 1.0; epsilon*=0.5; }
  if ( dparup[ks32]>1 ) { return 0; }//dparup[ks32] = 1.0; epsilon*=0.5; }
  if ( dpardown[ks12]<0 ) { return 0; }//dpardown[ks12] = 0.0; epsilon*=0.5; }
  if ( dpardown[ks13]<0 ) { return 0; }//dpardown[ks13] = 0.0; epsilon*=0.5; }
  if ( dpardown[ks32]<0 ) { return 0; }//;dpardown[ks32] = 0.0; epsilon*=0.5; }

  if ( fSinThetaForm==kSin2Theta )
    kSquared = false;
  else
    kSquared = true;

  m_BP->SetMNS( dparup[ks12], dparup[ks13], dparup[ks32], dparup[kdm12], dparup[kdm32], dparup[kcp], energy_GeV,  kSquared, m_iflux );
  //m_BP->SetMNS( dparup[ks12], dparup[ks13], dparup[ks32], dparup[kdm12], dparup[kdm32], dparup[kcp], energy_GeV,  kSquared );
  
  SetT2K( m_iflux );
  double P1 = m_BP->GetProb( m_iflux, m_ixsec );
  
  m_BP->SetMNS( dpardown[ks12], dpardown[ks13], dpardown[ks32], dpardown[kdm12], dpardown[kdm32], dpardown[kcp], energy_GeV,  kSquared, m_iflux );
  //m_BP->SetMNS( dpardown[ks12], dpardown[ks13], dpardown[ks32], dpardown[kdm12], dpardown[kdm32], dpardown[kcp], energy_GeV,  kSquared );
		
  SetT2K( m_iflux );

  double P2 = m_BP->GetProb( m_iflux, m_ixsec );
//   if ( par==ks32 && m_iflux==2 && m_ixsec==2 ) {
//     // check out change in numu retension due to s32
//     std::cout << "barger derivative " 
// 	      << "d(" << dparup[par] << " - " << dpardown[par] << " @ " << energy_GeV << " GeV ): " 
// 	      << P1 << "-" << P2 << "/" << 2*epsilon << " = " << (P1-P2)/(2*fabs(epsilon)) << std::endl;
//     std::cin.get();
//   }
  return (P1-P2)/(2*epsilon);
  
}

double WeightBargerTruthHist::CalculateDerivativeFromCache( int par, double energy_GeV ) {
  int index=0, bin=0;
  if ( !fGradCacheMade ) MakeGradCache();
  // Use pre-calculated oscillation probabilities.
  index = GetProbArrayIndex( m_iflux, m_ixsec );
  bin = EnergyBins->FindBin( energy_GeV );   // if constant bins (can optiize here if know bin formula)
  return GradOscProbs[index][par][bin];
}

double WeightBargerTruthHist::CalcEventsOsc( TruthBinInfo* truth_info, bool modifyhist ) {
  int ncomps = truth_info->GetNumberOfTruthHists();
  if ( !fCacheMade ) MakeCache();

  double nevents = 0;
  for (int comp=0; comp<ncomps; comp++) {
    TH1D* comp_hist = (TH1D*)truth_info->GetTruthHist( comp );
    if ( comp_hist->Integral()==0 ) continue;


    if ( abs(component_index[comp].mode)>30 ) {
      // NC
      if ( component_index[comp].flux==component_index[comp].xsec ) {
	nevents += comp_hist->Integral();
      }
      else {
	if ( modifyhist )
	  comp_hist->Reset(); // zero out the histogram
      }
    }
    else {
      // CC
      int index = GetProbArrayIndex( component_index[comp].flux, component_index[comp].xsec );
      for (int bin=1; bin<=EnergyBins->GetNbins(); bin++) {
	nevents += comp_hist->GetBinContent(bin)*OscProbs[index][bin]; /// GOOD PLACE FOR GPU KERNAL
	if ( modifyhist )
	  comp_hist->SetBinContent( bin, comp_hist->GetBinContent(bin)*OscProbs[index][bin] );
      }
    }
    
//     std::cout << "[" << comp << "] " << comp_hist->GetName() << ": index=" << index 
// 	      << " flux=" << component_index[comp].flux << " xsec=" << component_index[comp].xsec << " mode=" << component_index[comp].mode 
// 	      << " integral=" << comp_hist->Integral() << " (post-osc, mod=" << modifyhist << ")"
// 	      << std::endl;
  }
  return nevents;
}
