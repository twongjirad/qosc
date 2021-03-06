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

#include "WeightBargerOsc.hh"
#include <assert.h>
#include <iostream>
#include "BargerPropagator.h"
#include "TChain.h"
#include "TMath.h"
#include "RootVariableList.hh"
#include "ParameterManager.hh"
#include "ModelParameter.hh"

using namespace qosc;

WeightBargerOsc::WeightBargerOsc( SinThetaForms sin_theta_type, TChain* source, std::string nuE_GeV_var, std::string nuflux_var, std::string nuxsec_var, std::string mode_var, std::string weight_var ) 
  : WeightOscEvent( source )
{
  /** 
      This takes a source chain file (and its friend trees) along with the leaf names 
      for the neutrino energy in GeV (nuE_GeV_var), the neutrino flux flavor (nuflux_var)
      and the neutrino cross section flavor (nuxsec_var).
      It then utilizes a RootVariableList to create a hook into the ROOT chain in order 
      to access the value of these variables event by event.
  */

  m_BP = new BargerPropagator( false ); // This is Roger's oscillation calculator
  SetParamsSetFlag( false ); // Parameters have not been set yet
  SetSourceChain( source ); // The source chain is expected to have information on the flux and xsec flavor of the event

  m_nuE_GeV_var = nuE_GeV_var; // The name of the neutrino energy variable (in GeV)
  m_nuflux_var = nuflux_var; // Name of the neutrino flux flavor variable
  m_nuxsec_var =  nuxsec_var; // Name of the neutrino xsec variable
  m_mode_var = mode_var; // Mode number of interaction
  m_weight_var = weight_var; // Normalization factor (extracted by separate code and included as a friend tree)
 
  // Create instances of the variables here
  m_rootvars.SetChain( source );  
  m_rootvars.Add( m_nuE_GeV_var );
  m_rootvars.Add( m_nuflux_var );
  m_rootvars.Add( m_nuxsec_var );
  m_rootvars.Add( m_mode_var );
  if ( m_weight_var!="" )
    m_rootvars.Add( m_weight_var );

  fSinThetaForm = sin_theta_type;
  SetMapMode(false);
  m_POT = 1.0;
  fParsAssigned = false;
}

WeightBargerOsc::WeightBargerOsc( SinThetaForms sin_theta_type, TChain* source, std::string nuE_GeV_var, std::string nuflux_var, std::string nuxsec_var, std::string mode_var, std::string weight_var,
				  std::string s12name, std::string s13name, std::string s23name, std::string dm12name, std::string dm32name, std::string cpdname, ParameterManager* parman )
  : WeightOscEvent( source )
{
  /** 
      This takes a source chain file (and its friend trees) along with the leaf names 
      for the neutrino energy in GeV (nuE_GeV_var), the neutrino flux flavor (nuflux_var)
      and the neutrino cross section flavor (nuxsec_var).
      It then utilizes a RootVariableList to create a hook into the ROOT chain in order 
      to access the value of these variables event by event.
  */

  m_BP = new BargerPropagator( false ); // This is Roger's oscillation calculator
  SetParamsSetFlag( false ); // Parameters have not been set yet
  SetSourceChain( source ); // The source chain is expected to have information on the flux and xsec flavor of the event

  m_nuE_GeV_var = nuE_GeV_var; // The name of the neutrino energy variable (in GeV)
  m_nuflux_var = nuflux_var; // Name of the neutrino flux flavor variable
  m_nuxsec_var =  nuxsec_var; // Name of the neutrino xsec variable
  m_mode_var = mode_var; // Mode number of interaction
  m_weight_var = weight_var; // Normalization factor (extracted by separate code and included as a friend tree)
 
  // Create instances of the variables here
  m_rootvars.SetChain( source );  
  m_rootvars.Add( m_nuE_GeV_var );
  m_rootvars.Add( m_nuflux_var );
  m_rootvars.Add( m_nuxsec_var );
  m_rootvars.Add( m_mode_var );
  if ( m_weight_var!="" )
    m_rootvars.Add( m_weight_var );

  fSinThetaForm = sin_theta_type;
  SetMapMode(false);
  m_POT = 1.0;
  fParsAssigned = false;

  // assign pars
  AssignPars( s12name, s13name, s23name, dm12name, dm32name, cpdname, parman );
}

WeightBargerOsc::~WeightBargerOsc() {
  delete m_BP;
}

void WeightBargerOsc::SetMNS( double s12, double s13, double s23, double dm12, double dm23, double CP ) {
  SetMNS( s12, s13, s23, dm12, dm23, CP, fSinThetaForm );
}

void WeightBargerOsc::SetMNS( double s12, double s13, double s23, double dm12, double dm23, double CP, SinThetaForms sin_term_type ) {
  m_osc_pars[ks12] = s12;
  m_osc_pars[ks13] = s13;
  m_osc_pars[ks23] = s23;
  m_osc_pars[kdm12] = dm12;
  m_osc_pars[kdm32] = dm23;
  m_osc_pars[kcp] = CP;
  fParamsSet = true;
  fSinThetaForm = sin_term_type;
}

void WeightBargerOsc::PrintMNS() {
  std::cout << "------------------------" << std::endl;
  std::cout << "Setting PMNS Matrix: " << std::endl;
  if ( fParsAssigned ) {
    for (int i=0; i<kNumPars; i++) {
      std::cout << "  " << parnames[i] << " = " << m_osc_pars[i] << std::endl;
    }
  }
  else {  
    std::cout << "  sin2(2theta12) =" << m_osc_pars[ks12] << std::endl;
    std::cout << "  sin2(2theta13) =" << m_osc_pars[ks13] << std::endl;
    std::cout << "  sin2(2theta23) =" << m_osc_pars[ks23] << std::endl;
    std::cout << "  (dm12)^2 = " << m_osc_pars[kdm12] << " eV^2" << std::endl;
    std::cout << "  (dm23)^2 = " << m_osc_pars[kdm32] << " eV^2" << std::endl;
    std::cout << "  delta_CP = " << m_osc_pars[kcp]*180/3.14159 << " deg" << std::endl;
  }
  //std::cout << " Normal Hierarchy" << std::endl;
  std::cout << "------------------------" << std::endl;
}

void WeightBargerOsc::SetT2K( int flavortype ) {
  double baseline = 295.0; // km
  double density = 2.6; // g/cm3
  m_BP->propagateLinear( flavortype, baseline, density ); // probably initiates the oscillation machinery with the parameters for dealing with oscillations through uniform, dense matter
}

double WeightBargerOsc::CalculateWeight() {

  if (fParamsSet==false) {
    std::cout << "The oscillation parameters for the BargerPropagator instance has not been set. Cannot calculate weight." << std::endl;
    assert(false);
  }

  // Get Normalization weight
  double weight = m_POT;
  if ( m_weight_var!="" )
    weight *= m_rootvars.GetVariableValue( m_weight_var );

  // Get Neutrino Flavors
  int iflux = int( m_rootvars.GetVariableValue( m_nuflux_var ) );
  int ixsec = int( m_rootvars.GetVariableValue( m_nuxsec_var ) );

  // Get Mode (if NC or CC)
  int mode = int( m_rootvars.GetVariableValue( m_mode_var ) );
  if ( abs(mode)>30 ) {
    // no flavor-changing NC
    if (iflux!=ixsec ) 
      //weight = 0.0; //for debug
      return 0.0;
    else
      return weight;
  }

  // CC
  if (fMapMode==true) {
    return weight; // needed when making signal pdfs
  }

  
  // Get Energy
  double energy_GeV = m_rootvars.GetVariableValue(m_nuE_GeV_var);

  // SetMNS
  bool kSquared = false; // { true: sin^2( theta ), false: sin^2( 2 theta ) }
  if ( fSinThetaForm==kSin2Theta ) {
    kSquared = false;
  }
  else if ( fSinThetaForm==kSinTheta ) {
    kSquared = true;
  }    
  else
    assert(false);

  
  m_BP->SetMNS( m_osc_pars[ks12], m_osc_pars[ks13], m_osc_pars[ks23], m_osc_pars[kdm12], m_osc_pars[kdm32], m_osc_pars[kcp], energy_GeV,  kSquared, iflux );
  SetT2K(iflux);
  
  // neutrino types
  // 1:e 2:mu 3:tau   -1: e_bar -2: mu_bar -3: tau_bar
  int nu_in = iflux;
  int nu_out = ixsec;

  double prob = m_BP->GetProb( nu_in, nu_out );
  
  //   if (iflux==2 && ixsec==2) {
  //std::cout << "BP: P*w=" << prob*weight << "= p(" << prob << ") x w(" << weight << ") :  iflux=" << iflux << " ixsec=" << ixsec << " E=" << energy_GeV << " mode=" << mode << std::endl;
  //std::cin.get();
  //   }
  
  return prob*weight;
}

Weight* WeightBargerOsc::CloneWeight( TChain* source ) {

  WeightBargerOsc* clone = new WeightBargerOsc( fSinThetaForm, source, m_nuE_GeV_var, m_nuflux_var, m_nuxsec_var, m_mode_var, m_weight_var );
  clone->SetMNS( m_osc_pars[ks12], m_osc_pars[ks13], m_osc_pars[ks23], m_osc_pars[kdm12], m_osc_pars[kdm32], m_osc_pars[kcp], fSinThetaForm );
  return clone;
  
}

void WeightBargerOsc::AssignPars( std::string s12name, std::string s13name, std::string s23name, std::string dm12name, std::string dm32name, std::string cpdname, ParameterManager* parman ) {
  if ( parman->IsParameterDefined( s12name ) ) {
    parnames[ks12] = s12name;
    poscpars[ks12] = parman->GetParameter( s12name );
    parid[ks12] = parman->GetParameterIndex( s12name );
  }

  if ( parman->IsParameterDefined( s13name ) ) {
    parnames[ks13] = s13name;
    poscpars[ks13] = parman->GetParameter( s13name );
    parid[ks13] = parman->GetParameterIndex( s13name );
  }

  if ( parman->IsParameterDefined( s23name ) ) {
    parnames[ks23] = s23name;
    poscpars[ks23] = parman->GetParameter( s23name );
    parid[ks23] = parman->GetParameterIndex( s23name );
  }

  if ( parman->IsParameterDefined( dm12name ) ) {
    parnames[kdm12] = dm12name;
    poscpars[kdm12] = parman->GetParameter( dm12name );
    parid[kdm12] = parman->GetParameterIndex( dm12name );
  }

  if ( parman->IsParameterDefined( dm32name ) ) {
    parnames[kdm32] = dm32name;
    poscpars[kdm32] = parman->GetParameter( dm32name );
    parid[kdm32] = parman->GetParameterIndex( dm32name );
  }

  if ( parman->IsParameterDefined( cpdname ) ) {
    parnames[kcp] = cpdname;
    poscpars[kcp] = parman->GetParameter( cpdname );
    parid[kcp] = parman->GetParameterIndex( cpdname );
  }
  fParsAssigned = true;
}

void WeightBargerOsc::UpdateParameters( ParameterManager* parman ) {
  if ( parman==NULL ) {
    // used stored instances
    for (int i=0; i<kNumPars; i++) {
      m_osc_pars[i] = poscpars[i]->GetValue();
    }
  }
  else {
    // read out of passed in container
    for (int i=0; i<kNumPars; i++)
      m_osc_pars[i] = parman->GetParameterFromIndex( parid[i] )->GetValue();
  }
}
