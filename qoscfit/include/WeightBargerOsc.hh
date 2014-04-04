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
/** --------------------------------------------------------------------------------------------
 * \class WeightBargerOsc
 * \ingroup QoscFit
 * \brief Interface to R. Wendell's Prob3++ oscillation probability calculator
 * 
 * This class is a wrapper that provides an oscillation probability
 * weight derived from R. Wendell's Prob3++ class.
 * It also interfaces to a ROOT tree through the RootVariable.
 * ---------------------------------------------------------------------------------------------- */

#ifndef __WeightBargerOsc__
#define __WeightBargerOsc__

#include <string>
//#include "BargerPropagator.h"

#include "WeightOscEvent.hh"
#include "RootVariableList.hh"

#include "SinThetaForms.hh"

class TChain;
class BargerPropagator;

namespace qosc {

  class ParameterManager;
  class ModelParameter;
  
  class WeightBargerOsc : public WeightOscEvent {

  public:

    //typedef enum { kSin2Theta, kSinTheta } SinThetaForms;
    typedef enum { ks12=0, ks13, ks23, kdm12, kdm32, kcp } ParID;
    static const int kNumPars = 6;

    WeightBargerOsc( SinThetaForms sin_term_type, TChain* source, std::string nueE_GeV_var, std::string nuflux_var, std::string nuxsec_var, std::string mode_var, std::string weight_var );
    WeightBargerOsc( SinThetaForms sin_term_type, TChain* source, std::string nueE_GeV_var, std::string nuflux_var, std::string nuxsec_var, std::string mode_var, std::string weight_var,
		     std::string s12name, std::string s13name, std::string s23name, std::string dm12name, std::string dm32name, std::string cpdname, ParameterManager* parman );
    virtual ~WeightBargerOsc();

    virtual double CalculateWeight();
    virtual double CalculateWeight( double var ) { return CalculateWeight(); };
    virtual Weight* CloneWeight( TChain* source );    
    virtual void UpdateParameters( ParameterManager* );

    void SetMNS( double theta12, double theta13, double theta23, double dm12, double dm23, double CP );
    void SetMNS( double theta12, double theta13, double theta23, double dm12, double dm23, double CP, SinThetaForms atm_term_type );
    void PrintMNS();

    void SetT2K(int flavortype); //hmm, in the futre this should be changed to its own child class
    void SetPOT( double pot ) { m_POT = pot; };

    void SetSinThetaForm( SinThetaForms sin_term_type ) { fSinThetaForm=sin_term_type; };
    bool AreParsAssigned() { return fParsAssigned; };

  protected:

    BargerPropagator* m_BP;
    TChain* m_source_chain;
    RootVariableList m_rootvars;
    double m_POT;

    // Names of key variables needed by Barger Oscillator
    SinThetaForms fSinThetaForm;
    std::string m_nuE_GeV_var;
    std::string m_nuflux_var;
    std::string m_nuxsec_var;
    std::string m_mode_var;
    std::string m_weight_var;

    // PMNS Oscillation parameters
//     double m_s12; 
//     double m_s13;
//     double m_s23;
//     double m_dm12;
//     double m_dm23;
//     double m_CP;
    double m_osc_pars[6];

  protected:
    bool fParsAssigned;
    std::string parnames[6];
    int parid[6];
    ModelParameter* poscpars[6];
    void AssignPars( std::string s12name, std::string s13name, std::string s23name, std::string dm12name, std::string dm32name, std::string cpdname, ParameterManager* parman );
  };

}

#endif
