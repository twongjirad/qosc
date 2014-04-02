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
 * \class WeightBargerTruthHist
 * \ingroup QoscFit
 * \brief Class provides weights to oscillate expectation via TruthBinInfo class
 *
 * This class is a wrapper that provides an oscillation probability
 * weight derived from R. Wendell's Prob3++ class.
 * It takes as input the true energy histograms stored in a TruthBinInfo histogram.
 * ------------------------------------------------------------------------------------------- */

#ifndef __WeightBargerTruthHist__
#define __WeightBargerTruthHist__

#include "Weight.hh"
#include <string>
//#include "RootVariableList.hh"
#include "SinThetaForms.hh"

class TChain;
class TH1D;
class TAxis;
class BargerPropagator;

namespace qosc {

  class TruthBinInfo;
  //class HistRootVariable;

  class WeightBargerTruthHist : public Weight {

  public:

    //typedef enum { kSin2Theta, kSinTheta } SinThetaForms;
    typedef enum { ks12=0, ks13, ks32, kdm12, kdm32, kcp } ParID;
    static const int kNumPars = 6;

    WeightBargerTruthHist( TruthBinInfo* truth_bin_def );
    virtual ~WeightBargerTruthHist();

    void SetMNS( double theta12, double theta13, double theta23, double dm12, double dm23, double CP, SinThetaForms sin_term );
    void SetMNS( double theta12, double theta13, double theta23, double dm12, double dm23, double CP );
    double CalcEventsOsc( TruthBinInfo* truth_info, bool modifyhist=false ); // Takes an instance of TruthBinInfo and uses it to calculate the number of events after oscillation
    void PrintMNS();

    void SetMapMode( bool mapmode ) { fMapMode = mapmode; };
    void DontUseCache() { fUseCache = false; }; 
    void UseCache() { fUseCache = true; };
    void SetSinThetaForm( SinThetaForms sin_term ) { fSinThetaForm = sin_term; };
    SinThetaForms GetSinThetaForm() { return fSinThetaForm; };

  protected:
    // -------------------------------------------------------------------------------------
    // Configuration of oscillation probability calculation
    void Configure( TruthBinInfo* bin_def ); ///< reads truth bin configuration from 

    void BuildEnergyArray( TruthBinInfo* bin_def );
    void BuildComponentIndex( TruthBinInfo* bin_def );
    void BuildProbCache();
    // -------------------------------------------------------------------------------------

  protected:
    // -------------------------------------------------------------------------------------
    // Internal functions used to work with BargerPropagator
    void SetParameter( std::string parname,  double value );
    void SetIncomingNeutrino( std::string flavor );
    void SetOutgoingNeutrino( std::string flavor );
    void SetNEUTmode( int mode );
    void SetOscillationComp( int comp );
    void SetNeutrinoEnergy( double EnuGeV );
    double CalculateWeight() { return 0.; };
    double CalculateWeight( double variable );
    double CalculateDerivative( int par, double energy_GeV );
  
    double CalculateWeightFromCache( double energy_GeV );
    double CalculateWeightFromBarger( double energy_GeV );

    double CalculateDerivativeFromCache( int par, double energy_GeV );
    double CalculateDerivativeFromBarger( int par, double energy_GeV );
  

    void SetT2K(int flavortype); //hmm, in the futre this should be changed to its own child class
    // -------------------------------------------------------------------------------------

  protected:

    bool fParamsSet; /// flag that turns on when osc parameters set for the first time
    bool fMapMode; /// in map mode, oscillation weight is always one
    BargerPropagator* m_BP; /// oscillation weight engine.  later, can replace this to accomodate other oscillation probability models

    // variables used to set BargerPropagator
    int m_imode;
    int m_iflux;
    int m_ixsec;

    // Truth information configuration, set using the initial truthbininfo instance passed to constructor (see configure())
    int nEnergyBins;
    TAxis* EnergyBins; // energy of bins
    typedef struct component_def {
      int flux;
      int xsec;
      int mode;
    };
    component_def* component_index;

    // Probability Cache
    double* OscProbs[18]; // 18 flavor combinations: mu->mu, mu->e ... + CP conjugates
    double* GradOscProbs[18][6]; // Oscillation Probability Gradient for all 6 params

    // PMNS Oscillation parameters
    SinThetaForms fSinThetaForm;
    double m_s12; 
    double m_s13;
    double m_s23;
    double m_dm12;
    double m_dm23;
    double m_CP;

    // Cache flags
    bool fUseCache;
    bool fCalcGradient;
    int nFluxCombinations;
    bool fCacheMade;
    bool fGradCacheMade;

    // Interface to Cache
    // index comes from 1:e 2:mu 3:tau, -1:e-bar ...
    // index = 3*(abs( id )-1) + (abs(id)-1) + 9*( (1+sign(id))/2 )
    int GetProbArrayIndex( int incoming, int outgoing ); // uses internal m_iflux, m_ixsec
    void MakeCache();
    void MakeGradCache();

  protected:
    // Verbosity flag and its get/set function
    int fVerbose;
  public:
    void SetVerbose( int verbose ) { fVerbose = verbose; };
    int GetVerbose() { return fVerbose; };

  };

}

#endif
