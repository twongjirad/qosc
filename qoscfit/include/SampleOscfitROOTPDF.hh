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
* \class SampleOscfitROOTPDF
* \ingroup QoscFit
* \brief Implementation of the Sample abstract class for oscillation analysis framework
*
* This class has a concrete implementation of FillBins. 
* It can either build the histogram through 
*  (1) event-by-event reweighting of MC events or
*  (2) response functions (linear or cubic-spline response)
* Also includes the effect of oscillations. This is also implemented in an event-by-event reweighting
*   or through weight factors calculated from the information stored in the TruthBinInfo classes.
* -------------------------------------------------------------------------------------------*/

#ifndef __SampleOscfitROOTPDF__
#define __SampleOscfitROOTPDF__

#include <string>
#include <iostream>

#include "SampleROOTPDF.hh"
#include "SinThetaForms.hh"

class TFile;

namespace qosc {

  class Hist;
  class HistRootVariable;
  class ParameterManager;
  class WeightOscEvent;
  class WeightBargerTruthHist;
  class TruthBinInfo;
  class BinShiftParameter;

  class SampleOscfitROOTPDF : public SampleROOTPDF {
  
  public:
    SampleOscfitROOTPDF( std::string sampleName, std::string chainName );
    virtual ~SampleOscfitROOTPDF();

    virtual void FillBins( ParameterManager* systerms, bool fillnominal=false );
    virtual void CalculateFijValues( ParameterManager* parameters );
    virtual bool LoadSampleFromFile( TFile* file, double scalefactor=1.0 );
    virtual void WriteToFile();
    virtual void SetSampleHistogram( HistRootVariable* hist );
  
    // Oscillation Truth Bin interface [ Truth Bins for Oscillation Calculation ]
    int GetNumberOfTruthBinsPerTruthHist();
    int GetNumberOfTotalTruthBins();
    int GetNumberOfTruthHists();
    void GetListOfTruthHistNames( std::vector< std::string >& names );
    void GetListOfTruthHistDefinitions( std::vector< std::string >& hist_defs );
    double GetOscTruthBinExpectation( int irec, int itruth_global );
    std::string GetTruthBinHistName( int irec, int itruth_global );

    // Spline Truth Bin Interface [ Truth Bins for Spline/Response Function ]
    int GetNumberOfTotalTruthInfoBins( std::string bin_info_name );
    int GetNumberOfTruthInfoHists( std::string bin_info_name );
    int GetNumberOfTruthBinsPerTruthHistForBinInfo( std::string bin_info_name );
    void GetListOfTruthInfoHistNames( std::string bin_info_name, std::vector< std::string >& names );
    double GetTruthInfoBinExpectation( std::string bin_info_name, int irec, int itruth_global );
    std::string GetTruthBinInfoHistName( std::string bin_info_name, int irec, int itruth_global );

  public:  
    // Member for variables for saving histograms for output
    bool fSaveTruthModeHists;
    Hist** truthmode_histograms;
    void SetSaveTruthModeHistsFlag( bool savehists );
    void WriteTruthModeHists();
    void SetTemplateMode( bool ison );
    
    // OSCILLATION TOOLS
  public:
    void SetSinThetaForm( SinThetaForms sintheta ); /// going to move this to the manager level
    void UseOscEventReweight( SinThetaForms sinthetaform, std::string nuE_GeV_var, std::string nuflux_var, std::string nuxsec_var, std::string mode_var, std::string weight_var ); // move to manager
    void UseOscHistReweight(); // manager level
    void UpdateOscPars( ParameterManager* parameters );
    bool DoWeUseOscEventReweighting() { return fUseEventReweightOsc; };
    void SetOscBinInfoName( std::string bininfoname ) {  m_osc_bin_info_name=bininfoname; std::cout << "Set OscBinInfoName to " << m_osc_bin_info_name << std::endl; };
    void AddResponseBinInfoName( std::string bininfoname ) { m_response_bin_info_name.insert(bininfoname); m_response_bin_info_index[bininfoname] = m_response_bin_info_name.size()-1; }; 
    std::string GetOscBinInfoName()  { return m_osc_bin_info_name; };
    void GetResponseBinInfoNames( std::vector< std::string >& spline_info_objs ) { 
      for (  std::set< std::string >::iterator it=m_response_bin_info_name.begin(); it!=m_response_bin_info_name.end(); it++ ) spline_info_objs.push_back( *it );
    };
    int GetNumOfResponseBinInfoNames() { return m_response_bin_info_name.size(); };
    int GetResponseBinInfoIndex( std::string bin_info_name ) { return m_response_bin_info_index[bin_info_name]; };
  protected:
    WeightOscEvent* event_oscillator; ///< Responsible for oscillations using event by event reweight method
    WeightBargerTruthHist* oscillator; ///< Responsible for oscillations using histogram method (for speed gains)

    TruthBinInfo* m_truthbins; ///< A copy of the truth bin info instance. Each recon bin has one. Contains truth info needed to build the expectation.
    TruthBinInfo* m_truthbins_fij; ///< A copy of the truth bin info instance. Each recon bin has one. Contains truth info needed to build the Fijs.
    bool fUseEventReweightOsc; ///< flag if true, relies on event-by-event reweighting. if false, uses histogram reweighting
    bool fUseDifferentOscAndResponseBinning; ///< flag if true, means that the reponse function binning and the binning for oscillation reweighting is different
    bool fTruthToResponseMapMade;
    int* m_sample_truth_to_splinebin; /// dictionary from truth (oscillation binning) to response/spline binning
    //void MakeTruthBinToSplineBinMap(); 
    std::string m_osc_bin_info_name;
    std::set< std::string > m_response_bin_info_name;
    std::map< std::string, int > m_response_bin_info_index;
    int GetResponseTruthBin( int osc_truth_bin );
  
    Hist* original_histogram;
    Hist* predistortion_histogram;
    Hist* GetPredistortionHist() { return predistortion_histogram; };

    // FIJ FUNCTIONS: Reweighting by Response Functions
  protected:
    double* Fij;
    int m_fij_nelems;
    int m_fij_npars;
  public:
    void SetFij( int bin, int ipar, double fijvalue );
    double GetFij( int bin, int ipar );
    void GetBinFijs( int bin, double* fijvalues );

  protected:
    int* m_sampleid_table;
    void StoreResponseParSampleID( ParameterManager* );

    // BIN SHIFT PARAMETERS: Bin migration
  protected:
    BinShiftParameter** m_shiftpar_array;
    int m_nshiftpars;
    bool fCheckedForShiftPars;

  };

}

#endif
