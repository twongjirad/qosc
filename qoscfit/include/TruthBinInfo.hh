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
* \class TruthBinInfo 
* \ingroup QoscFit
* \brief Class containing information needed to define a set of a Sample instance.
*
* The class SampleSetupParser is responsible for building this class.
* We need to index the true energy bins. Indexing by number is desirable for looping through
*  the bins repeatedly.
* -------------------------------------------------------------------------------------------*/

#ifndef __TruthBinInfo__
#define __TruthBinInfo__

#include "UserBinInfo.hh"
#include <string>
#include <vector>
#include <map>

#include "TChain.h"
#include "TH1.h"
#include "RootVariableList.hh"
#include "SampleBinInfo.hh"

namespace qosc {

  class TruthBinInfo : public UserBinInfo {

  public:

    TruthBinInfo( std::string truth_name, std::string root_formula, SampleBinInfo* binfo, TChain* source_chain );
    virtual ~TruthBinInfo();

    virtual void Initialize(); ///< setup histograms and such
    virtual void SetupInstance( SampleBinInfo* binfo );

    virtual double GetBinReweight();
    virtual void FillBinInfo( double external_weight ); // important: fills histograms
    virtual UserBinInfo* Copy( std::string name ); // important: helps spread instances to many different bins
    virtual void Reset(); ///< Clears value of histograms

    void SetChain( TChain* chain );

    int GetDimensions() { return ndims; };
    TH1* GetTemplateHist() { return m_template_hist; };

  public:
    int ndims;
    std::string m_truth_name;
    std::string m_root_formula;
    TChain* m_source_chain;
    SampleBinInfo* m_bin_def;
    RootVariableList m_rootvars;
  
    TH1* m_template_hist;

    // Temporary Cut Lists. For each cut combination, there needs to be a histogram. 
    // We store cut definitions here, before allocating the array below;
    // These variables live here, because we use a recursive formula to build the list of cuts
    struct cut_set {
      int ncuts;
      std::string* names;
      std::string* formulas;
    };
    std::vector< cut_set* > m_cutset_list;
    //   std::vector< std::string* > m_cutlist_formulas;
    //   std::vector< int > m_cutlist_ncuts;
    //   std::vector< std::string* > m_cutlist_names;


    // Here the cuts are stored in sequence. Made from combining the cuts stored in the m_cutlist_formulas vector.
    // We need to be able to acces these histogram by names, but also by truth bin index.
    struct cut_data {
      std::string name;
      std::string formula;
      std::string histname;
      int id;
      int nbins;
      int ndims;
      int start_bin; // start of truth bin indexing (indexed from 0)
      int end_bin; //index of last bin (inclusive)
      TH1* hist;
    };
    int m_ncutlists;
    cut_data* m_cut_list;
    std::map< std::string, int > m_cut_lookup; // maps name to index

    //   std::string *m_cut_lists; // final combined cut lists
    //   std::string *m_hists_names; // list of hist names;
    //   TH1** m_hists_lists; // list of histogram instances
    bool fCutsStored; // were cuts placed in root var list?

    // Recursive cut maker. Only way to traverse combinations of multiple dimensions
    bool MakeCutlist( std::vector< std::string >& cutlist, std::vector< std::string >& namelist, int* index_tracker, int& ncutlists );
    std::string MakeCombinedCutName( std::string sub_cuts[], int nsubcuts );

  public:
    int GetNumberOfTruthHists() { return m_ncutlists; };
    int GetTotalTruthBins();
    int GetNumberOfTruthBinsPerHist();
    TH1* GetTruthHist( int index ); 
    TH1* GetTruthHist( std::string cutname );
    int GetCutIndex( std::string cutname );
    int GetCutIndex( std::string sub_cuts[], int nsubcuts );
    int GetCutTruthBinStart( std::string cutname );
    int GetCutTruthBinStart( std::string subcuts[], int nsubcuts );
    void GetListOfCutNames( std::vector< std::string >& cutnames );
    void GetListOfCutDefinitions( std::vector< std::string >& cutdefs );
    double GetTruthBinContentFromGlobalIndex( int global_truth_bin );
    std::string GetCutNameFromTruthBin( int bin_num );
    std::string GetCutDefinitionFromTruthBin( int bin_num );
    void Print();

    // Functions to go to and from the disk
  public:
    virtual void LoadBinInfoFromFile( TFile* file );
    virtual void Write();
  };
}
#endif
