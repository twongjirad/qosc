//-*- mode:c++; c-basic-offset:2;   -*-
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
* \class SampleManagerROOTPDF
* \ingroup AnalysisTools
* \brief Coordinates Sample Histograms built from ROOT Trees
*
* Responsibilities:
* (1) Act as container for Sample class instances. So dictionary routines, get/set stuff.
* (2) Organize the arrangement of bins. The AnalysisBase class will require that the bins
*     be linearized and indexed, ranging from 0 to N.
* (3) Coordinate the filling of the sample bins
*     By default, it calls the pure virtual Sample::Fill() function.
*     But the user can overload SampleManager::FillSamples() to coordinate sample fills.
*
* -------------------------------------------------------------------------------------------*/

#ifndef __SampleManagerROOTPDF__
#define __SampleManagerROOTPDF__

#include "SampleManager.hh"

#include "HistRootVariable.hh"
#include "HistCoordinator.hh"
#include "SampleROOTPDF.hh"

class TChain;

namespace qosc {

  class SampleManagerROOTPDF : public SampleManager {

  public:
    SampleManagerROOTPDF();
    virtual ~SampleManagerROOTPDF();

    virtual void RegisterSample( std::string name, SampleROOTPDF* asample ); ///< store a sample. bin definitions found in PDFRootVariablePDF
    void AddRootChain( std::string name, TChain* chain ); ///< store ROOT Chain in contrainer
    void AddRootChain( std::string name, void* chain ); ///< store ROOT Chain in contrainer (used for python bindings)
    virtual void FillSampleBins( ParameterManager* systerms, bool fillnominal=false ); ///< Required Code from SampleManager.hh
    virtual bool LoadSamplesFromFile( TFile* file, double scalefactor=1.0 ); ///< Load samples from file instead of building them from the chains
    void SetChainBounds( std::string chainname, int start, int end ); ///< only fill the sample bins using the entries found between start and end
    void RebuildHistogramsFromChains() { fBuildHistogramsFromChains = true; }; ///< force histograms to use the data chains to fill instead of other sources
    void AlwaysRebuildHistogramsFromChains() { fAlwaysRebuildHistogramsFromChains = true; }; ///< forces sample histograms to rebuilt by extracting info from ROOT chain

  protected:

    bool fBuildHistogramsFromChains;
    bool fAlwaysRebuildHistogramsFromChains;
    bool loadingFromFile;
    TFile* m_sourcefile;

    std::map< std::string, TChain* > m_chain_dict;
    typedef std::map< std::string, TChain* >::iterator ChainDictIter;
    ChainDictIter ChainDictBegin() { return m_chain_dict.begin(); };
    ChainDictIter ChainDictEnd() { return m_chain_dict.end(); };
  public:
    TChain* GetChain( std::string name ) { 
      ChainDictIter it = m_chain_dict.find( name );
      if ( it!=ChainDictEnd() ) return (*it).second;
      else return NULL;
    };
    bool IsChainDefined( std::string name ) { 
      if ( GetChain(name) ) return true;
      else return false;
    };

    std::map< std::string, int > m_chain_start;
    std::map< std::string, int > m_chain_end;

    void LoadSampleAndFijHistsFromFile( ParameterManager* parameters );
    bool DoSamplesHaveEventReweightParameter( ParameterManager* parameters );

  public:
    void Print();

    // for python bindings
    SampleROOTPDF* GetSample( std::string name ) { return dynamic_cast<SampleROOTPDF*>( SampleManager::GetSample(name) ); };

  };

}

#endif
