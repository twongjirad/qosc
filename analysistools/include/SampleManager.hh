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
* \class SampleManager
* \ingroup AnalysisTools
* \brief Class responsible for coordinating the samples
*
* Responsibilities:
* (1) Act as container for Sample class instances. So dictionary routines, get/set stuff.
* (2) Organize the arrangement of bins. The AnalysisBase class will require that the bins
*     be linearized and indexed, ranging from 0 to N.
* (3) Coordinate the filling of the sample bins
*     By default, it calls the pure virtual Sample::Fill() function.
*     But the user can overload SampleManager::FillSamples() to coordinate sample fills.
*
* Verbose levels:
*  0: quiet
*  1: summarize updates
*  2: more info
*  3: stop at each update of samples
* -------------------------------------------------------------------------------------------*/

#ifndef __SampleManager__
#define __SampleManager__

#include "Sample.hh"
#include "ParameterManager.hh"

namespace qosc {

  class SampleManager { 

  public:
    SampleManager();
    virtual ~SampleManager();

    void RegisterSample( std::string name, Sample* asample );
    void SetupBins(); ///< right now just calls set sample order using sample order in dictionary iterator

    virtual void FillSampleBins( ParameterManager* systerms, bool fillnominal=false ); ///< AnalysisBase expects bins to fill given Parameter pull values.
    virtual void SetSeed( int seed );
    virtual void MakeFakeDataSet( ParameterManager* systerms, bool includeStatVariation=true, bool includePullVariation=true );
    virtual void WriteHistograms();
    virtual bool LoadSamplesFromFile( TFile* file, double scalefactor=1.0 );
    void SetVerbose( int verbose ) { fVerbose = verbose; };
    int GetVerbose() { return fVerbose; };

    virtual void Print();
    virtual void PrintBins();

    int GetNumberOfBins();

      
    int GetNumberOfSampleBins( std::string samplename );
    double GetExpectedBinContent( int binindex );
    double GetDataBinContent( int binindex );
    void SetDataBinContent( int binindex, double data );

    void AddConfigurationTreeToFile( TFile* afile );
    void AddExpectedBinsToTree( TTree* atree );
    void AddObservedBinsToTree( TTree* atree );
    std::string GetExpectedBinBranchName( int global_index );
    std::string GetObservedBinBranchName( int global_index );

  public:  
    std::map< std::string, Sample* > m_sample_dict;
    typedef std::map< std::string, Sample* >::iterator SampleDictIter;
    SampleDictIter SampleDictBegin() { return m_sample_dict.begin(); };
    SampleDictIter SampleDictEnd() { return m_sample_dict.end(); };
  public:
    unsigned int GetNumberOfSamples() { return m_sample_dict.size(); };
    Sample* GetSample( std::string name );
    bool IsSampleDefined( std::string name ) {
      Sample* asample = GetSample( name );
      if ( asample ) return true;
      else return false;
    };

  protected:
    bool fSampleOrdered;
    std::vector< std::string > m_sample_order;
    std::map< std::string, int > m_sample_id;
    std::map< int, std::string > m_indexed_samples;
  public:
    typedef std::vector< std::string >::iterator SampleListIter;
    SampleListIter SampleListBegin() { return m_sample_order.begin(); };
    SampleListIter SampleListEnd() { return m_sample_order.end(); };
    void SetSampleOrder( std::string samplelist );
    Sample* GetSampleFromIndex( int index ) { return GetSample( m_indexed_samples[index] ); };
    int GetSampleIndex( std::string samplename ) { return m_sample_id[ samplename ]; };

  protected:
    int fTotalNumberOfBins; ///< Total Number of Bins
    std::map< int, std::string > m_bin_index_dict; ///< ties bin index to a sample name
  public:
    Sample* GetSampleFromBinIndex( int binnum ) {
      std::map< int, std::string >::iterator it = m_bin_index_dict.find( binnum );
      if ( it!=m_bin_index_dict.end() ) return GetSample((*it).second);
      else return NULL;
    };

  protected:


  protected:
    int fVerbose;

  };

}

#endif
