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
* \class Sample
* \ingroup AnalysisTools
* \brief Base class that represents the sample bins.
*
* What the AnalysisBase class require of this class:
*  (1) That the bins can be assigned an index ranging from 0 to N which the Analysis class will
*      pass to the minimizer.
*  (2) That the class will take care of it's own filling. (though there are gains to fill simultaneously across other samples.)
*  (3) That the sample can coordinate with systematic errors to fill the bins. The container class of systematic error terms,
*      SysTermManager, will be passed to the samples. The user will have to coordinate the sample and sys terms to know how
*      to work together. This doesn't seem too onerous of a task.
* Right now, I'm going to base this class around the assumption that the bins are filled using
* a ROOT tree. This is not a terrible assumption. But it would be nice to have some flexibility
* here. For instance, one might have a text file of bin data.
*
* In any case, I use my PDF library here, which skips the boiler plate of setting up
* access to different ROOT variables via the TTreeFormula class.
*
* However, a less complicated method would be to just simply setup chains to trees.
* Or for a shortcut, one can simply load previously generated histograms. This might be 
* what the data wants.
*
*
* --------------------------------------------------------------------------------------------*/


#ifndef __Sample__
#define __Sample__


#include <string>
#include <map>
#include "TRandom3.h"

#include "Hist.hh"
#include "ParameterManager.hh"

namespace qosc {

  class Sample {

  public:
    Sample( std::string m_sample_name );
    virtual ~Sample();

  protected:

    std::string m_sample_name;
    bool fSampleActive;

  public:

    std::string GetName() { return m_sample_name; };
    bool IsActive() { return fSampleActive; };
    void SetActiveFlag( bool active ) { fSampleActive = active; };
    virtual void FillBins( ParameterManager* systerms, bool fillnominal=false ) = 0; ///< you must implement this.
    virtual void MakeFakeDataSet( ParameterManager* systerms, bool includeStatVariation=true, bool includePullVariation=true ); ///< Fill the data hist using MC expectation
    virtual void MakeFakeDataSet( ParameterManager* systerms, bool includeStatVariation, std::map< std::string, double >& parvalues ); ///< Fill the data hist using MC expectation
    virtual void SetSeed( int seed );

    // Get Bin Functions
    double GetExpectedContentByGlobalIndex( int globalindex) { return GetContentByGlobalIndex( globalindex, histogram ); };
    double GetExpectedContentByLocalIndex( int localindex ) { return GetContentByLocalIndex( localindex, histogram ); };
    double GetExpectedContentByHistIndex( int histindex ) { return GetContentByHistIndex( histindex, histogram ); };
    double GetDataContentByGlobalIndex( int globalindex) { return GetContentByGlobalIndex( globalindex, datahist ); };
    double GetDataContentByLocalIndex( int localindex ) { return GetContentByLocalIndex( localindex, datahist ); };
    double GetDataContentByHistIndex( int histindex ) { return GetContentByHistIndex( histindex, datahist ); };

    // Set Bin Functions
    void SetExpectedContentByGlobalIndex( int globalindex, double value ) { SetContentByGlobalIndex( globalindex, value, histogram ); };
    void SetExpectedContentByLocalIndex( int localindex, double value ) { SetContentByLocalIndex( localindex, value, histogram ); };
    void SetExpectedContentByHistIndex( int histindex, double value ) { SetContentByHistIndex( histindex, value, histogram ); };
    void SetDataContentByGlobalIndex( int globalindex, double value ) { SetContentByGlobalIndex( globalindex, value, datahist ); };
    void SetDataContentByLocalIndex( int localindex, double value ) { SetContentByLocalIndex( localindex, value, datahist ); };
    void SetDataContentByHistIndex( int histindex, double value ) { SetContentByHistIndex( histindex, value, datahist ); };


    void WriteHistogram(); ///< Save the histogram and its components
    virtual bool LoadSampleFromFile( TFile* file, double scalefactor=1.0 ) { return true; }; // does nothing, you'll have to override this.
    virtual void PrintSampleInfo( bool printBinInfo=false );

  protected:
    double GetContentByGlobalIndex( int globalindex, Hist* hist);
    double GetContentByLocalIndex( int localindex, Hist* hist);
    double GetContentByHistIndex( int histindex, Hist* hist );

    void SetContentByGlobalIndex( int globalindex, double value, Hist* hist );
    void SetContentByLocalIndex( int localindex, double value, Hist* hist );
    void SetContentByHistIndex( int histindex, double value, Hist* hist );

  protected:
    Hist* histogram; ///< histogram which represent the bins values in the fit. For expectation. Has 1D and 2D support. (note: blind to ROOT tree integration)
    Hist* datahist;  ///< histogram for the data. Will always mirror the expectation histogram above.
  public:
    void SetSampleHistogram( Hist* hist );
    Hist* GetSampleHistogram() { return histogram; };
    Hist* GetDataHistogram() { return datahist; };
    int GetNumberOfHistogramBins() { return histogram->GetNumberOfBins(); };
  protected:
    void CheckHistogram(); ///< tests if histogram is OK.

  protected:
    int m_SampleID;
  public:
    int GetID() { return m_SampleID; };
    void SetID( int id ) { m_SampleID = id; };

  protected:
    int fNumberOfSampleBins;
    std::map< int, bool > m_histbin_active; ///< Set status of histogram bin. (allows one to toggle bins on and off).
    typedef std::map< int, bool >::iterator BinActiveDictIter;
  public:
    bool IsHistBinActive( int histbin ) {
      BinActiveDictIter it = m_histbin_active.find( histbin );
      if ( it!=m_histbin_active.end() )
	return (*it).second;
      else
	return true; // by default, if status not sepcifically set, we assume the bin to be active.
    }
    int GetNumberOfBins() { return fNumberOfSampleBins; };

  protected:
    std::map< int, int > m_local_to_hist; ///< maps sample's local bins (that is ordered set of active bins) to the bins in the PDF bins.
    std::map< int, int > m_global_to_local; ///< maps global analysis bin (assigned by external class) to local bin number
    std::map< int, int > m_local_to_global; ///< local bin number to global bin number
    typedef std::map<int,int>::iterator BinIndexDictIter;
  public:
    void IndexBins(); ///< pairs root bins to local bins
    void PairGlobalBinToLocalBin( int globalbin, int localbin ) {
      m_global_to_local[globalbin] = localbin;
      m_local_to_global[localbin] = globalbin;
    };
    int GetLocalIndexFromGlobal( int globalbin ) { return m_global_to_local[globalbin]; };


    // ------------------------------------------------------------
    // List of Parameters that apply to this sample
  protected:
    std::set< std::string > m_par_list;
  public:
    void AddSampleParameterName( std::string parname ) { 
      m_par_list.insert( parname );
    };
    bool DoesParameterApplyToSample( std::string parname ) {
      if ( m_par_list.find( parname )!=m_par_list.end() ) return true;
      return false;
    };
  
  
    //   // ------------------------------------------------------------
    //   // Fij histogram management
    // protected:
    //   std::map< std::string, Hist*> m_hist_fij_dict; ///< Fij histograms for each analysis sample: <parameter, Fij histogram*>
    //   typedef std::map< std::string, Hist*>::iterator FijDictIter;
    //   FijDictIter FijDictBegin() { return m_hist_fij_dict.begin(); };
    //   FijDictIter FijDictEnd() { return m_hist_fij_dict.end(); };
    // public:
    //   Hist* GetParameterFijHist( std::string parname ) { 
    //     FijDictIter it = m_hist_fij_dict.find( parname );
    //     if ( it!=FijDictEnd() ) return (*it).second;
    //     else { 
    //       std::cout << "GetParameterFijHist is returning null. looking for '" << parname << "'" << std::endl; 
    //       return NULL;
    //     }
    //   };
    //   virtual void DefineFijHistogram( std::string parname, FijParameter* par ); ///< Creates a Hist class ( simple Hist concrete type with histogram ) for the Fij
    //   bool IsFijHistDefinedForParameter( std::string parname ) {
    //     if ( GetParameterFijHist(parname) ) return true;
    //     else return false;
    //   };
    //   virtual void UpdateFijHistogram( std::string parname );
    //   void WriteFijHistograms();

    // ------------------------------------------------------------
    // Verbosity
  protected:
    int fVerbose;
  public:
    void SetVerbose( int verbose ) { fVerbose=verbose; };
    int GetVerbose() { return fVerbose; };

    // ------------------------------------------------------------
    // Random Number Generator
    // Used to create fake MC data
  protected:
    unsigned int m_seed; ///< random seed for generator
    TRandom3* m_rand_gen;
  public:
    void SetRandomGenSeed( unsigned int seed ) {
      if ( m_rand_gen ) delete m_rand_gen;
      m_rand_gen = new TRandom3( seed );
      m_seed = seed;
    };


  public:
    Double_t* GetExpectedBinAddressFromGlobalIndex( int global_index );
    Double_t* GetExpectedBinAddressFromLocalIndex( int local_index );
    Double_t* GetObservedBinAddressFromGlobalIndex( int global_index );
    Double_t* GetObservedBinAddressFromLocalIndex( int local_index );

    // ---------------------------------------------------------------------------
    // Expose Information about UserBinInfo Objects stored in Hist Histogram
    void GetListOfUserBinInfoNames( std::vector< std::string >& bin_info_names );
    UserBinInfo* GetMasterListBinInfoInstance( std::string bin_info_name );

  };

}
#endif
