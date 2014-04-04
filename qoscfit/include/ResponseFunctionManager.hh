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
/* --------------------------------------------------------------------------------------------
* \class ResponseFunctionManager
* \ingroup QoscFit
* \brief Container/Manager for Response Function objects
*
* Takes as input, a sample setup parser instance along with a fully initialize parameter manager.
* From the sample partser it extracts the TruthBinInfo instance for each define sample.
* From the sample manager, it grabs a list of the parameters.
* At a maximum, the number of response functions equals: (num truth bins) x ( num of pars )
* But this population of functions in this space will be very small. We will need a way
*  to optimize the loop through the response functions. Another class will do that.
*  This class is going to be pretty vanilla. 
* -------------------------------------------------------------------------------------------*/

#ifndef __ResponseFunctionManager__
#define __ResponseFunctionManager__

#include <string>
#include <map>
#include <vector>
#include <set>

namespace qosc {

  class SampleManager;
  class ParameterManager;
  class ResponseFunction;
  class SampleSetupParser;
  class ResponseAndOscBinMap;
  class TruthBinInfo;

  class ResponseFunctionManager {

  public:
    ResponseFunctionManager();
    ResponseFunctionManager( SampleManager* sampleman );
    virtual ~ResponseFunctionManager();

    void Initialize( SampleManager* sampleman );
    ResponseFunction* GetResponseFunction( int sampleid, int recon_bin_index, int bin_info_index, int spline_truth_bin_index );
    ResponseFunction* GetResponseFunction( int sampleid, int recon_bin_index, int truth_bin_index );
    ResponseFunction* GetResponseFunction( std::string samplename, int recon_bin_index, int truth_bin_index );
    void LoadResponseFunction( std::string samplename, int rec_bin, int truth_bin, ResponseFunction* response );
    void LoadResponseFunction( int sampleid, int rec_bin, int truth_bin, ResponseFunction* response );
    void LoadResponseFunction( int sampleid, int recon_bin_index, std::string bin_info_name, int spline_truth_bin_index, ResponseFunction* response );
    void LoadResponseFunction( std::string samplename, int recon_bin_index, std::string bin_info_name, int spline_truth_bin_index, ResponseFunction* response );
    void PrintShortReport();
    int GetNumberOfTruthBins( std::string sample );
    int GetNumberOfReconBins( std::string sample );
    bool IsSampleDefinedInManager( std::string sample );
    int GetSampleID( std::string sample );
    void Optimize(double spine_thresh=1.0e-6); ///< Removes unnecesary response parameters

  protected:
    void DestroyResponseFunctions( int sampleid );

    SampleManager* m_sampleman;
  
    int m_numsamples;

    // Data structures
    // A response function for each true energy bin 
    typedef struct sample_response_ptrs {
      std::string samplename;
      int sampleid;
      int n_truth_bins;   // for osc
      int n_rec_bins;     // for osc
      int n_total_funcs;  // n functions
      int n_bin_infos;
      int n_max_truthinfo;
      TruthBinInfo* oscInfoTemplate;
      TruthBinInfo** binInfoTemplates; // [n_bin_infos]
      std::map< std::string, int > spline_info_lookup;
      ResponseAndOscBinMap* oscToSplineMap;
      ResponseFunction** resp_ptrs; // Will become an array of Pointers using pointer math to access with target structure: [rec][n_bin_infos][true];
    };
    sample_response_ptrs* m_sample_func_list; // one for each sample AND bin info
    std::map< std::string, int > m_sample_lookup;

  protected:
    std::string m_truth_bin_info_name;
  public:
    void GetListOfSplineBinInfoNames( std::string sample, std::set< std::string >& bininfo );
    std::string GetTruthBinInfoName() { return m_truth_bin_info_name; };

  };

}

#endif
