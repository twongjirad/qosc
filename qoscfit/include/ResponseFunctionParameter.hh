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
/**
 * --------------------------------------------------------------------------------------------
 * \class ResponseFunctionParameter
 * \ingroup QoscFit
 * \brief Concrete implementation of a parameter for response functions
 *
 * -------------------------------------------------------------------------------------------*/

#ifndef __ResponseFunctionParameter__
#define __ResponseFunctionParameter__

#include "BasicParameter.hh"

#include <string>

namespace qosc {

  class ResponseFunction;
  class ResponseFunctionManager;
  class SampleManager;

  class ResponseFunctionParameter : public BasicParameter {

  public:
    ResponseFunctionParameter( std::string name, double mean, double sig, FitterParType t );
    virtual ~ResponseFunctionParameter();

    virtual std::string IsA() { return "ResponseFunctionParameter"; };
    int GetNumberOfTruthBins( std::string sample );
    int GetNumberOfReconBins( std::string sample );
    double GetResponse( std::string samplename, int recon_bin_index, int truth_bin_index, double x );
    double GetResponse( int sampleid, int recon_bin_index, int truth_bin_index, double x );
    int GetSampleID( std::string sample );

  protected:
    ResponseFunctionManager* m_funcMan;
  public:
    ResponseFunctionManager* GetResponseFunctionManager() { return m_funcMan; };
    void DefineSampleAndResponseStructure( SampleManager* sampleman );
    bool IsResponseLoaded( std::string samplename, int rec_bin, int truth_bin );
    ResponseFunction* GetResponseFunction( int sampleid, int recon_bin_index, int truth_bin_index );
    ResponseFunction* GetResponseFunction( std::string samplename, int recon_bin_index, int truth_bin_index );
    static std::string GenerateStandardResponseName( std::string samplename, int recon_bin_index, int truth_bin_index );

  };

  std::string generate_standard_response_name( std::string samplename, int recon_bin_index, int truth_bin_index );

}

#endif
