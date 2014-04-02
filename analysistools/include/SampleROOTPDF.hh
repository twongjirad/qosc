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
/* --------------------------------------------------------------------------------------------
* \class SampleROOTPDF
* \ingroup AnalysisTools
* \brief Class with functions needed for interface with ROOTPDF variables
* -------------------------------------------------------------------------------------------*/

#ifndef __SampleROOTPDF__
#define __SampleROOTPDF__


#include <string>
#include "Sample.hh"
#include "Hist.hh"

namespace qosc {

  class SampleROOTPDF : public Sample {
  
  public:
    SampleROOTPDF( std::string sampleName, std::string chainName );
    virtual ~SampleROOTPDF();

    virtual void PrintSampleInfo();

  protected:
    std::string m_chain_name;
  public:
    void FillBins( ParameterManager* systerms, bool fillnominal=false );  // This is an empty declartion. will assert if called.
    std::string GetChainName() { return m_chain_name; };
    bool DoesSampleUseChain( std::string chainname ) {
      if ( GetChainName()==chainname ) return true;
      else return false;
    };

  protected:
    double fScalingFactor;
  public:
    void SetScalingFactor( double scale ) { fScalingFactor = scale; };
    double GetScalingFactor() { return fScalingFactor; };

  };

}

#endif
