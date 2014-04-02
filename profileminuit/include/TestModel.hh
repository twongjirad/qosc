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
/**
 * -----------------------------------------------------------------------------------------------------------------------
 * \class TestModel
 *  \defgroup ProfileMinuit
 *  \brief This class implements a simple model that bins a gaussian. Works with tests/run_testmodel.cc
 *
 * -----------------------------------------------------------------------------------------------------------------------
 */

#ifndef __TestModel__
#define __TestModel__

#include "AnalysisClass.hh"
#include "TH1D.h"

namespace qosc {

  class ModelParameter;
  class ParameterManager;

  class TestModel : public AnalysisClass {

  public:

    TestModel( int numBins, double trueMean, double trueSigma, double trueSignalMax, double trueBGMax, double trueBGConstant, bool allMinuit=false );
    virtual ~TestModel();

    // ----------------------------------------------------------------
    // Implementation of pure virtual methods from Analysis Class
    virtual void UpdateModel() {} ;
    // ----------------------------------------------------------------

    int m_numBins;
    bool kUseExpectationForData;
    void MakeData( bool useExpectationForData=false );
    void InitParameters( double trueMean, double trueSigma, double trueSignalMax, double trueBGMax, double trueBGConstant, bool allminuit=false );
    TH1D* dataHist;

    // Parameters
    ParameterManager* m_parameters;

    double m_signal_mean;
    double m_signal_sigma;
    double m_signal_max;
    double m_bg_const;
    double m_bg_max;


    void SetSignalMean( double mean );
    double GetSignalMean();

    void SetSignalSigma( double sigma );
    double GetSignalSigma();

    void SetSignalMax( double max );
    double GetSignalMax();

    void SetBackgroundMax( double max );
    double GetBackgroundMax();

    void SetBackgroundConstant( double bgconst );
    double GetBackgroundConstant();

    // Bin Information
    double GetExpectation( int bin );
    double GetNominalExpectation( int bin );
    double GetDataBin( int bin );
    double GetFij( std::string parname, int bin );

    void Print();
  };

}

#endif
