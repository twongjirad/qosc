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

#include "TestModel.hh"
#include <cmath>
#include "TH1D.h"
#include "TRandom3.h"
#include "ParameterManager.hh"
#include "BasicParameter.hh"

using namespace qosc;

TestModel::TestModel( int numBins, double trueMean, double trueSigma, double trueSignalMax, double trueBGMax, double trueBGConstant, bool allminuit ) {
  m_numBins = numBins;
  m_parameters = new ParameterManager();
  
  InitParameters( trueMean, trueSigma, trueSignalMax, trueBGMax, trueBGConstant, allminuit );
  m_parameters->SetParameterValue( "signal_mean", 0. );
  m_parameters->SetParameterValue( "signal_sigma", 0. );
  m_parameters->SetParameterValue( "signal_max", 0. );
  m_parameters->SetParameterValue( "bg_const", 0.0 );
  m_parameters->SetParameterValue( "bg_max", 0.0 );

  MakeData( true );

}


TestModel::~TestModel() {
}

void TestModel::InitParameters( double trueMean, double trueSigma, double trueSignalMax, double trueBGMax, double trueBGConstant, bool allminuit ) {

  m_signal_mean = trueMean;
  m_signal_sigma = trueSigma;
  m_signal_max = trueSignalMax;
  m_bg_const = trueBGConstant;
  m_bg_max = trueBGMax;

  // Signal Parameters
  m_parameters->RegisterParameter( new BasicParameter( "signal_mean", 0.0, 1.0, kMinuitTerm ) ); /// Signal mean
  m_parameters->RegisterParameter( new BasicParameter( "signal_sigma", 0.0, 1.0, kMinuitTerm ) ); /// Signal sigma
  if ( allminuit )
    m_parameters->RegisterParameter( new BasicParameter( "signal_max", 0.0, 1.0, kMinuitTerm ) ); /// Signal Amplitude
  else
    m_parameters->RegisterParameter( new BasicParameter( "signal_max", 0.0, 1.0, kNewtonTerm ) ); /// Signal Amplitude


  // Background Parameters
  m_parameters->RegisterParameter( new BasicParameter( "bg_const", 0.0, 1.0, kMinuitTerm ) ); // Background constant
  if ( allminuit )
    m_parameters->RegisterParameter( new BasicParameter( "bg_max", 0.0, 1.0, kMinuitTerm ) ); // Background Amplitude
  else
    m_parameters->RegisterParameter( new BasicParameter( "bg_max", 0.0, 1.0, kNewtonTerm ) ); // Background Amplitude


  // Set Parameter Ordering
  m_parameters->DefineParameterOrdering( "signal_mean;signal_sigma;signal_max;bg_const;bg_max" );

  // Set Par Bounds
  m_parameters->GetParameter( "signal_mean" )->SetBounds( -1.0, 6.0 );
  m_parameters->GetParameter( "signal_sigma" )->SetBounds( -1.0, 6.0 );
  m_parameters->GetParameter( "signal_max" )->SetBounds( -1.0, 6.0 );
  m_parameters->GetParameter( "bg_const" )->SetBounds( -1.0, 6.0 );
  m_parameters->GetParameter( "bg_max" )->SetBounds( -1.0, 6.0 );

  m_parameters->Initialize();
}

void TestModel::MakeData( bool useExpectationForData ) {
  kUseExpectationForData=useExpectationForData;

  // assumes parameters have been loaded and initialized to true values
  TRandom3 rand(0);
  
  dataHist = new TH1D( "testmodel_data", "Pseudo-Data", m_numBins, 0, m_numBins );
  for (int n=0; n<m_numBins; n++) {
    double expectation = GetExpectation( n );
    if ( kUseExpectationForData ) {
      dataHist->SetBinContent( n+1, expectation );
    }
    else {
      dataHist->SetBinContent( n+1, rand.Poisson( expectation ) );
    }
  }
}

double TestModel::GetDataBin( int n ) {
  return dataHist->GetBinContent( n+1 );
}

double TestModel::GetExpectation( int n ) {
  double bg_max = m_bg_max*(1 + m_parameters->GetParameterValue( "bg_max" ));
  double bg_const = m_bg_const*(1 + m_parameters->GetParameterValue( "bg_const" ) );
  double signal_mean = m_signal_mean*(1+m_parameters->GetParameterValue( "signal_mean" ));
  double signal_sigma = m_signal_sigma*(1+m_parameters->GetParameterValue( "signal_sigma" ));
  double signal_max = m_signal_max*(1+m_parameters->GetParameterValue( "signal_max" ));

  double bg_expectation = bg_max*exp( -bg_const*double(n));
  double signal_expectation=signal_max*exp( -pow( (double(n)-signal_mean)/(signal_sigma),2) );
  double expectation = bg_expectation+signal_expectation;
  if ( expectation<0 ) expectation = 0.;
  return expectation;
}

double TestModel::GetNominalExpectation( int n ) {
  double bg_max = m_bg_max*(1 + 0.0);
  double bg_const = m_bg_const*(1 + m_parameters->GetParameterValue( "bg_const" ) );
  double signal_mean = m_signal_mean*(1+m_parameters->GetParameterValue( "signal_mean" ));
  double signal_sigma = m_signal_sigma*(1+m_parameters->GetParameterValue( "signal_sigma" ));
  double signal_max = m_signal_max*(1+0.0);

  double bg_expectation = bg_max*exp( -bg_const*double(n));
  double signal_expectation=signal_max*exp( -pow( (double(n)-signal_mean)/(signal_sigma),2) );
  double expectation = bg_expectation+signal_expectation;
  if ( expectation<0 ) expectation = 0.;
  return expectation;
}

double TestModel::GetFij( std::string parname, int bin ) {
  // keeping it simple
  if ( parname=="signal_max" || parname=="bg_max" )
    return 1.0; 
  else
    return 0.; // actually this is an error.
}

void TestModel::SetSignalMean( double mean ) {
  m_parameters->SetParameterValue( "signal_mean", mean );
}

double TestModel::GetSignalMean() {
  return m_parameters->GetParameterValue( "signal_mean" );
}

void TestModel::SetSignalSigma( double sigma ) {
  m_parameters->SetParameterValue( "signal_sigma", sigma );
}

double TestModel::GetSignalSigma() {
  return m_parameters->GetParameterValue( "signal_sigma" );
}

void TestModel::SetSignalMax( double max ) {
  m_parameters->SetParameterValue( "signal_max", max );
}

double TestModel::GetSignalMax() {
  return m_parameters->GetParameterValue( "signal_max" );
}

void TestModel::SetBackgroundMax( double max ) {
  m_parameters->SetParameterValue( "bg_max", max );
}

double TestModel::GetBackgroundMax() {
  return m_parameters->GetParameterValue( "bg_max" );
}

void TestModel::SetBackgroundConstant( double bgconst ) {
  m_parameters->SetParameterValue( "bg_const", bgconst );
}

double TestModel::GetBackgroundConstant() {
  return m_parameters->GetParameterValue( "bg_const" );
}

void TestModel::Print() {
  std::cout << "------------------------------------------------" << std::endl;
  std::cout << "Test model expectation." << std::endl;
  for (int n=0; n<m_numBins; n++) {
    std::cout << "[" << n << "] : " << GetExpectation( n ) << std::endl;
  }
}
