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

#include "ProfileMinuitI.hh"
#include <cstdlib>
#include <iostream>
#include "Fitter.hh"
#include "ProfileMinuit.hh"
#include "AnalysisClass.hh"
#include "ParameterManager.hh"

using namespace qosc;

ProfileMinuitI::ProfileMinuitI( Fitter* userFitter, AnalysisClass* userAnalysis ) 
  : FitterI( userFitter, userAnalysis )
{}

ProfileMinuitI::~ProfileMinuitI() {}

ProfileMinuit* ProfileMinuitI::RetypeFitter( Fitter* base_fitter ) {
  ProfileMinuit* fitter = dynamic_cast< ProfileMinuit* >( base_fitter );
  if ( fitter==NULL ) {
    std::cout << "ProfileMinuitI:  Passed and instance of the Fitter class not of type ProfileMinuit." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  return fitter;
}

void ProfileMinuitI::InitializeFitter() {
  ProfileMinuit* fitter = RetypeFitter( GetFitter() );
  AnalysisClass* base_analysis = GetAnalysisClass();
  fitter->SetFitterInterface( this );
  fitter->Initialize( UserInit_ReturnNumberOfBins(), UserInit_ReturnParameterManager() );
}


void ProfileMinuitI::UpdateFitter() {
  ProfileMinuit* fitter = RetypeFitter( GetFitter() );
  UserUpdate_PreCalculation();
  UserUpdate_LoadExpectationBins();
  UserUpdate_LoadObservedBins();
  UserUpdate_LoadFijValues();
  fitter->StoreLastUserPenalty( UserUpdate_ReturnUserPenalty() );
  UserUpdate_PostCalculation();
}

ProfileMinuit* ProfileMinuitI::GetFitter() {
  return (ProfileMinuit*)FitterI::GetFitter();
}
