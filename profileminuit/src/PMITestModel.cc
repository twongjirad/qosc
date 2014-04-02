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

#include "PMITestModel.hh"

using namespace qosc;

PMITestModel::PMITestModel( ProfileMinuit* mmfitter, TestModel* tm ) 
  : ProfileMinuitI( mmfitter, tm )
{
  
}

PMITestModel::~PMITestModel() {}

// Initialize
ParameterManager* PMITestModel::UserInit_ReturnParameterManager() {
  return GetModel()->m_parameters;
}

int PMITestModel::UserInit_ReturnNumberOfBins() {
  std::cout << GetModel()->m_numBins << std::endl;
  return GetModel()->m_numBins;
}

// UpdateFitter
void PMITestModel::UserUpdate_PreCalculation() {
}

void PMITestModel::UserUpdate_LoadExpectationBins() {
  TestModel* tm = GetModel();
  ProfileMinuit* mm = GetFitter();
  for (int n=0; n<tm->m_numBins; n++) {
    mm->SetBinNExpected( n, tm->GetNominalExpectation( n ) );
  }
}


void PMITestModel::UserUpdate_LoadObservedBins() {
  TestModel* tm = GetModel();
  ProfileMinuit* mm = GetFitter();
 
  for (int n=0; n<tm->m_numBins; n++) {
    mm->SetBinNObserved( n, tm->GetDataBin( n ) );
  }
}

void PMITestModel::UserUpdate_LoadFijValues() {
  TestModel* tm = GetModel();
  ProfileMinuit* mm = GetFitter();
  for (int n=0; n<tm->m_numBins; n++) {
    mm->SetFij( "signal_max", n, tm->GetFij( "signal_max", n ) );
    mm->SetFij( "bg_max", n, tm->GetFij( "bg_max", n ) );
  }
}

double PMITestModel::UserUpdate_ReturnUserPenalty() {
  return 0.;
}

void PMITestModel::UserUpdate_PostCalculation() {
}

void PMITestModel::UpdateModel() {
  //TestModel* tm = GetModel();
  ProfileMinuit* mm = GetFitter();
  mm->GetParametersFromMinuit();
  mm->GetParametersFromPullMinimizer();
}
