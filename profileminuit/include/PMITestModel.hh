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
 * \class ProfileMinuit
 *  \defgroup ProfileMinuit
 *  \brief Interface between MarginMinuit and TestModel
 *
 * ---------------------------------------------------------------------------------------------------------------------- */

#ifndef __PMITestModel__
#define __PMITestModel__

#include "ProfileMinuitI.hh"
#include "ProfileMinuit.hh"
#include "TestModel.hh"
#include "FitterI.hh"

namespace qosc {

  class PMITestModel : public ProfileMinuitI {

  public:

    PMITestModel( ProfileMinuit* mmfitter, TestModel* tm );
    virtual ~PMITestModel();

    // ----------------------------------------------------------------
    // Implementations of pure virtual function from ProfileMinuitI

    // Initialization
    virtual ParameterManager* UserInit_ReturnParameterManager();
    virtual int UserInit_ReturnNumberOfBins();  

    // UpdateFitter
    virtual void UserUpdate_PreCalculation();
    virtual void UserUpdate_LoadExpectationBins();
    virtual void UserUpdate_LoadObservedBins();
    virtual void UserUpdate_LoadFijValues();
    virtual double UserUpdate_ReturnUserPenalty();
    virtual void UserUpdate_PostCalculation();
  
    // ----------------------------------------------------------------
    // Implementations of pure virtual function from FitterI
    virtual void UpdateModel();

    // ----------------------------------------------------------------
    virtual TestModel* GetModel() { return (TestModel*)FitterI::GetAnalysisClass(); };
  };

}

#endif
