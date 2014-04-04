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
 * \class ProfileMinuitI
 *  \group ProfileMinuit
 *  \brief Abstract Interface class between ProfileMinuit Fitter class and some generic AnalysisClass
 *
 *  Meant to guide the user in writing the interface between the ProfileMinuit Fitting class and their analysis class.
 *  Implements all of the abstract FitterI classes and defines pure virtual methods that the user must implement.
 *  These methods tell the user what needs to be loaded into the fitting class.
 * -----------------------------------------------------------------------------------------------------------------------
 */

#ifndef __ProfileMinuitI__
#define __ProfileMinuitI__

#include "FitterI.hh"

namespace qosc {

  class Fitter;
  class AnalysisClass;
  class ProfileMinuit;
  class ParameterManager;

  class ProfileMinuitI : public FitterI {

  public:

    ProfileMinuitI( Fitter* , AnalysisClass* );
    virtual ~ProfileMinuitI();

    // --------------------------------------------------
    // Implemented pure virtual functions from FitterI
    virtual void UpdateFitter();
    // --------------------------------------------------

    // --------------------------------------------------
    // Implemented pure virtual functions from FitterI
    virtual ProfileMinuit* GetFitter();
    // --------------------------------------------------
  
    // --------------------------------------------------
    // Fitter Initialization.
    // Providing abtract functions that load the required information for the fitter.
    virtual void InitializeFitter(); 
    virtual ParameterManager* UserInit_ReturnParameterManager() = 0;
    virtual int UserInit_ReturnNumberOfBins() = 0;  
    // --------------------------------------------------

    /// ------------------------------------------------------------------------
    /// Called in UpdateFitter()
    // We are breaking up this function into explicit steps to help guide the user into what to implement
    virtual void UserUpdate_PreCalculation()  = 0;
    virtual void UserUpdate_LoadExpectationBins() = 0;
    virtual void UserUpdate_LoadObservedBins() = 0;
    virtual void UserUpdate_LoadFijValues() = 0;
    virtual double UserUpdate_ReturnUserPenalty()  = 0;
    virtual void UserUpdate_PostCalculation() = 0;
    /// ------------------------------------------------------------------------

  protected:

    ProfileMinuit* RetypeFitter( Fitter* base_fitter );

  };

}

#endif
