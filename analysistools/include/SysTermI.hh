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
* \class SysTermI
* \ingroup AnalysisTools
* \brief Interface between arbitrary MC eventweight generators and the EventReweightParameters
*
* -------------------------------------------------------------------------------------------*/

#ifndef __SysTermI__
#define __SysTermI__


#include <string>

namespace qosc {

  class SysTermI {
  
  public:

    enum SysTermType { kRootFormula, kRootFunction, kDataManager, kT2KReWeight }; ///< Types of Sys Error Types, Add it here.

    SysTermI( std::string name, SysTermType systype=SysTermI::kRootFormula );
    virtual ~SysTermI();

    std::string GetName() { return m_name; };
    SysTermType GetType() { return fType; };

    virtual double GetPlus1SigmaFracError() = 0;
    virtual double GetMinus1SigmaFracError() = 0;

    // Set Pull Value
    virtual void SetPullValue( double pull, bool flag=false )  {  m_pull_value = pull; };
    virtual double GetPullValue() { return m_pull_value; };

  protected:
    double m_pull_value;

  private:
  
    std::string m_name;
    SysTermType fType;

  };
}

#endif
