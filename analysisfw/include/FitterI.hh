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
/*
 * -----------------------------------------------------------------------------------------------------------------------
 * \class FitterI
 *  \group AnalysisYakuin
 *  \brief Abstract Interface class between AnalysisClass implementations and FitterClass implementations
 *
 *  Every realization will have to implement this.  Knows about the gory details of the fitter and the samples/parameters.
 *  Knows how to get what the fitter needs to do its thing.
 * -----------------------------------------------------------------------------------------------------------------------
 */


#ifndef __FitterI__
#define __FitterI__

namespace qosc {

  class Fitter;
  class AnalysisClass;

  class FitterI {

  public:
    FitterI( Fitter*, AnalysisClass* );
    virtual ~FitterI();

  public:
    virtual void UpdateFitter() = 0;
    virtual void UpdateModel() = 0;

  protected:
    AnalysisClass* m_userAnalysis;
  public:
    AnalysisClass* GetAnalysisClass() { return m_userAnalysis; };

  protected:
    Fitter* m_userFitter;
  public:
    Fitter* GetFitter() { return m_userFitter; };

  };

}

#endif
