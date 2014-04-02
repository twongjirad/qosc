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
 * \class FitYakuin
 *  \defgroup FitYakuin
 *  \brief 
 *
 *  Every realization will have to implement this.  Knows about the gory details of the fitter and the samples/parameters.
 *  Knows how to get what the fitter needs to do its thing.
 * -----------------------------------------------------------------------------------------------------------------------
 */

#ifndef __FitYakuin__
#define __FitYakuin__

namespace qosc {

  class Fitter;
  class AnalysisClass;
  class FitterI;

  class FitYakuin {

  public:

    FitYakuin( Fitter*, AnalysisClass*, FitterI* );
    virtual ~FitYakuin();

    void InitializeFitter();
    void RunFit();

    Fitter* GetFitter() { return m_theFitterClass; };
    AnalysisClass* GetAnalysis() { return m_theAnalysisClass; };

  protected:
  
    Fitter* m_theFitterClass; ///< Object that does the fitting
    AnalysisClass* m_theAnalysisClass; ///< Object that knows about the bins and parameters and how to manipulate them.
    FitterI* m_theFitterInterface; ///< Object that contains interface methods that knows how to pass information between the two type of classes.


  };

}

#endif
