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

#include "FitYakuin.hh"
#include "Fitter.hh"
#include "FitterI.hh"
#include "AnalysisClass.hh"

using namespace qosc;

FitYakuin::FitYakuin( Fitter* theFitter, AnalysisClass* theAnalysis, FitterI* theInterface ) {
  m_theFitterClass = theFitter;
  m_theAnalysisClass = theAnalysis;
  m_theFitterInterface = theInterface;
  
  InitializeFitter();

}


FitYakuin::~FitYakuin() {}

void FitYakuin::InitializeFitter() {
}

void FitYakuin::RunFit() {
}
