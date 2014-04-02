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
 * \class TestFit1
 *  \defgroup ProfileMinuit
 *  \brief Test implementation for the Iterative Solver. Works with tests/run_testfit1.cc.
 *
 * -----------------------------------------------------------------------------------------------------------------------
 */

#ifndef __TestFit1__
#define __TestFit1__

#include "IterativeSolver.hh"

namespace qosc {

  class TestFit1 : public IterativeSolver {
  
  public:
  
    TestFit1();
    virtual ~TestFit1();

    virtual void CalculateF( double* x, double* F );
    virtual void BuildJacobian( double* x, double** J );
  
  };

}

#endif
