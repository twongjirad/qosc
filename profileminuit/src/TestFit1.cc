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

#include "TestFit1.hh"
#include <cmath>
#include <assert.h>

using namespace qosc;

TestFit1::TestFit1() 
  : IterativeSolver( 2, IterativeSolver::kNewton )
{
  // constructor defines fit as having 2 parameters and to use the Newton method.
}

TestFit1::~TestFit1(){}

void TestFit1::CalculateF( double* x, double* F ) {
  // Calculate the Value
//   m_F[0] = x[0] - (2+x[0]*x[1]/2.0-x[1])/2.0;
//   m_F[1] = x[1] - (1.5 + cos(x[1])/2.0 - x[0])/2.0;

//   F[0] = x[0]*x[0]+x[1]*x[1];
//   F[1] = x[0]*x[0]-x[1]*x[1];

  // Example (3): from pdf stored in Papers
  // Solution is (1,1) with initial point (0.7,0.7)
  if ( IterativeSolver::fMethod==IterativeSolver::kNewton ) {
    F[0] = x[0]*x[0] - x[1]*x[1]*x[1]*x[1]; // x_1^2 - x_2^4
    F[1] = x[0] - x[1]*x[1]*x[1]; // x_1 - x_2^3
  }
  else if ( IterativeSolver::fMethod==IterativeSolver::kJacobi ) {
    // Haven't found an appropriate test equation/form
    F[0] = x[0] + pow(x[0],2)/pow(x[1],4) - 1;
    F[1] = x[1] + x[0]/x[1] - 1/x[1];
  }
  else {
    assert(false);
  }

}

void TestFit1::BuildJacobian( double* x, double** J ) {
//   m_J[0][0] = 2-x[1]/2.0;
//   m_J[0][1] = 1-x[0]/2;
//   m_J[1][0] = 1.0;
//   m_J[1][1] = 2.0 + sin(x[1])/2.0;

//   m_J[0][0] = 2*x[0];
//   m_J[0][1] = 2*x[1];
//   m_J[1][0] = 2*x[0];
//   m_J[1][1] = -2*x[1];

  // Example (3)
  J[0][0] = 2*x[0]; J[0][1] = -4*x[1]*x[1]*x[1];
  J[1][0] = 1;      J[1][1] = -3*x[1]*x[1];
    
}

