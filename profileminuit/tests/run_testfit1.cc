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

#include <iostream>
#include <string>

#include "TestFit1.hh"

using namespace qosc;

int main( int iarg, char** argv ) {

  /*
   *   Verbose levels:
   *    0: Quiet (at least base class will be)
   *    1: Print out summary of fit
   *    2: Print out more information about each step
   *    3: Print out maximum amount of info
   *    4: Pause after each interation
 */
  TestFit1* tf1 = new TestFit1();
  tf1->SetStoppingTolerance( 1.0e-6 );
  tf1->SetMethod( IterativeSolver::kNewton );
  //tf1->SetMethod( IterativeSolver::kJacobi );
  //tf1->SetVerbose( 0 );
  //tf1->SetVerbose( 1 );
  //tf1->SetVerbose( 2 );
  //tf1->SetVerbose( 3 );
  tf1->SetVerbose( 4 );

  double initial_point[2] = { 0.7, 0.7 };
  tf1->SetInitialPoint( initial_point );
  double solution[2];
  tf1->RunIteration( solution );
  
  std::cout << "Initial Point: (" << initial_point[0] << ", " << initial_point[1] << ")" << std::endl;
  std::cout << "Fitted Solution: (" << solution[0] << ", " << solution[1] << ")" << std::endl;
  std::cout << "For an initial value of (0.7, 0.7 ) the answer should be (1.0,1.0)" << std::endl;

  return 0;

};
