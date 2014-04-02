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
 * --------------------------------------------------------------------------
 *   \class IterativeSolver
 *   \ingroup ProfileMinuit
 *   \brief Iterative Solver Virtual Class
 *   Methods based on slides in IterativeMethodsLecture.pdf.
 *   
 *   When in NEWTON mode:
 *   Attempts to solve a non-linear equation F(x) = 0. Here x and F(x) is a vector of length n.
 *   F(x) may depend on any order function of x.
 *   User must specify how two things are to be calculated:
 *      (1) F(x), given x
 *      (2) The Jacobian, defined as J[i][j] = d(F_i)/d(x_j)
 *
 *   When in JACOBIAN mode:
 *   Attempts to solve a non-linear equation x(n+1) = F(x(n)).
 *   In order to help convergence, the update to x(n+1)  = x(n) + a*(F(x(n))-x(n)) = a*F(x(n)) + x(n)*(1-a)
 *       a is adjusted such that the level of change in x(n) is less than A (default=0.1)
 *       a is iteratively reduced by 0.5, starting at 1.0. If a becomes less than 1e-10, 
 *          then problem considered divergent.
 *   User must implement how to calculate F(x).
 *   User also has to implement empty Jacobian function.
 *   Note that the problem converges based on the form of F(x). 
 *      See lecture notes on convergence requirement (which requires knowledge of solution unfortunately).
 *
 *   To Do: Implement Broyden Quasi-Newton method to speed up fitter. Also allows higher-order expectation terms.
 *
 *   Verbose levels:
 *    0: Quiet (at least base class will be)
 *    1: Print out summary of fit
 *    2: Print out more information about each step
 *    3: Print out maximum amount of info
 *    4: Pause after each interation
 * --------------------------------------------------------------------------
*/

#ifndef __IterativeSolver__
#define __IterativeSolver__

#include "TMatrixD.h"
#include "TDecompLU.h"
#include <iostream>

class TRandom3;

namespace qosc {

  class IterativeSolver {

  public:

    typedef enum { kNewton, kJacobi } Method;

    IterativeSolver();
    IterativeSolver( int ndims, Method meth=kNewton, double tolerance=1.0e-6, int maxiters=1000 );
    virtual ~IterativeSolver();

    void InitializeSetup( int ndims, Method meth, double tolerance, int maxiters);
    void ClearVariables();

    void SetInitialPoint( double* x0 ); // what the user calls to seed the fit
    void SetMethod( Method meth ) { fMethod=meth; };
    void GetCurrentSolution( double* solution );
    virtual bool RunIteration( double* solution ); // what the user calls

    void SetStoppingTolerance( double tolerance ) { m_tolerance = tolerance; };
    void SetVerbose( int verbosity ) { m_verbose = verbosity; };
    void SetIterativeMethod( Method meth ) { fMethod = meth; };

    int GetVerbose() { return m_verbose; };
    int GetNdims() { return m_ndims; };
    void PrintF();
    void PrintX();

    // Functions only for Newton Method
    void Newton_PrintJacobian();
    void Newton_PrintInverseJacobian();
    void Newton_UseDampedMethod( bool use ) { fDamped = use; };
    void Newton_UseNoInverseMethod( bool use ) { fNoInverse = use; };

    // Functions only for Jacobian Method
    void Jacobian_SetStepLimitingFactor( double a ) { m_a_steplimitfactor = a;  };

    // optional public methods that the user can implement
    virtual bool IsSolutionValid( double* x );
    virtual void MakeValidSolution( double* x );
  
  protected:

    double* m_F; ///< function, F(x)
    double* m_x; ///< vector

    double** m_J; ///< Jacobian
    double** m_invJ; ///< Inverted Jacobian

    // User Implementation
    virtual void CalculateF( double* x, double* F ) = 0; // user must provide for a new matrix each iteration
    virtual void BuildJacobian( double* x, double** J ) = 0; // user must provide for a new jacobian each iteration using the passed vector (if using newton's method)

  private:

    // Common Methods: picks next step
    bool GetNextSolution( double* x, double* F, double* x_next );
    bool DoWeStop( double* x, double* F, double* x_next, double* F_next );

    // Newton Methods and Variables
    bool GetNextSolutionNewton( double* x, double* F, double* x_next );
    void Newton_InvertJacobian(); ///< Invert Jacobian for Newton
    virtual double Newton_CalculateFSquared( double* F ); ///< used in non-inversion methods
    bool Newton_GetStep( double* F, double* dx ); ///< get step size for newton
    double Newton_CalculateBasicStep( double* F, double* dx );
    double Newton_CalculateFSquaredStep( double* F, double* dx );
    bool Newton_CheckTermination( double* x, double* F, double* x_next, double* F_next );

    // Jacobi Methods and Variables
    bool GetNextSolutionJacobi( double* x, double* F, double* x_next );
    bool Jacobi_CheckTermination( double* x, double* F, double* x_next, double* F_next );

    //#private:
  protected:
  
    // common members
    int m_ndims;
    int m_max;
    int m_iter;
    int m_verbose;
    double m_tolerance;
    Method fMethod;
    TRandom3* m_numgen;

    // Newton members
    TDecompLU* LU;
    TMatrixD* m_JMatrix;
    TMatrixD* m_invJMatrix;
    bool fDamped;
    bool fNoInverse;
    double m_Fsquared;
    double* m_theta;  

    // Jacobian members
    double m_a_steplimitfactor;

  };

}

#endif
