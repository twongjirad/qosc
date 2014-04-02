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

#include "IterativeSolver.hh"
#include <iostream>
#include <cmath>
#include <assert.h>
#include "TDecompLU.h"
#include "TRandom3.h"

using namespace qosc;

IterativeSolver::IterativeSolver() {
  // Blank Constructor
}

IterativeSolver::IterativeSolver( int ndims, Method meth, double tolerance, int maxiters ) {
  InitializeSetup( ndims, meth, tolerance, maxiters );
}


IterativeSolver::~IterativeSolver(){
  ClearVariables();
}

void IterativeSolver::InitializeSetup( int ndims, Method meth, double tolerance, int maxiters ) {

  m_ndims = ndims;
  m_iter = 0;
  SetIterativeMethod( meth );
  SetVerbose( 0 );

  // Allocate common members
  m_x = new double[ndims];
  m_F = new double[ndims];
  m_tolerance = tolerance;
  m_max = maxiters;

  // Allocate the matrix (M) and Jacobian (J)
  fDamped = true;
  fNoInverse = false;
  m_theta = NULL; // used for no inverse calculate of step
  m_J = new double*[ndims]; // holds user jacobian values
  m_invJ = new double*[ndims]; // holds output of jacobian inverse values
  for (int n=0; n<ndims; n++) {
    m_J[n] = new double[ndims];
    m_invJ[n] = new double[ndims];
    m_x[n] = 0.;
    m_F[n] = 0.;
    for (int m=0; m<ndims; m++ ) {
      m_J[n][m] = 0.;
      m_invJ[n][m] = 0.;
    }
  }
  m_JMatrix = new TMatrixD( ndims, ndims ); // using this to get access to ROOT inversin methods
  m_invJMatrix = new TMatrixD( ndims, ndims ); // output of ROOT inverse
  LU = new TDecompLU( ndims ); // class that is going to do our inversions
  m_numgen = new TRandom3(1); // don't remember why this is here

  // Jacobian Initialization
  m_a_steplimitfactor = 0.1;

}


void IterativeSolver::ClearVariables() {
  
  // destroy common members
  delete [] m_x;
  delete [] m_F;

  // destory newton members
  for (int n=0; n<m_ndims; n++) {
    delete [] m_J[n];
    delete [] m_invJ[n];
  }
  delete [] m_J;
  delete [] m_invJ;
  delete LU;
  delete m_numgen;
  delete m_JMatrix;
  delete m_invJMatrix;
  if (m_theta) delete [] m_theta;

}


// ====================================================================================================
// OPTIONAL USER METHODO

bool IterativeSolver::IsSolutionValid( double* x ) {
  // This is the default. Implementing this can help convergence.
  // Iterator will shorten step until this is true
  return true;
}

void IterativeSolver::MakeValidSolution( double* x ) {
  // By default this does nothing.
  // But user can try to help convergence by modifiying solution here.
}

// ====================================================================================================
// COMMON METHODS

void IterativeSolver::SetInitialPoint( double* x0 ) {
  for (int n=0; n<m_ndims; n++) m_x[n] = x0[n];
}

void IterativeSolver::GetCurrentSolution( double* solution ) {
  for (int n=0; n<m_ndims; n++) solution[n] = m_x[n];
}

bool IterativeSolver::RunIteration( double* solution ) {

  // This must be common to both the Jacobi and Newton Methods
  if ( GetVerbose()>=1 ) {
    std::cout << "////////////////////////////////////////////////////////////////////////" << std::endl;
    std::cout << "IterativeSolver::RunIteration" << std::endl;
  }

  m_iter=0;
  double x_next[m_ndims];
  double F_next[m_ndims];
  bool converged = false;
  bool stepOK = true;
  bool terminate = false;
  double Fnorm = 0;
  bool error = false;

  while (m_iter<m_max && terminate==false ) {
    CalculateF( m_x, m_F ); // Goes out to user implementation

    stepOK = GetNextSolution( m_x, m_F, x_next ); // The iterative methods return the next step and their estimated error.
    CalculateF( x_next, F_next );

    // Barf
    if ( GetVerbose()>=2 ) {
      std::cout << "==============================================================================" << std::endl;
      std::cout << "[ Iteration #" << m_iter << " ]" << std::endl;
      PrintF();
      PrintX();
      std::cout << "Step dx: {";
      for (int n=0; n<m_ndims; n++) {
	std::cout <<  x_next[n]-m_x[n] << ", ";
      }
      std::cout << "}" << std::endl;
      if ( GetVerbose()>=3 ) {
	Newton_PrintJacobian();
	Newton_PrintInverseJacobian();
      }
    }

    if ( stepOK==false )  {
      // if calculation of next solution had an error. stop.
      terminate = true; 
      break;
    }
    else {
      // Update solution
      Fnorm = 0.;
      for (int n=0; n<m_ndims; n++) {
	m_x[n] = x_next[n];
	m_F[n] = F_next[n];
	Fnorm += m_F[n]*m_F[n];
	if ( m_x[n]!=m_x[n] ) error = true;
      }
      Fnorm = sqrt(Fnorm);
    }
    if ( Fnorm!=Fnorm ) error = true;
    
    terminate = DoWeStop( m_x, m_F, x_next, F_next );

    if ( this->GetVerbose()>=2 ) {
      std::cout << "Verbose: " << this->GetVerbose() << " " << m_verbose << " " << std::endl;
      std::cout << "[ End of iteration #" << m_iter << " ] " << std::endl;
      std::cout << "==============================================================================" << std::endl;
      if ( this->GetVerbose()>=4 || error ) {
	std::cout << "Press [ENTER] to continue." << std::endl;
	std::cin.get();
      }
    }

    if ( stepOK==true ) {
      // prep for next iteration
      m_iter++;
    }

  }//end of while loop
    
  for (int n=0; n<m_ndims; n++)
    solution[n] = m_x[n];
  if ( m_iter < m_max  // check if converged within iteration limit
       && stepOK==true ) // and last update was OK
    converged = true;
  
  if ( GetVerbose()>=1 ) {
    std::cout << "Solver stopped after " << m_iter << " of " << m_max << " iterations. ";
    std::cout << " status converged = " << converged << std::endl;
    std::cout << " |F(x)| = " << Fnorm << std::endl;
    std::cout << "{";
    for (int n=0; n<m_ndims; n++) std::cout << solution[n] << " ";
    std::cout << "}" << std::endl;
    std::cout << "////////////////////////////////////////////////////////////////////////" << std::endl;
  }
  return converged;
}

bool IterativeSolver::GetNextSolution( double* x, double* F, double* x_next ) {
  if ( fMethod==kNewton ) {
    return GetNextSolutionNewton( x, F, x_next );
  }
  else if ( fMethod==kJacobi ) {
    return GetNextSolutionJacobi( x, F, x_next );
  }
  else {
    std::cout << "IterativeSolver::GetNextSolution() -- ERROR. Unrecognized Iterative Method!" << std::endl;
    assert(false);
  }
}

bool IterativeSolver::DoWeStop( double* x, double* F, double* x_next, double* F_next ) {
  if ( fMethod==kNewton ) {
    return Newton_CheckTermination( x, F, x_next, F_next );
  }
  else if ( fMethod==kJacobi ) {
    return Jacobi_CheckTermination( x, F, x_next, F_next );
  }
  else {
    std::cout << "IterativeSolver::DoWeStop() -- ERROR. Unrecognized Iterative Method!" << std::endl;
    assert(false);
  }
}

void IterativeSolver::PrintF() {
  std::cout << "F: { ";
  for (int k=0; k<m_ndims; k++) {
    std::cout << m_F[k];
    if ( k!=m_ndims-1) std::cout << ", ";
  }
  std::cout << " }" << std::endl;
}

void IterativeSolver::PrintX() {
  std::cout << "X: { ";
  for (int k=0; k<m_ndims; k++) {
    std::cout << m_x[k];
    if ( k!=m_ndims-1) std::cout << ", ";
  }
  std::cout << " }" << std::endl;
}


// ====================================================================================================
// NEWTON METHODS

bool IterativeSolver::GetNextSolutionNewton( double* x, double* F, double* x_next ) {
  
  double delx[m_ndims]; // size of newton step
  
  BuildJacobian( x, m_J ); // Goes out to User's implementation
  if ( fNoInverse==false )
    Newton_InvertJacobian();

  bool stepOK = Newton_GetStep( F, delx ); // Implicitly uses Jacobian and Inverse of Jacobian
  for (int n=0; n<m_ndims; n++) {
    x_next[n] = m_x[n] - delx[n];
  }//end of n loop
  
  return stepOK;

}

bool IterativeSolver::Newton_CheckTermination( double* x, double* F, double* x_next, double* F_next ) {
  // looks like x and F are cruft, we have two copies because jacobi termination check needs both

  // calculate if we should terminate iteration
  double error = 0.;
  double Fnorm = 0.;
  double x_next_norm = 0.;
  double delx_norm = 0.;
  double delx[m_ndims];
  
  // calculate next step 
  bool stepOK = Newton_GetStep( F_next, delx ); 
  for (int n=0; n<m_ndims; n++) {
    x_next_norm += x_next[n]*x_next[n];
    delx_norm += delx[n]*delx[n];
    Fnorm += F_next[n]*F_next[n];
  }
  x_next_norm = sqrt(x_next_norm);
  delx_norm = sqrt(delx_norm);
  error = delx_norm/x_next_norm;
  Fnorm = sqrt(Fnorm);
  
  if (GetVerbose()>=2) {
    std::cout << "IterativeSolver::Newton_CheckTermination at iteration " << m_iter << std::endl;
    std::cout << "  Is x converging: |dx_next|/|x|=" << error << " < " << m_tolerance << " (" << (error<=m_tolerance)  << ") " << std::endl;
    std::cout << "  Is |F(x;k+1)| below tolerance: |F(x;k+1)|=" << Fnorm << " < " << m_tolerance << " (" << (Fnorm <=m_tolerance) << ")" << std::endl;
    std::cout << "  Was the last step OK: " << stepOK << std::endl;
    std::cout << "  |x_next|=" << x_next_norm
	      << ", |dx|=" << delx_norm
	      << std::endl;
    if ( fNoInverse )
      std::cout << "  F2 = " << Newton_CalculateFSquared( F_next ) << std::endl;
  }//if verbose

  // Determine termination
  if (error<=m_tolerance && Fnorm<=m_tolerance)
    return true;
  if (stepOK==false)
    return true;
  if ( Fnorm==0.0 )
    return true; // found solution

  return false;

}

void IterativeSolver::Newton_InvertJacobian() {
  // move the jacobian matrix into the m_Jrow array
  for (int row=0; row<m_ndims; row++) { 
    for (int col=0; col<m_ndims; col++) {
      (*m_JMatrix)[row][col] = m_J[row][col];
    }
  }

  LU->SetMatrix( (*m_JMatrix) );
  LU->Invert( *m_invJMatrix );

  for (int row=0; row<m_ndims; row++) { 
    for (int col=0; col<m_ndims; col++) {
      m_invJ[row][col] = (*m_invJMatrix)[row][col];
    }
  }
  //std::cout << "INVERTED" << std::endl;
  //PrintJacobian();
}

void IterativeSolver::Newton_PrintJacobian() {
  std::cout << "-------------------------------------------------------" << std::endl;
  std::cout << "Jacobian" << std::endl;
  std::cout << "{ ";
  for (int l=0; l<m_ndims; l++) {
    std::cout << "{ ";
    for (int k=0; k<m_ndims; k++) {
      std::cout << m_J[l][k];
      if ( k!=m_ndims-1 ) std::cout <<  ", ";
    }
    std::cout << "}";
    if ( l!=m_ndims-1 ) std::cout << ", ";
    std::cout << std::endl;
  }
  std::cout << "}" << std::endl;
  std::cout << "-------------------------------------------------------" << std::endl;
}

void IterativeSolver::Newton_PrintInverseJacobian() {
  std::cout << "-------------------------------------------------------" << std::endl;
  std::cout << "Inverse Jacobian" << std::endl;
  std::cout << "{ ";
  for (int l=0; l<m_ndims; l++) {
    std::cout << "{ ";
    for (int k=0; k<m_ndims; k++) {
      std::cout << m_invJ[l][k];
      if ( k!=m_ndims-1 ) std::cout <<  ", ";
    }
    std::cout << "}";
    if ( l!=m_ndims-1 ) std::cout << ", ";
    std::cout << std::endl;
  }
  std::cout << "}" << std::endl;
  std::cout << "-------------------------------------------------------" << std::endl;
}


bool IterativeSolver::Newton_GetStep( double* F, double* dx ) {

  bool stepOK = true;

  double dx_norm = 0.;
  if ( fNoInverse==false )
    dx_norm = Newton_CalculateBasicStep( F, dx );
  else 
    dx_norm = Newton_CalculateFSquaredStep( F, dx );
  
  if ( dx_norm!=dx_norm ) stepOK = false; // check for NAN
  if ( dx_norm<0 ) stepOK = false; // check for -inf

  // if regular newton method, return this step
  if ( fDamped==false ) return stepOK;

  // if F already 0? If so, we are done.
  double Fnorm = 0.;
  for (int n=0; n<m_ndims; n++)
    Fnorm += F[n]*F[n];
  if (Fnorm==0.0)
    return stepOK;

  // if damped, we adjust the step so that:
  //  |dx(k+1)| <= (1 - 0.5lambda)|dx(k)|
  //  
  double F_next[m_ndims];
  double x_next[m_ndims];
  
  double lambda = 2.; // set so that first iteration of next loop is 1.0
  double dxd[m_ndims]; // damped step
  double dxd_norm = 2.*dx_norm;
  
  do {
    lambda *= 0.5;
    if (lambda<1.0e-10)  {
      std::cout << "IterativeSolver::Newton_GetStep -- WARNING: could not dampen step:  |dx(k)|=" << dx_norm << "  |dx(k+1)|=" << dxd_norm << std::endl;
      return false;
    }
    // calculate next x with damping factor
    for (int n=0;n<m_ndims; n++) x_next[n] = m_x[n] - lambda*dx[n]; // calculate the next iteration
    // get the next value of F (goes out to user)
    CalculateF( x_next, F_next ); 
    // get the size of the next step
    dxd_norm = 0.;
    if ( fNoInverse==false )
      dxd_norm = Newton_CalculateBasicStep( F_next, dxd );
    else
      dxd_norm = Newton_CalculateFSquaredStep( F_next, dxd );
    // increase amount of damping if next step is larger than current step
    // so what it does is shorten the step size, until the next step will be smaller (forces step sizes to decrease monotonically)
    //if (IsSolutionValid( x_next )==false) std::cout << "invalid solution, shorten step" << std::endl;
    //} while ( dxd_norm >= (1-0.5*lambda)*dx_norm ); // 0.5 from quad rate of convergance
  } while ( dxd_norm >= dx_norm ); // looser
  
  for (int n=0; n<m_ndims; n++) dx[n] *= lambda; // return the damped step if we got out of the loop alive
  if (GetVerbose()>1) {
    std::cout << "IterativeSolver::Newton_GetStep() -- damped step returned, with damping factor " << lambda << ",  |dx(k)|=" << dx_norm << "  |dx(k+1)|=" << dxd_norm << " stepOK=" << stepOK << std::endl;
  }

  return stepOK;
}


double IterativeSolver::Newton_CalculateBasicStep( double* F, double* dx ) {
  /// !!!! This is an internal function
  /// !!!! note we've assumed:
  // (1) the jacobian has been inverted (its reference implicitly)
  // (2) that F(x) has been calculated for the current iteration
  //
  double dx_norm= 0.;
  for (int n=0; n<m_ndims; n++) {
    dx[n] = 0.;
    for (int k=0; k<m_ndims; k++) {
      dx[n] += m_invJ[n][k]*F[k]; 
    }
    dx_norm += dx[n]*dx[n];
  }//end of n loop
  dx_norm = sqrt(dx_norm);
  return dx_norm;
}

double IterativeSolver::Newton_CalculateFSquared( double* F ) {
  /// !!!! This is an internal function
  /// !!!! note we've assumed that F(x) has been calculated for the current iteration

  if ( m_theta==NULL ) { 
    m_theta = new double[m_ndims];
  }
  for (int n=0; n<m_ndims; n++) m_theta[n] = 1.0;
  
  double F2 = 0.;
  for (int n=0; n<m_ndims; n++) {
    F2 += sqrt(F[n]*F[n] + m_theta[n]*m_theta[n]) - m_theta[n];
  }
  return F2;
}

double IterativeSolver::Newton_CalculateFSquaredStep( double* F, double* dx ) {
  /// !!!! This is an internal function
  /// !!!! note we've assumed:
  // (1) the jacobian has been calculated

  // Calculate F-squared
  double F2 = Newton_CalculateFSquared( F );
  // implicitly, theta[m] has been defined

  // Calculate the gradient and the direction of deepest descent
  // Here we use the Jacobian
  double gradF2[m_ndims];
  double grad_norm = 0.;
  for (int n=0; n<m_ndims; n++) {
    gradF2[n] = 0.;
    for (int m=0; m<m_ndims; m++) {
      gradF2[n] += (F[m]/sqrt(F[m]*F[m] + m_theta[m]*m_theta[m]))*m_J[m][n];
    }
    grad_norm += gradF2[n]*gradF2[n];
  }
  grad_norm = sqrt(grad_norm);

  // Calculate Step
  double dx_norm = 0.;
  for (int n=0; n<m_ndims; n++) {
    dx[n] = F2*gradF2[n]/(grad_norm*grad_norm); //type2
    dx_norm += dx[n]*dx[n];
  }
  dx_norm = sqrt(dx_norm);
  
  return dx_norm;
}


// =============================================================================================================
// JACOBI METHODS

bool IterativeSolver::GetNextSolutionJacobi( double* x, double* F, double* x_next ) {
  // Remember by def:
  // 0 = F(x) = G(x) + x
  // We want to iterate x = G(x) until convergence
  double reducing_factor = 1.0;
  while ( reducing_factor>=1.0e-6 ) {

    double x_step_norm = 0.;
    double x_norm = 0.;
    double a = reducing_factor;
    for (int n=0; n<m_ndims; n++) {
      x_next[n] = a*F[n]+x[n]*(1-a);
      x_norm += x[n]*x[n];
      x_step_norm += (x_next[n]-x[n])*(x_next[n]-x[n]);
    }
    double relative_change = sqrt( x_step_norm/x_norm );

    if ( relative_change > m_a_steplimitfactor ) {
      reducing_factor *= 0.5;
    }
    else {
      break;
    }
  }
  std::cout << " jacobi solution made using reducing factor: " << reducing_factor << std::endl;
  return true;
}

bool IterativeSolver::Jacobi_CheckTermination( double* x, double* F, double* x_next, double* F_next ) {
  // For Jacobi-Iterator, check to see if converging

  double delF = 0.;
  double delX = 0.;
  double xnorm = 0.;
  double Fnorm = 0.;
  for (int n=0; n<m_ndims; n++) {
    delF += (F[n]-F_next[n])*(F[n]-F_next[n]);
    delX += (x[n]-x_next[n])*(x[n]-x_next[n]);
    xnorm += x[n]*x[n];
    Fnorm += F[n]*F[n];
  }
  delF = sqrt(delF);
  delX = sqrt(delX);
  Fnorm = sqrt(Fnorm);
  xnorm = sqrt(xnorm);
	       
  if ( GetVerbose()>1 ) {
    std::cout << "IterativeSolver::Jacobi_CheckTermination at iteration " << m_iter << std::endl;
    std::cout << "  Is x converging: |dx|/|x| = " << delX/xnorm << " < " << m_tolerance << " (" << ((delX/xnorm)<m_tolerance) << ")" << std::endl;
  }

  return false;
    
  if ( delX/xnorm < m_tolerance  )
    return true;
  if ( Fnorm==0.0 ) return true;

  return false;
}

void IterativeSolver::BuildJacobian( double* x, double** m_J ) {
  if ( fMethod==kNewton ) {
    std::cout << "Empty BuildJacobian called when running in Newton Minimization mode" << std::endl;
    std::cout << "You, the user, have to reimplement this." << std::endl;
    assert(false);
  }
}

// =============================================================================================================
