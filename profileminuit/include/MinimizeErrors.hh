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
 * -------------------------------------------------------------------------------------------
 *  \class MinimizeErrors
 *  \ingroup ProfileMinuit
 *  \brief Concrete Implementation of IterativeSolver tailored for Pull Term Minimization
 *
 *  Expects to be passed a set of N bins for both Observed (data) and
 *  Expected (MC) bins.
 *  Also K errors must be defined. For each error, a set of N Fij values
 *  must be passed to it.  Fij values is the change in the expected (MC) bin for a 1-sigma
 *  change in the systematic error value. Units of events/sigma.
 *
 * Extra tools:
 * Can define a sequential set of bins as a "sample". This allows one to ask for
 * the chi-squared value of those grouped bins.
 * Also, one may pass a TFile to which a TTree with the chi-squared values can be appended.
 *
 * To dos:
 * Allow easy overloading for the Minimization Function? Might be redundant to just creating
 * a new concrete implementation of IterativeSolver.
 *
 *   Verbose levels:
 *    0: Quiet (at least base class will be)
 *    1: Print out summary of fit
 *    2: Print out more information about each step
 *    3: Print out maximum amount of info
 *    4: Pause after each interation
 * -------------------------------------------------------------------------------------------
*/

#ifndef __MinimizeErrors__
#define __MinimizeErrors__

#include "IterativeSolver.hh"
#include <string>
#include <vector>
#include <map>

#include "TMatrixDSym.h"

class TFile;
class TTree;

namespace qosc {

  class SamplePartition;

  class MinimizeErrors : public IterativeSolver {

  public:

    MinimizeErrors();
    MinimizeErrors( int nerrors, int nbins, IterativeSolver::Method meth=IterativeSolver::kNewton );
    virtual ~MinimizeErrors();

    virtual bool RunIteration(double* solution);

    void DefineTrees( TFile* m_outputfile );
    void FillResult();
    void WriteResults();

    void SetFij( int errnum, int binnum, double value, std::string errorlabel="", std::string binlabel="" );
    void SetBinNObserved( int binnum, double observed );
    double GetBinNObserved( int binnum );
    void SetBinNExpected( int binnum, double expected );
    double GetBinNExpected( int binnum );
    void ImportCovarianceMatrix( TMatrixDSym* cov );
    void SetErrorCovariance( int errnum1, int errnum2, double cov ); 
    void SetBinLabel(int binnum, std::string binlabel );
    void SetErrorLabel( int errnum, std::string errorlabel );

    double CalculatePull(int n_bin, double* solution);
    double CalculatePartialPull(int n_bin, int err, double* solution);
    double CalculateDoublePartialPull(int n_bin, int err1, int err2, double* solution);
    double CalculateChiSquared(); /// Calculate chi-squared assuming zero-pull value
    double CalculateChiSquared( double* epsilons, double& error ); /// If given null, will used the stored value of x
    double CalculateChiSquared( double* epsilons, double& error, std::string sampleName ); /// If given null, will used the stored value of x
    double CalculatePullChiSquared( double* epsilons=NULL );

    double GetFractionalFij( int bin, int err ) { 
      if ( m_NExpected[bin] ) return m_Fij[bin][err];
      else return 0;
    };
    //double GetFij( int bin, int err ) { return m_Fij[bin][err]*m_NExpected[bin]; };
    double GetFij( int bin, int err ) { return m_Fij[bin][err]; };
    void PrintFij();
    void WriteFij();
  
    void PrintErrorList();
    void PrintErrorInfo( int err );
    void PrintBinList();

    void GetLinearizedSolution( double* solution );
    virtual void MakeValidSolution( double* solution );
    virtual bool IsSolutionValid( double* solution );
    void SetPullMin(int, double);
    void SetPullValue( int errid, double value );
    double GetPullValue( int errid ) { return IterativeSolver::m_x[errid]; };
    void DefineSample( std::string name, int firstbin, int lastbin );
    void FixParameter( int errnum ) { m_minimize_errorterm[errnum] = false; };
    void VaryParameter( int errnum ) { m_minimize_errorterm[errnum] = true; };
    void TermDoesNotPull( int errnum ) { m_term_pulls[errnum] = false; };
    void PenalizeParameter( int errnum, bool penalize ) { m_penalize_errorterm[errnum] = penalize; };

  protected:

    virtual void CalculateF( double* x, double* F );
    virtual void BuildJacobian( double* x, double** J );

    int m_nerrors;
    int m_nbins;

    double** m_Fij; /// array to store Fij
    double* m_NObserved;
    double* m_NExpected;
    double* m_binChiSquared;
    double* m_pull_min;
    // ------------------------------------------------
    // Important flags for each parameter. by default, all are set to true for each parameter.
    bool* m_minimize_errorterm; ///< flag to determine if chi-squared is minized with respect to this parameter
    bool* m_penalize_errorterm; ///< flag to determine if parameter included in penalty term
    bool* m_term_pulls; ///< flag that determines if a term modifies the expectation
    // ------------------------------------------------

    char** m_error_list;
    char** m_binlabel_list;

    TFile* m_outputfile;
    TTree* m_tree_fitresults;
    TTree* m_tree_errorterms;
    TTree* m_tree_binterms;

    // ----------------------------------------------------
    // Optimizations
    // (1) Pull Cache
    bool fCacheNeedsToBeRemade;
    double* m_cached_pulls; // [ nbins ]
    double* m_cached_pulls_x_source; // [ nerrors ]
    void MakePullCache( double* solution );
    //   // (2) Partial derivative cahce
    //   bool fPartialCacheNeedsToBeRemade;
    //   double* m_cached_partial_pulls; // [ nbins ][ nerrors ]
    //   void MakePartialPullCache();
    //   // (3) Double Partial derivative cache
    //   bool fDoublePartialCacheNeedsToBeRemade;
    //   double* m_cached_double_partial_pulls; // [ nbins ][ nerrors ][ nerrors ]
    //   void MakeDoublePartialPullCache();
    //   // ----------------------------------------------------

    double ChiSquared;
    double errChiSquared;
    int converged_flag;
    int minimized_flag;

    TMatrixD* linM;
    TMatrixD* linInvM;
    TDecompLU* linLU;
    TMatrixDSym* m_covariance;
    TMatrixDSym* m_inv_cov; ///< inverse of covariance matrix
    bool fNeedToInvertCov; 
    void InvertCovarianceMatrix();

    std::map< std::string, SamplePartition*  >  m_samplePartitions;

    // second order routines
  public:
    void SetF2ij( int binnum, double* f2array );
    void UseF2ij( bool use=true ) { fUseF2ij = use; }; // second order dependence on epsilons
  protected:
    double** m_F2ij; /// array to store second-order Fij
    bool fUseF2ij;

  public:

    TMatrixDSym* GetCovarianceMatrix() { return m_covariance; };
    TMatrixDSym* GetInverseCovarianceMatrix() { return m_inv_cov; };

  };

}

#endif
