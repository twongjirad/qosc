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

#include "MinimizeErrors.hh"
#include <iostream>
#include <cmath>
#include <assert.h>
#include "TFile.h"
#include "TTree.h"
#include "SamplePartition.hh"

using namespace std;
using namespace qosc;

MinimizeErrors::MinimizeErrors() {
}

MinimizeErrors::MinimizeErrors( int nerrors, int nbins, IterativeSolver::Method meth  ) : 
  IterativeSolver( nerrors, meth ) {

  m_nerrors = nerrors;
  m_nbins = nbins;
  converged_flag = 0;
  minimized_flag = 0;
  m_samplePartitions.clear();
  
  // allocate Fij, NObserved and Nexpected
  m_Fij = new double*[m_nbins];
  for (int n=0; n<m_nbins; n++) {
    m_Fij[n] = new double[m_nerrors];
    memset( m_Fij[n], 0, sizeof(double)*m_nerrors  );
  }
  m_NObserved = new double[m_nbins];
  m_NExpected = new double[m_nbins];
  m_binChiSquared = new double[m_nbins];
  m_binlabel_list = new char*[m_nbins];
  for (int n=0; n<m_nbins; n++) {
    m_binlabel_list[n] = new char[100];
    sprintf( m_binlabel_list[n], "bin %d", n );
  }
  m_error_list = new char*[m_nerrors];
  m_minimize_errorterm = new bool[m_nerrors];
  m_term_pulls = new bool[m_nerrors];
  m_penalize_errorterm = new bool[m_nerrors];
  for (int m=0; m<m_nerrors; m++) {
    m_error_list[m] = new char[200];
    m_minimize_errorterm[m] = true;
    m_penalize_errorterm[m] = true;
    m_term_pulls[m] = true;
    sprintf( m_error_list[m], "error %d", m );
  }
  m_pull_min = new double[m_nerrors];
  for (int k=0; k<m_nerrors; k++) m_pull_min[k] = 0.;

  // Allocate the Matrices
  linM = new TMatrixD( m_nerrors, m_nerrors );
  linInvM = new TMatrixD( m_nerrors, m_nerrors );
  linLU = new TDecompLU( m_nerrors );
  m_covariance = new TMatrixDSym( m_nerrors );
  // fill covariance with identity matrix by default
  for (int n=0; n<m_nerrors; n++) (*m_covariance)[n][n] = 1.0;
  for (int i=0; i<m_nerrors; i++) {
    for (int j=0; j<m_nerrors; j++) {
      if ( i==j ) continue;
      (*m_covariance)[i][j] = 0.;
    }
  }
  m_inv_cov = new TMatrixDSym( *m_covariance );
  fNeedToInvertCov = true;

  // second-order matrix
  m_F2ij = new double*[m_nbins];
  for (int n=0; n<m_nbins; n++) {
    m_F2ij[n] = new double[m_nerrors*m_nerrors];
    memset( m_F2ij[n], 0, sizeof(double)*m_nerrors*m_nerrors );
  }

  // Pull Cache
  m_cached_pulls = new double[m_nbins];
  m_cached_pulls_x_source = new double[m_nerrors];
  memset( m_cached_pulls, 1, sizeof(double)*m_nbins );
  memset( m_cached_pulls_x_source, -10, sizeof(double)*m_nerrors );

  m_verbose = 0;
  fCacheNeedsToBeRemade = true;
}

MinimizeErrors::~MinimizeErrors() {
  for (int n=0; n<m_nbins; n++) {
    delete [] m_Fij[n];
    delete [] m_F2ij[n];
  }
  for (int n=0; n<m_nerrors; n++)
    delete [] m_error_list[n];

  delete [] m_Fij;
  delete [] m_error_list;
  delete [] m_minimize_errorterm;
  delete [] m_penalize_errorterm;

  delete [] m_NObserved;
  delete [] m_NExpected;
  for (int n=0; n<m_nbins; n++)
    delete [] m_binlabel_list[n];
  delete [] m_binlabel_list;

}

void MinimizeErrors::SetBinLabel( int binnum, std::string binlabel ) {
  if (binlabel!="") sprintf(m_binlabel_list[binnum], binlabel.c_str());
}

void MinimizeErrors::SetErrorLabel( int errnum, std::string errorlabel ) {
  if (errorlabel!="") sprintf( m_error_list[errnum], errorlabel.c_str() );
}

void MinimizeErrors::SetFij( int errnum, int binnum, double value, std::string errorlabel, std::string binlabel ) {
  // Set the Fij values
  if (binnum<0 || binnum>=m_nbins ) {
    std::cout << "Bin number is out of range: " << binnum << "(<0 or >=" << m_nbins << ")" << std::endl;
    assert(false);
  }
  else if ( errnum<0 || errnum>=m_nerrors ) {
    std::cout << "Error number is out of range: " << errnum << "(<0 or >=" << m_nerrors << ")" << std::endl;
    assert(false);
  }
  if ( m_term_pulls[errnum] ) {
    m_Fij[binnum][errnum] = value;
  }
  else {
    m_Fij[binnum][errnum] = 0.0;
  }

  SetBinLabel( binnum, binlabel );
  SetErrorLabel( errnum, errorlabel );
}

void MinimizeErrors::SetF2ij( int binnum, double* f2ijarray ) {
  // Set the Fij values
  if (binnum<0 || binnum>=m_nbins ) {
    std::cout << "Bin number is out of range: " << binnum << "(<0 or >=" << m_nbins << ")" << std::endl;
    assert(false);
  }
  memcpy( m_F2ij[binnum], f2ijarray, sizeof(double)*m_nerrors*m_nerrors );
}


void MinimizeErrors::SetBinNObserved( int binnum, double observed ) {
  m_NObserved[binnum] = observed;
  // reset the flags
  minimized_flag = 0;
  converged_flag = 0;
}

double MinimizeErrors::GetBinNObserved( int binNum ) {
  return m_NObserved[binNum];
}

void MinimizeErrors::SetBinNExpected( int binnum, double expected ) {
  m_NExpected[binnum] = expected;
  // reset the flags
  minimized_flag = 0;
  converged_flag = 0;
}

double MinimizeErrors::GetBinNExpected( int binnum ) {
  return m_NExpected[binnum];
}

void MinimizeErrors::SetPullMin( int errnum, double pull_value_at_min ) {
  if ( errnum>=0 && errnum<=m_nerrors )
    m_pull_min[errnum] = pull_value_at_min;
}

void MinimizeErrors::SetPullValue( int errnum, double value ) {
  if ( errnum>=0 && errnum<=m_nerrors )
    m_x[errnum] = value;
}

void MinimizeErrors::CalculateF( double* x, double* F ) {

  bool error = false;
  if ( fNeedToInvertCov ) InvertCovarianceMatrix();

  for (int n=0; n<m_nerrors; n++) {
    F[n] = 0.;
    if ( !m_minimize_errorterm[n] ) continue; // if the parameter is held fixed for the pull method, the partial of the chi-squared is held at zero.

    // penalty portion: sum over all error values that pull
    for (int j=0; j<m_nerrors; j++)
      if ( m_term_pulls[j] )
	F[n] += (x[j]-m_pull_min[j])*(*m_inv_cov)[j][n] + (*m_inv_cov)[n][j]*(x[j]-m_pull_min[j]);

    // bin difference part (modified residual)
    for (int b=0; b<m_nbins; b++) {
      double bin_pull = CalculatePull( b, x );
      double dpull = CalculatePartialPull( b, n, x );
//       if ( m_NObserved[b]>0 ) 
// 	F[n] += 2.0*( m_Fij[b][n] - (m_NObserved[b]*m_Fij[b][n])/(bin_pull*m_NExpected[b])); // assumed sigma = 1.0
      //F[n] += 2.0*( m_NExpected[b]*dpull - m_NObserved[b]*(dpull/bin_pull) );
      F[n] += 2.0*( m_NExpected[b]*dpull );
      if ( m_NObserved[b]>0 ) {
	F[n] -= 2.0*m_NObserved[b]*(dpull/bin_pull);
      }
    }
    if ( F[n]!=F[n] ) {
      std::cout << "F[" << n << "] is NaN: E=" << m_NExpected[n] << " " << m_NObserved[n] << std::endl;
      error = true;
    }
  }

  if ( error ) {
    std::cout << "x= {";
    for (int i=0; i<m_nerrors; i++) {
      std::cout << x[i] << ", ";
      if ( m_minimize_errorterm[i] )
	m_x[i] = x[i];
    }
    std::cout << "}" << std::endl;
    PrintBinList();
    PrintErrorList();
    PrintFij();
    std::cout << "MinimizeErrors stopped because CalculteF calculated a NaN. Iteration #" << m_iter << std::endl;
    std::cin.get();
  }
}


void MinimizeErrors::BuildJacobian( double* x, double** J ) {

  // This is the second derivative of the chi-squared. Also the first derivative of F(n)
  if ( fNeedToInvertCov ) InvertCovarianceMatrix();

  double binpull = 0.;
  for (int l=0; l<m_nerrors; l++) {
    for (int k=0; k<m_nerrors; k++) {

      J[l][k] = 0.0;
      if ( !m_minimize_errorterm[k] || !m_minimize_errorterm[l] ) {
	if ( l==k ) J[l][k] = 1.0; // so i can invert this. this is a hack though.
	continue;
      }

      J[l][k] += (*m_inv_cov)[l][k] + (*m_inv_cov)[k][l];
      
      for (int n=0; n<m_nbins; n++) {
	binpull = CalculatePull( n, x );
	double dpull1 = CalculatePartialPull( n, l, x );
	double dpull2 = CalculatePartialPull( n, k, x );
	double ddpull = CalculateDoublePartialPull( n, l, k, x );
	//J[l][k] += 2.0*(m_NExpected[n]*ddpull - m_NObserved[n]*((binpull*ddpull - dpull1*dpull2)/pow( binpull, 2 ) ));
	J[l][k] += 2.0*(m_NExpected[n]*ddpull);
	if ( m_NObserved[n]> 0 )
	  J[l][k] -= 2.0*m_NObserved[n]*((binpull*ddpull - dpull1*dpull2)/pow( binpull, 2 ));
      }//end of n loop
    }//end of k
  }//end of l
  
}


void MinimizeErrors::PrintFij() {
  std::cout << "-------------------------------------------------------" << std::endl;
  std::cout << "Fij" << std::endl;
  for (int l=0; l<m_nerrors; l++) {
    std::cout << "Err " << l << " ";
    for (int k=0; k<m_nbins; k++) {
      std::cout << m_Fij[k][l] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << "-------------------------------------------------------" << std::endl;
}

#include "TH2D.h"

void MinimizeErrors::WriteFij() {
  TH2D fij("Fij", "Fij", m_nbins, 0, m_nbins, m_nerrors, 0, m_nerrors );
  for (int k=0; k<m_nbins; k++) {
    for (int l=0; l<m_nerrors; l++) {
      fij.SetBinContent( k+1, l+1, m_Fij[k][l] );
    }
  }
  fij.Write();
}

void MinimizeErrors::PrintErrorList() {
  std::cout << "--------------------------------------------------------------------------------" << std::endl;
  std::cout << " ERROR LIST: Total of " << m_nerrors << std::endl;
  for (int err=0; err<m_nerrors; err++)
    PrintErrorInfo( err );
  std::cout << "--------------------------------------------------------------------------------" << std::endl;
}

void MinimizeErrors::PrintErrorInfo( int err ) {
  std::cout << "   (" << err << ") " << m_error_list[err] << " = " << m_x[err] << " : central value=" << m_pull_min[err];
  if ( m_minimize_errorterm[err] )
    std::cout << " : Minimize Term ";
  if ( m_term_pulls[err] )
    std::cout << " : Term pulls ";
  if ( m_penalize_errorterm[err] )
    std::cout << " : Penalize Term";
  std::cout << std::endl;
}

void MinimizeErrors::PrintBinList() {
  std::cout << "--------------------------------------------------------------------------------" << std::endl;
  std::cout << " BIN LIST: Total of " << m_nbins << std::endl;
  double observed_total = 0.;
  double nominal_total = 0.;
  double pulled_total = 0.;
  for (int bin=0; bin<m_nbins; bin++) {
    double pull = CalculatePull(bin,m_x);
    double pulls[m_nerrors];
    memset(pulls,0,sizeof(double)*m_nerrors);
    for (int k=0; k<m_nerrors; k++) {
      if ( m_term_pulls[k] ) {
	//pull += m_Fij[bin][k]*m_x[k];
	//pulls[k] = m_Fij[bin][k]*m_x[k];
	pulls[k] = m_Fij[bin][k]*m_x[k];
      }
    }
    observed_total += m_NObserved[bin];
    nominal_total += m_NExpected[bin];
    pulled_total += m_NExpected[bin]*pull;
    std::cout << "   (" << bin << ") " << m_binlabel_list[bin] 
	      << " : Observed=" << m_NObserved[bin] 
	      << "  Expected=" << m_NExpected[bin] 
	      << "  Exp(pulled)=" << m_NExpected[bin]*pull
	      << "  Pull=" << pull 
	      << "  Diff=" << m_NExpected[bin]*pull-m_NObserved[bin]
	      << std::endl;
    std::cout << "  Fij[" << bin << "][k]*x = 1 + ( ";
    double pulltot = 1.;
    for (int k=0; k<m_nerrors; k++) {
      std::cout << pulls[k];
      if ( k<m_nerrors-1) std::cout << ", ";
      pulltot += pulls[k];
    }
    std::cout << ") = " << pulltot << std::endl;
  }
  std::cout << "--------------------------------------------------------------------------------" << std::endl;
  std::cout << " ( Total ) : Observed=" << observed_total 
	    << "  Nominal Expected=" << nominal_total
	    << "  Pulled Expected=" << pulled_total 
	    <<  std::endl;
  std::cout << "--------------------------------------------------------------------------------" << std::endl;
  std::cout << "Calculated using x={";
  for (int n=0; n<m_nerrors; n++) std::cout << m_x[n] << " ";
  std::cout << " }" << std::endl;
  std::cout << "--------------------------------------------------------------------------------" << std::endl;
}

void MinimizeErrors::DefineTrees( TFile* outputfile ) {

  m_outputfile = outputfile;
  m_outputfile->cd();

  char bin_def[100];
  m_tree_binterms = new TTree("BinTerms", "Fitter Bin Terms");
  m_tree_binterms->Branch("nBins", &m_nbins, "nBins/I");
  m_tree_binterms->Branch("BinTerm", bin_def, "BinTerm[100]/C");
  for (int n=0; n<m_nbins; n++) {
    sprintf( bin_def, m_binlabel_list[n] );
    m_tree_binterms->Fill();
  }
  m_tree_binterms->Write();

  char error_def[200];
  m_tree_errorterms = new TTree("ErrorTerms", "Fitter Error Terms");
  m_tree_errorterms->Branch("nErrors", &m_nerrors, "nErrors/I");
  m_tree_errorterms->Branch("ErrorTerm",error_def,"ErrorTerm[100]/C");
  for (int k=0; k<m_nerrors; k++) {
    sprintf( error_def, m_error_list[k] );
    m_tree_errorterms->Fill();
  }
  m_tree_errorterms->Write();

  m_tree_fitresults = new TTree("FitResults", "Chi-squared and pull values after each fit");
  m_tree_fitresults->Branch("ChiSquared", &ChiSquared,"ChiSquared/D");
  m_tree_fitresults->Branch("errChiSquared", &errChiSquared, "errChiSquared/D");
  m_tree_fitresults->Branch("minimized", &minimized_flag, "minimized/I");
  m_tree_fitresults->Branch("converged", &converged_flag, "converged/I");
  m_tree_fitresults->Branch("nBins", &m_nbins, "nBins/I" );
  m_tree_fitresults->Branch("nErrors", &m_nerrors, "nErrors/I" );
  // define the epsilon addresses
  for (int k=0; k<m_nerrors; k++) {
    char errdef[200];
    sprintf(errdef,"%s/D", m_error_list[k]);
    m_tree_fitresults->Branch( m_error_list[k], &m_x[k], errdef );
  }
  char bincs_def[200];
  sprintf( bincs_def, "binCS[%d]/D", m_nbins );
  //m_tree_fitresults->Branch( "binCS", m_binChiSquared, bincs_def );

}


void MinimizeErrors::FillResult() {
  m_tree_fitresults->Fill();
}

void MinimizeErrors::WriteResults() {
  m_tree_fitresults->Write();
}

bool MinimizeErrors::RunIteration(double* solution) {
  converged_flag = 0;

  // Run iterative method to find solution
  bool converged = IterativeSolver::RunIteration(solution);
  errChiSquared = 0.;
  double chi2 = CalculateChiSquared( solution, errChiSquared );
  // check that all the pulls are legal
  if ( converged && IsSolutionValid(solution)==false && fabs(errChiSquared)/fabs(chi2)>1.0e-4 ) {
    converged = false;
  }
  
  if (converged) converged_flag = 1;
  minimized_flag = 1;

  if ( GetVerbose()>=3 ) {
    std::cout << "Bin results: " << std::endl;
    for (int n=0; n<m_nbins; n++) {
      double mod_expectation = m_NExpected[n]*CalculatePull( n, solution );
      std::cout << "Bin [" << n << "]:"
		<< " Observed=" << m_NObserved[n] 
		<< " Expected=" << m_NExpected[n] 
		<< " E(1+fe)=" << mod_expectation 
		<< ", Res=" << mod_expectation-m_NObserved[n] 
		<< " pull=" << CalculatePull(n, solution)
		<< " bin chi2=" << m_binChiSquared[n] << std::endl;
    }
    std::cout << std::endl;
  }

  //std::cin.get();
  return converged;
}

void MinimizeErrors::GetLinearizedSolution( double* solution ) {
  
  for (int k=0; k<m_nerrors; k++) {
    for (int l=0; l<m_nerrors; l++) {
      (*linM)[k][l] = 0.;
      if (!m_minimize_errorterm[l] || !m_minimize_errorterm[k] ) {
	// zero out the solution for fixed error terms
	if ( l==k ) (*linM)[k][l] = 1.0;
	else (*linM)[k][l] = (*linM)[l][k] = 0.0;
	continue;
      }
      (*linM)[k][l] = (*m_inv_cov)[k][l];
      for (int n=0; n<m_nbins; n++) {
	if ( m_NExpected[n]>0 ) {
	  //(*linM)[k][l] += m_Fij[n][k]*m_Fij[n][l]*m_NObserved[n]/pow(m_NExpected[n],2);
	  (*linM)[k][l] += m_Fij[n][k]*m_Fij[n][l]*m_NObserved[n]; // if Fij are in fractional change per sigma
	  //(*linM)[k][l] += m_Fij[n][k]*m_Fij[n][l]*m_NObserved[n]*m_NExpected[n]/pow(m_NExpected[n],2); // in units of events per sigma
	}
      }
    }
  }
  
  linLU->SetMatrix( (*linM) );
  linLU->Invert( (*linInvM) );

  for (int k=0; k<m_nerrors; k++) {
    if ( !m_minimize_errorterm[k] ) continue; // leaves original value alone
    solution[k] = 0.;
    for (int l=0; l<m_nerrors; l++) {
      
      double v = 0.;

      for (int i=0; i<m_nerrors; i++) {
	//if ( !m_minimize_errorterm[i] ) continue;
	if ( m_term_pulls[i] )
	  v += (*m_inv_cov)[l][i]*m_pull_min[i];
      }
      
      for (int n=0; n<m_nbins; n++) {
	v += m_Fij[n][l]*(m_NObserved[n]-m_NExpected[n]); // if Fij are in fraction change
	//if ( m_NObserved[n] > 0 ) v += (m_NObserved[n]-m_NExpected[n])*(m_Fij[n][l]/m_NExpected[n]); // if Fij are in unit of events
      }
      solution[k] += (*linInvM)[k][l]*v;
    }
  }
  
  if (GetVerbose()>0) {
    std::cout << "MinimizeErrors::GetLinearizedSolution() x: { ";
    for (int k=0; k<m_nerrors; k++) {
      std::cout << solution[k] << ", ";
    }
    std::cout << "}" << std::endl;
  }

}

bool MinimizeErrors::IsSolutionValid( double* solution ) {
  double pull = 1.0;
  bool pullOK = false;
  while (pullOK==false) {
    pullOK = true;
    for (int n=0; n<m_nbins; n++) {
      pull = CalculatePull( n, solution );
      if (pull*m_NExpected[n]<0) {
	std::cout << "pull*En for bin #" << n << " is negative: pull=" << pull << " En=" << m_NExpected[n] << std::endl;
	return false;
      }
    }
  }
  return true;
}

void MinimizeErrors::MakeValidSolution( double* solution ) {

  // calculate pull[n], look for negative, don't start at a crap position
  int numhalves = 0;
  bool pullOK = false;
  while (pullOK==false && numhalves<20) {
    pullOK = IsSolutionValid(solution);
    if (pullOK==false) {
      // if pull bad... half solution
      if (GetVerbose()>0) std::cout << "MinimizeErrors::MakeValidSolution() bad initial solution, halving active parameters..." << std::endl;
      for (int k=0; k<m_nerrors; k++) {
	if ( !m_minimize_errorterm[k] ) continue;
	solution[k] *= 0.5;
      }
      numhalves++;
    }
  }
    
  if (GetVerbose()>0 && numhalves>=1) {
    std::cout << "Modifed x to make it a valid value, numhalves=" << numhalves << ": { ";
    for (int k=0; k<m_nerrors; k++) {
      std::cout << solution[k] << ", ";
    }
    std::cout << "}" << std::endl;
    if (GetVerbose()>2) {
      std::cout << "Enter to continue" << std::endl;
      std::cin.get();
    }
  }
}

double MinimizeErrors::CalculatePull(int n_bin, double* solution) {

  if ( solution==NULL ) solution = m_x;

  // Try using a cache. 
  // Check to see if pull values have changed.
  int ndiff = 0;
  ndiff = memcmp( m_cached_pulls_x_source, solution, sizeof(double)*m_nerrors );
  
  if ( ndiff!=0 ) {
    MakePullCache( solution );
  }
  
  return m_cached_pulls[n_bin];
}

double MinimizeErrors::CalculatePartialPull(int n_bin, int n_err,  double* solution) {
  if ( solution==NULL ) solution = m_x;
  if ( m_NExpected[n_bin]==0 ) return 0.0;  
  if ( !m_term_pulls[n_err] ) return 0.;
  //double dpull = (m_Fij[n_bin][n_err]/m_NExpected[n_bin]); // dividing by expectation because Fij given in number of events, not fraction.
  double dpull = m_Fij[n_bin][n_err]; // dividing by expectation because Fij given in number of events, not fraction.
  if ( fUseF2ij ) {
    for (int j=0; j<m_nerrors; j++) {
      if ( j!=n_err && m_term_pulls[j] )
	dpull += (*(m_F2ij[n_bin] + m_nerrors*n_err + j )/m_NExpected[n_bin])*solution[j];
    }
  }
  return dpull;
}

double MinimizeErrors::CalculateDoublePartialPull(int n_bin, int err1, int err2, double* solution) {
  if ( solution==NULL ) solution = m_x;
  if ( m_NExpected[n_bin]==0 ) return 0.0;
  if ( !m_minimize_errorterm[err1] || !m_minimize_errorterm[err2] ) return 0;
  if ( !fUseF2ij ) return 0;
  
  double ddpull = *(m_F2ij[n_bin] + m_nerrors*err1 + err2 )/m_NExpected[n_bin];
  return ddpull;
}

void MinimizeErrors::DefineSample( std::string name, int firstbin, int lastbin ) {

  // Note last bin is inclusive!
  std::map< std::string, SamplePartition* >::iterator it = m_samplePartitions.find( name );
  SamplePartition* sample = NULL;
  if ( it==m_samplePartitions.end() ) {
    sample = new SamplePartition();
    char branchname[200];
    sprintf(branchname, "%s_chi2",name.c_str());
    m_tree_fitresults->Branch( branchname, &sample->SampleChiSquared, std::string( std::string(branchname)+"/D").c_str() );
  }
  else {
    sample = (*it).second;
  }

  sample->name = name;
  sample->firstbin = firstbin;
  sample->lastbin = lastbin;
  sample->userange = true;

  m_samplePartitions[name] = sample;
  
}

// ------------------------------------------------------------------------------------------
// Pull Cache Functions

void MinimizeErrors::MakePullCache( double* solution ) {
//   if ( !m_cached_pulls ) {
//     m_cached_pulls = new double[m_nbins];    
//     memset( m_cached_pulls, 0, sizeof(double)*m_nbins );
//     m_cached_pulls_x_source = new double[m_nerrors];
//     memset( m_cached_pulls_x_source, 0, sizeof(double)*m_nerrors );
//   }

  memcpy( m_cached_pulls_x_source, solution, sizeof(double)*m_nerrors );

  for (int n=0; n<m_nbins; n++ ) {
    m_cached_pulls[n] = 1.0;
    if ( m_NExpected[n]==0 ) continue;

    for (int k=0;k<m_nerrors; k++) {
      //double epsilon = solution[k]-m_pull_min[err]
      //if ( !m_minimize_errorterm[k] ) continue; // if newton minmizer not minimizing this term
      if ( !m_term_pulls[k] ) continue; // if newton minmizer not minimizing this term
      //m_cached_pulls[n] += (m_Fij[n][k]/m_NExpected[n])*solution[k]; // if F is in change in bin
      m_cached_pulls[n] += m_Fij[n][k]*solution[k]; // if F is in fraction change
//       if ( fUseF2ij ) {
// 	for (int j=k+1; j<m_nerrors; j++) {
// 	  //if ( !m_minimize_errorterm[j] ) continue;
// 	  if ( !m_term_pulls[k] ) continue; // if newton minmizer not minimizing this term
// 	  m_cached_pulls[n] += *(m_F2ij[n] + m_nerrors*k + j )/m_NExpected[n]*solution[k]*solution[j];
// 	}
//       }
    }//end of error sum
    if ( m_cached_pulls[n]<0.0 ) m_cached_pulls[n] = 1.0e-10;
  } // end of bin loop

  fCacheNeedsToBeRemade = false;
  
}

// ------------------------------------------------------------------------------------------

void MinimizeErrors::ImportCovarianceMatrix( TMatrixDSym* cov ) {
  for (int i=0; i<m_nerrors; i++) {
    for (int j=i; j<m_nerrors; j++) {
      (*m_covariance)[i][j] = (*m_covariance)[j][i] = (*cov)[i][j];
    }
  }
  InvertCovarianceMatrix();
}

void MinimizeErrors::SetErrorCovariance( int errnum1, int errnum2, double cov ) {
  (*m_covariance)[errnum1][errnum2] = (*m_covariance)[errnum2][errnum1] = cov;
  fNeedToInvertCov = true;
}

void MinimizeErrors::InvertCovarianceMatrix() {
  if (m_inv_cov) delete m_inv_cov;
  m_inv_cov = new TMatrixDSym( *m_covariance );
  m_inv_cov->Invert();
  fNeedToInvertCov = false;
}

double MinimizeErrors::CalculateChiSquared() {
  double error = 0;
  ChiSquared = CalculateChiSquared( NULL, error );
  return ChiSquared;
}

double MinimizeErrors::CalculateChiSquared( double* epsilons, double& error ) {
  ChiSquared = CalculateChiSquared( epsilons, error, "__all__" );
  std::map< std::string, SamplePartition* >::iterator it;
  for (it=m_samplePartitions.begin(); it!=m_samplePartitions.end(); it++) {
    SamplePartition* samplepart = (*it).second;
    double sample_err = 0.;
    samplepart->SampleChiSquared = CalculateChiSquared( epsilons, sample_err, (*it).first );
  }
  return ChiSquared;
}


double MinimizeErrors::CalculateChiSquared( double* epsilons, double& error, std::string sampleName ) {

  // sometimes there are little tweaks that must be made to keep chi-squared from becoming NAN
  // the error introduced by these tweaks are stored in the error value
  error = 0.;
  
  if ( epsilons==NULL ) {
    epsilons = m_x;
  }

  double binpull = 0.0;
  double chi2 = 0;
  error = 0.;

  // A Sample Partion is essentially a subset of bins that have a name attached to it.
  // A partition may be specified by a range, [a,b], or a vector of bin numbers.
  // If this function is called with a sampleName found in the m_samplePartition dictionary, 
  //   the function will proceed to fill only the specified bins.
  // If the function is called with a sampleName of '__all__' or the sampleName is not found
  //  in the m_samplePartition dictionary, the ChiSquared for all the bins is calculated.

  SamplePartition* samplepart = NULL;
  std::map< std::string, SamplePartition* >::iterator it = m_samplePartitions.find( sampleName );
  if ( sampleName=="__all__" || it==m_samplePartitions.end() )
    samplepart = NULL;
  else {
    samplepart = (*it).second;
  }

  for (int n=0; n<m_nbins; n++) {

    if ( samplepart!=NULL ) {
      if ( samplepart->userange ) {
	if ( samplepart->firstbin>n || n>samplepart->lastbin )
	  continue;
      }
      else {
	// This is super slow. Will lead to N^2 time due to linear search.  Should used an ordered list and binary tree search, but oh well, I will likely not use this feature.
	bool found = false;
	for (unsigned int b=0; b<samplepart->binlist.size(); b++) {
	  if ( samplepart->binlist.at(b)==n ) {
	    found=true;
	    break;
	  }
	}
	if (found==false) continue;
      }
    }
    
    m_binChiSquared[n] = 0.;
    binpull = CalculatePull( n, epsilons );
    
    if ( binpull<0 ) std::cout << "MinimizeErrors::CalculateChiSquared() - WARNING!! bin " << n 
			       << " has negative pull = " << binpull 
			       << " expected = " << m_NExpected[n]
			       << " observed = " << m_NObserved[n]
			       << std::endl; // for debug. most of the time if CS==nan, gamma fell negative
    
    if ( binpull*m_NExpected[n]>=0 )  {
      m_binChiSquared[n] += 2.0*( m_NExpected[n]*binpull - m_NObserved[n] );
      if (m_NObserved[n]>0.0) {
	m_binChiSquared[n] += 2.0*m_NObserved[n]*log( m_NObserved[n] );
	
	if (binpull*m_NExpected[n]>0.0)
	  m_binChiSquared[n] += 2.0*(-m_NObserved[n]*log(m_NExpected[n]*binpull)); // where we can go nan!
	else
	  m_binChiSquared[n] += 2.0*(-m_NObserved[n]*log(1.0e-50)); // this is a hack, what to do?
      }
    }
    else {
      // we drive the modified expectation to very small: 1.0e-6;
      double mod_expectation = 1.0e-50;
      error += 2.0*( mod_expectation - m_NObserved[n] );
      if (m_NObserved[n]>0.0) {
	error += 2.0*m_NObserved[n]*(log( m_NObserved[n] )-log(mod_expectation));
      }	
    }
    chi2 += m_binChiSquared[n];

  }//end of loop over bins

  // pull terms chi-squared
  double chi2_pull = CalculatePullChiSquared( epsilons );
  
  if (samplepart)
    samplepart->SampleChiSquared = chi2;
  
  //if ( GetVerbose()>=2 )
  std::cout << "Chi-squared calculated for sample " << sampleName
	    << ": " << chi2+chi2_pull << " ( = [res, " << chi2 << "]+[pull, " << chi2_pull << "]),  error=" << error << std::endl;
  
  chi2+=chi2_pull;
  return chi2;
  
}

double MinimizeErrors::CalculatePullChiSquared( double* epsilons ) {
  // pull terms chi-squared, i.e. the Penalty Term
  double chi2_pull = 0.;
  // even if we don't minimize with respect to some term, we still need to penalize some terms

  if ( epsilons==NULL ) {
    epsilons = m_x;
  }

  for (int k=0; k<m_nerrors; k++) {    
    double term_penalty = 0.;
    if ( m_penalize_errorterm[k] ) {
      //if ( m_term_pulls[k] ) {
      for (int j=0; j<m_nerrors; j++) {
	//if ( m_term_pulls[j] )
	if ( m_penalize_errorterm[j] )
	  term_penalty += (epsilons[k]-m_pull_min[k])*(*m_inv_cov)[k][j]*(epsilons[j]-m_pull_min[j]);
      }
    }
    if ( GetVerbose()>=3 )
      std::cout << "MinimizeErrors::CalculatePullChiSquared: " << m_error_list[k] << " value=" << epsilons[k] << " penalty = " << term_penalty <<  " (" << m_penalize_errorterm[k] << ")" << std::endl;
    chi2_pull += term_penalty;
  }
    
    
  if ( GetVerbose()>=4 ) {
    m_covariance->Print();
    m_inv_cov->Print();
    std::cout << "MinimizeErrors::CalculatePullChiSquared - penalty=" << chi2_pull << std::endl;
  }

  return chi2_pull;
}

