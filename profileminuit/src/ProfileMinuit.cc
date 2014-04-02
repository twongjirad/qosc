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

#include "ProfileMinuit.hh"
#include <assert.h>
#include <iostream>
#include <cmath>
#include <cstdlib>

#include "TTree.h"
#include "TFitter.h"
#include "TRandom3.h"

#include "SamplePartition.hh"
#include "FitterI.hh"
#include "ModelParameter.hh"
#include "ProfileMinuitI.hh"

using namespace qosc;

ProfileMinuit* ProfileMinuit::gGlobalProfileMinuitInstance = NULL;

ProfileMinuit::ProfileMinuit() {
  
  SetFitterInterface( NULL );

  fNewtonSolverTolerance = 1.0e-9;
  fNewtonVerbose = 0;
  fMinuitVerbose = 0;
  fNbins = -1;
  SetVerbose(0);
  SetThisInstanceAsGlobal();

  fInitialized = false;

  m_parMan = NULL;
  fParsIndexed = false;

  m_Minuit = NULL;
  fMinuitInitialized = false;
  fNVariedNewtonPars = -1;
  m_minuit_float_par.clear();

  m_pull_min = NULL;
  fPullMinimizerInitialized = false;
  fNVariedMinuitPars = -1;
  fLastMinimizerCallConverged = false;
  m_tries = 0;
  fUseF2ij = false;
  
  fGiveMinuitGradient = false;
  m_grad_nelements = 0;
  m_grad_nbins = 0;
  m_grad_npars = 0;
  bin_gradients = 0;

  m_lastChiSquared = -1;
  m_lastPreFitChiSquared = -1;
  m_lastPreFitPullChiSquared = -1;
  m_lastFitterChiSquared = -1;
  m_lastFitterPullChiSquared = -1;
  m_lastUserPenalty = -1;
  m_converged_flag = -1;
  m_minuit_converged_flag = -1;
  m_newton_converged_flag = -1;

  SetMinuitFCN( &DefaultProfileMinuitFCN );
}

ProfileMinuit::~ProfileMinuit() {
  // Clear the parameter dictionary
  //delete m_parMan; // not owner, so do not destroy
  delete m_Minuit;
  delete m_pull_min;
 }


// ------------------------------------------------------------------------------------------
// Initilization

void ProfileMinuit::Initialize( int numbins, ParameterManager* parMan ) {
  if ( GetVerbose()>=1 ) {
    std::cout << "---------------------------------" << std::endl;
    std::cout << " ProfileMinuit::Initialize " << std::endl;
  }

  SetNbins( numbins );
  if ( fNbins<0 ) {
    std::cout << "---------------------------------" << std::endl;
    std::cout << " ProfileMinuit::Initialize " << std::endl;
    std::cout << " ERROR. The number of bins in the analysis must be set before intializing class." << std::endl;
    std::cout << " Somewhere, somehow, ProfileMinuit::SetNbins( n ) must be called" << std::endl;
    std::cout << "---------------------------------" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  m_parMan = parMan;
  if ( !m_parMan ) {
    std::cout << "---------------------------------" << std::endl;
    std::cout << " ProfileMinuit::Initialize " << std::endl;
    std::cout << " ERROR. The instance of the parameter manager is NULL." << std::endl;
    std::cout << "---------------------------------" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  if ( !GetParManager()->IsParameterListBuilt() ) {
    //IndexParameters(); // need to index all the parameters. this is to facilitate interface with newton's method solver and minuit.
    std::cout << "---------------------------------" << std::endl;
    std::cout << " ProfileMinuit::Initialize " << std::endl;
    std::cout << " ERROR. The parameter manager has not been initialized." << std::endl;
    std::cout << "---------------------------------" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  else {
    fParsIndexed = true;
  }

  InitializePullMinimizer(); // create pull minimizer class

  InitializeMinuit(); // create minuit class

  fInitialized = true;

  BuildGradientArray( fNbins, GetParManager()->NumberOfParameters() );

  fNVariedNewtonPars = GetParManager()->GetNumVariedNewtonPars();
  fNVariedMinuitPars = GetParManager()->GetNumVariedMinuitPars();

  UpdateExpectation();

  if ( GetVerbose()>=1 ) {
    Print();
    std::cout << " End of ProfileMinuit::Initialize " << std::endl;
    std::cout << "---------------------------------" << std::endl;
    if ( GetVerbose()>=1 ) {
      std::cout << "Press [ENTER] to continue." << std::endl;
      std::cin.get();
    }
  }
}

void ProfileMinuit::SetFitterInterface( ProfileMinuitI* fitter_interface ) {
  // Override Fitter:SetFitterInterface( ... ) so that can type-check interface
  ProfileMinuitI* userInterface = dynamic_cast< ProfileMinuitI* >( fitter_interface );
  if ( fitter_interface && userInterface==NULL ) {
    std::cout << "------------------------------------------" << std::endl;
    std::cout << "ProfileMinuit::SetFitterIterface" << std::endl;
    std::cout << " The interface class passed to the fitter, " << fitter_interface << ", does not inherit from ProfileMinuitI." << std::endl;
    std::cout << "------------------------------------------" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  Fitter::SetFitterInterface( userInterface );
}

void ProfileMinuit::IndexParameters() {
  GetParManager()->Initialize();
  fParsIndexed = true;
}

void ProfileMinuit::InitializeMinuit() {

  if ( GetNMinuitTerms()==0 ) {
    fMinuitInitialized = false;
    m_Minuit = NULL;
  }
  
  m_Minuit = new TFitter( GetNMinuitTerms() );
  double p_quiet = 0;
  double p_strategy = 1; // 0=economical, 1=default, 2=precise
  if ( GetVerbose()>0 )
    p_quiet = 0;
  m_Minuit->ExecuteCommand("SET PRINTOUT", &p_quiet, 1);
  m_Minuit->ExecuteCommand("SET STR", &p_strategy, 1);
  assert( m_fcn!=NULL );
  m_Minuit->SetFCN( m_fcn );
    
  PushParametersToMinuit();
  
  fMinuitInitialized = true;


}

void ProfileMinuit::SetMinuitVerbose( int verbose ) {
  fMinuitVerbose = verbose;
  if ( m_Minuit ) {
    double dverbose = fMinuitVerbose-1; // so that 0 is quiet for all verbose modes
    m_Minuit->ExecuteCommand("SET PRINTOUT", &dverbose, 1);
  }
}


void ProfileMinuit::InitializePullMinimizer() {
  // even if newton solver not used. we still use this class to store the bins.
  // not the best system for sure
  m_pull_min = new MinimizeErrors( GetNParameters(), fNbins, IterativeSolver::kNewton ); // later can handle Jacobian method if needed.
  m_pull_min->SetStoppingTolerance( fNewtonSolverTolerance );

  // Fij must be added
  fPullMinimizerInitialized = true;

  PushParStatusToPullMinimizer();
  
  // Add error labels
  assert( fParsIndexed );
  std::vector< std::string > pars;
//   GetParManager()->GetListOfParameterNames( pars );
//   for ( std::vector< std::string >::iterator it=pars.begin(); it!=pars.end(); it++ ) {
//     m_pull_min->SetErrorLabel( GetParManager()->GetParameterIndex( *it ), *it );

  for ( int ipar=0; ipar<GetParManager()->NumberOfParameters(); ipar++) { 
    m_pull_min->SetErrorLabel( ipar, GetParManager()->GetParameterNameFromIndex( ipar ) );

//     if ( GetParManager()->GetParameter(*it)->GetFitterParType()==kNewtonTerm ) {
//       // If Pull Term
//       m_pull_min->VaryParameter( GetParManager()->GetParameterIndex(*it) );
//     }
//     else {
//       m_pull_min->FixParameter( GetParManager()->GetParameterIndex(*it) );
//       m_pull_min->TermDoesNotPull( GetParManager()->GetParameterIndex( *it ) );
//     }
    
// //     if ( !GetParManager()->GetParameter(*it)->DoesTermPull() )
// //       m_pull_min->TermDoesNotPull( GetParManager()->GetParameterIndex( *it ) );
    
  }
  
  // Add correlation matrix
  m_pull_min->ImportCovarianceMatrix( GetParManager()->GetCovarianceMatrix() );
  if ( GetVerbose()>=2 ) {
    std::cout << "Loaded Covariance Matrix:" << std::endl;
    GetParManager()->GetCovarianceMatrix()->Print();
  }
}

void ProfileMinuit::PushParStatusToPullMinimizer() {

  for ( ParameterManager::ParDictIter it=GetParManager()->ParDictBegin(); it!=GetParManager()->ParDictEnd(); it++ ) {

    if ( (*it).second->GetFitterParType()==kNewtonTerm ) {
      // If Pull Term
      if ( (*it).second->IsVaried() ) {
	m_pull_min->VaryParameter( (*it).second->GetID() );
      }
      else
	m_pull_min->FixParameter( (*it).second->GetID() );
    }
    else {
      m_pull_min->FixParameter( (*it).second->GetID() );
      m_pull_min->TermDoesNotPull( (*it).second->GetID() );
    }
    m_pull_min->SetPullMin( GetParManager()->GetParameterIndex((*it).first), (*it).second->GetCentralValue() );
    //std::cout << GetParManager()->GetParameterIndex((*it).first) << " " << (*it).second->GetCentralValue() << std::endl;
  }
  
}


// ------------------------------------------------------------------------------------------
// Toggle Parameter Flags

void ProfileMinuit::FloatPar( std::string name ) {
 if ( IsParameterDefined( name ) ) {
    ModelParameter* par = GetParameter(name);
    par->VaryParameter();
    // Maybe this shouldn't be here. Should set the variation at the time of running.
    if ( par->GetFitterParType()==kMinuitTerm && fMinuitInitialized) {
      std::cout << "ProfileMinuit::FloatPar(" << name << ") - Releasing parameter (MinuitIndex=" << GetParManager()->GetMinuitIndex(name) << ")" << std::endl;
      m_Minuit->ReleaseParameter( GetParManager()->GetMinuitIndex(name) );
    }
    else if ( par->GetFitterParType()==kNewtonTerm && fPullMinimizerInitialized ) {
      GetPullMinimizer()->VaryParameter( GetParManager()->GetParameterIndex( name ) );
    }
    else {
      assert(false);
    }
  }
}

void ProfileMinuit::FixPar( std::string name ) {
  if ( IsParameterDefined( name ) ) {
    ModelParameter* par = GetParameter(name);
    par->FixParameter();
    par->TermDoesNotThrow();
    // Maybe this shouldn't be here. Should set the variation at the time of running.
    if ( par->GetFitterParType()==kMinuitTerm && fMinuitInitialized) {
      m_Minuit->FixParameter( GetParManager()->GetMinuitIndex(name) );
    }
    else if ( par->GetFitterParType()==kNewtonTerm && fPullMinimizerInitialized ) {
      //GetPullMinimizer()->FixParameter( GetParManager()->GetNewtonIndex(name) );
      GetPullMinimizer()->FixParameter( GetParManager()->GetParameterIndex( name ) );
    }
    else {
      assert(false);
    }
  }
}

void ProfileMinuit::FixAllPars() {
  if ( fMinuitInitialized && GetNMinuitTerms()>0) {
    std::vector< std::string > minuitpars;
    GetParManager()->GetListOfMinuitParNames( minuitpars );
    for ( std::vector< std::string >::iterator it=minuitpars.begin(); it!=minuitpars.end(); it++ ) {
      FixPar( *it );
    }
  }
}

void ProfileMinuit::FloatAllPars() {
  if ( fMinuitInitialized && GetNMinuitTerms()>0) {
    std::vector< std::string > minuitpars;
    GetParManager()->GetListOfMinuitParNames( minuitpars );
    for ( std::vector< std::string >::iterator it=minuitpars.begin(); it!=minuitpars.end(); it++ ) {
      FloatPar( *it );
    }
  }
}

// ------------------------------------------------------------------------------------------
// Get/Set Parameter Values

void ProfileMinuit::SetParameterValue( std::string name, double value ) {
  if (!IsParameterDefined(name) ) assert(false);
  ModelParameter* par = GetParameter(name);
  par->SetValue(value);
  if ( par->GetFitterParType()==kMinuitTerm && fMinuitInitialized ) {
    double step = par->GetLengthScale();
    m_Minuit->SetParameter( GetParManager()->GetMinuitIndex( name ), name.c_str(), par->GetValue(), step, par->GetLowBound(), par->GetHighBound() );
  }
  else if ( par->GetFitterParType()==kNewtonTerm && fPullMinimizerInitialized ) {
    std::cout << "Setting newton parameters in the fitter not yet defined (may not make sense)" << std::endl;
  }
}

void ProfileMinuit::GetParametersFromMinuit() {
  if ( fMinuitInitialized && GetNMinuitTerms()>0) {
    if ( GetVerbose()>=2 )
      std::cout << __PRETTY_FUNCTION__ << std::endl;
    std::vector< std::string > minuitpars;
    GetParManager()->GetListOfMinuitParNames( minuitpars );
    for ( std::vector< std::string >::iterator it=minuitpars.begin(); it!=minuitpars.end(); it++ ) {
      if ( GetVerbose()>=2 )
	std::cout << "Minuit set " << *it << ": " << m_Minuit->GetParameter( GetParManager()->GetMinuitIndex(*it)  ) << std::endl;
      GetParameter( *it )->SetValue( m_Minuit->GetParameter( GetParManager()->GetMinuitIndex(*it) ) );
    }
  }
}

void ProfileMinuit::GetParametersFromPullMinimizer() {
  if ( fPullMinimizerInitialized ) {
    for ( ParameterManager::ParListIter it=GetParManager()->ParListBegin(); it!=GetParManager()->ParListEnd(); it++ ) {
      GetParameter( *it )->SetValue( GetPullMinimizer()->GetPullValue( GetParameter(*it)->GetID() ) );
    }
  }
}

void ProfileMinuit::PushParametersToMinuit() {
  std::vector< std::string > minuitpars;
  GetParManager()->GetListOfMinuitParNames( minuitpars );
  for ( std::vector< std::string >::iterator it=minuitpars.begin(); it!=minuitpars.end(); it++ ) {
    std::string parname = *it;
    ModelParameter* par = GetParManager()->GetParameter(parname);
    double step = par->GetLengthScale(); // I don't know what this value really does.
//     if ( par->IsBounded() ) { 
//       step = 0.1*fabs(par->GetHighBound()-par->GetLowBound() );
//       std::cout << parname << " is Bounded." << std::endl;
//     }
    m_Minuit->SetParameter(GetParManager()->GetMinuitIndex(parname), parname.c_str(), par->GetValue(), step, par->GetLowBound(), par->GetHighBound() );
    

    if ( !par->IsVaried() ) m_Minuit->FixParameter( GetParManager()->GetMinuitIndex(parname) );
  }
}

void ProfileMinuit::PushParametersToPullMinimizer() {
  for ( ParameterManager::ParListIter it=GetParManager()->ParListBegin(); it!=GetParManager()->ParListEnd(); it++ ) {
    std::string parname = *it;
    ModelParameter* par = GetParManager()->GetParameter(parname);
    GetPullMinimizer()->SetPullValue( par->GetID(), par->GetValue() );
  }
}


// ------------------------------------------------------------------------------------------
// Misc Tools

void ProfileMinuit::AddFitResultBranchesToTree( TTree* tree ) {
  tree->Branch("ChiSquared", &m_lastChiSquared, "ChiSquared/D");
  tree->Branch("PreFitChiSquared", &m_lastPreFitChiSquared, "PreFitChiSquared/D");
  tree->Branch("PreFitPullChiSquared", &m_lastPreFitPullChiSquared, "PreFitPullChiSquared/D");
  tree->Branch("FitterChiSquared", &m_lastFitterChiSquared, "FitterChiSquared/D");
  tree->Branch("FitterPullChiSquared", &m_lastFitterPullChiSquared, "FitterPullChiSquared/D");
  tree->Branch("UserPenalty", &m_lastUserPenalty, "UserPenalty/D");
  tree->Branch("Converged", &m_converged_flag, "Converged/I" );
  tree->Branch("MinuitConverged", &m_minuit_converged_flag, "MinuitConverged/I" );
  tree->Branch("NewtonConverged", &m_newton_converged_flag, "NewtonConverged/I" );
}

void ProfileMinuit::Print() {

  std::cout << "===============================================" << std::endl;
  std::cout << "Margin Minuit (" << this << ")" << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << "Number of Bins: " << fNbins << std::endl;
  std::cout << "Number of Parameters: " << GetNParameters() << std::endl;
  std::cout << "  Minuit Parameters: " << fNVariedMinuitPars << " varied, " << GetNMinuitTerms() << " total" << std::endl;
  std::cout << "  Pull Terms: " << fNVariedNewtonPars << " varied, " << GetNNewtonTerms() << " total " << std::endl;
  if ( !fInitialized ) {
    std::cout << "===============================================" << std::endl;
    return;
  }
  
  std::cout << "-----------------------------------------------" << std::endl;
  if ( GetNMinuitTerms()>0 && fMinuitInitialized ) {
    std::cout << "Minuit Parameters [MinuitIndex/FitterIndex]" << std::endl;
    for ( int n=0; n<GetNMinuitTerms(); n++ ) {
      std::string parname = GetParManager()->GetParameterNameByMinuitIndex( n );
      std::cout << "  [" << n << " / " << GetParManager()->GetParameterIndex( parname ) << "] ";
      std::cout << parname;
      if ( GetParameter(parname)->IsVaried() ) std::cout << " [VARY] ";
      else std::cout << " [FIXED] ";
      std::cout << " = " << m_Minuit->GetParameter( GetParManager()->GetMinuitIndex( parname ) ) << std::endl;
    }
  }

  if ( fMinuitInitialized && fPullMinimizerInitialized ) 
    std::cout << "-----------------------------------------------" << std::endl;

  if ( fPullMinimizerInitialized ) {
    // Error list
    std::cout << "Newton Parameters [NewtonIndex/FitterIndex] [status] [Parname] [Info]" << std::endl;
    for (int n=0; n<GetNNewtonTerms(); n++) {
      std::string parname = GetParManager()->GetParameterNameByNewtonIndex( n );
      ModelParameter* par = GetParManager()->GetParameter( parname );
      std::cout << "  [" << n << " / " << GetParManager()->GetParameterIndex( parname ) << "] ";
      if ( par->IsVaried() ) std::cout << "[VARY] ";
      else std::cout << "[FIXED] ";
      std::cout << parname << " value=" << par->GetValue()  << " : "
		<< " central=" << par->GetCentralValue() 
		<< " value/sigma=" << par->GetValue()/par->GetSigma()
		<< " dist_from_central=" << ( par->GetValue()-par->GetCentralValue() )/par->GetSigma() 
		<< std::endl;
    }
    //m_pull_min->PrintBinList();
  }
  std::cout << "===============================================" << std::endl;
}

void ProfileMinuit::PrintNewtonFitterInfo() {
  m_pull_min->PrintErrorList();
  m_pull_min->PrintBinList();
  m_pull_min->PrintFij();
}

// ------------------------------------------------------------------------------------------
// Run Minimization

void ProfileMinuit::DoFit() {
  // returns true if converges
  // Run Minuit: this will pick a new set of osc parameters (and pull values)

  if ( !IsInitialized() ) {
    std::cout << "------------------------------------------" << std::endl;
    std::cout << "ProfileMinuit::RunFit" << std::endl;
    std::cout << " Class has not been properly initialized." << std::endl;
    std::cout << " Make sure the method ProfileMinuit::Initialize() is called" << std::endl;
    std::cout << " either in the main, but probably in any overriding intialize calls" << std::endl;
    std::cout << "------------------------------------------" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  if ( GetVerbose()>=2 ) {
    std::cout << "------------------------------------------" << std::endl;
    std::cout << "ProfileMinuit::RunFit" << std::endl;
  }

  fNVariedNewtonPars = GetParManager()->GetNumVariedNewtonPars();
  fNVariedMinuitPars = GetParManager()->GetNumVariedMinuitPars();

  GetPullMinimizer()->UseF2ij( fUseF2ij );

  // Store the central values of parameters into pull minimizer
  for ( ParameterManager::ParListIter it=GetParManager()->ParListBegin(); it!=GetParManager()->ParListEnd(); it++ ) {
    SetPullCentralValue( *it, GetParManager()->GetParameter( *it )->GetCentralValue() );
  }
  PushParametersToPullMinimizer();
  std::cout << "----- [ Initial pull values, X ] -----" << std::endl;
  GetPullMinimizer()->PrintX();
  std::cout << "--------------------------------------" << std::endl;

  // Pre-MINUIT
  PushParametersToMinuit();
  UpdateExpectation();
//   UserLoadExpectationBins( GetParManager() ); /// User must load the expected number of events in the bins
//   UserLoadObservedBins( GetParManager() ); /// Likewise, user must load the observed number of events in the bins
//   UserLoadFijValues( GetParManager() ); /// Finally, User must calculated Fij values
  double prefit_chi2 = GetChiSquared();
  double external_penalty = GetLastUserPenalty();
  m_lastPreFitPullChiSquared = GetPullChiSquared();
  StoreLastPrefitChiSquared( prefit_chi2+external_penalty );  
  RunPullMinimizer();

  bool minuitconverged;
  if ( GetNumVariedMinuitTerms()>0 )
    minuitconverged = RunMinuit();
  else
    minuitconverged = true;

  bool converged = minuitconverged & fLastMinimizerCallConverged;

  // store the convergence info (variable below are stored in a tree)
  ( converged ) ? m_converged_flag = 1 : m_converged_flag = 0;
  ( minuitconverged ) ? m_minuit_converged_flag = 1 : m_minuit_converged_flag = 0;
  ( fLastMinimizerCallConverged ) ? m_newton_converged_flag = 1 : m_newton_converged_flag = 0;

  // Get parameters from MINUIT and calculate last chi-squared
  GetParametersFromMinuit(); 
  UpdateExpectation();
//   UserLoadExpectationBins( GetParManager() ); /// User must load the expected number of events in the bins
//   UserLoadObservedBins( GetParManager() ); /// Likewise, user must load the observed number of events in the bins
//   UserLoadFijValues( GetParManager() ); /// Finally, User must calculated Fij values
  bool converges = RunPullMinimizer();

  m_lastFitterChiSquared = GetPullMinimizer()->CalculateChiSquared();
  m_lastFitterPullChiSquared = GetPullMinimizer()->CalculatePullChiSquared();
  m_lastUserPenalty = UserCalculateAdditionalPenalty( GetParManager() );
  m_lastChiSquared = m_lastFitterChiSquared+m_lastUserPenalty;


  if ( GetVerbose()>=1 ) {
    Print();
    std::cout << "End of ProfileMinuit::RunFit, status=";
    if ( converged ) std::cout << "converged" << std::endl;
    else {
      std::cout << "failed" << std::endl;
      std::cout << "  Minuit converged: " << minuitconverged << std::endl;
      std::cout << "  Newton converged: " << fLastMinimizerCallConverged << std::endl;
    }
    std::cout << "Final Chi-Squared: " << m_lastChiSquared << std::endl;
    std::cout << "------------------------------------------" << std::endl;
  }
  
}

bool ProfileMinuit::RunMinuit() {

  const char * interp[13];
  interp[ 0] = " 0: command executed normally";
  interp[ 1] = " 1: command is blank, ignored";
  interp[ 2] = " 2: command line unreadable, ignored";
  interp[ 3] = " 3: unknown command, ignored";
  interp[ 4] = " 4: abnormal termination (e.g., MIGRAD not converged)";
  interp[ 5] = " 5: command is a request to read PARAMETER definitions";
  interp[ 6] = " 6: 'SET INPUT' command";
  interp[ 7] = " 7: 'SET TITLE' command";
  interp[ 8] = " 8: 'SET COVAR' command";
  interp[ 9] = " 9: reserved";
  interp[10] = "10: END command";
  interp[11] = "11: EXIT or STOP command";
  interp[12] = "12: RETURN command";

  DefaultMinuitParValues();
  PushParametersToMinuit();

  // returns true if converges
  double miparms[2] = { 5000, 1000 };
  SetMinuitVerbose( fMinuitVerbose );

  if ( GetVerbose()>=2 ) 
    std::cout << "ProfileMinuit::RunMinuit [ max calls=" << miparms[0] << ", tolerance=" << miparms[1] << " ]" << std::endl;

  if ( fGiveMinuitGradient ) {
    double forcegrad[1];
    forcegrad[0] = 0.;
    GetMinuit()->ExecuteCommand("SET GRAD", forcegrad, 1 );
  }
  
  //int ret = m_Minuit->ExecuteCommand( "MIGRAD", miparms, 2 );
  int ret = m_Minuit->ExecuteCommand( "MINImize", NULL, 0 );
  //int ret = m_Minuit->ExecuteCommand( "MINImize", NULL, 0 );
  if ( GetVerbose()>=1 ) std::cout << "MINUIT return value: " << interp[ret] << std::endl;
  //std::cout << "==========================================================================" << std::endl;
  //std::cout << "MINUIT fit complete" << std::endl;
  //m_Minuit->PrintResults(1,0);
  //std::cout<< "==========================================================================" << std::endl;

  return (ret != 4);

}

bool ProfileMinuit::RunPullMinimizer() {
  SetNewtonSolverVerbose( fNewtonVerbose );
  
  if ( GetNNewtonTerms() == 0 || GetNumVariedNewtonTerms()==0 ) {
    if ( GetVerbose()>=2 ) std::cout << "ProfileMinuit::RunPullMinimizer: Skipped because zero Newton parameters" << std::endl;
    fLastMinimizerCallConverged = true;
    return true;
  }
  
  if ( GetVerbose()>=2 )
    std::cout << "ProfileMinuit::RunPullMinimizer" << std::endl;
  
  // we've assumed the Marginalizer has been properly setup.
  // (1) observed bin values have been stored
  // (2) the nominal expected bin values have been stored
  // (3) Fij values have been stored
  
  int fNPullTerms = GetNMinuitTerms()+GetNNewtonTerms();
  double pullValues[fNPullTerms];
  
  // first get a good initial value from which to start the iterations
  double linearSolution[fNPullTerms];
  
  // start with current values
  for (ParameterManager::ParListIter it=GetParManager()->ParListBegin(); it!=GetParManager()->ParListEnd(); it++) {
    if ( GetParManager()->GetParameter( *it )->GetFitterParType()!=kNewtonTerm )
      linearSolution[ GetParManager()->GetParameterIndex( *it ) ] = GetParManager()->GetParameterValue( *it );
    else
      linearSolution[ GetParManager()->GetParameterIndex( *it ) ] = 0.;
  }
  
  // solve for pulls analytically  
  if ( GetNNewtonTerms()>0 ) {
    if ( GetVerbose()>=2 )
      std::cout << "Calling MinimizerErrors::GetLinearizedSolution()" << std::endl;
    m_pull_min->GetLinearizedSolution( linearSolution );
  
    if ( GetNNewtonTerms()>0 )
      m_pull_min->MakeValidSolution( pullValues );
    for (int k=0; k<fNPullTerms; k++) pullValues[k] = linearSolution[k];
  }
  
  bool minimizerConverged = false;
  int maxtries = 100;
  
  TRandom3 rand(1);
  
  m_tries = 0;
  //while ( fMinimizeErrors==true && m_tries<maxtries && minimizerConverged==false) {
  while ( m_tries<maxtries && minimizerConverged==false) {
    if ( GetVerbose()>=2 ) {
      std::cout << "ProfileMinuit::RunPullMinimizer(). Newton Minimizer called. Attempt " << m_tries << " of " << maxtries << std::endl;      
    }
    
    // set initial point
    m_pull_min->SetInitialPoint( pullValues );

    double chierror = 0.;

    if ( GetNNewtonTerms()<1 ) break;
    
    // try to solve for epsilons given a set of Fij
    minimizerConverged = m_pull_min->RunIteration( pullValues );

    double ChiSquared = m_pull_min->CalculateChiSquared( pullValues, chierror );
    
    if ( minimizerConverged  ) {
      if ( ChiSquared==ChiSquared ) { // checking for NaN
	minimizerConverged = true;
      }
      else {
	minimizerConverged = false;
	do { 
	  // try another set of random initial pull terms. This time randomly set.
	  for (int k=0; k<fNPullTerms; k++)
	    pullValues[k] = 10.0*(rand.Uniform()-0.5);
	} while ( m_pull_min->IsSolutionValid( pullValues )==false );
	m_pull_min->MakeValidSolution( pullValues );
      }
    }
    else {
      // Reset pull values for another attempt
      do { 
	for (int k=0; k<fNPullTerms; k++)
	  pullValues[k] = 10.0*(rand.Uniform()-0.5);
      } while ( m_pull_min->IsSolutionValid( pullValues )==false );
      m_pull_min->MakeValidSolution( pullValues );
    }
    m_tries++;
    
    if ( GetVerbose()>=2 )  {
      std::cout << "ProfileMinuit::RunPullMinimizer()  Fit returned with Chi-Squared=" << ChiSquared << " (chi-error=" << chierror << ") status=";
      if ( minimizerConverged ) std::cout << "converged";
      else std::cout << "failed";
      std::cout << std::endl;
    }

  }//end of while not enough tries and minizer has not converged

  // pass parameter values found by newton minimizer back to parameters in parameter manager
  for ( ParameterManager::ParListIter it=GetParManager()->NewtonParListBegin(); it!=GetParManager()->NewtonParListEnd(); it++ ) {
    GetParameter( *it )->SetValue( pullValues[ GetParManager()->GetParameterIndex( *it ) ] );
  }
  
  if ( GetVerbose()>=2 ) {
    std::cout << "End of ProfileMinuit::RunPullMinimizer. Resulting parameter values:" << std::endl;
    for (int n=0; n<GetNParameters(); n++) {
      double sigma = GetParManager()->GetParameterFromIndex( n )->GetSigma();
      double value = GetParManager()->GetParameterFromIndex( n )->GetValue();
      std::cout << "[" << n << "] " << GetParManager()->GetParameterFromIndex( n )->GetName() << " value=: " << value << " value/sigma=" << value/sigma << std::endl;
    }
  }
  
  fLastMinimizerCallConverged = minimizerConverged;
  return minimizerConverged;
}
  
// ------------------------------------------------------------------------------------------
// Allowed interfaces to pull minimizer

double ProfileMinuit::GetChiSquared() {
  // Must use all parameters now, not just newton parameters
  // However, the indexing doesn't work for this yet
  int fNTerms = GetNMinuitTerms()+GetNNewtonTerms();
  double pullValues[ fNTerms ];
  std::vector< std::string > parlist; 
  for ( ParameterManager::ParListIter it=GetParManager()->ParListBegin(); it!=GetParManager()->ParListEnd(); it++ ) {
    pullValues[ GetParManager()->GetParameterIndex(*it) ] = GetParameter( *it )->GetValue();
  }
  
  double errors = 0;
  return m_pull_min->CalculateChiSquared( pullValues, errors );
}

double ProfileMinuit::GetPullChiSquared() {
  // Must use all parameters now, not just newton parameters
  // However, the indexing doesn't work for this yet
  int fNTerms = GetNMinuitTerms()+GetNNewtonTerms();
  double pullValues[ fNTerms ];
  std::vector< std::string > parlist; 
  for ( ParameterManager::ParListIter it=GetParManager()->ParListBegin(); it!=GetParManager()->ParListEnd(); it++ ) {
    pullValues[ GetParManager()->GetParameterIndex(*it) ] = GetParameter( *it )->GetValue();
  }
  
  return m_pull_min->CalculatePullChiSquared( pullValues );
}

void ProfileMinuit::SetBinLabel( int binnum, std::string name ) {
  m_pull_min->SetBinLabel( binnum, name );
}

void ProfileMinuit::SetBinNExpected( int binnum, double expected, double* gradient, int ngradpars ) {
  m_pull_min->SetBinNExpected( binnum, expected );
  if ( gradient!=NULL && ngradpars!=GetParManager()->NumberOfParameters() ) {
    std::cout << "WARNING - ProfileMinuit::SetBinExpected(...) was given " << std::endl;
    std::cout << "a number of gradient parameters not matching the number of parameters" << std::endl;
    std::cout << "given=" << ngradpars << " num minuit pars=" << GetParManager()->NumberOfParameters() << std::endl;
    exit(EXIT_FAILURE);
  }
  if ( gradient ) {
    memcpy( GetBinGradient( binnum ), gradient, sizeof(double)*ngradpars );
  }

}

void ProfileMinuit::SetBinNObserved( int binnum, double observed ) {
  m_pull_min->SetBinNObserved( binnum, observed );
}

void ProfileMinuit::SetFij( std::string pullName, int binNum, double fij ) {
  if ( !IsParameterDefined( pullName ) ) assert(false);
  m_pull_min->SetFij( GetParManager()->GetParameterIndex( pullName ), binNum, fij );
}

void ProfileMinuit::SetFij( int parIndex, int binNum, double fij ) {
  m_pull_min->SetFij( parIndex, binNum, fij );
}

void ProfileMinuit::PrintFij() {
  m_pull_min->PrintFij();
}

void ProfileMinuit::SetPullCentralValue( std::string parname, double minimum ) {
  if ( !IsParameterDefined( parname ) ) {
    std::cout << "ProfileMinuit::SetPullCentralValue - '" << parname << "' has not yet been defined." << std::endl;
    assert(false);
    exit(EXIT_FAILURE);
  }

  if ( !m_pull_min || !fPullMinimizerInitialized ) {
    std::cout << "ProfileMinuit::SetPullCentralValue - the Pull Minimizer has not yet been initilized yet." << std::endl;
    assert(false);
    exit( EXIT_FAILURE );
  }
  
  m_pull_min->SetPullMin( GetParManager()->GetParameterIndex( parname ), minimum );
  
}

// ------------------------------------------------------------------------------------------
// GRADIENT MANAGEMENT AND CALCULATIONS

void ProfileMinuit::SendGradientToMinuit( bool doit ) {
  fGiveMinuitGradient = doit;
}

void ProfileMinuit::BuildGradientArray( int nbins, int npars ) {
  m_grad_nbins = nbins;
  m_grad_npars = npars;
  m_grad_nelements = nbins*npars;
  bin_gradients = new double[ m_grad_nelements ];
  memset( bin_gradients, 0, sizeof(double)*m_grad_nelements );
}

void ProfileMinuit::DestroyGradientArray() {
  delete [] bin_gradients;
  bin_gradients = NULL;
}

double ProfileMinuit::GetGradientElement( int ibin, int ipar ) {
  return *(bin_gradients + m_grad_npars*ibin + ipar );
}

void ProfileMinuit::SetGradientElement( int ibin, int ipar, double grad ) {
  *(bin_gradients + m_grad_npars*ibin + ipar ) = grad;
}

double* ProfileMinuit::GetBinGradient( int ibin ) {
  return ( bin_gradients + m_grad_npars*ibin );
}

double ProfileMinuit::GetChiSquaredPartial( int parid ) {
  // calculate partial derivative with respect to parameter with id
  // really only use this MINUIT terms.
  double partialChi = 0.;
  for (int bin=0; bin<fNbins; bin++) {
    double dE = GetGradientElement( bin, parid );
    double E = GetPullMinimizer()->GetBinNExpected( bin );
    double O = GetPullMinimizer()->GetBinNObserved( bin );
    partialChi += 2*dE;
    if ( O>0 && E>0) {
      partialChi -= 2*dE*(O/E);
    }
  }
  
  return partialChi;
}

// ------------------------------------------------------------------------------------------

void ProfileMinuit::UpdateExpectation() {
  GetParametersFromMinuit(); /// asks the Minuit class what its parameter values are and loads them into the parameter objects kept in the ParameterManager class.
  PushParametersToPullMinimizer();
  //GetMMI()->UpdateModel();
  GetMMI()->UpdateFitter();
}

void ProfileMinuit::UserLoadExpectationBins( ParameterManager* parameters ) {
  std::cout << " ------------------------------------------------------------------------------------------ " << std::endl;
  std::cout << " MARGINMINUIT::UserLoadExpectationBins" << std::endl;
  std::cout << "  Since the program has called this function, something is wrong." << std::endl;
  std::cout << "  User must overload this function and use it to set the expected number of events in each bin." << std::endl;
  std::cout << "  The number of each bins should have been set already via SetNbins( int nbins )" << std::endl;
  std::cout << "  Set the expectation via: void SetBinNExpected( int binnum, double expected )" << std::endl;
  std::cout << "  Bin numbers start at 0." << std::endl;
  std::cout << " ------------------------------------------------------------------------------------------ " << std::endl;
  assert(false);
}

void ProfileMinuit::UserLoadObservedBins( ParameterManager* parameters ) {
  std::cout << " ------------------------------------------------------------------------------------------ " << std::endl;
  std::cout << " MARGINMINUIT::UserLoadObservedBins" << std::endl;
  std::cout << "  Since the program has called this function, something is wrong." << std::endl;
  std::cout << "  User must overload this function and use it to set the observed number of events in each bin." << std::endl;
  std::cout << "  The number of each bins is set by 'SetNbins( int nbins )'" << std::endl;
  std::cout << "  Set the observed number of events in a bin via 'SetBinNObserved( int binnum, double expected )'" << std::endl;
  std::cout << "  Bin numbers start at 0." << std::endl;
  std::cout << " ------------------------------------------------------------------------------------------ " << std::endl;
  assert(false);
}

void ProfileMinuit::UserLoadFijValues( ParameterManager* parameters ) {

  for (int n=0; n<GetNbins(); n++ ) {
    for ( ParameterManager::ParListIter it = GetParManager()->NewtonParListBegin(); it!=GetParManager()->NewtonParListEnd(); it++ ) {
      double Nshift = UserReturn1SigmaShiftInBinByNewtonParameter( n, *it, parameters->GetParameter( *it ) );
      double fij = 0;
      double Nexpected = GetPullMinimizer()->GetBinNExpected( n );
      if ( Nexpected!=0.0 ) fij = Nshift/Nexpected;
      SetFij( *it, n, fij );
    }
  }
}

double ProfileMinuit::UserReturn1SigmaShiftInBinByNewtonParameter( int binNum, std::string parname, ModelParameter* par ) {
  std::cout << " ------------------------------------------------------------------------------------------ " << std::endl;
  std::cout << " MARGINMINUIT::UserReturn1SigmaShiftInBinByNewtonParameter" << std::endl;
  std::cout << "  Since the program has called this function, something is wrong." << std::endl;
  std::cout << "  The user has two options for loading Fij values to the fitter." << std::endl;
  std::cout << "   (1) Overload 'double UserReturn1SigmaShiftInBinByNewtonParameter(...)'. The code will automatically put the value in the proper place." << std::endl;
  std::cout << "       This method is called by UserLoadFijValues(...)." << std::endl;
  std::cout << "   (2) Overload void UserLoadFijValues(...)" << std::endl;
  std::cout << "       For this method, the user must load a total number of such terms." << std::endl;
  std::cout << "       The total is equal to K x N where K is the number of Pull terms, " << std::endl;
  std::cout << "         while N is the number of bins." << std::endl;
  std::cout << "       The number of bins is can be retrieved via 'GetNbins()'. Bin numbers start at 0." << std::endl;
  std::cout << "       The user can get a list of pull term parameters (aka Newton parameters) via 'GetListOfNewtonParNames( std::vector< std::string >& )'" << std::endl;
  std::cout << "       Finally, set the Fij values via: SetFij( std::string parname, int binNum, double fij )" << std::endl;
  std::cout << " ------------------------------------------------------------------------------------------ " << std::endl;
  assert(false);  
}


// ------------------------------------------------------------------------------------------

void DefaultProfileMinuitFCN( int &ndim, double* grad_optional, double &result, double par[], int flag ) {
  ProfileMinuit* mm = ProfileMinuit::GetGlobalInstance();


  if ( mm->GetVerbose() ) {
    std::cout << "=========================================================================================================" << std::endl;
    std::cout << " DefaultProfileMinuitFCN" << std::endl;
    std::cout << "---------------------------------------------------------------------------------------------------------" << std::endl;
  }

  mm->UpdateExpectation();

  if ( mm->GetVerbose()>=3 ) mm->PrintNewtonFitterInfo();

  bool converges = mm->RunPullMinimizer();
  if (!converges) 
    std::cout << "!!! PullMinimizer did not converge" << std::endl;
  double postfit_chi2 = mm->GetChiSquared();
  double postfit_pull_chi2 = mm->GetPullChiSquared();
  double external_penalty = mm->GetLastUserPenalty();
  result = postfit_chi2+external_penalty;

  if ( flag==2 ) {
    ParameterManager* parMan = mm->GetParManager();
    std::cout << "Chi-2 gradient [flag=" << flag << "]" << std::endl;
    double chi2_grad[ parMan->NumberOfParameters() ];
    bool calcdone = mm->UserCalculateChi2Gradient( chi2_grad, parMan );
    if ( calcdone ) 
      memcpy( grad_optional, chi2_grad, sizeof(double)*parMan->NumberOfParameters() );
  }
  
  if ( mm->GetVerbose()>=1 ) {
    std::cout << "------------------------------------" << std::endl;
    std::cout << "DefaultProfileMinuitFunction[ flag=" << flag << "]" << std::endl;
    //std::cout << " prefit chi2: " << prefit_chi2 << std::endl;
    std::cout << " postfit_chi2: " << postfit_chi2 
	      << "  res=" << postfit_chi2-postfit_pull_chi2
	      << "  pull=" << postfit_pull_chi2
	      << std::endl;
    std::cout << " user external penalty: " << external_penalty << std::endl;
    //std::cout << " log( det|Jacobian| ): " << log(fisher_information) << std::endl;
    std::cout << " returned chi2: " << result << std::endl;
    std::cout << "---------------------------------------------------------------------------------------------------------" << std::endl;
    std::cout << " End of DefaultProfileMinuitFCN" << std::endl;
    std::cout << "=========================================================================================================" << std::endl;

    if ( mm->GetVerbose()>=3 ) {
      std::cout << "Hit [Enter] to continue." << std::endl;
      std::cin.get();
    }
  }
}

// ------------------------------------------------------------------------------------------
