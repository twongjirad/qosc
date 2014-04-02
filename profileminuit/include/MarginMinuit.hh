//-*- mode:c++; c-basic-offset:2;   -*-
#ifndef __MarginMinuit__
#define __MarginMinuit__

/*
 * -----------------------------------------------------------------------------------------------------------------------
 * \class MarginMinuit
 *  \defgroup MarginMinuit
 *  \brief Class responsible for managing the fit.
 *
 *  Verbose Levels (for MarginMinuit class. Different for Iterative solver and Minuit verbose levels)
 *   0: Quiet
 *   1: Print out summary of initialization and fitting
 *   2: Dump status of par values during each step
 *   3: Pause after each fit
 *   4: Pause at major steps in process
 *  For Newton Solver see MinimizeErrors class.
 *  For Minuit (note that verbose integer is offset by +1 from default MINUIT levels so that 0 is quiet for all pieces)
 *   0: Quiet
 *   1: Minimal output
 *   2: Normal output
 *   3: Additional output giving intermediate steps
 *   4: Maximum output, showing progress of minimizations
 * -----------------------------------------------------------------------------------------------------------------------
 */

#include <string>
#include <map>
#include <iostream>
#include <assert.h>
#include "MinimizeErrors.hh"
#include "TFitter.h"
#include "TTree.h"
#include "ModelParameter.hh"
#include "ParameterManager.hh"
#include "Fitter.hh"
#include "FitterI.hh"
#include "MarginMinuitI.hh"

typedef void (*MinuitFCN) ( int&, double*, double&, double[], int );

void DefaultMarginMinuitFCN( int&, double*, double&, double[], int );

class MarginMinuit : public Fitter {

public:
  MarginMinuit();
  virtual ~MarginMinuit();


  // Function Implementations From Fitter Base Class
  virtual void DoFit();
  virtual void SetFitterInterface( MarginMinuitI* interface ); /// overriding interface loader to type-check for MarginMinuitI*
  MarginMinuitI* GetMMI() { return (MarginMinuitI*)Fitter::GetFitterInterface(); }; /// way to get interface pointer typed as MarginMinuitI*

  // Class methods
public:
  void UpdateExpectation();
  void Initialize( int numbins, ParameterManager* );

  void SetNbins( int nbins ) { fNbins = nbins; };
  int GetNbins() { return fNbins; };

  void SetParValue( std::string name, double value );
  void SetParBounds( std::string name, double low, double high );

  void FloatPar( std::string name );
  void FixPar( std::string name );
  void FitParameter( std::string name, bool fit ) {
    if ( fit ) FloatPar( name );
    else FixPar( name );
  };
  void FixAllPars();
  void FloatAllPars();
  void ZeroParameterValues() {GetParManager()->ZeroParameterValues();}
  void ZeroNewtonParValues() {GetParManager()->ZeroNewtonParValues();}
  void ZeroMinuitParValues() {GetParManager()->ZeroMinuitParValues();}
  void DefaultMinuitParValues() {GetParManager()->DefaultMinuitParValues();}
  void UseSecondOrderFij( bool use=true ) { fUseF2ij = use; };
  bool DoWeUseF2ij() { return fUseF2ij; };

  // Called by user to run minimizers
  bool RunMinuit();
  bool RunPullMinimizer();
  void Print();
  void PrintNewtonFitterInfo();
  bool IsInitialized() { return fInitialized; };

  int GetNMinuitTerms() { return GetParManager()->NumberOfMinuitPars(); };
  int GetNNewtonTerms() { return GetParManager()->NumberOfNewtonPars(); };
  int GetNumVariedNewtonTerms() { return fNVariedNewtonPars; };
  int GetNumVariedMinuitTerms() { return fNVariedMinuitPars; };
  int GetNParameters() { return GetParManager()->NumberOfTotalPars(); };

  int GetVerbose() { return fVerbose; };
  void SetVerbose( int verbose ) { fVerbose=verbose; };
  void SetNewtonSolverVerbose( int verbose ) { fNewtonVerbose = verbose; if (m_pull_min) m_pull_min->SetVerbose(fNewtonVerbose); };
  int GetNewtonSolverVerbose() { if (m_pull_min) return m_pull_min->GetVerbose(); else return 0; };
  void SetMinuitVerbose( int verbose );

  void MinuitFunction( int&, double*, double&, double[], int ); // This is the default function we tell MINUIT to call.
  void SetPullCentralValue( std::string parname, double minimum );

  // ---------------------------------------------------------------
  // Needed User Functions due to lack of class that generically handles data and sets up analysis. To do.
public:
  virtual void UserLoadExpectationBins( ParameterManager* parameters );
  virtual void UserLoadObservedBins( ParameterManager* parameters );
  virtual void UserLoadFijValues( ParameterManager* parameters );
  virtual double UserReturn1SigmaShiftInBinByNewtonParameter( int binNum, std::string parName, ModelParameter* par );
  virtual double UserCalculateAdditionalPenalty( ParameterManager* parameters ) { return 0; };
  virtual double UserCalculateFisherInformationDet( ParameterManager* parameters ) { return 1.; };
  virtual void UserCalculateBinGradient( double* bin_gradient, ParameterManager* parameters ) { };
  virtual bool UserCalculateChi2Gradient( double* chi2_gradient, ParameterManager* parameters ) { return false; };
  virtual void UserPreCalculation( ParameterManager* parameters ) {};
  virtual void UserPostCalculation( ParameterManager* parameters ) {};
  // ---------------------------------------------------------------

  // ---------------------------------------------------------------
  // Because of the way we need to interface to MINUIT, we need a global pointer
public:
  static MarginMinuit* gGlobalMarginMinuitInstance;
  static MarginMinuit* GetGlobalInstance() { return gGlobalMarginMinuitInstance; };
  void SetThisInstanceAsGlobal() { MarginMinuit::gGlobalMarginMinuitInstance = this; };
  // ---------------------------------------------------------------

  // ----------------------------------------
  // Internal Parameters  
private:  
  bool fInitialized; ///makes sure we call the initializer for minuit
  int fVerbose;
  int fNewtonVerbose;
  int fMinuitVerbose;
  int fNbins;

  // ----------------------------------------
  // Parameter Manager
private:
  ParameterManager* m_parMan;
public:
  ParameterManager* GetParManager() { return m_parMan; };
  void SetParManager( ParameterManager* parman ) {
    if ( dynamic_cast< ParameterManager* >( parman )==false ) {
      std::cout << "MarginMinuit::SetParManager -- ERROR. Any Parameter Manager given to MarginMinuit must inherit from the ParameterManager class." << std::endl;
      assert(false);
    }
    if ( m_parMan ) delete m_parMan;
    m_parMan = parman;
  };
private:
  bool fParsIndexed;
protected:
  void IndexParameters(); /// called by initialize

  // ----------------------------------------
  // Get/Set
public:
  bool IsParameterDefined( std::string name ) {
    return GetParManager()->IsParameterDefined( name );
  };
  ModelParameter* GetParameter( std::string name ) {
    if ( IsParameterDefined( name ) ) return GetParManager()->GetParameter( name );
    else return NULL;
  };
  double GetParameterValue( std::string name ) {
    if (!IsParameterDefined(name) ) assert(false);
    return GetParameter( name )->GetValue();
  };
  void SetParameterValue( std::string name, double value );
  
  void AddParameter( ModelParameter* par ) {
    GetParManager()->RegisterParameter( par );
  };
  void GetListOfNewtonParNames( std::vector< std::string >& parlist ) { GetParManager()->GetListOfNewtonParNames( parlist ); };
  void GetListOfMinuitParNames( std::vector< std::string >& parlist ) { GetParManager()->GetListOfNewtonParNames( parlist ); };

  // --------------------------------------------------------
  // Parameter and Fitter Interface
public:
  void GetParametersFromMinuit(); ///< extracts minuit parameter values from the minuit class and stores them in our parameters
  void GetParametersFromPullMinimizer();
  void PushParametersToMinuit();
  void PushParametersToPullMinimizer();
  void PushParStatusToPullMinimizer();

  // --------------------------------------------------------
//   // Might be good to have such items
//   typedef std::vector< std::string >::iterator ParListIter;
//   ParListIter MinuitListBegin() { return m_minuit_pars.begin(); };
//   ParListIter MinuitListEnd()   { return m_minuit_pars.end(); };
//   ParListIter PullListBegin()   { return m_pull_pars.begin(); };
//   ParListIter PullListEnd()     { return m_pull_pars.end(); };
    
  // --------------------------------------------------------
  // Minuit
private:
  TFitter* m_Minuit;
  MinuitFCN m_fcn;
  bool fMinuitInitialized;
  int fNVariedMinuitPars;
  std::map< std::string, bool > m_minuit_float_par;
  void InitializeMinuit();
protected:
  void SetMinuitFCN( MinuitFCN fcn ) { m_fcn = fcn; };
  TFitter* GetMinuit() { return m_Minuit; };

  // --------------------------------------------------------
  // Pull Minimizer Methods
  // So, because of the inefficient way the code is structured,
  // this class acts as the interface between the calling analysis class and both the pull/minuit minimizers.
  // than the indexing of pull terms in the minimizer. because of this, we can't let the analysis class
  // talk to the pull minimizer directly, because the indexes for the pull terms probably won't match.
  // consequently, we have the following go between functions defining what and how an outside class can
  // get information from the pull minimizer.
private:
  MinimizeErrors* m_pull_min;
  bool fPullMinimizerInitialized;
  int fNVariedNewtonPars;
  bool fLastMinimizerCallConverged;
  double fNewtonSolverTolerance;
  int m_tries;
  void InitializePullMinimizer();
protected:
  bool fUseF2ij;
public:
  MinimizeErrors* GetPullMinimizer() { return m_pull_min; };
  double GetChiSquared();
  double GetPullChiSquared();
  void SetBinLabel( int binnum, std::string name );
  void SetBinNExpected( int binnum, double expected, double* gradient=NULL, int ngradpars=0 );
  void SetBinNObserved( int binnum, double observed );
  void SetFij( std::string pulltermName, int binNum, double fij );
  void SetFij( int parIndex, int binNum, double fij );
  void PrintFij();
  void PrintParameterInfo() { GetParManager()->Print(); };
  void SetNewtonSolverTolerance( double tolerance ) { fNewtonSolverTolerance = tolerance; };

  // --------------------------------------------------------

  // Functions for passing chi-squared gradient to MINUTI (for MINUIT terms only)
protected:
  bool fGiveMinuitGradient;
  int m_grad_nelements;
  int m_grad_nbins;
  int m_grad_npars;
  double* bin_gradients; // [ibin][ipar]
  void BuildGradientArray( int nbins, int npars );
  void DestroyGradientArray();
  double GetGradientElement( int ibin, int ipar );
  void SetGradientElement( int ibin, int ipar, double grad );
  double* GetBinGradient( int ibin );
  
public:
  virtual double GetChiSquaredPartial( int parid );
  bool IsGradientSentToMinuit() { return fGiveMinuitGradient; };
  void SendGradientToMinuit( bool doit );
  

  // Output via TTree
protected:
  double m_lastChiSquared;
  double m_lastPreFitChiSquared;
  double m_lastPreFitPullChiSquared;
  double m_lastFitterChiSquared;
  double m_lastFitterPullChiSquared;
  double m_lastUserPenalty;
  int m_converged_flag;

  int m_minuit_converged_flag;
  int m_newton_converged_flag;
public:
  virtual void StoreLastPrefitChiSquared( double chi ) { m_lastPreFitChiSquared = chi; };
  virtual void StoreLastUserPenalty( double penalty ) { m_lastUserPenalty = penalty; };
  virtual void AddFitResultBranchesToTree( TTree* tree );
  double GetLastUserPenalty() { return m_lastUserPenalty; };
  int DidLastCallConverge() { return m_converged_flag; };
  int DidLastMinuitCallConverged() { return m_minuit_converged_flag; };
  int DidLastPullMinimizerCallConverge() { return m_newton_converged_flag; };
  double GetPrefitChiSqaured() { return m_lastPreFitChiSquared; };
  double GetStoredChiSqaured() { return m_lastFitterChiSquared; };
  double GetStoredUserPenalty() { return m_lastUserPenalty; };
  double GetStoredPullChiSquared() { return m_lastFitterPullChiSquared; };

};

#endif
