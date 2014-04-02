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
 * \class ParameterManager
 *  \defgroup analysistools
 *  \brief Storage/Management of sys. error terms for fit.
 *
 * Note that we have a slightly expanded notion of what it means to be a systematic error term.
 * Actually, the word parameter is more fitting, but to change the name would require a lot of housekeeping 
 * we avoid for now, until the class has been rewritten.
 *
 * Responsibilities
 *  storage/registration/lookup for sys. terms
 *  management of covariance matrix between terms
 *  random parameter generation
 *  methods to get/set parameter vector
 *
 * Want the ability to also draw a random vector based on covariance matrix. 
 * -----------------------------------------------------------------------------------------------------------------------
 */

#ifndef __ParameterManager__
#define __ParameterManager__


#include <string>
#include <map>
#include <set>
#include <assert.h>

#include "TMatrixDSym.h"
#include "TMatrixD.h"
#include "TMatrixDSymEigen.h"
#include "TVector.h"
#include "TH2D.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TFile.h"

#include "ParamGen.hh"
#include "ModelParameter.hh"
#include "BasicParameter.hh"

namespace qosc {

  class ParameterManager {

  public:
    // ---------------------------------------------------------------------------
    // Constructors/Destructors

    ParameterManager();
    ParameterManager(std::string manager_name);
    virtual ~ParameterManager();

    // ---------------------------------------------------------------------------
    // Parameter Interface functions
  public:
    virtual void Initialize(); ///< Initializes Manager ... main job is to index parameters
    void RegisterParameter( ModelParameter* parterm ); ///< Register a systermatic error term
    void RegisterParameter( BasicParameter* parterm ); ///< Register a systermatic error term

    int  NumberOfMinuitPars() { return m_minuit_pars.size(); };
    int  NumberOfNewtonPars() { return m_newton_pars.size(); };
    int  NumberOfActivePars() { return m_parameter_list.size(); };
    int  NumberOfTotalPars() { return m_parameter_dict.size(); };
    int  NumberOfParameters() { return NumberOfActivePars(); };
    int  GetNumVariedNewtonPars();
    int  GetNumVariedMinuitPars();

    ModelParameter* GetParameter( std::string parname );
    std::string GetParameterTypeName( std::string parname );
    double GetParameterValue( std::string parname );
    void  SetParameterValue( std::string parname, double );
    bool  IsParameterDefined( std::string parname );
    void  SetParameterActiveStatus( std::string parname, bool active );
    bool  GetParameterActiveStatus( std::string parname );
    void  GetListOfParameterNames( std::vector< std::string >& parlist );
    void  GetListOfNewtonParNames( std::vector< std::string >& parlist );
    void  GetListOfMinuitParNames( std::vector< std::string >& parlist );
    void  GetParValueArray( double* );
    void  SetParValueByArray( double* parvalues );
    void  DefineParameterOrdering( std::string par_order ); ///< User can specify the ordering the systermatic error terms for whatever reasons. separate by ;:,
    void  SetToCentralValues();
    void  ZeroParameterValues();
    void  ZeroNewtonParValues();
    void  ZeroMinuitParValues();
    void  DefaultMinuitParValues();
    void  AddParameterValuesToTree( TTree* tree );
    void  AddListOfParametersToFile( TFile* file );

  protected:

    // ---------------------------------------------------------------------------
    // Instance Management

    std::string m_instance_name;
    int m_instanceID;
    void Start();

    // ---------------------------------------------------------------------------
    // Parameter Term Management

    // Dictionary
    //protected:
  public:
    std::map< std::string, ModelParameter* > m_parameter_dict; ///< Parameter stores important parameters and the Fij histograms associated with each sys term
    typedef std::map< std::string, ModelParameter* >::iterator ParDictIter;
  public:
    ParDictIter ParDictBegin() { return m_parameter_dict.begin(); };
    ParDictIter ParDictEnd() { return m_parameter_dict.end(); };

    // C omplete List (of active sys terms). Built when Initialize/Update is called.
  protected:
    std::vector< std::string > m_parameter_list;
    bool fParameterListReady;
  public:
    typedef std::vector< std::string >::iterator ParListIter;
    ParListIter ParListBegin() { 
      if  ( !fParameterListReady ) assert(false);
      return m_parameter_list.begin(); 
    };
    ParListIter ParListEnd() { 
      if ( !fParameterListReady ) assert(false);
      return m_parameter_list.end(); 
    };
    bool IsParameterListBuilt() { return fParameterListReady; };

    // --------------------------------------------------------
    // Parameter Indexing. Book keeping for the fitter.
  protected:
    std::map< std::string, int > m_parameter_index;
    std::map< int, std::string > m_parameter_from_index;
    ModelParameter** m_parameter_array; /// array of pointers. to avoid lookup table.
    std::string* m_parameter_name_array; /// array of strings
  public:
    int GetParameterIndex( std::string parname );
    ModelParameter* GetParameterFromIndex( int index );
    std::string GetParameterNameFromIndex( int index );

    // ---------------------------------------------------------------------------
    // Parameter group by type
  protected:
    // all of these get filled when initialize is called.
    std::vector< std::string > m_minuit_pars;
    std::vector< std::string > m_newton_pars;
    std::map< std::string, int > m_minuit_index; ///< id parameter number in the minuit class
    std::map< int, std::string > m_minuit_parname_fromindex; ///< id parameter number in the minuit class
    std::map< std::string, int > m_newton_index; ///< id parameter number in the minimizeerror class (the newton fitter)
    std::map< int, std::string > m_newton_parname_fromindex; ///< id parameter number in the minimizeerror class (the newton fitter)
  public:
    int GetMinuitIndex( std::string parname ) { return m_minuit_index[parname]; };
    int GetNewtonIndex( std::string parname ) { return m_newton_index[parname]; };
    std::string GetParameterNameByMinuitIndex( int index ) { return m_minuit_parname_fromindex[index]; };
    std::string GetParameterNameByNewtonIndex( int index ) { return m_newton_parname_fromindex[index]; };
  public:
    ParListIter NewtonParListBegin() { return m_newton_pars.begin(); };
    ParListIter NewtonParListEnd() { return m_newton_pars.end(); };
    ParListIter MinuitParListBegin() { return m_minuit_pars.begin(); };
    ParListIter MinuitParListEnd() { return m_minuit_pars.end(); };

    // ---------------------------------------------------------------------------
    // Parameter Term Group Management
    // -- we define groups in order to perform variations for subsets of systematice error terms
    // void DefineSysTermGroup( std::string groupname, std::string systerms );
    //std::map< std::string, std::set< std::string > > m_systerm_group_dict;

    // ---------------------------------------------------------------------------
    // Covariance Matrix Management
    // We primarily store covariance matrix information in the m_systerms_covariances dictionary.
    // We then build a covariance matrix once everything has been specified.
    // We do this so we don't have to rely on implementing the details of any input covariance matrices.
    // However, this strategy comes with the cost of extra boilerplate code for transferring cov matrix.
    // In future, might have a utility function to try and take care of this.
  protected:
    bool fCovMatrixReady;
    TMatrixDSym* m_cov_matrix;
    TMatrixDSym* m_invcov_matrix;
    TH2D* m_cov_matrix_hist;
    TMatrixD* m_cov_eigenmatrix;
    TMatrixD* m_cov_Teigenmatrix;
    TVectorD* m_cov_eigenvalues;
    void BuildCovarianceMatrix();

    typedef std::pair<std::string,std::string> ParameterPair;
    typedef std::map< ParameterPair, double > ParameterCovIter;
    std::map< ParameterPair, double > m_parameter_covariances;
    std::map< std::string, double > m_parameter_values_at_min; ///< likely obsolete
    std::set< std::string > m_parameter_wcov_term;
    bool DoesTermHaveCovariance( std::string parname ) { 
      std::set< std::string >::iterator it = m_parameter_wcov_term.find( parname );
      if ( it!=m_parameter_wcov_term.end() ) return true;
      else return false;
    };
  public:  
    void StoreParameterCentralValue( std::string term, double value_at_min );
    void SetParameterCovariance( std::string term1, std::string term2, double cov );
    double GetParameterCovariance( std::string term1, std::string term2 );
    TMatrixDSym* GetCovarianceMatrix() { return m_cov_matrix; };
    TH2D* GetCovarianceMatrixHist() { return m_cov_matrix_hist; };
    TMatrixDSym* GetInverseCovarianceMatrix() { return m_invcov_matrix; };

    // ---------------------------------------------------------------------------
    // Random Parameter term generation (needs to be rethought)
  protected:
    TRandom3* m_rand_gen;
    int m_rand_seed;
    ParamGen* m_param_gen;
  public:
    void GenerateRandomValues( double* values, bool rethrow_if_unphysical=false, int maxtries=1000 );
    void GenerateRandomValues( std::map< std::string, double>& values_dict, bool rethrow_if_unphysical=false, int maxtries=1000 );
    void SetRandomSeed( int seed ) {
      m_rand_seed = seed;
      delete m_rand_gen;
      m_rand_gen = new TRandom3( seed );
      if ( m_param_gen ) m_param_gen->SetSeed( seed );
    };

    // ---------------------------------------------------------------------------
    // Interface to TTree
    TTree* m_tree;
    bool fTreeBranchesSetup;
    double* m_parValueBranchBuffer;

    // ---------------------------------------------------------------------------
    // Misc and Utils
  public:
    void Print();
    void PrintParameterSummaryLine( ModelParameter* systerm );
    // ---------------------------------------------------------------------------


  private:

    static int NumInstances;

  };

}

#endif 
