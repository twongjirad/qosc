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

#include "ParameterManager.hh"

#include <sstream>
#include <iostream>
#include <assert.h>
#include <cmath>

#include "TDecompLU.h"

#include "SeparateString.hh"

using namespace qosc;

int ParameterManager::NumInstances = 0;

ParameterManager::ParameterManager() {
  Start();
}

ParameterManager::ParameterManager( std::string manager_name ) {
  Start();
  m_instance_name = manager_name; // override instance name
}

ParameterManager::~ParameterManager() {
  if ( m_cov_matrix ) delete m_cov_matrix;
  if ( m_invcov_matrix ) delete m_invcov_matrix;
  if ( m_cov_eigenmatrix ) delete m_cov_eigenmatrix;
  if ( m_cov_Teigenmatrix ) delete m_cov_Teigenmatrix;
  if ( m_cov_eigenvalues ) delete m_cov_eigenvalues;
  if ( m_parValueBranchBuffer ) delete [] m_parValueBranchBuffer;
  delete m_rand_gen;
}

void ParameterManager::Start() {
  m_instanceID = ParameterManager::NumInstances;
  ParameterManager::NumInstances++;
  std::stringstream ss;
  ss.str("");
  ss << "ParameterManagerInstance_" << m_instanceID;
  m_instance_name = ss.str();
  fParameterListReady = false;

  // Covariance Matrix Classes
  fCovMatrixReady = false;
  m_cov_matrix = NULL;
  m_cov_matrix_hist = NULL;
  m_invcov_matrix = NULL;
  m_cov_eigenmatrix = NULL;
  m_cov_Teigenmatrix = NULL;
  m_cov_eigenvalues = NULL;

  m_tree = NULL;
  fTreeBranchesSetup = false;
  m_parValueBranchBuffer = NULL;

  m_parameter_dict.clear();
  m_parameter_list.clear();

  // Random parameter value generator
  m_rand_seed = 0;
  m_rand_gen = new TRandom3( m_rand_seed );

  m_param_gen = NULL;  
}

// ------------------------------------------------------------------------------
// -------------------
// Parameter Access
// -------------------

ModelParameter* ParameterManager::GetParameter( std::string parname ) {
  ParDictIter it = m_parameter_dict.find( parname );
  if ( it!=ParDictEnd() ) return (*it).second;
  return NULL;
}

bool ParameterManager::IsParameterDefined( std::string parname ) {
  ModelParameter* parterm = GetParameter( parname );
  if ( parterm ) return true;
  return false;
}

std::string ParameterManager::GetParameterTypeName( std::string parname ) {
  assert( IsParameterDefined( parname ) );
  return GetParameter( parname )->IsA();
}

double ParameterManager::GetParameterValue( std::string parname ) {
  if ( !IsParameterDefined( parname ) ) {
    std::cout << "ParameterManager::GetParameterValue() - Error. Parameter=" << parname << " not defined." << std::endl;
    assert( false );
  }
  return GetParameter( parname )->GetValue();
}

void ParameterManager::SetParameterValue( std::string parname, double value ) {
  assert( IsParameterDefined( parname ) );
  return GetParameter( parname )->SetValue( value );
}


void ParameterManager::GetListOfParameterNames( std::vector< std::string >& parlist ) {
  for (ParListIter it=ParListBegin(); it!=ParListEnd(); it++) {
    parlist.push_back( (*it) );
  }
}

void ParameterManager::GetParValueArray( double* parvalues ) {
  if ( parvalues==NULL ) return;
  for ( ParDictIter it=ParDictBegin(); it!=ParDictEnd(); it++ ) {
    parvalues[ (*it).second->GetID() ] = (*it).second->GetValue();
  }
}

void ParameterManager::SetParValueByArray( double* parvalues ) {
  if ( parvalues==NULL ) assert(false);
  for ( ParDictIter it=ParDictBegin(); it!=ParDictEnd(); it++ ) {
    //std::cout << "set " << (*it).first << " to " << parvalues[ (*it).second->GetID() ] << std::endl;
    (*it).second->SetValue( parvalues[ (*it).second->GetID() ] );
  }
}

void ParameterManager::GetListOfNewtonParNames( std::vector< std::string >& parlist ) {
  for (ParListIter it=ParListBegin(); it!=ParListEnd(); it++) {
    if ( GetParameter(*it)->GetFitterParType()==kNewtonTerm )
      parlist.push_back( (*it) );
  }
}

void ParameterManager::GetListOfMinuitParNames( std::vector< std::string >& parlist ) {
  for (ParListIter it=ParListBegin(); it!=ParListEnd(); it++) {
    if ( GetParameter(*it)->GetFitterParType()==kMinuitTerm )
      parlist.push_back( (*it) );
  }
}

int ParameterManager::GetNumVariedNewtonPars() {
  int nvaried = 0;
  std::vector< std::string > parlist;
  GetListOfNewtonParNames( parlist );
  for (std::vector< std::string >::iterator it=parlist.begin(); it!=parlist.end(); it++) {
    if ( GetParameter( *it )->IsVaried() ) nvaried++;
  }
  return nvaried;
}

int ParameterManager::GetNumVariedMinuitPars() {
  int nvaried = 0;
  std::vector< std::string > parlist;
  GetListOfMinuitParNames( parlist );
  for (std::vector< std::string >::iterator it=parlist.begin(); it!=parlist.end(); it++) {
    if ( GetParameter( *it )->IsVaried() ) nvaried++;
  }
  return nvaried;
}

void ParameterManager::SetParameterActiveStatus( std::string parname, bool active ) {
  ModelParameter* parterm = GetParameter( parname );
  if ( parterm ) {
    parterm->SetActiveFlag( active );
    // This needs to reset the status of the manager. We need to rebuild the systematic error term list along with the covariance matrix.
  }
  else {
    std::cout << "ParameterManager::SetParameterActiveStatus - WARNING! Tried to set the status of a parameter that is not defined, " << parname << "." << std::endl;
  }
}

bool ParameterManager::GetParameterActiveStatus( std::string parname ) {
  ModelParameter* parterm = GetParameter( parname );
  if ( parterm ) return parterm->IsActive();
  else {
    std::cout << "ParameterManager::SetParameterActiveStatus - ERROR! Tried to get the status of a parameter that is not defined, " << parname << "." << std::endl;
    assert(false);
  }
}

int ParameterManager::GetParameterIndex( std::string sysname ) {
  std::map< std::string, int >::iterator it = m_parameter_index.find( sysname );
  if ( it!=m_parameter_index.end() ) {
    return (*it).second;
  }
  else {
    return -1;
  }
}

ModelParameter* ParameterManager::GetParameterFromIndex( int index ) {
  return m_parameter_array[index];
  return GetParameter( m_parameter_from_index[index] );
}


std::string ParameterManager::GetParameterNameFromIndex( int index ) {
  return m_parameter_name_array[index];
//   std::map< int, std::string >::iterator it = m_parameter_from_index.find( index );
//   if ( it!=m_parameter_from_index.end() ) {
//     return (*it).second;
//   }
  assert(false);
  return "";
}

// ------------------------------------------------------------------------------
// ---------------------------
// Parameter Covariance Storage
// ---------------------------

void ParameterManager::StoreParameterCentralValue( std::string term, double value_at_min ) {
  term = "systerm_"+term;
  m_parameter_values_at_min[term] = value_at_min;
}

void ParameterManager::SetParameterCovariance( std::string term1, std::string term2, double cov ) {
  // we ignore terms that are independent
  if (cov==0.0) return;
    
  ParameterCovIter it;
  ParameterPair p1 (term1, term2);
  ParameterPair p2 (term2, term1);
  if ( m_parameter_covariances.find( p1 )!=m_parameter_covariances.end() || m_parameter_covariances.find( p2 )!=m_parameter_covariances.end() ) return; // already stored
  else m_parameter_covariances[ p1 ] = cov;
  m_parameter_wcov_term.insert( term1 );
  m_parameter_wcov_term.insert( term2 );
}

double ParameterManager::GetParameterCovariance( std::string term1, std::string term2 ) {
  ParameterCovIter it;
  ParameterPair p1 (term1, term2);
  ParameterPair p2 (term2, term1);
  if ( m_parameter_covariances.find( p1 )!=m_parameter_covariances.end() ) return m_parameter_covariances[p1];
  else if ( m_parameter_covariances.find( p2 )!=m_parameter_covariances.end() ) return m_parameter_covariances[p2];
  else {
    //std::cout << term1 << ", " << term2 << " not found." << std::endl;
    return 0.;
  }
}

void ParameterManager::BuildCovarianceMatrix() {
  if ( m_cov_matrix ) {
    delete m_cov_matrix;
    m_cov_matrix = NULL;
  }
  if ( m_invcov_matrix ) {
    delete m_invcov_matrix;
    m_invcov_matrix = NULL;
  }
  int nterms = NumberOfActivePars();

  if ( nterms==0 ) return;

  std::cout << "ParameterManager::BuildCovarianceMatrix, n terms used in the fit=" << nterms << ", terms with cov=" << m_parameter_wcov_term.size() << std::endl;
  m_cov_matrix = new TMatrixDSym( nterms );
  m_invcov_matrix = new TMatrixDSym( nterms );
  std::string cov_mat_name = "cov_matrix_hist_"+m_instance_name;
  m_cov_matrix_hist = new TH2D( cov_mat_name.c_str(), "Covariance Matrix", nterms, 0, nterms, nterms, 0, nterms );
  std::set< std::string >::iterator it_covterms;

  for (int i=0; i<nterms; i++) {
    std::string termi = GetParameterNameFromIndex( i );
    it_covterms = m_parameter_wcov_term.find( termi );
    if ( it_covterms == m_parameter_wcov_term.end() ) {
      // no covariance terms.
      for (int j=0; j<nterms; j++) {
	if ( i!=j )
	  (*m_cov_matrix)[i][j] = (*m_cov_matrix)[j][i] = 0.;
	else
	  (*m_cov_matrix)[i][i] = 1.0; // this is to make sure the matrix is invertible

	m_cov_matrix_hist->SetBinContent( i+1, j+1, (*m_cov_matrix)[i][j] );
      }
    }
    else {
      // covariance for term defined
      for (int j=i; j<nterms; j++) {
	std::string termj = GetParameterNameFromIndex( j );
	if (termi=="" || termj=="" ) {
	  std::cout << "Sys Term not stored properly: termi=" << termi << ", termj=" << termj << std::endl;
	  assert(false);
	}
	(*m_cov_matrix)[i][j] = (*m_cov_matrix)[j][i] = GetParameterCovariance( termi, termj );
	m_cov_matrix_hist->SetBinContent( i+1, j+1, (*m_cov_matrix)[i][j] );
	m_cov_matrix_hist->SetBinContent( j+1, i+1, (*m_cov_matrix)[j][i] );
      }

      // sanity check
      if ( (*m_cov_matrix)[i][i]==0 ) {
	std::cout << "ParameterManager::BuildCovarianceMatrix - ERROR. No self correlation?" << std::endl;
	assert(false);
      }
    }
    m_cov_matrix_hist->GetXaxis()->SetBinLabel( i+1, GetParameterNameFromIndex( i ).c_str() );
    m_cov_matrix_hist->GetYaxis()->SetBinLabel( i+1, GetParameterNameFromIndex( i ).c_str() );
  }
  
  // should do tests of covariance matrix.
  // Find Inverse -- if invertible, disaster.
  TDecompLU* inverter = new TDecompLU( nterms );
  inverter->SetMatrix( *m_cov_matrix );
  inverter->Invert( (*(TMatrixD*)m_invcov_matrix) );
  delete inverter;

  // Find the eigen value and eigen vectors
  TMatrixDSymEigen* eigensolver = new TMatrixDSymEigen( *m_cov_matrix );
  if ( m_cov_eigenmatrix ) m_cov_eigenmatrix = NULL;
  if ( m_cov_Teigenmatrix ) m_cov_Teigenmatrix = NULL;
  if ( m_cov_eigenvalues ) m_cov_eigenvalues = NULL;
  m_cov_eigenmatrix = new TMatrixD( eigensolver->GetEigenVectors() );
  m_cov_Teigenmatrix = new TMatrixD( eigensolver->GetEigenVectors() ); m_cov_Teigenmatrix->T();
  m_cov_eigenvalues = new TVectorD( eigensolver->GetEigenValues() );
  delete eigensolver;

  // Setup Parameter generator. By definition, nominal is set at 1.0.
  TVectorD offset( nterms );
  for (int n=0; n<nterms; n++)
    offset(n) = 1. + GetParameterFromIndex( n )->GetCentralValue();
  m_param_gen = new ParamGen( offset, (*m_cov_matrix ) ); 
  m_param_gen->SetSeed( m_rand_seed );

  fCovMatrixReady = true;
}

// ------------------------------------------------------------------------------

void ParameterManager::RegisterParameter( ModelParameter* par ) {
  
  if ( !par ) {
    std::cout << "ParameterManager::RegisterParameter - Passed NULL Parameter instance!" << std::endl;
    assert(false);
  }

  if ( IsParameterDefined( par->GetName() ) == false ) {
    m_parameter_dict[ par->GetName() ] = par;
    fCovMatrixReady = false;
    fParameterListReady = false;
  }
  else {
    std::cout << "ParameterManager::RegisterParameter - WARNING! Reregistering parameter with name '" << par->GetName() << "'" << std::endl;
  }
}

void ParameterManager::RegisterParameter( BasicParameter* par ) {
  
  if ( !par ) {
    std::cout << "ParameterManager::RegisterParameter - Passed NULL Parameter instance!" << std::endl;
    assert(false);
  }

  RegisterParameter( dynamic_cast< ModelParameter* >( par ) );

}

void ParameterManager::DefineParameterOrdering( std::string par_order ) {
  // Reset the list and index map
  m_parameter_list.clear();
  m_parameter_index.clear();
  m_minuit_index.clear();
  m_newton_index.clear();
  m_minuit_parname_fromindex.clear();
  m_newton_parname_fromindex.clear();

  // must rebuild the covariance matrix
  fCovMatrixReady = false;

  int parIndex = 0;
  int minuitIndex = 0;
  int newtonIndex = 0;
  std::vector<std::string> parlist;
  SeparateString( par_order, parlist );
  for ( std::vector<std::string>::iterator it=parlist.begin(); it!=parlist.end(); it++ ) {
    std::string parname = *it;
    if ( GetParameterActiveStatus( parname )==true ) {
      m_parameter_list.push_back( parname );
      m_parameter_index[parname] = parIndex;
      m_parameter_from_index[parIndex] = parname;
      GetParameter( parname )->SetID( parIndex );
      if ( GetParameter( parname )->GetFitterParType()==kMinuitTerm ) {
	m_minuit_index[parname] = minuitIndex;
	m_minuit_parname_fromindex[minuitIndex] = parname;
	minuitIndex++;
      }
      else if ( GetParameter( parname )->GetFitterParType()==kNewtonTerm ) {
	m_newton_index[parname] = newtonIndex;
	m_newton_parname_fromindex[newtonIndex] = parname;
	newtonIndex++;
      }
      parIndex++;
    }
    else {
      std::cout << "ParameterManager::DefineParameterOrdering - WARNING. " << parname << " is not active. Skipping." << std::endl;
    }
  }
  fParameterListReady = true;
}

void ParameterManager::Initialize() {
  std::cout << "ParameterManager[" << m_instance_name << "]::Initialize()" << std::endl;
  // Give each systerm in the list an index ID. 
  // The user may choose the order the parameter herself.
  // However, by default we order by:
  //   first getting a list of all the types of systerm parameters.
  //   then we order by type.
  if ( !fParameterListReady ) {

    // if the parameter order has not been defined, we organize the terms by types by default
    std::map< std::string, std::set< std::string >* > par_types;
    typedef std::map< std::string, std::set< std::string >* >::iterator par_typelist;

    // Loop through dict of systematic terms and get the types
    m_minuit_pars.clear();
    m_newton_pars.clear();
    for (ParDictIter it=ParDictBegin(); it!=ParDictEnd(); it++) {
      ModelParameter* par = (*it).second;
      if ( par==NULL || par->IsA()=="" ) {
	std::cout << "ParameterManager::Initialize - ERROR! Blank name given to Parameter Type. Please rename." << std::endl;
	std::cout << "address: " << par << std::endl;
	assert(false);
      }
      if ( par_types.find( par->IsA() )==par_types.end() )
	par_types[par->IsA()] = new std::set< std::string >;

      par_types[ par->IsA() ]->insert( par->GetName() );
    }

    // Now loop through the types and order the parameters
    std::string parnames = "";
    for ( par_typelist it=par_types.begin(); it!=par_types.end(); it++ ) {
      
      for ( std::set< std::string >::iterator itsys = (*it).second->begin(); itsys!=(*it).second->end(); itsys++ ) {
	parnames += (*itsys)+";";
      }
    
      // we are done with this list. destroy it.
      delete par_types[ (*it).first ];
      par_types[ (*it).first ] = NULL;
    }

    // for good measure, empty the dict/map
    par_types.clear();
    
    // Now send the parameter list to the function that will do the indexing.
    DefineParameterOrdering( parnames );
  }
  
  
  // (2) Make list of Newton and Minuit terms
  m_minuit_pars.clear();
  m_newton_pars.clear();
  for (ParDictIter it=ParDictBegin(); it!=ParDictEnd(); it++) {
    ModelParameter* par = (*it).second;
    // store minuit and newton parameter types
    if ( par->GetFitterParType()==kNewtonTerm ) m_newton_pars.push_back( par->GetName() );
    else if ( par->GetFitterParType()==kMinuitTerm ) m_minuit_pars.push_back( par->GetName() );
    else assert(false);
  }

  // (3) Build ModelParameter Pointer Array
  m_parameter_array = new ModelParameter*[ NumberOfTotalPars() ];
  m_parameter_name_array = new std::string[ NumberOfTotalPars() ];
  memset( m_parameter_array, 0, sizeof(ModelParameter*)*NumberOfTotalPars() );
  for (int n=0; n<NumberOfTotalPars(); n++) {
    m_parameter_array[n] = GetParameter( m_parameter_from_index[n] );
    m_parameter_name_array[n] = m_parameter_from_index[n] ;
  }

  // (4) Build Covariance Matrix
  BuildCovarianceMatrix();

}

void ParameterManager::SetToCentralValues() {
  for ( ParDictIter it=ParDictBegin(); it!=ParDictEnd(); it++ ) {
    (*it).second->SetValue( (*it).second->GetCentralValue()  );
  }
}

void ParameterManager::ZeroParameterValues() {
  for ( ParDictIter it=ParDictBegin(); it!=ParDictEnd(); it++ ) {
    (*it).second->SetValue( 0.0 );
  }
}

void ParameterManager::ZeroNewtonParValues() {
  for ( ParDictIter it=ParDictBegin(); it!=ParDictEnd(); it++ ) {
    if ( (*it).second->GetFitterParType()==kNewtonTerm )
      (*it).second->SetValue( 0.0 );
  }  
}

void ParameterManager::ZeroMinuitParValues() {
  for ( ParDictIter it=ParDictBegin(); it!=ParDictEnd(); it++ ) {
    if ( (*it).second->GetFitterParType()==kMinuitTerm )
      (*it).second->SetValue( 0.0 );
  }  
}

void ParameterManager::DefaultMinuitParValues() {
  for ( ParDictIter it=ParDictBegin(); it!=ParDictEnd(); it++ ) {
    if ( (*it).second->GetFitterParType()==kMinuitTerm )
      (*it).second->GetDefault();
  }  
}

void ParameterManager::GenerateRandomValues( double* values, bool rethrow_if_unphysical, int maxtries ) {
  // This should probably just return a dict of values or an array.
  // Instead it sets the values, which is not always good.
  
  int nelements = NumberOfActivePars();
  if ( values==NULL )
    values = new double[nelements];

  int ntries = 0;
  while ( rethrow_if_unphysical && ntries<maxtries ) {

    // draw a random vector
    TMatrixD pulls(nelements,nelements);
    for (int i=0; i<nelements; i++) {
      for (int j=0; j<nelements; j++) {
	if (i==j) pulls[i][i] = m_rand_gen->Gaus( 0, sqrt(fabs((*m_cov_eigenvalues)[i])) );
	else pulls[i][j] = 0.;
      }
    }
    
    // Now get back to undiagonalized basis                                                                                                                                       
    TMatrixD A( pulls, TMatrixD::kMult, *m_cov_Teigenmatrix );
    TMatrixD B( *m_cov_eigenmatrix, TMatrixD::kMult, A );
    
    // grab the diagonals for the pull values                                                                                                                                     
    
    for (int i=0; i<nelements; i++) {
      values[i] = pulls[i][i];
      if ( rethrow_if_unphysical ) {
	rethrow_if_unphysical = false; // default is to exit loop
	if ( GetParameterFromIndex(i)->GetLowBound()>values[i] || values[i]>GetParameterFromIndex(i)->GetHighBound() ) {
	  rethrow_if_unphysical = true;
	  ntries++;
	  std::cout << "ParameterManager::Must rethrow parameters (tries=" << ntries << ")" << std::endl;
	}
      }
    }
  } // end of while loop
  
}


void ParameterManager::GenerateRandomValues( std::map< std::string, double > &value_dict, bool rethrow_if_unphysical, int maxtries ) {
  std::vector< double > values; 
  bool parsOK = false;
  int ntries = 0;
  while ( parsOK==false && ntries<maxtries) {

    m_param_gen->ThrowSet( values );
    parsOK = true;

    for (unsigned int i=0; i<values.size(); i++) {
      
      if ( GetParameterFromIndex(i)->DoesTermThrow() ) {
	double parval = values.at(i)-1.0;
	value_dict[ GetParameterNameFromIndex(i) ] = parval;
	if ( rethrow_if_unphysical && ( parval<GetParameterFromIndex(i)->GetLowBound() || parval>GetParameterFromIndex(i)->GetHighBound() ) ) {
	  parsOK = false;
	  std::cout << "ParameterManager::GenerateRandomValues: parameter=" << GetParameterNameFromIndex(i) << " was gen out of bounds: "
		    << GetParameterFromIndex(i)->GetLowBound() << " < " << parval << " < " << GetParameterFromIndex(i)->GetHighBound() << std::endl;
	}
      }
      else {
	value_dict[ GetParameterNameFromIndex(i) ] = GetParameterFromIndex(i)->GetValue();
      }
    }
    if ( parsOK==false ) {
      ntries++;
      std::cout << "ParameterManager::GenerateRandomValues - par values unphysical. rethrowing parameters. tries=" << ntries << " of " << maxtries << std::endl;
    }

  }
  return;
  
}

void ParameterManager::AddParameterValuesToTree( TTree* tree ) {
  m_tree = tree; // store pointer for reference. unlikely to use it. maybe use to check if already stored.

//   if ( m_parValueBranchBuffer ) delete [] m_parValueBranchBuffer;
//   m_parValueBranchBuffer = new double[GetNumActivePars()];

  char parname[100];
  char branchinfo[100];
  for ( ParListIter it=ParListBegin(); it!=ParListEnd(); it++ ) {
    sprintf( parname, "%s", (*it).c_str() );
    sprintf( branchinfo, "%s", std::string( (*it) + "/D" ).c_str() );
    //m_tree->Branch( parname, &m_parValueBranchBuffer[ GetParameterIndex( *it ) ], branchinfo );
    double* p_par = GetParameter( *it )->GetParameterAddress();
    m_tree->Branch( parname, p_par , branchinfo );
  }
  fTreeBranchesSetup = true;
}

// void ParameterManager::FillParameterValuesToTree() {
//   if ( !fTreeBranchesSetup ) {
//     std::cout << "ParameterManager::FillParameterValuesToTree() - WARNING! Tree has not been setup." << std::endl;
//     return;
//   }

//   for (int n=0; n<NumberOfActivePars(); n++) {
//     m_parValueBranchBuffer[ n ] = GetParameterFromIndex(n)->GetValue();
//   }

//   m_tree->Fill();
//  }

void ParameterManager::AddListOfParametersToFile( TFile* file ) {
  
  std::cout << "ParameterManager::AddListOfParametersToFile()" << std::endl;
  
  file->cd();
  TTree* tree = new TTree( "param_list", "Parameter List" );
  char name[100];
  char partypename[100];
  int active_status;
  int index;
  double sigma;

  tree->Branch( "ParName", name, "ParName[100]/C" );
  tree->Branch( "ParTypeName", partypename, "ParTypeName[100]/C" );
  tree->Branch( "Index", &index, "Index/I" );
  tree->Branch( "Sigma", &sigma, "Sigma/D" );
  tree->Branch( "ActiveStatus", &active_status, "ActiveStatus/I" );
  
  for ( ParDictIter it=ParDictBegin(); it!=ParDictEnd(); it++ ) {
    sprintf( name, "%s", (*it).first.c_str() );
    sprintf( partypename, "%s", (*it).second->IsA().c_str() );
    if ( GetParameterActiveStatus( (*it).first ) ) {
      active_status = 1;
      index = GetParameterIndex( (*it).first );
      sigma = (*it).second->GetSigma();
    }
    else {
      active_status = 0;
      index = -1;
      sigma = 0.;
    }
    tree->Fill();
  }

  tree->Write();
  delete tree;

}


void ParameterManager::Print() {
  std::cout << "-----------------------------------------------------------------------" << std::endl;
  
  std::cout << "************ ParameterManager[" << this << "]::Print() *****************" << std::endl;
  std::cout << "Instance Name: " << m_instance_name << std::endl;
  std::cout << "Number of registered terms: " << m_parameter_dict.size() << std::endl;
  std::cout << "Number of indexed/active terms: " << m_parameter_list.size() << std::endl;
  std::cout << "-----------------------------------------------------------------------" << std::endl;
  if ( m_parameter_list.size()>0 ) {
    std::cout << "Indexed/Active Terms" << std::endl;
    std::cout << "[index] Name : Parameter Type : Fitter Type : Value" << std::endl;
    std::cout << "-----------------------------------------------------------------------" << std::endl;
    for ( unsigned int index=0; index<m_parameter_list.size(); index++ ) {
      std::string parname = GetParameterNameFromIndex( index );
      std::cout << " [" << index << "] ";
      ModelParameter* parterm = m_parameter_dict[ parname ];
      PrintParameterSummaryLine( parterm );
    }
  }
  if ( m_parameter_list.size()!=m_parameter_dict.size() ) {
    std::cout << "Unindexed/Inactive Terms" << std::endl;
    std::cout << "------------------------" << std::endl;
    for ( ParDictIter it=ParDictBegin(); it!=ParDictEnd(); it++ ) {
      if ( m_parameter_index.find( (*it).first )!=m_parameter_index.end() ) continue;
      std::string parname = (*it).first;
      ModelParameter* parterm = (*it).second;
      PrintParameterSummaryLine( parterm );
    }
  }
  std::cout << "------------------------------------------------------" << std::endl;
  
}

void ParameterManager::PrintParameterSummaryLine( ModelParameter* par ) {
  std::cout << par->GetName() << " : " << par->IsA() << " : ";
  if ( par->GetFitterParType()==kMinuitTerm ) std::cout << " Minuit";
  else if ( par->GetFitterParType()==kNewtonTerm ) std::cout << " Newton";
  std::cout << " : Value=" << par->GetValue()
	    << " Mean=" << par->GetCentralValue() << " (" << 1+par->GetCentralValue() << ")"
	    << " Sigma=" << par->GetSigma() << "(" << (1+par->GetCentralValue())*par->GetSigma() << ")";
  //if ( par->GetFitterParType()==kMinuitTerm && par->IsActive() ) {
  if ( par->IsActive() ) {
    if ( par->IsVaried() )
      std::cout << " VARIED";
    else
      std::cout << " FIXED";
  }
  std::cout << std::endl;  
}
