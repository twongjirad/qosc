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

#include "HistRootVariable.hh"
#include <iostream>
#include <assert.h>
#include <cstdlib>
#include "TH1D.h"
#include "TH2D.h"
#include "TChain.h"
#include "Weight.hh"
#include "WeightRootVariables.hh"
#include "RootVariable.hh"
#include "RootVariableList.hh"
#include "RootVariableManager.hh"

using namespace qosc;

// --------------------------------------------------------
// 1D Histogram
// --------------------------------------------------------

HistRootVariable::HistRootVariable( std::string variable_name, std::string variable_formula, TChain* source_tree ) : Hist(variable_name) {
  // Specify 1D histogram using variable formula, no weight used, no selection used.
  CommonInitializationTasks( source_tree );
  fVariableDefinitionType = kTreeVarFormula;
  fWeightDefinitionType = kNoWeight;
  fSelectionType = kNoSelection;
  LoadVariable( variable_formula );
}

HistRootVariable::HistRootVariable( std::string variable_name, std::string variable_formula, std::string weight_formula, TChain* source_tree ) : Hist(variable_name) {
  // Specify 1D histogram using variable formula, weight formula, and no selection used.
  CommonInitializationTasks( source_tree );
  fVariableDefinitionType = kTreeVarFormula;
  fWeightDefinitionType = kWeightTreeFormula;
  fSelectionType = kNoSelection;
  LoadVariable( variable_formula );
  LoadWeights( weight_formula );
}

HistRootVariable::HistRootVariable( std::string variable_name, std::string variable_formula, std::string weight_formula, std::string selection_formula, TChain* source_tree ) : Hist(variable_name) {
  // Specify 1D histogram using variable formula, weight formula, and selection formula
  CommonInitializationTasks( source_tree );
  fVariableDefinitionType = kTreeVarFormula;
  fWeightDefinitionType = kWeightTreeFormula;
  fSelectionType = kSelectionFormula;
  LoadVariable( variable_formula );
  LoadWeights( weight_formula );
  LoadSelectionCuts( selection_formula );
}

HistRootVariable::HistRootVariable( std::string variable_name, 
				  double (*VariableFunction)(RootVariableList&), std::string function_variables, 
				  std::string weight_formula, std::string selection_formula, TChain* source_tree ) : Hist(variable_name) {
  // 1D histograms using variable function for filling, a weight formula, and selection formula
  CommonInitializationTasks( source_tree );
  fVariableDefinitionType = kVariableFunction;
  LoadVariableFunction( VariableFunction, function_variables );
  if ( weight_formula!="" ) {
    fWeightDefinitionType = kWeightTreeFormula;
    LoadWeights( weight_formula );
  }
  else {
    fWeightDefinitionType = kNoWeight;
  }

  if ( selection_formula!="" ) {
    fSelectionType = kSelectionFormula;
    LoadSelectionCuts( selection_formula );
  }
  else {
    fSelectionType = kNoSelection;
  }

}

HistRootVariable::HistRootVariable( std::string variable_name, double (*VariableFunction)(RootVariableList&), std::string function_variables, 
				  Weight* eventWeight, std::string selection_formula, TChain* source_tree ) : Hist(variable_name) {
  // 1D histograms using variable function for filling, a weight class, and selection formula
  CommonInitializationTasks( source_tree );
  fVariableDefinitionType = kVariableFunction;
  LoadVariableFunction( VariableFunction, function_variables );

  if ( eventWeight!=NULL ) {
    fWeightDefinitionType = kWeightClass;
    LoadWeightClass( "weightclass", eventWeight );
  }
  else 
    fWeightDefinitionType = kNoWeight;

  if ( selection_formula!="" ) {
    fSelectionType = kSelectionFormula;
    LoadSelectionCuts( selection_formula );
  }
  else
    fSelectionType = kNoSelection;

}

HistRootVariable::HistRootVariable( std::string variable_name, 
				  double (*VariableFunction)(RootVariableList&), std::string function_variables, 
				  Weight* eventWeight, 
				  double (*SelectionFunction)(RootVariableList&), std::string cut_function_variables,
				  TChain* source_tree ) : Hist(variable_name) {
  // 1D histograms using variable function for filling, a weight class, and selection function
  CommonInitializationTasks( source_tree );
  fVariableDefinitionType = kVariableFunction;
  LoadVariableFunction( VariableFunction, function_variables );

  if ( eventWeight!=NULL ) {
    fWeightDefinitionType = kWeightClass;
    LoadWeightClass( "weightclass", eventWeight );
  }
  else 
    fWeightDefinitionType = kNoWeight;

  if ( SelectionFunction!=NULL ) {
    m_cuts.SetChain( source_tree );
    fSelectionType = kSelectionFormula; // no difference in action between cut formula and cut function
    AddSelectionFunction( SelectionFunction, cut_function_variables );
  }
  else
    fSelectionType = kNoSelection;
  
}

HistRootVariable::HistRootVariable( std::string variable_name, std::string variable_formula, Weight* eventWeight, std::string selection_formula, TChain* source_tree ) : Hist(variable_name) {
  // 1D histogram using a formula for the fill variable, weight class, and selection formula
  CommonInitializationTasks( source_tree );
  fVariableDefinitionType = kTreeVarFormula;
  LoadVariable( variable_formula );
  if ( eventWeight!=NULL ) {
    fWeightDefinitionType = kWeightClass;
    LoadWeightClass( "weightclass", eventWeight );
  }
  else
    fWeightDefinitionType = kNoWeight;
  
  if ( selection_formula!="" ) {
    fSelectionType = kSelectionFormula;
    LoadSelectionCuts( selection_formula );
  }
  else
    fSelectionType = kNoSelection;
}

// --------------------------------------------------------
// 2D Histogram
// --------------------------------------------------------

// All Formulas Constructor //
HistRootVariable::HistRootVariable( std::string variable_name, 
				  std::string variable_formula_x, std::string variable_formula_y,
				  std::string weight_formula, std::string selection_formula, TChain* source_tree ) : Hist(variable_name) {
  CommonInitializationTasks( source_tree );
  fVariableDefinitionType = kTreeVarFormula;
  LoadVariable( variable_formula_x, "X" );
  LoadVariable( variable_formula_y, "Y" );

  if ( weight_formula!="" ) {
    fWeightDefinitionType = kWeightTreeFormula;
    LoadWeights( weight_formula );
  }
  else
    fWeightDefinitionType = kNoWeight;

  if ( selection_formula!="" ) {
    LoadSelectionCuts( selection_formula );
    fSelectionType = kSelectionFormula;
  }
  else
    fSelectionType = kNoSelection;
  
}

// Function Variables, Weight Formula, Selection Formula //
HistRootVariable::HistRootVariable( std::string variable_name, 
				  double (*VariableFunctionX)(RootVariableList&), std::string function_variables_X, 
				  double (*VariableFunctionY)(RootVariableList&), std::string function_variables_Y, 
				  std::string weight_formula, std::string selection_formula, TChain* source_tree ) : Hist(variable_name) {
  // 2D histograms using variable functions, weight formula, and selection formula
  CommonInitializationTasks( source_tree );
  fVariableDefinitionType = kVariableFunction;
  LoadVariableFunction( VariableFunctionX, function_variables_X, "X" );
  LoadVariableFunction( VariableFunctionY, function_variables_Y, "Y" );

  if ( weight_formula!="" ) {
    fWeightDefinitionType = kWeightTreeFormula;
    LoadWeights( weight_formula );
  }
  else {
    fWeightDefinitionType = kNoWeight;
  }

  if ( selection_formula!="" ) {
    fSelectionType = kSelectionFormula;
    LoadSelectionCuts( selection_formula );
  }
  else {
    fSelectionType = kNoSelection;
  }

}

// Variable Functions, Weight Class, Selection Formula //
HistRootVariable::HistRootVariable( std::string variable_name, 
				  double (*VariableFunctionX)(RootVariableList&), std::string function_variables_X, 
				  double (*VariableFunctionY)(RootVariableList&), std::string function_variables_Y, 
				  Weight* eventWeight, std::string selection_formula, TChain* source_tree ) : Hist(variable_name) {
  // 2D histograms using variable functions for filling, a weight class, and selection formula
  CommonInitializationTasks( source_tree );
  fVariableDefinitionType = kVariableFunction;
  LoadVariableFunction( VariableFunctionX, function_variables_X, "X" );
  LoadVariableFunction( VariableFunctionY, function_variables_Y, "Y" );

  if ( eventWeight!=NULL ) {
    fWeightDefinitionType = kWeightClass;
    LoadWeightClass( "weightclass", eventWeight );
  }
  else 
    fWeightDefinitionType = kNoWeight;

  if ( selection_formula!="" ) {
    fSelectionType = kSelectionFormula;
    LoadSelectionCuts( selection_formula );
  }
  else
    fSelectionType = kNoSelection;

}

// Variable Formulas, Weight Class, Selection Formula //
HistRootVariable::HistRootVariable( std::string variable_name, 
				  std::string variable_formula_X, std::string variable_formula_Y, 
				  Weight* eventWeight, std::string selection_formula, TChain* source_tree ) : Hist(variable_name) {
  // 2D histogram using a formula for the fill variable, weight class, and selection formula
  CommonInitializationTasks( source_tree );
  fVariableDefinitionType = kTreeVarFormula;
  LoadVariable( variable_formula_X, "X" );
  LoadVariable( variable_formula_Y, "Y" );

  if ( eventWeight!=NULL ) {
    fWeightDefinitionType = kWeightClass;
    LoadWeightClass( "weightclass", eventWeight );
  }
  else
    fWeightDefinitionType = kNoWeight;
  
  if ( selection_formula!="" ) {
    fSelectionType = kSelectionFormula;
    LoadSelectionCuts( selection_formula );
  }
  else
    fSelectionType = kNoSelection;
}

// Variable Functions, Weight Class, Selection Functions
HistRootVariable::HistRootVariable( std::string variable_name, 
				  double (*VariableFunctionX)(RootVariableList&), std::string function_variables_X,
				  double (*VariableFunctionY)(RootVariableList&), std::string function_variables_Y,
				  Weight* eventWeight, 
				  double (*SelectionFunction)(RootVariableList&), std::string cut_function_variables,
				  TChain* source_tree ) : Hist( variable_name ) {

  CommonInitializationTasks( source_tree );
  fVariableDefinitionType = kVariableFunction;
  LoadVariableFunction( VariableFunctionX, function_variables_X, "X" );
  LoadVariableFunction( VariableFunctionY, function_variables_Y, "Y" );

  if ( eventWeight!=NULL ) {
    fWeightDefinitionType = kWeightClass;
    LoadWeightClass( "weightclass", eventWeight );
  }
  else 
    fWeightDefinitionType = kNoWeight;

  if ( SelectionFunction!=NULL ) {
    m_cuts.SetChain( source_tree );
    fSelectionType = kSelectionFormula;
    AddSelectionFunction( SelectionFunction, cut_function_variables );
  }
  else 
    fSelectionType = kNoSelection;
}

// --------------------------------------------------------
// Copy Constructor
// --------------------------------------------------------
HistRootVariable::HistRootVariable( std::string pdf_copyname, HistRootVariable* origvar ) : Hist( pdf_copyname ) {
  // The copy constructor: we must replicate the state of the original. A bit of a mess.
  CommonInitializationTasks( origvar->GetSourceChain() );

  // Copy the Variable Type and the variable
  fVariableDefinitionType = origvar->GetVariableDefinitionType();
  if (fVariableDefinitionType==kTreeVarFormula) {
    LoadVariable( origvar->GetVariableFormula() );
  }
  else if (fVariableDefinitionType==kVariableFunction) {
    RootVariableFunction* orig_funcvar = (RootVariableFunction*)origvar->GetRootVariable();
    typedef double (*varfunction)(RootVariableList&);
    varfunction myvarfunc;
    myvarfunc = (varfunction) orig_funcvar->GetVariableFunction(); // Oh C, why do you let me get away with this?
    std::string func_vars = orig_funcvar->GetRootVariablesString();
    LoadVariableFunction( myvarfunc, func_vars );
  }

  // Set Weight Type
  fWeightDefinitionType = origvar->GetWeightDefinitionType();
  origvar->CopyWeightList( m_weight_dict ); // only copies pointers.
  
  // Set Selection Type: again, should have a container copy constructor
  fSelectionType = origvar->GetSelectionType();
  if ( fSelectionType == kSelectionFormula ) {
    m_cuts.SetChain( origvar->GetSourceChain() );
    std::vector< std::string > cut_names;
    origvar->GetSelectionRootVarList()->GetNameList( cut_names );
    std::string cut_vars = cut_names.at(0);
    for (unsigned int i=1; i<cut_names.size(); i++)
      cut_vars += ";"+cut_names.at(i);
    LoadSelectionCuts( cut_vars );
  }
  else {
    m_cuts.SetChain( origvar->GetSourceChain() );
  }

  // copy histogram if already defined
  if ( origvar->GetHistogramStatus()==Hist::kDefined ) {

    
    // Copy Histogram
    Hist::CopyHistogram( pdf_copyname, origvar );

    // Also if histogram has aready been defined, copy over any UserBinInfo for the histogram
    Hist::CopyBinInfo( origvar );
  }
  else {
    std::cout << "Warning, Histogram not defined in original variable."<< std::endl;
  }

}

HistRootVariable::~HistRootVariable() {
  //std::cout << "Calling HistRootVariable destructor" << std::endl;
  delete m_theRootVariable;
}


void HistRootVariable::CommonInitializationTasks( TChain* source_tree ) {
  fScanMode = false;
  fScanAllCuts = false;
  fNormalizeHist = true;
  fScanEntries = 100;
  //m_hist_pdf = NULL;
  m_source_tree = source_tree;
  return;
}

int HistRootVariable::GetCutResult( int entry ) {
  
  int selection_cut = 1;

  if ( fSelectionType== HistRootVariable::kNoSelection ) return selection_cut;

  if ( fSelectionType == HistRootVariable::kSelectionFormula ) {
    for ( RootVariableListIter it=m_cuts.Begin(); it!=m_cuts.End(); it++ ) {
      if ( m_cut_status[(*it).first]==kActive ) {
	RootVariable* cutvar = (*it).second;
	selection_cut *= int(cutvar->Value()+0.4);
      }
    }
  }

  if (fScanMode) {
    std::cout << "* Entry: " << entry << " * cut=" << selection_cut << " * ";
    if (fScanAllCuts) {
      for ( RootVariableListIter it=m_cuts.Begin(); it!=m_cuts.End(); it++ ) {
	std::string cutname = (*it).first;
	if ( m_cut_status[cutname]==kActive ) {
	  RootVariable* cutvar = m_cuts.GetVariable( cutname );
	  std::cout << cutname << " = " << cutvar->Value() << " * ";
	}
      }
    }
    std::cout << std::endl;
  }
  
  return selection_cut;
  
}

double HistRootVariable::GetWeight( int entry ) {
  if ( fScanMode ) std::cout << "* Entry: " << entry;
    
  double weight = 1.0;
  if ( fWeightDefinitionType==kNoWeight ) {
    if ( fScanMode ) std::cout << " * weight=1.0 (no weight) *" << std::endl;
    return weight;
  }


  std::map< std::string, Weight* >::iterator it_w;
  for (it_w=m_weight_dict.begin(); it_w!=m_weight_dict.end(); it_w++)
    if ( (*it_w).second!=NULL ) weight *= ((*it_w).second)->CalculateWeight();

  if ( fScanMode ) {
    std::cout << "* weight=" << weight << " * ";
    for (it_w=m_weight_dict.begin(); it_w!=m_weight_dict.end(); it_w++)
      std::cout << "* " << (*it_w).first << "=" << ((*it_w).second)->CalculateWeight();
    std::cout << " * " << std::endl;
  }

  return weight;
}

Weight* HistRootVariable::GetWeightClass( std::string weight_class_name ) {
  if ( m_weight_dict.find( weight_class_name )!=m_weight_dict.end() )
    return m_weight_dict[ weight_class_name ];
  else 
    return NULL;
}


double HistRootVariable::GetVariableValue( int entry, std::string var ) {
  double value = 0.0;
  if ( fVariableDefinitionType==kTreeVarFormula || fVariableDefinitionType==kVariableFunction) {
    // This is the price you pay for a strongly typed language (speed up possible if cast to double only)
    // Can I turn this into a macro?
    RootVariable* theVar = ( var=="X" ? m_theRootVariable : m_theRootVariableY );
    std::string vartype = theVar->GetType();
    
    if (vartype!="Char_t") {
      //value = static_cast< RootVariable<double>*>(m_theRootVariable)->Value();
      value = theVar->Value();
    }
    else {
      std::cout << "Root Variable List does not recognize this variable type: " << vartype << std::endl;
      assert(false);
    }
  }
  return value;
}

void HistRootVariable::ProcessEvent( int entry ) {

  if ( GetVerbose()>=2 ) {
    std::cout << "-------------------------------------------------------" << std::endl;
    std::cout << "Sample " << GetName() << ": entry(" << entry << ")" << std::endl;
    fScanMode = true;
    fScanAllCuts = true;
  }
  // selection portion  ------------------------------------------------------------------------
  int selection_cut = GetCutResult( entry );
  if (selection_cut==0) return;
    
  // weight calculation goes here ------------------------------------------------------------------------
  double weight = GetWeight( entry );
   if ( weight==0 ) return;

  // fill portion goes here ------------------------------------------------------------------------
  int dims = Hist::GetHistDimensions();
  double valueX = 0;
  double valueY = 0;
  if ( dims==1 ) {
    valueX = GetVariableValue( entry, "X" );
    ((TH1D*)GetHistogram())->Fill( valueX, weight );
    if ( fScanMode ) std::cout << "* Entry: " << entry << " * value=" << valueX << " * weight=" << weight << " *" << std::endl;
  }
  else if ( dims==2 ) {
    valueX = GetVariableValue( entry, "X" );
    valueY = GetVariableValue( entry, "Y" );
    ((TH2D*)GetHistogram())->Fill(valueX, valueY, weight);
  }
  else {
    std::cout << "HistRootVariable::ProcessEvent(): unsupported histogram dimension=" << dims << std::endl;
    assert(false);
  }
  
  // process fill bin info  ------------------------------------------------------------------------
  if ( dims==1 ) {
    int binfilled = ((TH1D*)GetHistogram())->FindBin( valueX );
    Hist::ProcessBinInfo( binfilled, weight );
  }
  else if ( dims==2 ) {
    int binX = ((TH2D*)GetHistogram())->GetXaxis()->FindBin( valueX );
    int binY = ((TH2D*)GetHistogram())->GetYaxis()->FindBin( valueY );
    if ( !Hist::fUseOverUnderFlow && binX>0 && binY>0 ) 
      Hist::ProcessBinInfo( binX-1, binY-1, weight );
    else
      Hist::ProcessBinInfo( binX, binY, weight );
  }
  // process fill hist info  ------------------------------------------------------------------------

  Hist::ProcessHistInfo( weight );

  // ------------------------------------------------------------------------ end of fill portion
  if ( GetVerbose()>=2 ) {
    fScanMode = false;
    fScanAllCuts = false;
    std::cout << "-------------------------------------------------------" << std::endl;
  }
  
}

void HistRootVariable::BuildPDF() {
  if ( fHistogramStatus==kUndefined ) {
    std::cout << "Cannot build PDF. Histogram still undefined" << std::endl;
    assert(false);
  }

  unsigned long bytes = 1;
  int entry = 0;
  while ( bytes>0 ) {
    bytes = m_source_tree->GetEntry( entry );
    if (bytes==0) break;

    if (entry%10000==0) std::cout << "HistRootVariable building PDF, entry " << entry << std::endl;
    if (fScanMode && entry>fScanEntries) break;

    ProcessEvent( entry );

    entry++;

  }//end of while loop over events

  if (fNormalizeHist)
    NormalizeHistogram();

  SetBuildStatus( false );  
}

void HistRootVariable::NormalizeHistogram() {
  Hist::SetIntegral( GetHistogram()->Integral() );
  if (Hist::GetIntegral()!=0.0)
    GetHistogram()->Scale(1.0/Hist::GetIntegral());
}

std::string HistRootVariable::GetVariableName() {
  return m_theRootVariable->GetVariableName();
}

std::string HistRootVariable::GetVariableFormula() {
  RootVariableFormula* rformula = dynamic_cast<RootVariableFormula*>(m_theRootVariable);
  if ( rformula==NULL ) {
    std::cout << "Tried to get a formula from an instance that is not of type RootVariableFormula." << std::endl;
    assert(false);
  }
  return rformula->GetFormulaName();
}

void HistRootVariable::AddSelection( std::string cut ) {
  m_cuts.Add( cut );
  m_cut_status[ cut ] = kActive;
  fSelectionType = kSelectionFormula;
}

void HistRootVariable::LoadVariable( std::string variable_name, std::string var ) {
  if ( fVariableDefinitionType==kTreeVarFormula ) {
    if ( var=="X" ) {
      m_variable_formula = variable_name;
      m_theRootVariable = RootVariableManager::GetTheRootVarManager()->MakeRootVariableFormula( variable_name, m_source_tree );
    }
    else if ( var=="Y" ) {
      m_variable_formulaY = variable_name;
      m_theRootVariableY = RootVariableManager::GetTheRootVarManager()->MakeRootVariableFormula( variable_name, m_source_tree );
    }
    else 
      assert(false);
    if ( var=="X" && m_theRootVariable->IsVariableAnArray() ) {
      std::cout << "---------------------------------------------------------------------" << std::endl;
      std::cout << "HistRootVariable, " << variable_name << ", has been created using an array ROOT tree data member.  " << std::endl;
      std::cout << "  This is not allowed yet." << std::endl;
      std::cout << "---------------------------------------------------------------------" << std::endl;
      assert(false);
    }
    if ( var=="Y" && m_theRootVariableY->IsVariableAnArray() ) {
      std::cout << "---------------------------------------------------------------------" << std::endl;
      std::cout << "HistRootVariable, " << variable_name << ", has been created using an array ROOT tree data member.  " << std::endl;
      std::cout << "  This is not allowed yet." << std::endl;
      std::cout << "---------------------------------------------------------------------" << std::endl;
      assert(false);
    }
  }
}

void HistRootVariable::LoadVariableFunction( double (*VariableFunction)(RootVariableList&), std::string function_variables, std::string var ) {
  // Load function variable
  // implicitly the variable name and source chain has been defined.  Should be OK, as this is protected
  if (GetSourceChain()==NULL) {
    std::cout << "Attempted to load a function variable for HistRootVariable that had a NULL tree." << std::endl;
    assert(false);
  }
  //m_theRootVariable = new RootVariableFunction( GetName(), VariableFunction, function_variables, GetSourceChain() );
  // We pass the creation and destruction of variables to the RootVarManager
  if (var=="X") m_theRootVariable = RootVariableManager::GetTheRootVarManager()->MakeRootVariableFunction( VariableFunction, function_variables, m_source_tree );
  else if ( var=="Y" ) m_theRootVariableY = RootVariableManager::GetTheRootVarManager()->MakeRootVariableFunction( VariableFunction, function_variables, m_source_tree );
}

void HistRootVariable::LoadWeights( std::string weights ) {
  // multiple weights can be defined by separating with a ;
  // the final weight will equal all product of all weights

  // we do this by defining a WeightRootVariables class
  WeightRootVariables* rootweight = new WeightRootVariables( GetSourceChain(), weights );
  m_weight_dict[weights] = rootweight;
  m_weight_status[weights] = true;

}

bool HistRootVariable::IsWeightClassDefined( std::string weightname ) {
  std::map< std::string, Weight* >::iterator it = m_weight_dict.find( weightname );
  if ( it!=m_weight_dict.end() ) return true;
  else return false;
}

void HistRootVariable::LoadSelectionCuts( std::string selection_formulas ) {
  // multiple weights can be defined by separating with a ;
  // the final weight will equal all product of all weights
  m_cuts.SetChain( m_source_tree );
  m_cuts.Add( selection_formulas );

  // check if multiplicity of all weights is 1: no arrays!
  std::vector< std::string > cut_names;
  m_cuts.GetNameList( cut_names );
  for (unsigned int i=0; i<cut_names.size(); i++) {
    RootVariable* cut_var = m_cuts.GetVariable( cut_names.at(i) );
    if ( cut_var->IsVariableAnArray() ) {
      std::cout << "The cut variable, " << cut_names.at(i) << ", is an array.  This is not allowed." << std::endl;
      assert(false);
    }
    m_str_selection_formulas.push_back( cut_names.at(i) );
    m_cut_status[ cut_names.at(i) ] = kActive;
  }
  
}

void HistRootVariable::AddSelectionFunction( double (*SelectionFunction)(RootVariableList&), std::string cut_function_variables ) {
  // Add a root variable function to the cut mix
  m_cuts.Add( SelectionFunction, cut_function_variables, m_source_tree ); // creates a root variable tree
  RootVariable* cut_var = m_cuts.GetVariable( SelectionFunction, cut_function_variables );
  if ( cut_var->IsVariableAnArray() ) {
    std::cout << "The cut variable, " << cut_var->GetVariableName() << ", is an array.  This is not allowed." << std::endl;
    assert(false);
  }
}

void HistRootVariable::SetWeightActive( std::string weight_name, bool isactive ) {
  if (isactive) m_weight_status[weight_name] = kActive;
  else m_weight_status[weight_name] = kInactive;
}

void HistRootVariable::SetCutActive( std::string cut_name, bool isactive ) {
  if (isactive) m_cut_status[cut_name] = kActive;
  else m_cut_status[cut_name] = kInactive;
}


void HistRootVariable::ClearHistogram() {
  if (m_hist_pdf) {
    m_hist_pdf->Reset();
  }
}

double HistRootVariable::GetProbability( double observable ) {

  // returns probability density

  if (fHistogramStatus==Hist::kUndefined)  {
    std::cout << "HistRootVariable::GetProbability: the histogram has not been defined yet" << std::endl;
    assert(false);
  }
  int bin = m_hist_pdf->FindBin( observable );
  //double prob = m_hist_pdf->Interpolate( observable );
  double prob = m_hist_pdf->GetBinContent(bin);
  double binwidth  = m_hist_pdf->GetBinWidth( bin );
  prob /= binwidth;
  if ( fNormalizeHist==false ) {
    if ( m_hist_pdf->Integral()>0 ) {
      prob /= m_hist_pdf->Integral();
    }
    else {
      std::cout << "Ill-defined PDF. Integral=" << m_hist_pdf->Integral() << std::endl;
      assert(false);
    }
  }
  
  return prob;
}

void HistRootVariable::CopyWeightList( std::map< std::string, Weight* >& weightCopy ) {
  // Copy our weight list to the argument list
  std::map< std::string, Weight* >::iterator it_w;
  for (it_w=m_weight_dict.begin(); it_w!=m_weight_dict.end(); it_w++) {
    weightCopy[(*it_w).first] = (*it_w).second;
  }
}

