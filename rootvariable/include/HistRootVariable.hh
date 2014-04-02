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
 * ------------------------------------------------------------------------------------------------
 * \class HistRootVariable
 * \ingroup PDFLibrary
 * \brief PDF which is designed to be built using data from a ROOT TTree/TChain
 *
 * Instances of this abstract class are suppose to be a representation of a quantity
 * derived from data in a ROOT Tree.
 * There are three components to specifying this class:
 *  (1) The variable: can eighter be a formula or a function (defined by a user made function pointer)
 *  (2) The weight: a weight class is defined via a formula, function or class to weight each event in the tree
 *  (3) The cut: a selection is defined via formula or function
 * 
 * For functions, one must defined a function with the following template:
 * double MyUserFunction( RootVariableList& );
 * A root variable list must be defined in conjuction with the function. At each event, the list
 * will be loaded with the variables so that the function can calculated the desired quantity.
 *
 * -----------------------------------------------------------------------------------------------
 */

#ifndef __HistRootVariable__
#define __HistRootVariable__

#include <string>
#include <map>
#include "Hist.hh"
#include "RootVariableList.hh"
#include "RootVariableFormula.hh"
#include "RootVariableFunction.hh"

class TH1D;
class TChain;


namespace qosc {

  class Weight;
  class RootVariable;

  class HistRootVariable : public Hist {

  public:
  
    // Constructor options, so many combinations
    // --------------------------------------------------------
    // 1D Histograms
    HistRootVariable( std::string variable_name, std::string variable_formula, TChain* source_tree ); // ( kTreeVarFormula, kNoWeight, kNoSelection): weight and selection has to be implented by user in the GetProbability function.
    HistRootVariable( std::string variable_name, std::string variable_formula, 
		      std::string weight_formula, TChain* source_tree ); // ( kTreeVarFormula, kWeightTreeFormula, kNoSelection): selection has to be implented by user in the GetProbability function. 
    HistRootVariable( std::string variable_name, std::string variable_formula,
		      std::string weight_formula, std::string selection_formula, TChain* source_tree ); // ( kTreeVarFormula, kWeightTreeFormula, kSelectionFormula)
    HistRootVariable( std::string variable_name, double (*VariableFunction)(RootVariableList&), std::string function_variables,
		      std::string weight_formula, std::string selection_formula, TChain* source_tree ); // ( kVariableFunction, kWeightTreeFormula, kSelectionFormula)
    HistRootVariable( std::string variable_name, std::string variable_formula, 
		      Weight* theEventWeight, std::string selection_formula, TChain* source_tree ); // ( kVariableFormula, kWeightClass, kSelectionFormula)
    HistRootVariable( std::string variable_name, double (*VariableFunction)(RootVariableList&), std::string function_variables,
		      Weight* theEventWeight, std::string selection_formula, TChain* source_tree ); // ( kVariableFunction, kWeightClass, kSelectionFormula)
    HistRootVariable( std::string variable_name, double (*VariableFunction)(RootVariableList&), std::string fill_function_variables,
		      Weight* theEventWeight, double (*SelectionFunction)(RootVariableList&), std::string cut_function_variables, TChain* source_tree ); // ( kVariableFunction, kWeightClass, kSelectionFunction)

    // --------------------------------------------------------
    // 2D Histograms
    HistRootVariable( std::string variable_name, 
		      std::string variable_formula_X, std::string variable_formula_Y,
		      std::string weight_formula,
		      std::string selection_formula, 
		      TChain* source_tree ); // ( All formula based )
    HistRootVariable( std::string variable_name, 
		      double (*VariableFunctionX)(RootVariableList&), std::string function_variables_X, double (*VariableFunctionY)(RootVariableList&), std::string function_variables_Y,
		      std::string weight_formula, 
		      std::string selection_formula, 
		      TChain* source_tree ); // (kVariableFunction , kWeightFormula, kSselectionFormula )
    HistRootVariable( std::string variable_name, 
		      double (*VariableFunctionX)(RootVariableList&), std::string function_variables_X, double (*VariableFunctionY)(RootVariableList&), std::string function_variables_Y,
		      Weight* theEventWeight, 
		      std::string selection_formula, 
		      TChain* source_tree ); // (kVariableFunction , kWeightClass, kSselectionFormula )
    HistRootVariable( std::string variable_name, 
		      std::string variable_formula_X, std::string variable_formula_Y, 
		      Weight* theEventWeight, 
		      std::string selection_formula, 
		      TChain* source_tree ); // ( kVariableFormula, kWeightClass, kSelectionFormula)
    HistRootVariable( std::string variable_name, 
		      double (*VariableFunctionX)(RootVariableList&), std::string function_variables_X,
		      double (*VariableFunctionY)(RootVariableList&), std::string function_variables_Y,
		      Weight* theEventWeight, 
		      double (*SelectionFunction)(RootVariableList&), std::string cut_function_variables,
		      TChain* source_tree ); // ( kVariableFunction, kWeightClass, kSelectionFunction)
  

    // --------------------------------------------------------
    // Copy Constructor
    HistRootVariable( std::string pdf_copyname, HistRootVariable* origvar ); ///< Copy constructor
  
    virtual ~HistRootVariable();
    virtual void BuildPDF();
    virtual void ProcessEvent( int entry );
    virtual double GetProbability( double observable ); ///< Required Probability Function
    virtual int GetCutResult( int entry ); ///< Get flag denoting if event passed all cuts
    virtual double GetWeight( int entry ); ///< Get weight of event
    virtual double GetVariableValue( int entry, std::string var="X" ); ///< Get Variable value
    void SetWeightActive( std::string weight_name, bool isactive);
    void SetCutActive( std::string cut_name, bool isactive);
    void SetNormalizeFlag( bool normme ) { fNormalizeHist = normme; };
    void NormalizeHistogram();
    void SetScanMode( bool scanmode ) { fScanMode = scanmode; };
    void SetScanAllCuts( bool scanmode ) { fScanAllCuts = scanmode; };
    void SetScanEntries( int entries ) { fScanEntries = entries; };
    void AddSelection( std::string cut );
    void ClearHistogram();

    // Here one can add adition weights and cuts
    void LoadWeightClass( std::string weight_class_name, Weight* weightclass ) { 
      m_weight_dict[weight_class_name] = weightclass; 
      fWeightDefinitionType = kWeightClass;
    };
    void AddSelectionFunction( double (*SelectionFunction)(RootVariableList&), std::string cut_variables );

    std::vector< std::string > GetCutFormulas() { return m_str_selection_formulas; }; // deep copy get
    std::string GetVariableName();
    std::string GetVariableFormula();
    TChain* GetSourceChain() { return m_source_tree; };


    // Status Flags: For checks and flow control
    enum VariableDefinitionType { kTreeVarFormula, kVariableFunction }; // Specifies the tree different ways I can specify a variable
    enum WeightDefinitionType { kNoWeight, kWeightClass, kWeightTreeFormula }; // The different ways I can specify the event weight
    //enum SelectionDefinition { kNoSelection, kSelectionFunction, kSelectionFormula };  // The different ways I can specify the selection function
    enum SelectionDefinition { kNoSelection, kSelectionFormula };  // The different ways I can specify the selection function
    enum WeightStatuses { kActive, kInactive };
    HistRootVariable::VariableDefinitionType GetVariableDefinitionType() { return fVariableDefinitionType; };
    HistRootVariable::WeightDefinitionType  GetWeightDefinitionType() { return fWeightDefinitionType; };
    HistRootVariable::SelectionDefinition GetSelectionType() { return fSelectionType; };


  protected:

    TChain* m_source_tree; // pointer to the ROOT tree form which we read the variable's values. We do not own it.  
    RootVariable* m_theRootVariable; ///< the RootVariable Object that we fill for each ROOT tree entry (may be formula or function)
    RootVariable* m_theRootVariableY; ///< the RootVariable Object used for 2D histograms. The Y variable.

    void CommonInitializationTasks( TChain* );
    void LoadVariable( std::string variable_name, std::string var="X" ); ///< creates the RootVariable as RootVariableFormula type
    void LoadVariableFunction( double (*VariableFunction)(RootVariableList&), std::string function_variables, std::string var="X" ); ///< creates the RootVariable as the RootVariableFunction type
    void LoadWeights( std::string weight_names ); // creates a list of RootVariables that will act as weights
    void LoadSelectionCuts( std::string selection_formulas ); // creates a list of RootVariables that will act as cuts
    RootVariable* GetRootVariable() { return m_theRootVariable; };
    RootVariableList* GetWeightRootVarList() { return &m_weights; };
    RootVariableList* GetSelectionRootVarList() { return &m_cuts; };
    void CopyWeightList( std::map< std::string, Weight* >& weightCopy );
    std::string GetCanonicalSelectionFunctionName() { return GetName()+"__selection_function__"; };

  public:
    bool IsWeightClassDefined( std::string weightname );
    Weight* GetWeightClass( std::string name );

    // Variable Formula
  protected:
    std::string m_variable_formula; ///< 1D variable, or X for 2D histograms
    std::string m_variable_formulaY; ///< 2D variable, or Y for 2D histograms

    // Weight Variables
    std::map< std::string, Weight* > m_weight_dict;
    std::vector< std::string > m_str_weight_formulas;
    std::map< std::string, int > m_weight_status;
    RootVariableList m_weights;

    // Select Formulas
    std::vector< std::string > m_str_selection_formulas;
    std::map< std::string, int > m_cut_status;
    RootVariableList m_cuts;
  
  private:


    VariableDefinitionType fVariableDefinitionType;
    WeightDefinitionType fWeightDefinitionType;
    SelectionDefinition fSelectionType;

    bool fScanMode;
    int fScanEntries;
    bool fScanAllCuts;
    bool fNormalizeHist;
  
  };
}

#endif
