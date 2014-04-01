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
#include "RootVariableFormula.hh"

#include <iostream>
#include <sstream>
#include <assert.h>
#include <cstdlib>
#include <cstdarg>

#include "TChain.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TTreeFormula.h"
#include "TTreeFormulaManager.h"
#include "TList.h"

#include "RootVariableManager.hh"

using namespace qosc;

RootVariableFormula::RootVariableFormula( std::string formula, TChain* source_tree )  :
  RootVariable() {

  // initialize internal members
  m_var_dim_nelems = NULL;
  m_num_specified_dims = 0;
  RootVariable::SetArrayFlag( false ); 
  RootVariable::SetChain( source_tree );
  LoadVariable( formula );

  lastCalculatedEntry = -1;
  useCachedValue = false;
  cachedValue = 0.;

}

RootVariableFormula::RootVariableFormula( std::string variable_name, std::string formula, TChain* source_tree )  :
  RootVariable( variable_name, source_tree )
{
  
  // initialize internal members
  m_var_dim_nelems = NULL;
  m_num_specified_dims = 0;
  RootVariable::SetArrayFlag( false ); 
  LoadVariable( formula );
}

RootVariableFormula::~RootVariableFormula() {
  delete m_var_dim_nelems;
  delete m_formula;
  delete [] m_input_indices;
}

void RootVariableFormula::LoadVariable( std::string formula ) {
  // I am relying on the user specifying the array indices in the title.
  // I believe this is necessary for ROOT to properly pack the data into the Tree.
  // So its not such a bad system. Though, it's not the most flexible.  

  // Make the Formula
  m_formula_name = formula;
  AllocateVariable();
  //std::cout << "load variable: " << formula << std::endl;

  // Parse the dimensions (use TTreeFormulaMananger is better)
  // Now the multiplicity of the answer takes on the leaf with the largest multiplicity
  int nleaves = m_formula->GetNcodes();
  int max_multiplicity = -1;
  int max_leaf = 0;
  for (int nleaf = 0; nleaf<nleaves; nleaf++) {
    TLeaf* formula_leaf = m_formula->GetLeaf( nleaf );
    if ( max_multiplicity < formula_leaf->GetNdata() ) {
      max_leaf = nleaf;
      max_multiplicity = formula_leaf->GetLenStatic();
    }
    //std::cout << "enable: " << formula_leaf->GetName() << std::endl;
    //RootVariable::GetChain()->SetBranchStatus( formula_leaf->GetBranch()->GetName(), 1 );
    //RootVarManager::GetTheRootVarManager()->RegisterBranch( (TChain*)formula_leaf->GetBranch()->GetTree(), formula_leaf->GetBranch()->GetName() );
    // register this branch with the root variable manager.
    RootVariableManager::GetTheRootVarManager()->RegisterBranch( GetChain(), formula_leaf->GetBranch()->GetName() ); // would like to avoid this so not circular dependency
  }

  //TLeaf* leafcount = m_formula->GetLeaf(0)->GetLeafCount(); // Sneaky, ROOT, sneaky.

  if ( nleaves>1 && max_multiplicity!=m_formula->GetNdata() ) {
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "WARNING!!!" << std::endl;
    std::cout << "leaf multiplicity from " << m_formula->GetLeaf(max_leaf)->GetName() << " does not match the formulas multiplicity" << std::endl;
    std::cout << " my max multiplicity=" << max_multiplicity << std::endl;
    std::cout << " formula Ndata = " << m_formula->GetNdata() << std::endl;
    std::cout << " this means my multiplicity algorithm is wrong. Attempts to access the array elements will then be off. may require a fix." << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
    //assert(false);
  }

  // Find the Leaf
//   std::cout << "Variable name: " << RootVariable<T>::GetVariableName()
// 	    << " (Leaf Name: " << m_formula->GetLeaf(max_leaf)->GetName() << ") "
// 	    << " multiplicity=" << max_multiplicity 
// 	    << " leafcount=" << leafcount
// 	    << std::endl;

  m_var_leaf = RootVariable::GetChain()->FindLeaf( m_formula->GetLeaf(max_leaf)->GetName() );
  if (m_var_leaf==NULL) {
    std::cout << "The specified variable, " << m_formula->GetLeaf(max_leaf)->GetName() << ", was not found in this tree!" << std::endl;
    assert( false );
  }
  m_var_branch = m_var_leaf->GetBranch();

  // Get the Leaf Array dimensions form branch title
  std::string titlename = m_var_branch->GetTitle(); // Get's the Branch's title.  This should include array info
  ParseArrayBrackets( titlename, m_var_dim_nelems, m_var_num_dims );
  // allocate array for indices, which the Value call will use later
  m_input_indices = new int[m_var_num_dims];

  // Determine if this variable is an array
  RootVariable::SetArrayFlag( false );
  int numdims = 0;
  for (int n=0; n<m_var_num_dims; n++)
    numdims += m_var_dim_nelems[n];
  int ndata = m_formula->GetNdata();
  //  std::cout << "ndata=" << ndata << std::endl;
//   if (leafcount) {
//     ndata = leafcount->GetMaximum();
//     std::cout << "leaf data: " << ndata << std::endl;
//   }

  if (ndata>1 && numdims>1) RootVariable::SetArrayFlag( true );  // THIS CAN HAVE BUGGY BEHAVIOR FOR VAR LENGTH ARRAYS
  //std::cout << ndata << " " << numdims << " " << RootVariable<T>::IsVariableAnArray() << std::endl;

  if ( RootVariable::IsVariableAnArray() ) {
    // Get number of specified array dims from variable
    // first split variable name into operations
    std::string operations = "();*+-/ ";
    std::vector< std::string > vars_in_formula;
    size_t start = 0;
    size_t end = 0;
    while ( start!=std::string::npos && end!=std::string::npos ) {
      end = formula.find_first_of( operations, start );
      std::string sub;
      if (end!=std::string::npos)
	sub = formula.substr( start, end-start );
      else
	sub = formula.substr( start );
      if (sub.size()>0) {
	// check to make sure its not just a constant
	std::string notletter = "[]1234567890.";
	size_t pos = sub.find_first_not_of( notletter );
	if (pos!=std::string::npos)
	  vars_in_formula.push_back( sub );
      }
      start = end+1;
    }
    
    // Now determine the specified dims
    m_num_specified_dims = 0;
    std::vector< std::string >::iterator it;
    for (it=vars_in_formula.begin(); it!=vars_in_formula.end(); it++) {
      std::string var = *it;
      if (var==m_formula_name) {
	int* indices_specified = NULL;
	int num_indices = 0;
	ParseArrayBrackets( var, indices_specified, num_indices );
	//std::cout << "defining variable: " << var << ": indices=" << num_indices << std::endl;
	m_num_specified_dims+=num_indices;
      }
//       else {
// 	std::cout << "other variable: " << var << std::endl;
//       }
    }

  }
  
  //std::cout << "Setup variable, " << m_formula_name << " := " << RootVariable<T>::GetVariableName() << ", multiplicity=" << m_formula->GetNdata() << ", " << std::endl;
  
}

void RootVariableFormula::ParseArrayBrackets( std::string array_name, int*& array_info , int& num_array_elems) {
  // This routine looks for bracket pairs and stores the contents into the vector passed to us.
  // e.g. array_name = "myarray[3][4]" returns array_info = { 3, 4 } and num_array_elemns = 2
  
  size_t str_pos = 0;
  bool findopen = true;
  bool findclose = false;
  size_t openpos, closepos;
  std::vector<int> array_indices;
  while (str_pos!=std::string::npos) {
    if (findopen) {
      str_pos = array_name.find("[", str_pos, 1);
      openpos = str_pos+1; //+1 to not include the bracket
      findopen = false; findclose = true;
    }
    else if (findclose) {
      str_pos = array_name.find("]",str_pos,1);
      closepos = str_pos;
      if ( closepos==std::string::npos ) {
	std::cout << "Error parsing variable name.  Mismatched array brackets!" << std::endl;
	assert(false);
      }
      size_t len = closepos-openpos;
      std::string str_array_index = array_name.substr( openpos, len );
      std::string numbers = "0123456789";
      bool anint = true;
      for (int n=0; n<str_array_index.size(); n++) {
	size_t foundnum = numbers.find( str_array_index[n], 0 );
	if (foundnum==std::string::npos) {
	  anint = false;
	  //std::cout << "This array parameter is not a number: " << str_array_index << std::endl;
	  //std::cout << "Is it a leaf in our tree? ";
	  TLeaf* par = RootVariable::GetChain()->FindLeaf( str_array_index.c_str() );
	  if (par) {
	    //std::cout << " Yes.";
	    //std::cout << " To be safe, we need to find its max value:";
	    TLeaf* leafcount = m_var_leaf->GetLeafCount(); // Sneaky, ROOT, sneaky.
	    int maxvalue = leafcount->GetMaximum();
	    //std::cout << " " << maxvalue << std::endl;
	    array_indices.push_back( maxvalue );
	    m_flexdims.push_back( str_array_index );
	    m_flexdim_max_elems[ str_array_index ] = maxvalue;
	    if ( array_indices.size()>1) {
	      std::cout << "Sorry, but only the first dimension can have a variable array length.  I don't know how to align the array otherwise." << std::endl;
	      std::cout << "If this somehow becomes a problem, then email taritree.wongjirad@gmail.com, to address it" << std::endl;
	      assert(false);
	    }
	    if ( m_flexdims.size()>1 ) {
	      std::cout << "Sorry, but only the first dimension can have a variable array length.  I don't know how to align the array otherwise." << std::endl;
	      std::cout << "If this somehow becomes a problem, then email taritree.wongjirad@gmail.com, to address it" << std::endl;
	      assert(false);
	    }
	  }
	  else {
	    std::cout << " No. Unacceptable. Must give a way to override this. Maybe specify the parameter." << std::endl;
	    assert(false);
	  }
	  break;
	}
      }
      if ( anint ) {
	int index = std::atoi( str_array_index.c_str() );
	array_indices.push_back( index );
      }
      
      findopen = true; findclose = false;
    }
    
  }//end of while loop
  
  // now make index array
  num_array_elems = array_indices.size();
  array_info = new int[ num_array_elems ];
  for (int n=0; n<num_array_elems; n++) {
    if (n==0 && m_flexdims.size()>0)
      array_info[n] = m_flexdim_max_elems[ m_flexdims.at(n) ];
    else
      array_info[n] = int(array_indices.at(n));
  }
}

void RootVariableFormula::AllocateVariable() {

  if (  RootVariable::GetChain()->GetTree() == NULL )  RootVariable::GetChain()->GetEntry(0); // If the chain has not loaded the tree, load the first entry and tree. do i need to load the aliases?
  try {
    m_formula = new TTreeFormula( RootVariable::GetVariableName().c_str(), m_formula_name.c_str(), RootVariable::GetChain() );
  }
  catch (int e) {
    std::cout << "Error=" << e << " occured while loading the TTree Formula \"" << m_formula_name << "\" for var=" << RootVariable::GetVariableName() 
	      << " built on TChain=" << RootVariable::GetChain() << " " << RootVariable::GetChain()->GetName() << std::endl;
    assert(false);
  }
  m_formula->SetQuickLoad( true ); /// Experimental...
  // Get leaf type and set it.
  TLeaf* leadFormulaLeaf = m_formula->GetLeaf( 0 );
  TLeaf* leadTreeLeaf = RootVariable::GetChain()->FindLeaf( leadFormulaLeaf->GetName() );
  RootVariable::SetType( leadTreeLeaf->GetTypeName() );
}

double RootVariableFormula::Value( int ndims, va_list args_list) {

  // We use the value cache to avoid unncessary disk access
  if ( useCachedValue==true ) {
    if ( GetChain()->GetReadEntry()==lastCalculatedEntry )
      return cachedValue;	
  }

  // Prep tree
  if ( RootVariable::HasTreeIDChanged() ) {
    m_formula->UpdateFormulaLeaves(); // will this cause a big slow down in performance?
  }
  m_formula->GetNdata(); // asking for the length of the array must set the internal dimension variable. 
  // if i dont call this, the full extent of the array is not loaded!

  // For single value data memebers
  if (RootVariable::IsVariableAnArray()==false)
    return m_formula->EvalInstance(0);

  if ( ndims!=m_var_num_dims ) {
    std::cout << "Calling for the value of " << RootVariable::GetVariableName() << " with the wrong dimensions." << std::endl;
    std::cout << "  ndims passed: " << ndims << std::endl;
    std::cout << "  dimensions according to tree: m_var_num_dims: " << m_var_num_dims << std::endl;
    assert(false);
  }
  
  // For arrays
  for (int n=m_num_specified_dims; n<m_var_num_dims; n++) {
    //std::cout << "instance=" << this << " var=" << m_formula_name << " dim " << n << " of " << m_var_num_dims;
    int indice = va_arg( args_list, int );
    //std::cout << " indice " << indice << std::endl;
    m_input_indices[n-m_num_specified_dims] = indice;
  }
  
  int index = 0;
  for (int i=m_num_specified_dims; i<m_var_num_dims; i++) {
    int offsetsize = 1;
    for (int j=i+1; j<m_var_num_dims; j++ )
      offsetsize *= m_var_dim_nelems[j];
    //std::cout << m_formula_name << ": index " << index << " from arg " << indices.at(i-m_num_specified_dims) << std::endl;
    index += offsetsize*m_input_indices[i-m_num_specified_dims];
  }
  //std::cout << "final index:  " << index << " " << m_formula->EvalInstance(index)<< " " << m_formula->GetNdata() << std::endl;

  if ( useCachedValue ) {
    lastCalculatedEntry = GetChain()->GetReadEntry();
    cachedValue = m_formula->EvalInstance(index);
  }
  return m_formula->EvalInstance(index);
}

std::string RootVariableFormula::PrintValue() {


  // Prep tree
  if ( RootVariable::HasTreeIDChanged() ) {
    m_formula->SetTree( RootVariable::GetChain()->GetTree() );
    m_formula->UpdateFormulaLeaves(); // will this cause a big slow down in performance?
  }
  m_formula->GetNdata(); // asking for the length of the array must set the internal dimension variable. 
  // if i dont call this, the full extent of the array is not loaded!
  // this is likely a bug on ROOT's part.

  
  return m_formula->PrintValue();
}
