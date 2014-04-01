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
 * \class RootVariableFormula
 * \ingroup rootvariable
 *
 * \brief This class provides access to ROOT tree data through the TTreeFormula class.
 *
 * Instances of this template class are suppose to be a representation of numbers in a ROOT tree.
 * Beware, this object can be kind of wonky as the kinks have yet to be all ironed out. 
 * Such problems will be addressed in the future.  Stick to single instance of variables
 * and the behavior will most likely be what you expected.
 *   e.g. have an array in the tree named, dir[100][3].  Specify only the instances you want
 *    such as dir[0][0], dir[0][1], etc. However, if you want the whole array, this is reliable.
 *  just specify 'dir' and the whole array should load correctly.
 * TClonedArray objects have worked, but be careful and validate.
 * 
 * \note This class could really use a good cleaning. (Used to be a template class and has cruft that shows that fact.)
 *
 * -----------------------------------------------------------------------------------------------
 */

#ifndef __RootVariableFormula__
#define __RootVariableFormula__

#include <string>
#include <cstdarg>
#include "RootVariable.hh"

class TChain;
class TTreeFormula;
class TTreeFormulaManager; // this is to help figure out dimensions and such (can substitute my dimensional parser)
class TLeaf;
class TBranch;

namespace qosc {

  class RootVariableFormula : public RootVariable {

  public:
  
    RootVariableFormula( std::string variable_name, std::string formula, TChain* source_tree );
    RootVariableFormula( std::string formula, TChain* source_tree );
    virtual ~RootVariableFormula();
    virtual double Value( int ndims, va_list args ); // gets the value of the variable or formula.  number of arguments depends on dimensions of data member
    virtual std::string PrintValue(); // prints the value as a string. note, this is the way to get text items from tree
    std::string GetFormulaName() { return m_formula_name; }; // gets the formula that defines this instance
    TTreeFormula* GetTreeFormula() { return m_formula; }; // get the TTreeFormula Object

  protected:
  
    virtual void LoadVariable( std::string variable_name ); ///< setups the formula

  private:

    void ParseArrayBrackets( std::string array_name, int*& array_info, int& num_array_elems ); ///< parses leaf name to determine dimension of variable -- not perfect
    void AllocateVariable(); ///< allocates the ttree formula instance

  private:

    TLeaf* m_var_leaf; // The TLeaf associated to this data memeber
    TBranch* m_var_branch; // The leaf's branch
    int m_user_array_dims; // number of dimensions of array
    int m_var_num_dims; // number of dimensions of variable array
    int* m_var_dim_nelems; // number of element in each dimension

    std::string m_formula_name;
    TTreeFormula* m_formula;
    TTreeFormulaManager* m_manager;

    // Items to handle the variable length index
    std::vector< std::string > m_flexdims; // variable name which controls flex dimension, order is left to right in array indices
    std::map< std::string, int > m_flexdim_max_elems;
    int m_num_specified_dims;

    // Temp space for indices;
    int* m_input_indices;

    // Caching variables
    int lastCalculatedEntry;
    double cachedValue;
    bool useCachedValue;

  };

}

#endif
