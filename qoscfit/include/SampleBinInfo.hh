#ifndef __SampleBinInfo__
#define __SampleBinInfo__
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
/** --------------------------------------------------------------------------------------------
* \class SampleBinInfo
* \ingroup QoscFit
* \brief Class containing the bin definition for a sample
*
* The class SampleSetupParser is responsible for building this class.
* -------------------------------------------------------------------------------------------*/

#include <string>
#include <vector>
#include <map>

namespace qosc {

  class SampleBinInfo {
  public:
    SampleBinInfo( std::string info_name );
    virtual ~SampleBinInfo();

    bool IsCompletelyDefined() { return true; };
  
    std::string name;
    std::string formula;
    std::string selection;
    std::string weight;
  
    int ndims;

    // Containers used to parse bin description of user
    int nXbinsets;
    std::vector<int> nbinsX;
    std::vector<double*> binedgesX;
    int nbinsX_total;
    double* xedges;
  
    int nYbinsets;
    std::vector<int> nbinsY;
    std::vector<double*> binedgesY;
    int nbinsY_total;
    double* yedges;
  
    // User can define a set of cuts for which a set of bins will be made
    int ncutsets;
    std::vector< std::string > cuts;
    std::vector< int > ncuts;
    std::map< std::string, std::string* > cutlists;
    std::map< std::string, std::string* > cutnames;

    void ProcessBinPattern( std::string command_line, size_t& pos, int& nbinsets, std::vector<int>& nbins, std::vector<double*>& binedges );
    void ProcessBinPatternX( std::string command_line, size_t& pos );
    void ProcessBinPatternY( std::string command_line, size_t& pos );
    void ProcessX();
    void ProcessY();
    void trim( std::string& s );
    void getNumberList( std::string command_line, size_t& pos, int& nnums, double*& numbers );

    void Print();
    
  };

}

#endif
