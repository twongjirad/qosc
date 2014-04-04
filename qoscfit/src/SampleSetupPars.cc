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

#include "SampleSetupPars.hh"

#include <iostream>
#include <assert.h>
#include <cstdlib>

#include "TChain.h"
#include "SampleBinInfo.hh"

using namespace qosc;

SampleSetupPars::SampleSetupPars( std::string sample_name ) { 
  name = sample_name; 
  ndims = nbinsX = nbinsY = -1;
  binedgesX = binedgesY = NULL;
  formula = selection = weight = "";
  bininfo_name_for_osc = "";
  bininfo_names_for_splines.clear();
}

SampleSetupPars::~SampleSetupPars() {
  delete binedgesX;
  delete binedgesY;
}

void SampleSetupPars::Print() {
  std::cout << "SampleSetupPars[" << name << "]" << std::endl;
  std::cout << " formula=" << formula
	    << " selection="<< selection
	    << " weight=" << weight 
	    << std::endl;
  std::cout << " ndims=" << ndims << " nbinsx=" << nbinsX << " nbinsy=" << nbinsY << std::endl;
  if ( binedgesX ) {
    std::cout << " BinsX: {";
    for (int i=0; i<=nbinsX; i++ ) std::cout << binedgesX[i] << ", ";
    std::cout << "}" << std::endl;
  }
  std::cout << " bin info list: ";
  for (std::vector< SampleBinInfo* >::iterator it=m_bininfo_list.begin(); it!=m_bininfo_list.end(); it++) {
    std::cout << " " << (*it)->name << " (" << (*it) << "), ";
  }
  std::cout << std::endl;
  if ( m_bininfo_list.size()==1 ) {
    std::cout << "  Oscillation Bin Info Definition: " << *m_bininfo_names.begin() << std::endl;
    std::cout << "  Spline/ResponseFunction Bin Info Definition: " << *m_bininfo_names.begin() << std::endl;
  }
  else {
    std::cout << "  Oscillation Bin Info Definition: " << bininfo_name_for_osc << std::endl;
    std::cout << "  Spline/ResponseFunction Bin Info Definition: ";
    for ( std::set< std::string >::iterator it=bininfo_names_for_splines.begin(); it!=bininfo_names_for_splines.end(); it++)
      std::cout << *it << " ";
    std::cout << std::endl;
  }
}
