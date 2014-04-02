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

#include "SampleBinInfo.hh"
#include <iostream>
#include <cmath>
#include <cstdlib>

//#define SAMPLEBININFO_DEBUG

using namespace qosc;

void SampleBinInfo::trim( std::string& s ) {
  size_t start = s.find_first_not_of(" ");
  size_t end = s.find_last_not_of(" \n\r\t");
  s = s.substr(start, end-start+1);
}

void SampleBinInfo::getNumberList( std::string command_line, size_t& pos, int& nnums, double*& numbers ) {
  // Parses pattern that looks like { Y_1, Y_2, Y_3, .... Y_N }
  // nnums = N
  // numbers = array of Y values
  
  bool found_close = false;
  nnums = 0;
  size_t next = command_line.find_first_not_of( " {[", pos )-1;
  std::vector< double > values;
  while ( !found_close ) {
    pos = next+1;
    next = command_line.find_first_of(",}",pos);
    std::string value = command_line.substr( pos, next-pos );
    trim(value);
    double fvalue = -1;
    if ( value!="" ) {
      fvalue = atof( value.c_str() );
      values.push_back( fvalue );
    }
    //std::cout << "getNumberList, value=" << fvalue << ", next up: " << command_line.substr( next, 1 ) << std::endl;
    if ( command_line.substr( next, 1 )=="}" )
      break;
  }
  pos = next+1;
  nnums = values.size();
  numbers = new double[nnums];
  for ( int i=0; i<nnums; i++ )
    numbers[i] = values.at(i);
}

SampleBinInfo::SampleBinInfo( std::string info_name ) { 
  name = info_name; 
  formula = selection = weight = "";
  ncutsets = 0;
  nXbinsets = 0;
  nYbinsets = 0;
  cuts.clear();
}

SampleBinInfo::~SampleBinInfo() {}

void SampleBinInfo::ProcessBinPatternX( std::string command_line, size_t& pos ) {
  ProcessBinPattern( command_line, pos, nXbinsets, nbinsX, binedgesX );
  ProcessX();
}

void SampleBinInfo::ProcessBinPatternY( std::string command_line, size_t& pos ) {
  ProcessBinPattern( command_line, pos, nYbinsets, nbinsY, binedgesY );
  ProcessY();
}

void SampleBinInfo::ProcessBinPattern( std::string command_line, size_t& pos, int& nbinsets, std::vector<int>& nbins, std::vector<double*>& binedges ) {
#ifdef SAMPLEBININFO_DEBUG
  std::cout << __PRETTY_FUNCTION__ << ": " << command_line << " startpos=" << pos << std::endl;
#endif

  size_t next = 0;
  size_t stop = command_line.find_first_of(";:",pos);
  if ( stop==std::string::npos )
    stop = command_line.length();
  while ( pos<stop ) {
    // first get the number of bins in the set
    next = command_line.find_first_of(" ", pos);
    int nbins_interval = atoi(command_line.substr(pos,next-pos).c_str());
    pos = next+1;

    // next get the min and max of the range
    double* binrange;// = new double[2];
    int nnums = 0;
#ifdef SAMPLEBININFO_DEBUG
    std::cout << " parsed number list: " << command_line.substr( pos, std::string::npos ) << std::endl;
#endif
    getNumberList( command_line, pos, nnums, binrange );
#ifdef SAMPLEBININFO_DEBUG
    std::cout << " nnums=" << nnums << std::endl;
#endif

    // store them
    nbins.push_back( nbins_interval );
    binedges.push_back( binrange );
    nbinsets = nbins.size();
#ifdef SAMPLEBININFO_DEBUG
    std::cout << " ProcessBinPattern: binning range #" << nbinsets << ": nbins=" << nbins_interval << " min=" << binrange[0] << " max=" << binrange[nnums-1] 
	      << " (pos=" << pos << ", npos=" << command_line.length() << ", stop=" << stop << ")" << std::endl;
#endif

    pos += 1;
  }
#ifdef SAMPLEBININFO_DEBUG
  std::cout << "End of ProcessBinPattern" << std::endl;
  std::cin.get();
#endif
}

void SampleBinInfo::ProcessX() {
  nbinsX_total = 0;
  for (int i=0; i<nXbinsets; i++) 
    nbinsX_total += nbinsX.at(i);
  xedges = new double[ nbinsX_total+1 ];
  int ibin = 0;
  double last_xmax = 0;
  for (int i=0; i<nXbinsets; i++) {
    double binstep = fabs(binedgesX.at(i)[1]-binedgesX.at(i)[0])/double( nbinsX.at(i) );
    for (int j=0; j<nbinsX.at(i); j++) {
      xedges[ibin] = j*binstep + binedgesX.at(i)[0];
      ibin++;
    }
    last_xmax = binedgesX.at(i)[1];
  }
  xedges[ibin] = last_xmax;
}

void SampleBinInfo::ProcessY() {
  nbinsY_total = 0;
  for (int i=0; i<nYbinsets; i++) 
    nbinsY_total += nbinsY.at(i);
  yedges = new double[ nbinsY_total+1 ];
  int ibin = 0;
  double last_xmax = 0;
  for (int i=0; i<nYbinsets; i++) {
    double binstep = fabs(binedgesY.at(i)[1]-binedgesY.at(i)[0])/double( nbinsY.at(i) );
    for (int j=0; j<nbinsY.at(i); j++) {
      yedges[ibin] = j*binstep + binedgesY.at(i)[0];
      ibin++;
    }
    last_xmax = binedgesY.at(i)[1];
  }
  yedges[ibin] = last_xmax;
}

void SampleBinInfo::Print() {
  std::cout << "SampleBinInfo[" << name << "]" << std::endl;
  std::cout << " formula=" << formula
	    << " selection="<< selection
	    << " weight=" << weight 
	    << std::endl;
  std::cout << " ndims=" << ndims << " nXbinsets=" << nXbinsets << " nYbinsets=" << nYbinsets << std::endl;
  std::cout << " X-bins: ";
  for (int i=0; i<nXbinsets; i++) {
    std::cout << "[(set #" << i+1 << ") nbins=" << nbinsX.at(i) << " {" << binedgesX.at(i)[0] << ", " << binedgesX.at(i)[1] << "}] ";
  }
  std::cout << std::endl;
  for (int i=0; i<ncutsets; i++) {
    std::cout << " -- cutlist #" << i+1 << " of " << ncutsets << ": " << cuts.at(i) << ", ncuts=" << ncuts.at(i) << " --" << std::endl;
    std::string* cutnamelist = cutnames[cuts.at(i)];
    if ( cutnamelist!=NULL) {
      for (int j=0; j<ncuts.at(i); j++)  {
	std::cout << "   (" << j+1 << ") " << cutnames[cuts.at(i)][j]  << " = '" << cutlists[cuts.at(i)][j] << "'" << std::endl;
      }
    }
    else {
      std::cout << "Cut list not given names!!" << std::endl;
      for (int j=0; j<ncuts.at(i); j++)  {
	std::cout << "   (" << j+1 << ") UNNAMED  = '" << cutlists[cuts.at(i)][j] << "'" << std::endl;
      }
    }
  }
}

