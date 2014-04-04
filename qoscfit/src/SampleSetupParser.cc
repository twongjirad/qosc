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


#include "SampleSetupParser.hh"
#include <cstdlib>
#include <iostream>
#include <assert.h>
#include <sstream>
#include <stdio.h>

#include "TChain.h"

#include "SampleSetupPars.hh"
#include "SampleBinInfo.hh"

#include "TruthBinInfo.hh"

#include "SampleManagerROOTPDF.hh"
#include "SampleOscfitROOTPDF.hh"

using namespace qosc;

std::string SampleSetupParser::__commands__[kNumCommands] 
= { "#", /// (0) Comment line
    "DEF_SAMPLE", /// (1) defines a type of observable
    "DEF_SAMPLE_BINS", // (2) defines binning for the sample's observable quantity
    "DEF_SAMPLE_VAR", // (3) defines the ROOT tree variable we are using to create the observable
    "DEF_SAMPLE_CUT",  // (4) defines the selection criteria for a sample (optional)
    "DEF_SAMPLE_WGT", // (5) defines the ROOT formula used select events for a sample (optional)
    "DEF_SAMPLE_OSC_BININFO", // (6) defines the INFO_BINS object that will serve for oscillation calculations. [ if not specified, and only 1 bin info assigned, will assume its for oscillations ]
    "DEF_SAMPLE_SPLINE_BININFO", // (7) defines the INFO_BINS object that will serve the spline binnings. [ if not specified, and only 1 bin info assigned, will assume its for oscillations ]
    "DEF_INFO_BINS", // (8) defines a new type of bin information to associate with each reconstructed bin
    "DEF_INFO_VAR", // (9) defines the ROOT formula to fill into the bin info histograms
    "DEF_INFO_CUTS", // (10) defines a set of mutually exclusive cuts. each cut has its own set of bins to fill
    "DEF_INFO_CUTS2", // (11) defines a set of mutually exclusive cuts. each cut has its own set of bins to fill
    "ADD_SAMPLE_INFO", // (12) assign the info instances to the sample definitions
    "DEF_INFO_CUTNAMES", // (13) defines a set of mutually exclusive cuts. each cut has its own set of bins to fill
    "DEF_CHAIN", // (14) Define the ROOT chains used to build the sample histograms
    "ADD_SAMPLE_CHAIN", // (15) Assign chain to a sample
    "OPT_SPLINE_MERGEBINS", // (16) Definings binning over which splines are merged
    "SET_PARSER_VERBOSITY", // (17) Set verbosity. Executes at position file was read
};
						   
std::set< std::string > SampleSetupParser::__command_set__( SampleSetupParser::__commands__, SampleSetupParser::__commands__+SampleSetupParser::kNumCommands );

// ------------------------------------------------------------------------------------------------------
// STRING PARSING TOOLS
void trim( std::string& s ) {
  size_t start = s.find_first_not_of(" ");
  size_t end = s.find_last_not_of(" \n\r\t");
  if ( start==std::string::npos ) return;
  s = s.substr(start, end-start+1);
}

void SampleSetupParser::getNumberList( std::string command_line, size_t& pos, int& nnums, double*& numbers ) {
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
    if ( fVerbosity>=2 )
      std::cout << "getNumberList, value=" << fvalue << ", next up: " << command_line.substr( next, 1 ) << std::endl;
    if ( command_line.substr( next, 1 )=="}" )
      break;
  }
  pos = next+1;
  nnums = values.size();
  numbers = new double[nnums];
  for ( int i=0; i<nnums; i++ )
    numbers[i] = values.at(i);
}

void SampleSetupParser::getStringList( std::string command_line, size_t& pos, int nnums, std::string *texts ) {
  size_t next = command_line.find_first_not_of( " {[", pos )-1;
  for ( int n=0; n<nnums; n++ ) {
    pos = next+1;
    next = command_line.find_first_of(",}",pos);
    std::string str = command_line.substr( pos, next-pos );
    trim(str);
    texts[n] = str;
  }
  pos = next+1;
}

void SampleSetupParser::parseStringList( std::string command_line, size_t& pos, int nnums, std::vector< std::string >& list ) {
  size_t next = command_line.find_first_not_of( " {[", pos )-1;
  for ( int n=0; n<nnums; n++ ) {
    pos = next+1;
    next = command_line.find_first_of(" ,}",pos);
    std::string str = command_line.substr( pos, next-pos );
    trim(str);
    list.push_back( str );
  }
  pos = next+1;
}

// ------------------------------------------------------------------------------------------------------


SampleSetupParser::SampleSetupParser( std::string setupfilename ) {
  m_setupfilename = setupfilename;
  m_input_file = NULL;
  SetVerbosity( 0 );
  OpenFile( m_setupfilename );
}


SampleSetupParser::~SampleSetupParser() {
  CloseFile();
}

void SampleSetupParser::OpenFile( std::string filename ) {
  if ( !m_input_file ) {
    m_input_file = new std::ifstream( filename.c_str() );
    if ( fVerbosity>=0 ) std::cout << "SampleSetupParser: Opened " << filename << std::endl;
  }
  else {
    std::cout << "SampleSetupParser: File already opened!" << std::endl;
    assert(false);
  }

  // Here should check file state
  if ( !m_input_file->good() ) {
    std::cout << "SampleSetupParser::OpenFile. Error opening " << filename << std::endl;
     assert(false);
     exit(EXIT_FAILURE);
  }
}

void SampleSetupParser::CloseFile() {
  if ( m_input_file ) {
    m_input_file->close();
    delete m_input_file;
  }
  m_input_file = NULL;
}

std::string SampleSetupParser::GetNextLine() {

  std::string sbuffer;

  char buffer[5000];

  (*m_input_file).getline(buffer,5000);
  sbuffer = buffer;
  while ( sbuffer=="" && !m_input_file->eof() && m_input_file->good() ) {
    (*m_input_file).getline(buffer,5000);
    sbuffer = buffer;
    if ( fVerbosity>=2) {
      std::cout << "buffer: "<< buffer << std::endl;
      std::cin.get();
    }
  }
  if ( fVerbosity>=2 ) std::cout << "GetNextLine: sbuffer: "<< sbuffer << std::endl;
  return sbuffer;
}

int SampleSetupParser::GetNextCommand( std::string& command, std::string& arg ) {

  std::string sbuffer = GetNextLine();
  trim( sbuffer );


  // CHECK FOR BLANK LINE
  if ( sbuffer=="" ) {
    command = "COMMENT_LINE";
    arg = "";
    return 0;
  }


  // CHECK FOR COMMENT SYMBOL
  std::string firstchar = sbuffer.substr(0,1);
  if ( firstchar=="#" ) {
    command = "COMMENT_LINE";
    arg = "";
    return 0;
  }

  // Find end of command
  size_t pos = sbuffer.find_first_of(" ");
  if ( pos==std::string::npos ) {
    // Got to the end of the string?
    command = sbuffer;
    return -1;
  }

  // extract command
  command = sbuffer.substr(0,pos);
  if ( SampleSetupParser::__command_set__.find( command )==SampleSetupParser::__command_set__.end() ) {
    // Not a valid command
    return -2;
  }

  // extract argument list
  arg = sbuffer.substr(pos+1,std::string::npos);
  if (fVerbosity>=2) 
    std::cout << command << " - " << arg << std::endl;
  trim(arg);

  return 0;
}

void SampleSetupParser::ParseCommand_DEF_SAMPLE( std::string command_line ) {
  trim(command_line);
  SampleSetupPars* sample = new SampleSetupPars( command_line );
  m_sample_setups[sample->name] = sample;
  if ( fVerbosity>=0 ) std::cout << "Sample Defined: " << sample->name << std::endl;
}

void SampleSetupParser::ParseCommand_DEF_SAMPLE_VAR( std::string command_line ) {
  size_t pos = 0;
  size_t next = command_line.find_first_of(" ");

  // Get the sample name and use it to get the SampleSetupPars instance stored in m_sample_setups  
  std::string samplename = command_line.substr(0,next);
  SampleSetupPars* sample = NULL;
  if ( m_sample_setups.find( samplename )==m_sample_setups.end() ) {
    std::cout << "SampleSetupParser::ParseCommand_DEF_SAMPLE_VAR. Could not find sample=" << samplename << std::endl;
    std::cout << " command: " << command_line << std::endl; 
    exit( EXIT_FAILURE );
  }
  else {
    sample = m_sample_setups[samplename];
  }
  
  // Get the variable formula
  pos = next+1;
  next = command_line.find_first_of(";:",pos);
  std::string sampleformula;
  if ( next==std::string::npos ) {
    // 1D
    sampleformula = command_line.substr( pos );
    trim(sampleformula);
    sample->formula = sampleformula;
  }
  else {
    // 2D
    sampleformula = command_line.substr( pos, next-pos );
    trim(sampleformula);
    sample->formula = sampleformula;

    sampleformula = command_line.substr( next+1 );
    trim(sampleformula);
    sample->formulaY = sampleformula;
  }

  if ( fVerbosity>=1 ) std::cout << "Parsed formula for sample='" << samplename << "': " << sampleformula << std::endl;
}

void SampleSetupParser::ParseCommand_DEF_SAMPLE_CUT( std::string command_line ) {
  size_t pos = 0;
  size_t next = command_line.find_first_of(" ");

  // Get the sample name and use it to get the SampleSetupPars instance stored in m_sample_setups  
  std::string samplename = command_line.substr(0,next);
  SampleSetupPars* sample = NULL;
  if ( m_sample_setups.find( samplename )==m_sample_setups.end() ) {
    std::cout << "SampleSetupParser::ParseCommand_DEF_SAMPLE_CUT. Could not find sample=" << samplename << std::endl;
    std::cout << " command: " << command_line << std::endl; 
    exit( EXIT_FAILURE );
  }
  else {
    sample = m_sample_setups[samplename];
  }
  
  // Get the variable formula
  pos = next+1;
  std::string sampleformula = command_line.substr( pos );
  trim(sampleformula);
  sample->selection = sampleformula;
  if ( fVerbosity>=1 ) std::cout << "Loaded cut formula for sample=" << sample->name << ": " << sample->selection << std::endl;
}

void SampleSetupParser::ParseCommand_DEF_SAMPLE_WGT( std::string command_line ) {
  size_t pos = 0;
  size_t next = command_line.find_first_of(" ");

  // Get the sample name and use it to get the SampleSetupPars instance stored in m_sample_setups  
  std::string samplename = command_line.substr(0,next);
  SampleSetupPars* sample = NULL;
  if ( m_sample_setups.find( samplename )==m_sample_setups.end() ) {
    std::cout << "SampleSetupParser::ParseCommand_DEF_SAMPLE_WGT. Could not find sample=" << samplename << std::endl;
    std::cout << " command: " << command_line << std::endl; 
    exit( EXIT_FAILURE );
  }
  else {
    sample = m_sample_setups[samplename];
  }
  
  // Get the variable formula
  pos = next+1;
  std::string sampleformula = command_line.substr( pos );
  trim(sampleformula);
  sample->weight = sampleformula;
  if ( fVerbosity>=1 ) std::cout << "Loaded eventweight formula for sample=" << sample->name << ": " << sample->weight << std::endl;
}

void SampleSetupParser::ParseCommand_DEF_SAMPLE_BINS( std::string command_line ) {
  size_t pos = 0;
  size_t next = command_line.find_first_of(" ");

  // Get the sample name and use it to get the SampleSetupPars instance stored in m_sample_setups
  std::string samplename = command_line.substr(0,next);
  SampleSetupPars* sample = NULL;
  if ( m_sample_setups.find( samplename )==m_sample_setups.end() ) {
    std::cout << "SampleSetupParser::ParseCommand_DEF_SAMPLE_BINS. Could not find sample=" << samplename << std::endl;
    std::cout << " command: " << command_line << std::endl; 
    exit( EXIT_FAILURE );
  }
  else {
    sample = m_sample_setups[samplename];
  }

  if (fVerbosity>=1 ) 
    std::cout << "DefineSampleBins: for sample=" << samplename << " (" << sample << ")" << std::endl;

  // Get the number of dimensions
  pos = next+1;
  next = command_line.find_first_of(" ",pos);
  sample->ndims = atoi( command_line.substr( pos, next-pos ).c_str() );

  if ( fVerbosity>=1 ) 
    std::cout << "DefineSampleBins: sample ndims=" << sample->ndims << std::endl;

  // ----- Get X binning --------
  pos = next+1;
  next = command_line.find_first_of(":;",pos);
  std::string binpattern;
  if ( next!=std::string::npos )
    binpattern = command_line.substr(pos,next-pos);
  else
    binpattern = command_line.substr(pos);
  trim(binpattern);
  SampleBinInfo* bininfo = new SampleBinInfo( samplename+"_xbins" );
  sample->bininfo = bininfo;
  sample->nbinsX = 0;
  size_t binstart = 0;
  bininfo->ProcessBinPatternX( binpattern, binstart );
  sample->nbinsX = bininfo->nbinsX_total;
  sample->binedgesX = new double[ bininfo->nbinsX_total+1 ];
  for (int i=0; i<=bininfo->nbinsX_total; i++)
    sample->binedgesX[i] = bininfo->xedges[i];
  
  if ( sample->ndims>1 ) {
    // ------- Get Y binning ------
    pos = next+1;
    binpattern = command_line.substr(pos);
    trim(binpattern);
    binstart = 0;
    sample->bininfo->ProcessBinPatternY( binpattern, binstart );
    sample->nbinsY = sample->bininfo->nbinsY_total;
    sample->binedgesY = new double[ sample->bininfo->nbinsY_total+1 ];
    for (int i=0; i<=sample->bininfo->nbinsY_total; i++)
      sample->binedgesY[i] = bininfo->yedges[i];
  }
  assert(sample->ndims<=2);

}

void SampleSetupParser::ParseCommand_DEF_INFO_BINS( std::string command_line ) {
  trim(command_line);
  size_t pos = 0;
  size_t next = pos;

  // get name of info
  next = command_line.find_first_of( " " );
  std::string info_name = command_line.substr(pos, next-pos);
  pos = next + 1;
  
  // define a new info par instance
  SampleBinInfo* bininfo = new SampleBinInfo( info_name );
  m_bininfo_dict[info_name] = bininfo;
  if ( fVerbosity>=1 ) std::cout << "Creating new SampleBinInfo instance: " << info_name << " (" << bininfo << ")" << std::endl;
  
  // get the ndims
  next = command_line.find_first_of(" ", pos);
  bininfo->ndims = atoi( command_line.substr(pos,next-pos).c_str() );
  pos = next+1;
  if ( fVerbosity>=2 ) std::cout << " ndims: " << bininfo->ndims << std::endl;

  // now get the X binning. specified as groupings of bins
  while ( pos<command_line.length() ) {
    // first get the number of bins in the set
    next = command_line.find_first_of(" ", pos);
    int nbins = atoi(command_line.substr(pos,next-pos).c_str());
    pos = next+1;

    // next get the min and max of the range
    double* binrange;// = new double[2];
    int nnums = 0;
    getNumberList( command_line, pos, nnums, binrange );

    // store them
    bininfo->nbinsX.push_back( nbins );
    bininfo->binedgesX.push_back( binrange );
    bininfo->nXbinsets = bininfo->nbinsX.size();
    if ( fVerbosity>=1 )
      std::cout << " binning range #" << bininfo->nXbinsets << ": nbins=" << nbins << " min=" << binrange[0] << " max=" << binrange[1] 
		<< " (pos=" << pos << ", npos=" << command_line.length() << ")" << std::endl;
    pos += 1;
  }
  bininfo->ProcessX();
  if ( fVerbosity>=1 ) std::cout << "Total number of bins for X: " << bininfo->nbinsX_total << std::endl;
  
  // haven't arrange for 2D just yet.
  assert( bininfo->ndims==1 );
  
}

void SampleSetupParser::ParseCommand_OPT_SPLINE_MERGEBINS( std::string command_line ) {
  trim(command_line);
  size_t pos = 0;
  size_t next = pos;

  // get name of info
  next = command_line.find_first_of( " " );
  std::string info_name = command_line.substr(pos, next-pos);
  pos = next + 1;
  
  // define a new info par instance
  SampleBinInfo* bininfo = new SampleBinInfo( info_name );
  m_splinebins_dict[info_name] = bininfo;
  if ( fVerbosity>=1 ) std::cout << "Creating new SampleBinInfo instance: " << info_name << " (" << bininfo << ")" << std::endl;
  
  // get the ndims
  next = command_line.find_first_of(" ", pos);
  bininfo->ndims = atoi( command_line.substr(pos,next-pos).c_str() );
  pos = next+1;
  if ( fVerbosity>=2 ) std::cout << " ndims: " << bininfo->ndims << std::endl;

  // now get the X binning. specified as groupings of bins
  while ( pos<command_line.length() ) {
    // first get the number of bins in the set
    next = command_line.find_first_of(" ", pos);
    int nbins = atoi(command_line.substr(pos,next-pos).c_str());
    pos = next+1;

    // next get the min and max of the range
    double* binrange;// = new double[2];
    int nnums = 0;
    getNumberList( command_line, pos, nnums, binrange );

    // store them
    bininfo->nbinsX.push_back( nbins );
    bininfo->binedgesX.push_back( binrange );
    bininfo->nXbinsets = bininfo->nbinsX.size();
    if ( fVerbosity>=1 )
      std::cout << " binning range #" << bininfo->nXbinsets << ": nbins=" << nbins << " min=" << binrange[0] << " max=" << binrange[1] 
		<< " (pos=" << pos << ", npos=" << command_line.length() << ")" << std::endl;
    pos += 1;
  }
  bininfo->ProcessX();
  if ( fVerbosity>=1 ) 
    std::cout << "Total number of bins for X: " << bininfo->nbinsX_total << std::endl;
  
  // haven't arrange for 2D just yet.
  assert( bininfo->ndims==1 );
  
}

void SampleSetupParser::ParseCommand_DEF_INFO_VAR( std::string command_line ) {
  size_t pos = 0;
  size_t next = command_line.find_first_of(" ");

  // Get the bin info ROOT variable formula
  std::string infoname = command_line.substr(0,next);
  SampleBinInfo* bininfo = NULL;
  if ( m_bininfo_dict.find( infoname )==m_bininfo_dict.end() ) {
    std::cout << "SampleSetupParser::ParseCommand_DEF_INFO_VAR. Could not find bininfo=" << infoname << std::endl;
    std::cout << " command: " << command_line << std::endl; 
    exit( EXIT_FAILURE );
  }
  else {
    bininfo = m_bininfo_dict[infoname];
  }
  
  // Get the variable formula
  pos = next+1;
  std::string sampleformula = command_line.substr( pos );
  trim(sampleformula);
  bininfo->formula = sampleformula;
  if ( fVerbosity>=2 ) std::cout << "BinInfo formula: " << sampleformula << std::endl;

}

void SampleSetupParser::ParseCommand_DEF_INFO_CUTS( std::string command_line ) {
  size_t pos = 0;
  size_t next = command_line.find_first_of(" ");

  // get name and instance of bininfo
  std::string infoname = command_line.substr(pos,next-pos);
  trim(infoname);
  SampleBinInfo* bininfo = m_bininfo_dict[ infoname ];
  pos = next+1;
  
  // get name for this set of cuts
  next = command_line.find_first_of(" ",pos);
  std::string cutname = command_line.substr(pos,next-pos);
  trim(cutname);
  pos = next+1;

  // get number of cuts
  next = command_line.find_first_of(" ",pos);
  int ncuts = atoi( command_line.substr(pos, next-pos ).c_str() );
  pos = next+1;

  // get list of cuts
  if ( fVerbosity>=1 ) std::cout << "Stored cut list for info=" << infoname << ": ncuts=" << ncuts << " cutname=" << cutname << std::endl;
  std::string* cutlist = new std::string[ncuts];
  getStringList( command_line, pos, ncuts, cutlist );

  // store them
  bininfo->cuts.push_back( cutname );
  bininfo->ncuts.push_back( ncuts );
  bininfo->cutlists[cutname] = cutlist;
  if ( fVerbosity>=2 )  {
    for (int i=0; i<ncuts; i++) {
      std::cout << " '" << bininfo->cutlists[ cutname ][i] << "', ";
    }
    std::cout << " " << bininfo->cuts.size() << std::endl;
  }
  bininfo->ncutsets = bininfo->cuts.size();

}

void SampleSetupParser::ParseCommand_DEF_INFO_CUTS2( std::string command_line ) {
  // Example: DEF_INFO_CUTS2 info_bin_name info_cut_set_name 2 {numu:iflux==2 && ixsec==2,numubar:iflux==-2}
  size_t pos = 0;
  size_t next = command_line.find_first_of(" ");

  // get name and instance of bininfo
  std::string infoname = command_line.substr(pos,next-pos);
  trim(infoname);
  SampleBinInfo* bininfo = m_bininfo_dict[ infoname ];
  pos = next+1;
  
  // get name for this set of cuts
  next = command_line.find_first_of(" ",pos);
  std::string cutname = command_line.substr(pos,next-pos);
  trim(cutname);
  pos = next+1;

  // get number of cuts
  next = command_line.find_first_of(" ",pos);
  int ncuts = atoi( command_line.substr(pos, next-pos ).c_str() );
  pos = next+1;

  // get list of cuts
  if ( fVerbosity>=1 ) std::cout << "Storing cut list for info=" << infoname << ": ncuts=" << ncuts << " cutname=" << cutname << std::endl;
  std::string* cutlist = new std::string[ncuts];
  getStringList( command_line, pos, ncuts, cutlist );

  // store them
  bininfo->cuts.push_back( cutname );
  bininfo->ncuts.push_back( ncuts );
  bininfo->cutlists[cutname] = new std::string[ncuts];
  bininfo->cutnames[cutname] = new std::string[ncuts];

  // get cuts out of list
  for (int i=0; i<ncuts; i++) {
    std::string cutdef = cutlist[i];
    if ( cutdef.find(":")==std::string::npos ) {
      std::cout << "ParseCommand_DEF_INFO_CUTS2: Illegal cut definition" << std::endl;
      std::cout << "must look like : { ... , cutname:cut formula, .... }" << std::endl;
      exit(EXIT_FAILURE);
    }
    std::string cutnom = cutdef.substr(0,cutdef.find(":"));
    std::string cutform = cutdef.substr(cutdef.find(":")+1,std::string::npos);
    bininfo->cutlists[cutname][i] = cutform;
    bininfo->cutnames[cutname][i] = cutnom;
  }
  delete [] cutlist;

  if ( fVerbosity>=2 )  {
    for (int i=0; i<ncuts; i++) {
      std::cout << " '" << bininfo->cutnames[ cutname ][i] << ":" << bininfo->cutlists[ cutname ][i] << "', ";
    }
    std::cout << " " << bininfo->cuts.size() << std::endl;
  }
  bininfo->ncutsets = bininfo->cuts.size();
  
}

void SampleSetupParser::ParseCommand_DEF_INFO_CUTNAMES( std::string command_line ) {
  size_t pos = 0;
  size_t next = command_line.find_first_of(" ");

  // get name and instance of bininfo
  std::string infoname = command_line.substr(pos,next-pos);
  trim(infoname);
  SampleBinInfo* bininfo = m_bininfo_dict[ infoname ];
  pos = next+1;
  
  // get name for this set of cuts
  next = command_line.find_first_of(" ",pos);
  std::string cutname = command_line.substr(pos,next-pos);
  trim(cutname);
  pos = next+1;

  // get number of cuts
  next = command_line.find_first_of(" ",pos);
  int ncuts = atoi( command_line.substr(pos, next-pos ).c_str() );
  pos = next+1;

  // get list of cut names
  if ( fVerbosity>=2 ) std::cout << "Stored cut list for info=" << infoname << ": ncuts=" << ncuts << " cutname=" << cutname << std::endl;
  std::string* cutlist = new std::string[ncuts];
  getStringList( command_line, pos, ncuts, cutlist );

  // store them
  bininfo->cutnames[cutname] =  cutlist;
  if ( fVerbosity>=2 ) {
    for (int i=0; i<ncuts; i++) {
      std::cout << " '" << bininfo->cutnames[cutname][i] << "', ";
    }
    std::cout << std::endl;
  }
  
}

void SampleSetupParser::ParseCommand_ADD_SAMPLE_INFO( std::string command_line ) {
  size_t pos = 0;
  size_t next = command_line.find_first_of(" ");
  
  // get sample name and instance
  std::string samplename = command_line.substr(0,next);
  trim(samplename);
  SampleSetupPars* sample = NULL;
  if ( m_sample_setups.find( samplename )==m_sample_setups.end() ) {
    std::cout << "SampleSetupParser::ParseCommand_ADD_SAMPLE_INFO. Could not find sample=" << samplename << std::endl;
    std::cout << " command: " << command_line << std::endl; 
    std::cout << " expecting: ADD_SAMPLE_INFO [sample name] [info name list]" << std::endl;
    exit( EXIT_FAILURE );
  }
  else {
    sample = m_sample_setups[samplename];
  }
  pos = next+1;
  
  // Get the list of bin info types to assign to each sample bin
  next = command_line.find_first_of( " ", pos );
  int ninfos = atoi( command_line.substr(pos,next-pos).c_str() );
  pos = next+1;

  std::string* infolist = new std::string[ninfos];
  getStringList( command_line, pos, ninfos, infolist );
  if ( fVerbosity>=0 ) std::cout << "Adding bin info classes to " << samplename << ": ";
  for (int i=0; i<ninfos; i++) {
    sample->m_bininfo_list.push_back( m_bininfo_dict[ infolist[i] ] );
    sample->m_bininfo_names.insert( sample->name+"_"+infolist[i] );
    if ( fVerbosity>=0 ) std::cout << infolist[i] << ", ";
  }
  if ( fVerbosity>=0 ) std::cout << std::endl;
  delete [] infolist;
}

void SampleSetupParser::ParseCommand_DEF_CHAIN( std::string command ) {
  size_t pos = 0;
  size_t next = 0;

  // get name of chain
  next = command.find_first_of(" ");
  std::string chainname = command.substr(pos,next-pos);
  pos = next+1;

  // get tree name
  next = command.find_first_of(" ", pos);
  std::string treename = command.substr(pos,next-pos);
  pos = next+1;
  
  // make TChain instance, add to dictionary
  TChain* chain = NULL;
  if ( m_chain_dict.find( chainname )==m_chain_dict.end() ) {
    chain = new TChain( treename.c_str() );
    m_chain_dict[chainname] = chain;
  }
  else {
    chain = m_chain_dict[chainname];
  }

  // get rootfile location
  next = command.find_first_of(" ",pos);
  std::string chain_file = command.substr(pos, next-pos);
  trim(chain_file);
  chain->Add( chain_file.c_str() );
  if ( fVerbosity>=0 )
    std::cout << "Loaded chain, " << treename << " (" << chain << "), from file: " << chain_file << ". entries " << chain->GetEntries() << std::endl;
}

void SampleSetupParser::ParseCommand_ADD_SAMPLE_CHAIN( std::string command_line ) {
  size_t pos = 0;
  size_t next = command_line.find_first_of(" ");

  // get sample name and instance
  std::string samplename = command_line.substr(0,next);
  trim(samplename);
  SampleSetupPars* sample = NULL;
  if ( m_sample_setups.find( samplename )==m_sample_setups.end() ) {
    std::cout << "SampleSetupParser::ParseCommand_ADD_SAMPLE_CHAIN. Could not find sample=" << samplename << std::endl;
    std::cout << " command: '" << command_line << "'" << std::endl; 
    exit( EXIT_FAILURE );
  }
  else {
    sample = m_sample_setups[samplename];
  }
  pos = next+1;
  
  // Get the name of the chain to add
  next = command_line.find_first_of( " ", pos );
  std::string chainname = command_line.substr(pos,next-pos);
  pos = next+1;

  // Add chain
  sample->chain = m_chain_dict[chainname];
  sample->chain_name = chainname;

}


void SampleSetupParser::ParseCommand_DEF_SAMPLE_OSC_BININFO( std::string command_line ) {
  trim(command_line);

  size_t pos = 0;
  size_t next = command_line.find_first_of(" ");

  // get sample name and instance
  std::string samplename = command_line.substr(0,next);
  trim(samplename);
  pos = next+1;

  // get the bin info name
  next = command_line.find_first_of(" ",pos);
  std::string bininfoname = command_line.substr(pos,next);
  trim(bininfoname);
  
  SampleSetupPars* sample_pars = NULL;
  if ( m_sample_setups.find( samplename )==m_sample_setups.end() ) {
    std::cout << "SampleSetupParser::ParseCommand_DEF_SAMPLE_OSC_BININFO. Could not find sample=" << samplename << std::endl;
    std::cout << " command: " << command_line << std::endl; 
    std::cout << " expect: DEF_SAMPLE_OSC_BININFO [sample name] [bin info name]" << std::endl;
    exit( EXIT_FAILURE );
  }
  else {
    sample_pars = m_sample_setups[samplename];
  }
  sample_pars->bininfo_name_for_osc = bininfoname;
  std::cout << "ParseCommand_DEF_SAMPLE_OSC_BININFO: " << sample_pars->bininfo_name_for_osc << std::endl;
}

void SampleSetupParser::ParseCommand_DEF_SAMPLE_SPLINE_BININFO( std::string command_line ) {
  trim(command_line);

  size_t pos = 0;
  size_t next = command_line.find_first_of(" ");

  // get sample name and instance
  std::string samplename = command_line.substr(0,next);
  trim(samplename);
  pos = next+1;

  // get the bin info names
  next = command_line.find_first_of(" ",pos);
  std::string bininfoname = command_line.substr(pos,next);
  trim(bininfoname);
  
  SampleSetupPars* sample_pars = NULL;
  if ( m_sample_setups.find( samplename )==m_sample_setups.end() ) {
    std::cout << "SampleSetupParser::ParseCommand_DEF_SAMPLE_SPLINE_BININFO. Could not find sample=" << samplename << std::endl;
    std::cout << " command: " << command_line << std::endl; 
    std::cout << " expecting: DEF_SAMPLE_SPLINE_BININFO [sample name] [bin info name]" << std::endl;
    exit( EXIT_FAILURE );
  }
  else {
    sample_pars = m_sample_setups[samplename];
  }
  std::cout << bininfoname << std::endl;
  std::cout << sample_pars << std::endl;
  sample_pars->bininfo_names_for_splines.insert(bininfoname);
}

void SampleSetupParser::ParseCommand_SET_PARSER_VERBOSITY( std::string command_line ) {
  trim(command_line);
  SetVerbosity( std::atoi( command_line.c_str() ) );
}

void SampleSetupParser::ParseSetup() {
  int result = 0;
  while ( !m_input_file->eof() && result==0 ) {
    std::string command, arg;
    result = GetNextCommand( command, arg );
    if ( fVerbosity>=1 ) 
      std::cout << result << ": command="<< command << " arg=" << arg << std::endl;

    if ( command=="COMMENT_LINE" )
      continue;

    if ( command=="DEF_SAMPLE" )
      ParseCommand_DEF_SAMPLE( arg );
    else if ( command=="DEF_SAMPLE_VAR" )
      ParseCommand_DEF_SAMPLE_VAR( arg );
    else if ( command=="DEF_SAMPLE_CUT" )
      ParseCommand_DEF_SAMPLE_CUT( arg );
    else if ( command=="DEF_SAMPLE_WGT" )
      ParseCommand_DEF_SAMPLE_WGT( arg );
    else if ( command=="DEF_SAMPLE_BINS" )
      ParseCommand_DEF_SAMPLE_BINS( arg );
    else if ( command=="DEF_INFO_BINS" )
      ParseCommand_DEF_INFO_BINS( arg );
    else if ( command=="DEF_INFO_VAR" )
      ParseCommand_DEF_INFO_VAR( arg );
    else if ( command=="DEF_INFO_CUTS" )
      ParseCommand_DEF_INFO_CUTS( arg );
    else if ( command=="DEF_INFO_CUTS2" )
      ParseCommand_DEF_INFO_CUTS2( arg );
    else if ( command=="DEF_INFO_CUTNAMES" )
      ParseCommand_DEF_INFO_CUTNAMES( arg );
    else if ( command=="ADD_SAMPLE_INFO")
      ParseCommand_ADD_SAMPLE_INFO( arg );
    else if ( command=="DEF_CHAIN" )
      ParseCommand_DEF_CHAIN( arg );
    else if ( command=="ADD_SAMPLE_CHAIN" )
      ParseCommand_ADD_SAMPLE_CHAIN( arg );
    else if ( command=="OPT_SPLINE_MERGEBINS" )
      ParseCommand_OPT_SPLINE_MERGEBINS( arg );
    else if ( command=="DEF_SAMPLE_OSC_BININFO" )
      ParseCommand_DEF_SAMPLE_OSC_BININFO( arg );
    else if ( command=="DEF_SAMPLE_SPLINE_BININFO" )
      ParseCommand_DEF_SAMPLE_SPLINE_BININFO( arg );
    else if ( command=="SET_PARSER_VERBOSITY" )
      ParseCommand_SET_PARSER_VERBOSITY( arg );
    else {
      std::cout << "Did not recognize command: " << command << " with arg=" << arg << std::endl;
      std::cout << "Possible commands include: " << std::endl;
      for ( std::set< std::string >::iterator it=__command_set__.begin(); it!=__command_set__.end(); it++ ) {
	std::cout << " " << *it << std::endl;
      }
      assert(false);
    }
  }

  if ( fVerbosity>=1 ) 
    PrintSetup();
}

void SampleSetupParser::PrintSetup() {
  
  std::cout << "====================================================" << std::endl;
  std::cout << "SampleSetupParser" << std::endl;
  std::cout << "----------------------------------------------------" << std::endl;
  for ( std::map< std::string, SampleSetupPars* >::iterator it=m_sample_setups.begin(); it!=m_sample_setups.end(); it++) {
    (*it).second->Print();
  }

  std::cout << "====================================================" << std::endl;
  
  for ( std::map< std::string, SampleBinInfo* >::iterator it=m_bininfo_dict.begin(); it!=m_bininfo_dict.end(); it++) {
    (*it).second->Print();
  }

  std::cout << "====================================================" << std::endl;
}


///-------------------------------------------------------------------------------------------------------------------

void SampleSetupParser::SetupSamples( SampleManagerROOTPDF* sampleMan ) {

  if ( fVerbosity>=0 )
    std::cout << "SampleSetupParser::SetupSamples for SampleMan=" << sampleMan << std::endl;
  
  // This method will take the loaded samples and true energy bins and load the given sample Manager with 
  // the necessary histograms

  // first load chains
  for ( std::map< std::string, TChain* >::iterator it_chain=m_chain_dict.begin(); it_chain!=m_chain_dict.end(); it_chain++ ) {
    sampleMan->AddRootChain( (*it_chain).first, (*it_chain).second );
  }

  // next create samples
  for ( std::map< std::string, SampleSetupPars* >::iterator it_sample = m_sample_setups.begin(); it_sample!=m_sample_setups.end(); it_sample++ ) {
    
    SampleSetupPars* setup = (*it_sample).second;
    if ( fVerbosity>=0 )
      std::cout << "Creating Sample: " << setup->name << " based on chain=" << setup->chain << std::endl;
    m_samples_defined.insert( setup->name );

    // First define the Root Variable PDF object
    HistRootVariable* hist;
    if ( setup->ndims==1 ) {
      // 1D Histogram
      if ( setup->selection!="" && setup->weight!="" ) {
	hist = new HistRootVariable( setup->name, setup->formula, setup->weight, setup->selection, setup->chain );
      }
      else if ( setup->selection!="" ) {
	hist = new HistRootVariable( setup->name, setup->formula, NULL, setup->selection, setup->chain );
      }
      else if ( setup->weight!="" ) {
	hist = new HistRootVariable( setup->name, setup->formula, setup->weight, setup->chain );
      }
						 
      // next define the histogram
      hist->SetHistogram( setup->name, setup->nbinsX, setup->binedgesX );
    }
    else if ( setup->ndims==2 ) {
      // 2D Histogram
      if ( setup->selection!="" && setup->weight!="" ) {
	hist = new HistRootVariable( setup->name, setup->formula, setup->formulaY, setup->weight, setup->selection, setup->chain );
      }
      else if ( setup->selection!="" ) {
	hist = new HistRootVariable( setup->name, setup->formula, setup->formulaY, NULL, setup->selection, setup->chain );
      }
      else if ( setup->weight!="" ) {
	hist = new HistRootVariable( setup->name, setup->formula, setup->formulaY, setup->weight, "", setup->chain );
      }
      hist->SetHistogram( setup->name, setup->nbinsX, setup->binedgesX, setup->nbinsY, setup->binedgesY );
    }
    else {
      assert(false);
    }
    
    // setup user bin information
    for ( std::vector< SampleBinInfo* >::iterator it=setup->m_bininfo_list.begin(); it!=setup->m_bininfo_list.end(); it++ ) {
      //if ( fVerbosity>=1 )
      std::cout << "Loading bin info, " << (*it)->name << ", into " << setup->name << std::endl;
      TruthBinInfo* binfo = new TruthBinInfo( (*it)->name, (*it)->formula, *it, setup->chain );
      hist->AddBinInfo( (*it)->name, binfo );
    }

    // now create Sample object [ Need to be able to generalize this? ]
    if ( fVerbosity>=1 )
      std::cout << "Creating SampleOscfitROOTPDF " << setup->name << std::endl;
    SampleOscfitROOTPDF* sample = new SampleOscfitROOTPDF( setup->name, setup->chain_name );
    std::cout << "Setting OscBinInfoName: " << setup->bininfo_name_for_osc << std::endl;
    sample->SetOscBinInfoName( setup->bininfo_name_for_osc );
    for ( std::set< std::string >::iterator it_bininfo=setup->bininfo_names_for_splines.begin(); it_bininfo!=setup->bininfo_names_for_splines.end(); it_bininfo++ ) {
      std::cout << "Inserting spline binfo name=" << *it_bininfo << " into sample" << std::endl;
      sample->AddResponseBinInfoName( *it_bininfo );
    }
    sample->SetSampleHistogram( hist );

    // add it to the sample manager
    sampleMan->RegisterSample( setup->name, sample );

  }//end of loop over samples


  sampleMan->SetupBins();

  
}

///-------------------------------------------------------------------------------------------------------------------

TruthBinInfo* SampleSetupParser::GetExampleTruthBinInfo( std::string sample ) {
  if ( m_truthinfo_dict.find( sample )!=m_truthinfo_dict.end() )
    return m_truthinfo_dict[ sample ];
  std::cout << "Did not find sample=" << sample << " in the TruthBinInfo template dictionary." << std::endl;
  assert(false);
}

///-------------------------------------------------------------------------------------------------------------------

void SampleSetupParser::SetVerbosity( int verbose ) {
  fVerbosity = verbose;
}

void SampleSetupParser::PrintVerbosity() {

  std::cout << "--------------------------------------------" << std::endl;
  std::cout << "SampleSetupParser::PrintVerbose" << std::endl;
  std::cout << "current verbose level" << fVerbosity << std::endl;
  std::cout << " [-1]: Quiet" << std::endl;
  std::cout << " [0]: Standard" << std::endl;
  std::cout << " [1]: Print Details of Progress" << std::endl;
  std::cout << " [2]: Puke Details of Parsing" << std::endl;
  std::cout << "--------------------------------------------" << std::endl;

}
