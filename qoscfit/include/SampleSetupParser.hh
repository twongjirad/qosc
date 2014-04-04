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
 * \class SampleSetupParser
 * \ingroup QoscFit
 * \brief Configures a SampleOscfitROOTPDF object based on setup stored in textfile 
 *
 * Guide to Verbosity
 * -1: Silent
 *  0: Default
 *  1: Details
 *  2: Parsing steps
 * For list of commands, see the source.
 * ---------------------------------------------------------------------------------------------- */


#ifndef __SampleSetupParser__
#define __SampleSetupParser__

#include <string>
#include <map>
#include <set>
#include <vector>
#include <fstream>
#include <iostream>

class TChain;

namespace qosc {

  class SampleManagerROOTPDF;
  class TruthBinInfo;
  class SampleBinInfo;
  class SampleSetupPars;

  class SampleSetupParser {

  public: 
    SampleSetupParser( std::string setupfilename );
    virtual ~SampleSetupParser();

  protected:
  public:
    std::string m_setupfilename;
    std::ifstream* m_input_file;

    void OpenFile( std::string filename );
    void CloseFile();
    int GetNextCommand( std::string& command, std::string& arg );
    std::string GetNextLine();
    void ParseSetup();

    void SetupSamples( SampleManagerROOTPDF* sampleMan ); ///< Loads sample manager with SampleOscfitROOTPDF instances as specified by setup file
    TruthBinInfo* GetExampleTruthBinInfo( std::string sample );
    bool WasSampleDefined( std::string samplename ) { return m_samples_defined.find( samplename )!=m_samples_defined.end(); };

  protected:

    // Commands
    static const int kNumCommands = 18;
    static std::string __commands__[kNumCommands];
    static std::set< std::string > __command_set__;

    void ParseCommand_DEF_CHAIN( std::string command_line );

    void ParseCommand_DEF_SAMPLE( std::string command_line );
    void ParseCommand_DEF_SAMPLE_BINS( std::string command_line );
    void ParseCommand_DEF_SAMPLE_VAR( std::string command_line );
    void ParseCommand_DEF_SAMPLE_CUT( std::string command_line );
    void ParseCommand_DEF_SAMPLE_WGT( std::string command_line );
    void ParseCommand_ADD_SAMPLE_INFO( std::string command_line );
    void ParseCommand_ADD_SAMPLE_CHAIN( std::string command_line );
    void ParseCommand_DEF_SAMPLE_OSC_BININFO( std::string command_line );
    void ParseCommand_DEF_SAMPLE_SPLINE_BININFO( std::string command_line );

    void ParseCommand_DEF_INFO_BINS( std::string command_line );
    void ParseCommand_DEF_INFO_VAR( std::string command_line );
    void ParseCommand_DEF_INFO_CUTS( std::string command_line );
    void ParseCommand_DEF_INFO_CUTS2( std::string command_line );
    void ParseCommand_DEF_INFO_CUTNAMES( std::string command_line );

    void ParseCommand_OPT_SPLINE_MERGEBINS( std::string command_line );

    void ParseCommand_SET_PARSER_VERBOSITY( std::string command_line );

  protected:
    void getNumberList( std::string command_line, size_t& pos, int& nnums, double*& numbers );
    void getStringList( std::string command_line, size_t& pos, int nnums, std::string *texts );
    void parseStringList( std::string command_line, size_t& pos, int nnums, std::vector< std::string >& list );
  public:
    std::set< std::string > m_samples_defined; //< Names of samples defined in the setup file
    std::map< std::string, TChain* > m_chain_dict; //< Container for all the chains loaded
    std::map< std::string, SampleBinInfo* > m_bininfo_dict; //< Stores the truth binning defined by DEF_INFO_BINS
    std::map< std::string, TruthBinInfo* > m_truthinfo_dict; //< Stores the truth bins object, stores the bins for each cut
    std::map< std::string, SampleSetupPars* > m_sample_setups; //< Stores the class which contains configuration parameters for a sample
    std::map< std::string, SampleBinInfo* > m_splinebins_dict; //< Stores the bin definition used to merge spline bins

  protected:
    int fVerbosity;
  public:
    void SetVerbosity( int verbose );
    void PrintVerbosity();
    void PrintSetup();

  };

}

#endif
