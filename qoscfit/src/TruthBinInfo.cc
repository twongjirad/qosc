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

#include "TruthBinInfo.hh"
#include <assert.h>
#include <cstdlib>
#include <iostream>
#include <sstream>

#include "TFile.h"

using namespace qosc;

TruthBinInfo::TruthBinInfo( std::string truth_name, std::string root_formula, SampleBinInfo* binfo, TChain* source_chain ) 
  : UserBinInfo( truth_name )
{
  m_truth_name = truth_name;
  m_root_formula = root_formula;
  m_source_chain = source_chain;
  m_bin_def = binfo;
  m_ncutlists = 0;
  SetChain( m_source_chain );
  SetupInstance( binfo );
  Initialize();
  //Print();
}


TruthBinInfo::~TruthBinInfo() {
  // Clean up.
  for ( std::vector< cut_set* >::iterator it=m_cutset_list.begin(); it!=m_cutset_list.end(); it++ ) {
    cut_set* cutset = (*it);
    delete [] cutset->names;
    delete [] cutset->formulas;
    delete cutset;
    cutset = NULL;
  }
  m_cutset_list.clear();

  for (int i=0; i<m_ncutlists; i++) {
    delete m_cut_list[i].hist;
  }
  delete [] m_cut_list;
}

void TruthBinInfo::SetChain( TChain* chain ) {
  m_source_chain = chain;
  if ( m_source_chain ) {
    m_rootvars.SetChain( chain );
    m_rootvars.Add( m_root_formula );
    if ( !fCutsStored ) {
      for (int cut=0; cut<m_ncutlists; cut++)
	m_rootvars.Add( m_cut_list[cut].formula );
    }
  }
}

double TruthBinInfo::GetBinReweight() {
  return 1.0;
}

void TruthBinInfo::FillBinInfo( double external_weight ) {
  for (int n=0; n<m_ncutlists; n++) {
    if ( m_rootvars.GetVariableValue( m_cut_list[n].formula )>0.9 ) {
      if ( ndims==1 ) {
	m_cut_list[n].hist->Fill( m_rootvars.GetVariableValue( m_root_formula ), external_weight );
// 	std::cout << "Fill hist: " << m_rootvars.GetVariableValue( m_root_formula ) << " " << external_weight << std::endl;
// 	std::cin.get();
      }
      else assert( false ); // not yet implemented
      //break; // mutually exclusive?
    }
  }
}

UserBinInfo* TruthBinInfo::Copy( std::string name ) {
  
  // create a copy with the above name

  // things to transfer...
  // (1) template histogram
  // (2) source chain
  // (3) cut variables
  // (4) fill variable formula(s)
  TruthBinInfo* clone = new TruthBinInfo( name, m_root_formula, m_bin_def, m_source_chain );


  return clone;
}


void TruthBinInfo::SetupInstance( SampleBinInfo* binfo ) {
  // This class reads in the bin configuration stored in the SampleBinInfo instance.
  // 
  if ( !binfo->IsCompletelyDefined() ) {
    std::cout << "Defining template instance of " << binfo->name << std::endl;
    std::cout << binfo->name << " is not completely defined." << std::endl;
    assert(false);
  }


  // create instance
  ndims = binfo->ndims;
    
  // set cuts
  for (int cut=0; cut<binfo->ncutsets; cut++) {

    // allocate pointer
    cut_set* cutset = new cut_set;

    // number of cuts in list
    int ncuts = binfo->ncuts.at(cut);
    cutset->ncuts = ncuts;
      
    // transfer cut expressions
    cutset->formulas = new std::string[ ncuts ];
    for (int i=0; i<ncuts; i++) cutset->formulas[i] = binfo->cutlists[ binfo->cuts.at(cut) ][i];

    // transfer cut names
    cutset->names = new std::string[ ncuts ];
    for (int i=0; i<ncuts; i++) cutset->names[i] = binfo->cutnames[ binfo->cuts.at(cut) ][i];
    m_cutset_list.push_back( cutset );
  }//end of loop over cuts

  // make template histogram
  TH1* truth_hist = NULL;
  if ( ndims==1 ) {
	
    // build name
    std::string histname = "h" + binfo->name + "_" + GetInstanceName();
      
    // setup bins
    int totalbins = 0;
    for (int i=0; i<binfo->nXbinsets; i++) totalbins += binfo->nbinsX.at(i);
    double* binedges = new double[ totalbins+1 ];
      
    // set first bin
    binedges[0] = binfo->binedgesX.at(0)[0]; 
    
    // loop through other edges
    int bin = 1;
    double lastbin = 0.;
    for (int i=0; i<binfo->nXbinsets; i++) {
      double binwidth = (binfo->binedgesX.at(i)[1]-binfo->binedgesX.at(i)[0])/double(binfo->nbinsX.at(i));
      for (int j=0; j<binfo->nbinsX.at(i); j++) {
	binedges[bin] = lastbin + binwidth;
	lastbin+=binwidth;
	bin++;
      }
    }//end of loop over bin sets
      
    TH1D* hist1d = new TH1D( histname.c_str(), histname.c_str(), totalbins, binedges );
    truth_hist = (TH1*)hist1d;
    delete [] binedges;
  }// if hist is 1D
  else {
    // 2D hist not yet ready
    std::cout << "not ready for truth::ndims=" << ndims << std::endl;
    assert(false);
  }
    
  // Set the template histogram
  m_template_hist = truth_hist;
    
}

void TruthBinInfo::Initialize() {
  // need to:
  //  (1) setup m_cut_lists/m_hists_names based on cutlists
  //  (2) setup m_hists_lists
  //  (3) setup root variable list
  int ncutsets = m_cutset_list.size();
  std::vector< std::string > cutlist;
  std::vector< std::string > namelist;
  int* index_tracker = new int[ncutsets];
  memset( index_tracker, 0, ncutsets*sizeof(int) );

  bool finished = false;
  while (!finished) {
    // for want to make combination of all unique combo of cuts.
    // using recursion-ish method. the index_tracker array tracks how far we've moved in the loop
    finished = MakeCutlist( cutlist, namelist, index_tracker, ncutsets ); // recursive nightmare
  }
  delete [] index_tracker;

  // OK, now make the array
  m_ncutlists = cutlist.size();
  m_cut_list = new cut_data[ m_ncutlists ];

  //std::cout << "TruthBinInfo[" << m_truth_name << ", " << GetInstanceName() << "] Making truth histograms (total=" << m_ncutlists <<")" << std::endl;
  int nbins = 0;
  for (int i=0; i<m_ncutlists; i++) {
    std::stringstream ss_histname;
    ss_histname << "h" << m_truth_name << "_" << namelist.at(i);// << "_id" << GetInstanceID();
    m_cut_list[i].name = namelist.at(i);
    m_cut_list[i].formula = cutlist.at(i);
    m_cut_list[i].histname = ss_histname.str();
    m_cut_list[i].hist = (TH1*)m_template_hist->Clone( m_cut_list[i].histname.c_str() );
    m_cut_list[i].hist->Reset();
    m_cut_list[i].hist->SetTitle(  m_cut_list[i].formula.c_str() );
    m_cut_list[i].id = i;
    m_cut_list[i].start_bin = nbins;
    m_cut_list[i].nbins = ((TH1D*)m_cut_list[i].hist)->GetNbinsX();
    m_cut_list[i].end_bin = nbins + ((TH1D*)m_cut_list[i].hist)->GetNbinsX() - 1;
    nbins += m_cut_list[i].nbins;

    m_cut_lookup[ namelist.at(i) ] = i;
    //std::cout << m_hists_lists[i]->GetName() << ": " << m_cut_lists[i];
    if ( m_source_chain ) {
      m_rootvars.Add( m_cut_list[i].formula );
      fCutsStored = true;
      //std::cout << "  added to root varlist";
    }
    //std::cout << std::endl;
  }
}

std::string TruthBinInfo::MakeCombinedCutName( std::string sub_cuts[], int nsubcuts ) {
  assert( nsubcuts==m_cutset_list.size() );
  std::string cutname = "";
  for (int i=0; i<nsubcuts; i++) {
    cutname += sub_cuts[i];
    if ( i<nsubcuts-1 ) cutname += "_";
  }
  return cutname;
}

bool TruthBinInfo::MakeCutlist( std::vector< std::string >& cutlist, std::vector< std::string >& namelist, int* index_tracker, int& ncutlists ) {

  
  // make a string with the current indicies.
  std::string cutformula = "";
  std::string name = "";
  for (int i=0; i<ncutlists; i++) {
    cutformula += "(" + m_cutset_list.at(i)->formulas[ index_tracker[i] ] + ")";
    name += m_cutset_list.at(i)->names[ index_tracker[i] ];
    if ( i<ncutlists-1 ) {
      cutformula += " && ";
      name += "_";
    }
  }
  cutlist.push_back( cutformula );
  namelist.push_back( name );

  // push forward the position of the cuts
  for (int i=ncutlists-1; i>=0; i-- ) {
    index_tracker[i]++;
    if ( index_tracker[i]<m_cutset_list.at(i)->ncuts ) {
      break;
    }
    else if ( index_tracker[i]==m_cutset_list.at(i)->ncuts ) {
      if ( i==0 ) return true; // we're done
      // this index is at the end, reset it.
      index_tracker[i] = 0;
    }
    else {
      assert(false); // oops
    }
  }
  
  return false; // not done
}


TH1* TruthBinInfo::GetTruthHist( int index ) { 
  return m_cut_list[index].hist;
}

TH1* TruthBinInfo::GetTruthHist( std::string cutname ) { 
  return m_cut_list[ m_cut_lookup[cutname] ].hist;
}

void TruthBinInfo::Print() {
  std::cout << "----------------------------------------------------------------------------" << std::endl;
  std::cout << "TruthBinInfo Instance: name=" << UserBinInfo::GetInstanceName() << " (" << this << ")" << std::endl;
  std::cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << std::endl;
  for (int n=0; n<GetNumberOfTruthHists(); n++) {
    std::cout << "[" << n << "] " 
	      << GetTruthHist(n)->GetName() 
	      << " (cut=" << m_cut_list[n].name << ") : " 
	      <<  m_cut_list[n].formula << " : "
	      << GetTruthHist(n)->Integral() << std::endl;
  }
  std::cout << "----------------------------------------------------------------------------" << std::endl;
}

int TruthBinInfo::GetCutIndex( std::string cutname ) {
  if ( m_cut_lookup.find( cutname )!=m_cut_lookup.end() ) return m_cut_lookup[cutname];
  std::cout << "TruthBinInfo::GetCutIndex did not find the cutname=" << cutname << std::endl;
  std::cout << " available cuts: " << std::endl;
  for ( std::map< std::string, int >::iterator it=m_cut_lookup.begin(); it!=m_cut_lookup.end(); it++ ) {
    std::cout << "  " << (*it).first << " [" << (*it).second << "]" << std::endl;
  }
  assert(false);
  return -1;
}

int TruthBinInfo::GetCutIndex( std::string subcuts[], int nsubcuts ) {
  std::string cutname = MakeCombinedCutName( subcuts, nsubcuts );
  return GetCutIndex( cutname );
}

int TruthBinInfo::GetCutTruthBinStart( std::string cutname ) {
  return m_cut_list[ GetCutIndex( cutname ) ].start_bin;
}

int TruthBinInfo::GetCutTruthBinStart( std::string subcuts[], int nsubcuts ) {
  return m_cut_list[ GetCutIndex( subcuts, nsubcuts ) ].start_bin;
}

int TruthBinInfo::GetTotalTruthBins() {
  return GetNumberOfTruthHists()*((TH1D*)m_template_hist)->GetNbinsX();
}

int TruthBinInfo::GetNumberOfTruthBinsPerHist() {
  return ((TH1D*)m_template_hist)->GetNbinsX();
}

void TruthBinInfo::GetListOfCutNames( std::vector< std::string >& cutnames ) {
  cutnames.clear();
  for (int cut=0; cut<m_ncutlists; cut++) {
    cutnames.push_back( m_cut_list[cut].name );
  }
}

void TruthBinInfo::GetListOfCutDefinitions( std::vector< std::string >& cutdefs ) {
  cutdefs.clear();
  for (int cut=0; cut<m_ncutlists; cut++) {
    cutdefs.push_back( m_cut_list[cut].formula );
  }
}

std::string TruthBinInfo::GetCutNameFromTruthBin( int bin_num ) {
  for (int cut=0; cut<m_ncutlists; cut++) {
    if ( m_cut_list[cut].start_bin<=bin_num && bin_num<=m_cut_list[cut].end_bin )
      return m_cut_list[cut].name;
  }
  assert(false);
}

std::string TruthBinInfo::GetCutDefinitionFromTruthBin( int bin_num ) {
  for (int cut=0; cut<m_ncutlists; cut++) {
    if ( m_cut_list[cut].start_bin<=bin_num && bin_num<=m_cut_list[cut].end_bin )
      return m_cut_list[cut].formula;
  }
  assert(false);
}

void TruthBinInfo::Reset() {
  for (int i=0; i<GetNumberOfTruthHists(); i++)
    GetTruthHist( i )->Reset();
}

double TruthBinInfo::GetTruthBinContentFromGlobalIndex( int global_truth_bin ) {
  int itruthhist = global_truth_bin/GetNumberOfTruthBinsPerHist();
  int itruth_bin = global_truth_bin%GetNumberOfTruthBinsPerHist();
  return ((TH1D*)GetTruthHist( itruthhist ))->GetBinContent( itruth_bin+1 );
}

void TruthBinInfo::LoadBinInfoFromFile( TFile* input ) {

  for (int i=0; i<GetNumberOfTruthHists(); i++) {
    TH1D* truthhist = (TH1D*)GetTruthHist( i );
    std::string name( truthhist->GetName() );
    size_t pos = name.rfind("__x");
    std::string stem = name.substr(0,pos);
    TH1D* from_file = (TH1D*)input->Get( stem.c_str() );
    if ( !from_file ) {
      std::cout << __PRETTY_FUNCTION__ << "[ " << GetInstanceName() << "]" << " could not load " << stem << std::endl;
      assert(false);
    }
    truthhist->Reset();
    truthhist->Add( from_file );
  }
  
}

void TruthBinInfo::Write() {
  for (int i=0; i<GetNumberOfTruthHists(); i++)
    GetTruthHist( i )->Write();
}
