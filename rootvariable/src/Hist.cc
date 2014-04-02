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

#include "Hist.hh"
#include <iostream>
#include <assert.h>
#include <sstream>
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"

#include "UserBinInfoList.hh"



using namespace qosc;

Hist::Hist( std::string name ) {
  m_hist_pdf = NULL;
  m_hist_pdf_stored = NULL;
  m_hist_pdf_2D = NULL;
  m_hist_pdf_2D_stored = NULL;
  fHistogramStatus = kUndefined;
  fUseOverUnderFlow = false;
  fNDims = -1;
  m_name = name;
  SetVerbose( 0 );
}

Hist::~Hist() {
  DestroyHistogram();
}

// --------------------------------------------------------------------------------------------------
// Setup 1D Histogram

void Hist::SetHistogram( TH1D* hist ) { 
  m_hist_pdf = hist; 
  SetHistogramStatus( kDefined );
  SetHistDimensions( 1 );
  SetupUserBinInfo();
}

void Hist::SetHistogram( std::string histname, int nbins, double xmin, double xmax ){
  // Define a new 1D histogram
  DestroyHistogram();
  SetHistogram( new TH1D( histname.c_str(), histname.c_str(), nbins, xmin, xmax ) ); 
}

void Hist::SetHistogram( std::string histname, int nbins, double* binedges ) {
  // Define a new 1D histogram
  DestroyHistogram();
  SetHistogram( new TH1D( histname.c_str(), histname.c_str(), nbins, binedges ) );
}

// --------------------------------------------------------------------------------------------------
// 2D Functions

void Hist::SetHistogram( std::string histname, int nbinsX, double xmin, double xmax, int nbinsY, double ymin, double ymax ){
  // Define a new 1D histogram
  if (m_hist_pdf_2D!=NULL) {
    ClearBinInfo();
    delete m_hist_pdf;
  }
  fNDims=2;
  m_hist_pdf_2D = new TH2D( histname.c_str(), histname.c_str(), nbinsX, xmin, xmax, nbinsY, ymin, ymax );
  fHistogramStatus = kDefined;
  for (int bx=0; bx<=nbinsX+1; bx++) {
    for (int by=0; by<=nbinsY+1; by++) {
      char binlabel[200];
      sprintf( binlabel, "%s_bin%d", GetName().c_str(), GetBinIndex(bx,by) );
      m_bininfo[GetBinIndex(bx,by)] = new UserBinInfoList( binlabel );
    }
  }
}

void Hist::SetHistogram( std::string histname, int nbinsX, double* binedgesX, int nbinsY, double* binedgesY ) {
  if (m_hist_pdf_2D!=NULL) {
    ClearBinInfo();
    delete m_hist_pdf_2D;
  }
  fNDims = 2;
  m_hist_pdf_2D = new TH2D( histname.c_str(), histname.c_str(), nbinsX, binedgesX, nbinsY, binedgesY );
  fHistogramStatus = kDefined;
  for (int bx=0; bx<=nbinsX+1; bx++) {
    for (int by=0; by<=nbinsY+1; by++) {
      m_bininfo[GetBinIndex(bx,by)] = new UserBinInfoList;
    }
  }
  std::cout << "Defined 2D Histogram for " << GetName() << ", nbinsX=" << nbinsX << ", nbinsY=" << nbinsY << std::endl;
}

// --------------------------------------------------------------------------------------------------
// Dimension Agnostic

TH1* Hist::GetHistogram() {
  int dims = GetHistDimensions();
  if ( dims==1 ) return m_hist_pdf;
  else if ( dims==2 ) return m_hist_pdf_2D;
  else return NULL;
}

void Hist::CopyHistogram( std::string histname, Hist* origpdf ) {
  CopyHistogram( histname, origpdf->GetHistogram(), origpdf->GetHistDimensions() );
}

void Hist::CopyHistogram( std::string histname, TH1* horig, int dims )  {
  // Warning. The UserBin Info classes are not getting copied.
  // This needs to get fixed.
  SetHistDimensions( dims );
  if (GetHistDimensions()==1) {
    int nbins = ((TH1D*)horig)->GetNbinsX();
    double binedges[nbins+1];
    for (int bin=0; bin<=nbins; bin++)
      binedges[bin] = horig->GetBinLowEdge(bin+1);
    SetHistogram(histname, nbins, binedges );
    for (int bin=0; bin<=nbins; bin++)
      GetHistogram()->SetBinContent( bin, ((TH1D*)horig)->GetBinContent(bin) );
  }
  else if ( GetHistDimensions()==2 ) {
    TH2D* hcopy = (TH2D*)horig->Clone( histname.c_str() );
    for (int bx=0; bx<=hcopy->GetNbinsX()+1; bx++) {
      for (int by=0; by<=hcopy->GetNbinsY()+1; by++) {
	hcopy->SetBinContent( bx, by, ((TH2D*)horig)->GetBinContent( bx, by ) );
      }
    }
    SetHistogram( hcopy );
  }
  else {
    std::cout << "Tried to copy histogram, but dimensions are invalid: " << GetHistDimensions() << std::endl;
    assert(false);
  }
}

// --------------------------------------------------------------------------------------------------
// Histogram Management

void Hist::DestroyHistogram() {
  ClearBinInfo();
  if ( m_hist_pdf ) delete m_hist_pdf;
  m_hist_pdf = NULL;
  if ( m_hist_pdf_2D ) delete m_hist_pdf_2D;
  m_hist_pdf_2D = NULL;
  if ( m_hist_pdf_stored ) delete m_hist_pdf_stored;
  if ( m_hist_pdf_2D_stored ) delete m_hist_pdf_2D_stored;
  m_hist_pdf_stored = NULL;
  m_hist_pdf_2D_stored = NULL;
  m_integral = 0;
  fNDims = -1;
  fHistogramStatus = kUndefined;
}

// --------------------------------------------------------------------------------------------------
// Bin Indexing

int Hist::GetTotalNumberOfBins() {
  // includes overflow and underflow
  if ( GetHistDimensions()==1 ) {
    return ((TH1D*)GetHistogram())->GetNbinsX()+2;
  }
  else if ( GetHistDimensions()==2 ) {
    return (((TH2D*)GetHistogram())->GetXaxis()->GetNbins()+2)*(((TH2D*)GetHistogram())->GetYaxis()->GetNbins()+2);
  }
  return 0;
}

int Hist::GetNumberOfBins() {
  if ( fUseOverUnderFlow ) return GetTotalNumberOfBins();
  
  if ( GetHistDimensions()==1 ) {
    return ((TH1D*)GetHistogram())->GetNbinsX();
  }
  else if ( GetHistDimensions()==2 ) {
    return (((TH2D*)GetHistogram())->GetXaxis()->GetNbins())*(((TH2D*)GetHistogram())->GetYaxis()->GetNbins());
  }
  return 0;
}


void Hist::GetNumberOfBinsXY( int& nbinsx, int& nbinsy ) {
  if ( GetHistDimensions()==1 ) {
    nbinsy = 0;
    if ( fUseOverUnderFlow ) nbinsx = m_hist_pdf->GetNbinsX()+2;
    else nbinsx = m_hist_pdf->GetNbinsX();
  }
  else if ( GetHistDimensions()==2 ) {
    if ( fUseOverUnderFlow ) {
      nbinsx = m_hist_pdf_2D->GetNbinsX()+2;
      nbinsy = m_hist_pdf_2D->GetNbinsY()+2;
    }
    else {
      nbinsx = m_hist_pdf_2D->GetNbinsX();
      nbinsy = m_hist_pdf_2D->GetNbinsY();
    }  
  }
  else {
    assert(false);
  }

}

void Hist::GetTotalNumberOfBinsXY( int& nbinsx, int& nbinsy ) {
  bool underflowflag = fUseOverUnderFlow;
  fUseOverUnderFlow = true;
  GetNumberOfBinsXY( nbinsx, nbinsy );
  fUseOverUnderFlow = underflowflag;
}

int Hist::GetBinIndex( int binx, int biny ) {

  assert( GetHistDimensions()!=-1 );

  if ( GetHistDimensions()==1 ) {
    return binx;
  }

  // below assumes 2D
  if ( biny==-1 ) return binx;

  int nbinsX, nbinsY;
  GetNumberOfBinsXY( nbinsX, nbinsY );
  
  // Check index range
  if ( ( binx<0 || binx>m_hist_pdf_2D->GetNbinsX()+1 ) || ( biny<0 || biny>m_hist_pdf_2D->GetNbinsY()+1) ) {
    std::cout << "Hist::GetBinIndex -- WARNING! Bin indices (" << binx << ", " << biny << ") out of range:"
	      << " x:[" << 0 << ", " << m_hist_pdf_2D->GetNbinsX()+1 << "] y:[" << 0 << ", " << m_hist_pdf_2D->GetNbinsY()+1 << "]" << std::endl;
    assert(false);
  }
  int index = biny*(nbinsX) + binx;
  return index;
}

void Hist::GetBinXYFromIndex( int binindex, int& binx, int& biny ) {

  int nbinsx, nbinsy;
  GetNumberOfBinsXY( nbinsx, nbinsy );

  if ( binindex<0 || binindex>(unsigned int)nbinsx*nbinsy ) {
    std::cout << "Bin Index out of range" << std::endl;
    std::cout << " given binindex=" << binindex << std::endl;
    std::cout << " nbinsx = " << nbinsx << std::endl;
    std::cout << " nbinsy = " << nbinsy << std::endl;
    std::cout << " totalbins: " << nbinsx*nbinsy << std::endl;
    assert(false);
  }
  biny = binindex/nbinsx;
  binx = binindex%nbinsx;
}

double Hist::GetBinContent( int binx, int biny ) {

  if ( GetHistDimensions()==1 ) {
    if ( fUseOverUnderFlow ) return m_hist_pdf->GetBinContent(  binx );
    else return m_hist_pdf->GetBinContent(  binx+1 );
  }
  else if ( GetHistDimensions()==2 ) {
    if (biny==-1) {
      int bx, by;
      GetBinXYFromIndex( binx, bx, by );
      if ( fUseOverUnderFlow ) return m_hist_pdf_2D->GetBinContent( bx, by );
      else return m_hist_pdf_2D->GetBinContent( bx+1, by+1 );
    }
    else {
      if ( fUseOverUnderFlow ) return m_hist_pdf_2D->GetBinContent( binx, biny );
      else return m_hist_pdf_2D->GetBinContent( binx+1, biny+1 );
    }
  }
  else {
    assert(false);
  }

}

void Hist::SetBinContent( int bin, double value ) {

  if ( GetHistDimensions()==1 ) {
    if ( fUseOverUnderFlow ) m_hist_pdf->SetBinContent(  bin, value );
    else m_hist_pdf->SetBinContent(  bin+1, value );
  }
  else if ( GetHistDimensions()==2 ) {
    int bx, by;
    GetBinXYFromIndex( bin, bx, by );
    if ( fUseOverUnderFlow ) m_hist_pdf_2D->SetBinContent( bx, by, value );
    else m_hist_pdf_2D->SetBinContent( bx+1, by+1, value );
  }
  else {
    assert(false);
  }
}

void Hist::SetBinContent( int binx, int biny, double value ) {
  if ( GetHistDimensions()==2 ) {
    if ( fUseOverUnderFlow ) m_hist_pdf_2D->SetBinContent( binx, biny, value );
    else m_hist_pdf_2D->SetBinContent( binx+1, biny+1, value );
  }
  else {
    assert(false);
  }
}

double* Hist::GetBinAddress( int binx, int biny ) {
  
  if ( GetHistDimensions()==1 ) {
    if ( !fUseOverUnderFlow ) ++binx; // skip the underflow bin
    return (&m_hist_pdf->fArray[binx]); // this is pointer math
  }
  else if ( GetHistDimensions()==2 ) {
    if (biny==-1) { // access 2D bin by unrolled index
      int bx, by;
      GetBinXYFromIndex( binx, bx, by );
      if ( !fUseOverUnderFlow ) {
	++bx;
	++by;
      }
      return &m_hist_pdf_2D->fArray[ m_hist_pdf_2D->GetBin( bx, by ) ];
    }
    else { // access 2D bin by ordered pair
      if ( !fUseOverUnderFlow ) { ++binx; ++biny; };
      return &m_hist_pdf_2D->fArray[ m_hist_pdf_2D->GetBin( binx, biny ) ];
    }
  }

  // should never get here
  assert(false);
  return NULL;
}


// --------------------------------------------------------------------------------------------------
// Bin Info Functions

void Hist::SetupUserBinInfo() {
  assert( GetHistDimensions()!=-1 );
  assert( GetHistogramStatus()!=kUndefined );
  for (int bin=0; bin<GetTotalNumberOfBins(); bin++) {
    std::stringstream ss;
    ss << GetName() << "_bin" << bin;
    m_bininfo[bin] = masterList.Copy( ss.str() );
  }
}

void Hist::ClearBinInfo() {
  // clear bin info lists
  if ( GetHistogramStatus()==Hist::kUndefined ) return;
  
  for (int bin=0; bin<GetTotalNumberOfBins(); bin++) {
    delete m_bininfo[bin];
    m_bininfo[bin] = NULL;
  }
  m_bininfo.clear();
}

void Hist::ProcessBinInfo( int binfilled , double weight) {
  UserBinInfoList* binlist = m_bininfo[binfilled];
  if ( !binlist ) {
    std::cout << GetName() << " failed to have a bin info list at bin " << binfilled << std::endl;
    return;
  }
  binlist->FillInfo( weight );
}

void Hist::ProcessBinInfo( int binX, int binY, double weight ) {
  int index = GetBinIndex(binX,binY);
  UserBinInfoList* binlist =   m_bininfo[index];
  if ( !binlist ) return;
  binlist->FillInfo( weight );
}

void Hist::AddBinInfo( std::string name, UserBinInfo* bininfo ) {

  masterList.AddInfo( name, bininfo );

  if ( GetHistogramStatus()==kUndefined ) {
    std::cout <<  "Hist::AddBinInfo: Warning - UserBinInfo '" << name << "' was added, but histogram has not been setup yet to define the bins." << std::endl;
    return;
  }

  // we make a copy for every bin
  for (unsigned int n=0; n<m_bininfo.size(); n++) {
    UserBinInfoList* infolist = m_bininfo[n];
    if ( !infolist ) {
      std::cout << "bin info list not yet initialized!" << std::endl;
      assert(false);
    }
    UserBinInfo* copy = bininfo->Copy( infolist->GetBinLabel()+"_"+name );
    //std::cout << " Adding " << name << " (" << copy->GetInstanceName() << ") to " << infolist->GetBinLabel() << std::endl;
    m_bininfo[n]->AddInfo( name, copy );
  }
  
}


UserBinInfoList* Hist::GetBinInfoList( int binX, int binY ) {
  if ( !fUseOverUnderFlow ) {
    binX++;
    if ( binY>-1 )
      binY++;
  }
  return GetBinInfoListInternal( binX, binY );
}

UserBinInfoList* Hist::GetBinInfoListInternal( int binX, int binY ) {
  // When dealing with bin info classes. we include over and underflow. beacause of this, we must take into account the over flow or underflow
  int binindex = GetBinIndex( binX, binY );

  std::map< int, UserBinInfoList* >::iterator it = m_bininfo.find( binindex );
  if ( it==m_bininfo.end() ) assert(false);
  UserBinInfoList* binlist = m_bininfo[binindex];
  return binlist;
}
  

UserBinInfo* Hist::GetBinInfo( std::string name, int binX, int binY ) {
  // operates for 1D and 2D functions. if 1D, ignore binY
  UserBinInfoList* binlist = GetBinInfoList( binX, binY );
  if ( !binlist ) return NULL;
  return binlist->GetInfo( name );
}

void Hist::CopyBinInfo( Hist* source ) {
  if ( GetHistogramStatus()==kUndefined ) return;
  assert( GetTotalNumberOfBins()==source->GetTotalNumberOfBins() );

  for (int bin=0; bin<GetTotalNumberOfBins(); bin++) {
    UserBinInfoList* infolist = source->GetBinInfoListInternal( bin );
    if ( infolist!=NULL ) {
      m_bininfo[bin]->CopyBinInfoClasses( infolist );
    }
  }
}


void Hist::WriteBinInfo() {
  if ( GetHistogramStatus()==kUndefined ) return;
  for ( unsigned int n=0; n<m_bininfo.size(); n++) {
    UserBinInfoList* infolist = GetBinInfoListInternal(n);
    if ( infolist ) infolist->Write();
  }
}

void Hist::ResetAllBinInfo() {
  for (int i=0; i<GetNumberOfBins(); i++) {
    UserBinInfoList* infolist = GetBinInfoList( i );
    for ( UserBinInfoListIter it=infolist->Begin(); it!=infolist->End(); it++ ) {
      (*it).second->Reset();
    }
  }
}

void Hist::ResetBinInfo( std::string name ) {
  for (int i=0; i<GetNumberOfBins(); i++) {
    UserBinInfo* info = GetBinInfo( name, i );
    info->Reset();
  }
}

// --------------------------------------------------------------------------------------------------
// Hist Info Functions

void Hist::AddHistInfo( std::string name, UserBinInfo* info ) {
  m_UserHistInfo.AddInfo( name, info );
}

UserBinInfo* Hist::GetHistInfo( std::string name ) {
  return m_UserHistInfo.GetInfo( name );
}

void Hist::ProcessHistInfo( double weight ) {
  m_UserHistInfo.FillInfo( weight );
}
  

void Hist::WriteHistInfo() {
  if ( GetHistogramStatus()==kUndefined ) return;
  m_UserHistInfo.Write();
}

// --------------------------------------------------------------------------------------------------
// Histogram State Storage

void Hist::StoreHistogram() {
  // we save the histogram for this pdf.
  TH1* myhist = GetHistogram();
  if ( !myhist ) {
    std::cout << "Hist Histogram, " << GetName() << ", was not initialized?! " << myhist << std::endl;
    assert(false);
  }
  std::string histname = std::string(myhist->GetName()) + "_stored";
  //std::cout << "Making clone: " << histname << std::endl;
  if ( GetHistDimensions()==1 ) {
    if ( m_hist_pdf==NULL ) assert(false);

    if ( m_hist_pdf_stored==NULL ) m_hist_pdf_stored = (TH1D*)m_hist_pdf->Clone( histname.c_str() );
    else m_hist_pdf_stored->Reset();

    for (int bin=0; bin<=m_hist_pdf->GetNbinsX(); bin++) {
      m_hist_pdf_stored->SetBinContent( bin, m_hist_pdf->GetBinContent(bin) );
    }
    //std::cout << "original=" << m_hist_pdf->Integral() << " - clone=" << m_hist_pdf_stored->Integral() << std::endl;
  }
  else if ( GetHistDimensions()==2 ) {
    if ( m_hist_pdf_2D==NULL ) assert(false);
    if ( m_hist_pdf_2D_stored==NULL ) {
      m_hist_pdf_2D_stored = (TH2D*)m_hist_pdf_2D->Clone( histname.c_str() );
    }
    else
      m_hist_pdf_2D_stored->Reset();
    for (int bx=0; bx<=m_hist_pdf_2D->GetNbinsX()+1; bx++)
      for (int by=0; by<=m_hist_pdf_2D->GetNbinsY()+1; by++)
	m_hist_pdf_2D_stored->SetBinContent( bx, by, m_hist_pdf_2D->GetBinContent( bx, by ) );
  }
}

TH1* Hist::GetStoredHistogram() {
  if ( GetHistDimensions()==1 && m_hist_pdf_stored ) return m_hist_pdf_stored;
  else if ( GetHistDimensions()==2 && m_hist_pdf_2D_stored ) return m_hist_pdf_2D_stored;
  return NULL;
}

// --------------------------------------------------------------------------------------------------
// PDF Update Functions

void Hist::UpdatePDF() {
  UpdatePDF( kReWeight );
}

void Hist::UpdatePDF( UpdateMode mode ) {
  // In future I might need to add a new mode, where I reweight from current or update reweight from stored.
  if ( mode!=kReWeight ) return;

  if ( GetVerbose()>0 ) {
    std::cout << "----------------------------------------------------------------------------------------" << std::endl;
    std::cout << "Hist::UpdatePDF() -- " << GetName() << " (mode=";
    if ( mode==kReWeight )
      std::cout << mode << " Reweight from stored.";
    std::cout << ")" << std::endl;
  }
  
  if ( GetHistDimensions()==1 ) {
    TH1D* hist = m_hist_pdf;
    if ( m_hist_pdf_stored!=NULL )
      hist = m_hist_pdf_stored;

    for (int bin=0; bin<=m_hist_pdf->GetNbinsX()+1; bin++) {
      double reweight = 1.0;
      UserBinInfoList* infolist = GetBinInfoListInternal(bin);
      if ( infolist==NULL ) continue;
      if ( GetVerbose() ) std::cout << " reweight with instances UserBinInfoList, " << infolist->GetBinLabel() << std::endl;
      UserBinInfoListIter it;
      for (it=infolist->Begin(); it!=infolist->End(); it++) {
	if ( GetVerbose() ) std::cout << "  reweighting factor from " << (*it).second->GetInstanceName() << ": " << (*it).second->GetBinReweight() << std::endl;
	reweight *= (*it).second->GetBinReweight();
      }
      if ( GetVerbose() ) std::cout << " [total reweight=" << reweight << "]" << std::endl;
      m_hist_pdf->SetBinContent( bin, reweight*hist->GetBinContent(bin) );
    }
  }
  else if ( GetHistDimensions()==2 ) {
    TH2D* hist = m_hist_pdf_2D;
    if ( m_hist_pdf_2D_stored!=NULL )
      hist = m_hist_pdf_2D_stored;

    for (int bx=0; bx<=m_hist_pdf_2D->GetNbinsX()+1; bx++) {
      for (int by=0; by<=m_hist_pdf_2D->GetNbinsY()+1; by++) {
	double reweight = 1.0;
	UserBinInfoList* infolist = GetBinInfoListInternal(bx,by);
	if ( infolist==NULL ) continue;
	if ( GetVerbose() ) std::cout << " reweight with instances UserBinInfoList, " << infolist->GetBinLabel() << std::endl;

	UserBinInfoListIter it;
	for (it=infolist->Begin(); it!=infolist->End(); it++) {
	  if ( GetVerbose() ) std::cout << "  reweighting factor from " << (*it).second->GetInstanceName() << ": " << (*it).second->GetBinReweight() << std::endl;
	  reweight *= (*it).second->GetBinReweight();
	}
	if ( GetVerbose() ) std::cout << " [total reweight=" << reweight << "]" << std::endl;
	m_hist_pdf_2D->SetBinContent( bx, by, reweight*hist->GetBinContent(bx,by) );
      }
    }
  }

  if ( GetVerbose() ) {
    std::cout << "End of Hist::UpdatePDF()" << std::endl;
    std::cout << "----------------------------------------------------------------------------------------" << std::endl;
  }
}

// --------------------------------------------------------------------------------------------------
// Load Histograms from File (allows for some time savings)

bool Hist::LoadHistogramsFromFile( TFile* file, double scalefactor ) {
  // this fills the histogram using one found in a file.
  // it must have the same name, number of bins else the program will stop.
  // for this to make sense, it must also have the same binning, but I am not going to try and check this and any other errors.
  int ndims = GetHistDimensions();
  if ( fHistogramStatus==kUndefined ) {
    std::cout << "Hist::LoadHistogramFromFile -- WARNING! The histogram has not been defined yet. Returning without doing anything." << std::endl;
    return false;
  }

  // Load the main histogram
  std::string histname = GetHistogram()->GetName();
  TH1* hist = (TH1*)file->Get( histname.c_str() );
  if ( hist==NULL ) {
    std::cout << "Hist::LoadHistogramFromFile -- WARNING! The histogram '" << histname << "' was not found in file. Tried to load for Hist=" << GetName() << std::endl;
  }
  else {
    if ( ndims==1 ) {
      TH1D* hist1d = (TH1D*)file->Get( histname.c_str() );
      if ( hist1d->GetNbinsX()!=((TH1D*)GetHistogram())->GetNbinsX() ) assert(false);
      ((TH1D*)GetHistogram())->Reset();
      ((TH1D*)GetHistogram())->Add( hist1d, scalefactor );
      if ( GetVerbose() ) std::cout << "Hist[" << GetName() << "]::LoadHistogramFromFile. Successfully loaded " << hist1d->GetName() << std::endl;
    }
    else if ( ndims==2 ) {
      TH2D* hist2d = (TH2D*)file->Get(histname.c_str() );
      if ( hist2d->GetXaxis()->GetNbins()!=((TH2D*)GetHistogram())->GetXaxis()->GetNbins()
	   || hist2d->GetYaxis()->GetNbins()!=((TH2D*)GetHistogram())->GetYaxis()->GetNbins() ) assert(false);
      ((TH2D*)GetHistogram())->Reset();
      ((TH2D*)GetHistogram())->Add(hist2d, scalefactor);
      if ( GetVerbose() ) std::cout << "Hist[" << GetName() << "]::LoadHistogramFromFile. Successfully loaded " << hist2d->GetName() << std::endl;
    }
  }

  // Load the stored histogram
  histname = std::string(GetHistogram()->GetName())+std::string("_stored");
  hist = (TH1*)file->Get( histname.c_str() );
  if ( hist==NULL ) {
    std::cout << "Hist::LoadHistogramFromFile -- WARNING! The histogram '" << histname << "' was not found in file. Tried to load for Hist=" << GetName() << std::endl;
  }
  else {
    StoreHistogram(); // makes a stored histogram if one hasn't already been made.
    if ( ndims==1 ) {
      TH1D* hist1d = (TH1D*)file->Get( histname.c_str() );
      if ( hist1d->GetNbinsX()!=m_hist_pdf_stored->GetNbinsX() ) assert(false);
      m_hist_pdf_stored->Reset();
      m_hist_pdf_stored->Add( hist1d, scalefactor );
      if ( GetVerbose() ) std::cout << "Hist[" << GetName() << "]::LoadHistogramFromFile. Successfully loaded " << hist1d->GetName() << std::endl;
    }
    else if ( ndims==2 ) {
      TH2D* hist2d = (TH2D*)file->Get(histname.c_str() );
      if ( hist2d->GetXaxis()->GetNbins()!=m_hist_pdf_2D_stored->GetXaxis()->GetNbins()
	   || hist2d->GetYaxis()->GetNbins()!=m_hist_pdf_2D_stored->GetYaxis()->GetNbins() ) assert(false);
      m_hist_pdf_2D_stored->Reset();
      m_hist_pdf_2D_stored->Add(hist2d, scalefactor);
      if ( GetVerbose() ) std::cout << "Hist[" << GetName() << "]::LoadHistogramFromFile. Successfully loaded " << hist2d->GetName() << std::endl;
    }
  }

  // Pass the file over to the UserBinInfo classes
  if ( GetVerbose() ) std::cout << "Load UserBinInfo objects from file (total: " << m_bininfo.size() << ")" << std::endl;
  for (unsigned int n=0; n<m_bininfo.size(); n++) {
    UserBinInfoList* infolist = m_bininfo[n];
    for ( UserBinInfoListIter it=infolist->Begin(); it!=infolist->End(); it++ ) {
      (*it).second->LoadBinInfoFromFile( file );
    }
  }

  // Pass the file over to the UserBinInfo classes
  if ( GetVerbose() ) std::cout << "Load HistInfo objects from file (total: " << m_bininfo.size() << ")" << std::endl;
  for ( UserBinInfoListIter it=GetHistInfoListBegin(); it!=GetHistInfoListEnd(); it++ ) {
    (*it).second->LoadBinInfoFromFile( file );
  }

  return true;
}

// --------------------------------------------------------------------------------------------------
// Misc. Tools

void Hist::Print() {
  std::cout << "========================================================" << std::endl;
  std::cout << "Hist Summary." << std::endl;
  std::cout << "Name: " << GetName() << " Address=" << this << std::endl;
  std::cout << "Histogram Status: ";
  ( GetHistogramStatus()==kDefined ) ? std::cout << "Defined" << std::endl : std::cout << "Undefined" << std::endl;
  if ( GetHistogramStatus()==kDefined ) {
    std::cout << "Histogram name: " << GetHistogram()->GetName() 
	      << ", Dimensions=" << GetHistDimensions() 
	      << ", Number of bins: " << GetNumberOfBins() 
	      << std::endl;
    ( fUseOverUnderFlow ) ? std::cout << "Indexes under- and overflow bins" << std::endl : std::cout << "Indices skip under and over-flow bins" << std::endl;
    std::cout << "Number of UserBinInfo class types associated to each bin: " << masterList.GetEntries() << std::endl;
    if ( masterList.GetEntries()>0 ) {
      for ( UserBinInfoListIter it=MasterBinInfoListBegin(); it!=MasterBinInfoListEnd(); it++ ) {
	std::cout << "- - - - -  - - - -  - - - -  - - - - " << std::endl;
	std::cout << "[" << (*it).first << " ] " << (*it).second->GetInstanceName() << std::endl;
	(*it).second->Print();
      }
      std::cout << "- - - - -  - - - -  - - - -  - - - - " << std::endl;
    }
  }
  std::cout << "========================================================" << std::endl;
}

void Hist::PrintBinInfoLists() {
  // dumps information on bin info list of histogram
  if ( GetHistogramStatus()==kUndefined ) return;

  int binx, biny;
  std::cout << "------------------------------------------------------------------------" << std::endl;
  std::cout << " Bin Info Lists for hist=" << GetName() << " (" << GetName() << ")" << std::endl;
  for (unsigned int n=0; n<m_bininfo.size(); n++) {
    std::cout << " -- bin " << n;
    if ( GetHistDimensions()==1 ) std::cout << " [x=" << n << "] --" << std::endl;
    else if ( GetHistDimensions()==2 ) {
      int nbinsx, nbinsy;
      GetTotalNumberOfBinsXY( nbinsx, nbinsy );
      binx = n/nbinsx;
      biny = n%nbinsx;
      std::cout << " [x=" << binx << ",y=" << biny << "] --" << std::endl;
    }
    
    UserBinInfoList* infolist = GetBinInfoListInternal(n);
    if ( infolist==NULL ) {
      std::cout << "  None Define" << std::endl;
      continue;
    }
    UserBinInfoListIter it;
    for (it = infolist->Begin(); it!=infolist->End(); it++) {
      std::cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - -" << std::endl;
      std::cout << "    " << (*it).first << " : " << (*it).second << std::endl;
      (*it).second->Print();
      std::cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - -" << std::endl;
    }
  }
  std::cout << "------------------------------------------------------------------------" << std::endl;
}
