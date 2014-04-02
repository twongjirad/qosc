#include "HistCoordinator.hh"
#include <assert.h>
#include <iostream>
#include "TChain.h"
#include "HistRootVariable.hh"
//#include "HistCondRootVar.hh"

using namespace qosc;

HistCoordinator::HistCoordinator( TChain* common_chain ) {
  m_source_chain = common_chain;
  fMaxEntries = -1;
  fStart = 0;
}

HistCoordinator::~HistCoordinator() {
  Clear();
}


void HistCoordinator::Add( HistRootVariable* apdf ) {
  m_PDFlist.insert( apdf );
}

// void HistCoordinator::Add( PDFCondRootVar* apdf ) {
//   m_CondPDFlist.insert( apdf );
// }

void HistCoordinator::Remove( HistRootVariable* apdf ) {
  //m_PDFlist.remove( apdf );
}

// void HistCoordinator::Remove( PDFCondRootVar* apdf ) {
//   //m_CondPDFlist.remove( apdf );
// }

void HistCoordinator::Clear() {
  m_PDFlist.clear();
  //m_CondPDFlist.clear();
}

void HistCoordinator::BuildPDFs() {
  if (DoAllPDFsShareCommonChain()==false) {
    std::cout << "HistCoordinator attempted to build PDFs in list but not all PDFs share the same ROOT tree source!" << std::endl;
    assert(false);
  }

  if (m_PDFlist.size()==0) return;

  
  std::set< HistRootVariable* >::iterator it;
  for (it=m_PDFlist.begin(); it!=m_PDFlist.end(); it++) {
    HistRootVariable* pdf = (*it);
    pdf->ClearHistogram();
  }
//   std::set< PDFCondRootVar* >::iterator it_maps;
//   for (it_maps=m_CondPDFlist.begin(); it_maps!=m_CondPDFlist.end(); it_maps++) {
//     PDFCondRootVar* pdf = (*it_maps);
//     pdf->ClearHistogram();
//   }

  int entry = fStart;
  bytesRead = m_source_chain->GetEntry( entry );
  if ( fStart!=0 ) std::cout << "HistCoordinator::BuildPDF - Starting from entry " << fStart << " (bytes=" << bytesRead << ")" << std::endl;
  if ( fMaxEntries>0 ) std::cout << "HistCoordinator::BuildPDF - Ending at entry " << fMaxEntries << std::endl;

  while (bytesRead>0 && ( fMaxEntries<0 || entry<fMaxEntries ) ) {
    if (entry%10000==0) std::cout << "HistCoordinator::BuildPDF: entry " << entry << " (bytes read: " << bytesRead << ")" << std::endl;

    ProcessEntry( entry );
    bytesRead = m_source_chain->GetEntry( ++entry );
  }

//   for (it_maps=m_CondPDFlist.begin(); it_maps!=m_CondPDFlist.end(); it_maps++) {
//     PDFCondRootVar* pdf = (*it_maps);
//     pdf->BuildMap();
//   }
  
}

int HistCoordinator::ProcessEntry( int entry ) {

  std::set< HistRootVariable* >::iterator it;
  for (it=m_PDFlist.begin(); it!=m_PDFlist.end(); it++) {
    HistRootVariable* pdf = (*it);
    pdf->ProcessEvent( entry );
  }
  
//   std::set< PDFCondRootVar* >::iterator it_maps;
//   for (it_maps=m_CondPDFlist.begin(); it_maps!=m_CondPDFlist.end(); it_maps++) {
//     PDFCondRootVar* pdf = (*it_maps);
//     pdf->ProcessEvent( entry );
//   }

  return 0;
}

void HistCoordinator::UpdatePDFs() {

  if (DoAllPDFsShareCommonChain()==false) {
    std::cout << "HistCoordinator attempted to update PDFs in list but not all PDFs share the same ROOT tree source!" << std::endl;
    assert(false);
  }

  if (m_PDFlist.size()==0) return;

  std::set< HistRootVariable* >::iterator it;
  for (it=m_PDFlist.begin(); it!=m_PDFlist.end(); it++) {
    HistRootVariable* pdf = (*it);
    if ( pdf->GetBuildStatus()==true )
      pdf->ClearHistogram();
  }

  int entry = 0;
  int bytes = m_source_chain->GetEntry( entry );
  while (bytes>0 && ( fMaxEntries<0 || entry<fMaxEntries ) ) {
    
    std::set< HistRootVariable* >::iterator it;
    for (it=m_PDFlist.begin(); it!=m_PDFlist.end(); it++) {
      HistRootVariable* pdf = (*it);
      if ( pdf->GetBuildStatus()==true )
	pdf->ProcessEvent(entry);
    }
    
    bytes = m_source_chain->GetEntry( entry++ );
  }

}

bool HistCoordinator::DoAllPDFsShareCommonChain() {
  std::set< HistRootVariable* >::iterator it;
  for (it=m_PDFlist.begin(); it!=m_PDFlist.end(); it++) {
    HistRootVariable* pdf = (*it);
    if ( m_source_chain!=pdf->GetSourceChain() ) {
      return false;
    }
  }  

//   std::set< PDFCondRootVar* >::iterator it_maps;
//   for (it_maps=m_CondPDFlist.begin(); it_maps!=m_CondPDFlist.end(); it_maps++) {
//     PDFCondRootVar* pdf = (*it_maps);    
//     if ( m_source_chain!=pdf->GetSourceChain() ) {
//       return false;
//     }
//   }  

  return true;
}

void HistCoordinator::WriteHistograms() {
  std::set< HistRootVariable* >::iterator it;
  for (it=m_PDFlist.begin(); it!=m_PDFlist.end(); it++) {
    HistRootVariable* pdf = (*it);
    TH1D* hist = (TH1D*)pdf->GetHistogram();
    hist->Write();
  }

//   std::set< PDFCondRootVar* >::iterator it_maps;
//   for (it_maps=m_CondPDFlist.begin(); it_maps!=m_CondPDFlist.end(); it_maps++) {
//     PDFCondRootVar* pdf = (*it_maps);    
//     pdf->WriteHistogram();
//   }  

}

void HistCoordinator::Print() {
  std::cout << "====================================================================" << std::endl;
  std::cout << "HistCoordinator for chain=" << m_source_chain << " " << m_source_chain->GetName() << std::endl; 
  std::cout << "Number of PDFs: " << GetNumberOfPDFs() << std::endl;
  int npdf = 1;
  for ( std::set< HistRootVariable* >::iterator it=m_PDFlist.begin(); it!=m_PDFlist.end(); it++) {
    std::cout << "(" << npdf << ") " << (*it)->GetName() << " chain=" << (*it)->GetSourceChain() << std::endl;
    npdf++;
  }
  std::cout << "====================================================================" << std::endl;
}
