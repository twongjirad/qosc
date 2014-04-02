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
 * ------------------------------------------------------------------------------------
 * \class HistCoordinator: 
 * \ingroup rootvariable
 * \brief Helps organize the filling of HistRootVariables from a ROOT TChain
 *
 * This is a specialized list which only takes HistRootVariable objects.  
 * Such objects need to loop over a ROOT chain to build their Hists.  
 * Objects registed with this class will have their PDFs built using one, coorindinated loop.
 * However, since disk access is expensive, we want to process all PDFs per event to minimize disk access.
 * Also to save time, only those branches which are being used will be accessed.
 *
 * ------------------------------------------------------------------------------------
 */

#ifndef __HistCoordinator__
#define __HistCoordinator__

#include <set>

class TChain;

namespace qosc {

  class HistRootVariable;
  //class PDFCondRootVar;

  class HistCoordinator {

  public:
    HistCoordinator( TChain* common_chain );
    virtual ~HistCoordinator();

    void Add( HistRootVariable* apdf ); ///< insert a PDF element
    //void Add( PDFCondRootVar* apdf ); ///< insert a PDF element
    void Remove( HistRootVariable* apdf ); ///< remove a PDF
    //void Remove( PDFCondRootVar* apdf ); ///< remove a PDF
    void Clear(); ///< clear the list
    void SetChain( TChain* chain ) { m_source_chain = chain; }; ///< Set the common chain
    void BuildPDFs(); ///< Rebuild all PDFs
    int ProcessEntry( int entry );
    void UpdatePDFs(); ///< Update all PDFs
    bool DoAllPDFsShareCommonChain();
    void WriteHistograms();
    void SetMaxEntries( int max ) { fMaxEntries = max; };
    void SetStartEntry( int nstart ) { fStart = nstart; };
    //int GetNumberOfPDFs() { return m_PDFlist.size()+m_CondPDFlist.size(); };
    int GetNumberOfPDFs() { return m_PDFlist.size(); };
    void Print();

  protected:

    TChain* m_source_chain;
    std::set< HistRootVariable* > m_PDFlist;
    //std::set< PDFCondRootVar* > m_CondPDFlist; ///< inelegant, but will work for now
    int fMaxEntries;
    int fStart;
    int bytesRead;

  private:

  };

}

#endif
