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
 * --------------------------------------------------------------------------------------------
 * \class Hist
 * \ingroup rootvariable
 * \brief Base class for Histogram classes. Extension of the ROOT histogram class.
 *
 *
 * [UserBinInfo Objects]
 * This class also has the ability to define UserBinInfo instances.
 * The purpose of this class is to allow the user to store information associated to each bin.
 * The user than can use this information to reweight bin-by-bin using the
 *   associated UserBinInfo instances.
 * If UpdatePDF( kReWeight ) is called, the UserBinInfo classes will be asked to provide a weight for 
 *  each bin.  If UpdatePDF( kNoReWeight ) is called, the userbininfo will be ignroed.
 * By default, the reweighting is compounded.  However, if SaveHistogram is called, a copy of
 *  the stored histogram is stored and all updates that use the user bin info weights will
 *  reweight the stored histogram contents in order to re-fill the histogram.
 *
 * Dimensions limited in principle to 3D, the max dimensions for ROOT histograms.
 * Higher dimensions will take some effort. But this is neutrino physics, there is barely enough stats
 * for 2D let alone N>2.
 *
 * [BIN INDEXING]
 * For setuping analyses, one must constantly be tracking bins. This leads to two choices that
 *  need to be made:
 *  (1) Choice of first index (0 or 1)
 *  (2) Include overflow and underflow in binning scheme
 * For setting up the analysis, I am of the opinion that using a 0-indexed scheme is more
 *  convenient. Of course, this means overriding ROOT's indexing scheme, which can lead
 *  to confusion.
 * Also, in setuping up an analysis, one must unroll the bins when passing them to a fitter.
 * Because, I must devise a bin numbering scheme anyway, I choose 0-indexing.
 * By default the overflow and underflow bins are IGNORED. But you can choose to include them
 *  by using UseOverUnderFlow( bool use );
 * Finally, just to make matters more complicated, however, when UserBinInfo classes are created,
 *  the overflow and underflow bins get their own instances.
 *
 * Binning scheme for 2D histogram: 
 *   ordered pair { (x,y) | x=[0,nbinsX-1], y[0,nbinsY-1] } 
 *   linear index { x | x[0, (nbinsX)*(nbinsY)-1 }
 * nbinsX and nbinsY includes overflow and underflow if fUseOverUnderFlow flag is set (default is NOT).
 * 
 * defined conversion from 1D to 2D via function GetBinXYFromIndex( index, &binx, &biny ) 
 *
 * --------------------------------------------------------------------------------------------
 */


#ifndef __Hist__
#define __Hist__

#include <map>
#include "UserBinInfoList.hh"
#include "UserBinInfo.hh"

#include "TH1D.h"
class TH2D;
class TFile;

namespace qosc {

  class Hist  {

  public:
    Hist( std::string name );
    virtual ~Hist();
  
    typedef enum { kReWeight, kNoReweight } UpdateMode; ///< kReWeight=Use UserBinInfo to rewight hist, kNoReweight=Ignore UserBinInfo

    // --------------------------------------------------------------------------------
    // Update Histogram using the UserBinInfo

    virtual void UpdatePDF(); ///< Feel free to overload the default method here, Update w/o reweight by default
    virtual void UpdatePDF( UpdateMode mode );

    // --------------------------------------------------------------------------------
    // Get/Set Hist information
    double GetBinContent( int binx, int biny=-1 );
    void SetBinContent( int binx, double value );
    void SetBinContent( int binx, int biny, double value );
    virtual double GetIntegral() { return m_integral; }; ///< Returns integral of histogram
    std::string GetName() { return m_name; };

    // --------------------------------------------------------------------------------
    // Misc
    void Print(); ///< summarizes the contents of the class
    void SetVerbose( int verbosity ) { m_verbose = verbosity; };
    int GetVerbose() { return m_verbose; };
    bool GetBuildStatus() { return fDoBuild; };
    void SetBuildStatus( bool dobuild ) { fDoBuild = dobuild; };

    // --------------------------------------------------------------------------------
    // Configuring the Histogram
    // Is there a better way than to double the number of functions?
    // Problem is that I'll have to double the functions in other parts of the program.

    // 1 Dimensional Functions
    void SetHistogram( TH1D* hist ); // Give the histogram to the PDF (this class thinks it owns it)
    void SetHistogram( std::string histname, int nbins, double xmin, double xmax );
    void SetHistogram( std::string histname, int nbins, double* binedges );

    // 2 Dimensional Functions
    void SetHistogram( std::string histname, int nbinsX, double xmin, double xmax, int nbinsY, double ymin, double ymax  );
    void SetHistogram( std::string histname, int nbinsX, double* binedgesX, int nbinsY, double* binedgesY );
    void SetHistogram( TH2D* hist ) { m_hist_pdf_2D = hist; fHistogramStatus = kDefined; fNDims=2; }; // Give the histogram to the PDF (this class thinks it owns it)
    void CopyHistogram( std::string copyname, TH2D* hist );
    void CopyHistogram( std::string copyname, Hist* origpdf );

    // Histogram Information
    enum HistogramStatuses { kUndefined, kDefined };
    HistogramStatuses GetHistogramStatus() { return fHistogramStatus; };
    int GetHistDimensions() { return fNDims; };
    void SetHistDimensions( int ndims ) { fNDims = ndims; };

    void CopyHistogram( std::string copyname, TH1* hist, int dims );
    TH1* GetHistogram();

    void StoreHistogram(); ///< Store a copy of the bin values to memory
    TH1* GetStoredHistogram(); ///< returns pointer to stored histogram
  
    // --------------------------------------------------------------------------------  
    // User Bin Info: method to store additional information for each bin besides the number of events there

  protected: 
    UserBinInfoList masterList;
  public:
    UserBinInfoListIter MasterBinInfoListBegin() { return masterList.Begin(); };
    UserBinInfoListIter MasterBinInfoListEnd() { return masterList.End(); };
    void AddBinInfo( std::string name, UserBinInfo* bininfo ); ///< When a user bin info class is added. Copies are made for each bin! Pay mind to your copy constructor.
    UserBinInfo* GetBinInfo( std::string name, int binX, int binY=-1 );
    UserBinInfo* GetMasterListBinInfo( std::string name ) { return masterList.GetInfo( name ); };
    UserBinInfoList* GetBinInfoList( int binX, int binY=-1 );
    UserBinInfoList* GetBinInfoListInternal( int binX, int binY=-1 );
    void CopyBinInfo( Hist* source ); ///< Deep Copy
    void PrintBinInfoLists(); // dumps information on bin info list of histogram
    void WriteBinInfo(); // save user info to file: invokes UserBinfo::Write() which can be overloaded by user
    void ResetAllBinInfo(); ///< loops through bin info and calls Reset
    void ResetBinInfo( std::string name ); ///< loops through each bin and resets only the userbininfo instance that matches 'name'

    // User Histogram Info. Same Idea as BinInfo, but this object will fill everything histogram filled, not just the bin
  protected:
    UserBinInfoList m_UserHistInfo; 
    void ProcessHistInfo( double weight=1.0 ); // have hist info objects process event
  public:
    void AddHistInfo( std::string name, UserBinInfo* info );
    UserBinInfo* GetHistInfo( std::string name );
    void WriteHistInfo();
    UserBinInfoListIter GetHistInfoListBegin() { return m_UserHistInfo.Begin(); };
    UserBinInfoListIter GetHistInfoListEnd() { return m_UserHistInfo.End(); };

    // --------------------------------------------------------------------------------  
    // Bin Indexing Option
  protected:
    bool fUseOverUnderFlow; // by default this this false
  public:
    void UseOverUnderFlow( bool use ) { fUseOverUnderFlow=use; };
    int GetBinIndex( int binx, int biny=-1 ); // returns unrolled index
    void GetBinXYFromIndex( int binindex, int& binx, int& biny );
    void GetNumberOfBinsXY( int& nbinsx, int& nbinsy ); ///< ignores overflow and underflow bins if fUseOverUnderFlow=false
    int GetNumberOfBins(); ///< ignores overflow and underflow bins if fUseOverUnderFlow=false

    // User Bin Information from File
    bool LoadHistogramsFromFile( TFile* file, double scalefactor=1.0 ); ///< load histogram and bin information data from file.  scalefactor used to scale the file's histograms.

    // Methods to access address of bins in a blind manner to details of histogram
  public:
    double* GetBinAddress( int binx, int biny=-1 ); // relies on binning convention

    // Important Settings
  protected:
    HistogramStatuses fHistogramStatus; 
    bool fDoBuild; ///< If false, will ignore commands to build the pdf

    int fNDims;
    TH1D* m_hist_pdf;
    TH2D* m_hist_pdf_2D;

    TH1D* m_hist_pdf_stored;
    TH2D* m_hist_pdf_2D_stored;
    void DestroyHistogram();

    double m_integral;
    void SetIntegral( double integral ) { m_integral = integral; };
    void SetHistogramStatus( HistogramStatuses status ) { fHistogramStatus = status; };
    int GetTotalNumberOfBins(); ///< note: includes under- and overflow by default
    void GetTotalNumberOfBinsXY(int& nbinsx, int& nbinsy ); ///< note: includes under- and overflow by default

    // Bin Info Functions
    void ClearBinInfo();
    void SetupUserBinInfo();
    void ProcessBinInfo( int binfilled, double weight=1.0 ); //< 1D
    void ProcessBinInfo( int binfilledX, int binfilledY, double weight=1.0 ); //< 2D

    std::map< int, UserBinInfoList* > m_bininfo; ///< bin info lists. one for each bin in the histogram.

  protected:
    int m_verbose;
    std::string m_name;
  };
}

#endif
