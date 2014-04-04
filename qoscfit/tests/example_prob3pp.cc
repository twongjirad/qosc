/*
  Tests of oscillation framework code and interface

 */

#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TH1D.h"
#include "TTree.h"
#include "TChain.h"

#include "BargerPropagator.h"

#include "HistRootVariable.hh"
#include "HistCoordinator.hh"
#include "WeightBargerOsc.hh"

using namespace qosc;

int main( int nargs, char** argv ) {
  
  // Make fake tree with neutrino info
  TFile* out = new TFile( "out.root", "RECREATE" );
  TTree* test = new TTree( "nutest", "Test tree" );
  int iflux, ixsec; // 1,2,3=e,mu,tau. <0 is anti-neutrino
  double Enu; // True energy
  double L=295; // km
  double weight=1.0; // other weight. Could be flux for example.
  double oscprob = 1.0;
  int mode = 1; 
  test->Branch( "Enu_GeV", &Enu, "Enu_GeV/D" );
  test->Branch( "L_km", &L, "L_km/D" );
  test->Branch( "iflux", &iflux, "iflux/I" );
  test->Branch( "ixsec", &ixsec, "ixsec/I" );
  test->Branch( "imode", &mode, "imode/I" );
  test->Branch( "weight", &weight, "weight/D" );
  test->Branch( "oscprob", &oscprob, "oscprob/D" );

  double meanEnu = 0.6; // GeV
  double sigmaEnu = 0.2; // GeV

  TRandom3 rand(1);
  int nevents = 1000;

  // in example set to numu-disappearance
  iflux = 2;
  ixsec = 2;

  // set osc pars
  double s12 = 0.312;
  double s13 = 0.0251;
  double s23 = 0.5;
  double dm12 = 7.5e-5;
  double dm32 = 2.5e-3;
  double cp = 0.0;
  double density = 2.6; // g/cm^2: matter density

  // Set up Prob3++
  // T: sin^2(x) F: sin^2(2x) 
  bool kSquared = true;
  BargerPropagator* bp = new BargerPropagator( false );
  
  for (int i=0; i<nevents; i++) {
    Enu = rand.Gaus(meanEnu, sigmaEnu);
    // setup oscillator
    bp->SetMNS( s12, s13, s23, dm12, dm32, cp, Enu, kSquared, iflux );
    bp->propagateLinear( iflux, L, density );
    oscprob = bp->GetProb( iflux, ixsec );
    test->Fill();
  }
  
  std::cout << "Number of events generated: " << test->GetEntries() << std::endl;

  test->Write();
  out->Close();

  // Load tree twice: once using machinery, other using more straight-forward techniques
  TChain* tc_test1 = new TChain( "nutest" );
  tc_test1->Add( "out.root" );

  TChain* tc_test2 = new TChain( "nutest" );
  tc_test2->Add( "out.root" );

  // Attach bracnhes to test1
  tc_test1->SetBranchAddress( "Enu_GeV", &Enu );
  tc_test1->SetBranchAddress( "L_km", &L );
  tc_test1->SetBranchAddress( "iflux", &iflux );
  tc_test1->SetBranchAddress( "ixsec", &ixsec );
  tc_test1->SetBranchAddress( "imode", &mode );
  tc_test1->SetBranchAddress( "weight", &weight );
  tc_test1->SetBranchAddress( "oscprob", &oscprob );

  out = new TFile( "out.root", "UPDATE" );

  // Load weight class: too much configured to T2K. ex: for SK would need L depedence.
  WeightBargerOsc* bpweight = new WeightBargerOsc( kSinTheta, tc_test2, "Enu_GeV", "iflux", "ixsec", "imode", "weight" );
  bpweight->SetMNS( s12, s13, s23, dm12, dm32, cp );

  // Now define a histogram
  TH1D* henu = new TH1D("hEnubins","",100,0,4);
  HistRootVariable* enu = new HistRootVariable( "histEnu", "Enu_GeV", tc_test2 );
  enu->LoadWeightClass( "osc", bpweight );
  enu->SetHistogram( henu );

  // We use the HistCoordinator class to hill the histogram
  HistCoordinator* coord = new HistCoordinator( tc_test2 );
  coord->Add( enu );

  // Output Setup and then Fill Histograms
  coord->Print();
  coord->BuildPDFs();

  // Build Comparison Histogram
  TH1D* henu_compare = (TH1D*)henu->Clone( "hEnubins_compare" );
  henu_compare->Reset();
  
  // Use ROOT technique
  tc_test1->Draw( "Enu_GeV>>hEnubins_compare", "weight*oscprob*(1==1)" );

  // Compare: Note offset. ROOT hists are by default indexed to 1. HistRootVariable is indexed to 0 by default (can be changed to include under and overflow using switch)
  int nbinsmatched = 0;
  for (int i=0; i<henu_compare->GetNbinsX(); i++) {
    if ( fabs( henu_compare->GetBinContent(i+1)-enu->GetBinContent(i) )>=1e-15 ) {
      std::cout << henu_compare->GetName() << " bin " << i << " disagrees: "
		<< " ROOT Draw method=" << henu_compare->GetBinContent(i+1) 
		<< " HistRootVariable method=" << enu->GetBinContent(i) << std::endl;
    }
    else
      nbinsmatched++;
  }

  henu_compare->Write();
  enu->GetHistogram()->Write();
  
  std::cout << nbinsmatched << " out of " << henu_compare->GetNbinsX() << " agree." << std::endl;
  std::cout << "FIN." << std::endl;
}
