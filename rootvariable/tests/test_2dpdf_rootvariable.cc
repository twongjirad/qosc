#include <iostream>
#include <string>
#include <cmath>

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TRandom3.h"

#include "HistCoordinator.hh"
#include "HistRootVariable.hh"

#include "RootVariableList.hh"

using namespace qosc;

double FillVar1( RootVariableList& vars ) {
  return vars.GetVariableValue( "var1" );
}

double FillVar2( RootVariableList& vars ) {
  return vars.GetVariableValue( "var2" );
}

double FillVar1FromArray( RootVariableList& vars ) {
  return vars.GetVariableValue( "array1", 1, 0 );
}

double FillVar2FromArray( RootVariableList& vars ) {
  return vars.GetVariableValue( "array1", 1, 1 );
}

double SelectPositive( RootVariableList& vars ) {
  if ( vars.GetVariableValue("cutvar1") > 0 )
    return 1;
  else
    return 0;
}

double SelectEven( RootVariableList& vars ) {
  if ( int(vars.GetVariableValue("cutvar2")+0.2)%2==0 )
    return 1.0;
  else
    return 0.;
}

int main( int argc, char** argv ) {


  // Function tests if the output of PDFRootVariables are correct.
  // Does this by filling a tree with gaussian values

  // -------------------------------------------------------------------------------------
  // First create test variable tree
  TFile* out = new TFile("output_test_2d_pdfrootvariable.root", "RECREATE" );

  // define tree and simple variables
  TTree* test_tree = new TTree( "test_tree", "Test Tree" );
  double var1, var2, array1[2];
  double cutvar1, cutvar2;
  test_tree->Branch( "var1", &var1, "var1/D" );
  test_tree->Branch( "var2", &var2, "var2/D" );
  test_tree->Branch( "array1", array1, "array1[2]/D" );
  test_tree->Branch( "cutvar1", &cutvar1, "cutvar1/D" );
  test_tree->Branch( "cutvar2", &cutvar2, "cutvar2/D" );

  // fill variables
  TRandom3 rand(1);
  const int num = 10000;
  for (int n=0; n<num; n++) {

    var1 = rand.Gaus( 0, 2.0 );
    var2 = rand.Gaus( 0, 5.0 );

    array1[0] = var1;
    array1[1] = fabs(var1);

    cutvar1 = rand.Gaus(0,1.0);
    cutvar2 = n;
    
    test_tree->Fill();
  }

  // write tree
  out->Write();
  out->Close();
  // -------------------------------------------------------------------------------------
  
  // Open tree and declare PDF objects
  out = new TFile("output_test_hists.root", "RECREATE" );

  // Now load Tree
  TChain* test = new TChain("test_tree");
  test->Add( "output_test_2d_pdfrootvariable.root" );
  std::cout << "Number of Entries in Tree: " << test->GetEntries() << std::endl;

  // Create reference histograms. Using ROOTs built in Draw functions
  TH2D* hall_ref = new TH2D("hall_ref", "All Entries", 100, -100, 100, 100, -100, 100 );
  test->Draw("var1:var2>>hall_ref");
  TH2D* heven_ref = new TH2D("heven_ref", "Even Throws", 100, -100, 100, 100, -100, 100 );
  test->Draw("var1:var2>>heven_ref", "cutvar2%2==0");
  TH2D* hboth_ref = new TH2D("hboth_ref", "Positive and Even Throws", 100, -100, 100, 100, -100, 100 );
  test->Draw("var1:var2>>hboth_ref", "cutvar1>0 && cutvar2%2==0" ); // positive and even
  
  std::cout << "Made Reference Histograms" << std::endl;

  // Create a PDFMaker for the tree
  HistCoordinator* pdfmaker = new HistCoordinator( test );
  
  // Now create a bunch of test pdfrootvariables
  HistRootVariable* hall_pdf = new HistRootVariable( "hall_pdf",
						   "var1","var2",
						   "", "", test );
  hall_pdf->SetHistogram( "hall_pdf", 100, -100, 100, 100, -100, 100 );
  pdfmaker->Add( hall_pdf );

  HistRootVariable* heven_pdf = new HistRootVariable( "heven_pdf",
						    "var1","var2",
						    "", "cutvar2%2==0", test );
  heven_pdf->SetHistogram( "heven_pdf", 100, -100, 100, 100, -100, 100 );
  pdfmaker->Add( heven_pdf );

  HistRootVariable* hboth_pdf = new HistRootVariable( "hboth_pdf",
						    &FillVar1, "var1", &FillVar2, "var2",
						    NULL,
						    &SelectPositive, "cutvar1",
						    test );
  hboth_pdf->AddSelectionFunction( &SelectEven, "cutvar2" );
  hboth_pdf->SetHistogram( "hboth_pdf", 100, -100, 100, 100, -100, 100 );
  pdfmaker->Add( hboth_pdf );

  pdfmaker->BuildPDFs();

  TH2D* hall_diff = new TH2D("hall_diff", "Diff in hall", 100, -100, 100, 100, -100, 100 );
  hall_diff->Add( hall_ref );
  hall_diff->Add( (TH2D*)hall_pdf->GetHistogram(), -1.0 );
  std::cout << "Integral of hall_ref: " << hall_ref->Integral() << std::endl;
  std::cout << "Integral of hall_pdf: " << ((TH2D*)hall_pdf->GetHistogram())->Integral() << std::endl;
  std::cout << "Integral of hall_diff: " << hall_diff->Integral() << std::endl;

  TH2D* heven_diff = new TH2D("heven_diff", "Diff in heven", 100, -100, 100, 100, -100, 100 );
  heven_diff->Add( heven_ref );
  heven_diff->Add( (TH2D*)heven_pdf->GetHistogram(), -1.0 );
  std::cout << "Integral of heven_ref: " << heven_ref->Integral() << std::endl;
  std::cout << "Integral of heven_pdf: " << ((TH2D*)heven_pdf->GetHistogram())->Integral() << std::endl;
  std::cout << "Integral of heven_diff: " << heven_diff->Integral() << std::endl;
  
  TH2D* hboth_diff = new TH2D("hboth_diff", "Diff in hboth", 100, -100, 100, 100, -100, 100 );
  hboth_diff->Add( hboth_ref );
  hboth_diff->Add( (TH2D*)hboth_pdf->GetHistogram(), -1.0 );
  std::cout << "Integral of hboth_ref: " << hboth_ref->Integral() << std::endl;
  std::cout << "Integral of hboth_pdf: " << ((TH2D*)hboth_pdf->GetHistogram())->Integral() << std::endl;
  std::cout << "Integral of hboth_diff: " << hboth_diff->Integral() << std::endl;

  out->Write();
  
}
