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

#include <iostream>
#include <string>
#include <cmath>

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TRandom3.h"

#include "HistCoordinator.hh"
#include "HistRootVariable.hh"

#include "RootVariableList.hh"

using namespace qosc;

double FillVar1( RootVariableList& vars ) {
  return vars.GetVariableValue( "var1" );
}

double FillVar1FromArray( RootVariableList& vars ) {
  return vars.GetVariableValue( "array1", 1, 0 );
}

double SelectPositive( RootVariableList& vars ) {
  if ( vars.GetVariableValue("cutvar1") > 0 )
    return 1;
  else
    return 0;
}

double SelectEven( RootVariableList& vars ) {
  if ( int(vars.GetVariableValue("cutvar2")+0.4)%2==0 )
    return 1.0;
  else
    return 0.;
}

int main( int argc, char** argv ) {


  // Function tests if the output of PDFRootVariables are correct.
  // Does this by filling a tree with gaussian values


  // First create test variable tree
  TFile* out = new TFile("output_test_pdfrootvariable.root", "RECREATE" );

  TTree* test_tree = new TTree( "test_tree", "Test Tree" );
  double var1, var2, array1[2];
  double cutvar1, cutvar2;
  test_tree->Branch( "var1", &var1, "var1/D" );
  test_tree->Branch( "var2", &var2, "var2/D" );
  test_tree->Branch( "array1", array1, "array1[2]/D" );
  test_tree->Branch( "cutvar1", &cutvar1, "cutvar1/D" );
  test_tree->Branch( "cutvar2", &cutvar2, "cutvar2/D" );

  TRandom3 rand(1);
  const int num = 10000;
  for (int n=0; n<num; n++) {

    var1 = rand.Gaus( 0, 20.0 );
    var2 = n;

    array1[0] = var1;
    array1[1] = fabs(var1);

    cutvar1 = var1;
    cutvar2 = n;
    
    test_tree->Fill();
  }

  out->Write();
  out->Close();
  
  out = new TFile("output_test_hists.root", "RECREATE" );

  // Now load Tree
  TChain* test = new TChain("test_tree");
  test->Add( "output_test_pdfrootvariable.root" );

  // Create reference histograms
  TH1D* hpos_ref = new TH1D("hpos_ref", "Positive", 100, -100, 100 );
  test->Draw("var1>>hpos_ref", "cutvar1>0");
  TH1D* heven_ref = new TH1D("heven_ref", "Even Throws", 100, -100, 100 );
  test->Draw("var1>>heven_ref", "cutvar2%2==0");
  TH1D* hboth_ref = new TH1D("hboth_ref", "Positive and Even Throws", 100, -100, 100 );
  test->Draw("var1>>hboth_ref", "cutvar1>0 && cutvar2%2==0" );
  
  // Create a PDFMaker for the tree
  HistCoordinator* pdfmaker = new HistCoordinator( test );
  
  // Now create a bunch of test pdfrootvariables
  
  // one based on variable function, weight, selection function
//   PDFRootVariable* hpos_pdf = new PDFRootVariable( "hpos_pdf",
// 						   "var1",
// 						   NULL,
// 						   "cutvar1>0",
// 						   test );


//   PDFRootVariable* hpos_pdf = new PDFRootVariable( "hpos_pdf",
// 						   "var1",
// 						   NULL,
// 						   "cutvar1>0",
// 						   test );

  HistRootVariable* hboth_pdf = new HistRootVariable( "hboth_pdf",
						     &FillVar1, "var1",
						     NULL,
						     &SelectPositive, "cutvar1",
						     test );
  hboth_pdf->AddSelectionFunction( &SelectEven, "cutvar2" );
  hboth_pdf->SetHistogram( "hboth_pdf", 100, -100, 100 );
  pdfmaker->Add( hboth_pdf );

  pdfmaker->BuildPDFs();
  
  TH1D* hboth_diff = new TH1D("hboth_diff", "Diff in hboth", 100, -100, 100 );
  hboth_diff->Add( hboth_ref );
  hboth_diff->Add( (TH1D*)hboth_pdf->GetHistogram(), -1.0 );

  std::cout << "Integral of hboth_ref: " << hboth_ref->Integral() << std::endl;
  std::cout << "Integral of hboth_pdf: " << ((TH1D*)hboth_pdf->GetHistogram())->Integral() << std::endl;
  std::cout << "Integral of hboth_diff: " << hboth_diff->Integral() << std::endl;
  out->Write();
  
}
