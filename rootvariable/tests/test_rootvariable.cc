#include "TTree.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TChain.h"
#include <iostream>
#include <string>
#include <cmath>

#include "RootVariable.hh"
#include "RootVariableFormula.hh"
#include "RootVariableFunction.hh"
#include "RootVariableList.hh"
#include "RootVariableManager.hh"

using namespace qosc;

double test_function( RootVariableList& vars ) {
  return vars.GetVariableValue( "y[1][2]" )*vars.GetVariableValue("x");
}

int main( int iarg, char** argv ) {

  int nsamples = 100;

  // Create a tree
  std::cout << "Creating tree" << std::endl;
  TFile* out = new TFile( "output_test_rootfile.root", "RECREATE" );

  TTree* tree = new TTree( "test", "Test Tree" );
  double x;
  double y[2][3];
  char s[50];
  tree->Branch( "x", &x, "x/D" );
  tree->Branch( "y", y, "y[2][3]/D" );
  tree->Branch( "s", s, "s[50]/C" );

  std::string muppets[7] = { "Kermit", "Ms. Piggy", "Animal", "Fozzy", "Gonzo", "BEAKER!", "Sweedish Chef" };

  TRandom3 rand(10);

  for (int i=0; i<nsamples; i++) {
    x = rand.Gaus(0,1);
    for (int a=0; a<2; a++) {
      for (int b=0; b<3; b++) {
	y[a][b] = rand.Gaus(1,2);
      }
    }
    sprintf( s, "%s", muppets[ int( rand.Uniform(0,7) ) ].c_str() );
    tree->Fill();
  }

  out->Write();

  out->Close();

  // Read in Tree two ways: RootVariable and Standard
  std::cout << "Verifying Output" << std::endl;
  
  TChain* standard = new TChain( "test" );
  standard->Add( "output_test_rootfile.root" );
  standard->SetBranchAddress( "x", &x );
  standard->SetBranchAddress( "y", y );
  standard->SetBranchAddress( "s", s );

  TChain* tc_rootvar = new TChain( "test" );
  tc_rootvar->Add( "output_test_rootfile.root" );
  RootVariableList* rootvars = new RootVariableList( "x:s:y[1][1]:x*y[0][2]/y[1][0]", tc_rootvar );
  RootVariable* rootvar_y =  new RootVariableFormula( "y", tc_rootvar );
  rootvars->Add( "myfunc", &test_function, "x:y[1][2]", tc_rootvar );

  RootVariableManager::GetTheRootVarManager()->DumpTheVariableCache();

  int xwrong = 0;
  int swrong = 0;
  int ywrong = 0;
  for (int i=0; i<nsamples; i++) {
    standard->GetEntry(i);
    tc_rootvar->GetEntry(i);

    double rootvar_x = rootvars->GetVariableValue( "x" );
    std::string rootvar_s = rootvars->GetVariable( "s" )->PrintValue();
    if ( rootvar_x!=x ) {
      std::cout << "Entry " << i << ": value of x disagrees. standard=" << x << " rootvar=" << rootvar_x <<  std::endl;
      xwrong++;
    }
    if ( rootvar_s!=std::string(s) ) {
      std::cout << "Entry " << i << ": value of s disagrees. standard=" << s << " rootvar=" << rootvar_s <<  std::endl;
      swrong++;
    }
    bool yok = true;
    int a=1;
    int b=1;
    if ( rootvars->GetVariableValue( "y[1][1]")!=y[a][b] ) {
      std::cout << "Entry " << i << ": value of y[" << a << "][" << b << "] disagrees. standard=" << y[a][b] << " rootvar=" << rootvars->GetVariableValue( "y[1][1]" ) <<  std::endl;
      yok = false;
    }
    // ==============================
    // Array axis in list not working
//     for (int a=0; a<2; a++) {
//       for (int b=0; b<3; b++) {
// 	if ( rootvars->GetVariableValue( "y", 2, a, b )!=y[a][b] ) {
// 	  std::cout << "Entry " << i << ": value of y[" << a << "][" << b << "] disagrees. standard=" << y[a][b] << " rootvar=" << rootvars->GetVariableValue( "y", 2, a, b ) <<  std::endl;
// 	  yok = false;
// 	}
//       }
//     }
    if ( yok==false )
      ywrong++;


    // ==============================
    // Array axis by variable
    a = 1;
    b = 2;
    double array_val = rootvar_y->Value( 2, a, b );
    if ( array_val!=y[a][b] ) {
      std::cout << "Entry " << i << ": value of y[" << a << "][" << b << "] disagrees. standard=" << y[a][b] << " rootvar=" << array_val << std::endl;
    }
    
    // Formula test: Formula here is good to 1e-14
    if ( std::fabs(rootvars->GetVariableValue( "x*y[0][2]/y[1][0]")-x*y[0][2]/y[1][0])>1e-14 ) {
      std::cout << "Entry " << i << ": value of formula disagrees. standard=" << x*y[0][2]/y[1][0] << " rootvar=" << rootvars->GetVariableValue( "x*y[0][2]/y[1][0]" ) <<  std::endl;
    }

    // ================
    // Function test
    if ( fabs(rootvars->GetVariableValue( "myfunc" )-x*y[1][2])>1e-14 ) {
      std::cout << "Entry " << i << ": value of formula disagrees. standard=" << x*y[1][2] << " rootvar=" << rootvars->GetVariableValue( "myfunc" ) <<  std::endl;
    }

    if ( i%10==0 ) {
      std::cout << "Sample Entry " << i << ": " << std::endl;
      std::cout << "  value of x:  standard=" << x << " rootvar=" << rootvar_x <<  std::endl;
      std::cout << "  value of s:  standard=" << s << " rootvar=" << rootvar_s <<  std::endl;
    }
  }

  std::cout << "Number of x mistakes: " << xwrong << std::endl;
  std::cout << "Number of s mistakes: " << swrong << std::endl;
  std::cout << "Number of y mistakes: " << ywrong << std::endl;

  return 0;
}
