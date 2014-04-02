#include <iostream>
#include <string>

#include "ProfileMinuit.hh"
#include "TestModel.hh"
#include "PMITestModel.hh"
#include "TFile.h"

using namespace qosc;

int main( int iarg, char** argv ) {

  /*
   *   Verbose levels:
   *    0: Quiet (at least base class will be)
   *    1: Print out summary of fit
   *    2: Print out more information about each step
   *    3: Print out maximum amount of info
   *    4: Pause after each interation
 */
  TFile* out = new TFile( "output_testmodel.root", "RECREATE" );


  //bool allMinuit = true;
  bool allMinuit = false;
  TestModel* model = new TestModel( 20, 10, 3, 20, 10, 0.1, allMinuit );


  ProfileMinuit* fitter = new ProfileMinuit();
  fitter->SetVerbose(1);
  fitter->SetNewtonSolverVerbose( 0 );
  fitter->SetMinuitVerbose( 0 );

  PMITestModel* interface = new PMITestModel( fitter, model );
  interface->InitializeFitter();

  //model->SetBackgroundMax( 1.2 );
  //model->SetSignalMax( 1.1 );
  model->SetSignalMean( 0.2 );
  model->SetSignalSigma( -0.2 );
  model->SetBackgroundConstant( 0.05 );

  model->m_parameters->Print();
  model->Print();
  std::cin.get();

  fitter->DoFit();

  out->Write();

  return 0;

};
