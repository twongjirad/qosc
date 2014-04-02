#include "SampleManagerROOTPDF.hh"
#include <iostream>
#include <assert.h>
#include "TChain.h"
#include "SampleROOTPDF.hh"

#include "ParameterManager.hh"
#include "RootVariableManager.hh"
#include "Weight.hh"
#include "EventReweightParameter.hh"
#include "SysTermI.hh"
//#include "ODHistInfo.hh"

using namespace qosc;

SampleManagerROOTPDF::SampleManagerROOTPDF()
  : SampleManager() 
{
  fBuildHistogramsFromChains = true;
  fAlwaysRebuildHistogramsFromChains = false;
  loadingFromFile = false;
}

SampleManagerROOTPDF::~SampleManagerROOTPDF() {
}

void SampleManagerROOTPDF::AddRootChain( std::string name, TChain* chain ) {
  if ( IsChainDefined( name ) ){
    std::cout << "SampleManagerROOTPDF::AddRootChain -- WARNING! Redefining chain with name: " << name << ". Potential Memory Leak." << std::endl;
  }
  m_chain_dict[name] = chain;
}

void SampleManagerROOTPDF::AddRootChain( std::string name, void* chain_ptr ) {
  TChain* chain = (TChain*)chain_ptr;
  AddRootChain( name, chain );
}

void SampleManagerROOTPDF::FillSampleBins( ParameterManager* parameters, bool fillnominal ) {

//   // ------------------------------------------------------------------------------------------------------
//   // First we check if we are to loading the samples from a TFile instance, and if so, if we've already done it.

//   if ( loadingFromFile ) {
//     //LoadSamplesFromFile( parameters ); //[obsolete?]
//   }

  // ------------------------------------------------------------------------------------------------------
  // Check if we need to loop through ROOT TChain entries to build the sample expectation histograms
  
  if ( fAlwaysRebuildHistogramsFromChains || fBuildHistogramsFromChains || DoSamplesHaveEventReweightParameter(parameters) ) {
    
    bool systermNeedsFutureLoop = false; // flag which will be set to true if any of the sys terms we encounter need this loop to be run again. (this is to force event-by-event reweighting)
    if ( GetVerbose()>=1 ) {
      std::cout << "----------------------------------------------------------------------" << std::endl;
      std::cout << "SampleManagerROOTPDF::FillSampleBins - Filling bins from TChains" << std::endl;
      std::cout << "----------------------------------------------------------------------" << std::endl;
    }

    // build pdfs for each chain.
    for ( ChainDictIter it = ChainDictBegin(); it!=ChainDictEnd(); it++) {
      // First get the ROOT chain we will use to fill.
      std::string chainname = (*it).first;
      TChain* chain = (*it).second;

      // We create a coordinator instance to manage all the histograms we will need to fill with this chain
      HistCoordinator* chainManager = new HistCoordinator( chain );
      
      // We bound the entries with which we use to fill the histograms. But only if the user has specified us to do so.
      if ( m_chain_start.find(chainname)!=m_chain_start.end() ) chainManager->SetStartEntry( m_chain_start[chainname] );
      if ( m_chain_end.find(chainname)!=m_chain_end.end() ) chainManager->SetMaxEntries( m_chain_end[chainname] );
    
      // Next loop through each of the sample histograms that need this chain.
      for ( SampleDictIter itsample = SampleDictBegin(); itsample!=SampleDictEnd(); itsample++ ) {
	std::string samplename = (*itsample).first;
	SampleROOTPDF* asample = dynamic_cast<SampleROOTPDF*>((*itsample).second);

	if ( !asample->DoesSampleUseChain(chainname) ) continue;
	if ( !asample->IsActive() ) continue;

	// Add the sample histogram to the manager.
	asample->GetSampleHistogram()->ResetAllBinInfo();
	chainManager->Add( (HistRootVariable*)asample->GetSampleHistogram() );

	// loop over hist info, and store pdf's
	// for ( UserBinInfoListIter it_info= asample->GetSampleHistogram()->GetHistInfoListBegin(); it_info!=asample->GetSampleHistogram()->GetHistInfoListEnd(); it_info++ ) {
	//((ODHistInfo*)(*it_info).second)->GetInfoPDF()->GetHistogram()->Reset();
	//chainManager->Add( ((ODHistInfo*)(*it_info).second)->GetInfoPDF() ); // note, this seems fishy. should not have something specific here.
	//}
      
	// Loop through the parameter manager and handle the systematic error terms:
	// (1) For EventReweightParameter classes: add the sys term's EventWeightGenerator to the PDF class
	std::vector< std::string > parlist;
	parameters->GetListOfParameterNames( parlist );
	for ( std::vector< std::string >::iterator itsys=parlist.begin(); itsys!=parlist.end(); itsys++) {
	  
	  //if ( !asample->DoesParameterApplyToSample( *itsys ) ) continue;
	  //if ( !parameters->DoesTermApplyToSample( *itsys ) ) continue;
	  if ( !parameters->GetParameter( *itsys )->IsActive() ) continue;
	  
	  
	  if ( parameters->GetParameterTypeName( *itsys )=="EventReweightParameter" ) {

	    // This type of parameter is using an event-by-event reweight method. 
	    // We need to give the sample histogram the class that generates the event-by-event weights. We only need to do this once.
	    // Maybe this should be assumed to be done correctly at the initialization of the sample.

	    // Get the Sample's histogram (represented as a HistRootVariable object)
	    HistRootVariable* samplePDF = dynamic_cast< HistRootVariable*>( asample->GetSampleHistogram() );
	    assert( samplePDF!=NULL );

	    // Cast the parameter to the EventReweightParameter type
	    EventReweightParameter* systerm = dynamic_cast<EventReweightParameter*>( parameters->GetParameter( *itsys ) );
	    assert( systerm!=NULL ); // This is not the right type!

	    // Get the event weight generator.
	    // This object should inherit from two different classes.
	    // Weight (from pdf lib) and SysTermI (from analysistools lib).
	    // Weight has methods needed for PDF objects to weight each event as it fills the histogram from a ROOT tree.
	    // SysTermI acts as the interface to an event weight generator outside of our package. It helps us get the weight we need.
	    // The two classes are redundant. Need to figure out a way to get rid of this structure.

	    Weight* eventweightgen = dynamic_cast< Weight* >( systerm->GetEventWeightGen( samplename ) );
	    SysTermI* systermgen = systerm->GetEventWeightGen( samplename );
	    if ( eventweightgen==NULL && systermgen==NULL ) {
	      std::cout << " No reweightgen defined for this sample. par=" << *itsys << " sample=" << samplename << std::endl;
	      continue;
	    }

	    std::cout << "Setup EventReweightParameter: sample=" << samplename << "( " << samplePDF->GetName() << ")" 
		      << " weightgen="<< eventweightgen << " systermgen=" << systermgen->GetName() << " " << systermgen << " value=" << systerm->GetValue() << std::endl;
	    assert( eventweightgen!=NULL );
	    
	    // Add it to the HistRootVariable object if it has not already been added
	    if ( !samplePDF->IsWeightClassDefined( systermgen->GetName() ) ) {
	      assert( eventweightgen!=NULL );
	      samplePDF->LoadWeightClass( systermgen->GetName(), eventweightgen );
	      std::cout << "Loading Weight Class " << systermgen->GetName() << " into sample=" << samplename << " (pdf rootvar=" << samplePDF->GetName() << ")" << std::endl;
	    }
	    
	    
	    // // associate this weight generator with any histinfo pdf's (we know about our analysis's implementation)
// 	    for ( UserBinInfoListIter it_info= asample->GetSampleHistogram()->GetHistInfoListBegin(); it_info!=asample->GetSampleHistogram()->GetHistInfoListEnd(); it_info++ ) {
// 	      HistRootVariable* fine_binning_pdf = ((ODHistInfo*)(*it_info).second)->GetInfoPDF();
// 	      if ( !fine_binning_pdf->IsWeightClassDefined( systermgen->GetName() ) ) {
// 		assert( eventweightgen!=NULL );
// 		fine_binning_pdf->LoadWeightClass( systermgen->GetName(), eventweightgen );
// 		std::cout << "Loading Weight Class " << systermgen->GetName() << " into pdfrootvariable=" << fine_binning_pdf->GetName() << std::endl;
// 	      }
// 	    }
	    
	    //[ now we have to set the weight value ]
	    systermgen->SetPullValue( systerm->GetValue(), true ); // true=reconfigure
	    
	    // Tell class we need to go through the event loop
	    systermNeedsFutureLoop = true;

	  }
	  else {
	    // do nothing
	  }

	  // OLD. USED TO WANT TO KNOW ABOUT ALL PARAMETERS. THE LESS I KNOW THE BETTER.
// 	  else if ( parameters->GetParameterTypeName( *itsys )=="UnifiedResponseCurveParameter" 
// 		    ||  parameters->GetParameterTypeName( *itsys )=="SpectrumDistortionParameter" ) {
// 	    // does nothing
// 	  }
// 	  else {
// 	    std::cout << "Unanticipated Sys Term Type: " << parameters->GetParameterTypeName( *itsys ) << std::endl;
// 	    assert(false);
// 	  }
	  
	}//end of loop over parameters that apply to the sample
      }//end of loop over sample(s) that apply to the chain
      
      // optimization. deactiveates all branches. this way we load a minimum of data from file.
      //RootVarManager::GetTheRootVarManager()->DeactivateUnusedBranches();
      //RootVarManager::GetTheRootVarManager()->DumpTheVariableCache();
      chainManager->Print();

      // tell the coordinator to build the pdfs (this loops over the chains and fills all the histograms we need.)
      chainManager->BuildPDFs();

      // forgot what this does.
      //RootVarManager::GetTheRootVarManager()->ResetManager();
    
      // Destroy the manager now that we are done with it.
      delete chainManager;


      // Now that in principle, needed data from ROOT file has been extracted, Store the original histogram values in order to use them for BinInfo reweighting.
      for ( SampleDictIter itsample = SampleDictBegin(); itsample!=SampleDictEnd(); itsample++ ) {
	std::string samplename = (*itsample).first;
        SampleROOTPDF* asample = dynamic_cast<SampleROOTPDF*>((*itsample).second);
	
	if ( !asample->DoesSampleUseChain(chainname) ) continue;
        if ( !asample->IsActive() ) continue;
	
	if ( asample->GetSampleHistogram() )
	  asample->GetSampleHistogram()->StoreHistogram();

	// loop over hist info, and store pdf's
	//for ( UserBinInfoListIter it_info= asample->GetSampleHistogram()->GetHistInfoListBegin(); it_info!=asample->GetSampleHistogram()->GetHistInfoListEnd(); it_info++ ) {
	//  ((ODHistInfo*)(*it_info).second)->GetInfoPDF()->StoreHistogram();
	//}
	
      }//end of loop over samples      
    }//end of loop over chains
    
    // If we don't need to run through this loop anymore, flag it here.
    if ( !systermNeedsFutureLoop )
      fBuildHistogramsFromChains = false; 
  }// if need to build loop from chains
  
  // ------------------------------------------------------------------------------------------------------
  // Finally, Proceed to fill the expectation bins

  SampleManager::FillSampleBins( parameters, fillnominal );

}

void SampleManagerROOTPDF::RegisterSample( std::string samplename, SampleROOTPDF* asample ) {
  std::string chainname = asample->GetChainName();
  if ( !GetChain( chainname ) ) {
    std::cout << "SampleManagerROOTPDF::RegisterSample -- ERROR. Chain with name='" << chainname << "' used to build Sample " << asample->GetName() << " was not found in the Sample Manager" << std::endl;
    exit(EXIT_FAILURE);
  }
  SampleManager::RegisterSample( samplename, asample );
}

void SampleManagerROOTPDF::Print() {
  std::cout << "====================================" << std::endl;
  std::cout << "SampleMangerROOTPDF" << std::endl;
  std::cout << "Number of Chains: " << m_chain_dict.size() << std::endl;
  int nchain = 0;
  for ( ChainDictIter it=ChainDictBegin(); it!=ChainDictEnd(); it++) {
    std::cout << "(" << nchain+1 << ") "<< (*it).first << " " << (*it).second << std::endl;
  }
  SampleManager::Print();
}

bool SampleManagerROOTPDF::LoadSamplesFromFile( TFile* file, double scalefactor ) {
  // In our special implementation. The sample manager not only builds the sample histograms, but also the Fij histograms for the systematic error terms.
  // so, to do this, we need a few extra switches

  // load the sample histograms
  bool loadedok = true;
  for ( SampleDictIter it=SampleDictBegin(); it!=SampleDictEnd(); it++ ) {
    (*it).second->SetVerbose( GetVerbose() );
    dynamic_cast< SampleROOTPDF* >( (*it).second )->LoadSampleFromFile( file, scalefactor );
  }
  
  // determine if we need to ask the samples to be remade
  if ( loadedok ) {
    fBuildHistogramsFromChains = false;
  }
  else 
    fBuildHistogramsFromChains = true;

  loadingFromFile = true;
  m_sourcefile = file;
  return true;
}

void SampleManagerROOTPDF::SetChainBounds( std::string chainname, int start, int end ) {
  m_chain_start[chainname] = start;
  m_chain_end[chainname] = end;
}

bool SampleManagerROOTPDF::DoSamplesHaveEventReweightParameter( ParameterManager* parameters ) {
  for ( ParameterManager::ParListIter it=parameters->ParListBegin(); it!=parameters->ParListEnd(); it++ )
    if ( parameters->GetParameterTypeName( *it )=="EventReweightParameter" ) return true;
  return false;
}
