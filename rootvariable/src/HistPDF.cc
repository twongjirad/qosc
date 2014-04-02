#include "HistPDF.hh"

using namespace qosc;

HistPDF::HistPDF( std::string pdfname ) 
  : Hist( pdfname )
{
}

HistPDF::HistPDF( std::string pdfname, TH1D* hist ) 
  : Hist( pdfname )
{

  Hist::SetHistogram( hist );

}

HistPDF::HistPDF( std::string pdfname, TH2D* hist ) 
  : Hist( pdfname )
{
  
  Hist::SetHistogram( hist );
  
}


HistPDF::~HistPDF(){}

double HistPDF::GetProbability( double value ) {
  double integral = 1.0;
  if ( Hist::GetHistogram()->Integral()!=0.0 ) integral = Hist::GetHistogram()->Integral();
  
  return Hist::GetHistogram()->GetBinContent( Hist::GetHistogram()->FindBin( value ) )/integral;
  
}
