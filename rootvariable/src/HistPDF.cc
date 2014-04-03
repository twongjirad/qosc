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

#include "HistPDF.hh"
#include "TH1D.h"

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
