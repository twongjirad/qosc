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

#include "OscTerm.hh"
#include "WeightBargerTruthHist.hh"

using namespace qosc;

OscTerm::OscTerm( std::string parname, double initial, double low, double high )
  : BasicParameter( parname, initial, 1.0, kMinuitTerm )
{
  SetBounds( low, high, 0.005*(high-low) );
  SetValue( initial );
  TermDoesNotPull();
  TermDoesNotThrow();  
}

OscTerm::~OscTerm() {

}

void OscTerm::SetValue( double value ) {
  ModelParameter::SetValue( value ); // this needs to be hidden from the user.
  ModelParameter::SetCentralValue( value ); // this needs to be hidden from the user.

  // I hate this so very much:
  //WeightBargerTruthHist::GetGlobalInstance()->SetParameter( GetName(), value ); // interface to the oscillation calculator, Sample histograms will reweight using this singleton.
  //std::cout << "OscTerm[" << GetName() << "]::SetValue( " << value << " )" << std::endl;
}
