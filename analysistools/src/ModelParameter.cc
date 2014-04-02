//-*- mode:c++; c-basic-offset:2;   -*-
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

#include "ModelParameter.hh"
#include <assert.h>

using namespace qosc;

ModelParameter::ModelParameter( ) {
  m_active_flag = true;
  SetID( -1 );
  SetBounds( -1e6, 1e6 );
  m_length_scale = 0.01;
  fBounded = false;
  m_ScalingFactor = 1.0;
  has_default = false;
  fDoesTermPull = true;
  fDoesTermThrow = true;
  SetValue( 0.0 );
  SetCentralValue( 0.0 );
  VaryParameter();
}

ModelParameter::ModelParameter( std::string name, double sig, FitterParType t) {
  m_parameterName = name;
  m_active_flag = true;
  SetID( -1 );
  SetFitterParType( t );
  m_ScalingFactor = 1.0;
  m_length_scale = 0.01;
  SetBounds( -1e6, 1e6 );
  fBounded = false;
  has_default = false;
  fDoesTermPull = true;
  fDoesTermThrow = true;
  VaryParameter();
  SetValue( 0.0 );
  SetSigma(sig);
  SetCentralValue( 0.0 );
}

ModelParameter::~ModelParameter() {
}

// std::string ModelParameter::IsA() {
//   std::cout << "Called ModelParameter[" << this << "]::IsA()" << std::endl;
//   assert(false);
// }

void ModelParameter::SetValue( double value ) {
  if ( fBounded ) {
    if ( value>GetHighBound() ) {
      m_parameter_value = GetHighBound();
      std::cout << "par[" << m_parameterName << "] value had to be clamped by upper bound: " << m_parameter_value << std::endl;
    }
    else if ( value<GetLowBound() ) {
      m_parameter_value = GetLowBound();
      std::cout << "par[" << m_parameterName << "] value had to be clamped by lower bound" << m_parameter_value << std::endl;
    }
    else m_parameter_value = value;
  }
  else 
    m_parameter_value = value;
}

void ModelParameter::SetBounds( double low, double high, double length_scale ) { 
  fBounded = true;
  if ( low<high ) {
    lowBound = low;
    highBound = high;
  }
  else {
    lowBound = high;
    highBound = low;
  }
  if ( length_scale<0 ) m_length_scale = 0.01*fabs( highBound-lowBound );
  else m_length_scale = fabs(length_scale);
}

void ModelParameter::SetDefault(double v) {
  m_default_val = v; 
  has_default = true; 
}


double ModelParameter::GetDefault() { 
  if (has_default)
    SetValue(m_default_val);
  return GetValue();
}
