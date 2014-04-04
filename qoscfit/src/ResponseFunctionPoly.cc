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

#include "ResponseFunctionPoly.hh"
#include <cmath>
#include <cstdlib>
#include <assert.h>

#include "TSpline.h"

#include "ResponseFunctionSpline.hh"

using namespace qosc;

ResponseFunctionPoly::ResponseFunctionPoly( double pars[], int nterms ) {
  m_nterms = nterms;
  params = NULL;
  if ( m_nterms>0 ) {
    params = new double[m_nterms];
    for (int i=0; i<m_nterms; i++) params[i] = pars[i];
    params[0] = 1.0; // hacky fix. but definitionally true.
  }
  m_weight = 0;
}

ResponseFunctionPoly::~ResponseFunctionPoly() {
  if ( m_nterms>0 ) {
    delete [] params;
  }
  m_nterms = 0;
}

double ResponseFunctionPoly::GetResponse( double x[] ) {
  double f = 0;
  for (int i=0; i<m_nterms; i++)
    f += params[i]*std::pow( x[0], i );
  return f;
}

double ResponseFunctionPoly::GetDerivative( double x[] ) {
  double f = 0;
  for (int i=1; i<m_nterms; i++)
    f += double(i)*params[i]*std::pow( x[0], i-1 );
  return f;
}

void ResponseFunctionPoly::Add( ResponseFunction* addition, double weight ) {
  if ( addition->IsA()=="Poly" ) {
    ResponseFunctionPoly* polyfunc = dynamic_cast<ResponseFunctionPoly*>(addition);
    if ( !polyfunc ) assert( false );

    double totalweight = weight+m_weight;
    if ( totalweight>0.) {
      for (int nterm=0; nterm<m_nterms; nterm++) {
      params[nterm] = (params[nterm]*m_weight+polyfunc->params[nterm]*weight)/totalweight;
      }
      m_weight = totalweight;
    }
  }
  else if ( addition->IsA()=="Spline" ) {
    // linearize for now
    ResponseFunctionSpline* splinefunc = dynamic_cast<ResponseFunctionSpline*>(addition);
    if ( !splinefunc ) assert( false );
    
    double totalweight = weight+m_weight;
    if ( totalweight>0.) {
      double inputx[1] = { 0 };
      params[0] = (params[0]*m_weight+splinefunc->GetResponse(inputx)*weight)/totalweight;
      params[1] = (params[1]*m_weight+splinefunc->GetDerivative(inputx)*weight)/totalweight;
    }
    m_weight = totalweight;
  }
  else 
    assert(false);
}
